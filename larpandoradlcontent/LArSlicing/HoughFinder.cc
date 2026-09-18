/**
 * @file   larpandoradlcontent/LArSlicing/HoughFinder.cc
 * @brief  Implementation of the Fast Hough Transform vertex finder
 */

#include <algorithm>
#include <cmath>
#include <numeric>

#include "larpandoradlcontent/LArSlicing/HoughFinder.h"

namespace lar_content
{

ThreeDVertexingHoughTransform::ThreeDVertexingHoughTransform(const std::vector<float> &thresholds, const float scalingFactor, const float tolerance, const int minVotes,
    const float nmsRadius, const int maxSeedClass, const bool useDynamicSeedClass) :
    m_thresholds(thresholds),
    m_scalingFactor(scalingFactor),
    m_tolerance(tolerance),
    m_minVotes(minVotes),
    m_maxSeedClass(maxSeedClass),
    m_nmsRadiusSq(nmsRadius * nmsRadius),
    m_useDynamicSeedClass(useDynamicSeedClass)
{
    float prev_t = 0.0f;
    for (const float t : m_thresholds)
    {
        m_binCenters.push_back((prev_t + t) / 2.0f);
        prev_t = t;
    }
    m_binCenters.push_back(0.0f);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexingHoughTransform::GetSeedAndVoterIndices(const std::vector<int> &predClasses, std::vector<int> &seedIndices, std::vector<int> &voterIndices) const
{
    const int numHits = predClasses.size();
    const int noiseClass = m_thresholds.size();

    int targetSeedClass = m_maxSeedClass;

    // Reserve the sizes for the vectors to avoid multiple allocations
    seedIndices.reserve(numHits);
    voterIndices.reserve(numHits);

    // Optional: Dynamically find the absolute lowest available class in this point cloud.
    // In some instances...there many not be a class 0 or 1 etc (i.e. true vertex is in gap)
    // so this should let it go as low as it needs to find seeds.
    //
    // TODO: Should there be a cut off? If m_maxSeedClass + X break and don't run?
    if (m_useDynamicSeedClass)
    {
        targetSeedClass = noiseClass - 1;
        for (int i = 0; i < numHits; ++i)
        {
            if (predClasses[i] < targetSeedClass)
                targetSeedClass = predClasses[i];
        }
    }

    for (int i = 0; i < numHits; ++i)
    {
        const int bestClass = predClasses[i];

        if (m_useDynamicSeedClass)
        {
            if (bestClass == targetSeedClass)
                seedIndices.push_back(i);
        }
        else
        {
            if (bestClass <= m_maxSeedClass)
                seedIndices.push_back(i);
        }

        // The last class is noise, so we don't want those voting.
        if (bestClass > targetSeedClass && bestClass != noiseClass)
            voterIndices.push_back(i);
    }

    // Resize the vectors to the actual number of seeds and voters found
    seedIndices.resize(seedIndices.size());
    voterIndices.resize(voterIndices.size());
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexingHoughTransform::CalculateSortScores(const std::vector<int> &seedIndices, const std::vector<int> &predClasses,
    const std::vector<float> &voteCounts, const std::vector<float> &logits, std::vector<float> &sortScores) const
{
    const int numSeeds = seedIndices.size();
    const int numClasses = m_thresholds.size() + 1;

    // Find the maximum votes across the slice
    float maxVotes = 0.0f;
    for (int s = 0; s < numSeeds; ++s)
        maxVotes = std::max(maxVotes, voteCounts[s]);

    // Define the relative strict threshold (within 3 votes of the max)
    const float strictThreshold = std::max(0.0f, maxVotes - 3.0f);

    for (int s = 0; s < numSeeds; ++s)
    {
        const int seedGlobalIdx = seedIndices[s];
        const int seedClass = predClasses[seedGlobalIdx];

        // Calculate Confidence
        const int logitOffset = seedGlobalIdx * numClasses;
        float maxLogit = -std::numeric_limits<float>::max();

        for (int c = 0; c < numClasses; ++c)
            maxLogit = std::max(maxLogit, logits[logitOffset + c]);

        float sumExp = 0.0f;
        for (int c = 0; c < numClasses; ++c)
            sumExp += std::exp(logits[logitOffset + c] - maxLogit);

        const float confidence = std::exp(logits[logitOffset + seedClass] - maxLogit) / sumExp;

        float score = ((m_maxSeedClass + 1) - seedClass) * 10000.0f;

        if (voteCounts[s] >= strictThreshold && voteCounts[s] >= 5.0f)
            score += 1000.0f + (confidence * 100.0f);
        else
            score += voteCounts[s];

        sortScores[s] = score;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<pandora::CartesianVector> ThreeDVertexingHoughTransform::Fit(const std::vector<pandora::CartesianVector> &hitPositions, const std::vector<float> &logits) const
{
    const int numHits = hitPositions.size();
    const int numClasses = m_thresholds.size() + 1;

    if (numHits == 0)
        return {};

    std::vector<float> predDists(numHits);
    std::vector<int> predClasses(numHits);
    std::vector<int> seedIndices;
    std::vector<int> voterIndices;

    for (int i = 0; i < numHits; ++i)
    {
        int bestClass = 0;
        float maxLogit = logits[i * numClasses];
        for (int c = 1; c < numClasses; ++c)
        {
            const float val = logits[i * numClasses + c];
            if (val > maxLogit)
            {
                maxLogit = val;
                bestClass = c;
            }
        }

        predDists[i] = m_binCenters[bestClass] * m_scalingFactor;
        predClasses[i] = bestClass;
    }

    // Separate hits into seeds and voters.
    this->GetSeedAndVoterIndices(predClasses, seedIndices, voterIndices);

    if (seedIndices.empty())
        return {};

    const int numSeeds = seedIndices.size();
    const int numVoters = voterIndices.size();
    std::vector<float> voteCounts(numSeeds, 0);
    std::vector<float> sortScores(numSeeds, 0);

    // Simple vectors for positions, to help with cache locality.
    std::vector<float> vX(numVoters), vY(numVoters), vZ(numVoters);
    std::vector<float> voterLowerBoundSq(numVoters), voterUpperBoundSq(numVoters);

    // Also store a weighting for each voter, based on the local density of hits around it.
    std::vector<float> voterWeights(numVoters, 1.0f);
    const float densityRadiusSq = 5.0f * 5.0f; // TODO: Param?

    // Precompute the voter positions and their distance bounds for voting.
    for (int v = 0; v < numVoters; ++v)
    {
        const int idx = voterIndices[v];
        const pandora::CartesianVector &pos = hitPositions[idx];

        vX[v] = pos.GetX();
        vY[v] = pos.GetY();
        vZ[v] = pos.GetZ();

        const float pd = predDists[idx];
        const float lb = std::max(0.0f, pd - m_tolerance);
        const float ub = pd + m_tolerance;
        voterLowerBoundSq[v] = lb * lb;
        voterUpperBoundSq[v] = ub * ub;

        int neighbors = 0;
        for (int v2 = 0; v2 < numVoters; ++v2)
        {
            const float dx = vX[v] - vX[v2];
            const float dy = vY[v] - vY[v2];
            const float dz = vZ[v] - vZ[v2];
            if ((dx * dx + dy * dy + dz * dz) < densityRadiusSq)
                neighbors++;
        }

        voterWeights[v] = 1.0f / static_cast<float>(std::max(1, neighbors));
    }

    // Voting loop - for each candidate, count how many voters are within the
    // predicted distance (with tolerance).
    for (int s = 0; s < numSeeds; ++s)
    {
        const int seedIdx = seedIndices[s];
        const float cX = hitPositions[seedIdx].GetX();
        const float cY = hitPositions[seedIdx].GetY();
        const float cZ = hitPositions[seedIdx].GetZ();

        float votes = 0.f;

        for (int v = 0; v < numVoters; ++v)
        {
            const float dx = cX - vX[v];
            const float dy = cY - vY[v];
            const float dz = cZ - vZ[v];
            const float distSq = (dx * dx) + (dy * dy) + (dz * dz);

            if (distSq >= voterLowerBoundSq[v] && distSq <= voterUpperBoundSq[v])
                votes += voterWeights[v];
        }
        voteCounts[s] = votes;
    }

    // Calculate sorting scores, to pick the best candidates first in NMS.
    this->CalculateSortScores(seedIndices, predClasses, voteCounts, logits, sortScores);

    std::vector<int> sortedCandIndices(numSeeds);
    std::iota(sortedCandIndices.begin(), sortedCandIndices.end(), 0);
    std::sort(sortedCandIndices.begin(), sortedCandIndices.end(), [&sortScores](int a, int b) { return sortScores[a] > sortScores[b]; });

    std::vector<pandora::CartesianVector> foundVertices;
    std::vector<bool> seedIsAvailable(numSeeds, true);
    std::vector<bool> voterIsAvailable(numVoters, true);

    // Non-Maximum Suppression Loop then final vertex creation
    for (const int candListIdx : sortedCandIndices)
    {
        if (!seedIsAvailable[candListIdx])
            continue;

        if (voteCounts[candListIdx] < m_minVotes)
            break;

        const int candGlobalIdx = seedIndices[candListIdx];
        const float cX = hitPositions[candGlobalIdx].GetX();
        const float cY = hitPositions[candGlobalIdx].GetY();
        const float cZ = hitPositions[candGlobalIdx].GetZ();

        int currentSupport = 0;
        std::vector<int> claimedVotersLocal;

        for (int v = 0; v < numVoters; ++v)
        {
            if (!voterIsAvailable[v])
                continue;

            const float dx = cX - vX[v];
            const float dy = cY - vY[v];
            const float dz = cZ - vZ[v];
            const float geomDistSq = (dx * dx) + (dy * dy) + (dz * dz);

            if (geomDistSq >= voterLowerBoundSq[v] && geomDistSq <= voterUpperBoundSq[v])
            {
                currentSupport++;
                claimedVotersLocal.push_back(v);
            }
        }

        if (currentSupport < m_minVotes)
            continue;

        const auto candPos = pandora::CartesianVector(cX, cY, cZ);
        foundVertices.push_back(candPos);

        // Consume the voters so they can't vote for subsequent nearby seeds
        for (const int localVoterIdx : claimedVotersLocal)
            voterIsAvailable[localVoterIdx] = false;

        for (int s = 0; s < numSeeds; ++s)
        {
            if (!seedIsAvailable[s])
                continue;

            const pandora::CartesianVector &otherCandPos = hitPositions[seedIndices[s]];
            if ((candPos - otherCandPos).GetMagnitudeSquared() < m_nmsRadiusSq)
            {
                seedIsAvailable[s] = false;
            }
        }
    }

    // INFO: Fallback logic in dynamic mode - Just use the best option.
    if (foundVertices.empty() && m_useDynamicSeedClass && numSeeds > 0)
    {
        const int bestCandListIdx = sortedCandIndices[0];
        const int bestGlobalIdx = seedIndices[bestCandListIdx];

        foundVertices.push_back(hitPositions[bestGlobalIdx]);
    }

    return foundVertices;
}

} // namespace lar_content
