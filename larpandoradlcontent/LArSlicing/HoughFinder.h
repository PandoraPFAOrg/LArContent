/**
 * @file   larpandoradlcontent/LArSlicing/HoughFinder.h
 * @brief  Header file for the Fast Hough Transform vertex finder
 */
#ifndef LAR_FAST_HOUGH_FINDER_H
#define LAR_FAST_HOUGH_FINDER_H 1

#include "Objects/CartesianVector.h"
#include <vector>

namespace lar_content
{

class ThreeDVertexingHoughTransform
{
public:
    /**
     * @brief Constructor
     *
     * @param thresholds List of distance bin edge thresholds
     * @param scalingFactor The scaling factor to convert network output back to physical cm
     * @param tolerance The distance tolerance (in cm) for a vote to count
     * @param minVotes Minimum number of votes required to keep a vertex seed
     * @param nmsRadius Radius (in cm) for Non-Maximum Suppression (NMS)
     * @param maxSeedClass Maximum distance class to consider as a seed for vertex finding
     * @param useDynamicSeedClass If true, dynamically finds the closest available class to use as seeds
     */
    ThreeDVertexingHoughTransform(const std::vector<float> &thresholds, const float scalingFactor = 400.0f, const float tolerance = 25.0f,
        const int minVotes = 3, const float nmsRadius = 35.0f, const int maxSeedClass = 2,
        const bool useDynamicSeedClass = false);

    /**
     * @brief Finds vertices based on hit positions and network distance logits
     *
     * @param hitPositions Vector of physical hit positions
     * @param logits Vector of network logits for each hit (flattened, size = N_hits * N_classes)
     *
     * @return Vector of found vertex positions
     */
    std::vector<pandora::CartesianVector> Fit(const std::vector<pandora::CartesianVector> &hitPositions, const std::vector<float> &logits) const;

private:
    /**
     * @brief Seed Vs Voter Selection Helper
     */
    void GetSeedAndVoterIndices(const std::vector<int> &predClasses, std::vector<int> &seedIndices, std::vector<int> &voterIndices) const;

    /**
     * @brief Score sorting helper
     */
    void CalculateSortScores(const std::vector<int> &seedIndices, const std::vector<int> &predClasses, const std::vector<float> &voteCounts,
        const std::vector<float> &logits, std::vector<float> &sortScores) const;

    std::vector<float> m_thresholds;    ///< The distance thresholds for each class bin edge. Must match the training thresholds.
    std::vector<float> m_binCenters;    ///< The pre-computed bin centers for each class, derived from the thresholds, used to convert class predictions to distance predictions.
    float m_scalingFactor;              ///< The scaling factor to convert the network output back to physical units, if needed.
    float m_tolerance;                  ///< The distance tolerance for a vote to count for a given seed.
    int m_minVotes;                     ///< The minimum number of votes required for a seed to be kept as a vertex.
    int m_maxSeedClass;                 ///< The maximum distance class to consider as a seed for vertex finding.
    bool m_useDynamicSeedClass;         ///< Whether to dynamically find the closest available class to use as seeds, if no seeds are found with the specified maxSeedClass.

    float m_nmsRadiusSq;                ///< Stored as squared radius for faster comparison
};

} // namespace lar_content

#endif // LAR_FAST_HOUGH_FINDER_H
