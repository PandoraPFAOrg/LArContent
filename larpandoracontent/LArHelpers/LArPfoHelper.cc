/**
 *  @file   larpandoracontent/LArHelpers/LArPfoHelper.cc
 *
 *  @brief  Implementation of the pfo helper class.
 *
 *  $Log: $
 */

#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArObjectHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

void LArPfoHelper::GetCoordinateVector(const ParticleFlowObject *const pPfo, const HitType &hitType, CartesianPointVector &coordinateVector)
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    for (const Cluster *const pCluster : clusterList)
        LArClusterHelper::GetCoordinateVector(pCluster, coordinateVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetCaloHits(const PfoList &pfoList, const HitType &hitType, CaloHitList &caloHitList)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
        LArPfoHelper::GetCaloHits(pPfo, hitType, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetCaloHits(const ParticleFlowObject *const pPfo, const HitType &hitType, CaloHitList &caloHitList)
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    for (const Cluster *const pCluster : clusterList)
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetIsolatedCaloHits(const PfoList &pfoList, const HitType &hitType, CaloHitList &caloHitList)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
        LArPfoHelper::GetIsolatedCaloHits(pPfo, hitType, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetIsolatedCaloHits(const ParticleFlowObject *const pPfo, const HitType &hitType, CaloHitList &caloHitList)
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    for (const Cluster *const pCluster : clusterList)
        caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllCaloHits2D(const ParticleFlowObject *const pPfo, CaloHitList &caloHitList)
{
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_W, caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetClusters(const PfoList &pfoList, const HitType &hitType, ClusterList &clusterList)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
        LArPfoHelper::GetClusters(pPfo, hitType, clusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetClusters(const ParticleFlowObject *const pPfo, const HitType &hitType, ClusterList &clusterList)
{
    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (hitType != LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        clusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArPfoHelper::GetNumberOfTwoDHits(const ParticleFlowObject *const pPfo)
{
    unsigned int totalHits(0);

    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
    for (const Cluster *const pCluster : clusterList)
        totalHits += pCluster->GetNCaloHits();

    return totalHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArPfoHelper::GetNumberOfThreeDHits(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    return total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetTwoDClusterList(const ParticleFlowObject *const pPfo, ClusterList &clusterList)
{
    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        clusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetThreeDClusterList(const ParticleFlowObject *const pPfo, ClusterList &clusterList)
{
    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        clusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllConnectedPfos(const PfoList &inputPfoList, PfoList &outputPfoList)
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
        LArPfoHelper::GetAllConnectedPfos(pPfo, outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllConnectedPfos(const ParticleFlowObject *const pPfo, PfoList &outputPfoList)
{
    if (outputPfoList.end() != std::find(outputPfoList.begin(), outputPfoList.end(), pPfo))
        return;

    outputPfoList.push_back(pPfo);
    LArPfoHelper::GetAllConnectedPfos(pPfo->GetParentPfoList(), outputPfoList);
    LArPfoHelper::GetAllConnectedPfos(pPfo->GetDaughterPfoList(), outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllDownstreamPfos(const PfoList &inputPfoList, PfoList &outputPfoList)
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
        LArPfoHelper::GetAllDownstreamPfos(pPfo, outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllDownstreamPfos(const ParticleFlowObject *const pPfo, PfoList &outputPfoList)
{
    if (outputPfoList.end() != std::find(outputPfoList.begin(), outputPfoList.end(), pPfo))
        return;

    outputPfoList.push_back(pPfo);
    LArPfoHelper::GetAllDownstreamPfos(pPfo->GetDaughterPfoList(), outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllDownstreamPfos(
    const pandora::ParticleFlowObject *const pPfo, pandora::PfoList &outputTrackPfoList, pandora::PfoList &outputLeadingShowerPfoList)
{
    if (LArPfoHelper::IsTrack(pPfo))
    {
        outputTrackPfoList.emplace_back(pPfo);
        for (const ParticleFlowObject *pChild : pPfo->GetDaughterPfoList())
        {
            if (std::find(outputTrackPfoList.begin(), outputTrackPfoList.end(), pChild) == outputTrackPfoList.end())
            {
                const int pdg{std::abs(pChild->GetParticleId())};
                if (pdg == E_MINUS)
                {
                    outputLeadingShowerPfoList.emplace_back(pChild);
                }
                else
                {
                    outputTrackPfoList.emplace_back(pChild);
                    LArPfoHelper::GetAllDownstreamPfos(pChild, outputTrackPfoList, outputLeadingShowerPfoList);
                }
            }
        }
    }
    else
    {
        outputLeadingShowerPfoList.emplace_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArPfoHelper::GetHierarchyTier(const ParticleFlowObject *const pPfo)
{
    const ParticleFlowObject *pParentPfo = pPfo;
    int tier(0);

    while (pParentPfo->GetParentPfoList().empty() == false)
    {
        if (1 != pParentPfo->GetParentPfoList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentPfo = *(pParentPfo->GetParentPfoList().begin());
        ++tier;
    }

    return tier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetTwoDLengthSquared(const ParticleFlowObject *const pPfo)
{
    if (!LArPfoHelper::IsTwoD(pPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float lengthSquared(0.f);

    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        lengthSquared += LArClusterHelper::GetLengthSquared(pCluster);
    }

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetThreeDLengthSquared(const ParticleFlowObject *const pPfo)
{
    if (!LArPfoHelper::IsThreeD(pPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float lengthSquared(0.f);

    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        lengthSquared += LArClusterHelper::GetLengthSquared(pCluster);
    }

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetClosestDistance(const ParticleFlowObject *const pPfo, const Cluster *const pCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float bestDistance(std::numeric_limits<float>::max());

    for (const Cluster *const pPfoCluster : clusterList)
    {
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster, pPfoCluster));

        if (thisDistance < bestDistance)
            bestDistance = thisDistance;
    }

    return bestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetThreeDSeparation(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2)
{
    ClusterList clusterList1, clusterList2;

    LArPfoHelper::GetClusters(pPfo1, TPC_3D, clusterList1);
    LArPfoHelper::GetClusters(pPfo2, TPC_3D, clusterList2);

    if (clusterList1.empty() || clusterList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return LArClusterHelper::GetClosestDistance(clusterList1, clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsTwoD(const ParticleFlowObject *const pPfo)
{
    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsThreeD(const ParticleFlowObject *const pPfo)
{
    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsTrack(const ParticleFlowObject *const pPfo)
{
    const int pdg(pPfo->GetParticleId());

    // muon, pion, proton, kaon
    return ((MU_MINUS == std::abs(pdg)) || (PI_PLUS == std::abs(pdg)) || (PROTON == std::abs(pdg)) || (K_PLUS == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsShower(const ParticleFlowObject *const pPfo)
{
    const int pdg(pPfo->GetParticleId());

    // electron, photon
    return ((E_MINUS == std::abs(pdg)) || (PHOTON == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetTrackScore(const ParticleFlowObject *const pPfo)
{
    float trackScore(-999.f);

    const PropertiesMap &metadata(pPfo->GetPropertiesMap());

    if (metadata.find("TrackScore") != metadata.end())
        trackScore = metadata.at("TrackScore");

    return trackScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArPfoHelper::GetPrimaryNeutrino(const ParticleFlowObject *const pPfo)
{
    try
    {
        const ParticleFlowObject *const pParentPfo = LArPfoHelper::GetParentNeutrino(pPfo);
        return pParentPfo->GetParticleId();
    }
    catch (const StatusCodeException &)
    {
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsFinalState(const ParticleFlowObject *const pPfo)
{
    if (pPfo->GetParentPfoList().empty() && !LArPfoHelper::IsNeutrino(pPfo))
        return true;

    if (LArPfoHelper::IsNeutrinoFinalState(pPfo))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsNeutrinoFinalState(const ParticleFlowObject *const pPfo)
{
    return ((pPfo->GetParentPfoList().size() == 1) && (LArPfoHelper::IsNeutrino(*(pPfo->GetParentPfoList().begin()))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsNeutrino(const ParticleFlowObject *const pPfo)
{
    const int absoluteParticleId(std::abs(pPfo->GetParticleId()));

    if ((NU_E == absoluteParticleId) || (NU_MU == absoluteParticleId) || (NU_TAU == absoluteParticleId))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsTestBeamFinalState(const ParticleFlowObject *const pPfo)
{
    return LArPfoHelper::IsTestBeam(LArPfoHelper::GetParentPfo(pPfo));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsTestBeam(const ParticleFlowObject *const pPfo)
{
    const PropertiesMap &properties(pPfo->GetPropertiesMap());
    const PropertiesMap::const_iterator iter(properties.find("IsTestBeam"));

    if (iter != properties.end())
        return ((iter->second > 0.f) ? true : false);

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetRecoNeutrinos(const PfoList *const pPfoList, PfoList &recoNeutrinos)
{
    if (!pPfoList)
        return;

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
            recoNeutrinos.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *LArPfoHelper::GetParentPfo(const ParticleFlowObject *const pPfo)
{
    const ParticleFlowObject *pParentPfo = pPfo;

    while (pParentPfo->GetParentPfoList().empty() == false)
    {
        pParentPfo = *(pParentPfo->GetParentPfoList().begin());
    }

    return pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *LArPfoHelper::GetParentNeutrino(const ParticleFlowObject *const pPfo)
{
    const ParticleFlowObject *const pParentPfo = LArPfoHelper::GetParentPfo(pPfo);

    if (!LArPfoHelper::IsNeutrino(pParentPfo))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *LArPfoHelper::GetVertex(const ParticleFlowObject *const pPfo)
{
    if (pPfo->GetVertexList().empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const Vertex *pVertex(nullptr);

    // ATTN : Test beam parent pfos contain an interaction and start vertex
    if (LArPfoHelper::IsTestBeam(pPfo) && pPfo->GetParentPfoList().empty())
    {
        pVertex = LArPfoHelper::GetVertexWithLabel(pPfo->GetVertexList(), VERTEX_START);
    }
    else
    {
        if (pPfo->GetVertexList().size() != 1)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        pVertex = *(pPfo->GetVertexList().begin());
    }

    if (VERTEX_3D != pVertex->GetVertexType())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return pVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *LArPfoHelper::GetTestBeamInteractionVertex(const ParticleFlowObject *const pPfo)
{
    if (pPfo->GetVertexList().empty() || !pPfo->GetParentPfoList().empty() || !LArPfoHelper::IsTestBeam(pPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const Vertex *pInteractionVertex(LArPfoHelper::GetVertexWithLabel(pPfo->GetVertexList(), VERTEX_INTERACTION));

    if (VERTEX_3D != pInteractionVertex->GetVertexType())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return pInteractionVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *LArPfoHelper::GetVertexWithLabel(const VertexList &vertexList, const VertexLabel vertexLabel)
{
    const Vertex *pTargetVertex(nullptr);

    for (const Vertex *pCandidateVertex : vertexList)
    {
        if (pCandidateVertex->GetVertexLabel() == vertexLabel)
        {
            if (!pTargetVertex)
            {
                pTargetVertex = pCandidateVertex;
            }
            else
            {
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
    }

    if (!pTargetVertex)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pTargetVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetSlidingFitTrajectory(const CartesianPointVector &pointVector, const CartesianVector &vertexPosition,
    const unsigned int layerWindow, const float layerPitch, LArTrackStateVector &trackStateVector, IntVector *const pIndexVector)
{
    LArPfoHelper::SlidingFitTrajectoryImpl(&pointVector, vertexPosition, layerWindow, layerPitch, trackStateVector, pIndexVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetSlidingFitTrajectory(const ParticleFlowObject *const pPfo, const Vertex *const pVertex,
    const unsigned int layerWindow, const float layerPitch, LArTrackStateVector &trackStateVector)
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList);
    LArPfoHelper::SlidingFitTrajectoryImpl(&caloHitList, pVertex->GetPosition(), layerWindow, layerPitch, trackStateVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetSlidingFitTrajectory(const CaloHitList *const pCaloHitList, const CartesianVector &vertexPosition, const unsigned int layerWindow,
    const float layerPitch, LArTrackStateVector &trackStateVector, IntVector *const pIndexVector, const bool return3DCaloHit)
{
    LArPfoHelper::SlidingFitTrajectoryImpl(pCaloHitList, vertexPosition, layerWindow, layerPitch, trackStateVector, pIndexVector, return3DCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArShowerPCA LArPfoHelper::GetPrincipalComponents(const CartesianPointVector &pointVector, const CartesianVector &vertexPosition)
{
    // Run the PCA analysis
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(pointVector, centroid, eigenValues, eigenVecs);

    // Require that principal eigenvalue should always be positive
    if (eigenValues.GetX() < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // By convention, principal axis should always point away from vertex
    const float testProjection(eigenVecs.at(0).GetDotProduct(vertexPosition - centroid));
    const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

    const CartesianVector primaryAxis(eigenVecs.at(0) * directionScaleFactor);
    const CartesianVector secondaryAxis(eigenVecs.at(1) * directionScaleFactor);
    const CartesianVector tertiaryAxis(eigenVecs.at(2) * directionScaleFactor);

    return LArShowerPCA(centroid, primaryAxis, secondaryAxis, tertiaryAxis, eigenValues);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArShowerPCA LArPfoHelper::GetPrincipalComponents(const ParticleFlowObject *const pPfo, const Vertex *const pVertex)
{
    CartesianPointVector pointVector;
    LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, pointVector);
    return LArPfoHelper::GetPrincipalComponents(pointVector, pVertex->GetPosition());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::SortByHitProjection(const LArTrackTrajectoryPoint &lhs, const LArTrackTrajectoryPoint &rhs)
{
    if (lhs.first != rhs.first)
        return (lhs.first < rhs.first);

    // ATTN Removed to support use with CartesianVector only (no CaloHit) input
    // if (lhs.second.GetCaloHit() && rhs.second.GetCaloHit())
    //     return (lhs.second.GetCaloHit()->GetInputEnergy() > rhs.second.GetCaloHit()->GetInputEnergy());

    const float dx(lhs.second.GetPosition().GetX() - rhs.second.GetPosition().GetX());
    const float dy(lhs.second.GetPosition().GetY() - rhs.second.GetPosition().GetY());
    const float dz(lhs.second.GetPosition().GetZ() - rhs.second.GetPosition().GetZ());
    return (dx + dy + dz > 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::SortByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    unsigned int nTwoDHitsLhs(0), nThreeDHitsLhs(0);
    float energyLhs(0.f);
    for (ClusterList::const_iterator iter = pLhs->GetClusterList().begin(), iterEnd = pLhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pClusterLhs = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pClusterLhs))
            nTwoDHitsLhs += pClusterLhs->GetNCaloHits();
        else
            nThreeDHitsLhs += pClusterLhs->GetNCaloHits();

        energyLhs += pClusterLhs->GetHadronicEnergy();
    }

    unsigned int nTwoDHitsRhs(0), nThreeDHitsRhs(0);
    float energyRhs(0.f);
    for (ClusterList::const_iterator iter = pRhs->GetClusterList().begin(), iterEnd = pRhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pClusterRhs = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pClusterRhs))
            nTwoDHitsRhs += pClusterRhs->GetNCaloHits();
        else
            nThreeDHitsRhs += pClusterRhs->GetNCaloHits();

        energyRhs += pClusterRhs->GetHadronicEnergy();
    }

    if (nTwoDHitsLhs != nTwoDHitsRhs)
        return (nTwoDHitsLhs > nTwoDHitsRhs);

    if (nThreeDHitsLhs != nThreeDHitsRhs)
        return (nThreeDHitsLhs > nThreeDHitsRhs);

    // ATTN Need an efficient (balance with well-motivated) tie-breaker here. Pfo length, for instance, is extremely slow.
    return (energyLhs > energyRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetBreadthFirstHierarchyRepresentation(const pandora::ParticleFlowObject *const pPfo, pandora::PfoList &pfoList)
{
    const ParticleFlowObject *pRoot{pPfo};
    PfoList parents{pRoot->GetParentPfoList()};
    while (!parents.empty())
    {
        if (parents.size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        pRoot = parents.front();
        parents = pRoot->GetParentPfoList();
    }
    PfoList queue;
    pfoList.emplace_back(pRoot);
    queue.emplace_back(pRoot);

    while (!queue.empty())
    {
        const PfoList &daughters{queue.front()->GetDaughterPfoList()};
        queue.pop_front();
        for (const ParticleFlowObject *pDaughter : daughters)
        {
            pfoList.emplace_back(pDaughter);
            queue.emplace_back(pDaughter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void LArPfoHelper::SlidingFitTrajectoryImpl(const T *const pT, const CartesianVector &vertexPosition, const unsigned int layerWindow,
    const float layerPitch, LArTrackStateVector &trackStateVector, IntVector *const pIndexVector, const bool return3DCaloHit)
{
    CartesianPointVector pointVector;

    for (const auto &nextPoint : *pT)
        pointVector.push_back(LArObjectHelper::TypeAdaptor::GetPosition(nextPoint));

    if (pointVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    std::sort(pointVector.begin(), pointVector.end(), LArClusterHelper::SortCoordinatesByPosition);

    LArTrackTrajectory trackTrajectory;
    IntVector indicesWithoutSpacePoints;
    if (pIndexVector)
        pIndexVector->clear();

    try
    {
        const ThreeDSlidingFitResult slidingFitResult(&pointVector, layerWindow, layerPitch);
        const CartesianVector minPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector maxPosition(slidingFitResult.GetGlobalMaxLayerPosition());

        if ((maxPosition - minPosition).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const CartesianVector seedPosition((maxPosition + minPosition) * 0.5f);
        const CartesianVector seedDirection((maxPosition - minPosition).GetUnitVector());

        const float scaleFactor((seedDirection.GetDotProduct(seedPosition - vertexPosition) > 0.f) ? +1.f : -1.f);

        int index(-1);
        for (const auto &nextPoint : *pT)
        {
            ++index;

            try
            {
                const float rL(slidingFitResult.GetLongitudinalDisplacement(LArObjectHelper::TypeAdaptor::GetPosition(nextPoint)));

                CartesianVector position(0.f, 0.f, 0.f);
                const StatusCode positionStatusCode(slidingFitResult.GetGlobalFitPosition(rL, position));

                if (positionStatusCode != STATUS_CODE_SUCCESS)
                    throw StatusCodeException(positionStatusCode);

                CartesianVector direction(0.f, 0.f, 0.f);
                const StatusCode directionStatusCode(slidingFitResult.GetGlobalFitDirection(rL, direction));

                if (directionStatusCode != STATUS_CODE_SUCCESS)
                    throw StatusCodeException(directionStatusCode);

                const float projection(seedDirection.GetDotProduct(position - seedPosition));

                trackTrajectory.push_back(LArTrackTrajectoryPoint(projection * scaleFactor,
                    LArTrackState(position, direction * scaleFactor, LArObjectHelper::TypeAdaptor::GetCaloHit(nextPoint, return3DCaloHit)), index));
            }
            catch (const StatusCodeException &statusCodeException1)
            {
                indicesWithoutSpacePoints.push_back(index);

                if (STATUS_CODE_FAILURE == statusCodeException1.GetStatusCode())
                    throw statusCodeException1;
            }
        }
    }
    catch (const StatusCodeException &statusCodeException2)
    {
        if (STATUS_CODE_FAILURE == statusCodeException2.GetStatusCode())
            throw statusCodeException2;
    }

    // Sort trajectory points by distance along track
    std::sort(trackTrajectory.begin(), trackTrajectory.end(), LArPfoHelper::SortByHitProjection);

    for (const LArTrackTrajectoryPoint &larTrackTrajectoryPoint : trackTrajectory)
    {
        trackStateVector.push_back(larTrackTrajectoryPoint.second);
        if (pIndexVector)
            pIndexVector->push_back(larTrackTrajectoryPoint.GetIndex());
    }

    // Store indices of spacepoints with no associated trajectory point at the end of pIndexVector
    if (pIndexVector)
    {
        for (const int index : indicesWithoutSpacePoints)
            pIndexVector->push_back(index);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template void LArPfoHelper::SlidingFitTrajectoryImpl(const CartesianPointVector *const, const CartesianVector &, const unsigned int,
    const float, LArTrackStateVector &, IntVector *const, const bool);
template void LArPfoHelper::SlidingFitTrajectoryImpl(
    const CaloHitList *const, const CartesianVector &, const unsigned int, const float, LArTrackStateVector &, IntVector *const, const bool);

} // namespace lar_content
