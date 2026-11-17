/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.cc
 *
 *  @brief  Implementation of the overshoot tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.h"

using namespace pandora;

namespace lar_content
{

OvershootTracksTool::OvershootTracksTool() :
    ThreeDKinkBaseTool(1),
    m_splitMode(true),
    m_maxVertexXSeparation(2.f),
    m_cosThetaCutForKinkSearch(0.94f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::GetIteratorListModifications(
    ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2)
        {
            if (iIter1 == iIter2)
                continue;

            try
            {
                const unsigned int nMatchedSamplingPoints1((*iIter1)->GetOverlapResult().GetNMatchedSamplingPoints());
                const unsigned int nMatchedSamplingPoints2((*iIter2)->GetOverlapResult().GetNMatchedSamplingPoints());
                IteratorList::const_iterator iIterA((nMatchedSamplingPoints1 >= nMatchedSamplingPoints2) ? iIter1 : iIter2);
                IteratorList::const_iterator iIterB((nMatchedSamplingPoints1 >= nMatchedSamplingPoints2) ? iIter2 : iIter1);

                Particle particle(*(*iIterA), *(*iIterB));
                const LArPointingCluster pointingClusterA1(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA1));
                const LArPointingCluster pointingClusterB1(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB1));
                const LArPointingCluster pointingClusterA2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA2));
                const LArPointingCluster pointingClusterB2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB2));

                LArPointingCluster::Vertex vertexA1, vertexB1, vertexA2, vertexB2;
                if (!LArPointingClusterHelper::TryGetClosestVerticesInX(pointingClusterA1, pointingClusterB1, vertexA1, vertexB1) ||
                    !LArPointingClusterHelper::TryGetClosestVerticesInX(pointingClusterA2, pointingClusterB2, vertexA2, vertexB2))
                    continue;

                if (!this->PassesVertexCuts(vertexA1, vertexB1) || !this->PassesVertexCuts(vertexA2, vertexB2))
                    continue;

                this->SetSplitPosition(vertexA1, vertexA2, vertexB1, vertexB2, particle);

                const bool isA1LowestInX(this->IsALowestInX(pointingClusterA1, pointingClusterB1));
                const bool isA2LowestInX(this->IsALowestInX(pointingClusterA2, pointingClusterB2));
                const bool isThreeDKink(m_majorityRulesMode ? true : this->IsThreeDKink(pAlgorithm, particle, isA1LowestInX, isA2LowestInX));

                if (isThreeDKink != m_splitMode)
                    continue;

                // Construct the modification object
                Modification modification;

                if (m_splitMode)
                {
                    modification.m_splitPositionMap[particle.m_pCommonCluster].push_back(particle.m_splitPosition);
                }
                else
                {
                    const bool vertex1AIsLowX(vertexA1.GetPosition().GetX() < vertexB1.GetPosition().GetX());
                    const Cluster *const pLowXCluster1(vertex1AIsLowX ? particle.m_pClusterA1 : particle.m_pClusterB1);
                    const Cluster *const pHighXCluster1(vertex1AIsLowX ? particle.m_pClusterB1 : particle.m_pClusterA1);
                    modification.m_clusterMergeMap[pLowXCluster1].push_back(pHighXCluster1);

                    const bool vertex2AIsLowX(vertexA2.GetPosition().GetX() < vertexB2.GetPosition().GetX());
                    const Cluster *const pLowXCluster2(vertex2AIsLowX ? particle.m_pClusterA2 : particle.m_pClusterB2);
                    const Cluster *const pHighXCluster2(vertex2AIsLowX ? particle.m_pClusterB2 : particle.m_pClusterA2);
                    modification.m_clusterMergeMap[pLowXCluster2].push_back(pHighXCluster2);
                }

                modification.m_affectedClusters.push_back(particle.m_pCommonCluster);
                modification.m_affectedClusters.push_back(particle.m_pClusterA1);
                modification.m_affectedClusters.push_back(particle.m_pClusterA2);
                modification.m_affectedClusters.push_back(particle.m_pClusterB1);
                modification.m_affectedClusters.push_back(particle.m_pClusterB2);

                modificationList.push_back(modification);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;

                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OvershootTracksTool::PassesVertexCuts(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB) const
{
    float longitudinalAB(-std::numeric_limits<float>::max()), transverseAB(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetImpactParameters(vertexA, vertexB, longitudinalAB, transverseAB);

    float longitudinalBA(-std::numeric_limits<float>::max()), transverseBA(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetImpactParameters(vertexB, vertexA, longitudinalBA, transverseBA);

    if (std::min(longitudinalAB, longitudinalBA) < m_minLongitudinalImpactParameter)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::SetSplitPosition(const LArPointingCluster::Vertex &vertexA1, const LArPointingCluster::Vertex &vertexA2,
    const LArPointingCluster::Vertex &vertexB1, const LArPointingCluster::Vertex &vertexB2, Particle &particle) const
{
    bool splitAtElementA(false), splitAtElementB(false);

    if (std::fabs(vertexA1.GetPosition().GetX() - vertexA2.GetPosition().GetX()) < m_maxVertexXSeparation)
    {
        splitAtElementA = true;
    }
    else if (std::fabs(vertexB1.GetPosition().GetX() - vertexB2.GetPosition().GetX()) < m_maxVertexXSeparation)
    {
        splitAtElementB = true;
    }

    if (!splitAtElementA && !splitAtElementB)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    particle.m_splitPosition1 = splitAtElementA ? vertexA1.GetPosition() : vertexB1.GetPosition();
    particle.m_splitPosition2 = splitAtElementA ? vertexA2.GetPosition() : vertexB2.GetPosition();

    CartesianVector splitPosition(0.f, 0.f, 0.f);
    float chiSquared(std::numeric_limits<float>::max());
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), LArClusterHelper::GetClusterHitType(particle.m_pClusterA1),
        LArClusterHelper::GetClusterHitType(particle.m_pClusterA2), particle.m_splitPosition1, particle.m_splitPosition2, splitPosition, chiSquared);

    particle.m_splitPosition = splitPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OvershootTracksTool::IsThreeDKink(
    ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle, const bool isA1LowestInX, const bool isA2LowestInX) const
{
    try
    {
        const TwoDSlidingFitResult &lowXFitResult1(isA1LowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA1)
                                                                 : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB1));
        const TwoDSlidingFitResult &highXFitResult1(isA1LowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB1)
                                                                  : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA1));
        const TwoDSlidingFitResult &lowXFitResult2(isA2LowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA2)
                                                                 : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB2));
        const TwoDSlidingFitResult &highXFitResult2(isA2LowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB2)
                                                                  : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA2));
        const TwoDSlidingFitResult &fitResultCommon3(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster));

        const float minusX(this->GetXSamplingPoint(particle.m_splitPosition, false, fitResultCommon3, lowXFitResult1, lowXFitResult2));
        const float plusX(this->GetXSamplingPoint(particle.m_splitPosition, true, fitResultCommon3, highXFitResult1, highXFitResult2));

        CartesianVector minus1(0.f, 0.f, 0.f), split1(particle.m_splitPosition1), plus1(0.f, 0.f, 0.f);
        CartesianVector minus2(0.f, 0.f, 0.f), split2(particle.m_splitPosition2), plus2(0.f, 0.f, 0.f);
        CartesianVector minus3(0.f, 0.f, 0.f), split3(particle.m_splitPosition), plus3(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != lowXFitResult1.GetGlobalFitPositionAtX(minusX, minus1)) ||
            (STATUS_CODE_SUCCESS != highXFitResult1.GetGlobalFitPositionAtX(plusX, plus1)) ||
            (STATUS_CODE_SUCCESS != lowXFitResult2.GetGlobalFitPositionAtX(minusX, minus2)) ||
            (STATUS_CODE_SUCCESS != highXFitResult2.GetGlobalFitPositionAtX(plusX, plus2)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon3.GetGlobalFitPositionAtX(minusX, minus3)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon3.GetGlobalFitPositionAtX(plusX, plus3)))
        {
            return true; // majority rules, by default
        }

        // Extract results
        const HitType hitType1(LArClusterHelper::GetClusterHitType(particle.m_pClusterA1));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(particle.m_pClusterA2));
        const HitType hitType3(LArClusterHelper::GetClusterHitType(particle.m_pCommonCluster));

        CartesianVector minus(0.f, 0.f, 0.f), split(0.f, 0.f, 0.f), plus(0.f, 0.f, 0.f);
        float chi2Minus(std::numeric_limits<float>::max()), chi2Split(std::numeric_limits<float>::max()),
            chi2Plus(std::numeric_limits<float>::max());
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, minus1, minus2, minus3, minus, chi2Minus);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, split1, split2, split3, split, chi2Split);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, plus1, plus2, plus3, plus, chi2Plus);

        // Apply final cuts
        const CartesianVector minusToSplit((split - minus).GetUnitVector());
        const CartesianVector splitToPlus((plus - split).GetUnitVector());
        const float dotProduct(minusToSplit.GetDotProduct(splitToPlus));

        if (dotProduct > m_cosThetaCutForKinkSearch)
            return false;
    }
    catch (StatusCodeException &)
    {
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

OvershootTracksTool::Particle::Particle(const TensorType::Element &elementA, const TensorType::Element &elementB) :
    m_splitPosition(0.f, 0.f, 0.f),
    m_splitPosition1(0.f, 0.f, 0.f),
    m_splitPosition2(0.f, 0.f, 0.f)
{
    const HitType commonView((elementA.GetClusterU() == elementB.GetClusterU()) ? TPC_VIEW_U
            : (elementA.GetClusterV() == elementB.GetClusterV())                ? TPC_VIEW_V
            : (elementA.GetClusterW() == elementB.GetClusterW())                ? TPC_VIEW_W
                                                                                : HIT_CUSTOM);

    if (HIT_CUSTOM == commonView)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_pCommonCluster = (TPC_VIEW_U == commonView) ? elementA.GetClusterU()
        : (TPC_VIEW_V == commonView)              ? elementA.GetClusterV()
                                                  : elementA.GetClusterW();
    m_pClusterA1 = (TPC_VIEW_U == commonView) ? elementA.GetClusterV() : elementA.GetClusterU();
    m_pClusterA2 = (TPC_VIEW_U == commonView) ? elementA.GetClusterW()
        : (TPC_VIEW_V == commonView)          ? elementA.GetClusterW()
                                              : elementA.GetClusterV();
    m_pClusterB1 = (TPC_VIEW_U == commonView) ? elementB.GetClusterV() : elementB.GetClusterU();
    m_pClusterB2 = (TPC_VIEW_U == commonView) ? elementB.GetClusterW()
        : (TPC_VIEW_V == commonView)          ? elementB.GetClusterW()
                                              : elementB.GetClusterV();

    if ((m_pClusterA1 == m_pClusterB1) || (m_pClusterA2 == m_pClusterB2))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OvershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SplitMode", m_splitMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxVertexXSeparation", m_maxVertexXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CosThetaCutForKinkSearch", m_cosThetaCutForKinkSearch));

    return ThreeDKinkBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
