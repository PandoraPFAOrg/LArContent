/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <cmath>

using namespace pandora;

namespace lar_content
{

CrossGapsAssociationAlgorithm::CrossGapsAssociationAlgorithm() :
    m_minClusterHits(10),
    m_minClusterLayers(6),
    m_slidingFitWindow(20),
    m_maxSamplingPoints(1000),
    m_sampleStepSize(0.5f),
    m_maxUnmatchedSampleRun(8),
    m_maxOnClusterDistance(1.5f),
    m_minMatchedSamplingPoints(10),
    m_minMatchedSamplingFraction(0.5f),
    m_crossTPCStepModifier(0.0f),
    m_crossTPCOnClusterDistanceModifier(0.0f),
    m_boostStartStep(0.0f),
    m_crossTPCMaxAdditionalSteps(1000.0f),
    m_gapTolerance(0.f),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN May want to opt-out completely if no gap information available
    // if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
    //     return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers)   
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    TwoDSlidingFitResultMap slidingFitResultMap;

    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
            slidingFitResultMap.insert(
                TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // ATTN This method assumes that clusters have been sorted by layer
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        const Cluster *const pInnerCluster = *iterI;
        TwoDSlidingFitResultMap::const_iterator fitIterI = slidingFitResultMap.find(pInnerCluster);

        if (slidingFitResultMap.end() == fitIterI)
            continue;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster = *iterJ;

            if (pInnerCluster == pOuterCluster)
                continue;

            TwoDSlidingFitResultMap::const_iterator fitIterJ = slidingFitResultMap.find(pOuterCluster);

            if (slidingFitResultMap.end() == fitIterJ)
                continue;

            if (!this->AreClustersAssociated(fitIterI->second, fitIterJ->second))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::AreClustersAssociated(const TwoDSlidingFitResult &innerFitResult, const TwoDSlidingFitResult &outerFitResult) const
{
    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetInnerPseudoLayer())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetOuterPseudoLayer())
        return false;
        
    bool isCrossingTPCCandidate = false;
    if (m_crossTPCStepModifier != 0.0) {
        
        std::set<unsigned int> innerVolIDs;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindLArTPCVolumeIds(innerFitResult, innerVolIDs));
        std::set<unsigned int> outerVolIDs;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindLArTPCVolumeIds(outerFitResult, outerVolIDs));

        std::vector<int> sharedClusterVolIDs;
    
        std::set_intersection(innerVolIDs.begin(), innerVolIDs.end(),
                          outerVolIDs.begin(), outerVolIDs.end(),
                          std::back_inserter(sharedClusterVolIDs));
    
        if (sharedClusterVolIDs.empty()) 
            isCrossingTPCCandidate = true; //set this to false to turn off the corrections for crossing the gap     
        
    }

    return (this->IsAssociated(innerFitResult.GetGlobalMaxLayerPosition(), innerFitResult.GetGlobalMaxLayerDirection(), outerFitResult, isCrossingTPCCandidate) &&
        this->IsAssociated(outerFitResult.GetGlobalMinLayerPosition(), outerFitResult.GetGlobalMinLayerDirection() * -1.f, innerFitResult, isCrossingTPCCandidate));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsAssociated(
    const CartesianVector &startPosition, const CartesianVector &startDirection, const TwoDSlidingFitResult &targetFitResult, const bool isCrossingTPCCandidate) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), hitType)};
    const float sampleStepSizeAdjusted{ratio * m_sampleStepSize};
    unsigned int nMatchedSamplingPoints(0), nUnmatchedSampleRun(0);

    bool shouldVisualize{false};
    if (m_visualize)
    {
        for (unsigned int iSample = 0; iSample < m_maxSamplingPoints; ++iSample)
        {
            const CartesianVector samplingPoint(startPosition + startDirection * static_cast<float>(iSample) * sampleStepSizeAdjusted);

            if (LArGeometryHelper::IsInGap(this->GetPandora(), samplingPoint, hitType, m_gapTolerance))
            {
                shouldVisualize = true;
                break;
            }
        }
    }
    
    unsigned int crossTPCGapAdditionalSteps = 0;
    if(isCrossingTPCCandidate) {
        if ( startDirection.GetOpeningAngle(CartesianVector(1.0, 0.0, 0.0)) < startDirection.GetOpeningAngle(CartesianVector(-1.0, 0.0, 0.0)) ) {
            crossTPCGapAdditionalSteps = static_cast<int>(std::round( m_crossTPCStepModifier / startDirection.GetCosOpeningAngle(CartesianVector(1.0, 0.0, 0.0))))  *      static_cast<int>(isCrossingTPCCandidate);
        }
        
        else {
            crossTPCGapAdditionalSteps = static_cast<int>(std::round( m_crossTPCStepModifier / startDirection.GetCosOpeningAngle(CartesianVector(-1.0, 0.0, 0.0))))  * static_cast<int>(isCrossingTPCCandidate);
        }
        if ( crossTPCGapAdditionalSteps > m_crossTPCMaxAdditionalSteps ) {
            crossTPCGapAdditionalSteps = m_crossTPCMaxAdditionalSteps;
        }
    }
    float numGapSteps = 0.0;
    for (unsigned int iSample = 0; iSample < m_maxSamplingPoints; ++iSample)
    {
    
        const CartesianVector samplingPoint(startPosition + startDirection * static_cast<float>(iSample) * sampleStepSizeAdjusted);

        if (LArGeometryHelper::IsInGap(this->GetPandora(), samplingPoint, hitType, m_gapTolerance))
        {
            if (shouldVisualize)
            {
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "", BLUE, 1));
            }
            numGapSteps = numGapSteps + 1.0;
            nUnmatchedSampleRun = 0; // ATTN Choose to also reset run when entering gap region
            continue;
        }

        if (this->IsNearCluster(samplingPoint, targetFitResult, numGapSteps, isCrossingTPCCandidate))
        {
            ++nMatchedSamplingPoints;
            nUnmatchedSampleRun = 0;
            if (shouldVisualize)
            {
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "", GREEN, 1));
            }
        }
        else if (++nUnmatchedSampleRun > m_maxUnmatchedSampleRun + crossTPCGapAdditionalSteps)
        {
            break;
        }
        else if (shouldVisualize)
        {
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "", RED, 1));
        }
    }
    
    if (shouldVisualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
        
    const float expectation(
        (targetFitResult.GetGlobalMaxLayerPosition() - targetFitResult.GetGlobalMinLayerPosition()).GetMagnitude() / sampleStepSizeAdjusted);
    const float matchedSamplingFraction(expectation > 0.f ? static_cast<float>(nMatchedSamplingPoints) / expectation : 0.f);

    if ((nMatchedSamplingPoints > m_minMatchedSamplingPoints) || (matchedSamplingFraction > m_minMatchedSamplingFraction))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsNearCluster(const CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult, const float numGapSteps, const bool isCrossingTPCCandidate) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), hitType)};
    const float maxOnClusterDistanceAdjusted{ratio * m_maxOnClusterDistance};

    float rL(std::numeric_limits<float>::max()), rT(std::numeric_limits<float>::max());
    targetFitResult.GetLocalPosition(samplingPoint, rL, rT);

    CartesianVector fitPosition(0.f, 0.f, 0.f);
    
    float additionalMaxOnClusterDistance = 0.0;
    if (isCrossingTPCCandidate) {
        additionalMaxOnClusterDistance = m_crossTPCOnClusterDistanceModifier * sqrt(std::max(numGapSteps - m_boostStartStep, 0.f));
    }

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPosition(rL, fitPosition))
    {
        if ((fitPosition - samplingPoint).GetMagnitudeSquared() < (maxOnClusterDistanceAdjusted + additionalMaxOnClusterDistance) * (maxOnClusterDistanceAdjusted + additionalMaxOnClusterDistance))
            return true;
    }

    CartesianVector fitPositionAtX(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPositionAtX(samplingPoint.GetX(), fitPositionAtX))
    {
        if ((fitPositionAtX - samplingPoint).GetMagnitudeSquared() < (maxOnClusterDistanceAdjusted + additionalMaxOnClusterDistance) * (maxOnClusterDistanceAdjusted + additionalMaxOnClusterDistance))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossGapsAssociationAlgorithm::FindLArTPCVolumeIds(const TwoDSlidingFitResult &fitResult, std::set<unsigned int> &volIds) const
{
    const Cluster *const pCluster{fitResult.GetCluster()};
    CaloHitList caloHits;
    LArClusterHelper::GetAllHits(pCluster, caloHits);
    
    for (const CaloHit *const pCaloHit : caloHits)
    {
        const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        if (pLArCaloHit == nullptr) {
            return STATUS_CODE_NOT_INITIALIZED;
        }
        volIds.insert(pLArCaloHit->GetLArTPCVolumeId());
    }
    
    return STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossGapsAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSamplingPoints", m_maxSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "CrossGapsAssociationAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxUnmatchedSampleRun", m_maxUnmatchedSampleRun));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxOnClusterDistance", m_maxOnClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingFraction", m_minMatchedSamplingFraction));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CrossTPCStepModifier", m_crossTPCStepModifier));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CrossTPCOnClusterDistanceModifier", m_crossTPCOnClusterDistanceModifier));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BoostStartStep", m_boostStartStep));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CrossTPCMaxAdditionalSteps", m_crossTPCMaxAdditionalSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GapTolerance", m_gapTolerance));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
