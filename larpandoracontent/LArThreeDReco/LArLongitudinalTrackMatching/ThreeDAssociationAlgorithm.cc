/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the 3D association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeDAssociationAlgorithm::ThreeDAssociationAlgorithm() :
    m_minClusterLength(1.10f),
    m_shortClusterLength(2.50f),
    m_minCosRelativeAngle(0.94f),
    m_maxGapDistanceSquared(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    m_clusterAttrMap.clear();

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster{ *iter };
        const float clusterLength{ LArClusterHelper::GetLength(pCluster) };

        if( clusterLength < m_minClusterLength )
            continue;

        this->PopulateClusterAttributes(pCluster);

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN This method assumes that clusters have been sorted with LArClusterHelper::SortByInnerPosition
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        const Cluster *const pInnerCluster{ *iterI };

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster{ *iterJ };

            if (pInnerCluster == pOuterCluster)
                continue;

            if (!this->AreClustersAssociated(pInnerCluster, pOuterCluster))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    CartesianVector currentClusterInnerCoordinate(0.f, 0.f, 0.f);
    CartesianVector currentClusterOuterCoordinate(0.f, 0.f, 0.f);
    CartesianVector testClusterInnerCoordinate(0.f, 0.f, 0.f);
    CartesianVector testClusterOuterCoordinate(0.f, 0.f, 0.f);

    if( m_clusterAttrMap.find(pCurrentCluster) == m_clusterAttrMap.end() )
    {
        LArClusterHelper::GetExtremalCoordinates(pCurrentCluster, currentClusterInnerCoordinate, currentClusterOuterCoordinate);
    }
    else
    {
        currentClusterInnerCoordinate = m_clusterAttrMap.at(pCurrentCluster).m_innerCoordinate;
        currentClusterOuterCoordinate = m_clusterAttrMap.at(pCurrentCluster).m_outerCoordinate;
    }
    
    if ( m_clusterAttrMap.find(pTestCluster) == m_clusterAttrMap.end() )
    {
        LArClusterHelper::GetExtremalCoordinates(pTestCluster, testClusterInnerCoordinate, testClusterOuterCoordinate);
    } 
    else 
    {
        testClusterInnerCoordinate = m_clusterAttrMap.at(pTestCluster).m_innerCoordinate;
        testClusterOuterCoordinate = m_clusterAttrMap.at(pTestCluster).m_outerCoordinate;
    }

    const CartesianVector &currentPosition{ isForward ? currentClusterOuterCoordinate : currentClusterInnerCoordinate };
    const CartesianVector &testPosition{ isForward ? testClusterOuterCoordinate : testClusterInnerCoordinate };

    // Find the extremal cluster by comparing the outermost coordinates (forward) or innermost coordinates (backward)
    // Pick the cluster with more hits if coordinates are the same
    if( isForward && ((LArClusterHelper::SortCoordinatesByPosition(currentPosition, testPosition) || 
                      ((currentPosition == testPosition) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster)))))
        return true;

    if( !isForward && ((LArClusterHelper::SortCoordinatesByPosition(testPosition, currentPosition) || 
                      ((currentPosition == testPosition) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster)))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDAssociationAlgorithm::AreClustersAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    if( (m_clusterAttrMap.find(pInnerCluster) == m_clusterAttrMap.end()) || (m_clusterAttrMap.find(pOuterCluster) == m_clusterAttrMap.end()) )
        throw pandora::StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const float closestDistance{ LArClusterHelper::GetClosestDistance(pInnerCluster, pOuterCluster) };

    if( closestDistance > m_maxGapDistanceSquared )
        return false;

    const CartesianVector &innerClusterInnerCoordinate{ m_clusterAttrMap.at(pInnerCluster).m_innerCoordinate };
    const CartesianVector &innerClusterOuterCoordinate{ m_clusterAttrMap.at(pInnerCluster).m_outerCoordinate };
    const CartesianVector &outerClusterInnerCoordinate{ m_clusterAttrMap.at(pOuterCluster).m_innerCoordinate };

    if( !LArClusterHelper::SortCoordinatesByPosition(innerClusterInnerCoordinate, outerClusterInnerCoordinate) )
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    const LArPcaHelper::EigenVectors &innerVectors{ m_clusterAttrMap.at(pInnerCluster).m_eigenVectors };
    const LArPcaHelper::EigenVectors &outerVectors{ m_clusterAttrMap.at(pOuterCluster).m_eigenVectors };

    const CartesianVector innerPrimaryAxis  { innerVectors.at(0).GetZ() > 0.f ? innerVectors.at(0) : innerVectors.at(0) * -1.f };
    const CartesianVector innerSecondaryAxis{ innerVectors.at(1).GetZ() > 0.f ? innerVectors.at(1) : innerVectors.at(1) * -1.f };
    const CartesianVector innerTertiaryAxis { innerVectors.at(2).GetZ() > 0.f ? innerVectors.at(2) : innerVectors.at(2) * -1.f };

    const CartesianVector outerPrimaryAxis  { outerVectors.at(0).GetZ() > 0.f ? outerVectors.at(0) : outerVectors.at(0) * -1.f };
    const CartesianVector outerSecondaryAxis{ outerVectors.at(1).GetZ() > 0.f ? outerVectors.at(1) : outerVectors.at(1) * -1.f };
    const CartesianVector outerTertiaryAxis { outerVectors.at(2).GetZ() > 0.f ? outerVectors.at(2) : outerVectors.at(2) * -1.f };

    const pandora::HitType clusterView = LArClusterHelper::GetClusterHitType(pInnerCluster);

    const float openingAngle{ innerPrimaryAxis.GetCosOpeningAngle(outerPrimaryAxis) };
    const double channelPitchAdjusted{ m_channelPitchTolerance * LArGeometryHelper::GetWirePitch(this->GetPandora(), clusterView) };

    // Bypass opening angle check if a short cluster extends a long one (or vice versa)
    if( (closestDistance > channelPitchAdjusted) ||
        ((m_clusterAttrMap.at(pInnerCluster).m_length <= m_shortClusterLength ) == (m_clusterAttrMap.at(pOuterCluster).m_length <= m_shortClusterLength))){

        if( openingAngle < m_minCosRelativeAngle )
            return false;
    }

    const bool isInnerLongest{ LArClusterHelper::SortByNHits(pInnerCluster, pOuterCluster) && 
        (m_clusterAttrMap.at(pInnerCluster).m_length > m_clusterAttrMap.at(pOuterCluster).m_length) };
    
    const CartesianVector &centroid{ isInnerLongest ? m_clusterAttrMap.at(pInnerCluster).m_centroid : m_clusterAttrMap.at(pOuterCluster).m_centroid };

    const CartesianVector &primaryAxis  { isInnerLongest ? innerPrimaryAxis   : outerPrimaryAxis   };
    const CartesianVector &secondaryAxis{ isInnerLongest ? innerSecondaryAxis : outerSecondaryAxis };
    const CartesianVector &tertiaryAxis { isInnerLongest ? innerTertiaryAxis  : outerTertiaryAxis  };

    const CartesianVector innerEndPrimaryAxis  { centroid + primaryAxis * primaryAxis.GetDotProduct(innerClusterOuterCoordinate - centroid) };
    const CartesianVector outerStartPrimaryAxis{ centroid + primaryAxis * primaryAxis.GetDotProduct(outerClusterInnerCoordinate - centroid) };

    const CartesianVector innerEndSecondaryAxis  { centroid + secondaryAxis * secondaryAxis.GetDotProduct(innerClusterOuterCoordinate - centroid) };
    const CartesianVector outerStartSecondaryAxis{ centroid + secondaryAxis * secondaryAxis.GetDotProduct(outerClusterInnerCoordinate - centroid) };
    
    const CartesianVector innerEndTertiaryAxis  { centroid + tertiaryAxis * tertiaryAxis.GetDotProduct(innerClusterOuterCoordinate - centroid) };
    const CartesianVector outerStartTertiaryAxis{ centroid + tertiaryAxis * tertiaryAxis.GetDotProduct(outerClusterInnerCoordinate - centroid) };

    if( (innerEndTertiaryAxis - outerStartTertiaryAxis).GetMagnitudeSquared() > m_shortClusterLength )
        return false;

    if( (innerEndSecondaryAxis - outerStartSecondaryAxis).GetMagnitudeSquared() > m_shortClusterLength )
        return false;

    if( (innerEndPrimaryAxis - outerStartPrimaryAxis).GetMagnitudeSquared() < m_maxGapDistanceSquared )
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShortClusterLength", m_shortClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapDistanceSquared", m_maxGapDistanceSquared));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChannelPitchTolerance", m_channelPitchTolerance));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::PopulateClusterAttributes( const pandora::Cluster *const pCluster ) const
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);   

    CartesianVector centroid (0.f, 0.f, 0.f);
    CartesianVector direction(0.f, 0.f, 0.f);

    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;

    LArPcaHelper::RunPca(clusterHits, centroid, eigenValues, eigenVectors);

    m_clusterAttrMap[pCluster].m_centroid = centroid;
    m_clusterAttrMap[pCluster].m_eigenVectors = eigenVectors;

    CartesianVector innerCoordinate(0.f, 0.f, 0.f);
    CartesianVector outerCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);

    const float clusterLength{ LArClusterHelper::GetLength(pCluster) };

    m_clusterAttrMap[pCluster].m_length = clusterLength;
    m_clusterAttrMap[pCluster].m_innerCoordinate = innerCoordinate;
    m_clusterAttrMap[pCluster].m_outerCoordinate = outerCoordinate;
}

} // namespace lar_content
