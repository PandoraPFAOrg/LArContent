/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDAssociationAlgorithm.h
 *
 *  @brief  Header file for the 3D association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREED_ASSOCIATION_ALGORITHM_H
#define LAR_THREED_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

namespace lar_content
{

/**
 *  @brief  ThreeDAssociationAlgorithm class
 */
class ThreeDAssociationAlgorithm : public ClusterAssociationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDAssociationAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  pInnerCluster address of the inner cluster
     *  @param  pOuterCluster address of the outer cluster
     *
     *  @return whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::Cluster *const pInnerCluster, const pandora::Cluster *const pOuterCluster) const;

    /**
     *  @brief Populate cluster attributes 
     *
     *  @param  pCluster address of the cluster
     */
    void PopulateClusterAttributes(const pandora::Cluster *const pCluster) const;

    float m_minClusterLength;            ///< minimum cluster length
    float m_shortClusterLength;          ///< length threshold for short cluster
    float m_minCosRelativeAngle;         ///< maximum allowed relative angle between associated clusters
    float m_maxGapDistanceSquared;       ///< maximum allowed distance (squared) between associated clusters
    float m_channelPitchTolerance;       ///< factor applied to channel pitch 

    // struct to cache cluster attributes
    struct ClusterAttr {
        float m_length;
        pandora::CartesianVector m_innerCoordinate = pandora::CartesianVector(0.f, 0.f, 0.f);
        pandora::CartesianVector m_outerCoordinate = pandora::CartesianVector(0.f, 0.f, 0.f);
        pandora::CartesianVector m_centroid = pandora::CartesianVector(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors m_eigenVectors;
    };

    typedef std::unordered_map<const pandora::Cluster*, ClusterAttr> ClusterAttrMap;
    mutable ClusterAttrMap m_clusterAttrMap;
};

} // namespace lar_content

#endif // #ifndef LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H
