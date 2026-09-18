/**
 *  @file   larpandoracontent/LArHelpers/LArClusterHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_HELPER_H
#define LAR_CLUSTER_HELPER_H 1

#include "Objects/Cluster.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  LArClusterHelper class
 */
class LArClusterHelper
{
public:
    typedef std::set<unsigned int> UIntSet;

    /**
     *  @brief  Get all hits (ordered and isoltated) associated with a two dimensional cluster
     *
     *  @param  pCluster the address of the cluster
     *  @param  caloHitList the output list of hits
     */
    static void GetAllHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get the hit type associated with a two dimensional cluster
     *
     *  @param  pCluster the address of the cluster
     *
     *  @return the cluster hit type
     */
    static pandora::HitType GetClusterHitType(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Divide an input cluster list into separate u, v and w lists (exception raised if alternative hit type encountered)
     *
     *  @param  inputClusters the input cluster list
     *  @param  clusterListU to receive the u clusters
     *  @param  clusterListV to receive the v clusters
     *  @param  clusterListW to receive the w clusters
     */
    static void GetClustersUVW(const pandora::ClusterList &inputClusters, pandora::ClusterList &clusterListU,
        pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW);

    /**
     *  @brief  Get the subset of clusters, from a provided list, that match the specified hit type
     *
     *  @param  inputClusters the input cluster list
     *  @param  hitType the specified hit type
     *  @param  clusterList to receive the clusters
     */
    static void GetClustersByHitType(const pandora::ClusterList &inputClusters, const pandora::HitType hitType, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get length squared of cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the length squared
     */
    static float GetLengthSquared(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get length of cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the length
     */
    static float GetLength(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get energy of cluster, based on length
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the energy
     */
    static float GetEnergyFromLength(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get number of layers spanned by cluster (1+Last-First)
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the layer span
     */
    static unsigned int GetLayerSpan(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return float
     */
    static float GetLayerOccupancy(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in a pair of clusters
     *
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     *
     *  @return float
     */
    static float GetLayerOccupancy(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between clusters in a pair of cluster lists
     *
     *  @param  clusterList1 the first cluster list
     *  @param  clusterList2 the second cluster list
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2);

    /**
     *  @brief  Get closest distance between a specified cluster and list of clusters
     *
     *  @param  pCluster address of the input cluster
     *  @param  clusterList list of input clusters
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster, const pandora::ClusterList &clusterList);

    /**
     *  @brief  Get closest distance between a pair of clusters
     *
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between a specified position and list of clusters
     *
     *  @param  position the position vector
     *  @param  clusterList list of input clusters
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::ClusterList &clusterList);

    /**
     *  @brief  Get closest distance between a specified position vector and the hits in a specified cluster
     *
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get closest distance between a specified position vector and the hits in a specified calo hit list
     *
     *  @param  position the position vector
     *  @param  caloHitList the list of calo hits
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get closest position in a list of clusters to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  clusterList list of input clusters
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::ClusterList &clusterList);

    /**
     *  @brief  Get closest position on a cluster to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get closest position of hits in a given calo hit list to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  caloHitList the list of calo hits
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get closest position of hits in a given calo hit list to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  caloHitList the list of calo hits
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::OrderedCaloHitList &caloHitList);

    /**
     *  @brief  Get pair of closest positions for a pair of clusters
     *
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  the closest position in the first cluster
     *  @param  the closest position in the second cluster
     */
    static void GetClosestPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianVector &position1, pandora::CartesianVector &position2);

    /**
     *  @brief  Get positions of the two most distant calo hits in a list of cluster (ordered by Z)
     *
     *  @param  clusterList the input cluster list
     *  @param  the inner extremal position
     *  @param  the outer extremal position
     */
    static void GetExtremalCoordinates(
        const pandora::ClusterList &clusterList, pandora::CartesianVector &innerCoordinate, pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief  Get positions of the two most distant calo hits in a cluster (ordered by Z)
     *
     *  @param  pCluster the input cluster
     *  @param  the inner extremal position
     *  @param  the outer extremal position
     */
    static void GetExtremalCoordinates(
        const pandora::Cluster *const pCluster, pandora::CartesianVector &innerCoordinate, pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief  Get positions of the two most distant calo hits in an ordered calo hit list (ordered by Z)
     *
     *  @param  orderedCaloHitList the ordered calo hit list
     *  @param  the inner extremal position
     *  @param  the outer extremal position
     */
    static void GetExtremalCoordinates(const pandora::OrderedCaloHitList &orderedCaloHitList, pandora::CartesianVector &innerCoordinate,
        pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief  Get positions of the two most distant points in a provided list (ordered by Z)
     *
     *  @param  coordinateVector the hit list
     *  @param  the inner extremal position
     *  @param  the outer extremal position
     */
    static void GetExtremalCoordinates(const pandora::CartesianPointVector &coordinateVector, pandora::CartesianVector &innerCoordinate,
        pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief  Get minimum and maximum X, Y and Z positions of the calo hits in a cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  the minimum positions (x,y,z)
     *  @param  the maximum positions (x,y,z)
     */
    static void GetClusterBoundingBox(
        const pandora::Cluster *const pCluster, pandora::CartesianVector &minimumCoordinate, pandora::CartesianVector &maximumCoordinate);

    /**
     *  @brief  Get vector of hit coordinates from an input cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  coordinateVector
     */
    static void GetCoordinateVector(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &coordinateVector);

    /**
     *  @brief  Get list of Calo hits from an input cluster that are contained in a bounding box.  The hits are sorted by position
     *
     *  @param  pCluster address of the cluster
     *  @param  lowerBound one opposing corner of the bounding box
     *  @param  upperBound the other opposing corner of the bounding box
     *  @param  caloHitList the CaloHitList to be filled
     */
    static void GetCaloHitListInBoundingBox(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBound,
        const pandora::CartesianVector &upperBound, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get the set of the daughter volumes that contains the cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  daughterVolumeIds output variable
     */
    static void GetDaughterVolumeIDs(const pandora::Cluster *const pCluster, UIntSet &daughterVolumeIds);

    /**
     *  @brief  Get average Z positions of the calo hits in a cluster in range xmin to xmax
     *
     *  @param  pCluster address of the cluster
     *  @param  xmin for range in x
     *  @param  xmax for range in x
     *  @param  averageZ to receive the average Z position
     *
     *  @return status code, faster than throwing in regular use-cases
     */
    static pandora::StatusCode GetAverageZ(const pandora::Cluster *const pCluster, const float xmin, const float xmax, float &averageZ);

    /**
     *  @brief  Determines if there is a hit between two hits.
     *
     *  @param  caloHits the collection of hits under consideration
     *  @param  pCaloHit1 the first calo hit
     *  @param  pCaloHit2 the second calo hit
     *
     *  @return true if there is a hit between the two hits, false otherwise
     */
    static bool HasBlockedPath(const pandora::CaloHitVector &caloHits, const pandora::CaloHit *const pCaloHit1, const pandora::CaloHit *const pCaloHit2);

    /**
     *  @brief  Constructs a list of hits from a cluster ordered along a sliding fit trajectory.
     *
     *  @param  pCluster the cluster whose hits are to be ordered
     *  @param  sfr the result of the sliding linear fit of the clusters
     *  @param  orderedHits the list in which to store the ordered hits
     */
    static void OrderHitsAlongTrajectory(const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &sfr, pandora::CaloHitList &orderedHits);

    /**
     *  @brief  Sort clusters by number of occupied layers, and by inner layer, then energy in event of a tie
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByNOccupiedLayers(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by number of hits, then layer span, then inner layer, then position, then pulse-height
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByNHits(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by layer span, then inner layer, then position, then pulse-height
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByLayerSpan(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by inner layer, then position, then pulse-height
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByInnerLayer(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by position, then pulse-height
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByPosition(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by pulse-height
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByPulseHeight(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort calo hits by their position (use Z, followed by X, followed by Y)
     *
     *  @param  pLhs address of first calo hit
     *  @param  pRhs address of second calo hit
     */
    static bool SortHitsByPosition(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs);

    /**
     *  @brief  Sort calo hits by their position (use X, followed by Z, followed by Y)
     *
     *  @param  pLhs address of first calo hit
     *  @param  pRhs address of second calo hit
     */
    static bool SortHitsByPositionInX(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs);

    /**
     *  @brief  Sort calo hits by their pulse height
     *
     *  @param  pLhs address of first calo hit
     *  @param  pRhs address of second calo hit
     */
    static bool SortHitsByPulseHeight(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs);

    /**
     *  @brief  Sort cartesian vectors by their position (use Z, followed by X, followed by Y)
     *
     *  @param  lhs first point
     *  @param  rhs second point
     */
    static bool SortCoordinatesByPosition(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_HELPER_H
