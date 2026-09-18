/**
 *  @file   larpandoracontent/LArReclustering/ShortTrackReclusteringAlgorithm.h
 *
 *  @brief  Tries to identify short tracks that are either lost, or under-clustered and recover missing hits.
 *          This algorithm looks at all three views of a PFO to try to identify if one or more views exhibit step changes in ADC values that
 *          might be indicative of a decay or inelastic interaction. If such a change occurs it looks for consistency across views and also
 *          examines nearby unclustered hits, or hits in nearby clusters to see if a more coherent clustering can be identified.
 *          If these various conditions are met, then reclustering is performed.
 *
 *  $Log: $
 */

#ifndef LAR_SHORT_TRACK_RECLUSTERING_ALGORITHM_H
#define LAR_SHORT_TRACK_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArReclustering/ThreeDReclusteringFigureOfMeritBaseTool.h"

namespace lar_content
{

/**
  *  @brief  ShortTrackReclusteringAlgorithm class
  */
class ShortTrackReclusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShortTrackReclusteringAlgorithm();

    /**
    *  @brief  Default destructor
    */
    ~ShortTrackReclusteringAlgorithm() = default;

private:
    typedef std::tuple<const pandora::CaloHit *, const pandora::CaloHit *, const pandora::CaloHit *> HitTriplet;
    typedef std::unordered_map<pandora::HitType, pandora::CaloHitSet> ViewToHitsMap;
    typedef std::unordered_map<pandora::HitType, pandora::ClusterList> ViewToClustersMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::Pfo *> ClusterToPfoMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::FloatVector> ClusterToAdcMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::CaloHitSet> ClusterToHitsMap;
    typedef std::unordered_map<const pandora::Pfo *, std::vector<HitTriplet>> PfoToHitTripletsMap;

    struct Partition
    {
        Partition(const pandora::Pfo *const pCurrentPfo, const HitTriplet &hitTriplet, const pandora::CaloHitList &hitsU,
            const pandora::CaloHitList &hitsV, const pandora::CaloHitList &hitsW) :
            m_pCurrentPfo(pCurrentPfo),
            m_hitTriplet(hitTriplet),
            m_hitsU(hitsU),
            m_hitsV(hitsV),
            m_hitsW(hitsW)
        {
        }

        const pandora::Pfo *m_pCurrentPfo;
        HitTriplet m_hitTriplet;
        pandora::CaloHitList m_hitsU, m_hitsV, m_hitsW;
    };
    typedef std::vector<Partition> PartitionVector;

    struct ProtoPfo
    {
        PandoraContentApi::ParticleFlowObject::Parameters m_pfoParameters;
        const pandora::Pfo *m_pOldPfo;
        std::unordered_map<pandora::HitType, const pandora::Cluster *> m_viewToClusterMap;
        std::unordered_map<pandora::HitType, pandora::CaloHitList> m_viewToHitsMap;
    };
    typedef std::vector<ProtoPfo> ProtoPfoVector;

    pandora::StatusCode Run();

    /**
     *  @brief  Helper function to get a list from Pandora content, and check that it is valid
     *
     *  @param  listName the name of the list to retrieve
     *  @param  pList the pointer to the list to be retrieved
     *
     *  @return true if the list is successfully retrieved and valid, false otherwise
     */
    template <typename T>
    bool GetList(const std::string &listName, const T *&pList) const;

    /**
     *  @brief  Collects the clusters currently assigned to PFOs, and maps them by view, and also maps clusters to their parent PFO
     *
     *  @param  pfoList the list of pfos to consider
     *  @param  viewToClustersMap the map in which to store the clusters, mapped by view
     *  @param  clusterToPfoMap the map in which to store the mapping of clusters to their parent PFO
     */
    void CollectClusters(const pandora::PfoList &pfoList, ViewToClustersMap &viewToClustersMap, ClusterToPfoMap &clusterToPfoMap) const;

    /**
     *  @brief  Loops over clusters and performs a sliding linear fit to get an ordered set of hits along the cluster trajectory, and the
     *          corresponding fit result, which are stored in maps for later use.
     *
     *  @param  viewToClustersMap the map of clusters mapped by view
     */
    void FitAndOrderClusters(const ViewToClustersMap &viewToClustersMap);

    /**
     *  @brief  Loops over clusters and looks for evidence of discontinuous changes in ADC values, and collects the corresponding hits
     *
     *  @param  clusterToPfoMap the map of clusters to their parent PFO
     *  @param  clusterToHitsMap the map in which to store the mapping of clusters to the hits associated with any identified ADC discontinuities
     */
    void FindAdcDiscontinuities(const ClusterToPfoMap &clusterToPfoMap, ClusterToHitsMap &clusterToHitsMap) const;

    /**
     *  @brief  Loops over the hits associated with identified ADC discontinuities and looks for corresponding hits in other views, to identify
     *  potential triplets of hits across views that are consistent with a common 3D position
     *
     *  @param  clusterToHitsMap the map of clusters to the hits associated with any identified ADC discontinuities
     *  @param  clusterToPfoMap the map of clusters to their parent PFO
     *  @param  pfoToHitTripletsMap the map in which to store the mapping of PFOs to triplets of hits across views that are consistent with
     *          a common 3D position
     */
    void MatchAdcDiscontinuities(
        const ClusterToHitsMap &clusterToHitsMap, const ClusterToPfoMap &clusterToPfoMap, PfoToHitTripletsMap &pfoToHitTripletsMap) const;

    /**
     *  @brief  Uses the identified discontinuity triplets to find coherent changes across all three views and proposes new partitions.
     *
     *  @param  pfoToHitTripletsMap the map of PFOs to triplets of hits across views that are consistent with a common 3D position
     *  @param  partitions the proposed partitions if coherent alternatives can be found
     */
    void PartitionDiscontinuities(const PfoToHitTripletsMap &pfoToHitTripletsMap, PartitionVector &partitions) const;

    /**
     *  @brief  Different partitions can have overlap. This function throws filters out the smaller partitions
     *
     *  @param[in,out]  partitions the proposed partitions, which are filtered in place
     */
    void FilterPartitions(PartitionVector &partitions) const;

    /**
     *  @brief  Applies the reclustering described by a given partition vector
     *
     *  @param  partitions the proposed reclustering partitions
     */
    void Recluster(const PartitionVector &partitions) const;

    /**
     *  @brief  Utility function that performs the partitioning of a given set of hits according to a provided hit about which the split
     *          should be performed. The split hit is then allocated to the partition that is most consistent with its ADC
     *
     *  @param  caloHits the set of hits to partition
     *  @param  pSplitHit the hit about which to partition the set of hits
     *  @param  cluster1Hits the vector in which to store the first set of partitioned hits
     *  @param  cluster2Hits the vector in which to store the second set of partitioned hits
     */
    void PartitionHits(const pandora::CaloHitList &caloHits, const pandora::CaloHit *const pSplitHit, pandora::CaloHitList &cluster1Hits,
        pandora::CaloHitList &cluster2Hits) const;

    /**
     *  @brief  Gets the hits associated with a cluster, ordered relative to a specified vertex
     *
     *  @param  clusterHits the hits associated with the cluster for which to retrieve the ordered hits
     *  @param  vertex the vertex relative to which to order the hits
     *  @param  orderedHits the vector in which to store the ordered hits
     */
    void OrderHitsRelativeToVertex(
        const pandora::CaloHitVector &clusterHits, const LArPointingCluster::Vertex &vertex, pandora::CaloHitVector &orderedHits) const;

    /**
     *  @brief  Gets the moving average of the ADC values for ordered sets of hits, window is backward looking, so considers the "current"i
     *          hit and the previous (window-1) hits.
     *
     *  @param  adcs an ordered set of ADC values for which to calculate the moving average
     *  @param  movingAdc the vector in which to store the moving average of the ADC values for the hits
     *  @param  movingVariance the vector in which to store the moving variance of the ADC values for the hits
     *  @param  window the size of the window over which to calculate the moving average
     */
    void GetAdcMovingAverage(const pandora::FloatVector &adcs, pandora::FloatVector &movingAdc, pandora::FloatVector &movingVariance,
        const size_t window = 3) const;

    /**
     *  @brief  Gets the median value for a vector. The input vector need not be sorted.
     *
     *  @param  values a set of values for which to calculate the median
     *  @return the median value for the input set of values
     */
    template <typename T>
    double GetMedian(const std::vector<T> &values) const;

    /**
     *  @brief  Gets the indices of hits in an ordered set of hits for which there is a step change in ADC values. This function looks for step changes in
     *          ADC, but with relatively stable ADC either side of the discontinuity, indicative of a decay or inelastic interaction. This method attempts
     *          to filter out Bragg peaks or regions of high volatility.
     *
     *  @param  hits an ordered set of hits for which to identify stable ADC discontinuities
     *  @param  discontinuities the vector in which to store the indices of hits for which there is a stable ADC discontinuity
     *  @param  window the size of the window over which to calculate the moving average and identify step changes
     */
    void GetStableAdcDiscontinuities(const pandora::CaloHitVector &hits, pandora::IntVector &discontinuities, const size_t window = 3) const;

    /**
     *  @brief  Gets the ADC values for a set of hits and normalizes relative to the median ADC
     *
     *  @param  hits a set of hits for which to get the normalized ADC values
     *  @param  normalizedAdc the vector in which to store the normalized ADC values for the hits
     */
    void NormalizeAdc(const pandora::CaloHitVector &hits, pandora::FloatVector &normalizedAdc) const;

    /**
     *  @brief  Identifies whether a set of hits exhibits a Bragg peak. This method considers the consensus view from computation of the linear slope,
     *          quadratic curvature, contrast ratio and monotonicity.
     *
     *  @param  hits the set of hits for which to identify whether there is a Bragg peak
     *  @param  start the index of the first hit in the set of hits to consider when looking for a Bragg peak
     *  @param  end the index of the last hit in the set of hits to consider when looking for a Bragg peak
     *
     *  @return true if there is evidence for a Bragg peak in the specified range of hits, false otherwise
     */
    bool IsBraggPeak(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the linear slope of the ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the linear slope score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the linear slope score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the linear slope score
     *
     *  @return a score indicative of the linear slope of the ADC values for the specified range of hits
     */
    float GetLinearSlopeScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the quadratic curvature of the ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the quadratic curvature score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the quadratic curvature score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the quadratic curvature score
     *
     *  @return a score indicative of the quadratic curvature of the ADC values for the specified range of hits
     */
    float GetQuadraticCurvatureScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the contrast in ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the contrast score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the contrast score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the contrast score
     *
     *  @return a score indicative of the contrast in ADC values for the specified range of hits
     */
    float GetContrastScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the monotonicity of the ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the monotonicity score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the monotonicity score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the monotonicity score
     *
     *  @return a score indicative of the monotonicity of the ADC values for the specified range of hits
     */
    float GetMonotonicityScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets the ratio of the median ADC value either side of the pivot.
     *
     *  @param  hits the set of hits for which to calculate the balance
     *  @param  pivot the index of the hit about which to calculate the balance
     *
     *  @return the ratio of the median ADC value either side of the pivot
     */
    float GetBalance(const pandora::CaloHitList &hits, const size_t pivot) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;            ///< Name of list of calo hits to consider during reclustering
    std::string m_pfoListName;                ///< Name of list of track-like pfos to consider for reclustering
    pandora::StringVector m_clusterListNames; ///< Names of lists of clusters to consider during reclustering
    float m_maxHitDiscrepancy;                ///< Maximum allowed discrepancy between 2D hit positions and triplet projections
    float m_maxHitDiscrepancySquared;         ///< Square of m_maxHitDiscrepancy, cached for efficiency
    float m_balanceThresholdLow;              ///< Threshold for the balance ratio below which a discontinuity is considered to be present
    float m_balanceThresholdHigh;             ///< Threshold for the balance ratio above which a discontinuity is considered to be present
    float m_braggLinearSlopeThreshold;        ///< Threshold for the linear slope score above which a Bragg peak is considered to be present
    float m_braggCurvatureThreshold; ///< Threshold for the quadratic curvature score above which a Bragg peak is considered to be present
    float m_braggContrastThreshold;  ///< Threshold for the contrast score above which a Bragg peak is considered to be present
    float m_braggMonotonicityThreshold; ///< Threshold for the monotonicity score above which a Bragg peak is considered to be present
    std::unordered_map<const pandora::Cluster *, TwoDSlidingFitResult> m_clusterToSFRMap;
    std::unordered_map<const pandora::Cluster *, pandora::CaloHitList> m_clusterToOrderedHitsMap;
};

} // namespace lar_content

#endif // #ifndef LAR_SHORT_TRACK_RECLUSTERING_ALGORITHM_H
