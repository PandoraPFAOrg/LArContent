/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.h
 *
 *  @brief  Header file for the track recovery algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
#define LAR_TRACK_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include <optional>
#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TrackRecoveryAlgorithm class
 */
class TrackRecoveryAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackRecoveryAlgorithm();

private:
    typedef std::unordered_map<pandora::HitType, const pandora::Cluster *> ViewToClusterMap;
    typedef std::unordered_map<const pandora::Pfo *, ViewToClusterMap> PfoToViewClusterMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::CaloHitList> ClusterToHitMap;
    typedef std::unordered_map<const pandora::Cluster *, const TwoDSlidingFitResult> ClusterToFitMap;
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::Cluster *> HitToClusterToMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::Pfo *> ClusterToPfoMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Performs a sliding linear fit to get an ordered set of hits along a cluster trajectory.
     *
     *  @param  pCluster the cluster for which hits should be ordered
     *  @param  orderedHits the list in which to store the ordered hits
     */
    std::optional<TwoDSlidingFitResult> FitAndOrderCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &orderedHits) const;

    /**
     *  @brief  Identifies hits that are missing from the PFO based on independent triplet matching. hitsA acts as the reference list to look
     *          up the hit triplets. Hits are deemed missing if two of the three hits in a triplet are present in the PFO.
     *
     *  @param  hitsA the list of hits in the reference view
     *  @param  hitsB the list of hits in the second view
     *  @param  hitsC the list of hits in the third view
     *  @param  unmatchedHitsB the set in which to store hits missing in only the second view
     *  @param  unmatchedHitsC the set in which to store hits missing in only the third view
     */
    void FindUnmatchedHits(const pandora::CaloHitList &hitsA, const pandora::CaloHitList &hitsB, const pandora::CaloHitList &hitsC,
        pandora::CaloHitSet &unmatchedHitsB, pandora::CaloHitSet &unmatchedHitsC) const;

    /**
     *  @brief  Identifies hits to merge from a collection of hits provisionally identified as being missing from a cluster.
     *          If a candidate hit lies along the existing set of hits, we look for a gap in the existing set that could be filled by this hit.
     *          If we find one, we consider this a valid addition. Alternatively, if a candidate hit lies outside the existing set of hits, we
     *          assume this is an extension of the track and consider this a valid addition. If there are no existing hits in the cluster, we
     *          assume the view was missing, and add all hits. The key assumption here is that the unmatched hits are already considered good
     *          candidates based on the triplet matching.
     *
     *  @param  pCluster the cluster for which to identify hits to merge
     *  @param  clusterHits the list of hits already associated to the cluster
     *  @param  unmatchedHits the set of hits not in the cluster, under consideration for merging
     *  @param  clusterToFitMap the map containing sliding fit results for clusters
     *  @param  mergeHits the list in which to store hits to merge
     */
    void IdentifyHitsToMerge(const pandora::Cluster *pCluster, const pandora::CaloHitList &clusterHits,
        const pandora::CaloHitSet &unmatchedHits, const ClusterToFitMap &clusterToFitMap, pandora::CaloHitList &mergeHits) const;

    /**
     *  @brief  Filters a list of hits to merge based on their perpendicular distance to the sliding fit for the cluster.
     *          This is a final quality cut to ensure we only merge hits that are consistent with the existing cluster.
     *
     *  @param  sfr the sliding fit result for the cluster
     *  @param[in,out]  mergeHits the list of hits provisionally identified for merging, which is filtered in place
     */
    void FilterHitsToMerge(const TwoDSlidingFitResult &sfr, pandora::CaloHitList &mergeHits) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputPfoListName;                ///< The name of the input PFO list containing track-like PFOs
    pandora::StringVector m_inputClusterListNames; ///< The list of cluster list names to consider when looking for clusters to recover
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
