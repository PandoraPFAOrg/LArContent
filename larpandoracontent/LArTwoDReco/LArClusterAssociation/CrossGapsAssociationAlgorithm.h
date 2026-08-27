/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h
 *
 *  @brief  Header file for the cross gaps association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSS_GAPS_ASSOCIATION_ALGORITHM_H
#define LAR_CROSS_GAPS_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CrossGapsAssociationAlgorithm class
 */
class CrossGapsAssociationAlgorithm : public ClusterAssociationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CrossGapsAssociationAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  innerFitResult two dimensional sliding fit result for the inner cluster
     *  @param  outerFitResult two dimensional sliding fit result for the outer cluster
     *
     *  @return boolean
     */
    bool AreClustersAssociated(const TwoDSlidingFitResult &innerFitResult, const TwoDSlidingFitResult &outerFitResult) const;

    /**
     *  @brief  Sample points along the extrapolation from a starting position to a target fit result to declare cluster association
     *
     *  @param  startPosition the start position
     *  @param  startDirection the start direction
     *  @param  targetFitResult the target fit result
     *
     *  @return boolean
     */
    bool IsAssociated(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection,
        const TwoDSlidingFitResult &targetFitResult, const bool isCrossingTPCCandidate) const;

    /**
     *  @brief  Whether a sampling point lies near a target 2d sliding fit result
     *
     *  @param  samplingPoint the sampling point
     *  @param  targetFitResult the target fit result
     *
     *  @return boolean
     */
    bool IsNearCluster(const pandora::CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult, const float numGapSteps, const bool isCrossingTPCCandidate) const;
    
    pandora::StatusCode FindLArTPCVolumeIds(const TwoDSlidingFitResult &fitResult, std::set<unsigned int> &volIds) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_minClusterHits;           ///< The minimum allowed number of hits in a clean cluster
    unsigned int m_minClusterLayers;         ///< The minimum allowed number of layers for a clean cluster
    unsigned int m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    unsigned int m_maxSamplingPoints;        ///< The maximum number of extension sampling points considered per association check
    float m_sampleStepSize;                  ///< The sampling step size used in association checks, units cm
    unsigned int m_maxUnmatchedSampleRun;    ///< The maximum run of unmatched (and non-gap) samples to consider before stopping
    float m_maxOnClusterDistance;            ///< The maximum distance between a sampling point and sliding fit to target cluster
    unsigned int m_minMatchedSamplingPoints; ///< Minimum number of matched sampling points to declare association
    float m_minMatchedSamplingFraction;      ///< Minimum ratio between matched sampling points and expectation to declare association
    float m_crossTPCStepModifier;            ///< Modifier to scale how much the angle between the startdirection of the sample points and the x-axis of the detector affects the m_maxUnmatchedSampleRun if crossing the APA
    float m_crossTPCOnClusterDistanceModifier; ///< Modifier to scale how much the number of samples that are in the APA gap affects the m_maxOnClusterDistance
    float m_boostStartStep;                ///< Minimum number of sample steps to cross the APA gap
    float m_crossTPCMaxAdditionalSteps;      ///< Maximum number of additional steps when sampling proccess crosses the APA gap to prevent circumstances where the value explodes due to angles close to parallel with the APA
    float m_gapTolerance;                    ///< The tolerance to use when querying whether a sampling point is in a gap, units cm
    bool m_visualize;                        ///< Whether to visualize sampling points/association diagnostics when enabled in settings
};

} // namespace lar_content

#endif // #ifndef LAR_CROSS_GAPS_ASSOCIATION_ALGORITHM_H
