/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeViewTransverseTracksAlgorithm.h
 *
 *  @brief  Header file for the three view transverse tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H
#define LAR_THREE_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

namespace lar_content
{

class TransverseTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewTransverseTracksAlgorithm class
 */
class ThreeViewTransverseTracksAlgorithm : public NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<TransverseOverlapResult>>
{
public:
    typedef NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<TransverseOverlapResult>> BaseAlgorithm;

    /**
     *  @brief  Default constructor
     */
    ThreeViewTransverseTracksAlgorithm();

private:
    typedef std::map<unsigned int, TransverseOverlapResult> FitSegmentToOverlapResultMap;
    typedef std::map<unsigned int, FitSegmentToOverlapResultMap> FitSegmentMatrix;
    typedef std::map<unsigned int, FitSegmentMatrix> FitSegmentTensor;

    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
        const pandora::Cluster *const pClusterW, TransverseOverlapResult &overlapResult);

    /**
     *  @brief  Get the number of matched points for three fit segments and accompanying sliding fit results
     *
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     *  @param  fitSegmentTensor the fit segment tensor
     */
    void GetFitSegmentTensor(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
        const TwoDSlidingFitResult &slidingFitResultW, FitSegmentTensor &fitSegmentTensor) const;

    /**
     *  @brief  Get the overlap result for three fit segments and the accompanying sliding fit results
     *
     *  @param  fitSegmentU the fit segment u
     *  @param  fitSegmentV the fit segment v
     *  @param  fitSegmentW the fit segment w
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     *  @param  transverseOverlapResult to receive the transverse overlap result
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode GetSegmentOverlap(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW,
        const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
        const TwoDSlidingFitResult &slidingFitResultW, TransverseOverlapResult &transverseOverlapResult) const;

    /**
     *  @brief  Get the best overlap result, by examining the fit segment tensor
     *
     *  @param  fitSegmentTensor the fit segment tensor
     *  @param  bestTransverseOverlapResult to receive the best transverse overlap result
     */
    void GetBestOverlapResult(const FitSegmentTensor &fitSegmentTensor, TransverseOverlapResult &bestTransverseOverlapResult) const;

    /**
     *  @brief  Get the best of the accumulated overlap results for the possible connected previous segments
     *
     *  @param  indexU the index u
     *  @param  indexV the index v
     *  @param  indexW the index w
     *  @param  sumResults the accumulated overlap results for all segment index combinations, indexed (u * nIndicesV + v) * nIndicesW + w
     *  @param  nIndicesV the number of v indices
     *  @param  nIndicesW the number of w indices
     *
     *  @return the best previous result, or an uninitialized result if there are none
     */
    TransverseOverlapResult GetBestPreviousOverlapResult(const unsigned int indexU, const unsigned int indexV, const unsigned int indexW,
        const TransverseOverlapResultVector &sumResults, const unsigned int nIndicesV, const unsigned int nIndicesW) const;

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<TransverseTensorTool *> TensorToolVector;
    TensorToolVector m_algorithmToolVector; ///< The algorithm tool vector

    unsigned int m_nMaxTensorToolRepeats;   ///< The maximum number of repeat loops over tensor tools
    unsigned int m_maxFitSegmentIndex;      ///< The maximum number of fit segments used when identifying best overlap result
    float m_pseudoChi2Cut;                  ///< The pseudo chi2 cut to identify matched sampling points
    float m_minSegmentMatchedFraction;      ///< The minimum segment matched sampling fraction to allow segment grouping
    unsigned int m_minSegmentMatchedPoints; ///< The minimum number of matched segment sampling points to allow segment grouping
    float m_minOverallMatchedFraction;      ///< The minimum matched sampling fraction to allow particle creation
    unsigned int m_minOverallMatchedPoints; ///< The minimum number of matched segment sampling points to allow particle creation
    float m_minSamplingPointsPerLayer;      ///< The minimum number of sampling points per layer to allow particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseTensorTool class
 */
class TransverseTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewTransverseTracksAlgorithm::MatchingType::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H
