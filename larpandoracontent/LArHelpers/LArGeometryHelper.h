/**
 *  @file   larpandoracontent/LArHelpers/LArGeometryHelper.h
 *
 *  @brief  Header file for the geometry helper class.
 *
 *  $Log: $
 */
#ifndef LAR_GEOMETRY_HELPER_H
#define LAR_GEOMETRY_HELPER_H 1

#include "Objects/Cluster.h"
#include "Pandora/PandoraEnumeratedTypes.h"
#include "Pandora/StatusCodes.h"

#include <unordered_map>

namespace pandora
{
class CartesianVector;
class Pandora;
} // namespace pandora

namespace lar_content
{

class TwoDSlidingFitResult;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArGeometryHelper class
 */
class LArGeometryHelper
{
public:
    typedef std::set<unsigned int> UIntSet;

    /**
     *  @brief  Struct to hold detector boundaries
     */
    struct DetectorBoundaries
    {
        std::pair<float, float> m_xBoundaries = {std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()}; ///< {low, high} x boundaries of the detector
        std::pair<float, float> m_yBoundaries = {std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()}; ///< {low, high} y boundaries of the detector
        std::pair<float, float> m_zBoundaries = {std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()}; ///< {low, high} z boundaries of the detector
    };

    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     */
    static float MergeTwoPositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const float position1, const float position2);

    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  direction1 the direction in the first view
     *  @param  direction2 the direction in the second view
     */
    static pandora::CartesianVector MergeTwoDirections(const pandora::Pandora &pandora, const pandora::HitType view1,
        const pandora::HitType view2, const pandora::CartesianVector &direction1, const pandora::CartesianVector &direction2);

    /**
     *  @brief  Merge 2D positions from two views to give 2D position in third view
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 output position in the third view
     *  @param  chi-squared
     */
    static void MergeTwoPositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &position3, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from two views to give 2D position in third view
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeTwoPositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &outputU,
        pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 2D positions for each view
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  view3 the third view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 the position in the third view
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeThreePositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::HitType view3, const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        const pandora::CartesianVector &position3, pandora::CartesianVector &outputU, pandora::CartesianVector &outputV,
        pandora::CartesianVector &outputW, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 2D positions for each view
     *
     *  @param  pandora the associated pandora instance
     *  @param  positionU input position in the U view
     *  @param  positionV input position in the V view
     *  @param  positionW input position in the W view
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeThreePositions(const pandora::Pandora &pandora, const pandora::CartesianVector &positionU,
        const pandora::CartesianVector &positionV, const pandora::CartesianVector &positionW, pandora::CartesianVector &outputU,
        pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from two views to give unified 3D position
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3D output position in 3D
     *  @param  chi-squared
     */
    static void MergeTwoPositions3D(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &position3D, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 3D position
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  view3 the third view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 the position in the third view
     *  @param  position3D output position in 3D
     *  @param  chi-squared
     */
    static void MergeThreePositions3D(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::HitType view3, const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        const pandora::CartesianVector &position3, pandora::CartesianVector &position3D, float &chiSquared);

    /**
     *  @brief  Project 3D position into a given 2D view
     *
     *  @param  pandora the associated pandora instance
     *  @param  position3D the position in 3D
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector ProjectPosition(
        const pandora::Pandora &pandora, const pandora::CartesianVector &position3D, const pandora::HitType view);

    /**
     *  @brief  Project 3D direction into a given 2D view
     *
     *  @param  pandora the associated pandora instance
     *  @param  direction3D the direction in 3D
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector ProjectDirection(
        const pandora::Pandora &pandora, const pandora::CartesianVector &direction3D, const pandora::HitType view);

    /**
     *  @brief  Return the wire pitch
     *
     *  @param  pandora the associated pandora instance
     *  @param  maxWirePitchDiscrepancy maximum allowed discrepancy between lar tpc wire pitch value in the w view
     */
    static float GetWireZPitch(const pandora::Pandora &pandora, const float maxWirePitchDiscrepancy = 0.01);

    /**
     *  @brief  Return the wire pitch
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     *  @param  maxWirePitchDiscrepancy maximum allowed discrepancy between lar tpc wire pitch values
     */
    static float GetWirePitch(const pandora::Pandora &pandora, const pandora::HitType view, const float maxWirePitchDiscrepancy = 0.01);

    /**
     *  @brief  Return the ratio of the wire pitch of the specified view to the minimum wire pitch for the detector
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     *
     *  @return The ratio of the specified view's wire pitch to the minimum wire pitch
     */
    static float GetWirePitchRatio(const pandora::Pandora &pandora, const pandora::HitType view);

    /**
     *  @brief  Return the wire axis (vector perpendicular to the wire direction and drift direction)
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector GetWireAxis(const pandora::Pandora &pandora, const pandora::HitType view);

    /**
     *  @brief  Whether a 2D test point lies in a registered gap with the associated hit type
     *
     *  @param  pandora the associated pandora instance
     *  @param  testPoint the test point
     *  @param  hitType the hit type
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsInGap(const pandora::Pandora &pandora, const pandora::CartesianVector &testPoint2D, const pandora::HitType hitType,
        const float gapTolerance = 0.f);

    /**
     *  @brief  Whether a 3D test point lies in a registered gap with the associated hit type
     *
     *  @param  pandora the associated pandora instance
     *  @param  testPoint the test point
     *  @param  hitType the hit type
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsInGap3D(const pandora::Pandora &pandora, const pandora::CartesianVector &testPoint3D, const pandora::HitType hitType,
        const float gapTolerance = 0.f);

    /**
     *  @brief  Whether there is a gap in a cluster (described via its sliding fit result) at a specified x sampling position
     *
     *  @param  pandora the associated pandora instance
     *  @param  xSample the x sampling position
     *  @param  slidingFitResult the sliding fit result for a cluster
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsXSamplingPointInGap(
        const pandora::Pandora &pandora, const float xSample, const TwoDSlidingFitResult &slidingFitResult, const float gapTolerance = 0.f);

    /**
     *  @brief  Calculate the total distance within a given 2D region that is composed of detector gaps
     *
     *  @param  pandora the associated pandora instance
     *  @param  minZ the start position in Z
     *  @param  maxZ the end position in Z
     *  @param  hitType the hit type
     */
    static float CalculateGapDeltaZ(const pandora::Pandora &pandora, const float minZ, const float maxZ, const pandora::HitType hitType);

    /**
     *  @brief  Find the sigmaUVW value for the detector geometry
     *
     *  @param  pandora the associated pandora instance
     *  @param  maxSigmaDiscrepancy maximum allowed discrepancy between lar tpc sigmaUVW values
     */
    static float GetSigmaUVW(const pandora::Pandora &pandora, const float maxSigmaDiscrepancy = 0.01);

    /**
     *  @brief  Return the set of common daughter volumes between two 2D clusters
     *
     *  @param  intersect the set of shared daughter volumes
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     */
    static void GetCommonDaughterVolumes(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, UIntSet &intersect);

    /**
     *  @brief  Get the detector boundaries of the pandora instance
     *
     *  @param  pandora the pandora instance
     */
    static DetectorBoundaries GetDetectorBoundaries(const pandora::Pandora &pandora);

    /**
     *  @brief  Obtain the list of 'line gap' boundaries in the detector,
                'line gaps' are those that span the drift coordinate
     *
     *  @param  pandora the pandora instance
     *  @param[out] detXGaps the set of x positions of the detector gaps
     */    
    static void GetDetectorXGaps(const pandora::Pandora &pandora, std::set<float> &detXGaps);

    /**
     *  @brief  Return whether an input point is within the bounds of the detector
     *
     *  @param detectorBoundaries the detector boundaries
     *  @param  position the input CartesianVector
     *
     *  @return Whether the input position is within the detector
     */
    static bool IsInDetector(const DetectorBoundaries &detectorBoundaries, const pandora::CartesianVector &position);

    /**
     *  @brief  Calculate a chi-squared value (including hit width) for a triplet of hits
     *
     *  @param  pandora the associated pandora instance
     *  @param  pCaloHitU the calo hit in the U view
     *  @param  pCaloHitV the calo hit in the V view
     *  @param  pCaloHitW the calo hit in the W view
     */
    static float CalculateChiSquared(const pandora::Pandora &pandora, const pandora::CaloHit *const pCaloHitU,
        const pandora::CaloHit *const pCaloHitV, const pandora::CaloHit *const pCaloHitW);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArGeometryHelper::GetWireZPitch(const pandora::Pandora &pandora, const float maxWirePitchWDiscrepancy)
{
    return LArGeometryHelper::GetWirePitch(pandora, pandora::TPC_VIEW_W, maxWirePitchWDiscrepancy);
}

} // namespace lar_content

#endif // #ifndef LAR_GEOMETRY_HELPER_H
