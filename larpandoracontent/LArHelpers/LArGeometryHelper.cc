/**
 *  @file   larpandoracontent/LArHelpers/LArGeometryHelper.cc
 *
 *  @brief  Implementation of the geometry helper class.
 *
 *  $Log: $
 */

#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"

#include "Geometry/DetectorGap.h"
#include "Geometry/LArTPC.h"

#include "Objects/CartesianVector.h"

#include "Pandora/Pandora.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "Plugins/LArTransformationPlugin.h"

using namespace pandora;

namespace lar_content
{

float LArGeometryHelper::MergeTwoPositions(const Pandora &pandora, const HitType view1, const HitType view2, const float position1, const float position2)
{
    if (view1 == view2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        return pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoW(position1, position2);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        return pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoW(position2, position1);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        return pandora.GetPlugins()->GetLArTransformationPlugin()->WUtoV(position1, position2);
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        return pandora.GetPlugins()->GetLArTransformationPlugin()->WUtoV(position2, position1);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        return pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoU(position1, position2);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        return pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoU(position2, position1);
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::MergeTwoDirections(
    const Pandora &pandora, const HitType view1, const HitType view2, const CartesianVector &direction1, const CartesianVector &direction2)
{
    // x components must be consistent
    if (direction1.GetX() * direction2.GetX() < 0.f)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // scale x components to a common value
    const float s1((std::fabs(direction1.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction2.GetX()) : 1.f);
    const float s2((std::fabs(direction2.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction1.GetX()) : 1.f);

    float pX(s1 * direction1.GetX()), pU(0.f), pV(0.f), pW(0.f);
    HitType newView(HIT_CUSTOM);

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        pU = s1 * direction1.GetZ();
        pV = s2 * direction2.GetZ();
        newView = TPC_VIEW_W;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        pV = s1 * direction1.GetZ();
        pU = s2 * direction2.GetZ();
        newView = TPC_VIEW_W;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        pW = s1 * direction1.GetZ();
        pU = s2 * direction2.GetZ();
        newView = TPC_VIEW_V;
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        pU = s1 * direction1.GetZ();
        pW = s2 * direction2.GetZ();
        newView = TPC_VIEW_V;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        pV = s1 * direction1.GetZ();
        pW = s2 * direction2.GetZ();
        newView = TPC_VIEW_U;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        pW = s1 * direction1.GetZ();
        pV = s2 * direction2.GetZ();
        newView = TPC_VIEW_U;
    }

    if (newView == TPC_VIEW_W)
        return CartesianVector(pX, 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoW(pU, pV)).GetUnitVector();

    if (newView == TPC_VIEW_U)
        return CartesianVector(pX, 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoU(pV, pW)).GetUnitVector();

    if (newView == TPC_VIEW_V)
        return CartesianVector(pX, 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->WUtoV(pW, pU)).GetUnitVector();

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions(const Pandora &pandora, const HitType view1, const HitType view2,
    const CartesianVector &position1, const CartesianVector &position2, CartesianVector &position3, float &chiSquared)
{
    if (view1 == view2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float X3((position1.GetX() + position2.GetX()) / 2.f);
    const float Z1(position1.GetZ());
    const float Z2(position2.GetZ());
    const float Z3(LArGeometryHelper::MergeTwoPositions(pandora, view1, view2, Z1, Z2));

    position3.SetValues(X3, 0.f, Z3);
    const float sigmaUVW(LArGeometryHelper::GetSigmaUVW(pandora));
    chiSquared = ((X3 - position1.GetX()) * (X3 - position1.GetX()) + (X3 - position2.GetX()) * (X3 - position2.GetX())) / (sigmaUVW * sigmaUVW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions(const Pandora &pandora, const HitType view1, const HitType view2, const CartesianVector &position1,
    const CartesianVector &position2, CartesianVector &outputU, CartesianVector &outputV, CartesianVector &outputW, float &chiSquared)
{
    if (view1 == view2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector position3(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(pandora, view1, view2, position1, position2, position3, chiSquared);

    float aveU(0.f), aveV(0.f), aveW(0.f);
    const float Z1(position1.GetZ()), Z2(position2.GetZ()), Z3(position3.GetZ()), aveX(position3.GetX());

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        aveU = Z1;
        aveV = Z2;
        aveW = Z3;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        aveV = Z1;
        aveW = Z2;
        aveU = Z3;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        aveW = Z1;
        aveU = Z2;
        aveV = Z3;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        aveV = Z1;
        aveU = Z2;
        aveW = Z3;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        aveW = Z1;
        aveV = Z2;
        aveU = Z3;
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        aveU = Z1;
        aveW = Z2;
        aveV = Z3;
    }

    outputU.SetValues(aveX, 0.f, aveU);
    outputV.SetValues(aveX, 0.f, aveV);
    outputW.SetValues(aveX, 0.f, aveW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions(const Pandora &pandora, const HitType view1, const HitType view2, const HitType view3,
    const CartesianVector &position1, const CartesianVector &position2, const CartesianVector &position3, CartesianVector &outputU,
    CartesianVector &outputV, CartesianVector &outputW, float &chiSquared)
{
    if ((view1 == view2) || (view2 == view3) || (view3 == view1))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position1, position2, position3, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position3, position1, position2, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position2, position3, position1, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position2, position1, position3, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position3, position2, position1, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position1, position3, position2, outputU, outputV, outputW, chiSquared);
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions(const Pandora &pandora, const CartesianVector &positionU, const CartesianVector &positionV,
    const CartesianVector &positionW, CartesianVector &outputU, CartesianVector &outputV, CartesianVector &outputW, float &chiSquared)
{
    const float YfromUV(pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoY(positionU.GetZ(), positionV.GetZ()));
    const float YfromUW(pandora.GetPlugins()->GetLArTransformationPlugin()->UWtoY(positionU.GetZ(), positionW.GetZ()));
    const float YfromVW(pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoY(positionV.GetZ(), positionW.GetZ()));

    const float ZfromUV(pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoZ(positionU.GetZ(), positionV.GetZ()));
    const float ZfromUW(pandora.GetPlugins()->GetLArTransformationPlugin()->UWtoZ(positionU.GetZ(), positionW.GetZ()));
    const float ZfromVW(pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoZ(positionV.GetZ(), positionW.GetZ()));

    // ATTN For detectors where w and z are equivalent, remain consistent with original treatment. TODO Use new treatment always.
    const bool useOldWZEquivalentTreatment(std::fabs(ZfromUW - ZfromVW) < std::numeric_limits<float>::epsilon());
    const float aveX((positionU.GetX() + positionV.GetX() + positionW.GetX()) / 3.f);
    const float aveY(useOldWZEquivalentTreatment ? YfromUV : (YfromUV + YfromUW + YfromVW) / 3.f);
    const float aveZ(useOldWZEquivalentTreatment ? (positionW.GetZ() + 2.f * ZfromUV) / 3.f : (ZfromUV + ZfromUW + ZfromVW) / 3.f);

    const float aveU(pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(aveY, aveZ));
    const float aveV(pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(aveY, aveZ));
    const float aveW(pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(aveY, aveZ));

    outputU.SetValues(aveX, 0.f, aveU);
    outputV.SetValues(aveX, 0.f, aveV);
    outputW.SetValues(aveX, 0.f, aveW);

    const float sigmaUVW(LArGeometryHelper::GetSigmaUVW(pandora));
    chiSquared = ((outputU.GetX() - positionU.GetX()) * (outputU.GetX() - positionU.GetX()) +
                     (outputV.GetX() - positionV.GetX()) * (outputV.GetX() - positionV.GetX()) +
                     (outputW.GetX() - positionW.GetX()) * (outputW.GetX() - positionW.GetX()) +
                     (outputU.GetZ() - positionU.GetZ()) * (outputU.GetZ() - positionU.GetZ()) +
                     (outputV.GetZ() - positionV.GetZ()) * (outputV.GetZ() - positionV.GetZ()) +
                     (outputW.GetZ() - positionW.GetZ()) * (outputW.GetZ() - positionW.GetZ())) /
        (sigmaUVW * sigmaUVW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions3D(const Pandora &pandora, const HitType view1, const HitType view2,
    const CartesianVector &position1, const CartesianVector &position2, CartesianVector &position3D, float &chiSquared)
{
    CartesianVector positionU(0.f, 0.f, 0.f), positionV(0.f, 0.f, 0.f), positionW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(pandora, view1, view2, position1, position2, positionU, positionV, positionW, chiSquared);

    position3D.SetValues(positionW.GetX(), pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoY(positionU.GetZ(), positionV.GetZ()),
        pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoZ(positionU.GetZ(), positionV.GetZ()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions3D(const Pandora &pandora, const HitType view1, const HitType view2, const HitType view3,
    const CartesianVector &position1, const CartesianVector &position2, const CartesianVector &position3, CartesianVector &position3D, float &chiSquared)
{
    CartesianVector positionU(0.f, 0.f, 0.f), positionV(0.f, 0.f, 0.f), positionW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeThreePositions(pandora, view1, view2, view3, position1, position2, position3, positionU, positionV, positionW, chiSquared);

    position3D.SetValues(positionW.GetX(), pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoY(positionU.GetZ(), positionV.GetZ()),
        pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoZ(positionU.GetZ(), positionV.GetZ()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::ProjectPosition(const Pandora &pandora, const CartesianVector &position3D, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        return CartesianVector(
            position3D.GetX(), 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(position3D.GetY(), position3D.GetZ()));
    }

    else if (view == TPC_VIEW_V)
    {
        return CartesianVector(
            position3D.GetX(), 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(position3D.GetY(), position3D.GetZ()));
    }

    else if (view == TPC_VIEW_W)
    {
        return CartesianVector(
            position3D.GetX(), 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(position3D.GetY(), position3D.GetZ()));
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::ProjectDirection(const Pandora &pandora, const CartesianVector &direction3D, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        return CartesianVector(
            direction3D.GetX(), 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(direction3D.GetY(), direction3D.GetZ()))
            .GetUnitVector();
    }

    else if (view == TPC_VIEW_V)
    {
        return CartesianVector(
            direction3D.GetX(), 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(direction3D.GetY(), direction3D.GetZ()))
            .GetUnitVector();
    }

    else if (view == TPC_VIEW_W)
    {
        return CartesianVector(
            direction3D.GetX(), 0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(direction3D.GetY(), direction3D.GetZ()))
            .GetUnitVector();
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::GetWirePitch(const Pandora &pandora, const HitType view, const float maxWirePitchDiscrepancy)
{
    if (view != TPC_VIEW_U && view != TPC_VIEW_V && view != TPC_VIEW_W && view != TPC_3D)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (view == TPC_3D)
    {
        const float pitchU{LArGeometryHelper::GetWirePitch(pandora, TPC_VIEW_U, maxWirePitchDiscrepancy)};
        const float pitchV{LArGeometryHelper::GetWirePitch(pandora, TPC_VIEW_V, maxWirePitchDiscrepancy)};
        const float pitchW{LArGeometryHelper::GetWirePitch(pandora, TPC_VIEW_W, maxWirePitchDiscrepancy)};
        return std::max({pitchU, pitchV, pitchW});
    }

    const LArTPCMap &larTPCMap(pandora.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "LArGeometryHelper::GetWirePitch - LArTPC description not registered with Pandora as required " << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    const float wirePitch(view == TPC_VIEW_U ? pFirstLArTPC->GetWirePitchU()
                                             : (view == TPC_VIEW_V ? pFirstLArTPC->GetWirePitchV() : pFirstLArTPC->GetWirePitchW()));

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        const float alternateWirePitch(
            view == TPC_VIEW_U ? pLArTPC->GetWirePitchU() : (view == TPC_VIEW_V ? pLArTPC->GetWirePitchV() : pLArTPC->GetWirePitchW()));

        if (std::fabs(wirePitch - alternateWirePitch) > maxWirePitchDiscrepancy)
        {
            std::cout << "LArGeometryHelper::GetWirePitch - LArTPC configuration not supported" << std::endl;
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }

    return wirePitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::GetWirePitchRatio(const Pandora &pandora, const HitType view)
{
    const float pitchU{LArGeometryHelper::GetWirePitch(pandora, TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(pandora, TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(pandora, TPC_VIEW_W)};
    const float pitchMin{std::min({pitchU, pitchV, pitchW})};
    const float pitchView{LArGeometryHelper::GetWirePitch(pandora, view)};

    return pitchView / pitchMin;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::GetWireAxis(const Pandora &pandora, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        return CartesianVector(0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(1.f, 0.f),
            pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(0.f, 1.f));
    }

    else if (view == TPC_VIEW_V)
    {
        return CartesianVector(0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(1.f, 0.f),
            pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(0.f, 1.f));
    }

    else if (view == TPC_VIEW_W)
    {
        return CartesianVector(0.f, pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(1.f, 0.f),
            pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(0.f, 1.f));
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::GetCommonDaughterVolumes(const Cluster *const pCluster1, const Cluster *const pCluster2, UIntSet &intersect)
{
    UIntSet daughterVolumeIds1, daughterVolumeIds2;

    LArClusterHelper::GetDaughterVolumeIDs(pCluster1, daughterVolumeIds1);
    LArClusterHelper::GetDaughterVolumeIDs(pCluster2, daughterVolumeIds2);

    std::set_intersection(daughterVolumeIds1.begin(), daughterVolumeIds1.end(), daughterVolumeIds2.begin(), daughterVolumeIds2.end(),
        std::inserter(intersect, intersect.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsInGap(const Pandora &pandora, const CartesianVector &testPoint2D, const HitType hitType, const float gapTolerance)
{
    // ATTN: input test point MUST be a 2D position vector
    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        if (pDetectorGap->IsInGap(testPoint2D, hitType, gapTolerance))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsInGap3D(const Pandora &pandora, const CartesianVector &testPoint3D, const HitType hitType, const float gapTolerance)
{
    const CartesianVector testPoint2D(LArGeometryHelper::ProjectPosition(pandora, testPoint3D, hitType));
    return LArGeometryHelper::IsInGap(pandora, testPoint2D, hitType, gapTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsXSamplingPointInGap(const Pandora &pandora, const float xSample, const TwoDSlidingFitResult &slidingFitResult, const float gapTolerance)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(slidingFitResult.GetCluster()));
    const CartesianVector minLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());

    const bool minLayerIsAtLowX(minLayerPosition.GetX() < maxLayerPosition.GetX());
    const CartesianVector &lowXCoordinate(minLayerIsAtLowX ? minLayerPosition : maxLayerPosition);
    const CartesianVector &highXCoordinate(minLayerIsAtLowX ? maxLayerPosition : minLayerPosition);

    if ((xSample > lowXCoordinate.GetX()) && (xSample < highXCoordinate.GetX()))
    {
        CartesianVector slidingFitPosition(0.f, 0.f, 0.f);

        if (STATUS_CODE_SUCCESS == slidingFitResult.GetGlobalFitPositionAtX(xSample, slidingFitPosition))
            return (LArGeometryHelper::IsInGap(pandora, slidingFitPosition, hitType, gapTolerance));
    }

    const CartesianVector lowXDirection(
        minLayerIsAtLowX ? slidingFitResult.GetGlobalMinLayerDirection() : slidingFitResult.GetGlobalMaxLayerDirection());
    const CartesianVector highXDirection(
        minLayerIsAtLowX ? slidingFitResult.GetGlobalMaxLayerDirection() : slidingFitResult.GetGlobalMinLayerDirection());

    const bool sampleIsNearerToLowX(std::fabs(xSample - lowXCoordinate.GetX()) < std::fabs(xSample - highXCoordinate.GetX()));
    const CartesianVector &startPosition(sampleIsNearerToLowX ? lowXCoordinate : highXCoordinate);
    const CartesianVector &startDirection(sampleIsNearerToLowX ? lowXDirection : highXDirection);

    if (std::fabs(startDirection.GetX()) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float pathLength((xSample - startPosition.GetX()) / startDirection.GetX());
    const CartesianVector samplingPoint(startPosition + startDirection * pathLength);

    return (LArGeometryHelper::IsInGap(pandora, samplingPoint, hitType, gapTolerance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateGapDeltaZ(const Pandora &pandora, const float minZ, const float maxZ, const HitType hitType)
{
    if (maxZ - minZ < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float gapDeltaZ(0.f);

    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap = dynamic_cast<const LineGap *>(pDetectorGap);

        if (!pLineGap)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const LineGapType lineGapType(pLineGap->GetLineGapType());

        if (!(((TPC_VIEW_U == hitType) && (TPC_WIRE_GAP_VIEW_U == lineGapType)) || ((TPC_VIEW_V == hitType) && (TPC_WIRE_GAP_VIEW_V == lineGapType)) ||
                ((TPC_VIEW_W == hitType) && (TPC_WIRE_GAP_VIEW_W == lineGapType))))
        {
            continue;
        }

        if ((pLineGap->GetLineStartZ() > maxZ) || (pLineGap->GetLineEndZ() < minZ))
            continue;

        const float gapMinZ(std::max(minZ, pLineGap->GetLineStartZ()));
        const float gapMaxZ(std::min(maxZ, pLineGap->GetLineEndZ()));

        if ((gapMaxZ - gapMinZ) > std::numeric_limits<float>::epsilon())
            gapDeltaZ += (gapMaxZ - gapMinZ);
    }

    return gapDeltaZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::GetSigmaUVW(const Pandora &pandora, const float maxSigmaDiscrepancy)
{
    const LArTPCMap &larTPCMap(pandora.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "LArGeometryHelper::GetSigmaUVW - LArTPC description not registered with Pandora as required " << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    const float sigmaUVW(pFirstLArTPC->GetSigmaUVW());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);

        if (std::fabs(sigmaUVW - pLArTPC->GetSigmaUVW()) > maxSigmaDiscrepancy)
        {
            std::cout << "LArGeometryHelper::GetSigmaUVW - Plugin does not support provided LArTPC configurations " << std::endl;
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }

    return sigmaUVW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArGeometryHelper::DetectorBoundaries LArGeometryHelper::GetDetectorBoundaries(const Pandora &pandora)
{
    DetectorBoundaries detectorBoundaries;

    const LArTPCMap &larTPCMap(pandora.GetGeometry()->GetLArTPCMap());

    for (const auto &entry : larTPCMap)
    {
        detectorBoundaries.m_xBoundaries.first =
            std::min(detectorBoundaries.m_xBoundaries.first, entry.second->GetCenterX() - (entry.second->GetWidthX() * 0.5f));
        detectorBoundaries.m_xBoundaries.second =
            std::max(detectorBoundaries.m_xBoundaries.second, entry.second->GetCenterX() + (entry.second->GetWidthX() * 0.5f));
        detectorBoundaries.m_yBoundaries.first =
            std::min(detectorBoundaries.m_yBoundaries.first, entry.second->GetCenterY() - (entry.second->GetWidthY() * 0.5f));
        detectorBoundaries.m_yBoundaries.second =
            std::max(detectorBoundaries.m_yBoundaries.second, entry.second->GetCenterY() + (entry.second->GetWidthY() * 0.5f));
        detectorBoundaries.m_zBoundaries.first =
            std::min(detectorBoundaries.m_zBoundaries.first, entry.second->GetCenterZ() - (entry.second->GetWidthZ() * 0.5f));
        detectorBoundaries.m_zBoundaries.second =
            std::max(detectorBoundaries.m_zBoundaries.second, entry.second->GetCenterZ() + (entry.second->GetWidthZ() * 0.5f));
    }

    return detectorBoundaries;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::GetDetectorXGaps(const Pandora &pandora, std::set<float> &detXGaps)
{
    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pDetectorGap));

        if (!pLineGap)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        
        if (pLineGap->GetLineGapType() == TPC_DRIFT_GAP)
        {
            detXGaps.insert(pLineGap->GetLineStartX());
            detXGaps.insert(pLineGap->GetLineEndX());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsInDetector(const DetectorBoundaries &detectorBoundaries, const CartesianVector &position)
{
    if ((position.GetX() < detectorBoundaries.m_xBoundaries.first) || (position.GetX() > detectorBoundaries.m_xBoundaries.second))
    {
        return false;
    }

    if ((position.GetY() < detectorBoundaries.m_yBoundaries.first) || (position.GetY() > detectorBoundaries.m_yBoundaries.second))
    {
        return false;
    }

    if ((position.GetZ() < detectorBoundaries.m_zBoundaries.first) || (position.GetZ() > detectorBoundaries.m_zBoundaries.second))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateChiSquared(
    const Pandora &pandora, const CaloHit *const pCaloHitU, const CaloHit *const pCaloHitV, const CaloHit *const pCaloHitW)
{
    const CartesianVector posU(pCaloHitU->GetPositionVector());
    const CartesianVector posV(pCaloHitV->GetPositionVector());
    const CartesianVector posW(pCaloHitW->GetPositionVector());

    const float y_uv(pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoY(posU.GetZ(), posV.GetZ()));
    const float y_uw(pandora.GetPlugins()->GetLArTransformationPlugin()->UWtoY(posU.GetZ(), posW.GetZ()));
    const float y_vw(pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoY(posV.GetZ(), posW.GetZ()));

    const float z_uv(pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoZ(posU.GetZ(), posV.GetZ()));
    const float z_uw(pandora.GetPlugins()->GetLArTransformationPlugin()->UWtoZ(posU.GetZ(), posW.GetZ()));
    const float z_vw(pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoZ(posV.GetZ(), posW.GetZ()));

    const float x((posU.GetX() + posV.GetX() + posW.GetX()) / 3.f);
    const float y((y_uv + y_uw + y_vw) / 3.f);
    const float z((z_uv + z_uw + z_vw) / 3.f);

    const float u(pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(y, z));
    const float v(pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(y, z));
    const float w(pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(y, z));

    const float sigma2_uvw{static_cast<float>(std::pow(LArGeometryHelper::GetSigmaUVW(pandora), 2))};
    const float sigma2_x_u{static_cast<float>(std::pow(0.5f * pCaloHitU->GetCellSize1(), 2))};
    const float sigma2_x_v{static_cast<float>(std::pow(0.5f * pCaloHitV->GetCellSize1(), 2))};
    const float sigma2_x_w{static_cast<float>(std::pow(0.5f * pCaloHitW->GetCellSize1(), 2))};

    float chiSquared = std::pow(x - posU.GetX(), 2) / sigma2_x_u + std::pow(x - posV.GetX(), 2) / sigma2_x_v + std::pow(x - posW.GetX(), 2) / sigma2_x_w;
    chiSquared += (std::pow(u - posU.GetZ(), 2) + std::pow(v - posV.GetZ(), 2) + std::pow(w - posW.GetZ(), 2)) / sigma2_uvw;

    return chiSquared;
}

} // namespace lar_content
