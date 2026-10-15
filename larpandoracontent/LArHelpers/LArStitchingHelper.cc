/**
 *  @file   larpandoracontent/LArHelpers/LArStitchingHelper.cc
 *
 *  @brief  Implementation of the helper class for multiple tpcs
 *
 *  $Log: $
 */

#include "Managers/GeometryManager.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

const LArTPC &LArStitchingHelper::FindClosestTPC(const Pandora &pandora, const LArTPC &inputTPC, const bool checkPositive)
{
    if (!MultiPandoraApi::GetPandoraInstanceMap().count(&pandora))
    {
        std::cout << "LArStitchingHelper::FindClosestTPC - functionality only available to primary/master Pandora instance " << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const LArTPCMap &larTPCMap(pandora.GetGeometry()->GetLArTPCMap());

    const LArTPC *pClosestTPC(nullptr);
    float closestSeparation(std::numeric_limits<float>::max());
    const float maxDisplacement(30.f); // TODO: 30cm should be fine, but can we do better than a hard-coded number here?

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC &checkTPC(*(mapEntry.second));

        if (&inputTPC == &checkTPC)
            continue;

        if (checkPositive != (checkTPC.GetCenterX() > inputTPC.GetCenterX()))
            continue;

        const float deltaX(std::fabs(checkTPC.GetCenterX() - inputTPC.GetCenterX()));
        const float deltaY(std::fabs(checkTPC.GetCenterY() - inputTPC.GetCenterY()));
        const float deltaZ(std::fabs(checkTPC.GetCenterZ() - inputTPC.GetCenterZ()));

        if (deltaY > maxDisplacement || deltaZ > maxDisplacement)
            continue;

        if (deltaX < closestSeparation)
        {
            closestSeparation = deltaX;
            pClosestTPC = &checkTPC;
        }
    }

    if (!pClosestTPC)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return (*pClosestTPC);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::CanTPCsBeStitched(const LArTPC &firstTPC, const LArTPC &secondTPC)
{
    if (&firstTPC == &secondTPC)
        return false;

    // ATTN: We assume that Pfos are crossing either an anode-anode boundary or a cathode-cathode boundary
    if (firstTPC.IsDriftInPositiveX() == secondTPC.IsDriftInPositiveX())
        return false;

    return LArStitchingHelper::AreTPCsAdjacent(firstTPC, secondTPC);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::AreTPCsAdjacent(const LArTPC &firstTPC, const LArTPC &secondTPC)
{
    // Check the relative positions of the centres of each drift volume
    const float maxDisplacement(30.f); // TODO: 30cm should be fine, but can we do better than a hard-coded number here?
    const float widthX(0.5f * (firstTPC.GetWidthX() + secondTPC.GetWidthX()));
    const float deltaX(std::fabs(firstTPC.GetCenterX() - secondTPC.GetCenterX()));
    const float deltaY(std::fabs(firstTPC.GetCenterY() - secondTPC.GetCenterY()));
    const float deltaZ(std::fabs(firstTPC.GetCenterZ() - secondTPC.GetCenterZ()));

    if (std::fabs(deltaX - widthX) > maxDisplacement || deltaY > maxDisplacement || deltaZ > maxDisplacement)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::AreTPCsAdjacent(const Pandora &pandora, const LArTPC &firstTPC, const LArTPC &secondTPC)
{
    // Check if first tpc is just upstream of second tpc
    try
    {
        const LArTPC &firstTPCCheck(LArStitchingHelper::FindClosestTPC(pandora, secondTPC, true));
        const LArTPC &secondTPCCheck(LArStitchingHelper::FindClosestTPC(pandora, firstTPC, false));

        if ((&firstTPCCheck == &firstTPC) && (&secondTPCCheck == &secondTPC))
            return true;
    }
    catch (pandora::StatusCodeException &)
    {
    }

    // Check if second tpc is just upstream of first tpc
    try
    {
        const LArTPC &firstTPCCheck(LArStitchingHelper::FindClosestTPC(pandora, secondTPC, false));
        const LArTPC &secondTPCCheck(LArStitchingHelper::FindClosestTPC(pandora, firstTPC, true));

        if ((&firstTPCCheck == &firstTPC) && (&secondTPCCheck == &secondTPC))
            return true;
    }
    catch (pandora::StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetTPCBoundaryCenterX(const LArTPC &firstTPC, const LArTPC &secondTPC)
{
    if (!LArStitchingHelper::AreTPCsAdjacent(firstTPC, secondTPC))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (firstTPC.GetCenterX() < secondTPC.GetCenterX())
    {
        return 0.5 * ((firstTPC.GetCenterX() + 0.5 * firstTPC.GetWidthX()) + (secondTPC.GetCenterX() - 0.5 * secondTPC.GetWidthX()));
    }
    else
    {
        return 0.5 * ((firstTPC.GetCenterX() - 0.5 * firstTPC.GetWidthX()) + (secondTPC.GetCenterX() + 0.5 * secondTPC.GetWidthX()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetTPCBoundaryWidthX(const LArTPC &firstTPC, const LArTPC &secondTPC)
{
    if (!LArStitchingHelper::AreTPCsAdjacent(firstTPC, secondTPC))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (firstTPC.GetCenterX() < secondTPC.GetCenterX())
    {
        return ((secondTPC.GetCenterX() - 0.5 * secondTPC.GetWidthX()) - (firstTPC.GetCenterX() + 0.5 * firstTPC.GetWidthX()));
    }
    else
    {
        return ((firstTPC.GetCenterX() - 0.5 * firstTPC.GetWidthX()) - (secondTPC.GetCenterX() + 0.5 * secondTPC.GetWidthX()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetTPCDisplacement(const LArTPC &firstTPC, const LArTPC &secondTPC)
{
    const float deltaX(firstTPC.GetCenterX() - secondTPC.GetCenterX());
    const float deltaY(firstTPC.GetCenterY() - secondTPC.GetCenterY());
    const float deltaZ(firstTPC.GetCenterZ() - secondTPC.GetCenterZ());

    return std::sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArStitchingHelper::GetClosestVertices(const LArTPC &larTPC1, const LArTPC &larTPC2, const LArPointingCluster &pointingCluster1,
    const LArPointingCluster &pointingCluster2, LArPointingCluster::Vertex &closestVertex1, LArPointingCluster::Vertex &closestVertex2)
{
    if (&larTPC1 == &larTPC2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Find the closest vertices based on relative X positions in tpc
    const float dxVolume(larTPC2.GetCenterX() - larTPC1.GetCenterX());
    const float dx1(pointingCluster1.GetOuterVertex().GetPosition().GetX() - pointingCluster1.GetInnerVertex().GetPosition().GetX());
    const float dx2(pointingCluster2.GetOuterVertex().GetPosition().GetX() - pointingCluster2.GetInnerVertex().GetPosition().GetX());

    // ATTN: the near vertices below come from the sign of each cluster's x extent, which is unreliable for clusters nearly
    // parallel to the boundary; where that choice fails, fall back to the pair of vertices closest in the y-z plane
    if (std::fabs(dx1) < std::numeric_limits<float>::epsilon() || std::fabs(dx2) < std::numeric_limits<float>::epsilon())
    {
        GetClosestVerticesInYZ(pointingCluster1, pointingCluster2, closestVertex1, closestVertex2);
        return;
    }

    const bool useInner1((dxVolume > 0.f) == (dx1 < 0.f)); // xVol1 - OUTER - INNER - | - xVol2  [xVol2-xVol1>0; xOuter1-xInner1<0]
    const bool useInner2((dxVolume > 0.f) == (dx2 > 0.f)); // xVol1 - | - INNER - OUTER - xVol2  [xVol2-xVol1>0; xOuter2-xInner2>0]

    // Confirm that these really are the closest pair of vertices by checking the other possible pairs
    const LArPointingCluster::Vertex &nearVertex1(useInner1 ? pointingCluster1.GetInnerVertex() : pointingCluster1.GetOuterVertex());
    const LArPointingCluster::Vertex &nearVertex2(useInner2 ? pointingCluster2.GetInnerVertex() : pointingCluster2.GetOuterVertex());

    const LArPointingCluster::Vertex &farVertex1(useInner1 ? pointingCluster1.GetOuterVertex() : pointingCluster1.GetInnerVertex());
    const LArPointingCluster::Vertex &farVertex2(useInner2 ? pointingCluster2.GetOuterVertex() : pointingCluster2.GetInnerVertex());

    const float dxNearNear(0.f);
    const float dyNearNear(nearVertex1.GetPosition().GetY() - nearVertex2.GetPosition().GetY());
    const float dzNearNear(nearVertex1.GetPosition().GetZ() - nearVertex2.GetPosition().GetZ());
    const float drNearNearSquared(dxNearNear * dxNearNear + dyNearNear * dyNearNear + dzNearNear * dzNearNear);

    const float dxNearFar(std::fabs(dx2));
    const float dyNearFar(nearVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());
    const float dzNearFar(nearVertex1.GetPosition().GetZ() - farVertex2.GetPosition().GetZ());
    const float drNearFarSquared(dxNearFar * dxNearFar + dyNearFar * dyNearFar + dzNearFar * dzNearFar);

    const float dxFarNear(std::fabs(dx1));
    const float dyFarNear(farVertex1.GetPosition().GetY() - nearVertex2.GetPosition().GetY());
    const float dzFarNear(farVertex1.GetPosition().GetZ() - nearVertex2.GetPosition().GetZ());
    const float drFarNearSquared(dxFarNear * dxFarNear + dyFarNear * dyFarNear + dzFarNear * dzFarNear);

    const float dxFarFar(std::fabs(dx1) + std::fabs(dx2));
    const float dyFarFar(farVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());
    const float dzFarFar(farVertex1.GetPosition().GetZ() - farVertex2.GetPosition().GetZ());
    const float drFarFarSquared(dxFarFar * dxFarFar + dyFarFar * dyFarFar + dzFarFar * dzFarFar);

    // ATTN: another pair of vertices is closer than the near vertices chosen from the x extents, which happens when those x extents
    // are unreliable for clusters nearly parallel to the boundary; fall back to the pair of vertices closest in the y-z plane
    if (drNearNearSquared > std::min(drFarFarSquared, std::min(drNearFarSquared, drFarNearSquared)))
    {
        GetClosestVerticesInYZ(pointingCluster1, pointingCluster2, closestVertex1, closestVertex2);
        return;
    }

    closestVertex1 = nearVertex1;
    closestVertex2 = nearVertex2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArStitchingHelper::GetClosestVerticesInYZ(const LArPointingCluster &pointingCluster1, const LArPointingCluster &pointingCluster2,
    LArPointingCluster::Vertex &closestVertex1, LArPointingCluster::Vertex &closestVertex2)
{
    float bestDrSquared(std::numeric_limits<float>::max());

    for (const LArPointingCluster::Vertex *const pVertex1 : {&pointingCluster1.GetInnerVertex(), &pointingCluster1.GetOuterVertex()})
    {
        for (const LArPointingCluster::Vertex *const pVertex2 : {&pointingCluster2.GetInnerVertex(), &pointingCluster2.GetOuterVertex()})
        {
            const float dy(pVertex1->GetPosition().GetY() - pVertex2->GetPosition().GetY());
            const float dz(pVertex1->GetPosition().GetZ() - pVertex2->GetPosition().GetZ());
            const float drSquared(dy * dy + dz * dz);

            if (drSquared < bestDrSquared)
            {
                bestDrSquared = drSquared;
                closestVertex1 = *pVertex1;
                closestVertex2 = *pVertex2;
            }
        }
    }

    if (std::numeric_limits<float>::max() == bestDrSquared)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::CalculateX0(const LArTPC &firstTPC, const LArTPC &secondTPC, const LArPointingCluster::Vertex &firstVertex,
    const LArPointingCluster::Vertex &secondVertex)
{
    if (&firstTPC == &secondTPC)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // ATTN: Assume that Pfos are crossing either an anode-anode boundary or a cathode-cathode boundary
    if (firstTPC.IsDriftInPositiveX() == secondTPC.IsDriftInPositiveX())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Assume that Pfos have opposite direction components in x, and have some direction component in the y-z plane
    const CartesianVector firstDirectionX(firstVertex.GetDirection().GetX(), 0.f, 0.f);
    const CartesianVector secondDirectionX(secondVertex.GetDirection().GetX(), 0.f, 0.f);

    if (std::fabs(firstDirectionX.GetX()) > 1.f - std::numeric_limits<float>::epsilon() ||
        std::fabs(secondDirectionX.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (-firstDirectionX.GetDotProduct(secondDirectionX) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Calculate displacement in x from relative displacement in y-z plane
    const CartesianVector firstPositionYZ(0.f, firstVertex.GetPosition().GetY(), firstVertex.GetPosition().GetZ());
    const CartesianVector firstDirectionYZ(0.f, firstVertex.GetDirection().GetY(), firstVertex.GetDirection().GetZ());

    const CartesianVector secondPositionYZ(0.f, secondVertex.GetPosition().GetY(), secondVertex.GetPosition().GetZ());
    const CartesianVector secondDirectionYZ(0.f, secondVertex.GetDirection().GetY(), secondVertex.GetDirection().GetZ());

    const float firstDirectionYZmag(firstDirectionYZ.GetMagnitude());
    const float secondDirectionYZmag(secondDirectionYZ.GetMagnitude());

    if (firstDirectionYZmag < std::numeric_limits<float>::epsilon() || secondDirectionYZmag < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float R1(firstDirectionYZ.GetUnitVector().GetDotProduct(firstPositionYZ - secondPositionYZ) / firstDirectionYZmag);
    const float X1(-1.f *
        firstDirectionX.GetUnitVector().GetDotProduct(secondVertex.GetPosition() - (firstVertex.GetPosition() - firstVertex.GetDirection() * R1)));

    const float R2(secondDirectionYZ.GetUnitVector().GetDotProduct(secondPositionYZ - firstPositionYZ) / secondDirectionYZmag);
    const float X2(-1.f *
        secondDirectionX.GetUnitVector().GetDotProduct(firstVertex.GetPosition() - (secondVertex.GetPosition() - secondVertex.GetDirection() * R2)));

    // ATTN: By convention, X0 is half the displacement in x (because both Pfos will be corrected)
    return (X1 + X2) * 0.25f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::SortTPCs(const pandora::LArTPC *const pLhs, const pandora::LArTPC *const pRhs)
{
    if (std::fabs(pLhs->GetCenterX() - pRhs->GetCenterX()) > std::numeric_limits<float>::epsilon())
        return (pLhs->GetCenterX() < pRhs->GetCenterX());

    if (std::fabs(pLhs->GetCenterY() - pRhs->GetCenterY()) > std::numeric_limits<float>::epsilon())
        return (pLhs->GetCenterY() < pRhs->GetCenterY());

    return (pLhs->GetCenterZ() < pRhs->GetCenterZ());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::HasPfoBeenStitched(const ParticleFlowObject *const pPfo)
{
    const PropertiesMap &properties(pPfo->GetPropertiesMap());
    const PropertiesMap::const_iterator iter(properties.find("X0"));

    if (iter != properties.end())
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetPfoX0(const ParticleFlowObject *const pPfo)
{
    // ATTN: If no stitch is present return a shift of 0.f
    if (!LArStitchingHelper::HasPfoBeenStitched(pPfo))
        return 0.f;

    // ATTN: HasPfoBeenStitched ensures X0 is present in the properties map
    return pPfo->GetPropertiesMap().at("X0");
}

} // namespace lar_content
