/**
 *  @file   larpandoracontent/LArHelpers/LArPointingClusterHelper.cc
 *
 *  @brief  Implementation of the pointing cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

float LArPointingClusterHelper::GetLengthSquared(const LArPointingCluster &pointingCluster)
{
    const LArPointingCluster::Vertex &innerVertex(pointingCluster.GetInnerVertex());
    const LArPointingCluster::Vertex &outerVertex(pointingCluster.GetOuterVertex());
    return (innerVertex.GetPosition() - outerVertex.GetPosition()).GetMagnitudeSquared();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPointingClusterHelper::GetLength(const LArPointingCluster &pointingCluster)
{
    return std::sqrt(LArPointingClusterHelper::GetLengthSquared(pointingCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::IsNode(const CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterVertex,
    const float minLongitudinalDistance, const float maxTransverseDistance)
{
    float rL(0.f), rT(0.f);
    LArPointingClusterHelper::GetImpactParameters(daughterVertex.GetPosition(), daughterVertex.GetDirection(), parentVertex, rL, rT);

    if (std::fabs(rL) > std::fabs(minLongitudinalDistance) || rT > maxTransverseDistance)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::IsEmission(const CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterVertex,
    const float minLongitudinalDistance, const float maxLongitudinalDistance, const float maxTransverseDistance, const float angularAllowance)
{
    float rL(0.f), rT(0.f);
    LArPointingClusterHelper::GetImpactParameters(daughterVertex.GetPosition(), daughterVertex.GetDirection(), parentVertex, rL, rT);

    if (std::fabs(rL) > std::fabs(minLongitudinalDistance) && (rL < 0 || rL > maxLongitudinalDistance))
        return false;

    const float tanSqTheta(std::pow(std::tan(M_PI * angularAllowance / 180.f), 2.0));

    if (rT * rT > maxTransverseDistance * maxTransverseDistance + rL * rL * tanSqTheta)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArPointingClusterHelper::GetProjectedPosition(const CartesianVector &vertexPosition,
    const CartesianVector &vertexDirection, const pandora::Cluster *const pCluster, const float projectionAngularAllowance)
{
    const CaloHit *pClosestCaloHit(nullptr);
    float closestDistanceSquared(std::numeric_limits<float>::max());
    const float minCosTheta(std::cos(M_PI * projectionAngularAllowance / 180.f));

    for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit : *layerEntry.second)
        {
            const CartesianVector hitProjection(pCaloHit->GetPositionVector() - vertexPosition);
            const float distanceSquared(hitProjection.GetMagnitudeSquared());

            if (distanceSquared > std::numeric_limits<float>::epsilon())
            {
                // TODO Try to give more weight to on-axis projections
                if (distanceSquared < closestDistanceSquared)
                {
                    if (-hitProjection.GetUnitVector().GetDotProduct(vertexDirection) > minCosTheta)
                    {
                        pClosestCaloHit = pCaloHit;
                        closestDistanceSquared = distanceSquared;
                    }
                }
            }
            else
            {
                return pCaloHit->GetPositionVector();
            }
        }
    }

    if (pClosestCaloHit)
        return pClosestCaloHit->GetPositionVector();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetClosestVertices(const bool useX, const bool useY, const bool useZ, const LArPointingCluster &pointingClusterI,
    const LArPointingCluster &pointingClusterJ, LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    if (!LArPointingClusterHelper::TryGetClosestVertices(useX, useY, useZ, pointingClusterI, pointingClusterJ, closestVertexI, closestVertexJ))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::TryGetClosestVertices(const bool useX, const bool useY, const bool useZ, const LArPointingCluster &pointingClusterI,
    const LArPointingCluster &pointingClusterJ, LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    if (pointingClusterI.GetCluster() == pointingClusterJ.GetCluster())
        return false;

    if (!useX && !useY && !useZ)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int useInnerI = 0; useInnerI < 2; ++useInnerI)
    {
        const LArPointingCluster::Vertex &vtxI(useInnerI == 1 ? pointingClusterI.GetInnerVertex() : pointingClusterI.GetOuterVertex());
        const LArPointingCluster::Vertex &endI(useInnerI == 0 ? pointingClusterI.GetInnerVertex() : pointingClusterI.GetOuterVertex());

        for (unsigned int useInnerJ = 0; useInnerJ < 2; ++useInnerJ)
        {
            const LArPointingCluster::Vertex &vtxJ(useInnerJ == 1 ? pointingClusterJ.GetInnerVertex() : pointingClusterJ.GetOuterVertex());
            const LArPointingCluster::Vertex &endJ(useInnerJ == 0 ? pointingClusterJ.GetInnerVertex() : pointingClusterJ.GetOuterVertex());

            const float vtxI_vtxJ_dx(useX ? (vtxI.GetPosition().GetX() - vtxJ.GetPosition().GetX()) : 0.f);
            const float vtxI_vtxJ_dy(useY ? (vtxI.GetPosition().GetY() - vtxJ.GetPosition().GetY()) : 0.f);
            const float vtxI_vtxJ_dz(useZ ? (vtxI.GetPosition().GetZ() - vtxJ.GetPosition().GetZ()) : 0.f);
            const float vtxI_vtxJ(vtxI_vtxJ_dx * vtxI_vtxJ_dx + vtxI_vtxJ_dy * vtxI_vtxJ_dy + vtxI_vtxJ_dz * vtxI_vtxJ_dz);

            const float vtxI_endJ_dx(useX ? (vtxI.GetPosition().GetX() - endJ.GetPosition().GetX()) : 0.f);
            const float vtxI_endJ_dy(useY ? (vtxI.GetPosition().GetY() - endJ.GetPosition().GetY()) : 0.f);
            const float vtxI_endJ_dz(useZ ? (vtxI.GetPosition().GetZ() - endJ.GetPosition().GetZ()) : 0.f);
            const float vtxI_endJ(vtxI_endJ_dx * vtxI_endJ_dx + vtxI_endJ_dy * vtxI_endJ_dy + vtxI_endJ_dz * vtxI_endJ_dz);

            const float endI_vtxJ_dx(useX ? (endI.GetPosition().GetX() - vtxJ.GetPosition().GetX()) : 0.f);
            const float endI_vtxJ_dy(useY ? (endI.GetPosition().GetY() - vtxJ.GetPosition().GetY()) : 0.f);
            const float endI_vtxJ_dz(useZ ? (endI.GetPosition().GetZ() - vtxJ.GetPosition().GetZ()) : 0.f);
            const float endI_vtxJ(endI_vtxJ_dx * endI_vtxJ_dx + endI_vtxJ_dy * endI_vtxJ_dy + endI_vtxJ_dz * endI_vtxJ_dz);

            const float endI_endJ_dx(useX ? (endI.GetPosition().GetX() - endJ.GetPosition().GetX()) : 0.f);
            const float endI_endJ_dy(useY ? (endI.GetPosition().GetY() - endJ.GetPosition().GetY()) : 0.f);
            const float endI_endJ_dz(useZ ? (endI.GetPosition().GetZ() - endJ.GetPosition().GetZ()) : 0.f);
            const float endI_endJ(endI_endJ_dx * endI_endJ_dx + endI_endJ_dy * endI_endJ_dy + endI_endJ_dz * endI_endJ_dz);

            if ((vtxI_vtxJ < std::min(vtxI_endJ, std::min(endI_vtxJ, endI_endJ))) && (endI_endJ > std::max(vtxI_endJ, std::max(endI_vtxJ, vtxI_vtxJ))))
            {
                closestVertexI = vtxI;
                closestVertexJ = vtxJ;
                return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::TryGetClosestVertices(const LArPointingCluster &pointingClusterI, const LArPointingCluster &pointingClusterJ,
    LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    return LArPointingClusterHelper::TryGetClosestVertices(true, true, true, pointingClusterI, pointingClusterJ, closestVertexI, closestVertexJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::TryGetClosestVerticesInX(const LArPointingCluster &pointingClusterI, const LArPointingCluster &pointingClusterJ,
    LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    return LArPointingClusterHelper::TryGetClosestVertices(true, false, false, pointingClusterI, pointingClusterJ, closestVertexI, closestVertexJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetClosestVertices(const LArPointingCluster &pointingClusterI, const LArPointingCluster &pointingClusterJ,
    LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    return LArPointingClusterHelper::GetClosestVertices(true, true, true, pointingClusterI, pointingClusterJ, closestVertexI, closestVertexJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetClosestVerticesInX(const LArPointingCluster &pointingClusterI, const LArPointingCluster &pointingClusterJ,
    LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    return LArPointingClusterHelper::GetClosestVertices(true, false, false, pointingClusterI, pointingClusterJ, closestVertexI, closestVertexJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetClosestVerticesInYZ(const LArPointingCluster &pointingClusterI,
    const LArPointingCluster &pointingClusterJ, LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ)
{
    return LArPointingClusterHelper::GetClosestVertices(false, true, true, pointingClusterI, pointingClusterJ, closestVertexI, closestVertexJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetImpactParametersInYZ(
    const LArPointingCluster::Vertex &initialVertex, const LArPointingCluster::Vertex &targetVertex, float &longitudinal, float &transverse)
{
    if (std::fabs(initialVertex.GetDirection().GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    const pandora::CartesianVector initialPosition(0.f, initialVertex.GetPosition().GetY(), initialVertex.GetPosition().GetZ());
    const pandora::CartesianVector initialDirection(0.f, initialVertex.GetDirection().GetY(), initialVertex.GetDirection().GetZ());
    const pandora::CartesianVector targetPosition(0.f, targetVertex.GetPosition().GetY(), targetVertex.GetPosition().GetZ());

    return LArPointingClusterHelper::GetImpactParameters(initialPosition, initialDirection.GetUnitVector(), targetPosition, longitudinal, transverse);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetImpactParameters(
    const LArPointingCluster::Vertex &pointingVertex, const LArPointingCluster::Vertex &targetVertex, float &longitudinal, float &transverse)
{
    return LArPointingClusterHelper::GetImpactParameters(
        pointingVertex.GetPosition(), pointingVertex.GetDirection(), targetVertex.GetPosition(), longitudinal, transverse);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetImpactParameters(
    const LArPointingCluster::Vertex &pointingVertex, const CartesianVector &targetPosition, float &longitudinal, float &transverse)
{
    return LArPointingClusterHelper::GetImpactParameters(
        pointingVertex.GetPosition(), pointingVertex.GetDirection(), targetPosition, longitudinal, transverse);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetImpactParameters(const CartesianVector &initialPosition, const CartesianVector &initialDirection,
    const CartesianVector &targetPosition, float &longitudinal, float &transverse)
{
    // sign convention for longitudinal distance:
    // -positive value means initial position is downstream of target position
    transverse = initialDirection.GetCrossProduct(targetPosition - initialPosition).GetMagnitude();
    longitudinal = -initialDirection.GetDotProduct(targetPosition - initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetAverageDirection(
    const LArPointingCluster::Vertex &firstVertex, const LArPointingCluster::Vertex &secondVertex, CartesianVector &averageDirection)
{
    const Cluster *const pFirstCluster(firstVertex.GetCluster());
    const Cluster *const pSecondCluster(secondVertex.GetCluster());

    if (pFirstCluster == pSecondCluster)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float energy1(pFirstCluster->GetHadronicEnergy());
    const float energy2(pSecondCluster->GetHadronicEnergy());

    if ((energy1 < std::numeric_limits<float>::epsilon()) || (energy2 < std::numeric_limits<float>::epsilon()))
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    averageDirection = (firstVertex.GetDirection() * energy1 + secondVertex.GetDirection() * energy2).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetIntersection(const LArPointingCluster::Vertex &firstVertex,
    const LArPointingCluster::Vertex &secondVertex, CartesianVector &intersectPosition, float &firstDisplacement, float &secondDisplacement)
{
    const Cluster *const pFirstCluster(firstVertex.GetCluster());
    const Cluster *const pSecondCluster(secondVertex.GetCluster());

    if (pFirstCluster == pSecondCluster)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    LArPointingClusterHelper::GetIntersection(firstVertex.GetPosition(), firstVertex.GetDirection(), secondVertex.GetPosition(),
        secondVertex.GetDirection(), intersectPosition, firstDisplacement, secondDisplacement);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetIntersection(const CartesianVector &a1, const CartesianVector &a2, const CartesianVector &b1,
    const CartesianVector &b2, CartesianVector &intersectPosition, float &firstDisplacement, float &secondDisplacement)
{
    // note: input lines are r = a1 + P * a2 and r = b1 + Q * b2
    const float cosTheta = a2.GetDotProduct(b2);

    // lines must be non-parallel
    if (1.f - std::fabs(cosTheta) < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // calculate the intersection (by minimising the distance between the lines)
    const float P = ((a2 - b2 * cosTheta).GetDotProduct(b1 - a1)) / (1.f - cosTheta * cosTheta);
    const float Q = ((a2 * cosTheta - b2).GetDotProduct(b1 - a1)) / (1.f - cosTheta * cosTheta);

    // position of intersection (or point of closest approach)
    intersectPosition = (a1 + a2 * P + b1 + b2 * Q) * 0.5f;

    // displacements of intersection from input vertices
    firstDisplacement = P;
    secondDisplacement = Q;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::GetIntersection(const LArPointingCluster::Vertex &vertexCluster, const Cluster *const pTargetCluster,
    CartesianVector &intersectPosition, float &displacementL, float &displacementT)
{
    displacementT = +std::numeric_limits<float>::max();
    displacementL = -std::numeric_limits<float>::max();

    float rL(0.f), rT(0.f);
    float figureOfMerit(std::numeric_limits<float>::max());
    float foundIntersection(false);

    const OrderedCaloHitList &orderedCaloHitList(pTargetCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList.begin(), iterEnd1 = orderedCaloHitList.end(); iter1 != iterEnd1; ++iter1)
    {
        for (CaloHitList::const_iterator iter2 = iter1->second->begin(), iterEnd2 = iter1->second->end(); iter2 != iterEnd2; ++iter2)
        {
            const CartesianVector &hitPosition = (*iter2)->GetPositionVector();

            LArPointingClusterHelper::GetImpactParameters(vertexCluster.GetPosition(), vertexCluster.GetDirection(), hitPosition, rL, rT);

            if (rT < figureOfMerit)
            {
                figureOfMerit = rT;

                displacementL = rL;
                displacementT = rT;
                intersectPosition = hitPosition;
                foundIntersection = true;
            }
        }
    }

    if (!foundIntersection)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex LArPointingClusterHelper::GetBestVertexEstimate(const LArPointingClusterVertexList &vertexList,
    const LArPointingClusterList &pointingClusterList, const float minLongitudinalDistance, const float maxLongitudinalDistance,
    const float maxTransverseDistance, const float angularAllowance)
{
    float bestAssociatedEnergy(0.f);
    LArPointingClusterVertexList::const_iterator bestVertexIter(vertexList.end());

    for (LArPointingClusterVertexList::const_iterator iter = vertexList.begin(), iterEnd = vertexList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &vertex(*iter);

        LArPointingClusterVertexList associatedVertices;
        LArPointingClusterHelper::CollectAssociatedClusters(vertex, pointingClusterList, minLongitudinalDistance, maxLongitudinalDistance,
            maxTransverseDistance, angularAllowance, associatedVertices);

        const float associatedEnergy(LArPointingClusterHelper::GetAssociatedEnergy(vertex, associatedVertices));

        if (associatedEnergy > bestAssociatedEnergy)
        {
            bestVertexIter = iter;
            bestAssociatedEnergy = associatedEnergy;
        }
    }

    if (vertexList.end() == bestVertexIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return (*bestVertexIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingClusterHelper::CollectAssociatedClusters(const LArPointingCluster::Vertex &vertex, const LArPointingClusterList &inputList,
    const float minLongitudinalDistance, const float maxLongitudinalDistance, const float maxTransverseDistance,
    const float angularAllowance, LArPointingClusterVertexList &outputList)
{
    for (LArPointingClusterList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster &pointingCluster = *iter;
        const LArPointingCluster::Vertex &innerVertex = pointingCluster.GetInnerVertex();
        const LArPointingCluster::Vertex &outerVertex = pointingCluster.GetOuterVertex();

        const float innerDistanceSquared = (innerVertex.GetPosition() - vertex.GetPosition()).GetMagnitudeSquared();
        const float outerDistanceSquared = (outerVertex.GetPosition() - vertex.GetPosition()).GetMagnitudeSquared();

        const LArPointingCluster::Vertex &chosenVertex((innerDistanceSquared < outerDistanceSquared) ? innerVertex : outerVertex);

        if (LArPointingClusterHelper::IsNode(vertex.GetPosition(), chosenVertex, minLongitudinalDistance, maxTransverseDistance) ||
            LArPointingClusterHelper::IsEmission(vertex.GetPosition(), chosenVertex, minLongitudinalDistance, maxLongitudinalDistance,
                maxTransverseDistance, angularAllowance))
        {
            outputList.push_back(chosenVertex);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPointingClusterHelper::GetAssociatedEnergy(const LArPointingCluster::Vertex &vertex, const LArPointingClusterVertexList &associatedVertices)
{
    float associatedEnergy(0.f);

    for (LArPointingClusterVertexList::const_iterator iter = associatedVertices.begin(), iterEnd = associatedVertices.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &clusterVertex(*iter);
        const Cluster *const pCluster(clusterVertex.GetCluster());

        const float clusterEnergy(LArClusterHelper::GetEnergyFromLength(pCluster));
        const float clusterLength(LArClusterHelper::GetLength(pCluster));
        const float deltaLength(clusterVertex.GetDirection().GetDotProduct(vertex.GetPosition() - clusterVertex.GetPosition()));

        if (deltaLength < std::numeric_limits<float>::epsilon())
        {
            associatedEnergy += clusterEnergy;
        }
        else if (deltaLength < clusterLength)
        {
            associatedEnergy += clusterEnergy * (1.f - (deltaLength / clusterLength));
        }
    }

    return associatedEnergy;
}

} // namespace lar_content
