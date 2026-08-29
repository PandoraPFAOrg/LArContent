/**
 *  @file   larpandoracontent/LArHelpers/LArClusterHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

void LArClusterHelper::GetAllHits(const Cluster *const pCluster, CaloHitList &caloHitList)
{
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    const CaloHitList &isolatedHitList{pCluster->GetIsolatedCaloHitList()};
    caloHitList.insert(caloHitList.end(), isolatedHitList.begin(), isolatedHitList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitType LArClusterHelper::GetClusterHitType(const Cluster *const pCluster)
{
    if (0 == pCluster->GetNCaloHits())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // TODO Find a way of confirming single hit-type clustering mode; currently only checked in ListPreparation algorithm
    // if (!pandora->GetSettings()->SingleHitTypeClusteringMode())
    //     throw StatusCodeException(STATUS_CODE_FAILURE);

    return (*(pCluster->GetOrderedCaloHitList().begin()->second->begin()))->GetHitType();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClustersUVW(const ClusterList &inputClusters, ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW)
{
    for (ClusterList::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(*iter));

        if (TPC_VIEW_U == hitType)
            clusterListU.push_back(*iter);
        else if (TPC_VIEW_V == hitType)
            clusterListV.push_back(*iter);
        else if (TPC_VIEW_W == hitType)
            clusterListW.push_back(*iter);
        else
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClustersByHitType(const ClusterList &inputClusters, const HitType hitType, ClusterList &clusterList)
{
    for (ClusterList::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        if (hitType == LArClusterHelper::GetClusterHitType(*iter))
            clusterList.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLengthSquared(const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // ATTN In 2D case, we will actually calculate the quadrature sum of deltaX and deltaU/V/W
    float minX(std::numeric_limits<float>::max()), maxX(-std::numeric_limits<float>::max());
    float minY(std::numeric_limits<float>::max()), maxY(-std::numeric_limits<float>::max());
    float minZ(std::numeric_limits<float>::max()), maxZ(-std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CartesianVector &hitPosition((*hitIter)->GetPositionVector());
            minX = std::min(hitPosition.GetX(), minX);
            maxX = std::max(hitPosition.GetX(), maxX);
            minY = std::min(hitPosition.GetY(), minY);
            maxY = std::max(hitPosition.GetY(), maxY);
            minZ = std::min(hitPosition.GetZ(), minZ);
            maxZ = std::max(hitPosition.GetZ(), maxZ);
        }
    }

    const float deltaX(maxX - minX), deltaY(maxY - minY), deltaZ(maxZ - minZ);
    return (deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLength(const Cluster *const pCluster)
{
    return std::sqrt(LArClusterHelper::GetLengthSquared(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetEnergyFromLength(const Cluster *const pCluster)
{
    const float dEdX(0.002f); // approximately 2 MeV/cm
    return (dEdX * LArClusterHelper::GetLength(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArClusterHelper::GetLayerSpan(const Cluster *const pCluster)
{
    return (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster)
{
    const unsigned int nOccupiedLayers(pCluster->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const unsigned int nOccupiedLayers(pCluster1->GetOrderedCaloHitList().size() + pCluster2->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + std::max(pCluster1->GetOuterPseudoLayer(), pCluster2->GetOuterPseudoLayer()) -
        std::min(pCluster1->GetInnerPseudoLayer(), pCluster2->GetInnerPseudoLayer()));

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const ClusterList &clusterList1, const ClusterList &clusterList2)
{
    if (clusterList1.empty() || clusterList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float closestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iterEnd1 = clusterList1.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster1, clusterList2));

        if (thisDistance < closestDistance)
            closestDistance = thisDistance;
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster, const ClusterList &clusterList)
{
    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float closestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pTestCluster = *iter;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster, pTestCluster));

        if (thisDistance < closestDistance)
            closestDistance = thisDistance;
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    LArClusterHelper::GetClosestPositions(pCluster1, pCluster2, closestPosition1, closestPosition2);

    return (closestPosition1 - closestPosition2).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const ClusterList &clusterList)
{
    return (position - LArClusterHelper::GetClosestPosition(position, clusterList)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const Cluster *const pCluster)
{
    return (position - LArClusterHelper::GetClosestPosition(position, pCluster)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const CaloHitList &caloHitList)
{
    return (position - LArClusterHelper::GetClosestPosition(position, caloHitList)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const ClusterList &clusterList)
{
    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());
    CartesianVector closestPosition(0.f, 0.f, 0.f);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pTestCluster = *iter;
        const CartesianVector thisPosition(LArClusterHelper::GetClosestPosition(position, pTestCluster));
        const float thisDistanceSquared((position - thisPosition).GetMagnitudeSquared());

        if (thisDistanceSquared < closestDistanceSquared)
        {
            distanceFound = true;
            closestDistanceSquared = thisDistanceSquared;
            closestPosition = thisPosition;
        }
    }

    if (distanceFound)
        return closestPosition;

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const Cluster *const pCluster)
{
    return LArClusterHelper::GetClosestPosition(position, pCluster->GetOrderedCaloHitList());
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const OrderedCaloHitList &caloHitList)
{
    const CaloHit *pClosestCaloHit(nullptr);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (const auto &entry : caloHitList)
    {
        for (const CaloHit *const pCaloHit : *entry.second)
        {
            const float distanceSquared((pCaloHit->GetPositionVector() - position).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                pClosestCaloHit = pCaloHit;
            }
        }
    }

    if (pClosestCaloHit)
        return pClosestCaloHit->GetPositionVector();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const CaloHitList &caloHitList)
{
    const CaloHit *pClosestCaloHit(nullptr);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distanceSquared((pCaloHit->GetPositionVector() - position).GetMagnitudeSquared());

        if (distanceSquared < closestDistanceSquared)
        {
            closestDistanceSquared = distanceSquared;
            pClosestCaloHit = pCaloHit;
        }
    }

    if (pClosestCaloHit)
        return pClosestCaloHit->GetPositionVector();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClosestPositions(
    const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianVector &outputPosition1, CartesianVector &outputPosition2)
{
    bool distanceFound(false);
    float minDistanceSquared(std::numeric_limits<float>::max());

    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    const OrderedCaloHitList &orderedCaloHitList1(pCluster1->GetOrderedCaloHitList());
    const OrderedCaloHitList &orderedCaloHitList2(pCluster2->GetOrderedCaloHitList());

    // Loop over hits in cluster 1
    for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList1.begin(), iter1End = orderedCaloHitList1.end(); iter1 != iter1End; ++iter1)
    {
        for (CaloHitList::const_iterator hitIter1 = iter1->second->begin(), hitIter1End = iter1->second->end(); hitIter1 != hitIter1End; ++hitIter1)
        {
            const CartesianVector &positionVector1((*hitIter1)->GetPositionVector());

            // Loop over hits in cluster 2
            for (OrderedCaloHitList::const_iterator iter2 = orderedCaloHitList2.begin(), iter2End = orderedCaloHitList2.end(); iter2 != iter2End; ++iter2)
            {
                for (CaloHitList::const_iterator hitIter2 = iter2->second->begin(), hitIter2End = iter2->second->end(); hitIter2 != hitIter2End; ++hitIter2)
                {
                    const CartesianVector &positionVector2((*hitIter2)->GetPositionVector());

                    const float distanceSquared((positionVector1 - positionVector2).GetMagnitudeSquared());

                    if (distanceSquared < minDistanceSquared)
                    {
                        minDistanceSquared = distanceSquared;
                        closestPosition1 = positionVector1;
                        closestPosition2 = positionVector2;
                        distanceFound = true;
                    }
                }
            }
        }
    }

    if (!distanceFound)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    outputPosition1 = closestPosition1;
    outputPosition2 = closestPosition2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClusterBoundingBox(const Cluster *const pCluster, CartesianVector &minimumCoordinate, CartesianVector &maximumCoordinate,
    bool useHitWidth)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float xmin(std::numeric_limits<float>::max());
    float ymin(std::numeric_limits<float>::max());
    float zmin(std::numeric_limits<float>::max());
    float xmax(-std::numeric_limits<float>::max());
    float ymax(-std::numeric_limits<float>::max());
    float zmax(-std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());
            const float halfWidth = useHitWidth ? 0.5f * pCaloHit->GetCellSize1() : 0.f;
            xmin = std::min(hit.GetX() - halfWidth, xmin);
            xmax = std::max(hit.GetX() + halfWidth, xmax);
            ymin = std::min(hit.GetY(), ymin);
            ymax = std::max(hit.GetY(), ymax);
            zmin = std::min(hit.GetZ(), zmin);
            zmax = std::max(hit.GetZ(), zmax);
        }
    }

    minimumCoordinate.SetValues(xmin, ymin, zmin);
    maximumCoordinate.SetValues(xmax, ymax, zmax);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArClusterHelper::GetAverageZ(const Cluster *const pCluster, const float xmin, const float xmax, float &averageZ)
{
    averageZ = std::numeric_limits<float>::max();

    if (xmin > xmax)
        return STATUS_CODE_INVALID_PARAMETER;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float zsum(0.f);
    int count(0);

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());

            if (hit.GetX() < xmin || hit.GetX() > xmax)
                continue;

            zsum += hit.GetZ();
            ++count;
        }
    }

    if (count == 0)
        return STATUS_CODE_NOT_FOUND;

    averageZ = zsum / static_cast<float>(count);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const ClusterList &clusterList, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
    OrderedCaloHitList orderedCaloHitList;

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(pCluster->GetOrderedCaloHitList()));
    }

    return LArClusterHelper::GetExtremalCoordinates(orderedCaloHitList, innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const Cluster *const pCluster, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
    return LArClusterHelper::GetExtremalCoordinates(pCluster->GetOrderedCaloHitList(), innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const OrderedCaloHitList &orderedCaloHitList, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointVector coordinateVector;

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;
            coordinateVector.push_back(pCaloHit->GetPositionVector());
        }
    }

    std::sort(coordinateVector.begin(), coordinateVector.end(), LArClusterHelper::SortCoordinatesByPosition);
    return LArClusterHelper::GetExtremalCoordinates(coordinateVector, innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const CartesianPointVector &coordinateVector, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
    if (coordinateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Find the extremal values of the X, Y and Z coordinates
    float xMin(+std::numeric_limits<float>::max());
    float yMin(+std::numeric_limits<float>::max());
    float zMin(+std::numeric_limits<float>::max());
    float xMax(-std::numeric_limits<float>::max());
    float yMax(-std::numeric_limits<float>::max());
    float zMax(-std::numeric_limits<float>::max());

    for (CartesianPointVector::const_iterator pIter = coordinateVector.begin(), pIterEnd = coordinateVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &pos = *pIter;
        xMin = std::min(pos.GetX(), xMin);
        xMax = std::max(pos.GetX(), xMax);
        yMin = std::min(pos.GetY(), yMin);
        yMax = std::max(pos.GetY(), yMax);
        zMin = std::min(pos.GetZ(), zMin);
        zMax = std::max(pos.GetZ(), zMax);
    }

    // Choose the coordinate with the greatest span (keeping any ties)
    const float xAve(0.5f * (xMin + xMax));
    const float yAve(0.5f * (yMin + yMax));
    const float zAve(0.5f * (zMin + zMax));

    const float xSpan(xMax - xMin);
    const float ySpan(yMax - yMin);
    const float zSpan(zMax - zMin);

    const bool useX((xSpan > std::numeric_limits<float>::epsilon()) && (xSpan + std::numeric_limits<float>::epsilon() > std::max(ySpan, zSpan)));
    const bool useY((ySpan > std::numeric_limits<float>::epsilon()) && (ySpan + std::numeric_limits<float>::epsilon() > std::max(zSpan, xSpan)));
    const bool useZ((zSpan > std::numeric_limits<float>::epsilon()) && (zSpan + std::numeric_limits<float>::epsilon() > std::max(xSpan, ySpan)));

    // Find the extremal hits separately for the chosen coordinates
    CartesianPointVector candidateVector;

    for (CartesianPointVector::const_iterator pIter = coordinateVector.begin(), pIterEnd = coordinateVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &pos = *pIter;

        if (useX)
        {
            if (((pos.GetX() - xMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetX() - xMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }

        if (useY)
        {
            if (((pos.GetY() - yMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetY() - yMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }

        if (useZ)
        {
            if (((pos.GetZ() - zMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetZ() - zMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }
    }

    // Finally, find the pair of hits that are separated by the greatest distance
    CartesianVector firstCoordinate(xAve, yAve, zAve);
    CartesianVector secondCoordinate(xAve, yAve, zAve);
    float maxDistanceSquared(+std::numeric_limits<float>::epsilon());

    for (CartesianPointVector::const_iterator iterI = candidateVector.begin(), iterEndI = candidateVector.end(); iterI != iterEndI; ++iterI)
    {
        const CartesianVector &posI = *iterI;

        for (CartesianPointVector::const_iterator iterJ = iterI, iterEndJ = candidateVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            const CartesianVector &posJ = *iterJ;

            const float distanceSquared((posI - posJ).GetMagnitudeSquared());

            if (distanceSquared > maxDistanceSquared)
            {
                maxDistanceSquared = distanceSquared;
                firstCoordinate = posI;
                secondCoordinate = posJ;
            }
        }
    }

    // Set the inner and outer coordinates (Check Z first, then X in the event of a tie)
    const float deltaZ(secondCoordinate.GetZ() - firstCoordinate.GetZ());
    const float deltaX(secondCoordinate.GetX() - firstCoordinate.GetX());

    if ((deltaZ > 0.f) || ((std::fabs(deltaZ) < std::numeric_limits<float>::epsilon()) && (deltaX > 0.f)))
    {
        innerCoordinate = firstCoordinate;
        outerCoordinate = secondCoordinate;
    }
    else
    {
        innerCoordinate = secondCoordinate;
        outerCoordinate = firstCoordinate;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetCoordinateVector(const Cluster *const pCluster, CartesianPointVector &coordinateVector)
{
    for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit : *layerEntry.second)
            coordinateVector.push_back(pCaloHit->GetPositionVector());
    }

    std::sort(coordinateVector.begin(), coordinateVector.end(), LArClusterHelper::SortCoordinatesByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetCaloHitListInBoundingBox(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBound,
    const pandora::CartesianVector &upperBound, pandora::CaloHitList &caloHitList)
{
    const bool useX(std::fabs(upperBound.GetX() - lowerBound.GetX()) > std::numeric_limits<float>::epsilon());
    const bool useY(std::fabs(upperBound.GetY() - lowerBound.GetY()) > std::numeric_limits<float>::epsilon());
    const bool useZ(std::fabs(upperBound.GetZ() - lowerBound.GetZ()) > std::numeric_limits<float>::epsilon());
    if (!useX && !useY && !useZ)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float minX(std::min(lowerBound.GetX(), upperBound.GetX()));
    const float maxX(std::max(lowerBound.GetX(), upperBound.GetX()));
    const float minY(std::min(lowerBound.GetY(), upperBound.GetY()));
    const float maxY(std::max(lowerBound.GetY(), upperBound.GetY()));
    const float minZ(std::min(lowerBound.GetZ(), upperBound.GetZ()));
    const float maxZ(std::max(lowerBound.GetZ(), upperBound.GetZ()));

    for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit : *layerEntry.second)
        {
            const CartesianVector &hitPosition = pCaloHit->GetPositionVector();
            if (useX &&
                (hitPosition.GetX() < minX - std::numeric_limits<float>::epsilon() || hitPosition.GetX() > maxX + std::numeric_limits<float>::epsilon()))
                continue;
            else if (useY &&
                (hitPosition.GetY() < minY - std::numeric_limits<float>::epsilon() || hitPosition.GetY() > maxY + std::numeric_limits<float>::epsilon()))
                continue;
            else if (useZ &&
                (hitPosition.GetZ() < minZ - std::numeric_limits<float>::epsilon() || hitPosition.GetZ() > maxZ + std::numeric_limits<float>::epsilon()))
                continue;

            caloHitList.push_back(pCaloHit);
        }
    }
    caloHitList.sort(LArClusterHelper::SortHitsByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetDaughterVolumeIDs(const Cluster *const pCluster, UIntSet &daughterVolumeIds)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit(*hIter);
            const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));

            if (pLArCaloHit)
                daughterVolumeIds.insert(pLArCaloHit->GetDaughterVolumeId());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::HasBlockedPath(const CaloHitVector &caloHits, const CaloHit *const pCaloHit1, const CaloHit *const pCaloHit2)
{
    // For each hit in the calo hit list, check if the line between the two hits passes through the hit
    const CartesianVector &pos1{pCaloHit1->GetPositionVector()};
    const CartesianVector &pos2{pCaloHit2->GetPositionVector()};
    for (const CaloHit *const pCaloHit : caloHits)
    {
        if (pCaloHit == pCaloHit1 || pCaloHit == pCaloHit2)
            continue;

        const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
        const double xmin{hitPosition.GetX() - 0.5 * pCaloHit->GetCellSize1()}, xmax{hitPosition.GetX() + 0.5 * pCaloHit->GetCellSize1()};
        const double zmin{hitPosition.GetZ() - 0.5 * pCaloHit->GetCellSize0()}, zmax{hitPosition.GetZ() + 0.5 * pCaloHit->GetCellSize0()};

        double entry{0.}, exit{1.};

        auto check_axis = [&](double p1, double p2, double minB, double maxB)
        {
            double t1{(minB - p1) / (p2 - p1)};
            double t2{(maxB - p1) / (p2 - p1)};

            if (t1 > t2)
                std::swap(t1, t2);

            entry = std::max(entry, t1);
            exit = std::min(exit, t2);

            return entry <= exit;
        };

        if (check_axis(pos1.GetX(), pos2.GetX(), xmin, xmax) && check_axis(pos1.GetZ(), pos2.GetZ(), zmin, zmax))
        {
            // We have an intervening hit
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::OrderHitsAlongTrajectory(const Cluster *const pCluster, const TwoDSlidingFitResult &sfr, CaloHitList &orderedHits)
{
    CaloHitList clusterHits;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHits);

    if (clusterHits.size() < 3)
    {
        orderedHits = clusterHits;
        return;
    }

    // Build a KD tree to allow quick identification of hits in the vicinity of a given position along the trajectory
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    float r{2.f};
    HitKDTree2D kdTree;
    HitKDNode2DList nodes;
    KDTreeBox boundingRegion{fill_and_bound_2d_kd_tree(clusterHits, nodes)};
    kdTree.build(nodes, boundingRegion);

    // We'll compare the size of the number of used hits to the number of hits in the cluster to determine when we're done
    const size_t nHits(clusterHits.size());
    CaloHitSet usedHits, prevSeed;
    // Start with the first hit in the list
    const CaloHit *pSeedHit{clusterHits.front()};
    bool seedIdentified{false};
    // Check that no other hit precedes this in the sliding fit trajectory, and if it does, update the seed hit until we find the true seed
    while (!seedIdentified)
    {
        seedIdentified = true;
        // Get sliding fit trajectory for the seed
        const CartesianVector &seedPos{pSeedHit->GetPositionVector()};
        float rL{0.f}, rT{0.f};
        CartesianVector fitDir(0, 0, 0);
        sfr.GetLocalPosition(seedPos, rL, rT);
        sfr.GetGlobalFitDirection(rL, fitDir);

        for (const CaloHit *pTestHit : clusterHits)
        {
            const CartesianVector &testPos{pTestHit->GetPositionVector()};
            const CartesianVector testDir{testPos - seedPos};
            const float lPos{testDir.GetDotProduct(fitDir)};

            // Check for upstream hit, and if found, select as new seed and repeat the process
            if (lPos < 0.f && !prevSeed.count(pTestHit))
            {
                prevSeed.insert(pSeedHit);
                pSeedHit = pTestHit;
                seedIdentified = false;
                break;
            }
        }
    }
    orderedHits.push_back(pSeedHit);
    usedHits.insert(pSeedHit);

    // Add hits along the trajectory
    while (usedHits.size() < nHits)
    {
        r = 2.f;
        // Get the local sliding fit trajectory
        const CartesianVector &seedPos{pSeedHit->GetPositionVector()};
        float rL{0.f}, rT{0.f};
        CartesianVector fitDir(0, 0, 0);
        sfr.GetLocalPosition(seedPos, rL, rT);
        sfr.GetGlobalFitDirection(rL, fitDir);

        // Find the nearby hits
        HitKDNode2DList found;
        bool ok{false};
        unsigned int nExpansions{0};
        do
        {
            // If we've failed to find usable hits after three expansions, just include all hits
            if (nExpansions > 2)
                r = std::numeric_limits<float>::max();
            KDTreeBox searchRegion(build_2d_kd_search_region(pSeedHit, r, r));
            kdTree.search(searchRegion, found);
            // Check if we found something other than the seed
            if (found.size() < 2)
            {
                r *= 2.f;
                ++nExpansions;
            }
            else
            {
                for (const auto &node : found)
                {
                    const CaloHit *const pTestHit{node.data};
                    // Ensure our hit region has usable hits
                    if (!usedHits.count(pTestHit))
                    {
                        ok = true;
                        break;
                    }
                }
                if (!ok)
                {
                    found.clear();
                    r *= 2.f;
                    ++nExpansions;
                }
            }
        } while (!ok);

        // Find the most upstream unused hit relative to the seed hit, and set this to be the new seed
        float minLPos(std::numeric_limits<float>::max());
        for (const auto &node : found)
        {
            const CaloHit *const pTestHit{node.data};
            if (usedHits.count(pTestHit))
                continue;
            const CartesianVector &testPos{pTestHit->GetPositionVector()};
            const CartesianVector testDir{testPos - seedPos};
            const float lPos{testDir.GetDotProduct(fitDir)};
            if (lPos < minLPos)
            {
                minLPos = lPos;
                pSeedHit = pTestHit;
            }
        }
        usedHits.insert(pSeedHit);
        orderedHits.emplace_back(pSeedHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNOccupiedLayers(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nOccupiedLayersLhs(pLhs->GetOrderedCaloHitList().size());
    const unsigned int nOccupiedLayersRhs(pRhs->GetOrderedCaloHitList().size());

    if (nOccupiedLayersLhs != nOccupiedLayersRhs)
        return (nOccupiedLayersLhs > nOccupiedLayersRhs);

    return SortByNHits(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNHits(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nHitsLhs(pLhs->GetNCaloHits());
    const unsigned int nHitsRhs(pRhs->GetNCaloHits());

    if (nHitsLhs != nHitsRhs)
        return (nHitsLhs > nHitsRhs);

    return SortByLayerSpan(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByLayerSpan(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return SortByInnerLayer(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByInnerLayer(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int innerLayerLhs(pLhs->GetInnerPseudoLayer());
    const unsigned int innerLayerRhs(pRhs->GetInnerPseudoLayer());

    if (innerLayerLhs != innerLayerRhs)
        return (innerLayerLhs < innerLayerRhs);

    return SortByPosition(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByPosition(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const CartesianVector deltaPositionIL(pRhs->GetCentroid(pRhs->GetInnerPseudoLayer()) - pLhs->GetCentroid(pLhs->GetInnerPseudoLayer()));

    if (std::fabs(deltaPositionIL.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionIL.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionIL.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionIL.GetX() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionIL.GetY()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionIL.GetY() > std::numeric_limits<float>::epsilon());

    const CartesianVector deltaPositionOL(pRhs->GetCentroid(pRhs->GetOuterPseudoLayer()) - pLhs->GetCentroid(pLhs->GetOuterPseudoLayer()));

    if (std::fabs(deltaPositionOL.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionOL.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionOL.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionOL.GetX() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionOL.GetY()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionOL.GetY() > std::numeric_limits<float>::epsilon());

    // Use pulse height to resolve ties
    return SortByPulseHeight(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByPulseHeight(const Cluster *const pLhs, const Cluster *const pRhs)
{
    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortHitsByPosition(const CaloHit *const pLhs, const CaloHit *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPositionVector() - pLhs->GetPositionVector());

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetY()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());

    // Use pulse height to resolve ties
    return SortHitsByPulseHeight(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortHitsByPositionInX(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPositionVector() - pLhs->GetPositionVector());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    return SortHitsByPosition(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortHitsByPulseHeight(const CaloHit *const pLhs, const CaloHit *const pRhs)
{
    // TODO: Think about the correct energy to use here (should they ever be different)
    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortCoordinatesByPosition(const CartesianVector &lhs, const CartesianVector &rhs)
{
    constexpr float epsilon(std::numeric_limits<float>::epsilon());

    const float deltaZ(rhs.GetZ() - lhs.GetZ());
    if (std::fabs(deltaZ) > epsilon)
        return (deltaZ > epsilon);

    const float deltaX(rhs.GetX() - lhs.GetX());
    if (std::fabs(deltaX) > epsilon)
        return (deltaX > epsilon);

    return (rhs.GetY() - lhs.GetY() > epsilon);
}

} // namespace lar_content
