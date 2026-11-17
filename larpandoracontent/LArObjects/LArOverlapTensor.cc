/**
 *  @file   larpandoracontent/LArObjects/LArOverlapTensor.cc
 *
 *  @brief  Implementation of the lar overlap tensor class.
 *
 *  $Log: $
 */

#include "Pandora/PandoraInputTypes.h"

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "Objects/Cluster.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArOverlapTensor.h"
#include "larpandoracontent/LArObjects/LArShowerOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

template <typename T>
const Cluster *OverlapTensor<T>::Element::GetCluster(const HitType hitType) const
{
    if ((hitType != TPC_VIEW_U) && (hitType != TPC_VIEW_V) && (hitType != TPC_VIEW_W))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    return (hitType == TPC_VIEW_U) ? m_pClusterU : (hitType == TPC_VIEW_V) ? m_pClusterV : m_pClusterW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetUnambiguousElements(const bool ignoreUnavailable, ElementList &elementList) const
{
    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
    {
        ClusterList clusterListU, clusterListV, clusterListW;
        this->GetConnectedClusters(iterU->first, ignoreUnavailable, nullptr, clusterListU, clusterListV, clusterListW);

        const Cluster *pClusterU(nullptr), *pClusterV(nullptr), *pClusterW(nullptr);
        if (!this->DefaultAmbiguityFunction(clusterListU, clusterListV, clusterListW, pClusterU, pClusterV, pClusterW))
            continue;

        // ATTN With HIT_CUSTOM definitions, it is possible to navigate from different U clusters to same combination
        if (iterU->first != pClusterU)
            continue;

        if (!pClusterU || !pClusterV || !pClusterW)
            continue;

        typename OverlapMatrix::const_iterator iterV = iterU->second.find(pClusterV);
        if (iterU->second.end() == iterV)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        typename OverlapList::const_iterator iterW = iterV->second.find(pClusterW);
        if (iterV->second.end() == iterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Element element(pClusterU, pClusterV, pClusterW, iterW->second);
        elementList.push_back(element);
    }

    std::sort(elementList.begin(), elementList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool OverlapTensor<T>::DefaultAmbiguityFunction(const ClusterList &clusterListU, const ClusterList &clusterListV,
    const ClusterList &clusterListW, const Cluster *&pClusterU, const Cluster *&pClusterV, const Cluster *&pClusterW) const
{
    if ((1 != clusterListU.size()) || (1 != clusterListV.size()) || (1 != clusterListW.size()))
        return false;

    pClusterU = *(clusterListU.begin());
    pClusterV = *(clusterListV.begin());
    pClusterW = *(clusterListW.begin());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedElements(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
    unsigned int &nU, unsigned int &nV, unsigned int &nW) const
{
    ClusterList clusterListU, clusterListV, clusterListW;
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, clusterListU, clusterListV, clusterListW);
    nU = clusterListU.size();
    nV = clusterListV.size();
    nW = clusterListW.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetSortedKeyClusters(ClusterVector &sortedKeyClusters) const
{
    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
        sortedKeyClusters.push_back(iterU->first);

    std::sort(sortedKeyClusters.begin(), sortedKeyClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::SetOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
    const pandora::Cluster *const pClusterW, const OverlapResult &overlapResult)
{
    OverlapList &overlapList = m_overlapTensor[pClusterU][pClusterV];
    typename OverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() != iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!overlapList.insert(typename OverlapList::value_type(pClusterW, overlapResult)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    ClusterList &navigationUV(m_clusterNavigationMapUV[pClusterU]);
    ClusterList &navigationVW(m_clusterNavigationMapVW[pClusterV]);
    ClusterList &navigationWU(m_clusterNavigationMapWU[pClusterW]);

    if (navigationUV.end() == std::find(navigationUV.begin(), navigationUV.end(), pClusterV))
        navigationUV.push_back(pClusterV);
    if (navigationVW.end() == std::find(navigationVW.begin(), navigationVW.end(), pClusterW))
        navigationVW.push_back(pClusterW);
    if (navigationWU.end() == std::find(navigationWU.begin(), navigationWU.end(), pClusterU))
        navigationWU.push_back(pClusterU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::ReplaceOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
    const pandora::Cluster *const pClusterW, const OverlapResult &overlapResult)
{
    typename TheTensor::iterator iterU = m_overlapTensor.find(pClusterU);

    if (m_overlapTensor.end() == iterU)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    typename OverlapMatrix::iterator iterV = iterU->second.find(pClusterV);

    if (iterU->second.end() == iterV)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    typename OverlapList::iterator iterW = iterV->second.find(pClusterW);

    if (iterV->second.end() == iterW)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    iterW->second = overlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::RemoveCluster(const pandora::Cluster *const pCluster)
{
    ClusterList additionalRemovals;

    if (m_clusterNavigationMapUV.erase(pCluster) > 0)
    {
        typename TheTensor::iterator iter = m_overlapTensor.find(pCluster);

        if (m_overlapTensor.end() != iter)
            m_overlapTensor.erase(iter);

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMapWU.begin(); navIter != m_clusterNavigationMapWU.end();)
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = std::find(thisIter->second.begin(), thisIter->second.end(), pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.push_back(thisIter->first);
        }
    }

    if (m_clusterNavigationMapVW.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            typename OverlapMatrix::iterator iter = iterU->second.find(pCluster);

            if (iterU->second.end() != iter)
                iterU->second.erase(iter);
        }

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMapUV.begin(); navIter != m_clusterNavigationMapUV.end();)
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = std::find(thisIter->second.begin(), thisIter->second.end(), pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.push_back(thisIter->first);
        }
    }

    if (m_clusterNavigationMapWU.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            for (typename OverlapMatrix::iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
            {
                typename OverlapList::iterator iter = iterV->second.find(pCluster);

                if (iterV->second.end() != iter)
                    iterV->second.erase(iter);
            }
        }

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMapVW.begin(); navIter != m_clusterNavigationMapVW.end();)
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = std::find(thisIter->second.begin(), thisIter->second.end(), pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.push_back(thisIter->first);
        }
    }

    additionalRemovals.sort(LArClusterHelper::SortByNHits);

    for (ClusterList::const_iterator iter = additionalRemovals.begin(), iterEnd = additionalRemovals.end(); iter != iterEnd; ++iter)
        this->RemoveCluster(*iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedElements(const Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
    ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    this->GetConnectedClusters(pCluster, ignoreUnavailable, &elementList, clusterListU, clusterListV, clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedClusters(const Cluster *const pCluster, const bool ignoreUnavailable, ElementList *const pElementList,
    ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    ClusterSet connectedClusters;
    this->ExploreConnections(pCluster, ignoreUnavailable, connectedClusters);

    // ATTN Now need to check that all clusters received are from fully available tensor elements
    if (pElementList)
        pElementList->clear();

    clusterListU.clear();
    clusterListV.clear();
    clusterListW.clear();

    ClusterSet addedU, addedV, addedW;

    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
    {
        if (!connectedClusters.count(iterU->first))
            continue;

        for (typename OverlapMatrix::const_iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
        {
            for (typename OverlapList::const_iterator iterW = iterV->second.begin(), iterWEnd = iterV->second.end(); iterW != iterWEnd; ++iterW)
            {
                if (ignoreUnavailable && (!iterU->first->IsAvailable() || !iterV->first->IsAvailable() || !iterW->first->IsAvailable()))
                    continue;

                if (pElementList)
                    pElementList->emplace_back(iterU->first, iterV->first, iterW->first, iterW->second);

                if (addedU.insert(iterU->first).second)
                    clusterListU.push_back(iterU->first);
                if (addedV.insert(iterV->first).second)
                    clusterListV.push_back(iterV->first);
                if (addedW.insert(iterW->first).second)
                    clusterListW.push_back(iterW->first);
            }
        }
    }

    if (pElementList)
        std::sort(pElementList->begin(), pElementList->end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::ExploreConnections(const Cluster *const pCluster, const bool ignoreUnavailable, ClusterSet &connectedClusters) const
{
    if (ignoreUnavailable && !pCluster->IsAvailable())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ClusterNavigationMap &navigationMap((TPC_VIEW_U == hitType) ? m_clusterNavigationMapUV
            : (TPC_VIEW_V == hitType)                                 ? m_clusterNavigationMapVW
                                                                      : m_clusterNavigationMapWU);

    if (!connectedClusters.insert(pCluster).second)
        return;

    ClusterNavigationMap::const_iterator iter = navigationMap.find(pCluster);

    if (navigationMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterList::const_iterator cIter = iter->second.begin(), cIterEnd = iter->second.end(); cIter != cIterEnd; ++cIter)
        this->ExploreConnections(*cIter, ignoreUnavailable, connectedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class OverlapTensor<float>;
template class OverlapTensor<TransverseOverlapResult>;
template class OverlapTensor<LongitudinalOverlapResult>;
template class OverlapTensor<FragmentOverlapResult>;
template class OverlapTensor<ShowerOverlapResult>;
template class OverlapTensor<DeltaRayOverlapResult>;

} // namespace lar_content
