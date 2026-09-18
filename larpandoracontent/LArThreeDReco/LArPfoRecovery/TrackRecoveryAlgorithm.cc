/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.cc
 *
 *  @brief  Implementation of the particle recovery algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

TrackRecoveryAlgorithm::TrackRecoveryAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::Run()
{
    // Get the list of track-like PFOs
    const PfoList *pPfoList(nullptr);
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList))
        return STATUS_CODE_SUCCESS;

    // Get maps between PFOs, views, clusters and hits
    PfoToViewClusterMap pfoToViewClusterMap;
    ClusterToPfoMap clusterToPfoMap;
    ClusterToHitMap clusterToHitMap;
    ClusterToFitMap clusterToFitMap;
    HitToClusterToMap hitToClusterMap;
    for (const Pfo *const pPfo : *pPfoList)
    {
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, TPC_3D})
        {
            ClusterList pfoClusters;
            LArPfoHelper::GetClusters(pPfo, view, pfoClusters);
            if (!pfoClusters.empty())
            {
                const Cluster *pCluster{pfoClusters.front()};
                pfoToViewClusterMap[pPfo][view] = pCluster;
                clusterToPfoMap[pCluster] = pPfo;
                if (view != TPC_3D)
                {
                    auto result = this->FitAndOrderCluster(pCluster, clusterToHitMap[pCluster]);
                    if (result)
                        clusterToFitMap.emplace(pCluster, *result);
                }
                else
                    LArClusterHelper::GetAllHits(pCluster, clusterToHitMap[pCluster]);
                for (const CaloHit *const pCaloHit : clusterToHitMap[pCluster])
                    hitToClusterMap[pCaloHit] = pCluster;
            }
        }
    }
    // Loop over all clusters to catch any that aren't associated to a PFO
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
        for (const Cluster *const pCluster : *pClusterList)
        {
            if (clusterToPfoMap.find(pCluster) == clusterToPfoMap.end())
            {
                LArClusterHelper::GetAllHits(pCluster, clusterToHitMap[pCluster]);
                for (const CaloHit *const pCaloHit : clusterToHitMap[pCluster])
                    hitToClusterMap[pCaloHit] = pCluster;
            }
        }
    }

    // Loop over the PFOs and project 3D hits from two views into the third view to look for unassociated hits
    std::unordered_map<const Pfo *, std::unordered_map<HitType, CaloHitList>> pfoToMergeHitsMap;
    for (const Pfo *const pPfo : *pPfoList)
    {
        const Cluster *const pClusterU{pfoToViewClusterMap[pPfo][TPC_VIEW_U]}, *const pClusterV{pfoToViewClusterMap[pPfo][TPC_VIEW_V]},
                                                                                   *const pClusterW{pfoToViewClusterMap[pPfo][TPC_VIEW_W]};
        const CaloHitList &hitsU{clusterToHitMap[pClusterU]}, &hitsV{clusterToHitMap[pClusterV]}, &hitsW{clusterToHitMap[pClusterW]};

        CaloHitSet unmatchedHitsU, unmatchedHitsV, unmatchedHitsW;
        CaloHitList mergeHitsU, mergeHitsV, mergeHitsW;
        this->FindUnmatchedHits(hitsU, hitsV, hitsW, unmatchedHitsV, unmatchedHitsW);
        this->FindUnmatchedHits(hitsV, hitsU, hitsW, unmatchedHitsU, unmatchedHitsW);
        this->FindUnmatchedHits(hitsW, hitsU, hitsV, unmatchedHitsU, unmatchedHitsV);

        // Construct 3D hits from the combinatorics of each view pair and project into the third view
        this->IdentifyHitsToMerge(pClusterU, hitsU, unmatchedHitsU, clusterToFitMap, mergeHitsU);
        this->IdentifyHitsToMerge(pClusterV, hitsV, unmatchedHitsV, clusterToFitMap, mergeHitsV);
        this->IdentifyHitsToMerge(pClusterW, hitsW, unmatchedHitsW, clusterToFitMap, mergeHitsW);

        auto filterMergeHits = [this](const CaloHitList &hits, CaloHitList &mergeHits, const HitType view)
        {
            CartesianPointVector newPos;
            for (const CaloHit *const pCaloHit : hits)
                newPos.emplace_back(pCaloHit->GetPositionVector());
            for (const CaloHit *const pCaloHit : mergeHits)
                newPos.emplace_back(pCaloHit->GetPositionVector());

            try
            {
                TwoDSlidingFitResult sfr(&newPos, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
                this->FilterHitsToMerge(sfr, mergeHits);
            }
            catch (const StatusCodeException &)
            {
            }
        };

        filterMergeHits(hitsU, mergeHitsU, TPC_VIEW_U);
        filterMergeHits(hitsV, mergeHitsV, TPC_VIEW_V);
        filterMergeHits(hitsW, mergeHitsW, TPC_VIEW_W);

        for (const CaloHit *const pCaloHit : mergeHitsU)
            pfoToMergeHitsMap[pPfo][TPC_VIEW_U].emplace_back(pCaloHit);
        for (const CaloHit *const pCaloHit : mergeHitsV)
            pfoToMergeHitsMap[pPfo][TPC_VIEW_V].emplace_back(pCaloHit);
        for (const CaloHit *const pCaloHit : mergeHitsW)
            pfoToMergeHitsMap[pPfo][TPC_VIEW_W].emplace_back(pCaloHit);
    }

    PfoVector sortedPfos;
    sortedPfos.reserve(pfoToMergeHitsMap.size());
    for (const auto &[pPfo, viewToMergeHitsMap] : pfoToMergeHitsMap)
        sortedPfos.emplace_back(pPfo);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo : sortedPfos)
    {
        auto &viewToMergeHitsMap{pfoToMergeHitsMap.at(pPfo)};
        for (auto &[view, mergeHits] : viewToMergeHitsMap)
        {
            const Cluster *pNewCluster{pfoToViewClusterMap[pPfo][view]};
            // Need to allow for the possibility that the view was entirely missing, in which case we need to create a new cluster to merge into
            for (const CaloHit *const pCaloHit : mergeHits)
            {
                if (hitToClusterMap.find(pCaloHit) != hitToClusterMap.end())
                {
                    const Cluster *pOldCluster{hitToClusterMap[pCaloHit]};
                    PandoraContentApi::RemoveFromCluster(*this, pOldCluster, pCaloHit);
                    CaloHitList &oldClusterHits{clusterToHitMap[pOldCluster]};
                    oldClusterHits.remove(pCaloHit);
                }
                PandoraContentApi::AddToCluster(*this, pNewCluster, pCaloHit);
                hitToClusterMap[pCaloHit] = pNewCluster;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::optional<TwoDSlidingFitResult> TrackRecoveryAlgorithm::FitAndOrderCluster(const Cluster *const pCluster, CaloHitList &orderedHits) const
{
    if (pCluster->GetNCaloHits() < 3)
    {
        LArClusterHelper::GetAllHits(pCluster, orderedHits);
        return std::nullopt;
    }
    const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
    try
    {
        TwoDSlidingFitResult sfr(pCluster, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
        LArClusterHelper::OrderHitsAlongTrajectory(pCluster, sfr, orderedHits);
        return sfr;
    }
    catch (const StatusCodeException &)
    {
        LArClusterHelper::GetAllHits(pCluster, orderedHits);
        return std::nullopt;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FindUnmatchedHits(const CaloHitList &hitsA, const CaloHitList &hitsB, const CaloHitList &hitsC,
    CaloHitSet &unmatchedHitsB, CaloHitSet &unmatchedHitsC) const
{
    if (hitsA.empty() || hitsB.empty() || hitsC.empty())
        return;
    const LArPlaneContextObject *pPlaneContextObject{
        dynamic_cast<const LArPlaneContextObject *>(PandoraContentApi::GetEventContextObject(*this, "PlaneContext"))};
    if (!pPlaneContextObject)
        return;
    const HitType viewB{hitsB.front()->GetHitType()}, viewC{hitsC.front()->GetHitType()};
    for (const CaloHit *const pCaloHitA : hitsA)
    {
        const LArPlaneContextObject::HitTriplet *pTriplet{pPlaneContextObject->GetHitTriplet(pCaloHitA)};
        if (!pTriplet)
            continue;
        const CaloHit *const pHitB{viewB == TPC_VIEW_W ? pTriplet->m_wHit : viewB == TPC_VIEW_V ? pTriplet->m_vHit : pTriplet->m_uHit};
        const CaloHit *const pHitC{viewC == TPC_VIEW_W ? pTriplet->m_wHit : viewC == TPC_VIEW_V ? pTriplet->m_vHit : pTriplet->m_uHit};
        auto bIterator{std::find(hitsB.begin(), hitsB.end(), pHitB)};
        auto cIterator{std::find(hitsC.begin(), hitsC.end(), pHitC)};
        // Only flag hits as unmatched if the other hits in the triplet are present in the PFO
        if (pHitB && bIterator == hitsB.end() && cIterator != hitsC.end())
            unmatchedHitsB.insert(pHitB);
        if (pHitC && cIterator == hitsC.end() && bIterator != hitsB.end())
            unmatchedHitsC.insert(pHitC);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::IdentifyHitsToMerge(const Cluster *pCluster, const CaloHitList &clusterHits, const CaloHitSet &unmatchedHits,
    const ClusterToFitMap &clusterToFitMap, CaloHitList &mergeHits) const
{
    CaloHitVector sortedHits;
    std::copy(unmatchedHits.begin(), unmatchedHits.end(), std::back_inserter(sortedHits));
    std::sort(sortedHits.begin(), sortedHits.end(), LArClusterHelper::SortHitsByPosition);

    if (pCluster && (clusterHits.size() > 1) && (clusterToFitMap.find(pCluster) != clusterToFitMap.end()))
    {
        float rLFront{0.f}, rLBack{0.f}, dummy{0.f};
        clusterToFitMap.at(pCluster).GetLocalPosition(clusterHits.front()->GetPositionVector(), rLFront, dummy);
        clusterToFitMap.at(pCluster).GetLocalPosition(clusterHits.back()->GetPositionVector(), rLBack, dummy);
        for (const CaloHit *const pCaloHit : sortedHits)
        {
            float rL{0.f};
            clusterToFitMap.at(pCluster).GetLocalPosition(pCaloHit->GetPositionVector(), rL, dummy);
            if (rLFront <= rL && rL <= rLBack)
            {
                // We're along the track, look for this hit acting as a blocking path
                for (auto iter = clusterHits.begin(); iter != std::prev(clusterHits.end()); ++iter)
                {
                    const CaloHit *const pHit1{*iter}, *const pHit2{*(std::next(iter))};
                    if (LArClusterHelper::HasBlockedPath({pCaloHit}, pHit1, pHit2))
                    {
                        mergeHits.emplace_back(pCaloHit);
                        break;
                    }
                }
            }
            else
            {
                // We're outside the track, extend it
                mergeHits.emplace_back(pCaloHit);
            }
        }
    }
    else
    {
        // Entire view is missing, or has less than 3 hits, so we should merge all unmatched hits
        for (const CaloHit *const pCaloHit : sortedHits)
            mergeHits.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FilterHitsToMerge(const TwoDSlidingFitResult &sfr, CaloHitList &mergeHits) const
{
    for (auto iter = mergeHits.begin(); iter != mergeHits.end();)
    {
        const CaloHit *const pCaloHit{*iter};
        const CartesianVector position{pCaloHit->GetPositionVector()};
        float rL(0.f), rT(0.f);
        sfr.GetLocalPosition(position, rL, rT);

        CartesianVector fitPosition(0, 0, 0);
        if (STATUS_CODE_SUCCESS != sfr.GetGlobalFitPosition(rL, fitPosition))
        {
            mergeHits.erase(iter++);
            continue;
        }

        float rTFit(0.f), dummy(0.f);
        sfr.GetLocalPosition(fitPosition, dummy, rTFit);

        CartesianVector fitDirection(0, 0, 0);
        if (STATUS_CODE_SUCCESS != sfr.GetGlobalFitDirection(rL, fitDirection))
        {
            mergeHits.erase(iter++);
            continue;
        }

        float dTdL(0.f);
        sfr.GetLocalDirection(fitDirection, dTdL);

        const float residual{(rT - rTFit) / std::sqrt(1.f + dTdL * dTdL)};
        if (std::abs(residual) > 0.5f)
            mergeHits.erase(iter++);
        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
