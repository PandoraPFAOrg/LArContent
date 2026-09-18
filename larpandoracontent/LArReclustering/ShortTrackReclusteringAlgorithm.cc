/**
 *  @file   larpandoracontent/LArReclustering/ShortTrackReclusteringAlgorithm.cc
 *
 *  @brief  Tries to identify short tracks that are either lost, or under-clustered and recover missing hits.
 *          This algorithm looks at all three views of a PFO to try to identify if one or more views exhibit step changes in ADC values that
 *          might be indicative of a decay or intelastic interaction. If such a change occurs it looks for consistency across views and also
 *          examines nearby unclustered hits, or hits in nearby clusters to see if a more cohenrent clustering can be identified.
 *          If these various conditions are met, then reclustering is performed.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArReclustering/ShortTrackReclusteringAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

ShortTrackReclusteringAlgorithm::ShortTrackReclusteringAlgorithm() :
    m_maxHitDiscrepancy{3.f},
    m_maxHitDiscrepancySquared{m_maxHitDiscrepancy * m_maxHitDiscrepancy},
    m_balanceThresholdLow{0.67f},
    m_balanceThresholdHigh{1.5f},
    m_braggLinearSlopeThreshold{0.16f},
    m_braggCurvatureThreshold{0.f},
    m_braggContrastThreshold{1.4f},
    m_braggMonotonicityThreshold{0.5f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShortTrackReclusteringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_IF(STATUS_CODE_SUCCESS, !this->GetList(m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList{nullptr};
    if (!this->GetList(m_pfoListName, pPfoList))
        return STATUS_CODE_SUCCESS;

    ViewToClustersMap viewToClustersMap;
    ClusterToPfoMap clusterToPfoMap;
    this->CollectClusters(*pPfoList, viewToClustersMap, clusterToPfoMap);

    this->FitAndOrderClusters(viewToClustersMap);

    // Loop over clusters, and look for evidence of discontinuous changes in ADC deposition and collect the corresponding hits.
    ClusterToHitsMap clusterToHitsMap;
    this->FindAdcDiscontinuities(clusterToPfoMap, clusterToHitsMap);

    // Find corresponding hits for discontinuity hits in other views.
    PfoToHitTripletsMap pfoToHitTripletsMap;
    this->MatchAdcDiscontinuities(clusterToHitsMap, clusterToPfoMap, pfoToHitTripletsMap);

    PartitionVector partitions;
    this->PartitionDiscontinuities(pfoToHitTripletsMap, partitions);
    this->FilterPartitions(partitions);

    if (!partitions.empty())
        this->Recluster(partitions);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool ShortTrackReclusteringAlgorithm::GetList(const std::string &listName, const T *&pList) const
{
    PandoraContentApi::GetList(*this, listName, pList);
    return pList && !pList->empty();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::CollectClusters(const PfoList &pfoList, ViewToClustersMap &viewToClustersMap, ClusterToPfoMap &clusterToPfoMap) const
{
    for (const Pfo *const pPfo : pfoList)
    {
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            ClusterList pfoClusterList;
            LArPfoHelper::GetClusters(pPfo, view, pfoClusterList);
            for (const Cluster *const pCluster : pfoClusterList)
            {
                viewToClustersMap[view].emplace_back(pCluster);
                clusterToPfoMap[pCluster] = pPfo;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::FitAndOrderClusters(const ViewToClustersMap &viewToClustersMap)
{
    for (const auto &[view, clusters] : viewToClustersMap)
    {
        for (const Cluster *const pCluster : clusters)
        {
            if (pCluster->GetNCaloHits() < 3)
            {
                m_clusterToOrderedHitsMap[pCluster] = CaloHitList();
                pCluster->GetOrderedCaloHitList().FillCaloHitList(m_clusterToOrderedHitsMap[pCluster]);
                continue;
            }
            try
            {
                m_clusterToSFRMap.emplace(pCluster, TwoDSlidingFitResult(pCluster, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view)));
                LArClusterHelper::OrderHitsAlongTrajectory(pCluster, m_clusterToSFRMap.at(pCluster), m_clusterToOrderedHitsMap[pCluster]);
            }
            catch (const StatusCodeException &)
            {
                m_clusterToOrderedHitsMap[pCluster] = CaloHitList();
                pCluster->GetOrderedCaloHitList().FillCaloHitList(m_clusterToOrderedHitsMap[pCluster]);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::FindAdcDiscontinuities(const ClusterToPfoMap &clusterToPfoMap, ClusterToHitsMap &clusterToHitsMap) const
{
    for (const auto &[pCluster, _] : clusterToPfoMap)
    {
        CaloHitList clusterHitList{m_clusterToOrderedHitsMap.at(pCluster)};
        // Can't perform the pointing cluster's sliding linear fit without at least 3 hits
        if (clusterHitList.size() < 3)
            continue;

        CaloHitVector clusterHits(clusterHitList.begin(), clusterHitList.end());
        IntVector discontinuities;
        this->GetStableAdcDiscontinuities(clusterHits, discontinuities);
        for (const int index : discontinuities)
            clusterToHitsMap[pCluster].insert(clusterHits.at(index));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::MatchAdcDiscontinuities(
    const ClusterToHitsMap &clusterToHitsMap, const ClusterToPfoMap &clusterToPfoMap, PfoToHitTripletsMap &pfoToHitTripletsMap) const
{
    ClusterVector sortedClusters;
    sortedClusters.reserve(clusterToHitsMap.size());
    for (const auto &[pCluster, _] : clusterToHitsMap)
        sortedClusters.emplace_back(pCluster);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {
        const auto &discontinuityHits{clusterToHitsMap.at(pCluster)};
        const Pfo *const pPfo{clusterToPfoMap.at(pCluster)};
        CaloHitList caloHits3D;
        CaloHitVector caloHits3Du, caloHits3Dv, caloHits3Dw;
        LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHits3D);
        // Split 3D hits according to which 2D view their parent belongs
        for (const CaloHit *const pCaloHit : caloHits3D)
        {
            const CaloHit *pParent{static_cast<const CaloHit *>(pCaloHit->GetParentAddress())};
            if (!pParent)
                continue;
            switch (pParent->GetHitType())
            {
                case TPC_VIEW_U:
                    caloHits3Du.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    caloHits3Dv.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    caloHits3Dw.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }

        for (const CaloHit *const pCaloHit : caloHits3D)
        {
            const CaloHit *pParent{static_cast<const CaloHit *>(pCaloHit->GetParentAddress())};
            if (!pParent)
                continue;
            if (discontinuityHits.find(pParent) == discontinuityHits.end())
                continue;

            CaloHitVector filteredHits3Du, filteredHits3Dv, filteredHits3Dw;
            if (pParent->GetHitType() != TPC_VIEW_U)
            {
                for (const CaloHit *const pCaloHitU : caloHits3Du)
                {
                    const CaloHit *pParentU{static_cast<const CaloHit *>(pCaloHitU->GetParentAddress())};
                    if (!pParentU)
                        continue;
                    const float xU{pCaloHitU->GetPositionVector().GetX()}, &x{pCaloHit->GetPositionVector().GetX()};
                    const float dxU{pCaloHitU->GetCellSize1()}, dx{pCaloHit->GetCellSize1()};
                    if (((xU - dxU) <= x && x <= (xU + dxU)) || ((x - dx) <= xU && xU <= (x + dx)))
                        filteredHits3Du.emplace_back(pCaloHitU);
                }
            }
            if (pParent->GetHitType() != TPC_VIEW_V)
            {
                for (const CaloHit *const pCaloHitV : caloHits3Dv)
                {
                    const CaloHit *pParentV{static_cast<const CaloHit *>(pCaloHitV->GetParentAddress())};
                    if (!pParentV)
                        continue;
                    const float xV{pCaloHitV->GetPositionVector().GetX()}, &x{pCaloHit->GetPositionVector().GetX()};
                    const float dxV{pCaloHitV->GetCellSize1()}, dx{pCaloHit->GetCellSize1()};
                    if (((xV - dxV) <= x && x <= (xV + dxV)) || ((x - dx) <= xV && xV <= (x + dx)))
                        filteredHits3Dv.emplace_back(pCaloHitV);
                }
            }
            if (pParent->GetHitType() != TPC_VIEW_W)
            {
                for (const CaloHit *const pCaloHitW : caloHits3Dw)
                {
                    const CaloHit *pParentW{static_cast<const CaloHit *>(pCaloHitW->GetParentAddress())};
                    if (!pParentW)
                        continue;
                    const float xW{pCaloHitW->GetPositionVector().GetX()}, &x{pCaloHit->GetPositionVector().GetX()};
                    const float dxW{pCaloHitW->GetCellSize1()}, dx{pCaloHit->GetCellSize1()};
                    if (((xW - dxW) <= x && x <= (xW + dxW)) || ((x - dx) <= xW && xW <= (x + dx)))
                        filteredHits3Dw.emplace_back(pCaloHitW);
                }
            }
            // Find the closest matching 3D hit in each view for each discontinuity hit
            Eigen::MatrixXf hitMatrixU(filteredHits3Du.size(), 3);
            LArEigenHelper::Vectorize3D(filteredHits3Du, hitMatrixU);
            Eigen::MatrixXf hitMatrixV(filteredHits3Dv.size(), 3);
            LArEigenHelper::Vectorize3D(filteredHits3Dv, hitMatrixV);
            Eigen::MatrixXf hitMatrixW(filteredHits3Dw.size(), 3);
            LArEigenHelper::Vectorize3D(filteredHits3Dw, hitMatrixW);

            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            Eigen::RowVectorXf row(3);
            row << pos.GetX(), pos.GetY(), pos.GetZ();

            const HitType view{pParent->GetHitType()};
            const CaloHit *pBestU{view == TPC_VIEW_U ? pCaloHit : nullptr}, *pBestV{view == TPC_VIEW_V ? pCaloHit : nullptr},
                *pBestW{view == TPC_VIEW_W ? pCaloHit : nullptr};
            if (view == TPC_VIEW_U || view == TPC_VIEW_V)
            {
                if (filteredHits3Dw.empty())
                    continue;
                Eigen::MatrixXf norms((hitMatrixW.rowwise() - row).array().pow(2).rowwise().sum());
                Eigen::Index index;
                norms.col(0).minCoeff(&index);
                pBestW = filteredHits3Dw.at(index);
            }
            if (view == TPC_VIEW_U || view == TPC_VIEW_W)
            {
                if (filteredHits3Dv.empty())
                    continue;
                Eigen::MatrixXf norms((hitMatrixV.rowwise() - row).array().pow(2).rowwise().sum());
                Eigen::Index index;
                norms.col(0).minCoeff(&index);
                pBestV = filteredHits3Dv.at(index);
            }
            if (view == TPC_VIEW_V || view == TPC_VIEW_W)
            {
                if (filteredHits3Du.empty())
                    continue;
                Eigen::MatrixXf norms((hitMatrixU.rowwise() - row).array().pow(2).rowwise().sum());
                Eigen::Index index;
                norms.col(0).minCoeff(&index);
                pBestU = filteredHits3Du.at(index);
            }

            // Calculate the distances between the found 3D hits and reject if the distances are too large
            const CartesianVector &posU{pBestU->GetPositionVector()}, &posV{pBestV->GetPositionVector()}, &posW{pBestW->GetPositionVector()};
            const float duv{(posU - posV).GetMagnitudeSquared()}, duw{(posU - posW).GetMagnitudeSquared()}, dvw{(posV - posW).GetMagnitudeSquared()};
            const CaloHit *selectedU{static_cast<const CaloHit *>(pBestU->GetParentAddress())},
                *selectedV{static_cast<const CaloHit *>(pBestV->GetParentAddress())},
                *selectedW{static_cast<const CaloHit *>(pBestW->GetParentAddress())};
            if (duv > m_maxHitDiscrepancySquared && duw > m_maxHitDiscrepancySquared)
                selectedU = nullptr;
            if (duv > m_maxHitDiscrepancySquared && dvw > m_maxHitDiscrepancySquared)
                selectedV = nullptr;
            if (duw > m_maxHitDiscrepancySquared && dvw > m_maxHitDiscrepancySquared)
                selectedW = nullptr;

            if (selectedU || selectedV || selectedW)
            {
                // We found candidates
                pfoToHitTripletsMap[pPfo].emplace_back(std::make_tuple(selectedU, selectedV, selectedW));
            }
            else
            {
                // We found no matches, but retain the single 2D discontinuity hit to see if we can recover something later
                selectedU = view == TPC_VIEW_U ? pParent : nullptr;
                selectedV = view == TPC_VIEW_V ? pParent : nullptr;
                selectedW = view == TPC_VIEW_W ? pParent : nullptr;
                pfoToHitTripletsMap[pPfo].emplace_back(std::make_tuple(selectedU, selectedV, selectedW));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::PartitionDiscontinuities(const PfoToHitTripletsMap &pfoToHitTripletsMap, PartitionVector &partitions) const
{
    PfoVector sortedPfos;
    sortedPfos.reserve(pfoToHitTripletsMap.size());
    for (const auto &[pPfo, _] : pfoToHitTripletsMap)
        sortedPfos.emplace_back(pPfo);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo : sortedPfos)
    {
        const auto &hitTriplets{pfoToHitTripletsMap.at(pPfo)};
        ClusterList pfoClusterList;
        LArPfoHelper::GetTwoDClusterList(pPfo, pfoClusterList);
        const Cluster *pClusterU{nullptr}, *pClusterV{nullptr}, *pClusterW{nullptr};
        for (const Cluster *const pCluster : pfoClusterList)
        {
            switch (LArClusterHelper::GetClusterHitType(pCluster))
            {
                case TPC_VIEW_U:
                {
                    pClusterU = pCluster;
                    break;
                }
                case TPC_VIEW_V:
                {
                    pClusterV = pCluster;
                    break;
                }
                case TPC_VIEW_W:
                {
                    pClusterW = pCluster;
                    break;
                }
                default:
                    break;
            }
        }

        for (const auto &[hitU, hitV, hitW] : hitTriplets)
        {
            const CaloHitList &orderedHitsU{pClusterU ? m_clusterToOrderedHitsMap.at(pClusterU) : CaloHitList()},
                &orderedHitsV{pClusterV ? m_clusterToOrderedHitsMap.at(pClusterV) : CaloHitList()},
                &orderedHitsW{pClusterW ? m_clusterToOrderedHitsMap.at(pClusterW) : CaloHitList()};
            size_t indexU{0}, indexV{0}, indexW{0};
            if (hitU && pClusterU)
            {
                auto it = std::find(orderedHitsU.begin(), orderedHitsU.end(), hitU);
                if (it != orderedHitsU.end())
                    indexU = std::distance(orderedHitsU.begin(), it);
            }
            if (hitV && pClusterV)
            {
                auto it = std::find(orderedHitsV.begin(), orderedHitsV.end(), hitV);
                if (it != orderedHitsV.end())
                    indexV = std::distance(orderedHitsV.begin(), it);
            }
            if (hitW && pClusterW)
            {
                auto it = std::find(orderedHitsW.begin(), orderedHitsW.end(), hitW);
                if (it != orderedHitsW.end())
                    indexW = std::distance(orderedHitsW.begin(), it);
            }

            // Don't split when clusters are too small
            int nSmall{0};
            nSmall += ((orderedHitsU.size() - indexU) < 3) || indexU < 3;
            nSmall += ((orderedHitsV.size() - indexV) < 3) || indexV < 3;
            nSmall += ((orderedHitsW.size() - indexW) < 3) || indexW < 3;
            if (nSmall >= 2)
                continue;

            const float balanceU{this->GetBalance(orderedHitsU, indexU)};
            const float balanceV{this->GetBalance(orderedHitsV, indexV)};
            const float balanceW{this->GetBalance(orderedHitsW, indexW)};
            int nStepUp{0}, nStepDown{0}, nMissing{0};
            nStepUp += balanceU > m_balanceThresholdHigh;
            nStepDown += balanceU && balanceU < m_balanceThresholdLow;
            nMissing += !balanceU;
            nStepUp += balanceV > m_balanceThresholdHigh;
            nStepDown += balanceV && balanceV < m_balanceThresholdLow;
            nMissing += !balanceV;
            nStepUp += balanceW > m_balanceThresholdHigh;
            nStepDown += balanceW && balanceW < m_balanceThresholdLow;
            nMissing += !balanceW;

            if (nStepUp >= 3 || nStepDown >= 3 || (nStepUp == 2 && nMissing == 1) || (nStepDown == 2 && nMissing == 1))
            {
                // We have a consistent discontinuous change in ADC across at least two views
                partitions.emplace_back(Partition(pPfo, std::make_tuple(hitU, hitV, hitW), orderedHitsU, orderedHitsV, orderedHitsW));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::FilterPartitions(PartitionVector &partitions) const
{
    for (auto iter1 = partitions.begin(); iter1 != partitions.end(); ++iter1)
    {
        const Pfo *const pPfo1{iter1->m_pCurrentPfo};
        for (auto iter2 = std::next(iter1); iter2 != partitions.end();)
        {
            const Pfo *const pPfo2{iter2->m_pCurrentPfo};
            if (pPfo2 == pPfo1)
                iter2 = partitions.erase(iter2);
            else
                ++iter2;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::Recluster(const PartitionVector &partitions) const
{
    std::string initialPfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Pfo>(*this, initialPfoListName));
    std::map<HitType, std::string> viewToClusterListNameMap;
    for (const std::string &clusterListName : m_clusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
        for (const Cluster *const pCluster : *pClusterList)
        {
            if (pCluster)
            {
                try
                {
                    viewToClusterListNameMap[LArClusterHelper::GetClusterHitType(pCluster)] = clusterListName;
                }
                catch (const StatusCodeException &)
                {
                    continue;
                }
            }
        }
    }

    const PfoList *pPfoList{nullptr};
    std::string tempPfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, tempPfoListName));

    ProtoPfoVector protoPfos;
    for (const auto &[pPfo, hitTriplet, hitsU, hitsV, hitsW] : partitions)
    {
        ProtoPfo protoPfo;
        protoPfo.m_pfoParameters.m_particleId = pPfo->GetParticleId();
        protoPfo.m_pfoParameters.m_charge = PdgTable::GetParticleCharge(pPfo->GetParticleId());
        protoPfo.m_pfoParameters.m_mass = PdgTable::GetParticleMass(pPfo->GetParticleId());
        protoPfo.m_pfoParameters.m_energy = 0.f;
        protoPfo.m_pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
        protoPfo.m_pOldPfo = pPfo;

        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const CaloHitList &hits{view == TPC_VIEW_U ? hitsU : (view == TPC_VIEW_V ? hitsV : hitsW)};
            const CaloHit *pHit{view == TPC_VIEW_U ? std::get<0>(hitTriplet) : (view == TPC_VIEW_V ? std::get<1>(hitTriplet) : std::get<2>(hitTriplet))};
            CaloHitList cluster1Hits, cluster2Hits;
            // Not all triplets are complete, so we need to check for nullptrs
            if (pHit)
            {
                const auto offset{std::distance(hits.begin(), std::find(hits.begin(), hits.end(), pHit))};
                // Only split the hits if the discontinuity hit is not at the start or end of the cluster
                if (offset != 0 && offset != static_cast<decltype(offset)>(hits.size() - 1))
                    this->PartitionHits(hits, pHit, cluster1Hits, cluster2Hits);
            }
            protoPfo.m_viewToHitsMap[view] = std::move(cluster2Hits);
        }
        protoPfos.emplace_back(std::move(protoPfo));
    }

    // Remove merge hits from their original PFOs, along with associated 3D hits
    for (auto &protoPfo : protoPfos)
    {
        const Pfo *const pPfo{protoPfo.m_pOldPfo};
        CaloHitList pfoHits3D;
        LArPfoHelper::GetCaloHits(pPfo, TPC_3D, pfoHits3D);

        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const CaloHitList &clusterHits{protoPfo.m_viewToHitsMap[view]};
            ClusterList clusterList;
            LArPfoHelper::GetClusters(pPfo, view, clusterList);
            for (const CaloHit *const pCaloHit : clusterHits)
                PandoraContentApi::RemoveFromCluster(*this, clusterList.front(), pCaloHit);

            ClusterList clusterList3D;
            LArPfoHelper::GetClusters(pPfo, TPC_3D, clusterList3D);
            for (const CaloHit *const pCaloHit : pfoHits3D)
            {
                const CaloHit *pParent{static_cast<const CaloHit *>(pCaloHit->GetParentAddress())};
                if (std::find(clusterHits.begin(), clusterHits.end(), pParent) != clusterHits.end())
                    PandoraContentApi::RemoveFromCluster(*this, clusterList3D.front(), pCaloHit);
            }
        }
    }

    // Create the new clusters and associate to the proto PFOs
    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        std::string newClusterListName;
        const ClusterList *pClusterList{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, newClusterListName));

        bool madeClusters{false};
        for (auto &protoPfo : protoPfos)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList = std::move(protoPfo.m_viewToHitsMap[view]);
            // Due to extremal split points, it's possible to have an empty list here
            if (parameters.m_caloHitList.empty())
                continue;

            const Cluster *pNewCluster(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
            protoPfo.m_pfoParameters.m_clusterList.emplace_back(pNewCluster);
            madeClusters = true;
        }
        if (madeClusters)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, viewToClusterListNameMap.at(view)));
    }

    // Create the new PFOs
    for (auto &protoPfo : protoPfos)
    {
        // Due to extremal split points, it's possible to have an empty list here
        if (protoPfo.m_pfoParameters.m_clusterList.empty())
            continue;
        const Pfo *pNewPfo{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, protoPfo.m_pfoParameters, pNewPfo));
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, tempPfoListName, m_pfoListName));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, initialPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::PartitionHits(
    const CaloHitList &caloHits, const CaloHit *const pSplitHit, CaloHitList &cluster1Hits, CaloHitList &cluster2Hits) const
{
    bool passedSplit{false};
    for (const CaloHit *const pCaloHit : caloHits)
    {
        if (pCaloHit == pSplitHit)
        {
            passedSplit = true;
            continue;
        }
        if (!passedSplit)
            cluster1Hits.emplace_back(pCaloHit);
        else
            cluster2Hits.emplace_back(pCaloHit);
    }

    FloatVector adcs1;
    for (const CaloHit *const pCaloHit : cluster1Hits)
        adcs1.emplace_back(pCaloHit->GetInputEnergy());
    const float median1{static_cast<float>(this->GetMedian(adcs1))};
    FloatVector adcs2;
    for (const CaloHit *const pCaloHit : cluster2Hits)
        adcs2.emplace_back(pCaloHit->GetInputEnergy());
    const float median2{static_cast<float>(this->GetMedian(adcs2))};
    const float splitAdc{pSplitHit->GetInputEnergy()};
    const float da1{std::abs(splitAdc - median1)}, da2{std::abs(splitAdc - median2)};
    // Add the split point hit to the cluster most consistent with its ADC
    if (da1 < da2)
        cluster1Hits.emplace_back(pSplitHit);
    else
        cluster2Hits.emplace_back(pSplitHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::OrderHitsRelativeToVertex(
    const CaloHitVector &clusterHits, const LArPointingCluster::Vertex &vertex, CaloHitVector &orderedHits) const
{
    const CartesianVector &position(vertex.GetPosition());
    const CartesianVector &direction(vertex.GetDirection());

    // Get the hit-vertex vector and project it onto the pointing cluster vector
    FloatVector projections;
    for (const CaloHit *const pCaloHit : clusterHits)
    {
        const CartesianVector hitDirection(pCaloHit->GetPositionVector() - position);
        projections.emplace_back(direction.GetDotProduct(hitDirection));
    }

    // Sort the hits according to their projections
    std::vector<size_t> hitIndices(clusterHits.size());
    std::iota(hitIndices.begin(), hitIndices.end(), 0);
    std::sort(hitIndices.begin(), hitIndices.end(), [&projections](size_t i, size_t j) { return projections[i] < projections[j]; });

    for (size_t i = 0; i < hitIndices.size(); ++i)
        orderedHits.emplace_back(clusterHits[hitIndices[i]]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::GetAdcMovingAverage(
    const FloatVector &adcs, FloatVector &movingAverage, FloatVector &movingVariance, const size_t window) const
{
    // Compute backward looking moving average and variance. Peforms rolling update for efficiency
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, window == 0);
    const size_t nHits{adcs.size()};
    movingAverage.resize(nHits);
    movingVariance.resize(nHits);

    FloatVector rollingWindow;
    float S{0.f}, V{0.f};
    size_t count{0};
    // initial window
    for (size_t i = 0; i < std::min(nHits, window); ++i)
    {
        rollingWindow.emplace_back(adcs[i]);
        S += adcs[i];
        V += adcs[i] * adcs[i];
        ++count;

        movingAverage[i] = static_cast<float>(this->GetMedian(rollingWindow));
        movingVariance[i] = count > 1 ? (V - S * S / count) / (count - 1) : 0.0; // sample variance
    }
    // rolling update
    for (size_t i = window; i < nHits; ++i)
    {
        rollingWindow.erase(rollingWindow.begin());
        rollingWindow.emplace_back(adcs[i]);
        const float o{adcs[i - window]};
        const float n{adcs[i]};

        S += n - o;
        V += n * n - o * o;

        movingAverage[i] = static_cast<float>(this->GetMedian(rollingWindow));
        movingVariance[i] = (V - S * S / count) / (count - 1); // sample variance
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
double ShortTrackReclusteringAlgorithm::GetMedian(const std::vector<T> &values) const
{
    std::vector<T> copy(values);
    const size_t mid{copy.size() / 2};

    if (copy.size() % 2 == 0)
    {
        std::nth_element(copy.begin(), copy.begin() + mid, copy.end());
        const double upper{copy[mid]};
        std::nth_element(copy.begin(), copy.begin() + mid - 1, copy.end());
        const double lower{copy[mid - 1]};
        return 0.5 * (lower + upper);
    }
    else
    {
        std::nth_element(copy.begin(), copy.begin() + mid, copy.end());
        return copy[mid];
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::GetStableAdcDiscontinuities(
    const pandora::CaloHitVector &hits, pandora::IntVector &discontinuities, const size_t window) const
{
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, window == 0);
    FloatVector normalizedAdc, movingAverage, movingVariance;
    this->NormalizeAdc(hits, normalizedAdc);
    this->GetAdcMovingAverage(normalizedAdc, movingAverage, movingVariance, window);

    // Look for step changes in ADC
    for (size_t i = window; i < hits.size(); ++i)
    {
        const float current{movingAverage[i]}, previous{movingAverage[i - 1] > 0 ? movingAverage[i - 1] : 1.f},
            previous2{((i >= 2) && (movingAverage[i - 2] > 0)) ? movingAverage[i - 2] : 1.f};
        const float ratio{std::max(current / previous, current / previous2)};

        if (ratio > 2.f)
        {
            // Possible discontinuity to more highly ionising particle, check local variance before step
            int nLow{0}, nHigh{0};
            for (size_t j = i - window; j < i; ++j)
            {
                if (movingVariance[j] < 0.3f) // Arbitrary threshold for now, but probably the right scale
                    ++nLow;
                else
                    ++nHigh;
            }
            if (nHigh >= nLow)
                continue;
            // Check for possible Bragg peak if the step is to more highly ionising particle
            const size_t start{i}, end{std::min(i + 2 * window - 1, hits.size() - 1)};
            if (this->IsBraggPeak(hits, start, end))
                continue;
            discontinuities.emplace_back(i);
        }
        else
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::NormalizeAdc(const pandora::CaloHitVector &hits, pandora::FloatVector &normalizedAdc) const
{
    normalizedAdc.clear();
    FloatVector adcs;
    for (const CaloHit *const pCaloHit : hits)
        adcs.emplace_back(pCaloHit->GetInputEnergy());
    // Find median via nth_element - generally faster than sorting (O(n) vs O(n log n))
    float median{static_cast<float>(this->GetMedian(adcs))};

    // This shouldn't be an issue, but just in case
    if (median < 1.)
        median = 1.;
    for (const CaloHit *const pCaloHit : hits)
        normalizedAdc.emplace_back(pCaloHit->GetInputEnergy() / median);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShortTrackReclusteringAlgorithm::IsBraggPeak(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    const float linearSlopeScore{this->GetLinearSlopeScore(hits, start, end)};
    // Curvature needs at least 5 hits for calculation
    const float curvatureScore{(end - start) >= 4 ? this->GetQuadraticCurvatureScore(hits, start, end) : 0.f};
    const float contrastScore{this->GetContrastScore(hits, start, end)};
    const float monotonicityScore{this->GetMonotonicityScore(hits, start, end)};

    int consensus{0};
    consensus += linearSlopeScore > 0.16f;
    consensus += curvatureScore > 0.f;
    consensus += contrastScore > 1.4f;
    consensus += monotonicityScore > 0.5f;

    return consensus >= 3;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetLinearSlopeScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    FloatVector adcs;
    for (size_t i = start; i <= end; ++i)
        adcs.emplace_back(hits[i]->GetInputEnergy());
    const size_t mid{adcs.size() / 2};
    std::nth_element(adcs.begin(), adcs.begin() + mid, adcs.end());
    const float median{adcs[mid]};

    const float dx{hits[end]->GetPositionVector().GetX() - hits[start]->GetPositionVector().GetX()};
    const float dz{hits[end]->GetPositionVector().GetZ() - hits[start]->GetPositionVector().GetZ()};
    const float dr{std::sqrt(dx * dx + dz * dz)};

    return (dr > 0) ? (hits[end]->GetInputEnergy() - hits[start]->GetInputEnergy()) / (median * dr) : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetQuadraticCurvatureScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    // Use standard three-point non-uniform finite difference for second derivative
    float curvature{0.f};
    for (size_t i = start + 1; i < end; ++i)
    {
        const float a_m{hits[i - 1]->GetInputEnergy()};
        const float a_0{hits[i]->GetInputEnergy()};
        const float a_p{hits[i + 1]->GetInputEnergy()};
        const float dx_m{hits[i]->GetPositionVector().GetX() - hits[i - 1]->GetPositionVector().GetX()};
        const float dz_m{hits[i]->GetPositionVector().GetZ() - hits[i - 1]->GetPositionVector().GetZ()};
        const float dr_m{std::sqrt(dx_m * dx_m + dz_m * dz_m)};
        const float dx_p{hits[i + 1]->GetPositionVector().GetX() - hits[i]->GetPositionVector().GetX()};
        const float dz_p{hits[i + 1]->GetPositionVector().GetZ() - hits[i]->GetPositionVector().GetZ()};
        const float dr_p{std::sqrt(dx_p * dx_p + dz_p * dz_p)};

        curvature += (2.f / (dr_m + dr_p)) * (((a_p - a_0) / dr_p) - ((a_0 - a_m) / dr_m));
    }

    return curvature / (end - start - 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetContrastScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    const size_t mid{(start + end) / 2};
    float muHead{0}, muTail{0};
    for (size_t i = start; i <= mid; ++i)
        muHead += hits[i]->GetInputEnergy();
    muHead /= (mid - start + 1);

    for (size_t i = mid + 1; i <= end; ++i)
        muTail += hits[i]->GetInputEnergy();
    muTail /= (end - mid);

    return (muHead > 0.f) ? muTail / muHead : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetMonotonicityScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    if (start == end)
        return 0.f;
    float monotonicity{0.f};
    for (size_t i = start + 1; i <= end; ++i)
    {
        const float delta{hits[i]->GetInputEnergy() - hits[i - 1]->GetInputEnergy()};
        // We only care about monotonically increasing cases here
        monotonicity += delta > 0.f ? 1.f : 0.f;
    }

    return monotonicity / (end - start);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetBalance(const pandora::CaloHitList &hits, const size_t pivot) const
{
    if (hits.empty() || pivot == 0 || pivot >= hits.size())
        return 0.f;

    FloatVector adcs;
    // Get pre-pivot median
    size_t i{0};
    auto it{hits.begin()};
    for (auto end = hits.end(); it != end && i < pivot; ++it, ++i)
        adcs.emplace_back((*it)->GetInputEnergy());
    size_t mid{adcs.size() / 2};
    std::nth_element(adcs.begin(), adcs.begin() + mid, adcs.end());
    const float medianDenom{!adcs.empty() ? adcs[mid] : 0.f};

    adcs.clear();

    // Get post-pivot median
    for (; it != hits.end(); ++it)
        adcs.emplace_back((*it)->GetInputEnergy());
    mid = adcs.size() / 2;
    std::nth_element(adcs.begin(), adcs.begin() + mid, adcs.end());
    const float medianNum{!adcs.empty() ? adcs[mid] : 0.f};

    return medianDenom > 0 ? medianNum / medianDenom : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShortTrackReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxHitDiscrepancy", m_maxHitDiscrepancy));
    m_maxHitDiscrepancySquared = m_maxHitDiscrepancy * m_maxHitDiscrepancy;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BalanceThresholdLow", m_balanceThresholdLow));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BalanceThresholdHigh", m_balanceThresholdHigh));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BraggLinearSlopeThreshold", m_braggLinearSlopeThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BraggCurvatureThreshold", m_braggCurvatureThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BraggContrastThreshold", m_braggContrastThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BraggMonotonicityThreshold", m_braggMonotonicityThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
