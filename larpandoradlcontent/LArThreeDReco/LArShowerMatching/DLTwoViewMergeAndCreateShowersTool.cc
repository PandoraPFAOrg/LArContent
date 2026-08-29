/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLTwoViewMergeAndCreateShowersTool.cc
 *
 *  @brief  Implementation of the merge and create showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLTwoViewMergeAndCreateShowersTool.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLTwoViewMergeAndCreateShowersTool::DLTwoViewMergeAndCreateShowersTool() :
    m_matchThreshold(0.5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoViewMergeAndCreateShowersTool::Run(DLMultiViewMatchingAlgorithm *const pAlgorithm,
    const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    DLMultiViewMatchingAlgorithm::ClusterGroupVector clusterGroupVector;
    pAlgorithm->GetConnectedGroups(clusterGroupVector);  

    bool madeParticles(false);
    for (const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup : clusterGroupVector)
    {
        int nZeroClusterViews(0), nOneClusterViews(0);

        for (const ClusterList &clusterList : {clusterGroup.m_clustersU, clusterGroup.m_clustersV, clusterGroup.m_clustersW})
        {
            if (clusterList.empty())
                ++nZeroClusterViews;
            else if (clusterList.size() == 1)
                ++nOneClusterViews;
        }

        if ((nZeroClusterViews != 1) || (nOneClusterViews == 2))
            continue;

        const HitType ignoreHitType(clusterGroup.m_clustersU.size() == 0 ? TPC_VIEW_U : clusterGroup.m_clustersV.size() == 0 ? TPC_VIEW_V : TPC_VIEW_W);
        const HitType hitType1(ignoreHitType == TPC_VIEW_U ? TPC_VIEW_V : ignoreHitType == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
        const HitType hitType2(ignoreHitType == TPC_VIEW_U ? TPC_VIEW_W : ignoreHitType == TPC_VIEW_V ? TPC_VIEW_U : TPC_VIEW_V);        

        const bool thisMadeParticles(this->CreateShower(pAlgorithm, clusterGroup, globalSimMatrix, hitType1, hitType2));
        madeParticles = madeParticles ? madeParticles : thisMadeParticles;
    }

    return madeParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoViewMergeAndCreateShowersTool::CreateShower(DLMultiViewMatchingAlgorithm *const pAlgorithm,
    const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix,
    const HitType hitType1, const HitType hitType2)
{
    const Cluster *pSeed1(nullptr), *pSeed2(nullptr);
    if (!this->FindSeed(clusterGroup, globalSimMatrix, hitType1, hitType2, pSeed1, pSeed2)) { return false; }

    this->MergeClusters(pAlgorithm, clusterGroup, globalSimMatrix, pSeed1, pSeed2);
    pAlgorithm->CreatePfo({pSeed1, pSeed2});

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoViewMergeAndCreateShowersTool::FindSeed(const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
    const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const HitType hitType1, const HitType hitType2,
    const Cluster *&pSeed1, const Cluster *&pSeed2)
{
    const ClusterList clusterList1(hitType1 == TPC_VIEW_U ? clusterGroup.m_clustersU : hitType1 == TPC_VIEW_V ? clusterGroup.m_clustersV : clusterGroup.m_clustersW);
    const ClusterList clusterList2(hitType2 == TPC_VIEW_U ? clusterGroup.m_clustersU : hitType2 == TPC_VIEW_V ? clusterGroup.m_clustersV : clusterGroup.m_clustersW);
    bool found(false);
    int highestHits(0);    
    
    for (const Cluster *const pCluster1 : clusterList1)
    {
        const auto iter1(globalSimMatrix.find(pCluster1));
        if (iter1 == globalSimMatrix.end()) { continue; }
        
        for (const Cluster *const pCluster2 : clusterList2)
        {
            const auto iter12(iter1->second.find(pCluster2));
            if (iter12 == iter1->second.end()) { continue; }
            if (iter12->second < m_matchThreshold) { continue; }
        
            const int nHits(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits());

            if (nHits > highestHits)
            {
                found = true;
                highestHits = nHits;
                pSeed1 = pCluster1;
                pSeed2 = pCluster2;
            }
        } 
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLTwoViewMergeAndCreateShowersTool::MergeClusters(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
    const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const Cluster *const pSeed1, const Cluster *const pSeed2)
{
    for (const Cluster *const pSeed : {pSeed1, pSeed2})
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pSeed));            
        const ClusterList &clusterList(hitType == TPC_VIEW_U ? clusterGroup.m_clustersU : hitType == TPC_VIEW_V ? clusterGroup.m_clustersV : clusterGroup.m_clustersW);
        const std::string clusterListName(hitType == TPC_VIEW_U ? m_clusterListNameU : hitType == TPC_VIEW_V ? m_clusterListNameV :m_clusterListNameW);
        
        for (const Cluster *const pClusterToMerge : clusterList)
        {
            if (pSeed == pClusterToMerge)
                continue;
            
            const auto iter1(globalSimMatrix.find(pSeed));
            if (iter1 == globalSimMatrix.end()) { continue; }
            const auto iter2(iter1->second.find(pClusterToMerge));
            if (iter2 == iter1->second.end()) { continue; }
            if (iter2->second < m_matchThreshold) { continue; }

            pAlgorithm->DeleteCluster(pClusterToMerge);            

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(
                *pAlgorithm, pSeed, pClusterToMerge, clusterListName, clusterListName));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoViewMergeAndCreateShowersTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    if (m_clusterListNameU.empty())
        m_clusterListNameU = "ClustersU";
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    if (m_clusterListNameV.empty())
        m_clusterListNameV = "ClustersV";
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));
    if (m_clusterListNameW.empty())
        m_clusterListNameW = "ClustersW";
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchThreshold", m_matchThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content

