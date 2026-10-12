/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLThreeViewMergeAndCreateShowersTool.cc
 *
 *  @brief  Implementation of the merge and create showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLThreeViewMergeAndCreateShowersTool.h"

using namespace pandora;

namespace lar_dl_content
{

DLThreeViewMergeAndCreateShowersTool::DLThreeViewMergeAndCreateShowersTool() :
    m_matchThreshold(0.5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLThreeViewMergeAndCreateShowersTool::Run(DLMultiViewMatchingAlgorithm *const pAlgorithm,
    const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    DLMultiViewMatchingAlgorithm::ClusterGroupVector clusterGroupVector;
    pAlgorithm->GetConnectedGroups(clusterGroupVector);  

    bool madeParticles(false);
    for (const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup : clusterGroupVector)
    {
        int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());                
            
        if ((nU * nV * nW) <= 1)
            continue;

        const bool thisMadeParticles(this->CreateShower(pAlgorithm, clusterGroup, globalSimMatrix));
        madeParticles = madeParticles || thisMadeParticles;
    }

    return madeParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLThreeViewMergeAndCreateShowersTool::CreateShower(DLMultiViewMatchingAlgorithm *const pAlgorithm,
    const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix)
{
    const Cluster *pSeedU(nullptr), *pSeedV(nullptr), *pSeedW(nullptr);
    if (!this->FindSeed(clusterGroup, globalSimMatrix, pSeedU, pSeedV, pSeedW)) { return false; }

    this->MergeClusters(pAlgorithm, clusterGroup, globalSimMatrix, pSeedU, pSeedV, pSeedW);
    pAlgorithm->CreatePfo({pSeedU, pSeedV, pSeedW});
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLThreeViewMergeAndCreateShowersTool::FindSeed(const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
    const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const Cluster *&pSeedU, const Cluster *&pSeedV, const Cluster *&pSeedW)
{
    bool found(false);
    int highestHits(0);
    
    for (const Cluster *const pClusterU : clusterGroup.m_clustersU)
    {
        const auto iterU(globalSimMatrix.find(pClusterU));
        if (iterU == globalSimMatrix.end()) { continue; }
        
        for (const Cluster *const pClusterV : clusterGroup.m_clustersV)
        {
            const auto iterUV(iterU->second.find(pClusterV));
            if (iterUV == iterU->second.end()) { continue; }
            if (iterUV->second < m_matchThreshold) { continue; }
            const auto iterV(globalSimMatrix.find(pClusterV));
            if (iterV == globalSimMatrix.end()) { continue; }
            
            for (const Cluster *const pClusterW : clusterGroup.m_clustersW)
            {
                const auto iterUW(iterU->second.find(pClusterW));
                if (iterUW == iterU->second.end()) { continue; }
                if (iterUW->second < m_matchThreshold) { continue; }
                const auto iterVW(iterV->second.find(pClusterW));
                if (iterVW == iterV->second.end()) { continue; }
                if (iterVW->second < m_matchThreshold) { continue; }                
        
                const int nHits(pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits());

                if (nHits > highestHits)
                {
                    found = true;
                    highestHits = nHits;
                    pSeedU = pClusterU;
                    pSeedV = pClusterV;
                    pSeedW = pClusterW;
                }
            }
        } 
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeViewMergeAndCreateShowersTool::MergeClusters(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
    const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const Cluster *const pSeedU, const Cluster *const pSeedV, const Cluster *const pSeedW)
{
    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const Cluster *const pSeed(hitType == TPC_VIEW_U ? pSeedU : hitType == TPC_VIEW_V ? pSeedV : pSeedW);
        const ClusterList &clusterList(hitType == TPC_VIEW_U ? clusterGroup.m_clustersU : hitType == TPC_VIEW_V ? clusterGroup.m_clustersV : clusterGroup.m_clustersW);
        const std::string clusterListName(hitType == TPC_VIEW_U ? m_clusterListNameU : hitType == TPC_VIEW_V ? m_clusterListNameV : m_clusterListNameW);
        
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

StatusCode DLThreeViewMergeAndCreateShowersTool::ReadSettings(const TiXmlHandle xmlHandle)
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

