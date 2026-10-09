/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/MultiViewMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the multi view matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Geometry/LArReadoutChannel.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLShowerHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLMultiViewMatchingAlgorithm::DLMultiViewMatchingAlgorithm() :
    m_trainingMode(false),
    m_trainingFileName("ShowerMatchingTraining.root"),
    m_trainingTreeName("trainingTree"),
    m_minNClusterHits(5),    
    m_nMaxRepeats(2),
    m_minWireOverlapFraction(0.8f),
    m_hitFeatureDim(14),
    m_polarRScaleFactor{1.f},
    m_cartesianXScaleFactor{1.f},
    m_cartesianZScaleFactor{1.f},
    m_matchThreshold(0.5f)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DLMultiViewMatchingAlgorithm::~DLMultiViewMatchingAlgorithm()
{
    if (m_trainingMode)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_trainingFileName, "UPDATE"));

        for ([[maybe_unused]] const HitType &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_view_data", "view", static_cast<int>(view)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                m_trainingTreeName + "_view_data", "pitch", lar_content::LArGeometryHelper::GetWirePitch(this->GetPandora(), view)));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName + "_view_data"));
        }
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName + "_view_data", m_trainingFileName, "UPDATE"));
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode DLMultiViewMatchingAlgorithm::Run()
{
    const ClusterList *pClusterListU(nullptr), *pClusterListV(nullptr), *pClusterListW(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetList(m_clusterListNameU, pClusterListU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetList(m_clusterListNameV, pClusterListV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetList(m_clusterListNameW, pClusterListW));
    
    ClusterExtentMap clusterExtentMap;
    this->PrepareClusters(pClusterListU, pClusterListV, pClusterListW, clusterExtentMap);

    this->FillNavigationMaps(clusterExtentMap);

    ClusterGroupVector clusterGroupVector;
    this->GetConnectedGroups(clusterGroupVector);

    if (m_trainingMode)
    {
        this->PrepareTrainingSample(clusterGroupVector);
    }
    else
    {
        SimilarityMatrix globalSimMatrix;
        this->FillGlobalSimMatrix(clusterGroupVector, globalSimMatrix);
        
        this->UpdateNavigationMaps(globalSimMatrix);
    
        bool repeat(true);        
        unsigned int repeatCounter(0);
        while (repeat && (repeatCounter < m_nMaxRepeats))
        {
            repeat = false;
            ++repeatCounter;            
            for (const auto &matchingTool : m_matchingToolVector)
            {
                const bool particlesMade(matchingTool->Run(this, globalSimMatrix));
                repeat = repeat || particlesMade;
            }
        }
    }

    this->CleanUp();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------  

template <typename T>
StatusCode DLMultiViewMatchingAlgorithm::GetList(const std::string listName, const T *&pList)
{
    pList = nullptr;

    if (PandoraContentApi::GetList(*this, listName, pList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    if ((!pList) || pList->empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

// //------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::PrepareClusters(const ClusterList *const pClusterListU, const ClusterList *const pClusterListV,
    const ClusterList *const pClusterListW, ClusterExtentMap &clusterExtentMap)
{
    for (const ClusterList *const pClusterList : {pClusterListU, pClusterListV, pClusterListW})
    {
        if (pClusterList->empty())
            continue;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterList->front()));
        ClusterList &filteredClusters((hitType == TPC_VIEW_U) ? m_filteredU : (hitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);
        
        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!pCluster->IsAvailable())
                continue;
            
            if (pCluster->GetNCaloHits() < m_minNClusterHits)
                continue;

            filteredClusters.emplace_back(pCluster);

            // Get drift-extent
            CartesianVector minCoord(0.f, 0.f, 0.f), maxCoord(0.f, 0.f, 0.f);
            LArClusterHelper::GetClusterBoundingBox(pCluster, minCoord, maxCoord, true);

            clusterExtentMap[pCluster] = std::pair(minCoord.GetX(), maxCoord.GetX());           
        }
    }    
}

// //------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::FillNavigationMaps(const ClusterExtentMap &clusterExtentMap)
{
    std::vector<std::pair<HitType, HitType>> viewCombinations({{TPC_VIEW_U, TPC_VIEW_V}, {TPC_VIEW_U, TPC_VIEW_W}, {TPC_VIEW_V, TPC_VIEW_W}});
    for (const auto& [hitType, projHitType] : viewCombinations)
    {
        const ClusterList &clusters((hitType == TPC_VIEW_U) ? m_filteredU : (hitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);
        const ClusterList &projClusters((projHitType == TPC_VIEW_U) ? m_filteredU : (projHitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);        
        NavigationMap &navigationMap((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);
        NavigationMap &projNavigationMap((projHitType == TPC_VIEW_U) ? m_navigationU : (projHitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);

        for (const Cluster *const pCluster : clusters)
        {
            for (const Cluster *const pProjCluster : projClusters)
            {
                if (this->DoClustersOverlap(pCluster, pProjCluster, clusterExtentMap))
                {
                    navigationMap[pCluster].emplace_back(pProjCluster);
                    projNavigationMap[pProjCluster].emplace_back(pCluster);
                }
            }
        }
    }
}

// //------------------------------------------------------------------------------------------------------------------------------------------  
    
bool DLMultiViewMatchingAlgorithm::DoClustersOverlap(const Cluster *const pCluster1, const Cluster *const pCluster2, const ClusterExtentMap &clusterExtentMap)
{
    const std::pair<float, float> &driftExtent1(clusterExtentMap.at(pCluster1));
    const std::pair<float, float> &driftExtent2(clusterExtentMap.at(pCluster2));
    const float minX(std::max(driftExtent1.first, driftExtent2.first));
    const float maxX(std::min(driftExtent1.second, driftExtent2.second));            
    const float driftSpan(maxX - minX);
            
    if (driftSpan < std::numeric_limits<float>::epsilon())
        return false;

    return this->DoClustersOverlapInWire(pCluster1, pCluster2, minX, maxX);
}

// //------------------------------------------------------------------------------------------------------------------------------------------      

bool DLMultiViewMatchingAlgorithm::DoClustersOverlapInWire(const Cluster *const pCluster1, const Cluster *const pCluster2, const float minX, const float maxX)
{
    CaloHitList caloHitList1, caloHitList2, filteredCaloHitList1, filteredCaloHitList2;
    LArClusterHelper::GetAllHits(pCluster1, caloHitList1);
    LArClusterHelper::GetAllHits(pCluster2, caloHitList2);
    
    for (unsigned int iList : {0, 1})
    {
        CaloHitList &caloHitList((iList == 0) ? caloHitList1 : caloHitList2);
        CaloHitList &filteredCaloHitList((iList == 0) ? filteredCaloHitList1 : filteredCaloHitList2);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if ((pCaloHit->GetPositionVector().GetX() < minX) || (pCaloHit->GetPositionVector().GetX() > maxX))
                continue;

            filteredCaloHitList.emplace_back(pCaloHit);
        }
    }

    // Search for overlap
    int overlapCount(0);
    const int nSamplingPoints1(filteredCaloHitList1.size());        
    const int nSamplingPoints2(filteredCaloHitList2.size());    
    for (const CaloHit *const pCaloHit1 : filteredCaloHitList1)
    {
        const LArCaloHit *const pLArHit1(dynamic_cast<const LArCaloHit *>(pCaloHit1));

        if (!pLArHit1)
            continue;
        
        for (const CaloHit *const pCaloHit2 : filteredCaloHitList2)
        {
            const LArCaloHit *const pLArHit2(dynamic_cast<const LArCaloHit *>(pCaloHit2));

            if (!pLArHit2)
                continue;
            
            // Readout volume should be the same
            if (LArGeometryHelper::GetReadoutVolume(this->GetPandora(), *pLArHit1).GetId() !=
                LArGeometryHelper::GetReadoutVolume(this->GetPandora(), *pLArHit2).GetId())
            {
                continue;
            }

            // pLAr1 get overlap in pLAr2 hit type...
            const LArReadoutChannel &readoutChannel1(LArGeometryHelper::GetReadoutChannel(this->GetPandora(), *pLArHit1));
            const LArReadoutChannel &readoutChannel2(LArGeometryHelper::GetReadoutChannel(this->GetPandora(), *pLArHit2));         
            const LArReadoutChannel::ChannelInterval &channelIntervalIn2(readoutChannel1.GetChannelInterval(pLArHit2->GetHitType()));
            
            if ((readoutChannel2.GetId() < channelIntervalIn2.first) || (readoutChannel2.GetId() > channelIntervalIn2.second))
                continue;

            overlapCount++;
            break;
        }
    }

    const float overlapFraction1(nSamplingPoints1 == 0 ? 0.f : float(overlapCount) / float(nSamplingPoints1));
    const float overlapFraction2(nSamplingPoints2 == 0 ? 0.f : float(overlapCount) / float(nSamplingPoints2));
    std::cout << "overlapFraction1: " << overlapFraction1 << std::endl;
    std::cout << "overlapFraction2: " << overlapFraction2 << std::endl;    
    
    return ((overlapFraction1 > m_minWireOverlapFraction) && (overlapFraction2 > m_minWireOverlapFraction));
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::GetConnectedGroups(ClusterGroupVector &clusterGroupVector)
{
    ClusterList usedU, usedV, usedW;
    
    for (const ClusterList &clusterList : {m_filteredU, m_filteredV, m_filteredW})
    {
        for (const Cluster *const pCluster : clusterList)
        {
            ClusterGroup clusterGroup;
            this->GetConnectedGroup(pCluster, clusterGroup, usedU, usedV, usedW);

            int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());
            
            if ((nU + nV + nW) <= 1)
                continue;

            clusterGroupVector.push_back(clusterGroup);  
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::GetConnectedGroup(const Cluster *const pCluster, ClusterGroup &clusterGroup, ClusterList &usedU, ClusterList &usedV, ClusterList &usedW)
{
    if (!pCluster->IsAvailable())
        return;
    
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    ClusterList &used((hitType == TPC_VIEW_U) ? usedU : (hitType == TPC_VIEW_V) ? usedV : usedW);
        
    if (std::find(used.begin(), used.end(), pCluster) != used.end())
        return;

    // Add cluster
    ClusterList &clusterGroupList((hitType == TPC_VIEW_U) ? clusterGroup.m_clustersU : (hitType == TPC_VIEW_V) ? clusterGroup.m_clustersV : clusterGroup.m_clustersW);
    clusterGroupList.emplace_back(pCluster);
    used.emplace_back(pCluster);
    
    // Now find its connections and add those
    const NavigationMap &navigation((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);
    const auto navIter(navigation.find(pCluster));

    if (navIter == navigation.end())
        return;
    
    for (const Cluster *const pNavCluster : navIter->second)
        this->GetConnectedGroup(pNavCluster, clusterGroup, usedU, usedV, usedW);
}

//------------------------------------------------------------------------------------------------------------------------------------------     
    
void DLMultiViewMatchingAlgorithm::FillGlobalSimMatrix(const ClusterGroupVector &clusterGroupVector, SimilarityMatrix &globalSimMatrix)
{
    // Get nu vertex projections
    const VertexList *pVertexList(nullptr);
    if (this->GetList(m_nuVertexListName, pVertexList) != STATUS_CODE_SUCCESS)
        return;
    const Vertex *const pVertex(pVertexList->front());
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, pVertex->GetVertexType() != VERTEX_3D);
    const std::map<HitType, CartesianVector> viewToVtxPos(
           {{TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U)},
            {TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V)},
            {TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W)}});

    std::set<float> detXGaps;
    LArGeometryHelper::GetDetectorXGaps(this->GetPandora(), detXGaps);

    // Add scores to the global similarity matrix
    for (const ClusterGroup &clusterGroup : clusterGroupVector)
    {
        this->PredictClusterSimilarityMatrix(clusterGroup.m_clustersU, clusterGroup.m_clustersV, clusterGroup.m_clustersW,
            viewToVtxPos, detXGaps, globalSimMatrix);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMultiViewMatchingAlgorithm::PredictClusterSimilarityMatrix(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW, 
    const std::map<HitType, CartesianVector> &viewToVtxPos, const std::set<float> &detXGaps, SimilarityMatrix &clusterSimMat)
{
    ClusterVector clusterVecU(clusterListU.begin(), clusterListU.end());
    std::sort(clusterVecU.begin(), clusterVecU.end(), LArClusterHelper::SortByNHits);
    ClusterVector clusterVecV(clusterListV.begin(), clusterListV.end());
    std::sort(clusterVecV.begin(), clusterVecV.end(), LArClusterHelper::SortByNHits);
    ClusterVector clusterVecW(clusterListW.begin(), clusterListW.end());
    std::sort(clusterVecW.begin(), clusterVecW.end(), LArClusterHelper::SortByNHits);
    int nClusters(clusterVecU.size() + clusterVecV.size() + clusterVecW.size());
    ClusterVector clusterVec;
    clusterVec.insert(clusterVec.end(), clusterVecU.begin(), clusterVecU.end());
    clusterVec.insert(clusterVec.end(), clusterVecV.begin(), clusterVecV.end());
    clusterVec.insert(clusterVec.end(), clusterVecW.begin(), clusterVecW.end());

    // Get our similarity matrix
    torch::InferenceMode guard;
    std::vector<torch::Tensor> tensorEncodedClusters;

    for (const Cluster *const pCluster : clusterVec)
    {
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);

        const HitType view(LArClusterHelper::GetClusterHitType(pCluster));

        std::vector<LArDLShowerHelper::HitFeatures> clusterFeatures;
        for (const CaloHit *const pCaloHit : clusterHits)
        {
            LArDLShowerHelper::HitFeatures hitFeatures;
            LArDLShowerHelper::CalculateHitFeatures(pCaloHit, detXGaps, viewToVtxPos.at(view), hitFeatures);
            clusterFeatures.emplace_back(hitFeatures);
        }

        torch::Tensor tensorCluster;
        this->MakeClusterTensor(clusterFeatures, view, nClusters, tensorCluster);
        torch::Tensor tensorEncodedCluster{m_modelEncoder.forward({tensorCluster}).toTensor()};
        tensorEncodedClusters.emplace_back(tensorEncodedCluster);
    }

    torch::Tensor tensorEncodedEvent{torch::cat(tensorEncodedClusters, 1)};
    tensorEncodedClusters.clear(); // Free memory
    torch::Tensor tensorAttnEvent{m_modelAttn.forward({tensorEncodedEvent}).toTensor()};
    tensorEncodedEvent = torch::Tensor(); // Free memory
    torch::Tensor tensorSimMat{m_modelSim.forward({tensorAttnEvent}).toTensor()};
    tensorAttnEvent = torch::Tensor();

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PopulateClusterSimilarityMatrix(tensorSimMat, clusterVec, clusterSimMat));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::MakeClusterTensor(const std::vector<LArDLShowerHelper::HitFeatures> &clusterFeatures, const HitType view,
    const int nClusters, torch::Tensor &tensorCluster) const
{
    // Global features
    const int nHits{static_cast<int>(clusterFeatures.size())};
    const float nHitsFeat{std::log(static_cast<float>(nHits))};
    const float nClustersFeat{std::log(static_cast<float>(nClusters))};
    const float zWidthFeat{LArGeometryHelper::GetWirePitch(this->GetPandora(), view) * m_cartesianZScaleFactor};

    // Fill tensor
    tensorCluster = torch::zeros({1, nHits, m_hitFeatureDim});
    auto accessor = tensorCluster.accessor<float, 3>();
    for (int i = 0; i < nHits; i++)
    {
        const LArDLShowerHelper::HitFeatures hitFeatures{clusterFeatures.at(i)};
        accessor[0][i][0] = hitFeatures.m_rRel * m_polarRScaleFactor;
        accessor[0][i][1] = hitFeatures.m_cosThetaRel;
        accessor[0][i][2] = hitFeatures.m_sinThetaRel;
        accessor[0][i][3] = hitFeatures.m_xRel * m_cartesianXScaleFactor;
        accessor[0][i][4] = hitFeatures.m_zRel * m_cartesianZScaleFactor;
        accessor[0][i][5] = hitFeatures.m_xWidth * m_cartesianXScaleFactor;
        accessor[0][i][6] = zWidthFeat;
        accessor[0][i][7] = hitFeatures.m_distToXGap * m_cartesianXScaleFactor;
        accessor[0][i][8] = std::log(hitFeatures.m_energy);
        accessor[0][i][9] = view == TPC_VIEW_U ? 1.f : 0.f;
        accessor[0][i][10] = view == TPC_VIEW_V ? 1.f : 0.f;
        accessor[0][i][11] = view == TPC_VIEW_W ? 1.f : 0.f;
        accessor[0][i][12] = nHitsFeat;
        accessor[0][i][13] = nClustersFeat;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMultiViewMatchingAlgorithm::PopulateClusterSimilarityMatrix(const torch::Tensor &tensorSimMat, const ClusterVector &clusterVector,
    SimilarityMatrix &clusterSimMat) const
{
    const int64_t nClusters{static_cast<int64_t>(clusterVector.size())};
    PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED, tensorSimMat.dim() != 3 || nClusters != tensorSimMat.size(-1) || nClusters != tensorSimMat.size(-2));
    
    auto accessor = tensorSimMat.accessor<float, 3>();
    auto iterI{clusterVector.begin()};
    for (int i = 0; i < nClusters; i++, iterI++)
    {
        auto iterJ{clusterVector.begin()};
        for (int j = 0; j < nClusters; j++, iterJ++)
            clusterSimMat[*iterI][*iterJ] = accessor[0][i][j];
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::UpdateNavigationMaps(const SimilarityMatrix &globalSimMatrix)
{
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        NavigationMap &navigationMap((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);

        for (auto &[pCluster, clusterList] : navigationMap)
        {
            for (auto iter = clusterList.begin(); iter != clusterList.end(); )
            {
                auto simIter1(globalSimMatrix.find(pCluster));
                if (simIter1 == globalSimMatrix.end()) {++iter; continue;}
                auto simIter2(simIter1->second.find(*iter));
                if (simIter2 == simIter1->second.end()) {++iter; continue;}                

                if (simIter2->second < m_matchThreshold)
                    iter = clusterList.erase(iter);
                else
                    ++iter;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::CreatePfo(const ClusterList &clusters)
{
    const PfoList *pPfoList(nullptr);
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = E_MINUS;
    pfoParameters.m_charge = PdgTable::GetParticleCharge(E_MINUS);
    pfoParameters.m_mass = PdgTable::GetParticleMass(E_MINUS);
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusters;
    
    const ParticleFlowObject *pPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }
}

// //------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::DeleteCluster(const Cluster *const pClusterToRemove)
{
    // Remove from filtered lists
    const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterToRemove));
    ClusterList &filteredClusters((hitType == TPC_VIEW_U) ? m_filteredU : (hitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);

    auto iter(std::find(filteredClusters.begin(), filteredClusters.end(), pClusterToRemove));
    if (iter != filteredClusters.end()) { filteredClusters.erase(iter); }

    // Remove from navigation maps
    NavigationMap &navigation((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);
    auto iterNav(navigation.find(pClusterToRemove));

    if (iterNav != navigation.end())
    {
        for (const Cluster *const pAssocCluster : iterNav->second)
        {
            const HitType assocHitType(LArClusterHelper::GetClusterHitType(pAssocCluster));
            NavigationMap &assocNavigation((assocHitType == TPC_VIEW_U) ? m_navigationU : (assocHitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);

            auto iterAssoc(assocNavigation.find(pAssocCluster));
            if (iterAssoc == assocNavigation.end()) { continue; }

            auto iterRemove(std::find(iterAssoc->second.begin(), iterAssoc->second.end(), pClusterToRemove));
            if (iterRemove != iterAssoc->second.end()) { iterAssoc->second.erase(iterRemove); }
        }

        navigation.erase(iterNav);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::PrepareTrainingSample(const ClusterGroupVector &clusterGroupVector)
{
    // Get nu vertex projections
    const VertexList *pVertexList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_nuVertexListName, pVertexList));
    const Vertex *const pVertex{pVertexList->front()};
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, pVertex->GetVertexType() != VERTEX_3D);
    const std::map<HitType, CartesianVector> viewToVtxPos(
           {{TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U)},
            {TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V)},
            {TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W)}});

    std::set<float> detXGaps;
    LArGeometryHelper::GetDetectorXGaps(this->GetPandora(), detXGaps);

    // Fill tree for each isolated group
    for (const ClusterGroup &clusterGroup : clusterGroupVector)
    {
        int currClusterID{0}, currMCID{0};
        std::map<const MCParticle *const, int> mcToID = {{nullptr, -1}}, mcToPDG = {{nullptr, 0}};
        std::vector<int> clusterView, clusterID;
        std::vector<int> mcID = {-1}, mcPDG = {0};
        std::vector<int> hitClusterID;
        std::vector<float> hitXRelPos, hitZRelPos;
        std::vector<float> hitRRelPos, hitSinThetaRelPos, hitCosThetaRelPos; // Polar coordinates
        std::vector<float> hitXWidth;
        std::vector<float> hitDistToXGap;
        std::vector<float> hitEnergy;
        std::vector<int> hitMCID;

        for (const ClusterList &clusterList : {clusterGroup.m_clustersU, clusterGroup.m_clustersV, clusterGroup.m_clustersW})
        {
            for (const Cluster *const pCluster : clusterList)
            {
                const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
                bool filled(false);

                CaloHitList clusterCaloHits;
                LArClusterHelper::GetAllHits(pCluster, clusterCaloHits);
                for (const CaloHit *const pCaloHit : clusterCaloHits)
                {
                    LArDLShowerHelper::HitFeatures hitFeatures;
                    LArDLShowerHelper::CalculateHitFeatures(pCaloHit, detXGaps, viewToVtxPos.at(view), hitFeatures);

                    try
                    {
                        const MCParticle *const pMainMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                        if (mcToID.find(pMainMC) == mcToID.end())
                        {
                            mcToID.insert({pMainMC, currMCID++});
                            mcToPDG.insert({pMainMC, pMainMC->GetParticleId()});
                            mcID.emplace_back(mcToID.at(pMainMC));
                            mcPDG.emplace_back(mcToPDG.at(pMainMC));
                        }

                        filled = true;
                        hitClusterID.emplace_back(currClusterID);
                        hitXRelPos.emplace_back(hitFeatures.m_xRel);
                        hitZRelPos.emplace_back(hitFeatures.m_zRel);
                        hitRRelPos.emplace_back(hitFeatures.m_rRel);
                        hitCosThetaRelPos.emplace_back(hitFeatures.m_cosThetaRel);
                        hitSinThetaRelPos.emplace_back(hitFeatures.m_sinThetaRel);
                        hitXWidth.emplace_back(hitFeatures.m_xWidth);
                        hitDistToXGap.emplace_back(hitFeatures.m_distToXGap);
                        hitEnergy.emplace_back(hitFeatures.m_energy);
                        hitMCID.emplace_back(mcToID.at(pMainMC));
                    }
                    catch(...) { continue; }
                }

                if (filled)
                {
                    clusterView.emplace_back(static_cast<int>(view));
                    clusterID.emplace_back(currClusterID);

                    currClusterID++;
                }
            }
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_id", &clusterID));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_view", &clusterView));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_id", &mcID));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_pdg", &mcPDG));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_cluster_id", &hitClusterID));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_rel_pos", &hitXRelPos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_z_rel_pos", &hitZRelPos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_r_rel_pos", &hitRRelPos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_ctheta_rel_pos", &hitCosThetaRelPos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_stheta_rel_pos", &hitSinThetaRelPos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_width", &hitXWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_gap_dist", &hitDistToXGap));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_energy", &hitEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_mc_id", &hitMCID));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::CleanUp()
{
    m_filteredU.clear();
    m_filteredV.clear();
    m_filteredW.clear();
    m_navigationU.clear();
    m_navigationV.clear();
    m_navigationW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode DLMultiViewMatchingAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "MatchingTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DLShowerMatchingTool *const pMatchingTool(dynamic_cast<DLShowerMatchingTool *>(*iter));

        if (!pMatchingTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_matchingToolVector.push_back(pMatchingTool);
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexListName", m_nuVertexListName));
    if (m_nuVertexListName.empty())
        m_nuVertexListName = "NeutrinoVertices3D";
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
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    if (m_outputPfoListName.empty())
        m_outputPfoListName = "ShowerParticles3D";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
    }
    else
    {
        std::string modelEncoderName, modelAttnName, modelSimName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelEncoderFileName", modelEncoderName));
        modelEncoderName = LArFileHelper::FindFileInPath(modelEncoderName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelEncoderName, m_modelEncoder);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelAttnFileName", modelAttnName));
        modelAttnName = LArFileHelper::FindFileInPath(modelAttnName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelAttnName, m_modelAttn);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelSimFileName", modelSimName));
        modelSimName = LArFileHelper::FindFileInPath(modelSimName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelSimName, m_modelSim);                
        
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "NMaxRepeats", m_nMaxRepeats));
               
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "HitFeatureDim", m_hitFeatureDim));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinNClusterHits", m_minNClusterHits));

     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinWireOverlapFraction", m_minWireOverlapFraction));
    
    const LArGeometryHelper::DetectorBoundaries detBounds{LArGeometryHelper::GetDetectorBoundaries(this->GetPandora())};
    double xLow{static_cast<double>(detBounds.m_xBoundaries.first)}, xHigh{static_cast<double>(detBounds.m_xBoundaries.second)};
    double yLow{static_cast<double>(detBounds.m_yBoundaries.first)}, yHigh{static_cast<double>(detBounds.m_yBoundaries.second)};
    double zLow{static_cast<double>(detBounds.m_zBoundaries.first)}, zHigh{static_cast<double>(detBounds.m_zBoundaries.second)};
    m_polarRScaleFactor = static_cast<float>(1. / std::sqrt(std::pow(xHigh - xLow, 2.) + std::pow(yHigh - yLow, 2.) + std::pow(zHigh - zLow, 2.)));
    m_cartesianXScaleFactor = static_cast<float>(1. / (xHigh - xLow));
    m_cartesianZScaleFactor = static_cast<float>(1. / (zHigh - zLow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MatchThreshold", m_matchThreshold));
    
    return STATUS_CODE_SUCCESS;
}   

// //------------------------------------------------------------------------------------------------------------------------------------------

template StatusCode DLMultiViewMatchingAlgorithm::GetList(const std::string, const ClusterList *&);
template StatusCode DLMultiViewMatchingAlgorithm::GetList(const std::string, const VertexList *&);
template StatusCode DLMultiViewMatchingAlgorithm::GetList(const std::string, const MCParticleList *&);

} // namespace lar_dl_content

