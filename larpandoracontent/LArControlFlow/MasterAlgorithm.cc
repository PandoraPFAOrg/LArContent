/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MasterAlgorithm::MasterAlgorithm() :
    m_workerInstancesInitialized(false),
    m_larCaloHitVersion(1),
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(false),
    m_visualizeOverallRecoStatus(false),
    m_shouldRemoveOutOfTimeHits(true),
    m_pSlicingWorkerInstance(nullptr),
    m_pSliceNuWorkerInstance(nullptr),
    m_pSliceCRWorkerInstance(nullptr),
    m_fullWidthCRWorkerWireGaps(true),
    m_passMCParticlesToWorkerInstances(false),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_inTimeMaxX0(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::ShiftPfoHierarchy(const ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const
{
    if (!pParentPfo->GetParentPfoList().empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PfoToLArTPCMap::const_iterator larTPCIter(pfoToLArTPCMap.find(pParentPfo));

    if (pfoToLArTPCMap.end() == larTPCIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    PfoList pfoList;
    LArPfoHelper::GetAllDownstreamPfos(pParentPfo, pfoList);

    if (m_visualizeOverallRecoStatus)
    {
        std::cout << "ShiftPfoHierarchy: x0 " << x0 << std::endl;
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "BeforeShiftCRPfos", GREEN));
    }

    for (const ParticleFlowObject *const pDaughterPfo : pfoList)
    {
        CaloHitList caloHitList;
        for (const Cluster *const pCluster : pDaughterPfo->GetClusterList())
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        }

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            PandoraContentApi::CaloHit::Metadata metadata;
            metadata.m_x0 = x0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
        }

        for (const Vertex *const pVertex : pDaughterPfo->GetVertexList())
        {
            PandoraContentApi::Vertex::Metadata metadata;
            metadata.m_x0 = x0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::AlterMetadata(*this, pVertex, metadata));
        }
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "AfterShiftCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::StitchPfos(
    const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete, PfoToLArTPCMap &pfoToLArTPCMap) const
{
    if (pPfoToEnlarge == pPfoToDelete)
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN Remove pfos from pfo to lar tpc map here, avoiding problems if stitching across multiple tpcs.
    pfoToLArTPCMap.erase(pPfoToEnlarge);
    pfoToLArTPCMap.erase(pPfoToDelete);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, m_recreatedPfoListName));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));

    for (const Vertex *const pDaughterVertex : daughterVertices)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, m_recreatedVertexListName));

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(PfoMopUpBaseAlgorithm::GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster, m_recreatedClusterListName, m_recreatedClusterListName));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Run()
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Reset());

    if (!m_workerInstancesInitialized)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->InitializeWorkerInstances());

    if (m_passMCParticlesToWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CopyMCParticles());

    if (m_pSliceNuWorkerInstance)
        PandoraApi::SetEventInformation(
            *m_pSliceNuWorkerInstance, this->GetPandora().GetRun(), this->GetPandora().GetSubrun(), this->GetPandora().GetEvent());

    PfoToFloatMap stitchedPfosToX0Map;
    VolumeIdToHitListMap volumeIdToHitListMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVolumeIdToHitListMap(volumeIdToHitListMap));

    if (m_shouldRunAllHitsCosmicReco)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayReconstruction(volumeIdToHitListMap));

        PfoToLArTPCMap pfoToLArTPCMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RecreateCosmicRayPfos(pfoToLArTPCMap));

        if (m_shouldRunStitching)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->StitchCosmicRayPfos(pfoToLArTPCMap, stitchedPfosToX0Map));
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        PfoList clearCosmicRayPfos, ambiguousPfos;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayHitRemoval(ambiguousPfos));
    }

    SliceVector sliceVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSlicing(volumeIdToHitListMap, sliceVector));

    if (m_shouldRunNeutrinoRecoOption || m_shouldRunCosmicRecoOption)
    {
        SliceHypotheses nuSliceHypotheses, crSliceHypotheses;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSliceReconstruction(sliceVector, nuSliceHypotheses, crSliceHypotheses));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SelectBestSliceHypotheses(nuSliceHypotheses, crSliceHypotheses));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::InitializeWorkerInstances()
{
    // ATTN Used to be in the regular Initialize callback, but detector gap list cannot be extracted in client app before the first event
    if (m_workerInstancesInitialized)
        return STATUS_CODE_ALREADY_INITIALIZED;

    try
    {
        const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
        const DetectorGapList &gapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

        for (const LArTPCMap::value_type &mapEntry : larTPCMap)
        {
            const unsigned int volumeId(mapEntry.second->GetLArTPCVolumeId());
            m_crWorkerInstances.push_back(
                this->CreateWorkerInstance(*(mapEntry.second), gapList, m_crSettingsFile, "CRWorkerInstance" + std::to_string(volumeId)));
        }

        if (m_shouldRunSlicing)
            m_pSlicingWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_slicingSettingsFile, "SlicingWorker");

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile, "SliceNuWorker");

        if (m_shouldRunCosmicRecoOption)
            m_pSliceCRWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile, "SliceCRWorker");
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "MasterAlgorithm: Exception during initialization of worker instances " << statusCodeException.ToString() << std::endl;
        return statusCodeException.GetStatusCode();
    }

    m_workerInstancesInitialized = true;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::CopyMCParticles() const
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    PandoraInstanceList pandoraWorkerInstances(m_crWorkerInstances);
    if (m_pSlicingWorkerInstance)
        pandoraWorkerInstances.push_back(m_pSlicingWorkerInstance);
    if (m_pSliceNuWorkerInstance)
        pandoraWorkerInstances.push_back(m_pSliceNuWorkerInstance);
    if (m_pSliceCRWorkerInstance)
        pandoraWorkerInstances.push_back(m_pSliceCRWorkerInstance);

    LArMCParticleFactory mcParticleFactory;

    for (const Pandora *const pPandoraWorker : pandoraWorkerInstances)
    {
        for (const MCParticle *const pMCParticle : *pMCParticleList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pPandoraWorker, pMCParticle, &mcParticleFactory));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const unsigned int nLArTPCs(larTPCMap.size());

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));

        if (!pLArCaloHit && (1 != nLArTPCs))
            return STATUS_CODE_INVALID_PARAMETER;

        const unsigned int volumeId(pLArCaloHit ? pLArCaloHit->GetLArTPCVolumeId() : 0);
        const LArTPC *const pLArTPC(larTPCMap.at(volumeId));

        LArTPCHitList &larTPCHitList(volumeIdToHitListMap[volumeId]);
        larTPCHitList.m_allHitList.push_back(pCaloHit);

        if (((pCaloHit->GetPositionVector().GetX() >= (pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX())) &&
                (pCaloHit->GetPositionVector().GetX() <= (pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX()))))
        {
            larTPCHitList.m_truncatedHitList.push_back(pCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    unsigned int workerCounter(0);

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());
        VolumeIdToHitListMap::const_iterator iter(volumeIdToHitListMap.find(larTPC.GetLArTPCVolumeId()));

        if (volumeIdToHitListMap.end() == iter)
            continue;

        for (const CaloHit *const pCaloHit : iter->second.m_allHitList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pCRWorker, pCaloHit));

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance " << ++workerCounter << " of " << m_crWorkerInstances.size() << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pCRWorker));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RecreateCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const PfoList *pCRPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pCRWorker, pCRPfos));

        PfoList newPfoList;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(*pCRPfos, newPfoList));

        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());

        for (const Pfo *const pNewPfo : newPfoList)
            pfoToLArTPCMap[pNewPfo] = &larTPC;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::StitchCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "RecreatedCRPfos", GREEN));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    for (StitchingBaseTool *const pStitchingTool : m_stitchingToolVector)
        pStitchingTool->Run(this, pRecreatedCRPfos, pfoToLArTPCMap, stitchedPfosToX0Map);

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "AfterStitchingCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::TagCosmicRayPfos(const PfoToFloatMap &stitchedPfosToX0Map, PfoList &clearCosmicRayPfos, PfoList &ambiguousPfos) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    PfoList nonStitchedParentCosmicRayPfos;
    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        if (!pPfo->GetParentPfoList().empty())
            continue;

        PfoToFloatMap::const_iterator pfoToX0Iter = stitchedPfosToX0Map.find(pPfo);
        const float x0Shift((pfoToX0Iter != stitchedPfosToX0Map.end()) ? pfoToX0Iter->second : 0.f);
        PfoList &targetList((std::fabs(x0Shift) > m_inTimeMaxX0) ? clearCosmicRayPfos : nonStitchedParentCosmicRayPfos);
        targetList.push_back(pPfo);
    }

    for (CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool : m_cosmicRayTaggingToolVector)
        pCosmicRayTaggingTool->FindAmbiguousPfos(nonStitchedParentCosmicRayPfos, ambiguousPfos, this);

    for (const Pfo *const pPfo : nonStitchedParentCosmicRayPfos)
    {
        const bool isClearCosmic(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo));

        if (isClearCosmic)
            clearCosmicRayPfos.push_back(pPfo);
    }

    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        const bool isClearCosmic(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo));
        PandoraContentApi::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["IsClearCosmic"] = (isClearCosmic ? 1.f : 0.f);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &clearCosmicRayPfos, "ClearCRPfos", RED));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &ambiguousPfos, "AmbiguousCRPfos", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayHitRemoval(const PfoList &ambiguousPfos) const
{
    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(ambiguousPfos, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        const ClusterList clusterList(pPfoToDelete->GetClusterList());
        const VertexList vertexList(pPfoToDelete->GetVertexList());

        // ATTN: If an ambiguous pfo has been stitched, reset the calo hit positions in preparation for subsequent algorithm chains
        if (LArStitchingHelper::HasPfoBeenStitched(pPfoToDelete))
        {
            CaloHitList caloHitList2D;
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_W, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_W, caloHitList2D);

            for (const CaloHit *const pCaloHit : caloHitList2D)
            {
                PandoraContentApi::CaloHit::Metadata metadata;
                metadata.m_x0 = 0.f;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
            }
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &vertexList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const
{
    for (const VolumeIdToHitListMap::value_type &mapEntry : volumeIdToHitListMap)
    {
        for (const CaloHit *const pCaloHit : (m_shouldRemoveOutOfTimeHits ? mapEntry.second.m_truncatedHitList : mapEntry.second.m_allHitList))
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            const HitType hitType(pCaloHit->GetHitType());
            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                continue;

            if (m_shouldRunSlicing)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSlicingWorkerInstance, pCaloHit));
            }
            else
            {
                if (sliceVector.empty())
                    sliceVector.push_back(CaloHitList());
                sliceVector.back().push_back(pCaloHit);
            }
        }
    }

    if (m_shouldRunSlicing)
    {
        if (m_printOverallRecoStatus)
            std::cout << "Running slicing worker instance" << std::endl;

        const PfoList *pSlicePfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSlicingWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSlicingWorkerInstance, pSlicePfos));

        if (m_visualizeOverallRecoStatus)
        {
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pSlicePfos, "OnePfoPerSlice", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        for (const Pfo *const pSlicePfo : *pSlicePfos)
        {
            sliceVector.push_back(CaloHitList());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_U, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_V, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_W, sliceVector.back());
        }
    }

    if (m_printOverallRecoStatus)
        std::cout << "Identified " << sliceVector.size() << " slice(s)" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSliceReconstruction(SliceVector &sliceVector, SliceHypotheses &nuSliceHypotheses, SliceHypotheses &crSliceHypotheses) const
{
    SliceVector selectedSliceVector;
    if (m_shouldRunSlicing && !m_sliceSelectionToolVector.empty())
    {
        SliceVector inputSliceVector(sliceVector);
        for (SliceSelectionBaseTool *const pSliceSelectionTool : m_sliceSelectionToolVector)
        {
            pSliceSelectionTool->SelectSlices(this, inputSliceVector, selectedSliceVector);
            inputSliceVector = selectedSliceVector;
        }
    }
    else
    {
        selectedSliceVector = std::move(sliceVector);
    }

    unsigned int sliceCounter(0);

    for (const CaloHitList &sliceHits : selectedSliceVector)
    {
        for (const CaloHit *const pSliceCaloHit : sliceHits)
        {
            // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
            const CaloHit *const pCaloHitInMaster(m_shouldRunSlicing ? static_cast<const CaloHit *>(pSliceCaloHit->GetParentAddress()) : pSliceCaloHit);

            if (m_shouldRunNeutrinoRecoOption)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceNuWorkerInstance, pCaloHitInMaster));

            if (m_shouldRunCosmicRecoOption)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceCRWorkerInstance, pCaloHitInMaster));
        }

        if (m_shouldRunNeutrinoRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running nu worker instance for slice " << (sliceCounter + 1) << " of " << selectedSliceVector.size() << std::endl;

            const PfoList *pSliceNuPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceNuWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceNuWorkerInstance, pSliceNuPfos));
            nuSliceHypotheses.push_back(*pSliceNuPfos);

            for (const ParticleFlowObject *const pPfo : *pSliceNuPfos)
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["SliceIndex"] = sliceCounter;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
        }

        if (m_shouldRunCosmicRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running cr worker instance for slice " << (sliceCounter + 1) << " of " << selectedSliceVector.size() << std::endl;

            const PfoList *pSliceCRPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceCRWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceCRWorkerInstance, pSliceCRPfos));
            crSliceHypotheses.push_back(*pSliceCRPfos);

            for (const ParticleFlowObject *const pPfo : *pSliceCRPfos)
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["SliceIndex"] = sliceCounter;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
        }

        ++sliceCounter;
    }

    // ATTN: If we swapped these objects at the start, be sure to swap them back in case we ever want to use sliceVector
    // after this function
    if (!(m_shouldRunSlicing && !m_sliceSelectionToolVector.empty()))
        sliceVector = std::move(selectedSliceVector);

    if (m_shouldRunNeutrinoRecoOption && m_shouldRunCosmicRecoOption && (nuSliceHypotheses.size() != crSliceHypotheses.size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::SelectBestSliceHypotheses(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses) const
{
    if (m_printOverallRecoStatus)
        std::cout << "Select best slice hypotheses" << std::endl;

    PfoList selectedSlicePfos;

    if (m_shouldPerformSliceId)
    {
        for (SliceIdBaseTool *const pSliceIdTool : m_sliceIdToolVector)
            pSliceIdTool->SelectOutputPfos(this, nuSliceHypotheses, crSliceHypotheses, selectedSlicePfos);
    }
    else if (m_shouldRunNeutrinoRecoOption != m_shouldRunCosmicRecoOption)
    {
        const SliceHypotheses &sliceHypotheses(m_shouldRunNeutrinoRecoOption ? nuSliceHypotheses : crSliceHypotheses);

        for (const PfoList &slice : sliceHypotheses)
            selectedSlicePfos.insert(selectedSlicePfos.end(), slice.begin(), slice.end());
    }

    PfoList newSlicePfoList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(selectedSlicePfos, newSlicePfoList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Reset()
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pCRWorker));

    if (m_pSlicingWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSlicingWorkerInstance));

    if (m_pSliceNuWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceNuWorkerInstance));

    if (m_pSliceCRWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceCRWorkerInstance));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const CaloHit *const pCaloHit) const
{
    const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
    if (pLArCaloHit == nullptr)
    {
        std::cout << "MasterAlgorithm: Could not cast CaloHit to LArCaloHit" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    LArCaloHitParameters parameters;
    pLArCaloHit->FillParameters(parameters);

    // For merged worker instances, remap the hit's volume ids to match the merged LArTPC's renumbered readout volumes
    if ((pPandora == m_pSlicingWorkerInstance) || (pPandora == m_pSliceNuWorkerInstance) || (pPandora == m_pSliceCRWorkerInstance))
    {
        const auto offsetIter(m_daughterVolumeIdOffsetMap.find(pLArCaloHit->GetLArTPCVolumeId()));

        if (m_daughterVolumeIdOffsetMap.end() == offsetIter)
            return STATUS_CODE_FAILURE;

        parameters.m_larTPCVolumeId = 0;
        parameters.m_daughterVolumeId = offsetIter->second + pLArCaloHit->GetDaughterVolumeId();
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, parameters, m_larCaloHitFactory));

    if (m_passMCParticlesToWorkerInstances)
    {
        MCParticleVector mcParticleVector;
        for (const auto &weightMapEntry : pLArCaloHit->GetMCParticleWeightMap())
            mcParticleVector.push_back(weightMapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetCaloHitToMCParticleRelationship(*pPandora, pLArCaloHit, pMCParticle, pLArCaloHit->GetMCParticleWeightMap().at(pMCParticle)));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const MCParticle *const pMCParticle, const LArMCParticleFactory *const pMCParticleFactory) const
{
    const LArMCParticle *const pLArMCParticle = dynamic_cast<const LArMCParticle *>(pMCParticle);

    if (!pLArMCParticle)
    {
        std::cout << "MasterAlgorithm::Copy - Expect to pass only LArMCParticles to Pandora worker instances." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    LArMCParticleParameters parameters;
    pLArMCParticle->FillParameters(parameters);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, parameters, *pMCParticleFactory));

    for (const MCParticle *const pDaughterMCParticle : pMCParticle->GetDaughterList())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pMCParticle, pDaughterMCParticle));

    for (const MCParticle *const pParentMCParticle : pMCParticle->GetParentList())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pParentMCParticle, pMCParticle));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const PfoList &inputPfoList, PfoList &newPfoList) const
{
    if (inputPfoList.empty())
        return STATUS_CODE_SUCCESS;

    // TODO if no pfo in input list is primary - raise exception

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Pfo *const pPfo : inputPfoList)
    {
        if (pPfo->GetParentPfoList().empty())
            this->Recreate(pPfo, nullptr, newPfoList);
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_recreatedClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_recreatedClusterListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_recreatedVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_recreatedVertexListName));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, m_recreatedPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, m_recreatedPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(static_cast<const CaloHit *>(pInputCaloHit->GetParentAddress()));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(static_cast<const CaloHit *>(pInputCaloHit->GetParentAddress()));

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit *>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit *>(pWorkerParentCaloHit->GetParentAddress()));
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit *>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit *>(pWorkerParentCaloHit->GetParentAddress()));
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *MasterAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const CaloHit *const pParentCaloHit) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void *>(pParentCaloHit);

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *MasterAlgorithm::CreateCluster(
    const Cluster *const pInputCluster, const CaloHitList &newCaloHitList, const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *MasterAlgorithm::CreateVertex(const Vertex *const pInputVertex) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition();
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *MasterAlgorithm::CreatePfo(
    const ParticleFlowObject *const pInputPfo, const ClusterList &newClusterList, const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;
    parameters.m_propertiesToAdd = pInputPfo->GetPropertiesMap();

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(
    const LArTPC &larTPC, const DetectorGapList &gapList, const std::string &settingsFile, const std::string &name) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora(name));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterCustomContent(pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraApi::SetEventInformation(*pPandora, this->GetPandora().GetRun(), this->GetPandora().GetSubrun(), this->GetPandora().GetEvent()));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The LArTPC
    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = larTPC.GetLArTPCVolumeId();
    larTPCParameters.m_centerX = larTPC.GetCenterX();
    larTPCParameters.m_centerY = larTPC.GetCenterY();
    larTPCParameters.m_centerZ = larTPC.GetCenterZ();
    larTPCParameters.m_widthX = larTPC.GetWidthX();
    larTPCParameters.m_widthY = larTPC.GetWidthY();
    larTPCParameters.m_widthZ = larTPC.GetWidthZ();
    larTPCParameters.m_wirePitchU = larTPC.GetWirePitchU();
    larTPCParameters.m_wirePitchV = larTPC.GetWirePitchV();
    larTPCParameters.m_wirePitchW = larTPC.GetWirePitchW();
    larTPCParameters.m_wireAngleU = larTPC.GetWireAngleU();
    larTPCParameters.m_wireAngleV = larTPC.GetWireAngleV();
    larTPCParameters.m_wireAngleW = larTPC.GetWireAngleW();
    larTPCParameters.m_sigmaUVW = larTPC.GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = larTPC.IsDriftInPositiveX();
    this->AppendReadoutVolumeParameters(larTPC, 0, larTPCParameters.m_readoutVolumeParametersVector);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    const float tpcMinX(larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX()), tpcMaxX(larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX());

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pGap));

        if (pLineGap &&
            (((pLineGap->GetLineEndX() >= tpcMinX) && (pLineGap->GetLineEndX() <= tpcMaxX)) ||
                ((pLineGap->GetLineStartX() >= tpcMinX) && (pLineGap->GetLineStartX() <= tpcMaxX))))
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            const LineGapType lineGapType(pLineGap->GetLineGapType());
            lineGapParameters.m_lineGapType = lineGapType;
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();

            if (m_fullWidthCRWorkerWireGaps &&
                ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W)))
            {
                lineGapParameters.m_lineStartX = -std::numeric_limits<float>::max();
                lineGapParameters.m_lineEndX = std::numeric_limits<float>::max();
            }

            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(
    const LArTPCMap &larTPCMap, const DetectorGapList &gapList, const std::string &settingsFile, const std::string &name) const
{
    if (larTPCMap.empty())
    {
        std::cout << "MasterAlgorithm::CreateWorkerInstance - no LArTPC details provided" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    // The Pandora instance
    const Pandora *const pPandora(new Pandora(name));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterCustomContent(pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraApi::SetEventInformation(*pPandora, this->GetPandora().GetRun(), this->GetPandora().GetSubrun(), this->GetPandora().GetEvent()));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The Parent LArTPC
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = 0;
    larTPCParameters.m_centerX = 0.5f * (parentMaxX + parentMinX);
    larTPCParameters.m_centerY = 0.5f * (parentMaxY + parentMinY);
    larTPCParameters.m_centerZ = 0.5f * (parentMaxZ + parentMinZ);
    larTPCParameters.m_widthX = parentMaxX - parentMinX;
    larTPCParameters.m_widthY = parentMaxY - parentMinY;
    larTPCParameters.m_widthZ = parentMaxZ - parentMinZ;
    larTPCParameters.m_wirePitchU = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchV = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchW = pFirstLArTPC->GetWirePitchW();
    larTPCParameters.m_wireAngleU = pFirstLArTPC->GetWireAngleU();
    larTPCParameters.m_wireAngleV = pFirstLArTPC->GetWireAngleV();
    larTPCParameters.m_wireAngleW = pFirstLArTPC->GetWireAngleW();
    larTPCParameters.m_sigmaUVW = pFirstLArTPC->GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = pFirstLArTPC->IsDriftInPositiveX();

    // Merge in readout volumes from every child TPC, keeping ids unique, and record the offset used for each original LArTPC so hits copied
    // into this worker can be remapped consistently in Copy()
    unsigned int idOffset(0);
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        // ATTN: This gets rebuilt on each call at startup (slicing), but it's always consistent, so it's fine to do this.
        m_daughterVolumeIdOffsetMap[pLArTPC->GetLArTPCVolumeId()] = idOffset;
        this->AppendReadoutVolumeParameters(*pLArTPC, idOffset, larTPCParameters.m_readoutVolumeParametersVector);
        idOffset += pLArTPC->GetReadoutVolumes().size();
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pGap));

        if (pLineGap)
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            lineGapParameters.m_lineGapType = pLineGap->GetLineGapType();
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();
            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::AppendReadoutVolumeParameters(const LArTPC &larTPC, const unsigned int idOffset,
    object_creation::LArReadoutVolumeParametersVector &readoutVolumeParametersVector) const
{
    for (const auto &[_, readoutVolume] : larTPC.GetReadoutVolumes())
    {
        object_creation::LArReadoutVolumeParameters readoutVolumeParams;
        readoutVolumeParams.m_id = idOffset + readoutVolume.GetId();
        readoutVolumeParams.m_center = readoutVolume.GetCenter();
        readoutVolumeParams.m_size = readoutVolume.GetSize();

        for (const LArReadoutUnit &readoutUnit : readoutVolume.GetReadoutUnits())
        {
            object_creation::LArReadoutUnitParameters unitParams;
            unitParams.m_id = readoutUnit.GetId();
            unitParams.m_view = readoutUnit.GetView();
            unitParams.m_referenceCoordinate = readoutUnit.GetReferenceCoordinate();
            unitParams.m_pitch = readoutUnit.GetPitch();
            unitParams.m_unitCenter = readoutUnit.GetUnitCenter();
            unitParams.m_unitSize = readoutUnit.GetUnitSize();

            for (const LArReadoutChannel &channel : readoutUnit.GetReadoutChannels())
            {
                object_creation::LArReadoutChannelParameters channelParams;
                channelParams.m_id = channel.GetId();
                channelParams.m_channelIntervalArray = channel.GetChannelIntervals();
                unitParams.m_channelParametersVector.push_back(channelParams);
            }
            readoutVolumeParams.m_readoutUnitParametersVector.push_back(unitParams);
        }
        readoutVolumeParametersVector.push_back(readoutVolumeParams);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RegisterCustomContent(const Pandora *const /*pPandora*/) const
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    ExternalSteeringParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        pExternalParameters = dynamic_cast<ExternalSteeringParameters *>(this->GetExternalParameters());
        if (!pExternalParameters)
            return STATUS_CODE_FAILURE;
    }

    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "SliceSelectionTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            SliceSelectionBaseTool *const pSliceSelectionTool(dynamic_cast<SliceSelectionBaseTool *>(pAlgorithmTool));
            if (!pSliceSelectionTool)
                return STATUS_CODE_INVALID_PARAMETER;
            m_sliceSelectionToolVector.push_back(pSliceSelectionTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldRunAllHitsCosmicReco,
            xmlHandle, "ShouldRunAllHitsCosmicReco", m_shouldRunAllHitsCosmicReco));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldRunStitching,
            xmlHandle, "ShouldRunStitching", m_shouldRunStitching));

    if (m_shouldRunStitching && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunStitching requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunStitching)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "StitchingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            StitchingBaseTool *const pStitchingTool(dynamic_cast<StitchingBaseTool *>(pAlgorithmTool));
            if (!pStitchingTool)
                return STATUS_CODE_INVALID_PARAMETER;
            m_stitchingToolVector.push_back(pStitchingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldRunCosmicHitRemoval,
            xmlHandle, "ShouldRunCosmicHitRemoval", m_shouldRunCosmicHitRemoval));

    if (m_shouldRunCosmicHitRemoval && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunCosmicHitRemoval requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "CosmicRayTaggingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool(dynamic_cast<CosmicRayTaggingBaseTool *>(pAlgorithmTool));
            if (!pCosmicRayTaggingTool)
                return STATUS_CODE_INVALID_PARAMETER;
            m_cosmicRayTaggingToolVector.push_back(pCosmicRayTaggingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldRunSlicing,
            xmlHandle, "ShouldRunSlicing", m_shouldRunSlicing));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldRunNeutrinoRecoOption,
            xmlHandle, "ShouldRunNeutrinoRecoOption", m_shouldRunNeutrinoRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldRunCosmicRecoOption,
            xmlHandle, "ShouldRunCosmicRecoOption", m_shouldRunCosmicRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_shouldPerformSliceId,
            xmlHandle, "ShouldPerformSliceId", m_shouldPerformSliceId));

    if (m_shouldPerformSliceId && (!m_shouldRunSlicing || !m_shouldRunNeutrinoRecoOption || !m_shouldRunCosmicRecoOption))
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldPerformSliceId requires ShouldRunSlicing and both neutrino and cosmic reconstruction options"
                  << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldPerformSliceId)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "SliceIdTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            SliceIdBaseTool *const pSliceIdIdTool(dynamic_cast<SliceIdBaseTool *>(pAlgorithmTool));
            if (!pSliceIdIdTool)
                return STATUS_CODE_INVALID_PARAMETER;
            m_sliceIdToolVector.push_back(pSliceIdIdTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() : pExternalParameters->m_printOverallRecoStatus,
            xmlHandle, "PrintOverallRecoStatus", m_printOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "VisualizeOverallRecoStatus", m_visualizeOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LArCaloHitVersion", m_larCaloHitVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldRemoveOutOfTimeHits", m_shouldRemoveOutOfTimeHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FullWidthCRWorkerWireGaps", m_fullWidthCRWorkerWireGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "PassMCParticlesToWorkerInstances", m_passMCParticlesToWorkerInstances));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CRSettingsFile", m_crSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuSettingsFile", m_nuSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SlicingSettingsFile", m_slicingSettingsFile));
    m_crSettingsFile = LArFileHelper::FindFileInPath(m_crSettingsFile, m_filePathEnvironmentVariable);
    m_nuSettingsFile = LArFileHelper::FindFileInPath(m_nuSettingsFile, m_filePathEnvironmentVariable);
    m_slicingSettingsFile = LArFileHelper::FindFileInPath(m_slicingSettingsFile, m_filePathEnvironmentVariable);

    if (m_passMCParticlesToWorkerInstances)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputMCParticleListName", m_inputMCParticleListName));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputHitListName", m_inputHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedPfoListName", m_recreatedPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedClusterListName", m_recreatedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedVertexListName", m_recreatedVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMaxX0", m_inTimeMaxX0));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const InputBool inputBool,
    const TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool)
{
    if (pExternalParameters && inputBool.IsInitialized())
    {
        outputBool = inputBool.Get();
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, xmlTag, outputBool));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
