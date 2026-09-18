/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoStitching/PfoStitchingAlgorithm.cc
 *
 *  @brief  Implementation of the track pfo stitching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Geometry/DetectorGap.h"
#include "Geometry/LArTPC.h"

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/PfoStitchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoStitchingAlgorithm::InferredLArTPC::InferredLArTPC(const PandoraApi::Geometry::LArTPC::Parameters &parameters) :
    LArTPC(parameters)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoStitchingAlgorithm::PfoStitchingAlgorithm() :
    m_visualise{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoStitchingAlgorithm::Run()
{
    PANDORA_RETURN_IF(STATUS_CODE_SUCCESS, m_stitchingToolVector.empty());

    PfoList sortedPfos;
    for (const std::string &listName : m_inputPfoListNames)
    {
        const PfoList *pPfos{nullptr};
        if (PandoraContentApi::GetList(*this, listName, pPfos) == STATUS_CODE_SUCCESS)
            sortedPfos.insert(sortedPfos.end(), pPfos->begin(), pPfos->end());
    }
    sortedPfos.sort(LArPfoHelper::SortByNHits);
    const PfoList *const pSortedPfoList{&sortedPfos};

    InferredLArTPCVector inferredTPCs;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->InferLArTPCs(inferredTPCs));
    if (inferredTPCs.size() == 1) // Might be in a per-TPC worker
        return STATUS_CODE_SUCCESS;

    PfoToLArTPCMap pfoToLArTPCMap;
    this->GetPfoToLArTPCMap(pSortedPfoList, inferredTPCs, pfoToLArTPCMap);
    if (pfoToLArTPCMap.empty())
        return STATUS_CODE_SUCCESS;

    PfoToFloatMap stitchedPfosToX0Map;
    for (StitchingBaseTool *const pStitchingTool : m_stitchingToolVector)
        pStitchingTool->Run(this, pSortedPfoList, pfoToLArTPCMap, stitchedPfosToX0Map);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoStitchingAlgorithm::InferLArTPCs(InferredLArTPCVector &inferredTPCs) const
{
    const LArTPC &parentLArTPC(this->GetPandora().GetGeometry()->GetLArTPC());
    const float parentMinX(parentLArTPC.GetCenterX() - 0.5f * parentLArTPC.GetWidthX());
    const float parentMaxX(parentLArTPC.GetCenterX() + 0.5f * parentLArTPC.GetWidthX());

    std::vector<std::pair<float, float>> driftGapXIntervals;
    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pDetectorGap));

        if (!pLineGap || (pLineGap->GetLineGapType() != TPC_DRIFT_GAP))
            continue;

        PANDORA_RETURN_IF(STATUS_CODE_INVALID_PARAMETER, pLineGap->GetLineStartX() > pLineGap->GetLineEndX());

        if (pLineGap->GetLineStartX() < parentMinX || pLineGap->GetLineEndX() > parentMaxX)
            continue;

        driftGapXIntervals.emplace_back(pLineGap->GetLineStartX(), pLineGap->GetLineEndX());
    }

    std::sort(driftGapXIntervals.begin(), driftGapXIntervals.end(),
        [](const std::pair<float, float> &lhs, const std::pair<float, float> &rhs) { return lhs.first < rhs.first; });

    std::vector<std::pair<float, float>> inferredTPCXIntervals;
    float inferredTPCMinX{parentMinX};
    for (const std::pair<float, float> &interval : driftGapXIntervals)
    {
        inferredTPCXIntervals.emplace_back(inferredTPCMinX, interval.first);
        inferredTPCMinX = interval.second;
    }
    inferredTPCXIntervals.emplace_back(inferredTPCMinX, parentMaxX);

    unsigned int volumeId{0};
    bool isDriftInPositiveX{true};
    for (const std::pair<float, float> &interval : inferredTPCXIntervals)
    {
        PandoraApi::Geometry::LArTPC::Parameters parameters;
        parameters.m_larTPCVolumeId = volumeId++;
        parameters.m_centerX = 0.5f * (interval.first + interval.second);
        parameters.m_centerY = parentLArTPC.GetCenterY();
        parameters.m_centerZ = parentLArTPC.GetCenterZ();
        parameters.m_widthX = interval.second - interval.first;
        parameters.m_widthY = parentLArTPC.GetWidthY();
        parameters.m_widthZ = parentLArTPC.GetWidthZ();
        parameters.m_wirePitchU = parentLArTPC.GetWirePitchU();
        parameters.m_wirePitchV = parentLArTPC.GetWirePitchV();
        parameters.m_wirePitchW = parentLArTPC.GetWirePitchW();
        parameters.m_wireAngleU = parentLArTPC.GetWireAngleU();
        parameters.m_wireAngleV = parentLArTPC.GetWireAngleV();
        parameters.m_wireAngleW = parentLArTPC.GetWireAngleW();
        parameters.m_sigmaUVW = parentLArTPC.GetSigmaUVW();
        parameters.m_isDriftInPositiveX = isDriftInPositiveX;
        isDriftInPositiveX = !isDriftInPositiveX;

        inferredTPCs.push_back(std::unique_ptr<InferredLArTPC>(new InferredLArTPC(parameters)));
    }

    PANDORA_RETURN_IF(STATUS_CODE_FAILURE, inferredTPCs.empty());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoStitchingAlgorithm::GetPfoToLArTPCMap(const PfoList *const pPfoList, const InferredLArTPCVector &inferredTPCs, PfoToLArTPCMap &pfoToLArTPCMap) const
{
    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
            continue;

        CartesianPointVector positions;
        LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, positions);

        // Find the LArTPC for PFOs contained in a single LArTPC
        // ATTN We should try relaxing this to stitch particles that span > 2 LArTPCs, would require dev in stitching tool & its helpers
        const LArTPC *pLArTPC{nullptr};
        for (const CartesianVector &pos : positions)
        {
            const LArTPC *pPosLArTPC{nullptr};
            for (const std::unique_ptr<InferredLArTPC> &pInferredLArTPC : inferredTPCs)
            {
                const float minX{pInferredLArTPC->GetCenterX() - 0.5f * pInferredLArTPC->GetWidthX()};
                const float maxX{pInferredLArTPC->GetCenterX() + 0.5f * pInferredLArTPC->GetWidthX()};

                if (pos.GetX() >= minX && pos.GetX() <= maxX)
                {
                    pPosLArTPC = pInferredLArTPC.get();
                    break;
                }
            }
            if (!pPosLArTPC) // 3D hit ended up in a gap or outside detector, ignoring it
                continue;

            if (pLArTPC && pLArTPC != pPosLArTPC) // Hits in > 1 LArTPC so PFO is not a stitching candidate, giving up
            {
                pLArTPC = nullptr;
                break;
            }

            pLArTPC = pPosLArTPC;
        }

        if (pLArTPC) // The Pfo is within a single LArTPC
            pfoToLArTPCMap.emplace(pPfo, pLArTPC);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string PfoStitchingAlgorithm::GetListName(const Cluster *const pCluster) const
{
    return this->GetListName<Cluster, ClusterList>(pCluster, m_daughterListNames);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string PfoStitchingAlgorithm::GetListName(const Vertex *const pVertex) const
{
    return this->GetListName<Vertex, VertexList>(pVertex, m_daughterListNames);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string PfoStitchingAlgorithm::GetListName(const ParticleFlowObject *const pPfo) const
{
    return this->GetListName<ParticleFlowObject, PfoList>(pPfo, m_inputPfoListNames);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TObject, typename TList>
const std::string PfoStitchingAlgorithm::GetListName(const TObject *const pObject, const StringVector &listNames) const
{
    for (const std::string &listName : listNames)
    {
        const TList *pObjectList{nullptr};
        if ((PandoraContentApi::GetList(*this, listName, pObjectList) == STATUS_CODE_SUCCESS) &&
            (pObjectList->end() != std::find(pObjectList->begin(), pObjectList->end(), pObject)))
        {
            return listName;
        }
    }
    PANDORA_THROW(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoStitchingAlgorithm::ShiftPfoHierarchy([[maybe_unused]] const ParticleFlowObject *const pParentPfo,
    [[maybe_unused]] const PfoToLArTPCMap &pfoToLArTPCMap, [[maybe_unused]] const float x0) const
{
    std::cout << "PfoStitchingAlgorithm: Stitching tools must not be configured to apply an x0 correction" << std::endl;
    PANDORA_THROW(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoStitchingAlgorithm::StitchPfos(
    const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete, PfoToLArTPCMap &pfoToLArTPCMap) const
{
    PANDORA_THROW_IF(STATUS_CODE_NOT_ALLOWED, pPfoToEnlarge == pPfoToDelete);

    if (m_visualise)
    {
        PfoList pfoToEnlargeList{pPfoToEnlarge};
        PfoList pfoToDeleteList;
        LArPfoHelper::GetAllDownstreamPfos(pPfoToDelete, pfoToDeleteList);
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoToEnlargeList, "PfoToEnlarge", GREEN, true, false));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoToDeleteList, "PfoToDelete", RED, true, false));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    pfoToLArTPCMap.erase(pPfoToEnlarge);
    pfoToLArTPCMap.erase(pPfoToDelete);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, this->GetListName(pPfoToDelete)));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));
    }

    for (const Vertex *const pDaughterVertex : daughterVertices)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, this->GetListName(pDaughterVertex)));
    }

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *const pParentCluster(PfoMopUpBaseAlgorithm::GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::MergeAndDeleteClusters(
                    *this, pParentCluster, pDaughterCluster, this->GetListName(pParentCluster), this->GetListName(pDaughterCluster)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoStitchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterListNames", m_daughterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "StitchingTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
    {
        StitchingBaseTool *const pStitchingTool(dynamic_cast<StitchingBaseTool *>(pAlgorithmTool));

        PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, !pStitchingTool);

        m_stitchingToolVector.push_back(pStitchingTool);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
