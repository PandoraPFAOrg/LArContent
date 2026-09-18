/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlTrackCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning based cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Objects/CaloHit.h"
#include "Objects/ParticleFlowObject.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArTrackShowerId/DlTrackCharacterisationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.cc"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.cc"
#include "larpandoracontent/LArHelpers/LArPfoHelper.cc"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlTrackCharacterisationAlgorithm::DlTrackCharacterisationAlgorithm() :
    m_trainingMode{false},
    m_minTrackHits{10},
    m_sequenceLength{128},
    m_maxChargeThreshold{10.f},
    m_detEdgeTolerance{2.5f},
    m_trainingFileName{""},
    m_trainingTreeName{""},
    m_detLimitMinX{std::numeric_limits<float>::max()},
    m_detLimitMaxX{std::numeric_limits<float>::lowest()},
    m_detLimitMinY{std::numeric_limits<float>::max()},
    m_detLimitMaxY{std::numeric_limits<float>::lowest()},
    m_detLimitMinZ{std::numeric_limits<float>::max()},
    m_detLimitMaxZ{std::numeric_limits<float>::lowest()}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlTrackCharacterisationAlgorithm::~DlTrackCharacterisationAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_trainingFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlTrackCharacterisationAlgorithm: Failed to save ROOT tree" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlTrackCharacterisationAlgorithm::Run()
{
    this->GetDetectorLimits();

    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlTrackCharacterisationAlgorithm::PrepareTrainingSample() const
{
    // Get features in U, V and W for all tracks
    PfoToTrackFeaturesMap pfoToTrackFeaturesMap;
    StatusCode status = this->GetAllTrackFeatures(pfoToTrackFeaturesMap);
    if (status == STATUS_CODE_FAILURE)
    {
        std::cout << "PfoList empty... moving on." << std::endl;
        return status;
    }

    for (auto const &[pPfo, trackFeatures] : pfoToTrackFeaturesMap)
    {
        std::map<HitType, std::vector<float>> viewX, viewWire, viewQ;

        for (auto const &[view, hitFeatures] : trackFeatures.GetViewToTrackHitFeaturesMap())
        {
            const unsigned int nHits(hitFeatures.size());
            std::vector<float> x, wire, q;
            x.reserve(nHits);
            wire.reserve(nHits);
            q.reserve(nHits);

            for (auto const &hit : hitFeatures)
            {
                // Set position coords as relative to the first hit
                x.emplace_back(hit.at(0) - hitFeatures.at(0).at(0));
                wire.emplace_back(hit.at(1) - hitFeatures.at(0).at(1));
                q.emplace_back(hit.at(2));
            }
            viewX[view] = x;
            viewWire[view] = wire;
            viewQ[view] = q;
        }
        for (auto &[view, x] : viewX)
        {
            const std::string varName(view == TPC_VIEW_U ? "uViewX" : (view == TPC_VIEW_V ? "vViewX" : "wViewX"));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, varName, &x));
        }
        for (auto &[view, wire] : viewWire)
        {
            const std::string varName(view == TPC_VIEW_U ? "uViewWire" : (view == TPC_VIEW_V ? "vViewWire" : "wViewWire"));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, varName, &wire));
        }
        for (auto &[view, q] : viewQ)
        {
            const std::string varName(view == TPC_VIEW_U ? "uViewQ" : (view == TPC_VIEW_V ? "vViewQ" : "wViewQ"));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, varName, &q));
        }

        const std::vector<float> auxillaryFeatures(trackFeatures.GetAuxillaryFeatures());
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTrackChildren",
            auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNTrkChildren))));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nShowerChildren",
            auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNShwChildren))));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTotalDescendants",
            auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNDescendants))));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nDescendantHits",
            auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNDescendantHits))));
        PANDORA_MONITORING_API(SetTreeVariable(
            this->GetPandora(), m_trainingTreeName, "recoTier", auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kRecoTier))));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "uViewNHits",
            static_cast<int>(auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNHitsU)))));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "vViewNHits",
            static_cast<int>(auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNHitsV)))));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "wViewNHits",
            static_cast<int>(auxillaryFeatures.at(static_cast<unsigned int>(AuxInputs::kNHitsW)))));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "truePdg", trackFeatures.GetTruePdg()));
        PANDORA_MONITORING_API(
            SetTreeVariable(this->GetPandora(), m_trainingTreeName, "isExiting", static_cast<int>(trackFeatures.GetExitingStatus())));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlTrackCharacterisationAlgorithm::Infer()
{
    PfoToTrackFeaturesMap pfoToTrackFeaturesMap;
    StatusCode status = this->GetAllTrackFeatures(pfoToTrackFeaturesMap);
    if (status == STATUS_CODE_FAILURE)
    {
        std::cout << "PfoList empty... moving on." << std::endl;
        return status;
    }

    for (auto const &[pPfo, trackFeatures] : pfoToTrackFeaturesMap)
    {
        // Check if we have at least one view above the minimum number of hits
        const bool uViewGood{trackFeatures.GetNHits(TPC_VIEW_U) >= m_minTrackHits};
        const bool vViewGood{trackFeatures.GetNHits(TPC_VIEW_V) >= m_minTrackHits};
        const bool wViewGood{trackFeatures.GetNHits(TPC_VIEW_W) >= m_minTrackHits};
        if (!uViewGood || !vViewGood || !wViewGood)
            continue;

        LArDLHelper::TorchInput sequenceInputU;
        LArDLHelper::TorchInput sequenceInputV;
        LArDLHelper::TorchInput sequenceInputW;
        const ViewToTrackHitFeaturesMap &trackHitFeatures(trackFeatures.GetViewToTrackHitFeaturesMap());
        this->CreateSequenceInput(trackHitFeatures.at(TPC_VIEW_U), sequenceInputU);
        this->CreateSequenceInput(trackHitFeatures.at(TPC_VIEW_V), sequenceInputV);
        this->CreateSequenceInput(trackHitFeatures.at(TPC_VIEW_W), sequenceInputW);

        LArDLHelper::TorchInput auxillaryInput;
        this->CreateAuxillaryInput(trackFeatures, auxillaryInput);

        LArDLHelper::TorchInputVector inputs;
        inputs.emplace_back(sequenceInputU);
        inputs.emplace_back(sequenceInputV);
        inputs.emplace_back(sequenceInputW);
        inputs.emplace_back(auxillaryInput);
        LArDLHelper::TorchOutput output;

        LArDLHelper::Forward(trackFeatures.GetExitingStatus() ? m_exitingModel : m_containedModel, inputs, output);

        this->SaveResultsToPfoMetadata(pPfo, output);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::SaveResultsToPfoMetadata(const ParticleFlowObject *pPfo, const LArDLHelper::TorchOutput &scores)
{
    // Due to the way PyTorch includes the softmax in the loss function we need to apply it here for PID scores
    const auto results = torch::softmax(scores, 1);
    const auto pidAccessor = results.accessor<float, 2>();
    const std::vector<std::string> pidClassNames{"muon", "pion", "proton", "kaon"};

    object_creation::ParticleFlowObject::Metadata newPfoMetadata;

    for (unsigned int p = 0; p < pidClassNames.size(); ++p)
    {
        const std::string metadataName{pidClassNames.at(p) + "_pid_score"};
        newPfoMetadata.m_propertiesToAdd[metadataName] = pidAccessor[0][p];
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, newPfoMetadata));
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DlTrackCharacterisationAlgorithm::GetAllTrackFeatures(PfoToTrackFeaturesMap &pfoToTrackFeaturesMap) const
{
    const PfoList *pTrackPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_trackPfoListName, pTrackPfoList));

    if (pTrackPfoList == nullptr || pTrackPfoList->empty())
        return STATUS_CODE_FAILURE;

    for (unsigned int p = 0; p < pTrackPfoList->size(); ++p)
    {
        PfoList::const_iterator iter = pTrackPfoList->begin();
        std::advance(iter, p);

        CaloHitList trackHits;
        LArPfoHelper::GetAllCaloHits2D(*iter, trackHits);
        if (trackHits.empty())
            continue;

        TrackFeatures trackFeatures = this->GetTrackFeatures(trackHits);
        this->GetTrackAuxillaryInfo(*iter, trackFeatures);

        // Use 3D hits for some auxillary information
        CaloHitList trackHits3D;
        LArPfoHelper::GetCaloHits(*iter, TPC_3D, trackHits3D);
        trackFeatures.SetExitingStatus(this->IsTrackExiting(trackHits3D));

        // Truth information
        if (m_trainingMode)
            trackFeatures.SetTruePdg(LArMCParticleHelper::GetMainMCParticle(*iter)->GetParticleId());

        // Deal with missing views
        trackFeatures.AddFeaturesToMissingViews();

        pfoToTrackFeaturesMap[*iter] = trackFeatures;
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlTrackCharacterisationAlgorithm::TrackFeatures DlTrackCharacterisationAlgorithm::GetTrackFeatures(const CaloHitList &caloHits) const
{
    TrackFeatures trackFeatures;
    // Add the hit information
    for (const CaloHit *const pCaloHit : caloHits)
        trackFeatures.AddHit(pCaloHit, pCaloHit->GetHitType(), m_maxChargeThreshold);

    return trackFeatures;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::GetTrackAuxillaryInfo(const ParticleFlowObject *pPfo, TrackFeatures &trackFeatures) const
{
    // Add the auxillary information
    const PfoList &childPfos{pPfo->GetDaughterPfoList()};
    unsigned int nTrackChildren{0}, nShowerChildren{0};
    for (const ParticleFlowObject *const pChildPfo : childPfos)
    {
        nTrackChildren += static_cast<float>(LArPfoHelper::IsTrack(pChildPfo));
        nShowerChildren += static_cast<float>(LArPfoHelper::IsShower(pChildPfo));
    }

    // Normalise the auxillary features
    trackFeatures.AddAuxillaryValue(nTrackChildren < 5 ? static_cast<float>(nTrackChildren) / 5.f : 1.f);
    trackFeatures.AddAuxillaryValue(nShowerChildren < 5 ? static_cast<float>(nShowerChildren) / 5.f : 1.f);

    // Descendents information
    PfoList allDescendantPfos;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, allDescendantPfos);
    // Pfo list includes the pfo in question so we need to subtract one to get the number of descendants
    const float nDescendants{(allDescendantPfos.size() - 1) < 5 ? static_cast<float>(allDescendantPfos.size() - 1) / 5.f : 1.f};
    trackFeatures.AddAuxillaryValue(nDescendants);

    unsigned int nDescendantHits{0};
    for (auto const *pDescendantPfo : allDescendantPfos)
    {
        if (pDescendantPfo == pPfo)
            continue;
        CaloHitList descendantHits;
        LArPfoHelper::GetAllCaloHits2D(pDescendantPfo, descendantHits);
        nDescendantHits += descendantHits.size();
    }
    trackFeatures.AddAuxillaryValue(nDescendantHits);

    const int recoTier{LArPfoHelper::GetHierarchyTier(pPfo)};
    trackFeatures.AddAuxillaryValue(recoTier < 5 ? static_cast<float>(recoTier) / 5.f : 1.f);

    trackFeatures.AddAuxillaryValue(trackFeatures.GetNHits(TPC_VIEW_U));
    trackFeatures.AddAuxillaryValue(trackFeatures.GetNHits(TPC_VIEW_V));
    trackFeatures.AddAuxillaryValue(trackFeatures.GetNHits(TPC_VIEW_W));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlTrackCharacterisationAlgorithm::IsTrackExiting(const CaloHitList &caloHits) const
{
    float xMin{std::numeric_limits<float>::max()}, xMax{std::numeric_limits<float>::lowest()};
    float yMin{std::numeric_limits<float>::max()}, yMax{std::numeric_limits<float>::lowest()};
    float zMin{std::numeric_limits<float>::max()}, zMax{std::numeric_limits<float>::lowest()};
    for (const CaloHit *pCaloHit : caloHits)
    {
        const CartesianVector &pos = pCaloHit->GetPositionVector();
        xMin = std::min(xMin, pos.GetX());
        xMax = std::max(xMax, pos.GetX());
        yMin = std::min(yMin, pos.GetY());
        yMax = std::max(yMax, pos.GetY());
        zMin = std::min(zMin, pos.GetZ());
        zMax = std::max(zMax, pos.GetZ());
    }

    if (xMin < m_detLimitMinX + m_detEdgeTolerance || xMax > m_detLimitMaxX - m_detEdgeTolerance)
        return true;
    if (yMin < m_detLimitMinY + m_detEdgeTolerance || yMax > m_detLimitMaxY - m_detEdgeTolerance)
        return true;
    if (zMin < m_detLimitMinZ + m_detEdgeTolerance || zMax > m_detLimitMaxZ - m_detEdgeTolerance)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::CreateSequenceInput(const TrackHitFeatures &trackHitFeatures, LArDLHelper::TorchInput &sequenceInput) const
{
    constexpr unsigned int nFeatures{3};
    LArDLHelper::InitialiseInput({1, m_sequenceLength, nFeatures}, sequenceInput);
    auto accessor = sequenceInput.accessor<float, 3>();

    // Do we need to truncate? If number of hits > m_sequenceLength we need ignore the first ones
    const unsigned int nHits(trackHitFeatures.size());
    const unsigned int firstHit{nHits <= m_sequenceLength ? 0 : nHits - m_sequenceLength};
    const unsigned int firstSeq{nHits <= m_sequenceLength ? m_sequenceLength - nHits : 0};

    if (nHits > 0)
    {
        for (unsigned int h = 0; h < m_sequenceLength; ++h)
        {
            if (h >= nHits)
                break;

            const std::vector<float> &hit{trackHitFeatures.at(firstHit + h)};
            for (unsigned int f = 0; f < nFeatures; ++f)
            {
                // Calculate positions relative to the first hit and normalise
                if (f < 2)
                    accessor[0][firstSeq + h][f] = (hit.at(f) - trackHitFeatures.at(0).at(f)) / 1000.f;
                else
                    accessor[0][firstSeq + h][f] = hit.at(f);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::CreateAuxillaryInput(const TrackFeatures &trackFeatures, LArDLHelper::TorchInput &auxillaryInput) const
{
    const std::vector<float> &auxFeatures{trackFeatures.GetAuxillaryFeatures()};
    LArDLHelper::InitialiseInput({1, static_cast<int>(auxFeatures.size())}, auxillaryInput);
    auto accessor = auxillaryInput.accessor<float, 2>();

    // We need to do some normalisation here for number of hits variables
    for (unsigned int f = 0; f < auxFeatures.size(); ++f)
    {
        if (f == static_cast<unsigned int>(AuxInputs::kNDescendantHits) || f == static_cast<unsigned int>(AuxInputs::kNHitsU) ||
            f == static_cast<unsigned int>(AuxInputs::kNHitsV) || f == static_cast<unsigned int>(AuxInputs::kNHitsW))
            accessor[0][f] = auxFeatures.at(f) > 0 ? std::log10(auxFeatures.at(f)) / 5.f : 0.f;
        else
            accessor[0][f] = auxFeatures.at(f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::GetDetectorLimits()
{
    const LArTPCMap &larTPCMap = this->GetPandora().GetGeometry()->GetLArTPCMap();
    for (auto const &[_, pTPC] : larTPCMap)
    {
        m_detLimitMinX = std::min(m_detLimitMinX, pTPC->GetCenterX() - 0.5f * pTPC->GetWidthX());
        m_detLimitMaxX = std::max(m_detLimitMaxX, pTPC->GetCenterX() + 0.5f * pTPC->GetWidthX());
        m_detLimitMinY = std::min(m_detLimitMinY, pTPC->GetCenterY() - 0.5f * pTPC->GetWidthY());
        m_detLimitMaxY = std::max(m_detLimitMaxY, pTPC->GetCenterY() + 0.5f * pTPC->GetWidthY());
        m_detLimitMinZ = std::min(m_detLimitMinZ, pTPC->GetCenterZ() - 0.5f * pTPC->GetWidthZ());
        m_detLimitMaxZ = std::max(m_detLimitMaxZ, pTPC->GetCenterZ() + 0.5f * pTPC->GetWidthZ());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlTrackCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinTrackHits", m_minTrackHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SequenceLength", m_sequenceLength));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ExitingModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_exitingModel);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ContainedModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_containedModel);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxChargeThreshold", m_maxChargeThreshold));

    if (m_maxChargeThreshold <= 0.f)
    {
        std::cerr << "DlTrackCharacterisationAlgorithm::ReadSettings: MaxChargeThreshold needs to be positive" << std::endl;
        return STATUS_CODE_FAILURE;
    }
    else
        return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::TrackFeatures::AddHit(const pandora::CaloHit *const pCaloHit, const HitType view, const float chargeThreshold = 1e9)
{
    const CartesianVector &pos{pCaloHit->GetPositionVector()};

    // Truncate (if required) and normalise the charge
    const float q{pCaloHit->GetMipEquivalentEnergy() <= chargeThreshold ? pCaloHit->GetMipEquivalentEnergy() / chargeThreshold : 1.f};

    m_viewToTrackHitFeatures[view].emplace_back(std::vector<float>({pos.GetX(), pos.GetZ(), q}));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::TrackFeatures::AddAuxillaryValue(const float value)
{
    m_auxillaryFeatures.emplace_back(value);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const DlTrackCharacterisationAlgorithm::ViewToTrackHitFeaturesMap &DlTrackCharacterisationAlgorithm::TrackFeatures::GetViewToTrackHitFeaturesMap() const
{
    return m_viewToTrackHitFeatures;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const DlTrackCharacterisationAlgorithm::TrackHitFeatures &DlTrackCharacterisationAlgorithm::TrackFeatures::GetAllHitFeatures(const HitType view) const
{
    return m_viewToTrackHitFeatures.at(view);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::vector<float> &DlTrackCharacterisationAlgorithm::TrackFeatures::GetHit(const unsigned int h, const HitType view) const
{
    return m_viewToTrackHitFeatures.at(view).at(h);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DlTrackCharacterisationAlgorithm::TrackFeatures::GetNHits(const HitType view) const
{
    return m_viewToTrackHitFeatures.count(view) ? m_viewToTrackHitFeatures.at(view).size() : 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::vector<float> &DlTrackCharacterisationAlgorithm::TrackFeatures::GetAuxillaryFeatures() const
{
    return m_auxillaryFeatures;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DlTrackCharacterisationAlgorithm::TrackFeatures::GetNAuxillaryFeatures() const
{
    return m_auxillaryFeatures.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::TrackFeatures::SetTruePdg(const int pdg)
{
    m_truePdg = pdg;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DlTrackCharacterisationAlgorithm::TrackFeatures::GetTruePdg() const
{
    return m_truePdg;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::TrackFeatures::SetExitingStatus(bool isExiting)
{
    m_isExiting = isExiting;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlTrackCharacterisationAlgorithm::TrackFeatures::GetExitingStatus() const
{
    return m_isExiting;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlTrackCharacterisationAlgorithm::TrackFeatures::AddFeaturesToMissingViews()
{
    if (!m_viewToTrackHitFeatures.count(TPC_VIEW_U))
        m_viewToTrackHitFeatures[TPC_VIEW_U] = TrackHitFeatures();
    if (!m_viewToTrackHitFeatures.count(TPC_VIEW_V))
        m_viewToTrackHitFeatures[TPC_VIEW_V] = TrackHitFeatures();
    if (!m_viewToTrackHitFeatures.count(TPC_VIEW_W))
        m_viewToTrackHitFeatures[TPC_VIEW_W] = TrackHitFeatures();
}

} // namespace lar_dl_content
