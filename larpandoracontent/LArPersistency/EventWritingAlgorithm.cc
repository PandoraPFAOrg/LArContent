/**
 *  @file   PandoraSDK/src/Persistency/EventWritingAlgorithm.cc
 *
 *  @brief  Implementation of the event writing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Persistency/BinaryFileWriter.h"
#include "Persistency/XmlFileWriter.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPersistency/EventWritingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EventWritingAlgorithm::EventWritingAlgorithm() :
    m_geometryFileType(UNKNOWN_FILE_TYPE),
    m_eventFileType(UNKNOWN_FILE_TYPE),
    m_pEventFileWriter(nullptr),
    m_pGeometryFileWriter(nullptr),
    m_shouldWriteGeometry(false),
    m_writtenGeometry(false),
    m_shouldWriteEvents(true),
    m_fileMajorVersion(2),
    m_fileMinorVersion(0),
    m_writtenEventGlobalHeader(false),
    m_shouldWriteMCRelationships(true),
    m_shouldWriteTrackRelationships(true),
    m_shouldOverwriteEventFile(false),
    m_shouldOverwriteGeometryFile(false),
    m_useLArCaloHits(true),
    m_larCaloHitVersion(1),
    m_useLArMCParticles(true),
    m_shouldFilterByNuanceCode(false),
    m_filterNuanceCode(0),
    m_shouldFilterByMCParticles(false),
    m_neutrinoInducedOnly(false),
    m_matchingMinPrimaryHits(15),
    m_nNonNeutrons(0),
    m_nMuons(0),
    m_nElectrons(0),
    m_nProtons(0),
    m_nPhotons(0),
    m_nChargedPions(0),
    m_shouldFilterByNeutrinoVertex(false),
    m_detectorHalfLengthX(-1.f),
    m_detectorHalfLengthY(-1.f),
    m_detectorHalfLengthZ(-1.f),
    m_coordinateOffsetX(0.f),
    m_coordinateOffsetY(0.f),
    m_coordinateOffsetZ(0.f),
    m_selectedBorderX(-1.f),
    m_selectedBorderY(-1.f),
    m_selectedBorderZ(-1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventWritingAlgorithm::~EventWritingAlgorithm()
{
    delete m_pEventFileWriter;
    delete m_pGeometryFileWriter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventWritingAlgorithm::Initialize()
{
    if (m_shouldWriteGeometry)
    {
        const FileMode fileMode(m_shouldOverwriteGeometryFile ? OVERWRITE : APPEND);

        if (BINARY == m_geometryFileType)
        {
            m_pGeometryFileWriter = new BinaryFileWriter(this->GetPandora(), m_geometryFileName, fileMode, m_fileMajorVersion, m_fileMinorVersion);
        }
        else if (XML == m_geometryFileType)
        {
            m_pGeometryFileWriter = new XmlFileWriter(this->GetPandora(), m_geometryFileName, fileMode, m_fileMajorVersion, m_fileMinorVersion);
        }
        else
        {
            return STATUS_CODE_FAILURE;
        }
    }

    if (m_shouldWriteEvents)
    {
        const FileMode fileMode(m_shouldOverwriteEventFile ? OVERWRITE : APPEND);

        if (BINARY == m_eventFileType)
        {
            m_pEventFileWriter = new BinaryFileWriter(this->GetPandora(), m_eventFileName, fileMode, m_fileMajorVersion, m_fileMinorVersion);
        }
        else if (XML == m_eventFileType)
        {
            m_pEventFileWriter = new XmlFileWriter(this->GetPandora(), m_eventFileName, fileMode, m_fileMajorVersion, m_fileMinorVersion);
        }
        else
        {
            return STATUS_CODE_FAILURE;
        }

        if (m_useLArCaloHits)
            m_pEventFileWriter->SetFactory(new LArHitFactory(m_larCaloHitVersion));

        if (m_useLArMCParticles)
            m_pEventFileWriter->SetFactory(new LArMCParticleFactory);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventWritingAlgorithm::Run()
{
    // ATTN Should complete geometry creation in LArSoft begin job, but some channel status service functionality unavailable at that point
    if (!m_writtenGeometry && m_pGeometryFileWriter && m_shouldWriteGeometry)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pGeometryFileWriter->WriteGeometry());
        m_writtenGeometry = true;
    }

    bool matchNuanceCode(!m_shouldFilterByNuanceCode || this->PassNuanceCodeFilter());
    bool matchParticles(!m_shouldFilterByMCParticles || this->PassMCParticleFilter());
    bool matchNeutrinoVertexPosition(!m_shouldFilterByNeutrinoVertex || this->PassNeutrinoVertexFilter());

    if (matchNuanceCode && matchParticles && matchNeutrinoVertexPosition && m_pEventFileWriter && m_shouldWriteEvents)
    {
        const CaloHitList *pCaloHitList = nullptr;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

        const TrackList *pTrackList = nullptr;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pTrackList));

        const MCParticleList *pMCParticleList = nullptr;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

        if (!m_writtenEventGlobalHeader)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileWriter->WriteGlobalHeader());
            m_writtenEventGlobalHeader = true;
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            m_pEventFileWriter->WriteEvent(*pCaloHitList, *pTrackList, *pMCParticleList, m_shouldWriteMCRelationships, m_shouldWriteTrackRelationships));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventWritingAlgorithm::PassNuanceCodeFilter() const
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    MCParticleVector mcNeutrinoList;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoList);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle *>(pMCNeutrino);

        if (pLArMCNeutrino && (pLArMCNeutrino->GetNuanceCode() == m_filterNuanceCode))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventWritingAlgorithm::PassMCParticleFilter() const
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::MCContributionMap mcParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcParticlesToGoodHitsMap);

    if (!m_neutrinoInducedOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, mcParticlesToGoodHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, mcParticlesToGoodHitsMap);
    }

    unsigned int nNonNeutrons(0), nMuons(0), nElectrons(0), nProtons(0), nPhotons(0), nChargedPions(0);

    MCParticleVector mcPrimaryVector;
    for (const auto &mapEntry : mcParticlesToGoodHitsMap)
        mcPrimaryVector.push_back(mapEntry.first);
    std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const unsigned int particleId(std::abs(pMCPrimary->GetParticleId()));
        if (NEUTRON != particleId)
            ++nNonNeutrons;

        if (MU_MINUS == particleId)
            ++nMuons;
        else if (E_MINUS == particleId)
            ++nElectrons;
        else if (PROTON == particleId)
            ++nProtons;
        else if (PHOTON == particleId)
            ++nPhotons;
        else if (PI_PLUS == particleId)
            ++nChargedPions;
    }

    if ((nNonNeutrons == m_nNonNeutrons) && (nMuons == m_nMuons) && (nElectrons == m_nElectrons) && (nProtons == m_nProtons) &&
        (nPhotons == m_nPhotons) && (nChargedPions == m_nChargedPions))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventWritingAlgorithm::PassNeutrinoVertexFilter() const
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        const CartesianVector &neutrinoInteractionVertex(pMCNeutrino->GetEndpoint());

        if ((neutrinoInteractionVertex.GetX() < (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) &&
            (neutrinoInteractionVertex.GetX() > (-m_coordinateOffsetX + m_selectedBorderX)) &&
            (neutrinoInteractionVertex.GetY() < (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) &&
            (neutrinoInteractionVertex.GetY() > (-m_coordinateOffsetY + m_selectedBorderY)) &&
            (neutrinoInteractionVertex.GetZ() < (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) &&
            (neutrinoInteractionVertex.GetZ() > (-m_coordinateOffsetZ + m_selectedBorderZ)))
        {
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventWritingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldWriteGeometry", m_shouldWriteGeometry));

    if (m_shouldWriteGeometry)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "GeometryFileName", m_geometryFileName));

        std::string fileExtension(m_geometryFileName.substr(m_geometryFileName.find_last_of(".")));
        std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

        if (std::string(".xml") == fileExtension)
        {
            m_geometryFileType = XML;
        }
        else if (std::string(".pndr") == fileExtension)
        {
            m_geometryFileType = BINARY;
        }
        else
        {
            std::cout << "EventReadingAlgorithm: Unknown geometry file type specified " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldWriteEvents", m_shouldWriteEvents));

    if (m_shouldWriteEvents)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "EventFileName", m_eventFileName));

        std::string fileExtension(m_eventFileName.substr(m_eventFileName.find_last_of(".")));
        std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

        if (std::string(".xml") == fileExtension)
        {
            m_eventFileType = XML;
        }
        else if (std::string(".pndr") == fileExtension)
        {
            m_eventFileType = BINARY;
        }
        else
        {
            std::cout << "EventReadingAlgorithm: Unknown event file type specified " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileMajorVersion", m_fileMajorVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileMinorVersion", m_fileMinorVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldWriteMCRelationships", m_shouldWriteMCRelationships));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldWriteTrackRelationships", m_shouldWriteTrackRelationships));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldOverwriteEventFile", m_shouldOverwriteEventFile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldOverwriteGeometryFile", m_shouldOverwriteGeometryFile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseLArCaloHits", m_useLArCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LArCaloHitVersion", m_larCaloHitVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseLArMCParticles", m_useLArMCParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldFilterByNuanceCode", m_shouldFilterByNuanceCode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldFilterByMCParticles", m_shouldFilterByMCParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldFilterByNeutrinoVertex", m_shouldFilterByNeutrinoVertex));

    if (m_shouldFilterByNuanceCode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FilterNuanceCode", m_filterNuanceCode));
    }

    if (m_shouldFilterByMCParticles)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoInducedOnly", m_neutrinoInducedOnly));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MatchingMinPrimaryHits", m_matchingMinPrimaryHits));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NNonNeutrons", m_nNonNeutrons));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NMuons", m_nMuons));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NElectrons", m_nElectrons));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NProtons", m_nProtons));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NPhotons", m_nPhotons));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NChargedPions", m_nChargedPions));
    }

    if (m_shouldFilterByNeutrinoVertex)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DetectorHalfLengthX", m_detectorHalfLengthX));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DetectorHalfLengthY", m_detectorHalfLengthY));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DetectorHalfLengthZ", m_detectorHalfLengthZ));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CoordinateOffsetX", m_coordinateOffsetX));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CoordinateOffsetY", m_coordinateOffsetY));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CoordinateOffsetZ", m_coordinateOffsetZ));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SelectedBorderX", m_selectedBorderX));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SelectedBorderY", m_selectedBorderY));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SelectedBorderZ", m_selectedBorderZ));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
