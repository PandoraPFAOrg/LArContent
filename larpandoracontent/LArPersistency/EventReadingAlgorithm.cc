/**
 *  @file   PandoraSDK/src/Persistency/EventReadingAlgorithm.cc
 *
 *  @brief  Implementation of the event reading algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Persistency/BinaryFileReader.h"
#include "Persistency/XmlFileReader.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPersistency/EventReadingAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

EventReadingAlgorithm::EventReadingAlgorithm() :
    m_skipToEvent(0),
    m_isEnhancedEventFile(false),
    m_useLArCaloHits(true),
    m_larCaloHitVersion(1),
    m_useLArMCParticles(true),
    m_larMCParticleVersion(2),
    m_pEventFileReader(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventReadingAlgorithm::~EventReadingAlgorithm()
{
    delete m_pEventFileReader;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::Initialize()
{
    if (!m_geometryFileName.empty())
    {
        const FileType geometryFileType(this->GetFileType(m_geometryFileName));

        if (BINARY == geometryFileType)
        {
            BinaryFileReader fileReader(this->GetPandora(), m_geometryFileName);
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, fileReader.ReadGeometry());
        }
        else if (XML == geometryFileType)
        {
            XmlFileReader fileReader(this->GetPandora(), m_geometryFileName);
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, fileReader.ReadGeometry());
        }
        else
        {
            return STATUS_CODE_FAILURE;
        }
    }

    if (!m_eventFileName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReplaceEventFileReader(m_eventFileName));

        if (m_isEnhancedEventFile)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->ReadGlobalHeader());

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->GoToEvent(m_skipToEvent));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::Run()
{
    if ((nullptr != m_pEventFileReader) && !m_eventFileName.empty())
    {
        try
        {
            m_pEventFileReader->ReadEvent();
        }
        catch (const StatusCodeException &)
        {
            this->MoveToNextEventFile();
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RepeatEventPreparation(*this));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventReadingAlgorithm::MoveToNextEventFile()
{
    if (m_eventFileNameVector.empty())
        throw StopProcessingException("All event files processed");

    m_eventFileName = m_eventFileNameVector.back();
    m_eventFileNameVector.pop_back();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReplaceEventFileReader(m_eventFileName));

    if (m_isEnhancedEventFile)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->ReadGlobalHeader());

    try
    {
        m_pEventFileReader->ReadEvent();
    }
    catch (const StatusCodeException &)
    {
        this->MoveToNextEventFile();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::ReplaceEventFileReader(const std::string &fileName)
{
    delete m_pEventFileReader;
    m_pEventFileReader = nullptr;

    std::cout << "EventReadingAlgorithm: Processing event file: " << fileName << std::endl;
    const FileType eventFileType(this->GetFileType(fileName));

    if (BINARY == eventFileType)
    {
        m_pEventFileReader = new BinaryFileReader(this->GetPandora(), fileName);
    }
    else if (XML == eventFileType)
    {
        m_pEventFileReader = new XmlFileReader(this->GetPandora(), fileName);
    }
    else
    {
        return STATUS_CODE_FAILURE;
    }

    if (m_useLArCaloHits)
        m_pEventFileReader->SetFactory(new LArHitFactory(m_larCaloHitVersion));

    if (m_useLArMCParticles)
        m_pEventFileReader->SetFactory(new LArMCParticleFactory(m_larMCParticleVersion));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FileType EventReadingAlgorithm::GetFileType(const std::string &fileName) const
{
    std::string fileExtension(fileName.substr(fileName.find_last_of(".")));
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

    if (std::string(".xml") == fileExtension)
    {
        return XML;
    }
    else if (std::string(".pndr") == fileExtension)
    {
        return BINARY;
    }
    else
    {
        std::cout << "EventReadingAlgorithm: Unknown file type specified " << fileName << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    ExternalEventReadingParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        pExternalParameters = dynamic_cast<ExternalEventReadingParameters *>(this->GetExternalParameters());

        if (!pExternalParameters)
            return STATUS_CODE_FAILURE;
    }

    if (pExternalParameters && !pExternalParameters->m_geometryFileName.empty())
    {
        m_geometryFileName = pExternalParameters->m_geometryFileName;
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GeometryFileName", m_geometryFileName));
    }

    if (pExternalParameters && !pExternalParameters->m_eventFileNameList.empty())
    {
        XmlHelper::TokenizeString(pExternalParameters->m_eventFileNameList, m_eventFileNameVector, ":");
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadVectorOfValues(xmlHandle, "EventFileNameList", m_eventFileNameVector));
    }

    if (!m_eventFileNameVector.empty())
    {
        std::reverse(m_eventFileNameVector.begin(), m_eventFileNameVector.end());
        m_eventFileName = m_eventFileNameVector.back();
        m_eventFileNameVector.pop_back();
    }

    if (pExternalParameters && pExternalParameters->m_skipToEvent.IsInitialized())
    {
        m_skipToEvent = pExternalParameters->m_skipToEvent.Get();
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SkipToEvent", m_skipToEvent));
    }

    if (m_geometryFileName.empty() && m_eventFileName.empty())
    {
        std::cout << "EventReadingAlgorithm - nothing to do; neither geometry nor event file specified." << std::endl;
        return STATUS_CODE_NOT_INITIALIZED;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IsEnhancedEventFile", m_isEnhancedEventFile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseLArCaloHits", m_useLArCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LArCaloHitVersion", m_larCaloHitVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LArMCParticleVersion", m_larMCParticleVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseLArMCParticles", m_useLArMCParticles));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
