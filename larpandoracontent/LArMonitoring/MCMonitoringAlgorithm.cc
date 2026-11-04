/**
 *  @file   larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

MCMonitoringAlgorithm::MCMonitoringAlgorithm() :
    m_caloHitListName(""),
    m_mcListName("Input"),
    m_visualise(false),
    m_colourByProcess(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCMonitoringAlgorithm::~MCMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::Run()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

    this->BuildMCHitMap();

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::BuildMCHitMap()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcListName, pMCParticleList));

    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_SUCCESS;

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    m_mcHitsMap.clear();
    for (const MCParticle *const pMC : *pMCParticleList)
        m_mcHitsMap[pMC] = CaloHitList();

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            m_mcHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }
    }

    std::map<MCProcess, Color> processColorMap{
        {MC_PROC_UNKNOWN, GRAY},
        {MC_PROC_PRIMARY, BLACK},
        {MC_PROC_COMPT,GREEN},
        {MC_PROC_PHOT,CYAN},
        {MC_PROC_E_IONI,ORANGE},
        {MC_PROC_CONV,MAGENTA},
        {MC_PROC_MU_IONI,ORANGE},
        {MC_PROC_HAD_IONI,ORANGE},
        {MC_PROC_PI_PLUS_INELASTIC,YELLOW}
    };
    for (const auto &[pMC, hits] : m_mcHitsMap)
    {
        const MCParticleList &parentList{pMC->GetParentList()};
        std::string desc{""};
        if (parentList.size() == 1)
        {
            std::ostringstream oss;
            oss << parentList.front();
            desc += std::to_string(parentList.front()->GetParticleId()) + " (" + oss.str() + ")" + " -> ";
            if (parentList.front()->GetParticleId() == 111)
            {
                std::cout << "MCParticle: PDG " << pMC->GetParticleId() << " with parent PDG 111 has " << hits.size() << " hits." << std::endl;
            }
            else if (std::abs(parentList.front()->GetParticleId()) == MU_MINUS)
            {
                if (std::abs(pMC->GetParticleId()) == E_MINUS)
                {
                    const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
                    const MCProcess process{pLArMC ? pLArMC->GetProcess() : MC_PROC_UNKNOWN};
                    if (process == MC_PROC_DECAY)
                    {
                        std::cout << "Found electron from muon with process: " << process << std::endl;
                        desc += " Decay";
                    }
                }
            }
        }
        desc += std::to_string(pMC->GetParticleId());

        if (!hits.empty())
        {
            if (m_colourByProcess)
            {
                const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
                const MCProcess process{pLArMC ? pLArMC->GetProcess() : MC_PROC_UNKNOWN};
                Color color{processColorMap.count(process) ? processColorMap.at(process) : GRAY};
                if (process == MC_PROC_COMPT)
                {
                    if (!pLArMC->GetDaughterList().empty())
                    {
                        color = RED;
                    }
                }

                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, desc, color));
            }
            else
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, desc, AUTOITER));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ColourByProcess", m_colourByProcess));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCListName", m_mcListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
