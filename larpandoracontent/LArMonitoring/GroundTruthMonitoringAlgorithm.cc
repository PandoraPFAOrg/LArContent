/**
 *  @file   larpandoracontent/LArMonitoring/GroundTruthMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArMonitoring/GroundTruthMonitoringAlgorithm.h"
#include "larpandoracontent/LArUtility/RollUp.h"

using namespace pandora;

namespace lar_content
{

GroundTruthMonitoringAlgorithm::GroundTruthMonitoringAlgorithm() :
    m_caloHitListName("CaloHitList2D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GroundTruthMonitoringAlgorithm::Run()
{
    LArMCParticleHelper::MCContributionMap uMCHitsMap, vMCHitsMap, wMCHitsMap;
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    const std::map<HitType, float> lengthThresholds{
        {TPC_VIEW_U, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)},
        {TPC_VIEW_V, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)},
        {TPC_VIEW_W, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)}};
    m_rollUp = RollUpper(std::make_unique<RollUpEMAndAmbiguousDeltaRayHitsPolicy>(0, lengthThresholds));

    CaloHitList uHits, vHits, wHits;
    this->PartitionViews(*pCaloHitList, uHits, vHits, wHits);
    this->MakeMCToHitsMap(uHits, uMCHitsMap);
    this->MakeMCToHitsMap(vHits, vMCHitsMap);
    this->MakeMCToHitsMap(wHits, wMCHitsMap);

    for (const auto &mcMap : {uMCHitsMap, vMCHitsMap, wMCHitsMap})
       this->Visualize(mcMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GroundTruthMonitoringAlgorithm::MakeMCToHitsMap(const CaloHitList &caloHitList, LArMCParticleHelper::MCContributionMap &mcMap) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pMC{m_rollUp.RollUpCaloHit(pCaloHit)};
            mcMap[pMC].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GroundTruthMonitoringAlgorithm::PartitionViews(const CaloHitList &hits2D, CaloHitList &uHits, CaloHitList &vHits, CaloHitList &wHits) const
{
    for (const CaloHit *pCaloHit : hits2D)
    {
        const HitType view{pCaloHit->GetHitType()};
        switch (view)
        {
            case HitType::TPC_VIEW_U:
                uHits.emplace_back(pCaloHit);
                break;
            case HitType::TPC_VIEW_V:
                vHits.emplace_back(pCaloHit);
                break;
            case HitType::TPC_VIEW_W:
                wHits.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GroundTruthMonitoringAlgorithm::Visualize(const LArMCParticleHelper::MCContributionMap &mcMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    std::map<HitType, std::string> views{{TPC_VIEW_U, "U"}, {TPC_VIEW_V, "V"}, {TPC_VIEW_W, "W"}};

    // Organise by tier
    std::map<int, MCParticleList> tierToMCMap;
    for (const auto &[pMC, caloHits] : mcMap)
    {
        const int tier{LArMCParticleHelper::GetHierarchyTier(pMC)};
        tierToMCMap[tier].emplace_back(pMC);
    }

    // Display hits
    for (const auto &[tier, mcParticles] : tierToMCMap)
    {
        for (const auto &pMC : mcParticles)
        {
            const CaloHitList caloHits{mcMap.at(pMC)};
            const HitType view{caloHits.front()->GetHitType()};
            const int pdg{pMC->GetParticleId()};
            std::string description{views.at(view) + std::to_string(tier) + " - " + std::to_string(pdg)};
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, description, AUTOITER));
        }
    }

    // Display true neutrino vertex
    for (const auto &[pMC, caloHits] : mcMap)
    {
        const MCParticle *pParentMC{LArMCParticleHelper::GetParentMCParticle(pMC)};
        const int parentPdg{pParentMC ? std::abs(pParentMC->GetParticleId()) : 0};
        if (parentPdg == 12 || parentPdg == 14 || parentPdg == 16)
        {
            const HitType view{caloHits.front()->GetHitType()};
            const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
            float x{0.f}, u{0.f}, v{0.f}, w{0.f};
            LArVertexHelper::GetTrueVertexPosition(pParentMC->GetVertex(), transform, x, u, v, w);
            const CartesianVector nuVertex(x, 0.f, view == TPC_VIEW_U ? u : view == TPC_VIEW_V ? v : w);
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &nuVertex, "nu vtx", BLACK, 2));
            break;
        }
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GroundTruthMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
