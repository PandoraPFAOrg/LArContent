/**
 *  @file   src/BasicVertexingMetrics.cc
 *
 *  @brief  Implementation of the basic vertexing metrics.
 *
 *  $Log: $
 */

#include <cmath>
#include <limits>
#include <map>
#include <string>

#include "Objects/CartesianVector.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArSlicing/BasicVertexingMetrics.h"

#define DEBUG_MODE 0
#if DEBUG_MODE
#define HEP_EVD_PANDORA_HELPERS 1
#include "hep_evd.h"
#endif

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

BasicVertexingMetrics::BasicVertexingMetrics() :
    m_caloHitListName(""),
    m_inputVertexListName(""),
    m_treeName{"threeDVertexMetrics"},
    m_fileName{"threeDVertexMetrics.root"}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

BasicVertexingMetrics::~BasicVertexingMetrics()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "BasicVertexingMetrics: Unable to write to ROOT tree" << std::endl;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode BasicVertexingMetrics::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    auto statusCode(PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (statusCode != STATUS_CODE_SUCCESS || nullptr == pCaloHitList || pCaloHitList->empty())
    {
        std::cout << "BasicVertexingMetrics::WriteMetrics - No CaloHits in this slice." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    const VertexList *pInputVertexList{nullptr};
    statusCode = PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList);
    bool noVertexFound(nullptr == pInputVertexList || pInputVertexList->empty());

    // INFO: Allow for the possibility of no vertices in this slice, as this is a valid case.
    if (noVertexFound)
    {
        std::cout << "BasicVertexingMetrics::WriteMetrics - No vertices in this slice." << std::endl;

        // INFO: Add a singular fake entry, to ensure that the tree is filled for this slice, even if no vertices were found.
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position =
            CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
        pInputVertexList = new VertexList({pVertex});
    }

    std::map<const MCParticle *, CaloHitList> neutrinoToHitMap;
    std::map<const CaloHit *, const MCParticle *> hitToNeutrinoMap;
    unsigned int failureCount{0};

    // First, get a map of every hit to its neutrino, and the other way around.
    for (const auto pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);

            if (!pMCParticle)
                throw std::runtime_error("CaloHit with no associated MCParticle");

            const MCParticle *const pNeutrino(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

            if (nullptr == pNeutrino || !LArMCParticleHelper::IsNeutrino(pNeutrino))
                continue;

            neutrinoToHitMap[pNeutrino].push_back(pCaloHit);
            hitToNeutrinoMap[pCaloHit] = pNeutrino;
        }
        catch (StatusCodeException &e)
        {
            ++failureCount;
            continue;
        }
    }

    std::cout << "BasicVertexingMetrics::WriteMetrics - Found " << neutrinoToHitMap.size() << " neutrinos in this slice." << std::endl;
    std::cout << "BasicVertexingMetrics::WriteMetrics - Failed to find neutrino for " << failureCount << " hits, out of "
              << pCaloHitList->size() << std::endl;

    unsigned int mostHits{0};
    float mostEnergetic{0.f};

    for (const auto &neutrinoIt : neutrinoToHitMap)
    {
        const MCParticle *const pNeutrino(neutrinoIt.first);
        const CaloHitList allRecoHitsForNeutrino(neutrinoIt.second);

        if (allRecoHitsForNeutrino.size() > mostHits)
            mostHits = allRecoHitsForNeutrino.size();

        float totalEnergy{0.f};
        for (const auto &hit : allRecoHitsForNeutrino)
            totalEnergy += hit->GetInputEnergy();

        if (totalEnergy > mostEnergetic)
            mostEnergetic = totalEnergy;
    }

#if DEBUG_MODE
    // DEBUG: Add visualization of the semantic labels to EVD.
    HepEVD::Hits hitsToVis;

    int evdHitIdx{0};
    for (const auto pCaloHit : *pCaloHitList)
    {
        if (nullptr == pCaloHit)
            continue;

        const auto x = pCaloHit->GetPositionVector().GetX();
        const auto y = pCaloHit->GetPositionVector().GetY();
        const auto z = pCaloHit->GetPositionVector().GetZ();
        const auto e = pCaloHit->GetInputEnergy();

        HepEVD::Hit *evdHit = new HepEVD::Hit({x, y, z}, e);

        hitsToVis.push_back(evdHit);

        evdHitIdx++;
    }

    HepEVD::setHepEVDGeometry(this->GetPandora().GetGeometry());
    HepEVD::getServer()->addHits(hitsToVis);

    HepEVD::Markers pointsToVis;
#endif

    // With that, we can start to write over the found vertices, and compare them to the MC truth.
    for (const auto &foundVertex : *pInputVertexList)
    {
        const CartesianVector &foundVertexPos = foundVertex->GetPosition();

        float closestDistance{std::numeric_limits<float>::max()};
        for (const auto &neutrinoIt : neutrinoToHitMap)
        {
            const MCParticle *const pNeutrino(neutrinoIt.first);
            const auto &trueVertex(pNeutrino->GetVertex());
            const float distSq = foundVertexPos.GetDistanceSquared(trueVertex);
            const float dist = std::sqrt(distSq);

            if (dist < closestDistance)
                closestDistance = dist;
        }

#if DEBUG_MODE
        HepEVD::Point *recoVtx = new HepEVD::Point({foundVertexPos.GetX(), foundVertexPos.GetY(), foundVertexPos.GetZ()});
        recoVtx->setColour("red");
        pointsToVis.push_back(*recoVtx);
#endif

        for (const auto &neutrinoIt : neutrinoToHitMap)
        {
            const MCParticle *const pNeutrino(neutrinoIt.first);
            const CaloHitList allRecoHitsForNeutrino(neutrinoIt.second);

            const bool hasMostHits(allRecoHitsForNeutrino.size() >= mostHits);

            float totalEnergy{0.f};
            for (const auto &hit : allRecoHitsForNeutrino)
                totalEnergy += hit->GetInputEnergy();

            const bool hasMostEnergy(totalEnergy >= mostEnergetic);

            const float dist = std::sqrt(foundVertexPos.GetDistanceSquared(pNeutrino->GetVertex()));
            const bool isClosestVertex(dist <= closestDistance);

            const auto parentInstance = MultiPandoraApi::GetPrimaryPandoraInstance(&(*this).GetPandora());

#if DEBUG_MODE
            HepEVD::Point *mcVtx = new HepEVD::Point({pNeutrino->GetVertex().GetX(), pNeutrino->GetVertex().GetY(), pNeutrino->GetVertex().GetZ()});

            if (hasMostHits)
                mcVtx->setColour("blue");
            else if (hasMostEnergy)
                mcVtx->setColour("cyan");
            else if (isClosestVertex)
                mcVtx->setColour("magenta");
            else
                mcVtx->setColour("green");

            pointsToVis.push_back(*mcVtx);
#endif

            // Now, we can write out the reco and MC truth values for this vertex.
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "run", static_cast<int>(parentInstance->GetRun())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "subrun", static_cast<int>(parentInstance->GetSubrun())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "event", static_cast<int>(parentInstance->GetEvent())));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuEnergy", pNeutrino->GetEnergy()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuPDG", pNeutrino->GetParticleId()));
            PANDORA_MONITORING_API(
                SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNumNusInSlice", static_cast<int>(neutrinoToHitMap.size())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxX", pNeutrino->GetVertex().GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxY", pNeutrino->GetVertex().GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxZ", pNeutrino->GetVertex().GetZ()));
            PANDORA_MONITORING_API(
                SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNumHitsInSlice", static_cast<int>(allRecoHitsForNeutrino.size())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "biggestContributor", static_cast<int>(hasMostHits)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mostEnergetic", static_cast<int>(hasMostEnergy)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "closestVertex", static_cast<int>(isClosestVertex)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "noVertexInSlice", static_cast<int>(noVertexFound)));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVtxX", foundVertexPos.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVtxY", foundVertexPos.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVtxZ", foundVertexPos.GetZ()));
            PANDORA_MONITORING_API(
                SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNumSliceNumHits", static_cast<int>(pCaloHitList->size())));

            PANDORA_MONITORING_API(
                SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDx", foundVertexPos.GetX() - pNeutrino->GetVertex().GetX()));
            PANDORA_MONITORING_API(
                SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDy", foundVertexPos.GetY() - pNeutrino->GetVertex().GetY()));
            PANDORA_MONITORING_API(
                SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDz", foundVertexPos.GetZ() - pNeutrino->GetVertex().GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDr", dist));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
    }

#if DEBUG_MODE
    HepEVD::getServer()->addMarkers(pointsToVis);
    HepEVD::saveState("FoundVertices");
    HepEVD::startServer();
#endif

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode BasicVertexingMetrics::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
