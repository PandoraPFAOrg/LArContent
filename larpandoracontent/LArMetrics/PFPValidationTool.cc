/**
 *  @file   larpandoracontent/LArMetrics/PFPValidationTool.cc
 *
 *  @brief  Implementation of the pfp validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArMetrics/PFPValidationTool.h"

using namespace pandora;

namespace lar_content
{

PFPValidationTool::PFPValidationTool() :
    m_pNuVertexList(nullptr),
    m_nuVertexListName("NeutrinoVertices3D"),
    m_michelPDG(777)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PFPValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu,
    const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Get reco neutrino vertex (will handle if not found)
    m_pNuVertexList = nullptr;
    PandoraContentApi::GetList(*pAlgorithm, m_nuVertexListName, m_pNuVertexList);

    PFPTreeVars pfpTreeVars;
    pfpTreeVars.m_run = this->GetPandora().GetRun();
    pfpTreeVars.m_subrun = this->GetPandora().GetSubrun();
    pfpTreeVars.m_event = this->GetPandora().GetEvent();

    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMC(targetMC.at(i));
        const Pfo *const pBestMatch(bestRecoMatch.at(i));

        this->GetMCParticleInfo(pMCNu, pMC, pfpTreeVars);

        if (pBestMatch)
        {
            this->GetRecoParticleInfo(pMC, pBestMatch, pfpTreeVars);
        }
        else
        {
            pfpTreeVars.m_recoVertexX.push_back(m_invalidLargeFloat);
            pfpTreeVars.m_recoVertexY.push_back(m_invalidLargeFloat);
            pfpTreeVars.m_recoVertexZ.push_back(m_invalidLargeFloat);
            pfpTreeVars.m_vertexAcc.push_back(m_invalidLargeFloat);
            pfpTreeVars.m_recoLength.push_back(m_invalidSmallFloat);
            pfpTreeVars.m_recoDisplacement.push_back(m_invalidSmallFloat);
            pfpTreeVars.m_isTrack.push_back(m_invalidInt);
            pfpTreeVars.m_isShower.push_back(m_invalidInt);
        }

        this->GetMatchingInfo(mcMatchesVec, pMC, pBestMatch, pfpTreeVars);
        this->GetAltMatchInfo(mcMatchesVec, targetMC, pMC, pBestMatch, pfpTreeVars);
    }

    this->FillTree(pfpTreeVars);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetMCParticleInfo(const MCParticle *const pMCNu, const MCParticle *const pMCTarget, PFPTreeVars &pfpTreeVars)
{
    // Energy
    pfpTreeVars.m_trueEnergy.push_back(pMCTarget->GetEnergy());
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pMCTarget));
    if (!pLArMCParticle)
    {
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
    pfpTreeVars.m_trueVisEnergy.push_back(pLArMCParticle->GetVisibleEnergy());

    // Vertex
    const CartesianVector &trueVertex((pMCTarget->GetParticleId() == PHOTON) ? pMCTarget->GetEndpoint() : pMCTarget->GetVertex());
    const CartesianVector &trueNuVertex(pMCNu->GetVertex());
    pfpTreeVars.m_trueVertexX.push_back(trueVertex.GetX());
    pfpTreeVars.m_trueVertexY.push_back(trueVertex.GetY());
    pfpTreeVars.m_trueVertexZ.push_back(trueVertex.GetZ());
    pfpTreeVars.m_trueDisplacement.push_back((trueNuVertex - trueVertex).GetMagnitude());

    // Endpoint
    const CartesianVector &trueEnd(pMCTarget->GetEndpoint());
    pfpTreeVars.m_trueEndX.push_back(trueEnd.GetX());
    pfpTreeVars.m_trueEndY.push_back(trueEnd.GetY());
    pfpTreeVars.m_trueEndZ.push_back(trueEnd.GetZ());

    // Length
    pfpTreeVars.m_trueLength.push_back((trueVertex - trueEnd).GetMagnitude());

    // Initial direction
    const CartesianVector &mcMom(pMCTarget->GetMomentum());
    const float momMagSq(mcMom.GetMagnitudeSquared());
    if (momMagSq < std::numeric_limits<float>::epsilon())
    {
        pfpTreeVars.m_trueThetaXZ.push_back(m_invalidAngle);
        pfpTreeVars.m_trueThetaYZ.push_back(m_invalidAngle);
        pfpTreeVars.m_trueDirX.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_trueDirY.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_trueDirZ.push_back(m_invalidLargeFloat);
    }
    else
    {
        pfpTreeVars.m_trueThetaXZ.push_back(atan2(mcMom.GetX(), mcMom.GetZ()));
        pfpTreeVars.m_trueThetaYZ.push_back(asin(mcMom.GetY() / mcMom.GetMagnitude()));
        const CartesianVector trueDir(mcMom.GetUnitVector());
        pfpTreeVars.m_trueDirX.push_back(trueDir.GetX());
        pfpTreeVars.m_trueDirY.push_back(trueDir.GetY());
        pfpTreeVars.m_trueDirZ.push_back(trueDir.GetZ());
    }

    // End direction
    const CartesianVector trueEndDir(pLArMCParticle->GetEndDirection());
    if (momMagSq < std::numeric_limits<float>::epsilon())
    {
        pfpTreeVars.m_trueEndDirX.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_trueEndDirY.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_trueEndDirZ.push_back(m_invalidLargeFloat);
    }
    else
    {
        pfpTreeVars.m_trueEndDirX.push_back(trueEndDir.GetX());
        pfpTreeVars.m_trueEndDirY.push_back(trueEndDir.GetY());
        pfpTreeVars.m_trueEndDirZ.push_back(trueEndDir.GetZ());
    }

    // TruePDG
    int truePDG(pMCTarget->GetParticleId());
    if (abs(pMCTarget->GetParticleId()) == PHOTON) // 111 = photon from pi0
    {
        if (pMCTarget->GetParentList().front()->GetParticleId() == PI_ZERO)
            truePDG = PI_ZERO;
    }
    else if (std::abs(pMCTarget->GetParticleId()) == E_MINUS) // 777 = michel electron
    {
        const MCParticle *const pMCParent(pMCTarget->GetParentList().front());

        if ((abs(pMCParent->GetParticleId()) == MU_MINUS) || (abs(pMCParent->GetParticleId()) == PI_PLUS))
        {
            if (((pMCParent->GetEndpoint() - pMCTarget->GetVertex()).GetMagnitude() < m_maxMichelSep) && (LArMCParticleHelper::IsDecay(pMCTarget)))
            {
                truePDG = m_michelPDG;
            }
        }
    }

    pfpTreeVars.m_truePDG.push_back(truePDG);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetRecoParticleInfo(const MCParticle *const pMCTarget, const Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars)
{
    try
    {
        const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(pBestMatch));
        pfpTreeVars.m_recoVertexX.push_back(pRecoVertex->GetPosition().GetX());
        pfpTreeVars.m_recoVertexY.push_back(pRecoVertex->GetPosition().GetY());
        pfpTreeVars.m_recoVertexZ.push_back(pRecoVertex->GetPosition().GetZ());

        // Signed vertexAcc
        float vertexAcc(m_invalidLargeFloat);
        try
        {
            const CartesianVector &trueVertex((pMCTarget->GetParticleId() == PHOTON) ? pMCTarget->GetEndpoint() : pMCTarget->GetVertex());
            vertexAcc = (pRecoVertex->GetPosition() - trueVertex).GetMagnitude();
            const float sign(
                (vertexAcc < std::numeric_limits<float>::epsilon() || pMCTarget->GetMomentum().GetMagnitude() < std::numeric_limits<float>::epsilon()) ? 1.f
                : (pRecoVertex->GetPosition() - trueVertex).GetOpeningAngle(pMCTarget->GetMomentum()) < (M_PI * 0.5) ? 1.f
                : -1.f);
            vertexAcc *= sign;
        }
        catch (StatusCodeException &) { vertexAcc = m_invalidLargeFloat; }

        pfpTreeVars.m_vertexAcc.push_back(vertexAcc);

        try
        {
            pfpTreeVars.m_recoLength.push_back(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pBestMatch)));
        }
        catch (...)
        {
            pfpTreeVars.m_recoLength.push_back(m_invalidSmallFloat);
        }

        if (m_pNuVertexList && !m_pNuVertexList->empty())
        {
            pfpTreeVars.m_recoDisplacement.push_back((m_pNuVertexList->front()->GetPosition() - pRecoVertex->GetPosition()).GetMagnitude());
        }
        else
        {
            pfpTreeVars.m_recoDisplacement.push_back(m_invalidSmallFloat);
        }
    }
    catch (...)
    {
        pfpTreeVars.m_recoVertexX.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_recoVertexY.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_recoVertexZ.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_vertexAcc.push_back(m_invalidLargeFloat);
        pfpTreeVars.m_recoLength.push_back(m_invalidSmallFloat);
        pfpTreeVars.m_recoDisplacement.push_back(m_invalidSmallFloat);
    }

    // PID
    pfpTreeVars.m_isTrack.push_back(LArPfoHelper::IsTrack(pBestMatch));
    pfpTreeVars.m_isShower.push_back(LArPfoHelper::IsShower(pBestMatch));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetMatchingInfo(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticle *const pMCTarget,
    const Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars)
{
    // Need to loop over matches to find this match :(
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() != pMCTarget)
            continue;

        const CaloHitList &mcHits(mcMatches.GetMC()->GetCaloHits());
        pfpTreeVars.m_nMCHits2D.push_back(mcHits.size());

        for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            // Get vectors to fill
            IntVector &viewNMCHits(hitType == TPC_VIEW_U ? pfpTreeVars.m_nMCHitsU
                    : hitType == TPC_VIEW_V              ? pfpTreeVars.m_nMCHitsV
                                                         : pfpTreeVars.m_nMCHitsW);
            IntVector &viewNPfoHits(hitType == TPC_VIEW_U ? pfpTreeVars.m_nPfoHitsU
                    : hitType == TPC_VIEW_V               ? pfpTreeVars.m_nPfoHitsV
                                                          : pfpTreeVars.m_nPfoHitsW);
            FloatVector &viewCompleteness(hitType == TPC_VIEW_U ? pfpTreeVars.m_completenessU
                    : hitType == TPC_VIEW_V                     ? pfpTreeVars.m_completenessV
                                                                : pfpTreeVars.m_completenessW);
            FloatVector &viewCompletenessADC(hitType == TPC_VIEW_U ? pfpTreeVars.m_completenessADCU
                    : hitType == TPC_VIEW_V                        ? pfpTreeVars.m_completenessADCV
                                                                   : pfpTreeVars.m_completenessADCW);
            FloatVector &viewPurity(hitType == TPC_VIEW_U ? pfpTreeVars.m_purityU
                    : hitType == TPC_VIEW_V               ? pfpTreeVars.m_purityV
                                                          : pfpTreeVars.m_purityW);
            FloatVector &viewPurityADC(hitType == TPC_VIEW_U ? pfpTreeVars.m_purityADCU
                    : hitType == TPC_VIEW_V                  ? pfpTreeVars.m_purityADCV
                                                             : pfpTreeVars.m_purityADCW);

            // Calculate matching vars
            float totalEnergy(0.f); // just for function call
            CaloHitVector viewMCHits;
            this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);
            viewNMCHits.push_back(viewMCHits.size());

            if (pBestMatch)
            {
                CaloHitList viewPfoHits;
                LArPfoHelper::GetCaloHits(pBestMatch, hitType, viewPfoHits);
                viewNPfoHits.push_back(viewPfoHits.size());

                // Find the best match
                const LArHierarchyHelper::RecoHierarchy::Node *pRecoMatchNode(nullptr);
                for (const auto pRecoNode : mcMatches.GetRecoMatches())
                {
                    if (pRecoNode->GetRecoParticles().front() == pBestMatch)
                        pRecoMatchNode = pRecoNode;
                }

                viewCompleteness.push_back(mcMatches.GetCompleteness(pRecoMatchNode, hitType, false));
                viewCompletenessADC.push_back(mcMatches.GetCompleteness(pRecoMatchNode, hitType, true));
                viewPurity.push_back(mcMatches.GetPurity(pRecoMatchNode, hitType, false));
                viewPurityADC.push_back(mcMatches.GetPurity(pRecoMatchNode, hitType, true));
            }
            else
            {
                viewNPfoHits.push_back(0);
                viewCompleteness.push_back(0);
                viewCompletenessADC.push_back(0);
                viewPurity.push_back(0);
                viewPurityADC.push_back(0);
            }
        }

        if (pBestMatch)
        {
            pfpTreeVars.m_hasMatch.push_back(1);
            pfpTreeVars.m_nPfoHits2D.push_back(LArPfoHelper::GetNumberOfTwoDHits(pBestMatch));
            pfpTreeVars.m_nPfoHits3D.push_back(LArPfoHelper::GetNumberOfThreeDHits(pBestMatch));

            // Find the best match
            const LArHierarchyHelper::RecoHierarchy::Node *pRecoMatchNode(nullptr);
            for (const auto pRecoNode : mcMatches.GetRecoMatches())
            {
                if (pRecoNode->GetRecoParticles().front() == pBestMatch)
                    pRecoMatchNode = pRecoNode;
            }

            pfpTreeVars.m_completeness.push_back(mcMatches.GetCompleteness(pRecoMatchNode, false));
            pfpTreeVars.m_completenessADC.push_back(mcMatches.GetCompleteness(pRecoMatchNode, true));
            pfpTreeVars.m_purity.push_back(mcMatches.GetPurity(pRecoMatchNode, false));
            pfpTreeVars.m_purityADC.push_back(mcMatches.GetPurity(pRecoMatchNode, true));
        }
        else
        {
            pfpTreeVars.m_hasMatch.push_back(0);
            pfpTreeVars.m_nPfoHits2D.push_back(0);
            pfpTreeVars.m_nPfoHits3D.push_back(0);
            pfpTreeVars.m_completeness.push_back(0.f);
            pfpTreeVars.m_completenessADC.push_back(0.f);
            pfpTreeVars.m_purity.push_back(0.f);
            pfpTreeVars.m_purityADC.push_back(0.f);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetAltMatchInfo(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticleVector &targetMC,
    const MCParticle *const pMCTarget, const Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars)
{
    // First find MCNode, sorry for looping over matches :(
    const LArHierarchyHelper::MCHierarchy::Node *pMCTargetNode(nullptr);
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() != pMCTarget)
            continue;

        pMCTargetNode = mcMatches.GetMC();
        break;
    }

    if (!pMCTargetNode)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // Find the highest alternative completeness
    float altCompleteness(m_invalidSmallFloat), altPurity(m_invalidSmallFloat);
    int altMCIndex(m_invalidInt);

    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        for (const LArHierarchyHelper::RecoHierarchy::Node *const pRecoNode : mcMatches.GetRecoMatches())
        {
            if (pBestMatch)
            {
                if (pBestMatch == pRecoNode->GetRecoParticles().front())
                    continue;
            }

            float thisCompleteness(m_invalidSmallFloat), thisPurity(m_invalidSmallFloat);
            this->GetAltMetrics(pMCTargetNode, pRecoNode, thisCompleteness, thisPurity);

            if ((thisCompleteness > altCompleteness) && (thisCompleteness > std::numeric_limits<float>::epsilon()))
            {
                const MCParticle *const pThisMC(mcMatches.GetMC()->GetMCParticles().front());
                const auto thisMCIter(std::find(targetMC.begin(), targetMC.end(), pThisMC));

                if (thisMCIter == targetMC.end())
                    continue;

                const int thisMCIndex(thisMCIter - targetMC.begin());

                altCompleteness = thisCompleteness;
                altPurity = thisPurity;
                altMCIndex = thisMCIndex;
            }
        }
    }

    // Tell us about the alt match
    int isUpstreamHierarchy(m_invalidInt);
    int isSameMC(m_invalidInt);
    if (altMCIndex != m_invalidInt)
    {
        // Is upstream hierarchy
        const MCParticle *pSeed(pMCTarget);
        while (pSeed->GetParentList().size() != 0)
        {
            if (pSeed->GetParentList().front() == targetMC.at(altMCIndex))
                isUpstreamHierarchy = 1;

            pSeed = pSeed->GetParentList().front();
        }

        if (isUpstreamHierarchy != 1)
            isUpstreamHierarchy = 0;

        // Is same MC
        isSameMC = (targetMC.at(altMCIndex) == pMCTarget) ? 1 : 0;
    }

    pfpTreeVars.m_altCompleteness.push_back(altCompleteness);
    pfpTreeVars.m_altPurity.push_back(altPurity);
    pfpTreeVars.m_altPDG.push_back(altMCIndex == m_invalidInt ? m_invalidInt : targetMC.at(altMCIndex)->GetParticleId());
    pfpTreeVars.m_altIsUpstreamHierarchy.push_back(isUpstreamHierarchy);
    pfpTreeVars.m_altIsSameMC.push_back(isSameMC);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetAltMetrics(const LArHierarchyHelper::MCHierarchy::Node *const pMCNode,
    const LArHierarchyHelper::RecoHierarchy::Node *const pRecoNode, float &completeness, float &purity)
{
    const CaloHitList &mcHits(pMCNode->GetCaloHits());
    const CaloHitList &recoHits(pRecoNode->GetCaloHits());

    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    completeness = mcHits.size() == 0 ? m_invalidSmallFloat : static_cast<float>(intersection.size()) / mcHits.size();
    purity = recoHits.size() == 0 ? m_invalidSmallFloat : static_cast<float>(intersection.size()) / recoHits.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::FillTree(PFPTreeVars &pfpTreeVars)
{
    IntVector &truePDG = pfpTreeVars.m_truePDG;
    FloatVector &trueEnergy = pfpTreeVars.m_trueEnergy;
    FloatVector &trueVisEnergy = pfpTreeVars.m_trueVisEnergy;
    FloatVector &trueThetaXZ = pfpTreeVars.m_trueThetaXZ;
    FloatVector &trueThetaYZ = pfpTreeVars.m_trueThetaYZ;
    IntVector &isTrack = pfpTreeVars.m_isTrack;
    IntVector &isShower = pfpTreeVars.m_isShower;
    IntVector &hasMatch = pfpTreeVars.m_hasMatch;
    IntVector &nMCHitsU = pfpTreeVars.m_nMCHitsU;
    IntVector &nMCHitsV = pfpTreeVars.m_nMCHitsV;
    IntVector &nMCHitsW = pfpTreeVars.m_nMCHitsW;
    IntVector &nMCHits2D = pfpTreeVars.m_nMCHits2D;
    IntVector &nPfoHitsU = pfpTreeVars.m_nPfoHitsU;
    IntVector &nPfoHitsV = pfpTreeVars.m_nPfoHitsV;
    IntVector &nPfoHitsW = pfpTreeVars.m_nPfoHitsW;
    IntVector &nPfoHits2D = pfpTreeVars.m_nPfoHits2D;
    IntVector &nPfoHits3D = pfpTreeVars.m_nPfoHits3D;
    FloatVector &completeness = pfpTreeVars.m_completeness;
    FloatVector &completenessU = pfpTreeVars.m_completenessU;
    FloatVector &completenessV = pfpTreeVars.m_completenessV;
    FloatVector &completenessW = pfpTreeVars.m_completenessW;
    FloatVector &completenessADC = pfpTreeVars.m_completenessADC;
    FloatVector &completenessADCU = pfpTreeVars.m_completenessADCU;
    FloatVector &completenessADCV = pfpTreeVars.m_completenessADCV;
    FloatVector &completenessADCW = pfpTreeVars.m_completenessADCW;
    FloatVector &purity = pfpTreeVars.m_purity;
    FloatVector &purityU = pfpTreeVars.m_purityU;
    FloatVector &purityV = pfpTreeVars.m_purityV;
    FloatVector &purityW = pfpTreeVars.m_purityW;
    FloatVector &purityADC = pfpTreeVars.m_purityADC;
    FloatVector &purityADCU = pfpTreeVars.m_purityADCU;
    FloatVector &purityADCV = pfpTreeVars.m_purityADCV;
    FloatVector &purityADCW = pfpTreeVars.m_purityADCW;
    FloatVector &altCompleteness = pfpTreeVars.m_altCompleteness;
    FloatVector &altPurity = pfpTreeVars.m_altPurity;
    IntVector &altPDG = pfpTreeVars.m_altPDG;
    IntVector &altIsUpstreamHierarchy = pfpTreeVars.m_altIsUpstreamHierarchy;
    IntVector &altIsSameMC = pfpTreeVars.m_altIsSameMC;
    FloatVector &trueVertexX = pfpTreeVars.m_trueVertexX;
    FloatVector &trueVertexY = pfpTreeVars.m_trueVertexY;
    FloatVector &trueVertexZ = pfpTreeVars.m_trueVertexZ;
    FloatVector &trueEndX = pfpTreeVars.m_trueEndX;
    FloatVector &trueEndY = pfpTreeVars.m_trueEndY;
    FloatVector &trueEndZ = pfpTreeVars.m_trueEndZ;
    FloatVector &trueDirX = pfpTreeVars.m_trueDirX;
    FloatVector &trueDirY = pfpTreeVars.m_trueDirY;
    FloatVector &trueDirZ = pfpTreeVars.m_trueDirZ;
    FloatVector &trueEndDirX = pfpTreeVars.m_trueEndDirX;
    FloatVector &trueEndDirY = pfpTreeVars.m_trueEndDirY;
    FloatVector &trueEndDirZ = pfpTreeVars.m_trueEndDirZ;
    FloatVector &trueLength = pfpTreeVars.m_trueLength;
    FloatVector &trueDisplacement = pfpTreeVars.m_trueDisplacement;
    FloatVector &recoVertexX = pfpTreeVars.m_recoVertexX;
    FloatVector &recoVertexY = pfpTreeVars.m_recoVertexY;
    FloatVector &recoVertexZ = pfpTreeVars.m_recoVertexZ;
    FloatVector &vertexAcc = pfpTreeVars.m_vertexAcc;
    FloatVector &recoLength = pfpTreeVars.m_recoLength;
    FloatVector &recoDisplacement = pfpTreeVars.m_recoDisplacement;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "Run", pfpTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "Subrun", pfpTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "Event", pfpTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TruePDG", &truePDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueEnergy", &trueEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueVisEnergy", &trueVisEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueThetaXZ", &trueThetaXZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueThetaYZ", &trueThetaYZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_IsTrack", &isTrack));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_IsShower", &isShower));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_HasMatch", &hasMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHitsU", &nMCHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHitsV", &nMCHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHitsW", &nMCHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHits2D", &nMCHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHitsU", &nPfoHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHitsV", &nPfoHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHitsW", &nPfoHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHits2D", &nPfoHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHits3D", &nPfoHits3D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Completeness", &completeness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessU", &completenessU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessV", &completenessV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessW", &completenessW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessADC", &completenessADC));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessADCU", &completenessADCU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessADCV", &completenessADCV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessADCW", &completenessADCW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Purity", &purity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityU", &purityU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityV", &purityV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityW", &purityW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityADC", &purityADC));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityADCU", &purityADCU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityADCV", &purityADCV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityADCW", &purityADCW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "ALT_Completeness", &altCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "ALT_Purity", &altPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "ALT_PDG", &altPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "ALT_IsUpstreamHierarchy", &altIsUpstreamHierarchy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "ALT_IsSameMC", &altIsSameMC));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_VertexX", &trueVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_VertexY", &trueVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_VertexZ", &trueVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndX", &trueEndX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndY", &trueEndY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndZ", &trueEndZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_DirX", &trueDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_DirY", &trueDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_DirZ", &trueDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndDirX", &trueEndDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndDirY", &trueEndDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndDirZ", &trueEndDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_Length", &trueLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_Displacement", &trueDisplacement));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexX", &recoVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexY", &recoVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexZ", &recoVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexAcc", &vertexAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Length", &recoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Displacement", &recoDisplacement));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "PFPTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PFPValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexListName", m_nuVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MichelPDG", m_michelPDG));

    return BaseValidationTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
