/**
 *  @file   larpandoracontent/LArHelpers/LArHierarchyHelper.cc
 *
 *  @brief  Implementation of the lar hierarchy helper class.
 *
 *  $Log: $
 */

#include "Pandora/PdgTable.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

#include <numeric>
#include <unordered_set>

namespace lar_content
{

using namespace pandora;

LArHierarchyHelper::FoldingParameters::FoldingParameters() :
    m_foldToLeadingShowers{false},
    m_foldToTier{false},
    m_foldDynamic{false},
    m_cosAngleTolerance{0.9962f},
    m_tier{1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::FoldingParameters::FoldingParameters(const bool foldDynamic, const float cosAngleTolerance) :
    m_foldToLeadingShowers{false},
    m_foldToTier{false},
    m_foldDynamic{foldDynamic},
    m_cosAngleTolerance{cosAngleTolerance},
    m_tier{1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::FoldingParameters::FoldingParameters(const int foldingTier) :
    m_foldToLeadingShowers{false},
    m_foldToTier{true},
    m_foldDynamic{false},
    m_cosAngleTolerance{0.9962f},
    m_tier{foldingTier}
{
    if (m_tier < 1)
    {
        std::cout << "LArHierarchyHelper: Error - attempting to fold to non-positive tier" << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::QualityCuts::QualityCuts() :
    m_minPurity{0.8f},
    m_minCompleteness{0.65f},
    m_selectRecoHits{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::QualityCuts::QualityCuts(const float minPurity, const float minCompleteness, const bool selectRecoHits) :
    m_minPurity{minPurity},
    m_minCompleteness{minCompleteness},
    m_selectRecoHits{selectRecoHits}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::MCHierarchy() :
    m_nextNodeId{1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::MCHierarchy(const ReconstructabilityCriteria &recoCriteria) :
    m_recoCriteria(recoCriteria),
    m_nextNodeId{1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::~MCHierarchy()
{
    for (const auto &[pRoot, nodeVector] : m_interactions)
    {
        (void)pRoot;
        for (const Node *pNode : nodeVector)
            delete pNode;
    }
    m_interactions.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::FillHierarchy(const MCParticleList &mcParticleList, const CaloHitList &caloHitList, const FoldingParameters &foldParameters)
{
    const auto predicate = [](const MCParticle *pMCParticle) { return std::abs(pMCParticle->GetParticleId()) == NEUTRON; };
    m_mcToHitsMap.clear();
    for (const CaloHit *pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            m_mcToHitsMap[pMCParticle].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }

    MCParticleList rootNodes;
    for (const MCParticle *pMCParticle : mcParticleList)
    {
        const MCParticleList &parentList{pMCParticle->GetParentList()};
        if (parentList.empty())
        {
            rootNodes.emplace_back(pMCParticle);
        }
    }

    for (const MCParticle *pRoot : rootNodes)
    {
        MCParticleSet primarySet;
        LArHierarchyHelper::GetMCPrimaries(pRoot, primarySet);
        MCParticleList primaries(primarySet.begin(), primarySet.end());
        primaries.sort(LArMCParticleHelper::SortByMomentum);
        if (m_recoCriteria.m_removeNeutrons)
            primaries.erase(std::remove_if(primaries.begin(), primaries.end(), predicate), primaries.end());
        if (foldParameters.m_foldToTier && foldParameters.m_tier == 1)
        {
            for (const MCParticle *pPrimary : primaries)
            {
                MCParticleList allParticles{pPrimary};
                if (!m_recoCriteria.m_removeNeutrons)
                {
                    LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles);
                }
                else
                {
                    // Collect track-like and shower-like particles together, but throw out neutrons and descendents
                    MCParticleList dummy;
                    LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles, allParticles, dummy);
                }
                CaloHitList allHits;
                for (const MCParticle *pMCParticle : allParticles)
                {
                    // Not all MC particles will have hits
                    if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                    {
                        const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                        allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                    }
                }
                m_interactions[pRoot].emplace_back(new Node(*this, allParticles, allHits));
            }
        }
        else if (foldParameters.m_foldToLeadingShowers)
        {
            for (const MCParticle *pPrimary : primaries)
            {
                MCParticleList allParticles{pPrimary};
                int pdg{std::abs(pPrimary->GetParticleId())};
                const bool isShower{pdg == E_MINUS || pdg == PHOTON};
                const bool isNeutron{pdg == NEUTRON};
                if (isShower || (isNeutron && !m_recoCriteria.m_removeNeutrons))
                    LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles);
                CaloHitList allHits;
                for (const MCParticle *pMCParticle : allParticles)
                {
                    // ATTN - Not all MC particles will have hits
                    if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                    {
                        const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                        allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                    }
                }
                Node *pNode{new Node(*this, allParticles, allHits)};
                m_interactions[pRoot].emplace_back(pNode);
                if (!(isShower || isNeutron))
                {
                    // Find the children of this particle and recursively add them to the hierarchy
                    const MCParticleList &children{pPrimary->GetDaughterList()};
                    for (const MCParticle *pChild : children)
                        pNode->FillHierarchy(pChild, foldParameters);
                }
            }
        }
        else if (foldParameters.m_foldDynamic)
        {
            for (const MCParticle *pPrimary : primaries)
            {
                MCParticleList leadingParticles, childParticles;
                this->InterpretHierarchy(pPrimary, leadingParticles, childParticles, foldParameters.m_cosAngleTolerance);
                CaloHitList allHits;
                for (const MCParticle *pMCParticle : leadingParticles)
                {
                    // ATTN - Not all MC particles will have hits
                    if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                    {
                        const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                        allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                    }
                }

                Node *pNode{new Node(*this, leadingParticles, allHits)};
                m_interactions[pRoot].emplace_back(pNode);
                for (const MCParticle *pChild : childParticles)
                    pNode->FillHierarchy(pChild, foldParameters);
            }
        }
        else
        {
            // Unfolded and folded to tier > 1 have the same behaviour for primaries
            for (const MCParticle *pPrimary : primaries)
            {
                MCParticleList allParticles{pPrimary};
                CaloHitList allHits;
                for (const MCParticle *pMCParticle : allParticles)
                {
                    // ATTN - Not all MC particles will have hits
                    if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                    {
                        const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                        allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                    }
                }
                Node *pNode{new Node(*this, allParticles, allHits)};
                m_interactions[pRoot].emplace_back(pNode);
                // Find the children of this particle and recursively add them to the hierarchy
                const MCParticleList &children{pPrimary->GetDaughterList()};
                for (const MCParticle *pChild : children)
                    pNode->FillHierarchy(pChild, foldParameters);
            }
        }

        Node *pLeadingLepton{nullptr};
        float leadingLeptonEnergy{-std::numeric_limits<float>::max()};
        for (const Node *pNode : m_interactions[pRoot])
        {
            const MCParticle *pMC{pNode->GetLeadingMCParticle()};
            if (pMC)
            {
                const int pdg{std::abs(pMC->GetParticleId())};
                if ((pdg == MU_MINUS || pdg == E_MINUS || pdg == TAU_MINUS) && pMC->GetEnergy() > leadingLeptonEnergy)
                {
                    pLeadingLepton = const_cast<Node *>(pNode);
                    leadingLeptonEnergy = pMC->GetEnergy();
                }
            }
        }
        if (pLeadingLepton)
            pLeadingLepton->SetLeadingLepton();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArHierarchyHelper::MCHierarchy::NodeVector &LArHierarchyHelper::MCHierarchy::GetInteractions(const pandora::MCParticle *pRoot) const
{
    if (m_interactions.find(pRoot) == m_interactions.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return m_interactions.at(pRoot);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::GetRootMCParticles(MCParticleList &rootMCParticles) const
{
    for (auto iter = m_interactions.begin(); iter != m_interactions.end(); ++iter)
        rootMCParticles.emplace_back(iter->first);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::InterpretHierarchy(
    const MCParticle *const pRoot, MCParticleList &leadingParticles, MCParticleList &childParticles, const float cosAngleTolerance) const
{
    leadingParticles.emplace_back(pRoot);
    MCParticleList foldCandidates, childCandidates;
    const MCParticleList &children{pRoot->GetDaughterList()};
    for (const MCParticle *pMCParticle : children)
    {
        const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
        if (!pLArMCParticle)
            continue;
        if (LArMCParticleHelper::IsInelasticScatter(pMCParticle) || LArMCParticleHelper::IsElasticScatter(pMCParticle))
        {
            // Elastic and inelastic scattering can either lead to folding, distinct nodes or disposable hits, all other processes
            // are either distinct nodes, or disposable
            if (pMCParticle->GetParticleId() == pRoot->GetParticleId())
                foldCandidates.emplace_back(pMCParticle);
            else
                childCandidates.emplace_back(pMCParticle);
        }
        else if (!m_recoCriteria.m_removeNeutrons || (m_recoCriteria.m_removeNeutrons && pLArMCParticle->GetProcess() != MC_PROC_N_CAPTURE))
        {
            // Non-scattering process particles become leading candidates unless it's neutron capture and we're removing neutrons
            childCandidates.emplace_back(pMCParticle);
        }
    }
    const MCParticle *pBestFoldCandidate{nullptr};
    float bestDp{std::numeric_limits<float>::max()};
    for (const MCParticle *pMCParticle : foldCandidates)
    {
        if (foldCandidates.size() == 1)
        {
            // No alternative options, so this is either the best folding option by default, or a sufficiently large scatter to
            // treat as a new particle for reconstruction purposes
            if (LArMCParticleHelper::AreTopologicallyContinuous(pRoot, pMCParticle, cosAngleTolerance))
                pBestFoldCandidate = pMCParticle;
            else
                childCandidates.emplace_back(pMCParticle);
        }
        else
        {
            // Assess which, if any, of the children might be a continuation of the trajectory, otherwise move to child candidates
            if (LArMCParticleHelper::AreTopologicallyContinuous(pRoot, pMCParticle, cosAngleTolerance))
            {
                const float dp{pRoot->GetMomentum().GetMagnitude() - pMCParticle->GetMomentum().GetMagnitude()};
                if (dp < bestDp)
                {
                    pBestFoldCandidate = pMCParticle;
                    bestDp = dp;
                }
            }
            else
            {
                childCandidates.emplace_back(pMCParticle);
            }
        }
    }
    if (pBestFoldCandidate)
    {
        leadingParticles.emplace_back(pBestFoldCandidate);
        // Having found a particle to fold back at this level, continue to explore its downstream hierarchy for further folding
        // opportunities and make their respective children leading particles for the folded node we are creating
        this->CollectContinuations(pBestFoldCandidate, leadingParticles, childCandidates, cosAngleTolerance);
    }
    for (const MCParticle *pMCParticle : childCandidates)
    {
        // Consider if the child particle will produce enough downstream hits to warrant inclusion
        if (this->IsReconstructable(pMCParticle))
            childParticles.emplace_back(pMCParticle);
        else
        {
            MCParticleList localHierarchy{pMCParticle};
            CaloHitList localHits;
            LArMCParticleHelper::GetAllDescendentMCParticles(pMCParticle, localHierarchy);
            for (const MCParticle *pLocalMCParticle : localHierarchy)
            {
                if (m_mcToHitsMap.find(pLocalMCParticle) != m_mcToHitsMap.end())
                {
                    const CaloHitList &caloHits(m_mcToHitsMap.at(pLocalMCParticle));
                    localHits.insert(localHits.begin(), caloHits.begin(), caloHits.end());
                }
            }
            if (this->IsReconstructable(localHits))
                childParticles.emplace_back(pMCParticle);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::CollectContinuations(
    const MCParticle *pRoot, MCParticleList &continuingParticles, MCParticleList &childParticles, const float cosAngleTolerance) const
{
    const MCParticleList &children{pRoot->GetDaughterList()};
    MCParticleList foldCandidates;
    for (const MCParticle *pMCParticle : children)
    {
        const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
        if (!pLArMCParticle)
            continue;
        // Only elastic and inelastic scattering can lead to folding
        if (LArMCParticleHelper::IsInelasticScatter(pMCParticle) || LArMCParticleHelper::IsElasticScatter(pMCParticle))
        {
            if (pMCParticle->GetParticleId() == pRoot->GetParticleId())
                foldCandidates.emplace_back(pMCParticle);
        }
        else if (!m_recoCriteria.m_removeNeutrons || (m_recoCriteria.m_removeNeutrons && pLArMCParticle->GetProcess() != MC_PROC_N_CAPTURE))
        {
            // Non-scattering process particles become leading candidates unless it's neutron capture and we're removing neutrons
            childParticles.emplace_back(pMCParticle);
        }
    }
    const MCParticle *pBestFoldCandidate{nullptr};
    float bestDp{std::numeric_limits<float>::max()};
    for (const MCParticle *pMCParticle : foldCandidates)
    {
        if (foldCandidates.size() == 1)
        {
            // No alternative options, so this is either the best folding option by default, or a sufficiently large scatter to
            // treat as a new particle for reconstruction purposes
            if (LArMCParticleHelper::AreTopologicallyContinuous(pRoot, pMCParticle, cosAngleTolerance))
                pBestFoldCandidate = pMCParticle;
        }
        else
        {
            // Assess which, if any, of the children might be a continuation of the trajectory, otherwise move to child candidates
            if (LArMCParticleHelper::AreTopologicallyContinuous(pRoot, pMCParticle, cosAngleTolerance))
            {
                const float dp{pRoot->GetMomentum().GetMagnitude() - pMCParticle->GetMomentum().GetMagnitude()};
                if (dp < bestDp)
                {
                    pBestFoldCandidate = pMCParticle;
                    bestDp = dp;
                }
            }
        }
    }
    if (pBestFoldCandidate)
    {
        continuingParticles.emplace_back(pBestFoldCandidate);
        const MCParticleList &newLeadingParticles{pBestFoldCandidate->GetDaughterList()};
        // We need to add the children as child particles to ensure these sub-hierarchies are explored...
        childParticles.insert(childParticles.begin(), newLeadingParticles.begin(), newLeadingParticles.end());
        // but this current best fold candidate may have been added to the child particles by previously, so remove it
        const auto iter{std::find(childParticles.begin(), childParticles.end(), pBestFoldCandidate)};
        if (iter != childParticles.end())
            childParticles.erase(iter);
        // Having found a particle to fold back at this level, continue to explore its downstream hierarchy for further folding
        // opportunities and make their respective children child particles for the folded node we are creating
        LArHierarchyHelper::MCHierarchy::CollectContinuations(pBestFoldCandidate, continuingParticles, childParticles, cosAngleTolerance);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::GetFlattenedNodes(const MCParticle *const pRoot, NodeVector &nodeVector) const
{
    NodeList queue;
    for (const Node *pNode : m_interactions.at(pRoot))
    {
        nodeVector.emplace_back(pNode);
        queue.emplace_back(pNode);
    }
    while (!queue.empty())
    {
        const NodeVector &children{queue.front()->GetChildren()};
        queue.pop_front();
        for (const Node *pChild : children)
        {
            nodeVector.emplace_back(pChild);
            queue.emplace_back(pChild);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::RegisterNode(const Node *pNode)
{
    m_nodeToIdMap.insert(std::make_pair(pNode, m_nextNodeId));
    ++m_nextNodeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::MCHierarchy::ToString() const
{
    std::string str;
    for (const auto &[pRoot, nodeVector] : m_interactions)
    {
        const LArMCParticle *const pLArRoot{dynamic_cast<const LArMCParticle *const>(pRoot)};
        if (pLArRoot)
            str += "=== MC Interaction : PDG " + std::to_string(pLArRoot->GetParticleId()) +
                " Energy: " + std::to_string(pLArRoot->GetEnergy()) + " Nuance: " + std::to_string(pLArRoot->GetNuanceCode()) + "\n";
        else
            str += "=== MC Interaction : PDG " + std::to_string(pRoot->GetParticleId()) + " Energy: " + std::to_string(pRoot->GetEnergy()) + "\n";
        for (const Node *pNode : nodeVector)
            str += "   " + pNode->ToString("") + "\n";
    }

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::IsReconstructable(const pandora::MCParticle *pMCParticle) const
{
    if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
    {
        unsigned int nHitsU{0}, nHitsV{0}, nHitsW{0};
        for (const CaloHit *pCaloHit : m_mcToHitsMap.at(pMCParticle))
        {
            const HitType view{pCaloHit->GetHitType()};
            if (view == TPC_VIEW_U)
                ++nHitsU;
            else if (view == TPC_VIEW_V)
                ++nHitsV;
            else if (view == TPC_VIEW_W)
                ++nHitsW;
        }
        const unsigned int nHits{nHitsU + nHitsV + nHitsW};
        unsigned int nGoodViews{0};
        nGoodViews += nHitsU >= m_recoCriteria.m_minHitsForGoodView ? 1 : 0;
        nGoodViews += nHitsV >= m_recoCriteria.m_minHitsForGoodView ? 1 : 0;
        nGoodViews += nHitsW >= m_recoCriteria.m_minHitsForGoodView ? 1 : 0;

        return nHits >= m_recoCriteria.m_minHits && nGoodViews >= m_recoCriteria.m_minGoodViews;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::IsReconstructable(const CaloHitList &caloHits) const
{
    unsigned int nHitsU{0}, nHitsV{0}, nHitsW{0};
    for (const CaloHit *pCaloHit : caloHits)
    {
        const HitType view{pCaloHit->GetHitType()};
        if (view == TPC_VIEW_U)
            ++nHitsU;
        else if (view == TPC_VIEW_V)
            ++nHitsV;
        else if (view == TPC_VIEW_W)
            ++nHitsW;
    }
    const unsigned int nHits{nHitsU + nHitsV + nHitsW};
    unsigned int nGoodViews{0};
    nGoodViews += nHitsU >= m_recoCriteria.m_minHitsForGoodView ? 1 : 0;
    nGoodViews += nHitsV >= m_recoCriteria.m_minHitsForGoodView ? 1 : 0;
    nGoodViews += nHitsW >= m_recoCriteria.m_minHitsForGoodView ? 1 : 0;

    return nHits >= m_recoCriteria.m_minHits && nGoodViews >= m_recoCriteria.m_minGoodViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::Node::Node(MCHierarchy &hierarchy, const MCParticle *pMCParticle, const int tier) :
    m_hierarchy(hierarchy),
    m_mainParticle(pMCParticle),
    m_tier{tier},
    m_pdg{0},
    m_isLeadingLepton{false}
{
    if (pMCParticle)
    {
        m_pdg = pMCParticle->GetParticleId();
        m_mcParticles.emplace_back(pMCParticle);
    }
    m_hierarchy.RegisterNode(this);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::Node::Node(MCHierarchy &hierarchy, const MCParticleList &mcParticleList, const CaloHitList &caloHitList, const int tier) :
    m_hierarchy(hierarchy),
    m_mcParticles(mcParticleList),
    m_caloHits(caloHitList),
    m_mainParticle(nullptr),
    m_tier{tier},
    m_pdg{0},
    m_isLeadingLepton{false}
{
    if (!mcParticleList.empty())
    {
        m_mainParticle = mcParticleList.front();
        m_pdg = m_mainParticle->GetParticleId();
    }
    m_mcParticles.sort(LArMCParticleHelper::SortByMomentum);
    m_caloHits.sort();
    m_hierarchy.RegisterNode(this);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::Node::~Node()
{
    m_mcParticles.clear();
    m_caloHits.clear();
    for (const Node *node : m_children)
        delete node;
    m_children.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::Node::FillHierarchy(const MCParticle *pRoot, const FoldingParameters &foldParameters)
{
    if (foldParameters.m_foldDynamic)
    {
        MCParticleList leadingParticles, childParticles;
        m_hierarchy.InterpretHierarchy(pRoot, leadingParticles, childParticles, foldParameters.m_cosAngleTolerance);
        CaloHitList allHits;
        for (const MCParticle *pMCParticle : leadingParticles)
        {
            // ATTN - Not all MC particles will have hits
            if (m_hierarchy.m_mcToHitsMap.find(pMCParticle) != m_hierarchy.m_mcToHitsMap.end())
            {
                const CaloHitList &caloHits(m_hierarchy.m_mcToHitsMap.at(pMCParticle));
                allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
            }
        }

        Node *pNode{new Node(m_hierarchy, leadingParticles, allHits, this->m_tier + 1)};
        m_children.emplace_back(pNode);
        for (const MCParticle *pChild : childParticles)
            pNode->FillHierarchy(pChild, foldParameters);
    }
    else
    {
        MCParticleList allParticles{pRoot};
        const int pdg{std::abs(pRoot->GetParticleId())};
        const bool isShower{pdg == E_MINUS || pdg == PHOTON};
        const bool isNeutron{pdg == NEUTRON};

        if (foldParameters.m_foldToTier && LArMCParticleHelper::GetHierarchyTier(pRoot) >= foldParameters.m_tier)
            LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles);
        else if (foldParameters.m_foldToLeadingShowers && (isShower || (isNeutron && !m_hierarchy.m_recoCriteria.m_removeNeutrons)))
            LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles);
        else if (m_hierarchy.m_recoCriteria.m_removeNeutrons && isNeutron)
            return;

        CaloHitList allHits;
        for (const MCParticle *pMCParticle : allParticles)
        {
            // ATTN - Not all MC particles will have hits
            if (m_hierarchy.m_mcToHitsMap.find(pMCParticle) != m_hierarchy.m_mcToHitsMap.end())
            {
                const CaloHitList &caloHits(m_hierarchy.m_mcToHitsMap.at(pMCParticle));
                allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
            }
        }

        if (!allParticles.empty())
        {
            const bool hasChildren{(foldParameters.m_foldToTier && LArMCParticleHelper::GetHierarchyTier(pRoot) < foldParameters.m_tier) ||
                (!foldParameters.m_foldToTier && !foldParameters.m_foldToLeadingShowers) ||
                (foldParameters.m_foldToLeadingShowers && !(isShower || isNeutron))};
            // Only add the node if it either has children, or is a leaf node with hits
            if (hasChildren || (!hasChildren && !allHits.empty()))
            {
                Node *pNode{new Node(m_hierarchy, allParticles, allHits, this->m_tier + 1)};
                m_children.emplace_back(pNode);
                if (hasChildren)
                {
                    // Find the children of this particle and recursively add them to the hierarchy
                    const MCParticleList &children{pRoot->GetDaughterList()};
                    for (const MCParticle *pChild : children)
                        pNode->FillHierarchy(pChild, foldParameters);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::Node::FillFlat(const MCParticle *pRoot)
{
    MCParticleList allParticles{pRoot};
    if (!m_hierarchy.m_recoCriteria.m_removeNeutrons)
    {
        LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles);
    }
    else
    {
        MCParticleList neutrons;
        LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles, allParticles, neutrons);
    }
    CaloHitList allHits;
    for (const MCParticle *pMCParticle : allParticles)
    {
        // ATTN - Not all MC particles will have hits
        if (m_hierarchy.m_mcToHitsMap.find(pMCParticle) != m_hierarchy.m_mcToHitsMap.end())
        {
            const CaloHitList &caloHits(m_hierarchy.m_mcToHitsMap.at(pMCParticle));
            allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
        }
    }
    if (!allParticles.empty())
    {
        Node *pNode{new Node(m_hierarchy, allParticles, allHits, this->m_tier + 1)};
        m_children.emplace_back(pNode);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArHierarchyHelper::MCHierarchy::Node::GetId() const
{
    return m_hierarchy.m_nodeToIdMap.at(this);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::Node::IsReconstructable() const
{
    const bool enoughHits{m_caloHits.size() >= m_hierarchy.m_recoCriteria.m_minHits};
    if (!enoughHits)
        return false;
    bool enoughGoodViews{false};
    unsigned int nHitsU{0}, nHitsV{0}, nHitsW{0};
    for (const CaloHit *pCaloHit : m_caloHits)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                ++nHitsU;
                break;
            case TPC_VIEW_V:
                ++nHitsV;
                break;
            case TPC_VIEW_W:
                ++nHitsW;
                break;
            default:
                break;
        }
        unsigned int nGoodViews{0};
        if (nHitsU >= m_hierarchy.m_recoCriteria.m_minHitsForGoodView)
            ++nGoodViews;
        if (nHitsV >= m_hierarchy.m_recoCriteria.m_minHitsForGoodView)
            ++nGoodViews;
        if (nHitsW >= m_hierarchy.m_recoCriteria.m_minHitsForGoodView)
            ++nGoodViews;
        if (nGoodViews >= m_hierarchy.m_recoCriteria.m_minGoodViews)
        {
            enoughGoodViews = true;
            break;
        }
    }

    return enoughGoodViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::Node::IsTestBeamParticle() const
{
    if (m_mainParticle)
        return LArMCParticleHelper::IsBeamParticle(m_mainParticle);
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::Node::IsCosmicRay() const
{
    if (m_mainParticle)
        return LArMCParticleHelper::IsCosmicRay(m_mainParticle);
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::MCHierarchy::Node::ToString(const std::string &prefix) const
{
    std::string str(prefix + "PDG: " + std::to_string(m_pdg) + " Energy: " + std::to_string(m_mainParticle ? m_mainParticle->GetEnergy() : 0) +
        " Hits: " + std::to_string(m_caloHits.size()) + "\n");
    for (const Node *pChild : m_children)
        str += pChild->ToString(prefix + "   ");

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria::ReconstructabilityCriteria() :
    m_minHits{30},
    m_minHitsForGoodView{10},
    m_minGoodViews{2},
    m_removeNeutrons{true}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria::ReconstructabilityCriteria(const ReconstructabilityCriteria &obj) :
    m_minHits{obj.m_minHits},
    m_minHitsForGoodView{obj.m_minHitsForGoodView},
    m_minGoodViews{obj.m_minGoodViews},
    m_removeNeutrons{obj.m_removeNeutrons}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria::ReconstructabilityCriteria(
    const unsigned int minHits, const unsigned int minHitsForGoodView, const unsigned int minGoodViews, const bool removeNeutrons) :
    m_minHits{minHits},
    m_minHitsForGoodView{minHitsForGoodView},
    m_minGoodViews{minGoodViews},
    m_removeNeutrons{removeNeutrons}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::RecoHierarchy()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::~RecoHierarchy()
{
    for (const auto &[pRoot, nodeVector] : m_interactions)
    {
        for (const Node *pNode : nodeVector)
            delete pNode;
    }
    m_interactions.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::FillHierarchy(const PfoList &pfoList, const FoldingParameters &foldParameters)
{
    PfoList rootNodes;
    for (const ParticleFlowObject *pPfo : pfoList)
    {
        const PfoList &parentList{pPfo->GetParentPfoList()};
        if (parentList.empty())
        {
            rootNodes.emplace_back(pPfo);
        }
    }

    for (const ParticleFlowObject *const pRoot : rootNodes)
    {
        PfoSet primarySet;
        LArHierarchyHelper::GetRecoPrimaries(pRoot, primarySet);
        PfoList primaries(primarySet.begin(), primarySet.end());
        primaries.sort(LArPfoHelper::SortByNHits);
        if (foldParameters.m_foldToTier && foldParameters.m_tier == 1)
        {
            for (const ParticleFlowObject *pPrimary : primaries)
            {
                PfoList allParticles;
                // ATTN - pPrimary gets added to the list of downstream PFOs, not just the child PFOs
                LArPfoHelper::GetAllDownstreamPfos(pPrimary, allParticles);
                CaloHitList allHits;
                for (const ParticleFlowObject *pPfo : allParticles)
                    LArPfoHelper::GetAllCaloHits2D(pPfo, allHits);
                m_interactions[pRoot].emplace_back(new Node(*this, allParticles, allHits));
            }
        }
        else if (foldParameters.m_foldToLeadingShowers)
        {
            for (const ParticleFlowObject *pPrimary : primaries)
            {
                PfoList allParticles;
                int pdg{std::abs(pPrimary->GetParticleId())};
                const bool isShower{pdg == E_MINUS};
                // ATTN - pPrimary gets added to the list of downstream PFOs, not just the child PFOs
                if (isShower)
                    LArPfoHelper::GetAllDownstreamPfos(pPrimary, allParticles);
                else
                    allParticles.emplace_back(pPrimary);

                CaloHitList allHits;
                for (const ParticleFlowObject *pPfo : allParticles)
                    LArPfoHelper::GetAllCaloHits2D(pPfo, allHits);
                Node *pNode{new Node(*this, allParticles, allHits)};
                m_interactions[pRoot].emplace_back(pNode);
                if (!isShower)
                {
                    // Find the children of this particle and recursively add them to the hierarchy
                    const PfoList &children{pPrimary->GetDaughterPfoList()};
                    for (const ParticleFlowObject *pChild : children)
                        pNode->FillHierarchy(pChild, foldParameters);
                }
            }
        }
        else
        {
            // Dynamic fold, Unfolded and fold to tier > 1 have the same behaviour for primaries
            for (const ParticleFlowObject *pPrimary : primaries)
            {
                PfoList allParticles{pPrimary};
                CaloHitList allHits;
                for (const ParticleFlowObject *pPfo : allParticles)
                    LArPfoHelper::GetAllCaloHits2D(pPfo, allHits);
                Node *pNode{new Node(*this, allParticles, allHits)};
                m_interactions[pRoot].emplace_back(pNode);
                // Find the children of this particle and recursively add them to the hierarchy
                const PfoList &children{pPrimary->GetDaughterPfoList()};
                for (const ParticleFlowObject *pChild : children)
                    pNode->FillHierarchy(pChild, foldParameters.m_foldToLeadingShowers);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::RecoHierarchy::GetInteractions(const pandora::ParticleFlowObject *pRoot) const
{
    if (m_interactions.find(pRoot) == m_interactions.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return m_interactions.at(pRoot);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::GetRootPfos(PfoList &rootPfos) const
{
    for (auto iter = m_interactions.begin(); iter != m_interactions.end(); ++iter)
        rootPfos.emplace_back(iter->first);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::GetFlattenedNodes(const ParticleFlowObject *pRoot, NodeVector &nodeVector) const
{
    NodeList queue;
    for (const Node *pNode : m_interactions.at(pRoot))
    {
        nodeVector.emplace_back(pNode);
        queue.emplace_back(pNode);
    }
    while (!queue.empty())
    {
        const NodeVector &children{queue.front()->GetChildren()};
        queue.pop_front();
        for (const Node *pChild : children)
        {
            nodeVector.emplace_back(pChild);
            queue.emplace_back(pChild);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::RecoHierarchy::ToString() const
{
    std::string str;
    for (const auto &[pRoot, nodeVector] : m_interactions)
    {
        str += "=== Reco Interaction : PDG " + std::to_string(pRoot->GetParticleId()) + "\n";
        for (const Node *pNode : nodeVector)
            str += "   " + pNode->ToString("") + "\n";
    }

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::Node::Node(const RecoHierarchy &hierarchy, const ParticleFlowObject *pPfo, const int tier) :
    m_hierarchy{hierarchy},
    m_mainPfo{pPfo},
    m_tier{tier},
    m_pdg{0}
{
    if (pPfo)
    {
        m_pdg = pPfo->GetParticleId();
        m_pfos.emplace_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::Node::Node(const RecoHierarchy &hierarchy, const PfoList &pfoList, const CaloHitList &caloHitList, const int tier) :
    m_hierarchy(hierarchy),
    m_mainPfo{nullptr},
    m_tier{tier},
    m_pdg{0}
{
    if (!pfoList.empty())
    {
        m_mainPfo = pfoList.front();
        m_pdg = pfoList.front()->GetParticleId();
    }
    m_pfos = pfoList;
    m_pfos.sort(LArPfoHelper::SortByNHits);
    m_caloHits = caloHitList;
    m_caloHits.sort();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::Node::~Node()
{
    m_pfos.clear();
    m_caloHits.clear();
    for (const Node *node : m_children)
        delete node;
    m_children.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::Node::FillHierarchy(const ParticleFlowObject *pRoot, const FoldingParameters &foldParameters)
{
    PfoList allParticles;
    int pdg{std::abs(pRoot->GetParticleId())};
    const bool isShower{pdg == E_MINUS};
    if (foldParameters.m_foldToTier && LArPfoHelper::GetHierarchyTier(pRoot) >= foldParameters.m_tier)
        LArPfoHelper::GetAllDownstreamPfos(pRoot, allParticles);
    else if (foldParameters.m_foldToLeadingShowers && isShower)
        LArPfoHelper::GetAllDownstreamPfos(pRoot, allParticles);
    else
        allParticles.emplace_back(pRoot);

    CaloHitList allHits;
    for (const ParticleFlowObject *pPfo : allParticles)
        LArPfoHelper::GetAllCaloHits2D(pPfo, allHits);
    const bool hasChildren{(foldParameters.m_foldToTier && LArPfoHelper::GetHierarchyTier(pRoot) < foldParameters.m_tier) ||
        (!foldParameters.m_foldToTier && !foldParameters.m_foldToLeadingShowers) || (foldParameters.m_foldToLeadingShowers && !isShower)};

    if (hasChildren || (!hasChildren && !allHits.empty()))
    {
        Node *pNode{new Node(m_hierarchy, allParticles, allHits, m_tier + 1)};
        m_children.emplace_back(pNode);

        if (hasChildren)
        {
            const PfoList &children{pRoot->GetDaughterPfoList()};
            for (const ParticleFlowObject *pChild : children)
                pNode->FillHierarchy(pChild, foldParameters);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::Node::FillFlat(const ParticleFlowObject *pRoot)
{
    PfoList allParticles;
    LArPfoHelper::GetAllDownstreamPfos(pRoot, allParticles);
    CaloHitList allHits;
    for (const ParticleFlowObject *pPfo : allParticles)
        LArPfoHelper::GetAllCaloHits2D(pPfo, allHits);
    Node *pNode{new Node(m_hierarchy, allParticles, allHits, m_tier + 1)};
    m_children.emplace_back(pNode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList &LArHierarchyHelper::RecoHierarchy::Node::GetRecoParticles() const
{
    return m_pfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &LArHierarchyHelper::RecoHierarchy::Node::GetCaloHits() const
{
    return m_caloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArHierarchyHelper::RecoHierarchy::Node::GetParticleId() const
{
    return m_pdg;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::RecoHierarchy::Node::ToString(const std::string &prefix) const
{
    std::string str(
        prefix + "PDG: " + std::to_string(m_pdg) + " Tier: " + std::to_string(m_tier) + " Hits: " + std::to_string(m_caloHits.size()) + "\n");
    for (const Node *pChild : m_children)
        str += pChild->ToString(prefix + "   ");

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCMatches::MCMatches(const MCHierarchy::Node *pMCParticle) :
    m_pMCParticle{pMCParticle}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCMatches::AddRecoMatch(const RecoHierarchy::Node *pReco, const int nSharedHits, const CaloHitList &selectedRecoHits)
{
    m_recoNodes.emplace_back(pReco);
    m_sharedHits.emplace_back(nSharedHits);
    m_selectedRecoHitsMap.insert(std::make_pair(pReco, selectedRecoHits));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList LArHierarchyHelper::MCMatches::GetSelectedRecoHits(const RecoHierarchy::Node *pReco) const
{
    auto iter{m_selectedRecoHitsMap.find(pReco)};
    if (iter != m_selectedRecoHitsMap.end())
        return iter->second;
    else
        return pReco->GetCaloHits();
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MCMatches::GetSharedHits(const RecoHierarchy::Node *pReco) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    int index = iter - m_recoNodes.begin();

    return static_cast<int>(m_sharedHits[index]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetPurity(const RecoHierarchy::Node *pReco, const bool adcWeighted) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CaloHitList recoHits = LArHierarchyHelper::MCMatches::GetSelectedRecoHits(pReco);
    const CaloHitList &mcHits{m_pMCParticle->GetCaloHits()};
    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetPurity(intersection, recoHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetPurity(const RecoHierarchy::Node *pReco, const HitType view, const bool adcWeighted) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CaloHitList recoHits;
    for (const CaloHit *pCaloHit : LArHierarchyHelper::MCMatches::GetSelectedRecoHits(pReco))
        if (pCaloHit->GetHitType() == view)
            recoHits.emplace_back(pCaloHit);
    CaloHitList mcHits;
    for (const CaloHit *pCaloHit : m_pMCParticle->GetCaloHits())
        if (pCaloHit->GetHitType() == view)
            mcHits.emplace_back(pCaloHit);

    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetPurity(intersection, recoHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetCompleteness(const RecoHierarchy::Node *pReco, const bool adcWeighted) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CaloHitList recoHits = LArHierarchyHelper::MCMatches::GetSelectedRecoHits(pReco);
    const CaloHitList &mcHits{m_pMCParticle->GetCaloHits()};
    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetCompleteness(intersection, mcHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetCompleteness(const RecoHierarchy::Node *pReco, const HitType view, const bool adcWeighted) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CaloHitList recoHits;
    for (const CaloHit *pCaloHit : LArHierarchyHelper::MCMatches::GetSelectedRecoHits(pReco))
        if (pCaloHit->GetHitType() == view)
            recoHits.emplace_back(pCaloHit);
    CaloHitList mcHits;
    for (const CaloHit *pCaloHit : m_pMCParticle->GetCaloHits())
        if (pCaloHit->GetHitType() == view)
            mcHits.emplace_back(pCaloHit);

    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetCompleteness(intersection, mcHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetPurity(const CaloHitVector &intersection, const CaloHitList &recoHits, const bool adcWeighted) const
{
    float purity{0.f};
    if (!intersection.empty())
    {
        if (adcWeighted)
        {
            float adcSum{0.f};
            for (const CaloHit *pCaloHit : recoHits)
                adcSum += pCaloHit->GetInputEnergy();
            if (adcSum > std::numeric_limits<float>::epsilon())
            {
                for (const CaloHit *pCaloHit : intersection)
                    purity += pCaloHit->GetInputEnergy();
                purity /= adcSum;
            }
        }
        else
        {
            purity = intersection.size() / static_cast<float>(recoHits.size());
        }
    }

    return purity;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetCompleteness(const CaloHitVector &intersection, const CaloHitList &mcHits, const bool adcWeighted) const
{
    float completeness{0.f};
    if (!intersection.empty())
    {
        if (adcWeighted)
        {
            float adcSum{0.f};
            for (const CaloHit *pCaloHit : mcHits)
                adcSum += pCaloHit->GetInputEnergy();
            if (adcSum > std::numeric_limits<float>::epsilon())
            {
                for (const CaloHit *pCaloHit : intersection)
                    completeness += pCaloHit->GetInputEnergy();
                completeness /= adcSum;
            }
        }
        else
        {
            completeness = intersection.size() / static_cast<float>(mcHits.size());
        }
    }

    return completeness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCMatches::IsQuality(const LArHierarchyHelper::QualityCuts &qualityCuts) const
{
    if (m_recoNodes.empty())
        return false;

    int nAboveThreshold{0};
    if (m_recoNodes.size() != 1)
    {
        for (const LArHierarchyHelper::RecoHierarchy::Node *const pNode : m_recoNodes)
            if (this->GetCompleteness(pNode) > 0.1f)
                ++nAboveThreshold;
        if (nAboveThreshold != 1)
            return false;
    }

    if (this->GetPurity(m_recoNodes.front()) < qualityCuts.m_minPurity)
        return false;

    if (this->GetCompleteness(m_recoNodes.front()) < qualityCuts.m_minCompleteness)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MatchInfo::MatchInfo(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy) :
    MatchInfo(mcHierarchy, recoHierarchy, QualityCuts())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MatchInfo::MatchInfo(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy, const QualityCuts &qualityCuts) :
    m_mcHierarchy{mcHierarchy},
    m_recoHierarchy{recoHierarchy},
    m_qualityCuts{qualityCuts}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchInfo::Match()
{
    MCParticleList rootMCParticles;
    m_mcHierarchy.GetRootMCParticles(rootMCParticles);
    PfoList rootPfos;
    m_recoHierarchy.GetRootPfos(rootPfos);
    std::map<const MCHierarchy::Node *, MCMatches> mcToMatchMap;

    for (const MCParticle *const pRootMC : rootMCParticles)
    {
        MCHierarchy::NodeVector mcNodes;
        m_mcHierarchy.GetFlattenedNodes(pRootMC, mcNodes);
        std::sort(mcNodes.begin(), mcNodes.end(), [](const MCHierarchy::Node *lhs, const MCHierarchy::Node *rhs)
            { return lhs->GetCaloHits().size() > rhs->GetCaloHits().size(); });

        // Get all of the hits from the MC nodes (for selecting reco hits)
        CaloHitList allMCHits;
        std::unordered_map<const CaloHit*, const MCHierarchy::Node*> hitToMCNodeMap;
        std::vector<bool> isReconstructable;
        isReconstructable.reserve(mcNodes.size());
        std::unordered_set<intptr_t> mcHitIds;

        for (const MCHierarchy::Node *pMCNode : mcNodes)
        {
            const CaloHitList &mcHits{pMCNode->GetCaloHits()};
            allMCHits.insert(allMCHits.begin(), mcHits.begin(), mcHits.end());

            const bool reconstructable{pMCNode->IsReconstructable()};
            isReconstructable.push_back(reconstructable);

            if (reconstructable)
                for (const CaloHit *pMCHit : mcHits)
                    hitToMCNodeMap[pMCHit] = pMCNode;
        }

        // Populate a map of MC hit IDs for fast lookup if we're applying the
        // reco hit selection quality cut.
        if (m_qualityCuts.m_selectRecoHits)
        {
            mcHitIds.reserve(allMCHits.size());
            for (const CaloHit *pMCHit : allMCHits)
                mcHitIds.insert(reinterpret_cast<intptr_t>(pMCHit->GetParentAddress()));
        }

        for (const ParticleFlowObject *const pRootPfo : rootPfos)
        {
            RecoHierarchy::NodeVector recoNodes;
            m_recoHierarchy.GetFlattenedNodes(pRootPfo, recoNodes);

            std::sort(recoNodes.begin(), recoNodes.end(), [](const RecoHierarchy::Node *lhs, const RecoHierarchy::Node *rhs)
                { return lhs->GetCaloHits().size() > rhs->GetCaloHits().size(); });

            for (const RecoHierarchy::Node *pRecoNode : recoNodes)
            {
                // Get the selected list of reco hits that overlap with all of the MC hits
                // or just use all of the hits in the reco node
                const CaloHitList selectedRecoHits = (m_qualityCuts.m_selectRecoHits == true)
                    ? LArHierarchyHelper::MatchInfo::GetSelectedRecoHits(pRecoNode, allMCHits, &mcHitIds)
                    : pRecoNode->GetCaloHits();

                // Count shared hits using the map instead of O(N^2) set_intersection
                std::unordered_map<const MCHierarchy::Node*, size_t> sharedCounts;
                for (const CaloHit *pRecoHit : selectedRecoHits)
                {
                    auto it = hitToMCNodeMap.find(pRecoHit);
                    if (it != hitToMCNodeMap.end())
                        sharedCounts[it->second]++;
                }

                const MCHierarchy::Node *pBestNode{nullptr};
                size_t bestSharedHits{0};

                for (auto const& [pMCNode, count] : sharedCounts)
                {
                    if (count > bestSharedHits)
                    {
                        bestSharedHits = count;
                        pBestNode = pMCNode;
                    }
                }

                if (pBestNode)
                {
                    auto iter{mcToMatchMap.find(pBestNode)};
                    if (iter != mcToMatchMap.end())
                    {
                        MCMatches &match(iter->second);
                        match.AddRecoMatch(pRecoNode, static_cast<int>(bestSharedHits), selectedRecoHits);
                    }
                    else
                    {
                        MCMatches match(pBestNode);
                        match.AddRecoMatch(pRecoNode, static_cast<int>(bestSharedHits), selectedRecoHits);
                        mcToMatchMap.insert(std::make_pair(pBestNode, match));
                    }
                }
                else
                {
                    m_unmatchedReco.emplace_back(pRecoNode);
                }
            }
        }
    }

    for (auto [pMCNode, matches] : mcToMatchMap)
    {
        // We need to figure out which MC interaction hierarchy the matches belongs to
        for (const MCParticle *const pRootMC : rootMCParticles)
        {
            MCHierarchy::NodeVector mcNodes;
            m_mcHierarchy.GetFlattenedNodes(pRootMC, mcNodes);
            if (std::find(mcNodes.begin(), mcNodes.end(), pMCNode) != mcNodes.end())
            {
                m_matches[pRootMC].emplace_back(matches);
                break;
            }
        }
    }

    const auto predicate = [](const MCMatches &lhs, const MCMatches &rhs)
    { return lhs.GetMC()->GetCaloHits().size() > rhs.GetMC()->GetCaloHits().size(); };

    for (const MCParticle *const pRootMC : rootMCParticles)
    {
        std::sort(m_matches[pRootMC].begin(), m_matches[pRootMC].end(), predicate);

        MCHierarchy::NodeVector mcNodes;
        m_mcHierarchy.GetFlattenedNodes(pRootMC, mcNodes);

        for (const MCHierarchy::Node *pMCNode : mcNodes)
        {
            if (pMCNode->IsReconstructable() && mcToMatchMap.find(pMCNode) == mcToMatchMap.end())
            {
                MCMatches match(pMCNode);
                m_matches[pRootMC].emplace_back(match);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNMCNodes(const MCParticle *const pRoot) const
{
    if (m_matches.find(pRoot) == m_matches.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return static_cast<unsigned int>(m_matches.at(pRoot).size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNNeutrinoMCNodes(const MCParticle *const pRoot) const
{
    if (m_matches.find(pRoot) == m_matches.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    unsigned int nNodes{0};
    for (const MCMatches &match : m_matches.at(pRoot))
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (!(pNode->IsCosmicRay() || pNode->IsTestBeamParticle()))
            ++nNodes;
    }

    return nNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNCosmicRayMCNodes(const MCParticle *const pRoot) const
{
    if (m_matches.find(pRoot) == m_matches.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    unsigned int nNodes{0};
    for (const MCMatches &match : m_matches.at(pRoot))
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsCosmicRay())
            ++nNodes;
    }

    return nNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNTestBeamMCNodes(const MCParticle *const pRoot) const
{
    if (m_matches.find(pRoot) == m_matches.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    unsigned int nNodes{0};
    for (const MCMatches &match : m_matches.at(pRoot))
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsTestBeamParticle())
            ++nNodes;
    }

    return nNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList LArHierarchyHelper::MatchInfo::GetSelectedRecoHits(
    const RecoHierarchy::Node *pRecoNode, const CaloHitList &allMCHits, const std::unordered_set<intptr_t> *pMCHitIds) const
{
    // Select all of the reco node hit Ids that overlap with the allMCHits Ids
    CaloHitList selectedHits;
    if (!pRecoNode)
        return selectedHits;

    const std::unordered_set<intptr_t> *pIdsToUse = pMCHitIds;
    std::unordered_set<intptr_t> mcHitIds;

    // Build a map of MC hit IDs, for fast lookup, if one wasn't provided.
    if (!pIdsToUse)
    {
        mcHitIds.reserve(allMCHits.size());
        for (const CaloHit *pMCHit : allMCHits)
            mcHitIds.insert(reinterpret_cast<intptr_t>(pMCHit->GetParentAddress()));
        pIdsToUse = &mcHitIds;
    }

    for (const CaloHit *pRecoHit : pRecoNode->GetCaloHits())
    {
        if (pIdsToUse->find(reinterpret_cast<intptr_t>(pRecoHit->GetParentAddress())) != pIdsToUse->end())
            selectedHits.emplace_back(pRecoHit);
    }
    return selectedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchInfo::Print(const MCHierarchy &mcHierarchy) const
{
    const std::streamsize originalPrecision{std::cout.precision()};
    const std::ios_base::fmtflags originalFormat{std::cout.flags()};

    MCParticleList rootMCParticles;
    mcHierarchy.GetRootMCParticles(rootMCParticles);

    for (const MCParticle *const pRootMC : rootMCParticles)
    {
        const LArHierarchyHelper::MCMatchesVector &matches{this->GetMatches(pRootMC)};

        MCParticleList primaries;
        for (const LArHierarchyHelper::MCMatches &match : matches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pMCNode{match.GetMC()};
            if (pMCNode->GetHierarchyTier() == 1)
            {
                const MCParticle *const pLeadingMC{pMCNode->GetLeadingMCParticle()};
                primaries.emplace_back(pLeadingMC);
            }
        }
        if (primaries.size() == 0)
            continue;
        primaries.sort(LArMCParticleHelper::SortByMomentum);
        const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primaries)};

        const LArMCParticle *const pLArRoot{dynamic_cast<const LArMCParticle *const>(pRootMC)};
        if (pLArRoot)
            std::cout << "=== MC Interaction : PDG " << std::to_string(pLArRoot->GetParticleId())
                      << " Energy: " << std::to_string(pLArRoot->GetEnergy()) << " Type: " << descriptor.ToString() << std::endl;
        else
            std::cout << "=== MC Interaction : PDG " << std::to_string(pRootMC->GetParticleId())
                      << " Energy: " << std::to_string(pRootMC->GetEnergy()) << " Type: " << descriptor.ToString() << std::endl;

        unsigned int nNeutrinoMCParticles{this->GetNNeutrinoMCNodes(pRootMC)}, nNeutrinoRecoParticles{0};
        unsigned int nCosmicMCParticles{this->GetNCosmicRayMCNodes(pRootMC)}, nCosmicRecoParticles{0};
        unsigned int nTestBeamMCParticles{this->GetNTestBeamMCNodes(pRootMC)}, nTestBeamRecoParticles{0};
        std::cout << "   === Matches ===" << std::endl;
        std::cout << std::fixed << std::setprecision(2);
        for (const MCMatches &match : m_matches.at(pRootMC))
        {
            const MCHierarchy::Node *pMCNode{match.GetMC()};
            const int pdg{pMCNode->GetParticleId()};
            const size_t mcHits{pMCNode->GetCaloHits().size()};
            const std::string tag{pMCNode->IsTestBeamParticle() ? "(Beam) " : pMCNode->IsCosmicRay() ? "(Cosmic) " : ""};
            std::cout << "   MC " << tag << pdg << " hits " << mcHits << std::endl;
            const RecoHierarchy::NodeVector &nodeVector{match.GetRecoMatches()};

            for (const RecoHierarchy::Node *pRecoNode : nodeVector)
            {
                const CaloHitList recoHits = match.GetSelectedRecoHits(pRecoNode);
                const unsigned int nRecoHits{static_cast<unsigned int>(recoHits.size())};
                const unsigned int sharedHits{match.GetSharedHits(pRecoNode)};
                const float purity{match.GetPurity(pRecoNode)};
                const float completeness{match.GetCompleteness(pRecoNode)};
                if (completeness > 0.1f)
                    std::cout << "   Matched " << sharedHits << " out of " << nRecoHits << " with purity " << purity << " and completeness "
                              << completeness << std::endl;
                else
                    std::cout << "   (Below threshold) " << sharedHits << " out of " << nRecoHits << " with purity " << purity
                              << " and completeness " << completeness << std::endl;
            }
            if (nodeVector.empty())
            {
                std::cout << "      Unmatched" << std::endl;
            }
            else if (match.IsQuality(this->GetQualityCuts()))
            {
                if (pMCNode->IsTestBeamParticle())
                    ++nTestBeamRecoParticles;
                else if (pMCNode->IsCosmicRay())
                    ++nCosmicRecoParticles;
                else
                    ++nNeutrinoRecoParticles;
            }
        }

        if (LArMCParticleHelper::IsNeutrino(pRootMC))
        {
            std::cout << "   Neutrino Interaction Summary:" << std::endl;
            if (nNeutrinoMCParticles)
            {
                std::cout << "   Good final state particles: " << nNeutrinoRecoParticles << " of " << nNeutrinoMCParticles << " : "
                          << (100 * nNeutrinoRecoParticles / static_cast<float>(nNeutrinoMCParticles)) << "%" << std::endl;
            }
        }
        else if (LArMCParticleHelper::IsCosmicRay(pRootMC))
        {
            std::cout << "   Cosmic Ray Interaction Summary:" << std::endl;
            std::cout << std::fixed << std::setprecision(1);
            if (nCosmicMCParticles)
            {
                std::cout << "   Good cosmics: " << nCosmicRecoParticles << " of " << nCosmicMCParticles << " : "
                          << (100 * nCosmicRecoParticles / static_cast<float>(nCosmicMCParticles)) << "%" << std::endl;
            }
        }
        else if (LArMCParticleHelper::IsBeamParticle(pRootMC))
        {
            std::cout << "   Test Beam Interaction Summary:" << std::endl;
            std::cout << std::fixed << std::setprecision(1);
            if (nTestBeamMCParticles)
            {
                std::cout << "   Good test beam particles: " << nTestBeamRecoParticles << " of " << nTestBeamMCParticles << " : "
                          << (100 * nTestBeamRecoParticles / static_cast<float>(nTestBeamMCParticles)) << "%" << std::endl;
            }
            if (nCosmicMCParticles)
            {
                std::cout << "   Matched cosmics: " << nCosmicRecoParticles << " of " << nCosmicMCParticles << " : "
                          << (100 * nCosmicRecoParticles / static_cast<float>(nCosmicMCParticles)) << "%" << std::endl;
            }
        }
        if (!this->GetUnmatchedReco().empty())
            std::cout << "   Unmatched reco: " << this->GetUnmatchedReco().size() << std::endl;
    }

    std::cout.precision(originalPrecision);
    std::cout.flags(originalFormat);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchInfo::GetRootMCParticles(MCParticleList &rootMCParticles) const
{
    for (auto iter = m_matches.begin(); iter != m_matches.end(); ++iter)
        rootMCParticles.emplace_back(iter->first);

    rootMCParticles.sort(LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::FillMCHierarchy(
    const MCParticleList &mcParticleList, const CaloHitList &caloHitList, const FoldingParameters &foldParameters, MCHierarchy &hierarchy)
{
    hierarchy.FillHierarchy(mcParticleList, caloHitList, foldParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::FillRecoHierarchy(const PfoList &pfoList, const FoldingParameters &foldParameters, RecoHierarchy &hierarchy)
{
    hierarchy.FillHierarchy(pfoList, foldParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchHierarchies(MatchInfo &matchInfo)
{
    matchInfo.Match();
}

//------------------------------------------------------------------------------------------------------------------------------------------
// private

void LArHierarchyHelper::GetMCPrimaries(const MCParticle *pRoot, MCParticleSet &primaries)
{
    try
    {
        MCParticleList visible;
        LArMCParticleHelper::GetFirstVisibleMCParticles(pRoot, visible);
        // Check if we still need to use sets here
        for (const MCParticle *const pMCParticle : visible)
            primaries.insert(pMCParticle);
    }
    catch (const StatusCodeException &)
    {
        if (pRoot->GetParticleId() != 111 && pRoot->GetParticleId() < 1e9)
            std::cout << "LArHierarchyHelper::MCHierarchy::FillHierarchy: MC particle with PDG code " << pRoot->GetParticleId()
                      << " at address " << pRoot << " has no associated primary particle" << std::endl;
    }
}

void LArHierarchyHelper::GetRecoPrimaries(const ParticleFlowObject *pRoot, PfoSet &primaries)
{
    if (LArPfoHelper::IsNeutrino(pRoot))
    {
        const PfoList &children{pRoot->GetDaughterPfoList()};
        // Check if we still need to use sets here
        for (const ParticleFlowObject *const pPfo : children)
            primaries.insert(pPfo);
    }
    else
    {
        // Might want different handling here for test beam and cosmic ray particles, but for now, just treat
        // a non-neutrino root node as something to be stored in the primaries set directly
        primaries.insert(pRoot);
    }
}

} // namespace lar_content
