/**
 *  @file   larpandoracontent/LArMetrics/HierarchyValidationTool.cc
 *
 *  @brief  Implementation of the hierarchy validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArMetrics/HierarchyValidationTool.h"

using namespace pandora;

namespace lar_content
{

HierarchyValidationTool::HierarchyValidationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu,
    [[maybe_unused]] const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    HierarchyTreeVars hierarchyTreeVars;
    hierarchyTreeVars.m_run = this->GetPandora().GetRun();
    hierarchyTreeVars.m_subrun = this->GetPandora().GetSubrun();
    hierarchyTreeVars.m_event = this->GetPandora().GetEvent();

    // Build true hierarchy out of 'target' particles
    Hierarchy hierarchy;
    this->BuildVisibleHierarchy(pMCNu, pMCNu, targetMC, 1, hierarchy);

    // Fill tree variables for each particle
    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMC(targetMC.at(i));
        const Pfo *const pBestMatch(bestRecoMatch.at(i));

        this->FillTrueVariables(pMC, targetMC, hierarchy, hierarchyTreeVars);

        if (pBestMatch)
        {
            this->FillRecoVariables(pBestMatch, bestRecoMatch, hierarchyTreeVars);
        }
        else
        {
            hierarchyTreeVars.m_recoTier.push_back(m_invalidInt);
            hierarchyTreeVars.m_recoParentIndex.push_back(m_invalidInt);
        }
    }

    this->FillTree(hierarchyTreeVars);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationTool::BuildVisibleHierarchy(const MCParticle *const pMCParticle, const MCParticle *const pMCParent,
    const MCParticleVector &targetMC, const int childTier, Hierarchy &hierarchy)
{
    for (const MCParticle *const pMCChild : pMCParticle->GetDaughterList())
    {
        // If child is not a target, then leapfrog
        if (std::find(targetMC.begin(), targetMC.end(), pMCChild) == targetMC.end())
        {
            this->BuildVisibleHierarchy(pMCChild, pMCParent, targetMC, childTier, hierarchy);
        }
        else
        {
            hierarchy.insert(std::make_pair(pMCChild, std::make_pair(pMCParent, childTier)));
            this->BuildVisibleHierarchy(pMCChild, pMCChild, targetMC, (childTier + 1), hierarchy);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationTool::FillTrueVariables(
    const MCParticle *const pMC, const MCParticleVector &targetMC, const Hierarchy &hierarchy, HierarchyTreeVars &hierarchyTreeVars)
{
    const auto hierarchyIter(hierarchy.find(pMC));

    if (hierarchyIter == hierarchy.end())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const int trueTier(hierarchyIter->second.second);
    hierarchyTreeVars.m_trueTier.push_back(trueTier);

    // Find index of parent
    if (trueTier == 1)
    {
        // Set neutrino index to be -1
        hierarchyTreeVars.m_trueParentIndex.push_back(-1);
    }
    else
    {
        const MCParticle *const pMCParent(hierarchyIter->second.first);
        const auto targetMCIter(std::find(targetMC.begin(), targetMC.end(), pMCParent));

        if (targetMCIter == targetMC.end())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        hierarchyTreeVars.m_trueParentIndex.push_back(targetMCIter - targetMC.begin());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationTool::FillRecoVariables(const Pfo *const pBestMatch, const PfoVector &bestMatches, HierarchyTreeVars &hierarchyTreeVars)
{
    const int recoTier(LArPfoHelper::GetHierarchyTier(pBestMatch));
    hierarchyTreeVars.m_recoTier.push_back(recoTier);

    // Find index of parent
    if (recoTier == 1)
    {
        // Set neutrino index to be -1
        hierarchyTreeVars.m_recoParentIndex.push_back(-1);
    }
    else
    {
        const Pfo *const pParentPfo(pBestMatch->GetParentPfoList().front());
        const auto bestMatchesIter(std::find(bestMatches.begin(), bestMatches.end(), pParentPfo));

        // If reco parent is not in the best-match list, the particle has been split
        // I don't want this to feed into a hierarchy building 'failure' so mark as -1 and remove from downstream metrics
        if (bestMatchesIter == bestMatches.end())
            hierarchyTreeVars.m_recoParentIndex.push_back(-1);
        else
            hierarchyTreeVars.m_recoParentIndex.push_back(bestMatchesIter - bestMatches.begin());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationTool::FillTree([[maybe_unused]] HierarchyTreeVars &hierarchyTreeVars)
{
#ifdef MONITORING
    IntVector &trueTier = hierarchyTreeVars.m_trueTier;
    IntVector &trueParentIndex = hierarchyTreeVars.m_trueParentIndex;
    IntVector &recoTier = hierarchyTreeVars.m_recoTier;
    IntVector &recoParentIndex = hierarchyTreeVars.m_recoParentIndex;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HierarchyTree", "MC_HierarchyTier", &trueTier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HierarchyTree", "MC_ParentIndex", &trueParentIndex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HierarchyTree", "BM_HierarchyTier", &recoTier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HierarchyTree", "BM_ParentIndex", &recoParentIndex));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "HierarchyTree"));
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return BaseValidationTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
