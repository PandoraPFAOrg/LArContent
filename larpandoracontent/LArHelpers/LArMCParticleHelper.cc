/**
 *  @file   larpandoracontent/LArHelpers/LArMCParticleHelper.cc
 *
 *  @brief  Implementation of the lar monte carlo particle helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"

#include "Pandora/PdgTable.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <cstdlib>

namespace lar_content
{

using namespace pandora;

LArMCParticleHelper::PrimaryParameters::PrimaryParameters() :
    m_minPrimaryGoodHits(15),
    m_minHitsForGoodView(5),
    m_minPrimaryGoodViews(2),
    m_selectInputHits(true),
    m_maxPhotonPropagation(2.5f),
    m_minHitSharingFraction(0.9f),
    m_foldBackHierarchy(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::DoesPrimaryMeetCriteria(const MCParticle *const pMCParticle, std::function<bool(const MCParticle *const)> fCriteria)
{
    try
    {
        const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
        return fCriteria(pPrimaryMCParticle);
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::DoesLeadingMeetCriteria(const MCParticle *const pMCParticle, std::function<bool(const MCParticle *const)> fCriteria)
{
    try
    {
        const MCParticle *const pLeadingMCParticle = LArMCParticleHelper::GetLeadingMCParticle(pMCParticle);
        return fCriteria(pLeadingMCParticle);
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBeamNeutrinoFinalState(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && LArMCParticleHelper::IsNeutrino(pParentMCParticle));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsTriggeredBeamParticle(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && (nuance == 2001));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBeamParticle(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && ((nuance == 2000) || (nuance == 2001)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsLeadingBeamParticle(const MCParticle *const pMCParticle)
{
    // ATTN: Only the parent triggered beam particle has nuance code 2001
    const int parentNuance(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));
    return (LArMCParticleHelper::IsLeading(pMCParticle) && (parentNuance == 2001 || parentNuance == 2000));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsCosmicRay(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) &&
        ((nuance == 3000) || ((nuance == 0) && !LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsEM(const MCParticle *const pMCParticle)
{
    const int pdg{std::abs(pMCParticle->GetParticleId())};
    return (pdg == E_MINUS || pdg == PHOTON);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArMCParticleHelper::GetNuanceCode(const MCParticle *const pMCParticle)
{
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pMCParticle));
    if (pLArMCParticle)
        return pLArMCParticle->GetNuanceCode();

    std::cout << "LArMCParticleHelper::GetNuanceCode - Error: Can't cast to LArMCParticle" << std::endl;
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrino(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    if ((nuance == 0) || (nuance == 2000) || (nuance == 2001) || (nuance == 3000))
        return false;

    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));
    return ((NU_E == absoluteParticleId) || (NU_MU == absoluteParticleId) || (NU_TAU == absoluteParticleId));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsCCInteraction(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParent(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

    if (!LArMCParticleHelper::IsNeutrino(pParent))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pParent));
    if (!pLArMCParticle)
    {
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return pLArMCParticle->GetIsCC();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsPrimary(const pandora::MCParticle *const pMCParticle)
{
    try
    {
        return (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle);
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsLeading(const pandora::MCParticle *const pMCParticle)
{
    try
    {
        return (LArMCParticleHelper::GetLeadingMCParticle(pMCParticle) == pMCParticle);
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArMCParticleHelper::GetHierarchyTier(const pandora::MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;
    int tier(0);

    while (pParentMCParticle->GetParentList().empty() == false)
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        ++tier;
    }

    return tier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsVisible(const MCParticle *const pMCParticle)
{
    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));

    if ((E_MINUS == absoluteParticleId) || (MU_MINUS == absoluteParticleId) || (PI_PLUS == absoluteParticleId) || (K_PLUS == absoluteParticleId) ||
        (SIGMA_MINUS == absoluteParticleId) || (SIGMA_PLUS == absoluteParticleId) || (HYPERON_MINUS == absoluteParticleId) ||
        (PROTON == absoluteParticleId) || (PHOTON == absoluteParticleId) || (NEUTRON == absoluteParticleId))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTrueNeutrinos(const MCParticleList *const pMCParticleList, MCParticleVector &trueNeutrinos)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCParticle))
            trueNeutrinos.push_back(pMCParticle);
    }

    std::sort(trueNeutrinos.begin(), trueNeutrinos.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTrueTestBeamParticles(const MCParticleList *const pMCParticleList, MCParticleVector &trueTestBeamParticles)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsTriggeredBeamParticle(pMCParticle))
            trueTestBeamParticles.push_back(pMCParticle);
    }

    std::sort(trueTestBeamParticles.begin(), trueTestBeamParticles.end(), LArMCParticleHelper::SortByMomentum);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::GetTrueVertex(const MCParticleList *const pMCParticleList, CartesianVector &trueVertex)
{
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    if (!primaries.empty())
    {
        const MCParticle *primary{primaries.front()};
        const MCParticleList &parents{primary->GetParentList()};
        if (parents.size() == 1)
        {
            const MCParticle *trueNeutrino{parents.front()};
            if (LArMCParticleHelper::IsNeutrino(trueNeutrino))
            {
                trueVertex = primaries.front()->GetVertex();
                return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetParentMCParticle(const MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;

    while (pParentMCParticle->GetParentList().empty() == false)
        pParentMCParticle = LArMCParticleHelper::GetNextParentMCParticle(pParentMCParticle);

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetNextParentMCParticle(const MCParticle *const pMCParticle)
{
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, pMCParticle->GetParentList().size() != 1);
    return *(pMCParticle->GetParentList().begin());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetAllDescendentMCParticles(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &descendentMCParticleList)
{
    for (const MCParticle *pDaughterMCParticle : pMCParticle->GetDaughterList())
    {
        if (std::find(descendentMCParticleList.begin(), descendentMCParticleList.end(), pDaughterMCParticle) == descendentMCParticleList.end())
        {
            descendentMCParticleList.emplace_back(pDaughterMCParticle);
            LArMCParticleHelper::GetAllDescendentMCParticles(pDaughterMCParticle, descendentMCParticleList);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetAllDescendentMCParticles(const MCParticle *const pMCParticle, MCParticleList &descendentTrackParticles,
    MCParticleList &leadingShowerParticles, MCParticleList &leadingNeutrons)
{
    for (const MCParticle *pDaughterMCParticle : pMCParticle->GetDaughterList())
    {
        if (std::find(descendentTrackParticles.begin(), descendentTrackParticles.end(), pDaughterMCParticle) == descendentTrackParticles.end())
        {
            const int pdg{std::abs(pDaughterMCParticle->GetParticleId())};
            if (pdg == E_MINUS || pdg == PHOTON)
            {
                leadingShowerParticles.emplace_back(pDaughterMCParticle);
            }
            else if (pdg == NEUTRON)
            {
                leadingNeutrons.emplace_back(pDaughterMCParticle);
            }
            else
            {
                descendentTrackParticles.emplace_back(pDaughterMCParticle);
                LArMCParticleHelper::GetAllDescendentMCParticles(pDaughterMCParticle, descendentTrackParticles, leadingShowerParticles, leadingNeutrons);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetAllAncestorMCParticles(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &ancestorMCParticleList)
{
    const MCParticleList &parentMCParticleList = pMCParticle->GetParentList();
    if (parentMCParticleList.empty())
        return;
    if (parentMCParticleList.size() != 1)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const MCParticle *pParentMCParticle = *parentMCParticleList.begin();
    if (std::find(ancestorMCParticleList.begin(), ancestorMCParticleList.end(), pParentMCParticle) == ancestorMCParticleList.end())
    {
        ancestorMCParticleList.push_back(pParentMCParticle);
        LArMCParticleHelper::GetAllAncestorMCParticles(pParentMCParticle, ancestorMCParticleList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPrimaryMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcPrimaryVector)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsPrimary(pMCParticle))
            mcPrimaryVector.push_back(pMCParticle);
    }

    std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetLeadingMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcLeadingVector)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const bool isBeamParticle(LArMCParticleHelper::IsBeamParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));

        if ((isBeamParticle && LArMCParticleHelper::IsLeading(pMCParticle)) || (!isBeamParticle && LArMCParticleHelper::IsPrimary(pMCParticle)))
        {
            mcLeadingVector.push_back(pMCParticle);
        }
    }

    std::sort(mcLeadingVector.begin(), mcLeadingVector.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetFirstVisibleMCParticles(const MCParticle *const pRoot, MCParticleList &visibleParticleList)
{
    if (LArMCParticleHelper::IsVisible(pRoot))
    {
        visibleParticleList.emplace_back(pRoot);
    }
    else
    {
        for (const MCParticle *const pMCParticle : pRoot->GetDaughterList())
            LArMCParticleHelper::GetFirstVisibleMCParticles(pMCParticle, visibleParticleList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetPrimaryMCParticle(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcVector.push_back(pParentMCParticle);
    }

    // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
    for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(), iterEnd = mcVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        if (LArMCParticleHelper::IsVisible(pNextParticle))
            return pNextParticle;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetLeadingMCParticle(const MCParticle *const pMCParticle, const int hierarchyTierLimit)
{
    // ATTN: If not beam particle return primary particle
    const bool isBeamParticle(LArMCParticleHelper::IsBeamParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));

    if (!isBeamParticle)
        return LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);

    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcVector.push_back(pParentMCParticle);
    }

    int hierarchyTier(0);
    const MCParticle *pLeadingMCParticle(nullptr);

    // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
    for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(), iterEnd = mcVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        // ATTN: Squash any invisible particles (e.g. pizero)
        if (!LArMCParticleHelper::IsVisible(pNextParticle))
            continue;

        if (hierarchyTier <= hierarchyTierLimit)
            pLeadingMCParticle = pNextParticle;

        hierarchyTier++;
    }

    if (!pLeadingMCParticle)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pLeadingMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCPrimaryMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        try
        {
            const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
            mcPrimaryMap[pMCParticle] = pPrimaryMCParticle;
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCLeadingMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcLeadingMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        try
        {
            const MCParticle *const pLeadingMCParticle = LArMCParticleHelper::GetLeadingMCParticle(pMCParticle);
            mcLeadingMap[pMCParticle] = pLeadingMCParticle;
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCToSelfMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcToSelfMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        mcToSelfMap[pMCParticle] = pMCParticle;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCToHitsMap(const CaloHitList &caloHitList2D, MCContributionMap &mcToHitsMap, const bool allowSharing)
{
    for (const CaloHit *const pCaloHit : caloHitList2D)
    {
        try
        {
            if (!allowSharing)
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                mcToHitsMap[pMC].emplace_back(pCaloHit);
            }
            else
            {
                const MCParticleWeightMap &mcWeightMap{pCaloHit->GetMCParticleWeightMap()};

                for (const auto &[pMC, _] : mcWeightMap)
                    mcToHitsMap[pMC].emplace_back(pCaloHit);
            }
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCToHitsMap(const CaloHitList *const pCaloHitList2D, const MCParticleList *const pMCParticleList,
    LArMCParticleHelper::MCContributionMap &mcToHitsMap)
{
    PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList2D, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCParticleList &mcHierarchy)
{
    try
    {
        for (const auto &[mc, hits] : mcToHitsMap)
        {
            mcHierarchy.emplace_back(mc);
            LArMCParticleHelper::GetAllAncestorMCParticles(mc, mcHierarchy);
        }
    }
    catch (const StatusCodeException &e)
    {
        throw;
    }

    // Move the neutrino to the front of the list
    auto pivot{
        std::find_if(mcHierarchy.begin(), mcHierarchy.end(), [](const MCParticle *mc) -> bool { return LArMCParticleHelper::IsNeutrino(mc); })};
    if (pivot != mcHierarchy.end())
        std::rotate(mcHierarchy.begin(), pivot, std::next(pivot));
    else
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetMainMCParticle(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
    return MCParticleHelper::GetMainMCParticle(&clusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::SortByMomentum(const MCParticle *const pLhs, const MCParticle *const pRhs)
{
    // Sort by momentum (prefer higher momentum)
    const float momentumLhs(pLhs->GetMomentum().GetMagnitudeSquared());
    const float momentumRhs(pRhs->GetMomentum().GetMagnitudeSquared());

    if (std::fabs(momentumLhs - momentumRhs) > std::numeric_limits<float>::epsilon())
        return (momentumLhs > momentumRhs);

    // Sort by energy (prefer lighter particles)
    if (std::fabs(pLhs->GetEnergy() - pRhs->GetEnergy()) > std::numeric_limits<float>::epsilon())
        return (pLhs->GetEnergy() < pRhs->GetEnergy());

    // Sort by PDG code (prefer smaller numbers)
    if (pLhs->GetParticleId() != pRhs->GetParticleId())
        return (pLhs->GetParticleId() < pRhs->GetParticleId());

    // Sort by vertex position (tie-breaker)
    const float positionLhs(pLhs->GetVertex().GetMagnitudeSquared());
    const float positionRhs(pRhs->GetVertex().GetMagnitudeSquared());

    return (positionLhs < positionRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const MCRelationMap &mcToTargetMCMap,
    CaloHitToMCMap &hitToMCMap, MCContributionMap &mcToTrueHitListMap)
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const MCParticle *pTargetParticle(pHitParticle);

            // ATTN Do not map back to target if mc to primary mc map or mc to self map not provided
            if (!mcToTargetMCMap.empty())
            {
                MCRelationMap::const_iterator mcIter = mcToTargetMCMap.find(pHitParticle);

                if (mcToTargetMCMap.end() == mcIter)
                    continue;

                pTargetParticle = mcIter->second;
            }

            mcToTrueHitListMap[pTargetParticle].push_back(pCaloHit);
            hitToMCMap[pCaloHit] = pTargetParticle;
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectReconstructableMCParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList,
    const PrimaryParameters &parameters, std::function<bool(const MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Obtain map: [mc particle -> target mc particle]
    LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
    parameters.m_foldBackHierarchy ? LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToTargetMCMap)
                                   : LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToTargetMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    // Unless selectInputHits == false
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToTargetMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    // Obtain maps: [hit -> target mc particle], [target mc particle -> list of hits]
    CaloHitToMCMap trueHitToTargetMCMap;
    MCContributionMap targetMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToTargetMCMap, trueHitToTargetMCMap, targetMCToTrueHitListMap);

    // Obtain vector: target mc particles
    MCParticleVector targetMCVector;
    if (parameters.m_foldBackHierarchy)
    {
        LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, targetMCVector);
    }
    else
    {
        std::copy(pMCParticleList->begin(), pMCParticleList->end(), std::back_inserter(targetMCVector));
    }

    // Select MCParticles matching criteria
    MCParticleVector candidateTargets;
    LArMCParticleHelper::SelectParticlesMatchingCriteria(targetMCVector, fCriteria, candidateTargets, parameters, false);

    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(candidateTargets, targetMCToTrueHitListMap, mcToTargetMCMap, parameters, selectedMCParticlesToHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectReconstructableTestBeamHierarchyMCParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList,
    const PrimaryParameters &parameters, std::function<bool(const MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Obtain map: [mc particle -> target mc particle]
    LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
    parameters.m_foldBackHierarchy ? LArMCParticleHelper::GetMCLeadingMap(pMCParticleList, mcToTargetMCMap)
                                   : LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToTargetMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    // Unless selectInputHits == false
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToTargetMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    // Obtain maps: [hit -> target mc particle], [target mc particle -> list of hits]
    CaloHitToMCMap trueHitToTargetMCMap;
    MCContributionMap targetMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToTargetMCMap, trueHitToTargetMCMap, targetMCToTrueHitListMap);

    // Obtain vector: target mc particles
    MCParticleVector targetMCVector;
    if (parameters.m_foldBackHierarchy)
    {
        LArMCParticleHelper::GetLeadingMCParticleList(pMCParticleList, targetMCVector);
    }
    else
    {
        std::copy(pMCParticleList->begin(), pMCParticleList->end(), std::back_inserter(targetMCVector));
    }

    // Select MCParticles matching criteria
    MCParticleVector candidateTargets;
    LArMCParticleHelper::SelectParticlesMatchingCriteria(targetMCVector, fCriteria, candidateTargets, parameters, true);

    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(candidateTargets, targetMCToTrueHitListMap, mcToTargetMCMap, parameters, selectedMCParticlesToHitsMap);

    // ATTN: The parent particle must be in the hierarchy map, event if not reconstructable
    for (const MCParticle *const pMCParticle : candidateTargets)
    {
        if (!LArMCParticleHelper::IsBeamParticle(pMCParticle))
            continue;

        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
        if (selectedMCParticlesToHitsMap.find(pParentMCParticle) == selectedMCParticlesToHitsMap.end())
        {
            CaloHitList caloHitList;
            selectedMCParticlesToHitsMap.insert(MCContributionMap::value_type(pParentMCParticle, caloHitList));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap,
    PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy)
{
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(
        pfoList, MCContributionMapVector({selectedMCParticleToHitsMap}), pfoToReconstructable2DHitsMap, foldBackHierarchy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const PfoList &pfoList,
    const MCContributionMap &selectedMCParticleToHitsMap, PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy)
{
    LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(
        pfoList, MCContributionMapVector({selectedMCParticleToHitsMap}), pfoToReconstructable2DHitsMap, foldBackHierarchy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        CaloHitList pfoHitList;
        LArMCParticleHelper::CollectReconstructable2DHits(pPfo, selectedMCParticleToHitsMaps, pfoHitList, foldBackHierarchy);

        if (!pfoToReconstructable2DHitsMap.insert(PfoContributionMap::value_type(pPfo, pfoHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const PfoList &pfoList,
    const MCContributionMapVector &selectedMCParticleToHitsMaps, PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        CaloHitList pfoHitList;
        LArMCParticleHelper::CollectReconstructableTestBeamHierarchy2DHits(pPfo, selectedMCParticleToHitsMaps, pfoHitList, foldBackHierarchy);

        if (!pfoToReconstructable2DHitsMap.insert(PfoContributionMap::value_type(pPfo, pfoHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(const PfoContributionMap &pfoToReconstructable2DHitsMap,
    const MCContributionMapVector &selectedMCParticleToHitsMaps, PfoToMCParticleHitSharingMap &pfoToMCParticleHitSharingMap,
    MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap)
{
    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoToReconstructable2DHitsMap)
        sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : sortedPfos)
    {
        for (const MCContributionMap &mcParticleToHitsMap : selectedMCParticleToHitsMaps)
        {
            MCParticleVector sortedMCParticles;
            for (const auto &mapEntry : mcParticleToHitsMap)
                sortedMCParticles.push_back(mapEntry.first);
            std::sort(sortedMCParticles.begin(), sortedMCParticles.end(), PointerLessThan<MCParticle>());

            for (const MCParticle *const pMCParticle : sortedMCParticles)
            {
                // Add map entries for this Pfo & MCParticle if required
                if (pfoToMCParticleHitSharingMap.find(pPfo) == pfoToMCParticleHitSharingMap.end())
                    if (!pfoToMCParticleHitSharingMap.insert(PfoToMCParticleHitSharingMap::value_type(pPfo, MCParticleToSharedHitsVector())).second)
                        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT); // ATTN maybe overkill

                if (mcParticleToPfoHitSharingMap.find(pMCParticle) == mcParticleToPfoHitSharingMap.end())
                    if (!mcParticleToPfoHitSharingMap.insert(MCParticleToPfoHitSharingMap::value_type(pMCParticle, PfoToSharedHitsVector())).second)
                        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // Check this Pfo & MCParticle pairing hasn't already been checked
                MCParticleToSharedHitsVector &mcHitPairs(pfoToMCParticleHitSharingMap.at(pPfo));
                PfoToSharedHitsVector &pfoHitPairs(mcParticleToPfoHitSharingMap.at(pMCParticle));

                if (std::any_of(mcHitPairs.begin(), mcHitPairs.end(),
                        [&](const MCParticleCaloHitListPair &pair) { return (pair.first == pMCParticle); }))
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                if (std::any_of(pfoHitPairs.begin(), pfoHitPairs.end(), [&](const PfoCaloHitListPair &pair) { return (pair.first == pPfo); }))
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // Add records to maps if there are any shared hits
                const CaloHitList sharedHits(
                    LArMCParticleHelper::GetSharedHits(pfoToReconstructable2DHitsMap.at(pPfo), mcParticleToHitsMap.at(pMCParticle)));

                if (!sharedHits.empty())
                {
                    mcHitPairs.push_back(MCParticleCaloHitListPair(pMCParticle, sharedHits));
                    pfoHitPairs.push_back(PfoCaloHitListPair(pPfo, sharedHits));

                    std::sort(mcHitPairs.begin(), mcHitPairs.end(),
                        [](const MCParticleCaloHitListPair &a, const MCParticleCaloHitListPair &b) -> bool
                        {
                            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size()
                                                                         : LArMCParticleHelper::SortByMomentum(a.first, b.first));
                        });

                    std::sort(pfoHitPairs.begin(), pfoHitPairs.end(),
                        [](const PfoCaloHitListPair &a, const PfoCaloHitListPair &b) -> bool {
                            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size()
                                                                         : LArPfoHelper::SortByNHits(a.first, b.first));
                        });
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetClusterToReconstructable2DHitsMap(const pandora::ClusterList &clusterList,
    const MCContributionMap &selectedMCToHitsMap, ClusterContributionMap &clusterToReconstructable2DHitsMap)
{
    LArMCParticleHelper::GetClusterToReconstructable2DHitsMap(clusterList, MCContributionMapVector({selectedMCToHitsMap}), clusterToReconstructable2DHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetClusterToReconstructable2DHitsMap(const pandora::ClusterList &clusterList,
    const MCContributionMapVector &selectedMCToHitsMaps, ClusterContributionMap &clusterToReconstructable2DHitsMap)
{
    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList caloHitList;
        LArMCParticleHelper::CollectReconstructable2DHits(pCluster, selectedMCToHitsMaps, caloHitList);

        if (!clusterToReconstructable2DHitsMap.insert(ClusterContributionMap::value_type(pCluster, caloHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBremsstrahlung(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_E_BREM:
        case MC_PROC_MU_BREM:
        case MC_PROC_HAD_BREM:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsCapture(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_MU_MINUS_CAPTURE_AT_REST:
        case MC_PROC_N_CAPTURE:
        case MC_PROC_CHIPS_NUCLEAR_CAPTURE_AT_REST:
        case MC_PROC_HAD_FRITIOF_CAPTURE_AT_REST:
        case MC_PROC_HAD_BERTINI_CAPTURE_AT_REST:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsDecay(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_DECAY:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsElasticScatter(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_COULOMB_SCAT:
        case MC_PROC_NEUTRON_INELASTIC:
        case MC_PROC_HAD_ELASTIC:
        case MC_PROC_RAYLEIGH:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsInelasticScatter(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_COMPT:
        case MC_PROC_PHOTON_INELASTIC:
        case MC_PROC_PROTON_INELASTIC:
        case MC_PROC_PI_PLUS_INELASTIC:
        case MC_PROC_PI_MINUS_INELASTIC:
        case MC_PROC_ION_INELASTIC:
        case MC_PROC_HE3_INELASTIC:
        case MC_PROC_ALPHA_INELASTIC:
        case MC_PROC_ANTI_HE3_INELASTIC:
        case MC_PROC_ANTI_ALPHA_INELASTIC:
        case MC_PROC_ANTI_DEUTERON_INELASTIC:
        case MC_PROC_ANTI_NEUTRON_INELASTIC:
        case MC_PROC_ANTI_PROTON_INELASTIC:
        case MC_PROC_ANTI_TRITON_INELASTIC:
        case MC_PROC_DEUTERON_INELASTIC:
        case MC_PROC_KAON_PLUS_INELASTIC:
        case MC_PROC_KAON_MINUS_INELASTIC:
        case MC_PROC_LAMBDA_INELASTIC:
        case MC_PROC_TRITON_INELASTIC:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsIonisation(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_E_IONI:
        case MC_PROC_MU_IONI:
        case MC_PROC_HAD_IONI:
        case MC_PROC_ION_IONI:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNuclear(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_ELECTRON_NUCLEAR:
        case MC_PROC_PHOTON_NUCLEAR:
        case MC_PROC_MU_NUCLEAR:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsPairProduction(const MCParticle *const pMCParticle)
{
    const LArMCParticle *pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
    if (!pLArMCParticle)
        return false;

    switch (pLArMCParticle->GetProcess())
    {
        case MC_PROC_MU_PAIR_PROD:
        case MC_PROC_HAD_PAIR_PROD:
            return true;
        default:
            return false;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::CausesShower(const MCParticle *const pMC, int nDescendantElectrons)
{
    if (nDescendantElectrons > 1)
        return true;

    if (std::abs(pMC->GetParticleId()) == E_MINUS)
        nDescendantElectrons++; // Including the parent particle, ie. the first in the recursion, as a descendant

    for (const MCParticle *pChildMC : pMC->GetDaughterList())
    {
        if (LArMCParticleHelper::CausesShower(pChildMC, nDescendantElectrons))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArMCParticleHelper::GetSharedHits(const CaloHitList &hitListA, const CaloHitList &hitListB)
{
    CaloHitList sharedHits;

    for (const CaloHit *const pCaloHit : hitListA)
    {
        if (std::find(hitListB.begin(), hitListB.end(), pCaloHit) != hitListB.end())
            sharedHits.push_back(pCaloHit);
    }

    return sharedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::AreTopologicallyContinuous(const MCParticle *const pMCParent, const MCParticle *const pMCChild, const float cosAngleTolerance)
{
    CartesianVector childDirection{pMCChild->GetEndpoint() - pMCChild->GetVertex()};
    if (childDirection.GetMagnitude() < std::numeric_limits<float>::epsilon())
        return true;
    childDirection = childDirection.GetUnitVector();

    const MCParticle *pMCUpstream{pMCParent};
    while (true)
    {
        CartesianVector parentDirection{pMCUpstream->GetEndpoint() - pMCUpstream->GetVertex()};
        if (parentDirection.GetMagnitude() > std::numeric_limits<float>::epsilon())
        {
            parentDirection = parentDirection.GetUnitVector();
            return parentDirection.GetDotProduct(childDirection) >= cosAngleTolerance;
        }
        else
        {
            const MCParticleList &parentList{pMCUpstream->GetParentList()};
            const size_t size{parentList.size()};
            if (size == 1)
                pMCUpstream = parentList.front();
            else if (size == 0)
                return true;
            else
                return false;
        }
    }

    return false;
}

// private
//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructable2DHits(const ParticleFlowObject *const pPfo,
    const MCContributionMapVector &selectedMCParticleToHitsMaps, CaloHitList &reconstructableCaloHitList2D, const bool foldBackHierarchy)
{

    PfoList pfoList;

    // If foldBackHierarchy collect all 2D calo hits in pfo hierarchy
    // else collect hits directly belonging to pfo
    if (foldBackHierarchy)
    {
        LArPfoHelper::GetAllDownstreamPfos(pPfo, pfoList);
    }
    else
    {
        pfoList.push_back(pPfo);
    }

    LArMCParticleHelper::CollectReconstructable2DHits(pfoList, selectedMCParticleToHitsMaps, reconstructableCaloHitList2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructableTestBeamHierarchy2DHits(const ParticleFlowObject *const pPfo,
    const MCContributionMapVector &selectedMCParticleToHitsMaps, CaloHitList &reconstructableCaloHitList2D, const bool foldBackHierarchy)
{

    PfoList pfoList;
    // If foldBackHierarchy collect all 2D calo hits in pfo hierarchy
    // else collect hits directly belonging to pfo
    if (foldBackHierarchy)
    {
        // ATTN: Only collect downstream pfos for daughter test beam particles & cosmics
        if (pPfo->GetParentPfoList().empty() && LArPfoHelper::IsTestBeam(pPfo))
        {
            pfoList.push_back(pPfo);
        }
        else
        {
            LArPfoHelper::GetAllDownstreamPfos(pPfo, pfoList);
        }
    }
    else
    {
        pfoList.push_back(pPfo);
    }

    LArMCParticleHelper::CollectReconstructable2DHits(pfoList, selectedMCParticleToHitsMaps, reconstructableCaloHitList2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructable2DHits(
    const PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps, CaloHitList &reconstructableCaloHitList2D)
{
    CaloHitList caloHitList2D;
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_U, caloHitList2D);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_V, caloHitList2D);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_W, caloHitList2D);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_U, caloHitList2D); // TODO check isolated usage throughout
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_V, caloHitList2D);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_W, caloHitList2D);

    // Filter for only reconstructable hits
    for (const CaloHit *const pCaloHit : caloHitList2D)
    {
        bool isTargetHit(false);
        for (const MCContributionMap &mcParticleToHitsMap : selectedMCParticleToHitsMaps)
        {
            // ATTN This map is unordered, but this does not impact search for specific target hit
            for (const MCContributionMap::value_type &mapEntry : mcParticleToHitsMap)
            {
                if (std::find(mapEntry.second.begin(), mapEntry.second.end(), pCaloHit) != mapEntry.second.end())
                {
                    isTargetHit = true;
                    break;
                }
            }
            if (isTargetHit)
                break;
        }

        if (isTargetHit)
            reconstructableCaloHitList2D.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructable2DHits(const pandora::Cluster *const pCluster,
    const MCContributionMapVector &selectedMCToHitsMaps, pandora::CaloHitList &reconstructableCaloHitList2D)
{
    const CaloHitList &isolatedCaloHitList{pCluster->GetIsolatedCaloHitList()};
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    for (const CaloHit *pCaloHit : isolatedCaloHitList)
        caloHitList.push_back(pCaloHit);

    // Filter for only reconstructable hits
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        bool isTargetHit{false};
        for (const MCContributionMap &mcParticleToHitsMap : selectedMCToHitsMaps)
        { // ATTN This map is unordered, but this does not impact search for specific target hit
            for (const MCContributionMap::value_type &mapEntry : mcParticleToHitsMap)
            {
                if (std::find(mapEntry.second.begin(), mapEntry.second.end(), pCaloHit) != mapEntry.second.end())
                {
                    isTargetHit = true;
                    break;
                }
            }
            if (isTargetHit)
                break;
        }

        if (isTargetHit)
            reconstructableCaloHitList2D.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
    CaloHitList &selectedCaloHitList, const bool selectInputHits, const float maxPhotonPropagation)
{
    if (!selectInputHits)
    {
        selectedCaloHitList.insert(selectedCaloHitList.end(), pCaloHitList->begin(), pCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToTargetMCMap.find(pHitParticle);

            if (mcToTargetMCMap.end() == mcIter)
                continue;

            // ATTN With folding on or off, still require primary particle to review hierarchy details
            const MCParticle *const pPrimaryParticle = LArMCParticleHelper::GetPrimaryMCParticle(pHitParticle);

            if (PassMCParticleChecks(pPrimaryParticle, pPrimaryParticle, pHitParticle, maxPhotonPropagation))
                selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsDescendentOf(const MCParticle *const pMCParticle, const int pdg, const bool isChargeSensitive)
{
    const MCParticle *pCurrentParticle = pMCParticle;
    while (!pCurrentParticle->GetParentList().empty())
    {
        if (pCurrentParticle->GetParentList().size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const MCParticle *pParent = *(pCurrentParticle->GetParentList().begin());
        const bool found{isChargeSensitive ? pParent->GetParticleId() == pdg : std::abs(pParent->GetParticleId()) == std::abs(pdg)};
        if (found)
            return true;
        pCurrentParticle = pParent;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetBreadthFirstHierarchyRepresentation(const MCParticle *const pMCParticle, MCParticleList &mcParticleList)
{
    const MCParticle *const pRoot{LArMCParticleHelper::GetParentMCParticle(pMCParticle)};
    MCParticleList queue;
    mcParticleList.emplace_back(pRoot);
    queue.emplace_back(pRoot);

    while (!queue.empty())
    {
        const MCParticleList &daughters{queue.front()->GetDaughterList()};
        queue.pop_front();
        for (const MCParticle *pDaughter : daughters)
        {
            mcParticleList.emplace_back(pDaughter);
            queue.emplace_back(pDaughter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectParticlesByHitCount(const MCParticleVector &candidateTargets, const MCContributionMap &mcToTrueHitListMap,
    const MCRelationMap &mcToTargetMCMap, const PrimaryParameters &parameters, MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Apply restrictions on the number of good hits associated with the MCParticles
    for (const MCParticle *const pMCTarget : candidateTargets)
    {
        MCContributionMap::const_iterator trueHitsIter = mcToTrueHitListMap.find(pMCTarget);
        if (mcToTrueHitListMap.end() == trueHitsIter)
            continue;

        const CaloHitList &caloHitList(trueHitsIter->second);

        // Remove shared hits where target particle deposits below threshold energy fraction
        CaloHitList goodCaloHitList;
        LArMCParticleHelper::SelectGoodCaloHits(
            &caloHitList, mcToTargetMCMap, goodCaloHitList, parameters.m_selectInputHits, parameters.m_minHitSharingFraction);

        if (goodCaloHitList.size() < parameters.m_minPrimaryGoodHits)
            continue;

        unsigned int nGoodViews(0);
        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, goodCaloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, goodCaloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, goodCaloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (nGoodViews < parameters.m_minPrimaryGoodViews)
            continue;

        if (!selectedMCParticlesToHitsMap.insert(MCContributionMap::value_type(pMCTarget, caloHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
    CaloHitList &selectedGoodCaloHitList, const bool selectInputHits, const float minHitSharingFraction)
{
    if (!selectInputHits)
    {
        selectedGoodCaloHitList.insert(selectedGoodCaloHitList.end(), pSelectedCaloHitList->begin(), pSelectedCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pSelectedCaloHitList)
    {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap())
            mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        MCParticleWeightMap targetWeightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToTargetMCMap.find(pMCParticle);

            if (mcToTargetMCMap.end() != mcIter)
                targetWeightMap[mcIter->second] += weight;
        }

        MCParticleVector mcTargetVector;
        for (const auto &mapEntry : targetWeightMap)
            mcTargetVector.push_back(mapEntry.first);
        std::sort(mcTargetVector.begin(), mcTargetVector.end(), PointerLessThan<MCParticle>());

        const MCParticle *pBestTargetParticle(nullptr);
        float bestTargetWeight(0.f), targetWeightSum(0.f);

        for (const MCParticle *const pTargetMCParticle : mcTargetVector)
        {
            const float targetWeight(targetWeightMap.at(pTargetMCParticle));
            targetWeightSum += targetWeight;

            if (targetWeight > bestTargetWeight)
            {
                bestTargetWeight = targetWeight;
                pBestTargetParticle = pTargetMCParticle;
            }
        }

        if (!pBestTargetParticle || (targetWeightSum < std::numeric_limits<float>::epsilon()) || ((bestTargetWeight / targetWeightSum) < minHitSharingFraction))
            continue;

        selectedGoodCaloHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectParticlesMatchingCriteria(const MCParticleVector &inputMCParticles,
    std::function<bool(const MCParticle *const)> fCriteria, MCParticleVector &selectedParticles, const PrimaryParameters &parameters, const bool isTestBeam)
{
    for (const MCParticle *const pMCParticle : inputMCParticles)
    {
        if (parameters.m_foldBackHierarchy)
        {
            if (!fCriteria(pMCParticle))
                continue;
        }
        else
        {
            if (isTestBeam)
            {
                if (!LArMCParticleHelper::DoesLeadingMeetCriteria(pMCParticle, fCriteria))
                    continue;
            }
            else
            {
                if (!LArMCParticleHelper::DoesPrimaryMeetCriteria(pMCParticle, fCriteria))
                    continue;
            }
        }
        selectedParticles.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
    const MCParticle *const pHitMCParticle, const float maxPhotonPropagation)
{
    if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
        return false;

    if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != pOriginalPrimary->GetParticleId()) &&
        (E_MINUS != std::abs(pOriginalPrimary->GetParticleId())))
    {
        if ((pThisMCParticle->GetEndpoint() - pThisMCParticle->GetVertex()).GetMagnitude() > maxPhotonPropagation)
            return false;
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (PassMCParticleChecks(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle, maxPhotonPropagation))
            return true;
    }

    return false;
}

} // namespace lar_content
