/**
 *  @file   larpandoracontent/LArHelpers/LArMCParticleHelper.h
 *
 *  @brief  Header file for the lar monte carlo particle helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_PARTICLE_HELPER_H
#define LAR_MC_PARTICLE_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include <functional>
#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  LArMCParticleHelper class
 */
class LArMCParticleHelper
{
public:
    typedef std::unordered_map<const pandora::MCParticle *, const pandora::MCParticle *> MCRelationMap;
    typedef std::unordered_map<const pandora::MCParticle *, int> MCParticleIntMap;

    typedef std::unordered_map<const pandora::MCParticle *, const pandora::ParticleFlowObject *> MCToPfoMap;

    typedef std::unordered_map<const pandora::CaloHit *, const pandora::MCParticle *> CaloHitToMCMap;
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::ParticleFlowObject *> CaloHitToPfoMap;

    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCContributionMap;
    typedef std::vector<MCContributionMap> MCContributionMapVector;

    typedef std::unordered_map<const pandora::Cluster *, pandora::CaloHitList> ClusterContributionMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, pandora::CaloHitList> PfoContributionMap;
    typedef std::unordered_map<const pandora::MCParticle *, PfoContributionMap> MCToPfoMatchingMap;

    typedef std::pair<const pandora::MCParticle *, pandora::CaloHitList> MCParticleCaloHitListPair;
    typedef std::pair<const pandora::ParticleFlowObject *, pandora::CaloHitList> PfoCaloHitListPair;
    typedef std::pair<const pandora::Cluster *, pandora::CaloHitList> ClusterCaloHitListPair;

    typedef std::vector<MCParticleCaloHitListPair> MCParticleToSharedHitsVector;
    typedef std::vector<PfoCaloHitListPair> PfoToSharedHitsVector;

    typedef std::map<const pandora::ParticleFlowObject *, MCParticleToSharedHitsVector> PfoToMCParticleHitSharingMap;
    typedef std::map<const pandora::MCParticle *, PfoToSharedHitsVector> MCParticleToPfoHitSharingMap;

    /**
     *  @brief   PrimaryParameters class
     */
    class PrimaryParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        PrimaryParameters();

        unsigned int m_minPrimaryGoodHits;  ///< the minimum number of primary good Hits
        unsigned int m_minHitsForGoodView;  ///< the minimum number of Hits for a good view
        unsigned int m_minPrimaryGoodViews; ///< the minimum number of primary good views
        bool m_selectInputHits;             ///< whether to select input hits
        float m_maxPhotonPropagation;       ///< the maximum photon propagation length
        float m_minHitSharingFraction;      ///< the minimum Hit sharing fraction
        bool m_foldBackHierarchy; ///< whether to fold the hierarchy back to the primary (neutrino) or leading particles (test beam)
    };

    /**
     *  @brief  Returns true if passed particle whose primary meets the passed criteria
     *
     *  @param  pMCParticle the input mc particle
     *  @param  fCriteria the given criteria
     */
    static bool DoesPrimaryMeetCriteria(const pandora::MCParticle *const pMCParticle, std::function<bool(const pandora::MCParticle *const)> fCriteria);

    /**
     *  @brief  Returns true if passed particle whose leading meets the passed criteria
     *
     *  @param  pMCParticle the input mc particle
     *  @param  fCriteria the given criteria
     */
    static bool DoesLeadingMeetCriteria(const pandora::MCParticle *const pMCParticle, std::function<bool(const pandora::MCParticle *const)> fCriteria);

    /**
     *  @brief  Returns true if passed a primary neutrino final state MCParticle
     */
    static bool IsBeamNeutrinoFinalState(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a primary triggered beam MCParticle
     */
    static bool IsTriggeredBeamParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a primary beam MCParticle
     */
    static bool IsBeamParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a leading beam MCParticle
     */
    static bool IsLeadingBeamParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Return true if passed a primary cosmic ray MCParticle
     */
    static bool IsCosmicRay(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Return true if passed a photon or electorn/position
     */
    static bool IsEM(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the nuance code of an MCParticle
     */
    static unsigned int GetNuanceCode(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is a neutrino or antineutrino
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return boolean
     */
    static bool IsNeutrino(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is part of a CC interaction
     *
     *  @param[in]  pMCParticle the input mc particle
     *
     *  @return boolean
     */
    static bool IsCCInteraction(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a provided mc particle matches the implemented definition of being primary
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return boolean
     */
    static bool IsPrimary(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a provided mc particle matches the implemented definition of being leading
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return boolean
     */
    static bool IsLeading(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Determine the position in the hierarchy for the MCParticle
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return integer
     */
    static int GetHierarchyTier(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is visible (i.e. long-lived charged particle)
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return boolean
     */
    static bool IsVisible(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get neutrino MC particles from an input MC particle list
     *
     *  @param  pMCParticleList the input MC particle list
     *  @param  trueNeutrinos to receive the vector of neutrino MC particles
     */
    static void GetTrueNeutrinos(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &trueNeutrinos);

    /**
     *  @brief  Get triggered test beam MC particles from an input MC particle list
     *
     *  @param  pMCParticleList the input MC particle list
     *  @param  trueTestBeamParticles to receive the vector of neutrino MC particles
     */
    static void GetTrueTestBeamParticles(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &trueTestBeamParticles);

    /**
     *  @brief  Retrieve the true neutrino vertex.
     *
     *  @param  pMCParticleList The list of MC particles to search for the true vertex
     *  @param  trueVertex The output true vertex position
     *
     *  @returns true if a vertex is found, false otherwise
     */
    static bool GetTrueVertex(const pandora::MCParticleList *const pMCParticleList, pandora::CartesianVector &trueVertex);

    /**
     *  @brief  Get the first visible MC particles given a root particle. For example, given a neutrino this would return the primaries (the visible
     *          final state particles or the first visible descendents of invisible final state particles - note photons and neutrons are considered
     *          visible for this purpose).
     *
     *  @param  pRoot the input mc particle
     *  @param  visibleParticleList the output list of visible particles (if pRoot is visible this will contain only pRoot)
     */
    static void GetFirstVisibleMCParticles(const pandora::MCParticle *const pRoot, pandora::MCParticleList &visibleParticleList);

    /**
     *  @brief  Get the primary parent mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetPrimaryMCParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the leading particle in the hierarchy, for use at ProtoDUNE
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetLeadingMCParticle(const pandora::MCParticle *const pMCParticle, const int hierarchyTierLimit = 1);

    /**
     *  @brief  Get vector of primary MC particles from an input list of MC particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcPrimaryVector the output mc particle vector
     */
    static void GetPrimaryMCParticleList(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &mcPrimaryVector);

    /**
     *  @brief  Get vector of leading MC particles from an input list of MC particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcLeadingVector the output mc particle vector
     */
    static void GetLeadingMCParticleList(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &mcLeadingVector);

    /**
     *  @brief  Get the parent mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the parent mc particle
     */
    static const pandora::MCParticle *GetParentMCParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the next parent mc particle one tier below the current mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the parent mc particle one tier away
     */
    static const pandora::MCParticle *GetNextParentMCParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get all descendent mc particles
     *
     *  @param  pMCParticle the input mc particle
     *  @param  descendentMCParticleList the output descendent mc particle list
     */
    static void GetAllDescendentMCParticles(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &descendentMCParticleList);

    /**
     *  @brief  Get all descendent mc particles, separated into track-like, shower-like and neutron branches.
     *          This method collects together all track-like particles (i.e. not electron, photon or neutron) downstream of the root
     *          particle, stopping at a leading shower or neutron and then storing that leading shower or neutron in a separate list.
     *
     *  @param  pMCParticle the input mc particle
     *  @param  descendentTrackParticles the output list of descendent track-like particles
     *  @param  leadingShowerParticles the output list of leading shower particles
     *  @param  leadingNeutrons the output list of leading neutrons
     */
    static void GetAllDescendentMCParticles(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &descendentTrackParticles,
        pandora::MCParticleList &leadingShowerParticles, pandora::MCParticleList &leadingNeutrons);

    /**
     *  @brief  Get all ancestor mc particles
     *
     *  @param  pMCParticle the input mc particle
     *  @param  ancestorMCParticleList the output ancestor mc particle list
     */
    static void GetAllAncestorMCParticles(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &ancestorMCParticleList);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) and their primary parent mc particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcPrimaryMap the output mapping between mc particles and their parents
     */
    static void GetMCPrimaryMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) and their leading parent mc particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcLeadingMap the output mapping between mc particles and their leading parent
     */
    static void GetMCLeadingMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcLeadingMap);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) to themselves (to be used when not folding particles to their primaries)
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcToSelfMap the output mapping between mc particles and themselves
     */
    static void GetMCToSelfMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcToSelfMap);

    /*
     *  @brief  Retrieve the map from MC to calo hits for all particles
     *
     *  @param  caloHitList The calo hit list for which MC particle matches are to be determined
     *  @param  mcToHitsMap The map to populate
     *  @param  allowSharing Whether to allow hits to be shared between multiple MC particles
     **/
    static void GetMCToHitsMap(const pandora::CaloHitList &caloHitList, MCContributionMap &mcToHitsMap, const bool allowSharing = false);

    /*
     *  @brief  Retrieve the map from MC to calo hits for reconstructable particles
     *
     *  @param  pCaloHitList2D The calo hit list for which MC particle matches are to be determined
     *  @param  pMCParticleList The MC particle list for which calo hit matches are to be determined
     *  @param  mcToHitsMap The map to populate
     **/
    static void GetMCToHitsMap(const pandora::CaloHitList *const pCaloHitList2S, const pandora::MCParticleList *const pMCParticleList,
        LArMCParticleHelper::MCContributionMap &mcToHitsMap);

    /*
     *  @brief  Construct a list of the MC particles from the MC to calo hits map, completing the interaction hierarchy with the invisible
     *          upstream particles.
     *
     *  @param  mcToHitsMap The map of reconstructible MC particles to calo hits
     *  @param  mcHierarchy The output list of MC particles representing the interaction
     *
     *  @return The StatusCode resulting from the function
     **/
    static void CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, pandora::MCParticleList &mcHierarchy);

    /**
     *  @brief  Find the mc particle making the largest contribution to 2D clusters in a specified pfo
     *
     *  @param  pPfo address of the pfo to examine
     *
     *  @return address of the main mc particle
     */
    static const pandora::MCParticle *GetMainMCParticle(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Sort mc particles by their momentum
     *
     *  @param  pLhs address of first mc particle
     *  @param  pRhs address of second mc particle
     */
    static bool SortByMomentum(const pandora::MCParticle *const pLhs, const pandora::MCParticle *const pRhs);

    /**
     *  @brief  Match calo hits to their parent particles
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  mcToTargetMCMap the mc particle to target (primary or self) mc particle map
     *  @param  hitToMCMap output mapping between calo hits and their main MC particle
     *  @param  mcToTrueHitListMap output mapping between MC particles and their associated hits
     */
    static void GetMCParticleToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const MCRelationMap &mcToTargetMCMap,
        CaloHitToMCMap &hitToMCMap, MCContributionMap &mcToTrueHitListMap);

    /**
     *  @brief  Select target, reconstructable mc particles that match given criteria.
     *
     *  @param  pMCParticleList the address of the list of MCParticles
     *  @param  pCaloHitList the address of the list of CaloHits
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected mcparticles to their hits
     */
    static void SelectReconstructableMCParticles(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList,
        const PrimaryParameters &parameters, std::function<bool(const pandora::MCParticle *const)> fCriteria,
        MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Select target, reconstructable mc particles in the relevant hierarchy that match given criteria.
     *
     *  @param  pMCParticleList the address of the list of MCParticles
     *  @param  pCaloHitList the address of the list of CaloHits
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected mcparticles to their hits
     */
    static void SelectReconstructableTestBeamHierarchyMCParticles(const pandora::MCParticleList *pMCParticleList,
        const pandora::CaloHitList *pCaloHitList, const PrimaryParameters &parameters,
        std::function<bool(const pandora::MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Get mapping from Pfo to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMap the input mapping from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     *  @param  foldBackHierarchy whether to fold the particle hierarchy back to the primaries
     */
    static void GetPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap,
        PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy);

    /**
     *  @brief  Get mapping from Pfo in reconstructed test beam hierarchy to reconstructable 2D hits (=good hits belonging to a selected
     *          reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMap the input mapping from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     *  @param  foldBackHierarchy whether to fold the particle hierarchy back to the leading particles
     */
    static void GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList,
        const MCContributionMap &selectedMCParticleToHitsMap, PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy);

    /**
     *  @brief  Get mapping from Pfo to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMaps the input vector of mappings from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     *  @param  foldToHierarchy whether to fold the particle hierarchy back to the primaries
     */
    static void GetPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy);

    /**
     *  @brief  Get mapping from Pfo in reconstructed test beam hierarchy to reconstructable 2D hits (=good hits belonging to a selected
     *          reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMaps the input vector of mappings from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     *  @param  foldBackHierarchy whether to fold the particle hierarchy back to the leading particles
     */
    static void GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList,
        const MCContributionMapVector &selectedMCParticleToHitsMaps, PfoContributionMap &pfoToReconstructable2DHitsMap, const bool foldBackHierarchy);

    /**
     *  @brief  Get mapping from cluster to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  clusterList the input list of clusters
     *  @param  selectedMCToHitsMap the input mapping from selected reconstructable MCParticles to their hits
     *  @param  clusterToReconstructable2DHitsMap the output mapping from clusters to their reconstructable 2D hits
     */
    static void GetClusterToReconstructable2DHitsMap(const pandora::ClusterList &clusterList, const MCContributionMap &selectedMCToHitsMap,
        ClusterContributionMap &clusterToReconstructable2DHitsMap);

    /**
     *  @brief  Get mapping from cluster to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  clusterList the input list of clusters
     *  @param  selectedMCToHitsMaps the input vector of mappings from selected reconstructable MCParticles to their hits
     *  @param  clusterToReconstructable2DHitsMap the output mapping from cluster to their reconstructable 2D hits
     */
    static void GetClusterToReconstructable2DHitsMap(const pandora::ClusterList &clusterList,
        const MCContributionMapVector &selectedMCToHitsMaps, ClusterContributionMap &clusterToReconstructable2DHitsMap);

    /**
     *  @brief  Get the mappings from Pfo -> pair (reconstructable MCparticles, number of reconstructable 2D hits shared with Pfo)
     *          reconstructable MCParticle -> pair (Pfo, number of reconstructable 2D hits shared with MCParticle)
     *
     *  @param  pfoToReconstructable2DHitsMap the input mapping from Pfos to reconstructable 2D hits
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  pfoToMCParticleHitSharingMap the output mapping from Pfos to selected reconstructable MCParticles and the number hits shared
     *  @param  mcParticleToPfoHitSharingMap the output mapping from selected reconstructable MCParticles to Pfos and the number hits shared
     */
    static void GetPfoMCParticleHitSharingMaps(const PfoContributionMap &pfoToReconstructable2DHitsMap,
        const MCContributionMapVector &selectedMCParticleToHitsMaps, PfoToMCParticleHitSharingMap &pfoToMCParticleHitSharingMap,
        MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap);

    /**
     *  @brief  Select a subset of calo hits representing those that represent "reconstructable" regions of the event
     *
     *  @param  pCaloHitList the address of the input calo hit list
     *  @param  mcToTargetMCMap the mc particle to target (primary or self) mc particle map
     *  @param  selectedCaloHitList to receive the populated selected calo hit list
     *  @param  selectInputHits whether to select input hits
     *  @param  maxPhotonPropagation the maximum photon propagation length
     */
    static void SelectCaloHits(const pandora::CaloHitList *const pCaloHitList, const MCRelationMap &mcToTargetMCMap,
        pandora::CaloHitList &selectedCaloHitList, const bool selectInputHits, const float maxPhotonPropagation);

    /**
     *  @brief  Determine if the MC particle is a descendent of a particle with the given PDG code.
     *
     *  @param  pMCParticle the descendent particle
     *  @param  pdg the PDG code of the ancestor particle
     *  @param  isChargeSensitive whether or not to consider the sign of the PDG code when looking for the ancestor (default: false)
     *
     *  @return true if the MC particle has an ancestor with the matching PDG code, false otherwise
     */
    static bool IsDescendentOf(const pandora::MCParticle *const pMCParticle, const int pdg, const bool isChargeSensitive = false);

    /**
     *  @brief  Retrieve a linearised representation of the MC particle hierarchy in breadth first order. This iterates over the MC
     *          hierarchy in a manor that sees primaries at the front of the list, with progressively deeper tiers later in the list.
     *          This is useful for some visualisation cases.
     *
     *  @param  pMCParticle an MC particle in the hierarchy - can be any particle
     *  @param  mcParticleList the output MC particle list
     */
    static void GetBreadthFirstHierarchyRepresentation(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &mcParticleList);

    /**
     *  @brief  Filter an input vector of MCParticles to ensure they have sufficient good hits to be reconstructable
     *
     *  @param  candidateTargets candidate reconstructable MCParticles
     *  @param  mcToTrueHitListMap mapping from candidates reconstructable MCParticles to their true hits
     *  @param  mcToTargetMCMap the mc particle to target (primary or self) mc particle map
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected mcparticles to their hits
     */
    static void SelectParticlesByHitCount(const pandora::MCParticleVector &candidateTargets, const MCContributionMap &mcToTrueHitListMap,
        const MCRelationMap &mcToTargetMCMap, const PrimaryParameters &parameters, MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Get the hits in the intersection of two hit lists
     *
     *  @param  hitListA an input hit list
     *  @param  hitListB another input hit list
     *
     *  @return The hits that are found in both hitListA and hitListB
     */
    static pandora::CaloHitList GetSharedHits(const pandora::CaloHitList &hitListA, const pandora::CaloHitList &hitListB);

    /*
     *  @brief  Check whether or not an MC particle comes from a Bremsstrahlung process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from Bremsstrahlung
     */
    static bool IsBremsstrahlung(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle comes from a capture process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from capture
     */
    static bool IsCapture(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle comes from a decay process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from decay
     */
    static bool IsDecay(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle came from an elastic scattering process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from an elastic scatter
     */
    static bool IsElasticScatter(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle came from an inelastic scattering process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from an inelastic scatter
     */
    static bool IsInelasticScatter(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle comes from an ionisation process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from ionisation
     */
    static bool IsIonisation(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle comes from a nuclear interaction process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from nuclear interaction
     */
    static bool IsNuclear(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether or not an MC particle comes from a pair production process
     *
     *  @param  pMCParticle The MC particle to consider
     *
     *  @return  Whether or not the MC particle came from pair production
     */
    static bool IsPairProduction(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Check whether an MC particle's decay tree contains more than one electron, indicating an EM shower.
     *
     *  @param  pMC                   The MC particle that defines the root of the subtree to check
     *  @param  nDescendantElectrons  The running count of electrons seen so far,
     *                                this is a recursive parameter that should be 0 at the top level
     *
     *  @return  Whether the MC particle's subtree contains more than one electron
     */
    static bool CausesShower(const pandora::MCParticle *const pMC, int nDescendantElectrons = 0);

    /**
     *  @brief  Determine if two MC particles are topologically continuous within a given tolerance.
     *          If the parent does not travel any distance, a travelling parent is sought and the comparison made between this and the child.
     *          If no travelling parent can be found, the particles are treated as continuous.
     *
     *  @param  pMCParent The parent MC particle
     *  @param  pMCChild The child MC particle
     *  @param  cosAngleTolerance The cosine of the maximum acceptable angular deviation
     *
     *  @return True if the partiles are topologically continous, false otherwise.
     */
    static bool AreTopologicallyContinuous(
        const pandora::MCParticle *const pMCParent, const pandora::MCParticle *const pMCChild, const float cosAngleTolerance);

private:
    /**
     *  @brief  For a given Pfo, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pPfo the input pfo
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     *  @param  foldBackHierarchy whether to fold the particle hierarchy back to primaries
     */
    static void CollectReconstructable2DHits(const pandora::ParticleFlowObject *const pPfo, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        pandora::CaloHitList &reconstructableCaloHitList2D, const bool foldBackHierarchy);

    /**
     *  @brief  For a given Pfo, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *          and belong in the test beam particle interaction hierarchy
     *
     *  @param  pPfo the input pfo
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     *  @param  foldBackHierarchy whether to fold the particle hierarchy back to leading particles
     */
    static void CollectReconstructableTestBeamHierarchy2DHits(const pandora::ParticleFlowObject *const pPfo,
        const MCContributionMapVector &selectedMCParticleToHitsMaps, pandora::CaloHitList &reconstructableCaloHitList2D, const bool foldBackHierarchy);

    /**
     *  @brief  For a given Pfo list, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input pfo list
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     */
    static void CollectReconstructable2DHits(const pandora::PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        pandora::CaloHitList &reconstructableCaloHitList2D);

    /**
     *  @brief  For a given cluster, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pCluster the input cluster
     *  @param  selectedMCToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     */
    static void CollectReconstructable2DHits(const pandora::Cluster *const pCluster,
        const MCContributionMapVector &selectedMCParticleToHitsMaps, pandora::CaloHitList &reconstructableCaloHitList2D);

    /**
     *  @brief  Apply further selection criteria to end up with a collection of "good" calo hits that can be use to define whether
     *          a target mc particle is reconstructable.
     *
     *  @param  pSelectedCaloHitList the address of the calo hit list (typically already been through some selection procedure)
     *  @param  mcToTargetMCMap the mc particle to target (primary or self) mc particle map
     *  @param  selectedGoodCaloHitList to receive the populated good selected calo hit list
     *  @param  selectInputHits whether to select input hits
     *  @param  minHitSharingFraction the minimum Hit sharing fraction
     */
    static void SelectGoodCaloHits(const pandora::CaloHitList *const pSelectedCaloHitList, const MCRelationMap &mcToTargetMCMap,
        pandora::CaloHitList &selectedGoodCaloHitList, const bool selectInputHits, const float minHitSharingFraction);

    /**
     *  @brief  Select mc particles matching given criteria from an input list
     *
     *  @param  inputMCParticles input vector of MCParticles
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     *  @param  selectedParticles the output vector of particles selected
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  isTestBeam whether the mc particles correspond to the test beam case or the neutrino case
     */
    static void SelectParticlesMatchingCriteria(const pandora::MCParticleVector &inputMCParticles,
        std::function<bool(const pandora::MCParticle *const)> fCriteria, pandora::MCParticleVector &selectedParticles,
        const PrimaryParameters &parameters, const bool isTestBeam);

    /**
     *  @brief  Whether it is possible to navigate from a primary mc particle to a downstream mc particle without "passing through" a neutron
     *
     *  @param  pOriginalPrimary the address of the original primary mc particle
     *  @param  pThisMCParticle the address of the current mc particle in the primary decay chain
     *  @param  pHitMCParticle the address of the mc particle associated to a calo hit
     *  @param  maxPhotonPropagation the maximum photon propagation length
     *
     *  @return boolean
     */
    static bool PassMCParticleChecks(const pandora::MCParticle *const pOriginalPrimary, const pandora::MCParticle *const pThisMCParticle,
        const pandora::MCParticle *const pHitMCParticle, const float maxPhotonPropagation);
};

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_HELPER_H
