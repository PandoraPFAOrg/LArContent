/**
 *  @file   larpandoracontent/LArHelpers/LArHierarchyHelper.h
 *
 *  @brief  Header file for the lar hierarchy helper class.
 *
 *  $Log: $
 */
#ifndef LAR_HIERARCHY_HELPER_H
#define LAR_HIERARCHY_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

namespace lar_content
{

/**
 *  @brief  LArHierarchyHelper class
 */
class LArHierarchyHelper
{
public:
    /**
     *  @brief  FoldingParameters class
     */
    class FoldingParameters
    {
    public:
        /**
         *  @brief  Default constructor
         */
        FoldingParameters();

        /**
         *  @brief  Constructor
         *
         *  @param  foldDynamic Whether or not to apply dynamic folding to the hierarchy
         */
        FoldingParameters(const bool foldDynamic, const float cosAngleTolerance = 0.9962f);

        /**
         *  @brief  Constructor.
         *
         *  If folding back to tier 2, any MC particle/PFO ("particles") at tier 1 will be allocated their own node. At tier 2, the
         *  particles will be allocated as the main particle for a node and all of their children will also be incorprated into the node
         *
         *  @param  foldingTier The tier at which level child particles should be folded back (> 0)
         */
        FoldingParameters(const int foldingTier);

        bool m_foldToLeadingShowers; ///< Whether or not to fold shower children to the leading shower particle
        bool m_foldToTier;           ///< Whether or not to apply folding based on particle tier
        bool m_foldDynamic;          ///< Whether or not to use process and topological information to make folding decisions
        float m_cosAngleTolerance;   ///< Cosine of the maximum angle at which topologies can be considered continuous
        int m_tier;                  ///< If folding to a tier, the tier to be combined with its child particles
    };

    /**
     *  @brief  QualityCuts class
     */
    class QualityCuts
    {
    public:
        /**
         *  @brief Default constructor
         */
        QualityCuts();

        /**
         *  @brief Constructor
         *
         *  @param  minPurity The minimum purity for a cut to be considered good
         *  @param  minCompleteness The minimum completeness for a cut to be considered good
	 *  @param  selectRecoHits Whether to only use reco hits that overlap with the MC particle hits
         */
        QualityCuts(const float minPurity, const float minCompleteness, const bool selectRecoHits = false);

        const float m_minPurity;       ///< The minimum purity for a match to be considered good
        const float m_minCompleteness; ///< The minimum completeness for a match to be considered good
        const bool m_selectRecoHits;   ///< Whether to only use reco hits that overlap with the MC particle hits
    };

    /**
     *  @brief  MCHierarchy class
     */
    class MCHierarchy
    {
    public:
        /**
         *  @brief   ReconstructabilityCriteria class
         */
        class ReconstructabilityCriteria
        {
        public:
            /**
             *  @brief  Default constructor
             */
            ReconstructabilityCriteria();

            /**
             *  @brief  Copy constructor
             */
            ReconstructabilityCriteria(const ReconstructabilityCriteria &obj);

            /**
             *  @brief  Constructor
             *
             *  @param  minHits The total minimum number of hits for a particle to be considered reconstructable
             *  @param  minHitsForGoodView The number of hits within a view for a particle to be considered reconstructable
             *  @param  minGoodViews The minimum number of good views for a particle to be considered reconstructable
             *  @param removeNeutrons Whether to remove neutrons and downstream particles from consideration
             */
            ReconstructabilityCriteria(const unsigned int minHits, const unsigned int minHitsForGoodView, const unsigned int minGoodViews,
                const bool removeNeutrons);

            const unsigned int m_minHits;            ///< the minimum number of primary good Hits
            const unsigned int m_minHitsForGoodView; ///< the minimum number of Hits for a good view
            const unsigned int m_minGoodViews;       ///< the minimum number of primary good views
            const bool m_removeNeutrons;             ///< whether to remove neutrons and their downstream particles
        };

        class Node;
        typedef std::vector<const Node *> NodeVector;
        typedef std::list<const Node *> NodeList;
        typedef std::map<const pandora::MCParticle *, NodeVector> MCNodeVectorMap;

        /**
         *  @brief  Node class
         */
        class Node
        {
        public:
            /**
             *  @brief  Create a node with a primary MC particle
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  pMCParticle The primary MC particle with which this node should be created
             *  @param  tier The tier that should be assigned to this node
             */
            Node(MCHierarchy &hierarchy, const pandora::MCParticle *pMCParticle, const int tier = 1);

            /**
             *  @brief  Create a node from a list of MC particles
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  mcParticleList The MC particle list with which this node should be created
             *  @param  caloHitList The CaloHit list with which this node should be created
             *  @param  tier The tier that should be assigned to this node
             */
            Node(MCHierarchy &hierarchy, const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList, const int tier = 1);

            /**
             *  @brief Destructor
             */
            virtual ~Node();

            /**
             *  @brief  Return whether or not this node should be considered reconstructable
             *
             *  @return true if reconstructable, false otherwise
             */
            bool IsReconstructable() const;

            /**
             *  @brief  Recursively fill the hierarchy based on the criteria established for this MCHierarchy
             *
             *  @param  pRoot The MC particle acting as the root of the current branch of the hierarchy
             *  @param  foldParameters The folding parameters
             */
            void FillHierarchy(const pandora::MCParticle *pRoot, const FoldingParameters &foldParameters);

            /**
             *  @brief  Fill this node by folding all descendent particles to this node
             *
             *  @param  pRoot The MC particle acting as the root of the current branch of the hierarchy
             */
            void FillFlat(const pandora::MCParticle *pRoot);

            /**
             *  @brief  Return the vector of children for this node
             *
             *  @return The vector of children
             */
            const NodeVector &GetChildren() const;

            /**
             *  @brief  Retrieve the unique ID of this node
             *
             *  @return The unique ID of this node
             */
            int GetId() const;

            /**
             *  @brief  Retrieve the leading MC particle associated with this node
             *
             *  @return The main MC particle associated with this node
             */
            const pandora::MCParticle *GetLeadingMCParticle() const;

            /**
             *  @brief  Retrieve the MC particles associated with this node
             *
             *  @return The MC particles associated with this node
             */
            const pandora::MCParticleList &GetMCParticles() const;

            /**
             *  @brief  Retrieve the CaloHits associated with this node
             *
             *  @return The list of CaloHits associated with this node
             */
            const pandora::CaloHitList &GetCaloHits() const;

            /**
             *  @brief  Retrieve the PDG code for the leading particle in this node
             *
             *  @return The PDG code for the leading particle in this node
             */
            int GetParticleId() const;

            /**
             *  @brief  Retrieve the hierarchy tier of this node
             *
             *  @return The hierarchy tier of this node
             */
            int GetHierarchyTier() const;

            /**
             *  @brief  Check if this is a particle induced by a neutrino interaction
             *
             *  @return Whether or not this is neutrino induced
             */
            bool IsNeutrinoInduced() const;

            /**
             *  @brief  Check if this is a test beam particle
             *
             *  @return Whether or not this is a test beam particle
             */
            bool IsTestBeamParticle() const;

            /**
             *  @brief  Check if this is a cosmic ray particle
             *
             *  @return Whether or not this is a cosmic ray
             */
            bool IsCosmicRay() const;

            /**
             *  @brief  Returns whether or not this particle is the leading lepton in the event
             *
             *  @return Whether or not this is the leading lepton
             */
            bool IsLeadingLepton() const;

            /**
             *  @brief  Produce a string representation of the hierarchy
             *
             *  @return The string representation of the hierarchy
             */
            const std::string ToString(const std::string &prefix) const;

        private:
            /**
             *  @brief  Tags the particle as the leading lepton
             */
            void SetLeadingLepton();

            MCHierarchy &m_hierarchy;                  ///< The parent MC hierarchy
            pandora::MCParticleList m_mcParticles;     ///< The list of MC particles of which this node is composed
            pandora::CaloHitList m_caloHits;           ///< The list of calo hits of which this node is composed
            NodeVector m_children;                     ///< The child nodes of this node
            const pandora::MCParticle *m_mainParticle; ///< The leading MC particle for this node
            int m_tier;                                ///< The hierarchy tier for this node
            int m_pdg;                                 ///< The PDG code of the leading MC particle for this node
            bool m_isLeadingLepton;                    ///< Whether or not this node is the leading lepton

            friend class MCHierarchy;
        };

        /**
         *  @brief  Default constructor
         */
        MCHierarchy();

        /**
         *  @brief  Construct a new MCHierarchy object using specified reconstructability criteria
         *
         *  @param  recoCriteria The reconstructability criteria to be applied
         */
        MCHierarchy(const ReconstructabilityCriteria &recoCriteria);

        /**
         *  @brief Destructor
         */
        virtual ~MCHierarchy();

        /**
         *  @brief  Creates an MC hierarchy representation. Without folding this will be a mirror image of the standard MCParticle
         *          relationships. However, with folding options selected the hierarchy structure will group together MC particles into
         *          nodes based on the folding requirements.
         *
         *          If only folding back to primaries, the hierarchy will be relatively flat, with a top-level neutrino or test beam
         *          particle, if appropriate, and then a set of leaf nodes, one for each primary particles also containing the MC particles
         *          (and corresponding hits) from daughter particles.
         *
         *          If only folding back to leading shower particles, the hierarchy will largely mirror the standard MCParticle hierarchy,
         *          but, when a shower particle is reached (for this purpose an electron or photon), this particle and all daughter
         *          particles will be represented by a single leaf node.
         *
         *          If folding back to both primary and leading shower particles the hierarchy will again be rather flat, but in this case,
         *          if a primary track-like particle (i.e. not an electron or photon) has a downstream shower particle then all downstream
         *          particles above the shower-like particle will be folded into the primary node, but a new, daughter leaf node will be
         *          created for the shower-like particle and all of its daughters, and a parent-child relationship will be formed between
         *          the primary node and shower node.
         *
         *  @param  mcParticleList The list of MC particles with which to fill the hierarchy
         *  @param  caloHitList The list of hits with which to fill the hierarchy
         *  @param  foldParameters The folding parameters to use for the hierarchy
         */
        void FillHierarchy(const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList, const FoldingParameters &foldParameters);

        /**
         *  @brief  Interpret the hierarchy below a particular particle to determine if and how it should be folded. Folded particles are
         *          added to the leadingParticles list and child particles are added to the childParticles list.
         *
         *  @param  pRoot The root of the hierarchy to interpret
         *  @param  leadingParticles The output list of particles that should be folded into the root particle
         *  @param  childParticles The output list of particles that should be considered children of the folded particle
         *  @param  cosAngleTolerance The cosine of the maximum angle for which trajectories are considered continuous
         */
        void InterpretHierarchy(const pandora::MCParticle *const pRoot, pandora::MCParticleList &leadingParticles,
            pandora::MCParticleList &childParticles, const float cosAngleTolerance) const;

        /**
         *  @brief  Retrieve the neutrino at the root of the hierarchy if it exists
         *
         *  @return The address of the incident neutrino (nullptr if it doesn't exist)
         */
        const pandora::MCParticle *GetNeutrino() const;

        /**
         *  @brief  Retrieve the root nodes in this hierarchy
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The primary nodes in the requested interaction hierarchy
         */
        const NodeVector &GetInteractions(const pandora::MCParticle *pRoot) const;

        /**
         *  @brief  Retrieve the root MC particles of the interaction hierarchies
         *
         *  @param  rootMCParticles The output list of root MC particles
         */
        void GetRootMCParticles(pandora::MCParticleList &rootMCParticles) const;

        /**
         *  @brief  Retrieve a flat vector of the ndoes in the hierarchy
         *
         *  @param  pRoot The root MC particle for an interaction
         *  @param  nodeVector The output vector for the nodes in the hierarchy in breadth first order
         */
        void GetFlattenedNodes(const pandora::MCParticle *const pRoot, NodeVector &nodeVector) const;

        /**
         *  @brief  Register a node with the hierarchy
         *
         *  @param  pNode  The node to register
         */
        void RegisterNode(const Node *pNode);

        /**
         *  @brief  Produce a string representation of the hierarchy
         *
         *  @return The string representation of the hierarchy
         */
        const std::string ToString() const;

    private:
        /**
         *  @brief  Identify downstream particles that represent continuations of the parent particle from a reconstruction perspective
         *
         *  @param  pRoot The root MC particle
         *  @param  continuingParticles An output list of the particles identified as continuations
         *  @param  childParticles An output list of the particles identified as child particles given any continuations
         *  @param  cosAngleTolerance The cosine of the maximum angle for which trajectories are considered continuous
         */
        void CollectContinuations(const pandora::MCParticle *pRoot, pandora::MCParticleList &continuingParticles,
            pandora::MCParticleList &childParticles, const float cosAngleTolerance) const;

        /**
         *  @brief  Checks if an individual particle meets reconstructability criteria
         *
         *  @param  pMCParticle  The MC particle to assess
         *
         *  @return Whether or not the MC particle meets reconstructability criteria
         */
        bool IsReconstructable(const pandora::MCParticle *pMCParticle) const;

        /**
         *  @brief  Checks if a set of hits meet reconstructability criteria
         *
         *  @param  caloHits  The calo hits to assess
         *
         *  @return Whether or not the hits meet reconstructability criteria
         */
        bool IsReconstructable(const pandora::CaloHitList &caloHits) const;

        MCNodeVectorMap m_interactions;            ///< Map from incident particles (e.g. neutrino) to primaries
        ReconstructabilityCriteria m_recoCriteria; ///< The criteria used to determine if the node is reconstructable
        std::map<const pandora::MCParticle *, pandora::CaloHitList> m_mcToHitsMap; ///< The map between MC particles and calo hits
        std::map<const Node *, int> m_nodeToIdMap;                                 ///< A map from nodes to unique ids
        int m_nextNodeId;                                                          ///< The ID to use for the next node
    };

    /**
     *  @brief   RecoHierarchy class
     */
    class RecoHierarchy
    {
    public:
        class Node;
        typedef std::vector<const Node *> NodeVector;
        typedef std::list<const Node *> NodeList;
        typedef std::map<const pandora::ParticleFlowObject *, NodeVector> RecoNodeVectorMap;

        /**
         *  @brief  Node class
         */
        class Node
        {
        public:
            /**
             *  @brief  Create a node with a primary PFO
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  pPfo The primary PFO with which this node should be created
             *  @param  tier The tier that should be assigned to this node
             */
            Node(const RecoHierarchy &hierarchy, const pandora::ParticleFlowObject *pPfo, const int tier = 1);

            /**
             *  @brief  Create a node from a list of PFOs
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  pfoList The PFO list with which this node should be created
             *  @parasm caloHitList The CaloHit list with which this node should be created
             *  @param  tier The tier that should be assigned to this node
             */
            Node(const RecoHierarchy &hierarchy, const pandora::PfoList &pfoList, const pandora::CaloHitList &caloHitList, const int tier = 1);

            /**
             *  @brief Destructor
             */
            virtual ~Node();

            /**
             *  @brief  Recursively fill the hierarchy based on the criteria established for this RecoHierarchy
             *
             *  @param  pRoot The PFO acting as the root of the current branch of the hierarchy
             *  @param  foldParameters The folding parameters
             */
            void FillHierarchy(const pandora::ParticleFlowObject *pRoot, const FoldingParameters &foldParameters);

            /**
             *  @brief  Fill this node by folding all descendent particles to this node
             *
             *  @param  pRoot The PFO acting as the root of the current branch of the hierarchy
             */
            void FillFlat(const pandora::ParticleFlowObject *pRoot);

            /**
             *  @brief  Return the vector of children for this node
             *
             *  @return The vector of children
             */
            const NodeVector &GetChildren() const;

            /**
             *  @brief  Retrieve the PFOs associated with this node
             *
             *  @return The PFOs associated with this node
             */
            const pandora::PfoList &GetRecoParticles() const;

            /**
             *  @brief  Retrieve the leading reco particle for this node
             *
             *  return  The leading reco particle for this node
             */
            const pandora::ParticleFlowObject *GetLeadingPfo() const;

            /**
             *  @brief  Retrieve the CaloHits associated with this node
             *
             *  @return The list of CaloHits associated with this node
             */
            const pandora::CaloHitList &GetCaloHits() const;

            /**
             *  @brief  Retrieve the PDG code for the leading particle in this node
             *          Note, for reco objects the PDG codes represent tracks (muon PDG) and showers (electron PDG)
             *
             *  @return The PDG code for the leading particle in this node
             */
            int GetParticleId() const;

            /**
             *  @brief  Retrieve the hierarchy tier of this node
             *
             *  @return The hierarchy tier of this node
             */
            int GetHierarchyTier() const;

            /**
             *  @brief  Produce a string representation of the hierarchy
             *
             *  @return The string representation of the hierarchy
             */
            const std::string ToString(const std::string &prefix) const;

        private:
            const RecoHierarchy &m_hierarchy;             ///< The parent reco hierarchy
            pandora::PfoList m_pfos;                      ///< The list of PFOs of which this node is composed
            pandora::CaloHitList m_caloHits;              ///< The list of calo hits of which this node is composed
            NodeVector m_children;                        ///< The child nodes of this node
            const pandora::ParticleFlowObject *m_mainPfo; ///< The leading particle flow object for this node
            int m_tier;                                   ///< The hierarchy tier for this node
            int m_pdg;                                    ///< The particle ID (track = muon, shower = electron)
        };

        /**
         *  @brief  Default constructor
         */
        RecoHierarchy();

        /**
         *  @brief Destructor
         */
        virtual ~RecoHierarchy();

        /**
         *  @brief  Creates a reconstructed hierarchy representation. Without folding this will be a mirror image of the standard
         *          ParticleFlowObject (PFO) relationships. However, with folding options selected the hierarchy structure will group
         *          together PFOs into nodes based on the folding requirements.
         *
         *          If only folding back to primaries, the hierarchy will be relatively flat, with a top-level neutrino or test beam
         *          particle, if appropriate, and then a set of leaf nodes, one for each primary particles also containing the PFOs (and
         *          corresponding hits) from daughter particles.
         *
         *          If only folding back to leading shower particles, the hierarchy will largely mirror the standard PFO hierarchy, but,
         *          when a shower particle is reached (based on the track/shower characterisation), this particle and all daughter particles
         *          will be represented by a single leaf node.
         *
         *          If folding back to both primary and leading shower particles the hierarchy will again be rather flat, but in this case,
         *          if a primary track-like particle has a downstream shower particle then all downstream particles above the shower-like
         *          particle will be folded into the primary node, but a new, daughter leaf node will be created for the shower-like
         *          particle and all of its daughters, and a parent-child relationship will be formed between the primary node and shower
         *          node.
         *
         *  @param  pfoList The list of PFOs with which to fill the hierarchy
         *  @param  foldParameters The folding parameters to use for the hierarchy
         */
        void FillHierarchy(const pandora::PfoList &pfoList, const FoldingParameters &foldParameters);

        /**
         *  @brief  Retrieve the root nodes in the hierarchy for a given interaction
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The root nodes in this hierarchy
         */
        const NodeVector &GetInteractions(const pandora::ParticleFlowObject *pRoot) const;

        /**
         *  @brief  Retrieve the root particle flow objects of the interaction hierarchies
         *
         *  @param  rootPfos The output list of root particle flow objects
         */
        void GetRootPfos(pandora::PfoList &rootPfos) const;

        /**
         *  @brief  Retrieve a flat vector of the nodes in the hierarchy
         *
         *  @param  pRoot The root particle flow object for an interaction
         *  @param  nodeVector The output vector for the nodes in the hierarchy in breadth first order
         */
        void GetFlattenedNodes(const pandora::ParticleFlowObject *const pRoot, NodeVector &nodeVector) const;

        /**
         *  @brief  Retrieve the neutrino at the root of the hierarchy if it exists
         *
         *  @return The address of the incident neutrino (nullptr if it doesn't exist)
         */
        const pandora::ParticleFlowObject *GetNeutrino() const;

        /**
         *  @brief  Produce a string representation of the hierarchy
         *
         *  @return The string representation of the hierarchy
         */
        const std::string ToString() const;

    private:
        RecoNodeVectorMap m_interactions; ///< Map from the root PFO (e.g. neutrino) to primaries
    };

    typedef std::map<const RecoHierarchy::Node *, const pandora::CaloHitList> SelectedRecoHitsMap;

    /**
     *  @brief  MCMatches class
     */
    class MCMatches
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pMCParticle The MCParticle being matched
         */
        MCMatches(const MCHierarchy::Node *pMCParticle);

        /**
         *  @brief  Add a reconstructed node as a match for this MC node
         *
         *  @param  pReco The reconstructed node that matches this MC node
         *  @param  nSharedHits The number of hits shared between reco and MC nodes
	 *  @param  selectedRecoHits The reco node hits that are selected from the MC particle hits
         */
        void AddRecoMatch(const RecoHierarchy::Node *pReco, const int nSharedHits, const pandora::CaloHitList &selectedRecoHits);

        /**
         *  @brief  Retrieve the MC node
         *
         *  @return The MC node
         */
        const MCHierarchy::Node *GetMC() const;

        /**
         *  @brief  Retrieve the vector of matched reco nodes
         *
         *  @return The vector of matched reco nodes
         */
        const RecoHierarchy::NodeVector &GetRecoMatches() const;

        /**
	 *  @brief  Retrieve the selected hits for the given reco node
	 *
	 *  @param  pReco The reco node pointer
	 *
	 *  @return The list of selected calo reco hits
	 */
        const pandora::CaloHitList GetSelectedRecoHits(const RecoHierarchy::Node *pReco) const;

        /**
         *  @brief  Retrieve the number of shared hits in the match
         *
         *  @param  pReco The reco node to consider
         *
         *  @return The number of shared hits
         */
        unsigned int GetSharedHits(const RecoHierarchy::Node *pReco) const;

        /**
         *  @brief  Retrieve the purity of the match
         *
         *  @param  pReco The reco node to consider
         *  @param  adcWeighted Whether or not to weight purity according to the charge contribution
         *
         *  @return The purity of the match
         */
        float GetPurity(const RecoHierarchy::Node *pReco, const bool adcWeighted = false) const;

        /**
         *  @brief  Retrieve the purity of the match
         *
         *  @param  pReco The reco node to consider
         *  @param  view The view for which purity should be calculated
         *  @param  adcWeighted Whether or not to weight purity according to the charge contribution
         *
         *  @return The purity of the match
         */
        float GetPurity(const RecoHierarchy::Node *pReco, const pandora::HitType view, const bool adcWeighted = false) const;

        /**
         *  @brief  Retrieve the completeness of the match
         *
         *  @param  pReco The reco node to consider
         *  @param  adcWeighted Whether or not to weight completeness according to the charge contribution
         *
         *  @return The completeness of the match
         */
        float GetCompleteness(const RecoHierarchy::Node *pReco, const bool adcWeighted = false) const;

        /**
         *  @brief  Retrieve the completeness of the match
         *
         *  @param  pReco The reco node to consider
         *  @param  view The view for which purity should be calculated
         *  @param  adcWeighted Whether or not to weight completeness according to the charge contribution
         *
         *  @return The completeness of the match
         */
        float GetCompleteness(const RecoHierarchy::Node *pReco, const pandora::HitType view, const bool adcWeighted = false) const;

        /**
         *  @brief  Get the number of reco nodes matched (both above and below quality cut thresholds) to the MC node
         *
         *  @return The number of reco nodes matched to the MC node
         */
        size_t GetNRecoMatches() const;

        /**
         *  @brief  Get whether this match passes quality cuts
         *
         *  @param  qualityCuts The quality cuts to pass
         *
         *  @return Whether or not this match passes quality cuts
         */
        bool IsQuality(const QualityCuts &qualityCuts) const;

    private:
        /**
         *  @brief  Core purity calculation given intersecting hits and reco hits
         *
         *  @param  intersection The intersecting reco and MC hits
         *  @param  recoHits The reco hits
         *  @param  adcWeighted Whether or not to weight purity according to the charge contribution
         *
         *  @return The purity of the match
         */
        float GetPurity(const pandora::CaloHitVector &intersection, const pandora::CaloHitList &recoHits, const bool adcWeighted) const;

        /**
         *  @brief  Core completeness calculation given intersecting hits and MC hits
         *
         *  @param  intersection The intersecting reco and MC hits
         *  @param  mcHits The MC hits
         *  @param  adcWeighted Whether or not to weight completeness according to the charge contribution
         *
         *  @return The completeness of the match
         */
        float GetCompleteness(const pandora::CaloHitVector &intersection, const pandora::CaloHitList &mcHits, const bool adcWeighted) const;

        const MCHierarchy::Node *m_pMCParticle;    ///< MC node associated with any matches
        RecoHierarchy::NodeVector m_recoNodes;     ///< Matched reco nodes
        pandora::IntVector m_sharedHits;           ///< Number of shared hits for each match
        SelectedRecoHitsMap m_selectedRecoHitsMap; ///< Map storing the selected CaloHits for a given reco node
    };

    typedef std::vector<MCMatches> MCMatchesVector;
    typedef std::map<const pandora::MCParticle *, MCMatchesVector> InteractionInfo;

    /**
     *  @brief  MatchInfo class
     */
    class MatchInfo
    {
    public:
        /**
         *  @brief  Default constructor
         *
         *  @param  mcHierarchy The MC hierarchy
         *  @param  recoHierarchy The reco hierarchy
         */
        MatchInfo(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy);

        /**
         *  @brief  Constructor
         *
         *  @param  mcHierarchy The MC hierarchy
         *  @param  recoHierarchy The reco hierarchy
         *  @param  qualityCuts The quality cuts to be applied to matched nodes
         */
        MatchInfo(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy, const QualityCuts &qualityCuts);

        /**
         *  @brief  Match the nodes in the MC and reco hierarchies.
         *
         */
        void Match();

        /**
         *  @brief  Retrieve the vector of matches (this will include null matches - i.e. MC nodes with no corresponding reco)
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The vector of matches
         */
        const MCMatchesVector &GetMatches(const pandora::MCParticle *const pRoot) const;

        /**
         *  @brief  Retrieve the vector of unmatched reco nodes
         */
        const RecoHierarchy::NodeVector &GetUnmatchedReco() const;

        /**
         *  @brief  Retrieve the number of MC nodes available to match
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The number of MC nodes available to match
         */
        unsigned int GetNMCNodes(const pandora::MCParticle *const pRoot) const;

        /**
         *  @brief  Retrieve the number of neutrino interaction derived MC nodes available to match
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The number of MC nodes available to match
         */
        unsigned int GetNNeutrinoMCNodes(const pandora::MCParticle *const pRoot) const;

        /**
         *  @brief  Retrieve the number of cosmic ray derived MC nodes available to match
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The number of MC nodes available to match
         */
        unsigned int GetNCosmicRayMCNodes(const pandora::MCParticle *const pRoot) const;

        /**
         *  @brief  Retrieve the number of test beam derived MC nodes available to match
         *
         *  @param  pRoot The root of the interaction hierarchy
         *
         *  @return The number of MC nodes available to match
         */
        unsigned int GetNTestBeamMCNodes(const pandora::MCParticle *const pRoot) const;

        /**
         *  @brief  Retrieve the MC hierarchy used for the matching
         *
         *  @return The MCHierarchy used for matching
         */
        const MCHierarchy &GetMCHierarchy() const;

        /**
         *  @brief  Retrieve the reco hierarchy used for the matching
         *
         *  @return The RecoHierarchy used for matching
         */
        const RecoHierarchy &GetRecoHierarchy() const;

        /**
         *  @brief  Retrieve the root MC particles of the interaction hierarchies
         *
         *  @param  rootMCParticles The output list of root MC particles
         */
        void GetRootMCParticles(pandora::MCParticleList &rootMCParticles) const;

        /**
         *  @brief  Retrieve the quality cuts for matching
         *
         *  @return The quality cuts
         */
        const QualityCuts &GetQualityCuts() const;

        /**
	 *  @brief  Get the selected reco node hits that are contained in the allMCHits list
	 *
	 *  @param  pRecoNode The reco node pointer
	 *  @param  allMCHits The CaloHit list of all of the MC particles that will be used for the selection
     *  @param  pMCHitIds Optional map for looking up IDs of MC. Made if missing.
	 *
	 *  @return The selected reco calo hits
	 */
        const pandora::CaloHitList GetSelectedRecoHits(const RecoHierarchy::Node *pRecoNode, const pandora::CaloHitList &allMCHits,
            const std::unordered_set<intptr_t> *pMCHitIds = nullptr) const;

        /**
         *  @brief  Prints information about which reco nodes are matched to the MC nodes, information about hit sharing, purity and
         *          completeness.
         *
         *  @param  mcHierarchy The MC hierarchy
         */
        void Print(const MCHierarchy &mcHierarchy) const;

    private:
        const MCHierarchy &m_mcHierarchy;          ///< The MC hierarchy for the matching procedure
        const RecoHierarchy &m_recoHierarchy;      ///< The Reco hierarchy for the matching procedure
        InteractionInfo m_matches;                 ///< The map between an interaction and the vector of good matches from MC to reco
        RecoHierarchy::NodeVector m_unmatchedReco; ///< The vector of unmatched reco nodes
        QualityCuts m_qualityCuts;                 ///< The quality cuts to be applied to matches
    };

    /**
     *  @brief  Fill an MC hierarchy based on the specified folding criteria (see MCHierarchy::FillHierarchy for details)
     *
     *  @param  mcParticleList The MCParticle list to use to fill this hierarchy
     *  @param  caloHitList The list of CaloHits to use to fill this hierarchy
     *  @param  foldParameters The folding parameters to use for the hierarchy
     *  @param  hierarchy The output MC hierarchy
     */
    static void FillMCHierarchy(const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList,
        const FoldingParameters &foldParameters, MCHierarchy &hierarchy);

    /**
     *  @brief  Fill a reconstructed hierarchy based on the specified folding criteria (see RecoHierarchy::FillHierarchy for details)
     *
     *  @param  pfoList The ParticleFlowObject list to use to fill this hierarchy
     *  @param  foldParameters The folding parameters to use for the hierarchy
     *  @param  hierarchy The output reconstructed hierarchy
     */
    static void FillRecoHierarchy(const pandora::PfoList &pfoList, const FoldingParameters &foldParameters, RecoHierarchy &hierarchy);

    /**
     *  @brief  Finds the matches between reconstructed and MC hierarchies.
     *
     *  @param  matchInfo The output match information
     */
    static void MatchHierarchies(MatchInfo &matchInfo);

private:
    typedef std::set<const pandora::MCParticle *> MCParticleSet;
    typedef std::set<const pandora::ParticleFlowObject *> PfoSet;

    /**
     *  @brief  Retrieves the primary MC particles from a list and returns the root (neutrino) for hierarchy, if it exists.
     *
     *  @param  pRoot The root MC particle (e.g. neutrino) for which primaries should be collected
     *  @param  primaries The output set of primary MC particles
     */
    static void GetMCPrimaries(const pandora::MCParticle *pRoot, MCParticleSet &primaries);

    /**
     *  @brief  Retrieves the primary PFOs from a list and returns the root (neutrino) for hierarchy, if it exists.
     *
     *  @param  pRoot The root particle flow object (e.g. neutrino) for which primaries should be collected
     *  @param  primaries The output set of primary PFOs
     */
    static void GetRecoPrimaries(const pandora::ParticleFlowObject *pRoot, PfoSet &primaries);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy::NodeVector &LArHierarchyHelper::MCHierarchy::Node::GetChildren() const
{
    return m_children;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticleList &LArHierarchyHelper::MCHierarchy::Node::GetMCParticles() const
{
    return m_mcParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArHierarchyHelper::MCHierarchy::Node::GetCaloHits() const
{
    return m_caloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArHierarchyHelper::MCHierarchy::Node::GetLeadingMCParticle() const
{
    return m_mainParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArHierarchyHelper::MCHierarchy::Node::GetParticleId() const
{
    return m_pdg;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArHierarchyHelper::MCHierarchy::Node::GetHierarchyTier() const
{
    return m_tier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArHierarchyHelper::MCHierarchy::Node::IsNeutrinoInduced() const
{
    return !(LArHierarchyHelper::MCHierarchy::Node::IsTestBeamParticle() || LArHierarchyHelper::MCHierarchy::Node::IsCosmicRay());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArHierarchyHelper::MCHierarchy::Node::IsLeadingLepton() const
{
    return m_isLeadingLepton;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArHierarchyHelper::MCHierarchy::Node::SetLeadingLepton()
{
    m_isLeadingLepton = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::RecoHierarchy::Node::GetChildren() const
{
    return m_children;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *LArHierarchyHelper::RecoHierarchy::Node::GetLeadingPfo() const
{
    return m_mainPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArHierarchyHelper::RecoHierarchy::Node::GetHierarchyTier() const
{
    return m_tier;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy::Node *LArHierarchyHelper::MCMatches::GetMC() const
{
    return m_pMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::MCMatches::GetRecoMatches() const
{
    return m_recoNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline size_t LArHierarchyHelper::MCMatches::GetNRecoMatches() const
{
    return m_recoNodes.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCMatchesVector &LArHierarchyHelper::MatchInfo::GetMatches(const pandora::MCParticle *const pRoot) const
{
    if (m_matches.find(pRoot) == m_matches.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return m_matches.at(pRoot);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::MatchInfo::GetUnmatchedReco() const
{
    return m_unmatchedReco;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy &LArHierarchyHelper::MatchInfo::GetMCHierarchy() const
{
    return m_mcHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy &LArHierarchyHelper::MatchInfo::GetRecoHierarchy() const
{
    return m_recoHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::QualityCuts &LArHierarchyHelper::MatchInfo::GetQualityCuts() const
{
    return m_qualityCuts;
}

} // namespace lar_content

#endif // #ifndef LAR_HIERARCHY_HELPER_H
