/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MASTER_ALGORITHM_H
#define LAR_MASTER_ALGORITHM_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingBaseTool.h"
#include "larpandoracontent/LArControlFlow/SliceIdBaseTool.h"
#include "larpandoracontent/LArControlFlow/SliceSelectionBaseTool.h"

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h"

#include <unordered_map>

namespace lar_content
{

class LArMCParticleFactory;

typedef std::vector<pandora::CaloHitList> SliceVector;
typedef std::vector<pandora::PfoList> SliceHypotheses;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MasterAlgorithm class
 */
class MasterAlgorithm : public pandora::ExternallyConfiguredAlgorithm, public StitchingPfoOperations
{
public:
    /**
     *  @brief  Default constructor
     */
    MasterAlgorithm();

    /**
     *  @brief  External steering parameters class
     */
    class ExternalSteeringParameters : public pandora::ExternalParameters
    {
    public:
        pandora::InputBool m_shouldRunAllHitsCosmicReco;  ///< Whether to run all hits cosmic-ray reconstruction
        pandora::InputBool m_shouldRunStitching;          ///< Whether to stitch cosmic-ray muons crossing between volumes
        pandora::InputBool m_shouldRunCosmicHitRemoval;   ///< Whether to remove hits from tagged cosmic-rays
        pandora::InputBool m_shouldRunSlicing;            ///< Whether to slice events into separate regions for processing
        pandora::InputBool m_shouldRunNeutrinoRecoOption; ///< Whether to run neutrino reconstruction for each slice
        pandora::InputBool m_shouldRunCosmicRecoOption;   ///< Whether to run cosmic-ray reconstruction for each slice
        pandora::InputBool m_shouldPerformSliceId;        ///< Whether to identify slices and select most appropriate pfos
        pandora::InputBool m_printOverallRecoStatus;      ///< Whether to print current operation status messages
    };

    /**
     *  @brief  Shift a Pfo hierarchy by a specified x0 value
     *
     *  @param  pParentPfo the address of the parent pfo
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  x0 the x0 correction relative to the input pfo
     */
    void ShiftPfoHierarchy(const pandora::ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const override;

    /**
     *  @brief  Stitch together a pair of pfos
     *
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     */
    void StitchPfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete,
        PfoToLArTPCMap &pfoToLArTPCMap) const override;

protected:
    /**
     *  @brief  LArTPCHitList class
     */
    class LArTPCHitList
    {
    public:
        pandora::CaloHitList m_allHitList;       ///< The list of all hits originating from a given LArTPC
        pandora::CaloHitList m_truncatedHitList; ///< The list of hits confined within LArTPC boundaries for given beam t0 value
    };

    typedef std::map<unsigned int, LArTPCHitList> VolumeIdToHitListMap;

    pandora::StatusCode Run() override;

    /**
     *  @brief  Initialize pandora worker instances
     */
    pandora::StatusCode InitializeWorkerInstances();

    /**
     *  @brief  Copy mc particles in the named input list to all pandora worker instances
     */
    pandora::StatusCode CopyMCParticles() const;

    /**
     *  @brief  Get the mapping from lar tpc volume id to lists of all hits, and truncated hits
     *
     *  @param  volumeIdToHitListMap to receive the populated volume id to hit list map
     */
    pandora::StatusCode GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const;

    /**
     *  @brief  Run the cosmic-ray reconstruction worker instances
     *
     *  @param  volumeIdToHitListMap the volume id to hit list map
     */
    pandora::StatusCode RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const;

    /**
     *  @brief  Recreate cosmic-ray pfos (created by worker instances) in the master instance
     *
     *  @param  pfoToLArTPCMap to receive the populated pfo to lar tpc map
     */
    pandora::StatusCode RecreateCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const;

    /**
     *  @brief  Stitch together cosmic-ray pfos crossing between adjacent lar tpcs
     *
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedPfosToX0Map to receive the map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     */
    pandora::StatusCode StitchCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) const;

    /**
     *  @brief  Tag clear, unambiguous cosmic-ray pfos
     *
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     *  @param  clearCosmicRayPfos to receive the list of clear cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous cosmic-ray pfos for further analysis
     */
    pandora::StatusCode TagCosmicRayPfos(
        const PfoToFloatMap &stitchedPfosToX0Map, pandora::PfoList &clearCosmicRayPfos, pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Run cosmic-ray hit removal, freeing hits in ambiguous pfos for further processing
     *
     *  @param  ambiguousPfos the list of ambiguous cosmic-ray pfos
     */
    pandora::StatusCode RunCosmicRayHitRemoval(const pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Run the event slicing procedures, dividing available hits up into distinct 3D regions
     *
     *  @param  volumeIdToHitListMap the volume id to hit list map
     *  @param  sliceVector to receive the populated slice vector
     */
    pandora::StatusCode RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const;

    /**
     *  @brief  Process each slice under different reconstruction hypotheses
     *
     *  @param  sliceVector the slice vector
     *  @param  nuSliceHypotheses to receive the vector of slice neutrino hypotheses
     *  @param  crSliceHypotheses to receive the vector of slice cosmic-ray hypotheses
     */
    pandora::StatusCode RunSliceReconstruction(SliceVector &sliceVector, SliceHypotheses &nuSliceHypotheses, SliceHypotheses &crSliceHypotheses) const;

    /**
     *  @brief  Examine slice hypotheses to identify the most appropriate to provide in final event output
     *
     *  @param  nuSliceHypotheses the vector of slice neutrino hypotheses
     *  @param  crSliceHypotheses the vector of slice cosmic-ray hypotheses
     */
    pandora::StatusCode SelectBestSliceHypotheses(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses) const;

    /**
     *  @brief  Reset all worker instances
     */
    pandora::StatusCode Reset() override;

    /**
     *  @brief  Copy a specified calo hit to the provided pandora instance
     *
     *  @param  pPandora the address of the target pandora instance
     *  @param  pCaloHit the address of the calo hit
     */
    pandora::StatusCode Copy(const pandora::Pandora *const pPandora, const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Copy a specified mc particle to the provided pandora instance
     *
     *  @param  pPandora the address of the target pandora instance
     *  @param  pMCParticle the address of the mc particle
     *  @param  pMCParticleFactory the address of the mc particle factory, allowing decoration of instances with information beyond that expected by sdk
     */
    pandora::StatusCode Copy(const pandora::Pandora *const pPandora, const pandora::MCParticle *const pMCParticle,
        const LArMCParticleFactory *const pMCParticleFactory) const;

    /**
     *  @brief  Recreate a specified list of pfos in the current pandora instance
     *
     *  @param  inputPfoList the input pfo list
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::PfoList &inputPfoList, pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Recreate a specified pfo in the current pandora instance
     *
     *  @param  pInputPfo the input pfo
     *  @param  pNewParentPfo the new parent of the new output pfo (nullptr if none)
     *  @param  newPfoList to receive the list of new pfos
     */
    virtual pandora::StatusCode Recreate(const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ParticleFlowObject *const pNewParentPfo, pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Create a new calo hit in the current pandora instance, based upon the provided input calo hit
     *
     *  @param  pInputCaloHit the address of the input calo hit
     *  @param  pParentCaloHit the address of the parent calo hit
     *
     *  @return the address of the new calo hit
     */
    const pandora::CaloHit *CreateCaloHit(const pandora::CaloHit *const pInputCaloHit, const pandora::CaloHit *const pParentCaloHit) const;

    /**
     *  @brief  Create a new cluster in the current pandora instance, based upon the provided input cluster
     *
     *  @param  pInputCluster the address of the input cluster
     *  @param  newCaloHitList the list of calo hits for the new cluster
     *  @param  newIsolatedCaloHitList the list of isolated calo hits for the new cluster
     *
     *  @return the address of the new cluster
     */
    const pandora::Cluster *CreateCluster(const pandora::Cluster *const pInputCluster, const pandora::CaloHitList &newCaloHitList,
        const pandora::CaloHitList &newIsolatedCaloHitList) const;

    /**
     *  @brief  Create a new vertex in the current pandora instance, based upon the provided input vertex
     *
     *  @param  pInputVertex the address of the input vertex
     *
     *  @return the address of the new vertex
     */
    const pandora::Vertex *CreateVertex(const pandora::Vertex *const pInputVertex) const;

    /**
     *  @brief  Create a new pfo in the current pandora instance, based upon the provided input pfo
     *
     *  @param  pInputPfo the address of the input pfo
     *  @param  newClusterList the list of clusters for the new pfo
     *  @param  newVertexList the list of vertices for the new pfo
     *
     *  @return the address of the new pfo
     */
    const pandora::ParticleFlowObject *CreatePfo(const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ClusterList &newClusterList, const pandora::VertexList &newVertexList) const;

    /**
     *  @brief  Create a pandora worker instance to handle a single LArTPC
     *
     *  @param  larTPC the lar tpc
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *  @param  name the pandora instance name
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPC &larTPC, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile, const std::string &name) const;

    /**
     *  @brief  Create a pandora worker instance to handle a number of LArTPCs
     *
     *  @param  larTPCMap the lar tpc map
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *  @param  name the pandora instance name
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPCMap &larTPCMap, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile, const std::string &name) const;

    /**
     *  @brief  Append the readout volume parameters from a given LArTPC into a readout volume parameters vector, offsetting each readout volume's
     *          id to maintain uniqueness when merging multiple LArTPCs into one worker instance
     *
     *  @param  larTPC the source LArTPC
     *  @param  idOffset the offset to apply to each readout volume's id
     *  @param  readoutVolumeParametersVector the vector to append to
     */
    void AppendReadoutVolumeParameters(const pandora::LArTPC &larTPC, const unsigned int idOffset,
        object_creation::LArReadoutVolumeParametersVector &readoutVolumeParametersVector) const;

    /**
     *  @brief  Register custom content, such as algorithms or algorithm tools, with a specified pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     */
    virtual pandora::StatusCode RegisterCustomContent(const pandora::Pandora *const pPandora) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    /**
     *  @brief  Read settings from external steering parameters block, if present, otherwise from xml as standard
     *
     *  @param  pExternalParameters the address of the external parameters (if present)
     *  @param  inputBool the input boolean value, from the external parameters (if present)
     *  @param  xmlHandle the xml handle
     *  @param  xmlTag the xml tag
     *  @param  outputBool to receive the output boolean value
     */
    pandora::StatusCode ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters,
        const pandora::InputBool inputBool, const pandora::TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool);

    bool m_workerInstancesInitialized; ///< Whether all worker instances have been initialized

    unsigned int m_larCaloHitVersion; ///< The LArCaloHit version for LArCaloHitFactory

    bool m_shouldRunAllHitsCosmicReco;  ///< Whether to run all hits cosmic-ray reconstruction
    bool m_shouldRunStitching;          ///< Whether to stitch cosmic-ray muons crossing between volumes
    bool m_shouldRunCosmicHitRemoval;   ///< Whether to remove hits from tagged cosmic-rays
    bool m_shouldRunSlicing;            ///< Whether to slice events into separate regions for processing
    bool m_shouldRunNeutrinoRecoOption; ///< Whether to run neutrino reconstruction for each slice
    bool m_shouldRunCosmicRecoOption;   ///< Whether to run cosmic-ray reconstruction for each slice
    bool m_shouldPerformSliceId;        ///< Whether to identify slices and select most appropriate pfos
    bool m_printOverallRecoStatus;      ///< Whether to print current operation status messages
    bool m_visualizeOverallRecoStatus;  ///< Whether to display results of current operations
    bool m_shouldRemoveOutOfTimeHits;   ///< Whether to remove out of time hits

    PandoraInstanceList m_crWorkerInstances;          ///< The list of cosmic-ray reconstruction worker instances
    const pandora::Pandora *m_pSlicingWorkerInstance; ///< The slicing worker instance
    const pandora::Pandora *m_pSliceNuWorkerInstance; ///< The per-slice neutrino reconstruction worker instance
    const pandora::Pandora *m_pSliceCRWorkerInstance; ///< The per-slice cosmic-ray reconstruction worker instance
    mutable std::map<unsigned int, unsigned int> m_daughterVolumeIdOffsetMap; ///< Maps original LArTPC volume id -> id offset for merged instances

    bool m_fullWidthCRWorkerWireGaps;        ///< Whether wire-type line gaps in cosmic-ray worker instances should cover all drift time
    bool m_passMCParticlesToWorkerInstances; ///< Whether to pass mc particle details (and links to calo hits) to worker instances

    typedef std::vector<StitchingBaseTool *> StitchingToolVector;
    typedef std::vector<CosmicRayTaggingBaseTool *> CosmicRayTaggingToolVector;
    typedef std::vector<SliceIdBaseTool *> SliceIdToolVector;
    typedef std::vector<SliceSelectionBaseTool *> SliceSelectionToolVector;

    StitchingToolVector m_stitchingToolVector;               ///< The stitching tool vector
    CosmicRayTaggingToolVector m_cosmicRayTaggingToolVector; ///< The cosmic-ray tagging tool vector
    SliceIdToolVector m_sliceIdToolVector;                   ///< The slice id tool vector
    SliceSelectionToolVector m_sliceSelectionToolVector;     ///< The slice selection tool vector

    std::string m_filePathEnvironmentVariable; ///< The environment variable providing a list of paths to xml files
    std::string m_crSettingsFile;              ///< The cosmic-ray reconstruction settings file
    std::string m_nuSettingsFile;              ///< The neutrino reconstruction settings file
    std::string m_slicingSettingsFile;         ///< The slicing settings file

    std::string m_inputMCParticleListName;  ///< The input mc particle list name
    std::string m_inputHitListName;         ///< The input hit list name
    std::string m_recreatedPfoListName;     ///< The output recreated pfo list name
    std::string m_recreatedClusterListName; ///< The output recreated cluster list name
    std::string m_recreatedVertexListName;  ///< The output recreated vertex list name

    float m_inTimeMaxX0;                   ///< Cut on X0 to determine whether particle is clear cosmic ray
    LArCaloHitFactory m_larCaloHitFactory; ///< Factory for creating LArCaloHits during hit copying
};

} // namespace lar_content

#endif // #ifndef LAR_MASTER_ALGORITHM_H
