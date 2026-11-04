/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

namespace lar_dl_content
{

/**
 *  @brief  DlHitTrackShowerIdAlgorithm class
 */
class DlHitTrackShowerIdAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlHitTrackShowerIdAlgorithm();

    virtual ~DlHitTrackShowerIdAlgorithm();

private:
    typedef std::unordered_map<const pandora::MCParticle *, const pandora::MCParticle *> MCFoldingMap;
    typedef std::vector<size_t> IndexVector;

    enum Category
    {
        UNINITIALISED = 0,
        MIP = 1,
        HIP = 2,
        EM_TRACK = 3,
        EM_SHOWER = 4,
        DIFFUSE = 5
    };

    struct Bounds
    {
        float xMin = std::numeric_limits<float>::max();
        float xMax = std::numeric_limits<float>::lowest();
        float zMin = std::numeric_limits<float>::max();
        float zMax = std::numeric_limits<float>::lowest();
        float rMax = 0;
        float xRange = 1.f;
        float zRange = 1.f;
    };

    pandora::StatusCode Run();

    /**
     *  @brief  Prepare training sample for the network
     *
     *  @return The StatusCode indicating success or failure in producing a training sample
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Run network inference
     *
     *  @return The StatusCode indicating success or failure in running inference
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Fold EM showers to their leading particle.
     *
     *  @param mcHitMap The map of MC particles to their hits
     *  @param mcFoldingMap The map of MC particles to their hits
     *  @param leadingHitMap The map of leading MC particles to their consolidated hits
     */
    void FoldToLeading(const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, const MCFoldingMap &mcFoldingMap,
        lar_content::LArMCParticleHelper::MCContributionMap &leadingHitMap) const;

    /**
     *  @brief  Consolidate instances of hits from MC particles where folding is appropriate.
     *          Delta rays often own hits along muon or pion tracks, and these should be
     *          owned by the parent particle instead.
     *
     *  @param mcHitMap The map of MC particles to their hits
     *  @param instanceHitMap The map of MC particles to their consolidated hits
     */
    void ConsolidateInstances(const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap,
        lar_content::LArMCParticleHelper::MCContributionMap &instanceHitMap) const;

    /**
     *  @brief  Determine if an MC particle is diffuse
     *
     *  @param pLArMC The MC particle to check
     *  @param mcToHitsMap The map of MC particles to their hits, used to determine if the particle is diffuse
     *  @return true if the MC particle is diffuse, false otherwise
     */
    bool IsDiffuse(const lar_content::LArMCParticle *const pLArMC, const lar_content::LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /**
     *  @brief  Determine if an MC particle is highly ionising
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is highly ionising, false otherwise
     */
    bool IsHip(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is minimally ionising
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is minimally ionising, false otherwise
     */
    bool IsMip(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is a shower-like EM particle (essentially, does it Bremm?)
     *
     *  @param pLArMC The MC particle to check
     *  @param mcHitMap The map of MC particles to their hits, used to determine if the particle undergoes Bremmstrahlung
     *  @return true if the MC particle is a shower, false otherwise
     */
    bool IsShowerLikeEM(const lar_content::LArMCParticle *const pLArMC, const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap) const;

    /**
     *  @brief  Determine if an MC particle is track-like EM particle (effectively, does it only ionise? In practice, this function assumes we've
     *          already ruled out diffuse and shower-like EM particles, so this is effectively a PDG check)
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is track-like EM, false otherwise
     */
    bool IsTrackLikeEM(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle undergoes Bremmstrahlung, with charge deposition in a Bremm branch.
     *
     *  @param pLArMC The MC particle to check
     *  @param mcHitMap The map of MC particles to their hits, used to determine if there is charge deposition along a Bremm branch
     *  @param onBremmBranch Whether we're in a branch that has undergone Bremmstrahlung.
     *  @return true if the MC particle undergoes Bremmstrahlung, false otherwise
     */
    bool UndergoesBremmstrahlung(const lar_content::LArMCParticle *const pLArMC, const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap,
        bool onBremmBranch) const;

    /**
     *  @brief  Determine if all hits in a particle are a result of non-Bremmstrahlung induced Compton scatters.
     *
     *  @param pLArMC The MC particle to check
     *  @param mcHitMap The map of MC particles to their hits, used to determine if there is charge deposition due to Compton scattering
     *  @param onBremmBranch Whether we're in a branch that has undergone Bremmstrahlung.
     *  @return true if all MC particle hits are due to non-Bremmstrahlung induced Compton scatters, false otherwise
     */
    bool AllHitsCompton(const lar_content::LArMCParticle *const pLArMC, const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap,
        bool onBremmBranch) const;

    /**
     *  @brief  Choose the owner of a child MC particle based on its parent and the folding map
     *
     *  @param pMC The MC particle
     *  @param mcFoldingMap The map of MC particles to their folded counterparts
     *  @param mcHitMap The map of MC particles to their hits
     *  @param particleOwnedHits The list of hits owned by the particle
     *  @param parentOwnedHits The list of hits owned by the parent particle
     *
     *  @return The parent MC particle if there is allocation to a parent, nullptr otherwise
     */
    const lar_content::LArMCParticle *AllocateHitOwner(const pandora::MCParticle *const pMC, const MCFoldingMap &mcFoldingMap,
        const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, pandora::CaloHitList &particleOwnedHits,
        pandora::CaloHitList &parentOwnedHits) const;

    /**
     *  @brief  Choose the owner of a child MC particle based on its parent and the folding map
     *
     *  @param pMC The MC particle
     *  @param mcHitMap The map of MC particles to their hits
     *  @param particleOwnedHits The list of hits owned by the particle
     *  @param parentOwnedHits The list of hits owned by the parent particle
     *
     *  @return The parent MC particle if there is allocation to a parent, nullptr otherwise
     */
    const lar_content::LArMCParticle *AllocateHitOwner(const pandora::MCParticle *const pMC,
        const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, pandora::CaloHitList &particleOwnedHits,
        pandora::CaloHitList &parentOwnedHits) const;

    /**
     *  @brief  Get planar coordinates of vertex in given view
     *
     *  @param  view The view into which the vertex should be projected
     *  @param  vx The x position of the vertex
     *  @param  vz The z position of the vertex
     *
     *  @return The status code
     */
    pandora::StatusCode GetVertexPlanarCoordinates(pandora::HitType view, float &vx, float &vz) const;

    /**
     *  @brief  Get extrema of hit coordinates in the specified hit list
     *
     *  @param  caloHitList The list from which extrema should be extracted
     *  @param  vx The x position of the vertex
     *  @param  vz The z position of the vertex
     *  @param  bounds The output structure of extremal coordinate information
     */
    void GetCoordinateExtrema(const pandora::CaloHitList &caloHitList, [[maybe_unused]]const float vx, [[maybe_unused]]const float vz,
        Bounds &bounds) const;

    /**
     *  @brief  Populate input vectors for the deep learning network from the specified hit list
     *
     *  @param  caloHitList The list from which to populate the input vectors
     *  @param  bounds The coordinate bounds structure
     *  @param  vx The x position of the vertex
     *  @param  vz The z position of the vertex
     *  @param  x_rel The vertex relative x coordinate vector
     *  @param  z_rel The vertex relative z coordinate vector
     *  @param  x_abs The absolute x coordinate vector
     *  @param  z_abs The absolute z coordinate vector
     *  @param  rr The radial coordinate vector
     *  @param  cosTheta The cosine of the theta angle vector
     *  @param  sinTheta The sine of the theta angle vector
     *  @param  wirePitch The vector of wire pitches
     *  @param  wireAngle The vector of wire angles
     *  @param  adc The ADC value vector
     *  @param  width The width vector
     *  @param  sortedCaloHitList The output sorted calo hit list
     */
    void PopulateInputVectors(const pandora::CaloHitList &caloHitList, const Bounds &bounds, const float vx, const float vz, pandora::FloatVector &x_rel,
        pandora::FloatVector &z_rel, pandora::FloatVector &x_abs, pandora::FloatVector &z_abs, pandora::FloatVector &rr, pandora::FloatVector &cosTheta,
        pandora::FloatVector &sinTheta, pandora::FloatVector &wirePitch, pandora::FloatVector &wireAngle, pandora::FloatVector &adc, pandora::FloatVector &width,
        pandora::CaloHitVector &sortedCaloHitList) const;

    void TraverseChildren(const pandora::MCParticle *const pRoot, const lar_content::LArMCParticleHelper::MCContributionMap &mcToHitsMap, const std::string &tab) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< Name of input calo hit list
    bool m_trainingMode; ///< Training mode
    int m_maxSeqLen; ///< Maximum sequence length when using polar coordinates
    std::string m_rootFileName; ///< Name of the ROOT file to save the training sample
    std::string m_rootTreeName; ///< Name of the ROOT tree to save the training sample
    std::string m_vertexListName; ///< Name of the vertex list to use for vertex-relative coordinates
    LArDLHelper::TorchModel m_model; ///< The model
};

} // namespace lar_dl_content

#endif // LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
