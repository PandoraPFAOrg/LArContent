/**
 *  @file   larpandoracontent/LArControlFlow/CosmicRayTaggingTool.h
 *
 *  @brief  Header file for the cosmic-ray tagging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TAGGING_TOOL_H
#define LAR_COSMIC_RAY_TAGGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingBaseTool.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CosmicRayTaggingTool class
 */
class CosmicRayTaggingTool : public CosmicRayTaggingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRayTaggingTool();

    pandora::StatusCode Initialize();
    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm);

protected:
    /**
     *  @brief  Class to encapsulate the logic required determine if a Pfo should or shouldn't be tagged as a cosmic ray
     */
    class CRCandidate
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pandora the relevant pandora instance
         *  @param  pPfo the address of the candidate pfo
         *  @param  slice the slice id
         */
        CRCandidate(const pandora::Pandora &pandora, const pandora::ParticleFlowObject *const pPfo, const unsigned int sliceId);

        const pandora::ParticleFlowObject *const m_pPfo; ///< Address of the candidate Pfo
        unsigned int m_sliceId;                          ///< Slice ID
        bool m_canFit;                                   ///< If there are a sufficient number of 3D hits to perform a fitting
        pandora::CartesianVector m_endPoint1;            ///< First fitted end point in 3D
        pandora::CartesianVector m_endPoint2;            ///< Second fitted end point in 3D
        double m_length;                                 ///< Straight line length of the linear fit
        double m_curvature;                              ///< Measure of the curvature of the track
        double m_theta;                                  ///< Direction made with vertical

    private:
        /**
         *  @brief  Calculate all variables which require a fit
         *
         *  @param  slidingFitResult the three dimensional sliding fit result
         */
        void CalculateFitVariables(const ThreeDSlidingFitResult &slidingFitResult);
    };

    typedef std::list<CRCandidate> CRCandidateList;

    /**
     *  @brief  Get the 3D calo hit cluster associated with a given Pfo, and check if it has sufficient hits
     *
     *  @param  pPfo input Pfo
     *  @param  pCluster3D to receive the address of the 3D cluster
     *
     *  @return whether the Pfo has sufficient hits?
     */
    bool GetValid3DCluster(const pandora::ParticleFlowObject *const pPfo, const pandora::Cluster *&pCluster3D) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, pandora::PfoList> PfoToPfoListMap;

    /**
     *  @brief  Get mapping between Pfos that are associated with it other by pointing
     *
     *  @param  parentCosmicRayPfos input list of Pfos
     *  @param  pfoAssociationsMap to receive the output mapping between associated Pfos
     */
    void GetPfoAssociations(const pandora::PfoList &parentCosmicRayPfos, PfoToPfoListMap &pfoAssociationMap) const;

    /**
     *  @brief  Check whethe two Pfo endpoints are associated by distance of closest approach
     *
     *  @param  endPoint1 position vector of an endpoint of Pfo 1
     *  @param  endDir1 direction vector of an endpoint of Pfo 1
     *  @param  endPoint2 position vector of an endpoint of Pfo 2
     *  @param  endDir2 direction vector of an endpoint of Pfos
     *
     *  @return whether the Pfos are associated
     */
    bool CheckAssociation(const pandora::CartesianVector &endPoint1, const pandora::CartesianVector &endDir1,
        const pandora::CartesianVector &endPoint2, const pandora::CartesianVector &endDir2) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, unsigned int> PfoToSliceIdMap;

    /**
     *  @brief  Break the event up into slices of associated Pfos
     *
     *  @param  parentCosmicRayPfos input list of Pfos
     *  @param  pfoAssociationMap mapping between Pfos and other associated Pfos
     *  @param  pfoToSliceIdMap to receive the mapping between Pfos and their slice ID
     */
    void SliceEvent(const pandora::PfoList &parentCosmicRayPfos, const PfoToPfoListMap &pfoAssociationMap, PfoToSliceIdMap &pfoToSliceIdMap) const;

    /**
     *  @brief  Fill a slice iteratively using Pfo associations
     *
     *  @param  pPfo Pfo to add to the slice
     *  @param  pfoAssociationMap mapping between Pfos and other associated Pfos
     *  @param  slice the slice to add Pfos to
     */
    void FillSlice(const pandora::ParticleFlowObject *const pPfo, const PfoToPfoListMap &pfoAssociationMap, pandora::PfoList &slice) const;

    /**
     *  @brief  Make a list of CRCandidates
     *
     *  @param  parentCosmicRayPfos input list of Pfos
     *  @param  pfoToSliceIdMap input mapping between Pfos and their slice id
     *  @param  candidates to receive the output list of CRCandidates
     */
    void GetCRCandidates(const pandora::PfoList &parentCosmicRayPfos, const PfoToSliceIdMap &pfoToSliceIdMap, CRCandidateList &candidates) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, bool> PfoToBoolMap;

    /**
     *  @brief  Check if each candidate is "in time"
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToInTimeMap output mapping between candidates Pfos and if they are in time
     */
    void CheckIfInTime(const CRCandidateList &candidates, PfoToBoolMap &pfoToInTimeMap) const;

    /**
     *  @brief  Check if each candidate is "contained" (contained = no associations to Y or Z detector faces, but the Pfo could still enter or exit by and X-face)
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsContainedMap output mapping between candidates Pfos and if they are contained
     */
    void CheckIfContained(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsContainedMap) const;

    /**
     *  @brief  Check if each candidate is "top to bottom"
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsTopToBottomMap output mapping between candidates Pfos and if they are top to bottom
     */
    void CheckIfTopToBottom(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsTopToBottomMap) const;

    typedef std::set<unsigned int> UIntSet;
    typedef std::unordered_map<int, bool> IntBoolMap;

    /**
     *  @brief  Get the slice indices which contain a likely neutrino Pfo
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToInTimeMap input map between Pfo and if in time
     *  @param  pfoToIsContainedMap input map between Pfo and if contained
     *  @param  neutrinoSliceSet output set of slice indices containing a likely neutrino Pfo
     */
    void GetNeutrinoSlices(const CRCandidateList &candidates, const PfoToBoolMap &pfoToInTimeMap, const PfoToBoolMap &pfoToIsContainedMap,
        UIntSet &neutrinoSliceSet) const;

    /**
     *  @brief  Tag Pfos which are likely to be a CR muon
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToInTimeMap input map between Pfo and if in time
     *  @param  pfoToIsTopToBottomMap input mapping between candidate Pfos and if they are top to bottom
     *  @param  neutrinoSliceSet input set of slice indices containing a likely neutrino Pfo
     *  @param  pfoToIsLikelyCRMuonMap to receive the output mapping between Pfos and a boolean deciding if they are likely a CR muon
     */
    void TagCRMuons(const CRCandidateList &candidates, const PfoToBoolMap &pfoToInTimeMap, const PfoToBoolMap &pfoToIsTopToBottomMap,
        const UIntSet &neutrinoSliceSet, PfoToBoolMap &pfoToIsLikelyCRMuonMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::pair<const ThreeDSlidingFitResult, const ThreeDSlidingFitResult> SlidingFitPair;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, SlidingFitPair> PfoToSlidingFitsMap;
    typedef std::vector<pandora::PfoList> SliceList;

    /**
     *  @brief  Choose a set of cuts using a keyword - "cautious" = remove as few neutrinos as possible
     *          "nominal" = optimised to maximise CR removal whilst preserving neutrinos
     *          "aggressive" = remove CR muons and allow more neutrinos to be tagged
     */
    std::string m_cutMode;

    float m_angularUncertainty;    ///< The uncertainty in degrees for the angle of a Pfo
    float m_positionalUncertainty; ///< The uncertainty in cm for the position of Pfo endpoint in 3D
    float m_maxAssociationDist; ///< The maximum distance from endpoint to point of closest approach, typically a multiple of LAr radiation length

    unsigned int m_minimumHits; ///< The minimum number of hits for a Pfo to be considered

    float m_inTimeMargin; ///< The maximum distance outside of the physical detector volume that a Pfo may be to still be considered in time
    float m_inTimeMaxX0;  ///< The maximum pfo x0 (determined from shifted vertex) to allow pfo to still be considered in time
    float m_marginY;      ///< The minimum distance from a detector Y-face for a Pfo to be associated
    float m_marginZ;      ///< The minimum distance from a detector Z-face for a Pfo to be associated
    float m_maxNeutrinoCosTheta; ///< The maximum cos(theta) that a Pfo can have to be classified as a likely neutrino
    float m_minCosmicCosTheta;   ///< The minimum cos(theta) that a Pfo can have to be classified as a likely CR muon
    float m_maxCosmicCurvature;  ///< The maximum curvature that a Pfo can have to be classified as a likely CR muon

    float m_face_Xa; ///< Anode      X face
    float m_face_Xc; ///< Cathode    X face
    float m_face_Yb; ///< Bottom     Y face
    float m_face_Yt; ///< Top        Y face
    float m_face_Zu; ///< Upstream   Z face
    float m_face_Zd; ///< Downstream Z face
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TAGGING_TOOL_H
