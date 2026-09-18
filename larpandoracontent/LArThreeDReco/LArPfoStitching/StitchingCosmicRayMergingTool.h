/**
 *  @file   LArContent/LArThreeDReco/LArPfoStitching/StitchingCosmicRayMergingTool.h
 *
 *  @brief  Header file for the stitching cosmic ray merging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H
#define LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  StitchingCosmicRayMergingTool class
 */
class StitchingCosmicRayMergingTool : public StitchingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    StitchingCosmicRayMergingTool();

    void RunStitching(const pandora::Algorithm *const pAlgorithm, const StitchingPfoOperations *const pStitchingOperations,
        const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) override;

    /**
     *  @brief  PfoAssociation class
     */
    class PfoAssociation
    {
    public:
        /**
         *  @brief  Vertex enumeration
         */
        enum VertexType
        {
            UNDEFINED = 0,
            INNER = 1,
            OUTER = 2
        };

        /**
         *  @brief  Constructor
         *
         *  @param  parent the parent vertex type
         *  @param  daughter the daughter vertex type
         *  @param  fom the figure of merit
         */
        PfoAssociation(const VertexType parent, const VertexType daughter, const float fom);

        /**
         *  @brief  Get parent
         *
         *  @return the parent
         */
        VertexType GetParent() const;

        /**
         *  @brief  Get daughter
         *
         *  @return the daughter
         */
        VertexType GetDaughter() const;

        /**
         *  @brief  Get figure of merit
         *
         *  @return the figure of merit
         */
        float GetFigureOfMerit() const;

    private:
        VertexType m_parent;   ///< The parent vertex type
        VertexType m_daughter; ///< The daughter vertex type
        float m_fom;           ///< The figure of merit
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    /**
     *  @brief  Select Pfos from the input list of Pfos
     *
     *  @param  pInputPfoList the input list of Pfos
     *  @param  pfoToLArTPCMap the input mapping between Pfos and tpc
     *  @param  outputPfoList the output list of Pfos
     */
    void SelectPfos(const pandora::PfoList *pInputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, pandora::PfoList &outputPfoList) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, LArPointingCluster> ThreeDPointingClusterMap;

    /**
     *  @brief  Build a 3D pointing cluster for each Pfo
     *
     *  @param  inputPfoList the input list of Pfos
     *  @param  pfoToLArTPCMap the input mapping between Pfos and tpc
     *  @param  pointingClusterMap  the mapping between Pfos and their corresponding 3D pointing clusters
     */
    void BuildPointingClusterMaps(
        const pandora::PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, ThreeDPointingClusterMap &pointingClusterMap) const;

    typedef std::unordered_map<const pandora::LArTPC *, pandora::PfoList> LArTPCToPfoMap;

    /**
     *  @brief  Build a list of Pfos for each tpc
     *
     *  @param  inputPfoList the input list of Pfos
     *  @param  pfoToLArTPCMap the input mapping between Pfos and tpc
     *  @param  larTPCToPfoMap the output mapping between tpc and Pfos
     */
    void BuildTPCMaps(const pandora::PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, LArTPCToPfoMap &larTPCToPfoMap) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, PfoAssociation> PfoAssociationMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, PfoAssociationMap> PfoAssociationMatrix;

    /**
     *  @brief  Create associations between Pfos using 3D pointing clusters
     *
     *  @param  larTPCToPfoMap the input mapping between tpc and Pfos
     *  @param  pointingClusterMap the input mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoAssociationMatrix the output matrix of associations between Pfos
     */
    void CreatePfoMatches(const LArTPCToPfoMap &larTPCToPfoMap, const ThreeDPointingClusterMap &pointingClusterMap,
        PfoAssociationMatrix &pfoAssociationMatrix) const;

    /**
     *  @brief  Create associations between Pfos using 3D pointing clusters
     *
     *  @param  larTPC1 the tpc description for the first Pfo
     *  @param  larTPC2 the tpc description for the second Pfo
     *  @param  pPfo1 the first Pfo
     *  @param  pPfo2 the second Pfo
     *  @param  pointingClusterMap the input mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoAssociationMatrix the output matrix of associations between Pfos
     */
    void CreatePfoMatches(const pandora::LArTPC &larTPC1, const pandora::LArTPC &larTPC2, const pandora::ParticleFlowObject *const pPfo1,
        const pandora::ParticleFlowObject *const pPfo2, const ThreeDPointingClusterMap &pointingClusterMap,
        PfoAssociationMatrix &pfoAssociationMatrix) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, pandora::PfoList> PfoMergeMap;

    /**
     *  @brief  Select the best associations between Pfos; create a mapping between associated Pfos, handling any ambiguities
     *
     *  @param  pfoAssociationMatrix the input list of all associations between Pfos
     *  @param  pfoSelectedMatches the output list of good associations between pfos (candidates for merging)
     */
    void SelectPfoMatches(const PfoAssociationMatrix &pfoAssociationMatrix, PfoMergeMap &pfoSelectedMatches) const;

    /**
     *  @brief  Create an initial map of Pfo merges to be made
     *
     *  @param  pfoMatches the input list of associated Pfos
     *  @param  pfoMerges the output mapping between Pfos to enlarge and Pfos to delete
     */
    void SelectPfoMerges(const PfoMergeMap &pfoMatches, PfoMergeMap &pfoMerges) const;

    /**
     *  @brief  Collect up associations between Pfos
     *
     *  @param  pSeedPfo the seed Pfo (will be enlarged after merging)
     *  @param  pCurrentPfo the target Pfo (will be deleted after merging)
     *  @param  pfoMerges the list of associations between Pfos
     *  @param  vetoSet the set of associated Pfos that have already been considered
     *  @param  associatedList the output list of associated Pfos
     */
    void CollectAssociatedPfos(const pandora::ParticleFlowObject *const pSeedPfo, const pandora::ParticleFlowObject *const pCurrentPfo,
        const PfoMergeMap &pfoMerges, const pandora::PfoSet &vetoSet, pandora::PfoList &associatedList) const;

    /**
     *  @brief  Identify the vertex Pfo and then re-order the map of merges so that the vertex Pfo will be enlarged
     *
     *  @param  pfoToLArTPCMap the mapping between pfos and tpc
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  inputPfoMerges the input map of Pfo merges
     *  @param  outputPfoMerges the re-ordered map of Pfo merges
     */
    void OrderPfoMerges(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
        const PfoMergeMap &inputPfoMerges, PfoMergeMap &outputPfoMerges) const;

    /**
     *  @brief  Apply X0 corrections, and then stitch together Pfos
     *
     *  @param  pAlgorithm the address of the parent stitching algorithm
     *  @param  pointingClusterMap  the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoMerges the input map of Pfo merges
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     */
    void StitchPfos(const pandora::Algorithm *const pAlgorithm, const StitchingPfoOperations *const pStitchingOperations,
        const ThreeDPointingClusterMap &pointingClusterMap, const PfoMergeMap &pfoMerges, PfoToLArTPCMap &pfoToLArTPCMap,
        PfoToFloatMap &stitchedPfosToX0Map) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, LArPointingCluster::Vertex> PfoToPointingVertexMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, PfoToPointingVertexMap> PfoToPointingVertexMatrix;

    /**
     *  @brief  Shift a pfo given its pfo stitching pair
     *
     *  @param  pPfoToShift the pfo of the stitching pair to shift
     *  @param  pMatchedPfo the pfo of the stitching pair to remain stationary
     *  @param  x0 the distance by which pPfoToShift is to be shifted (direction of shift is determined in method)
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  pfoToPointingVertexMatrix the map [pfo -> map [matched pfo -> pfo stitching vertex]]
     */
    void ShiftPfo(const pandora::Algorithm *const pAlgorithm, const StitchingPfoOperations *const pStitchingOperations,
        const pandora::ParticleFlowObject *const pPfoToShift, const pandora::ParticleFlowObject *const pMatchedPfo, const float x0,
        const PfoToLArTPCMap &pfoToLArTPCMap, const PfoToPointingVertexMatrix &pfoToPointingVertexMatrix) const;

    /**
     *  @brief  Calculate x0 shift for a group of associated Pfos
     *
     *  @param  pfoToLArTPCMap the mapping between pfos and tpc
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoVector the vector of parent Pfos to stitch together
     *  @param  x0 the output x0 value
     *  @param  pfoToPointingVertexMatrix map of pfo to a map of matched pfo and the corresponding pointing vertex used in stitching
     *
     *  @return bool true if all x0 contributions were consistent, false otherwise
     */
    bool CalculateX0(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
        const pandora::PfoVector &pfoVector, float &x0, PfoToPointingVertexMatrix &pfoToPointingVertexMatrix) const;

    bool m_useXcoordinate;
    bool m_alwaysApplyT0Calculation;
    int m_halfWindowLayers;
    float m_minLengthSquared;
    float m_minCosRelativeAngle;
    float m_relaxMinLongitudinalDisplacement; ///< The minimum value of the longitudinal impact parameter for association if both verticies fall in the detector gap
    float m_minLongitudinalDisplacementX;
    float m_maxLongitudinalDisplacementX;
    float m_maxTransverseDisplacement;
    float m_relaxCosRelativeAngle;
    float m_relaxTransverseDisplacement;
    unsigned int m_minNCaloHits3D;
    float m_maxX0FractionalDeviation; ///< The maximum allowed fractional difference of an X0 contribution for matches to be stitched
    float m_boundaryToleranceWidth;   ///< The distance from the APA/CPA boundary inside which the deviation consideration is ignored
    bool m_onlyFinalStatePfos;        ///< Select only final state PFOs
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H
