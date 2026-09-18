/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.h
 *
 *  @brief  Header file for the stitching tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_BASE_TOOL_H
#define LAR_STITCHING_BASE_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h"

namespace lar_content
{

/**
 *  @brief  StitchingBaseTool class
 */
class StitchingBaseTool : public pandora::AlgorithmTool
{
public:
    typedef std::unordered_map<const pandora::ParticleFlowObject *, const pandora::LArTPC *> PfoToLArTPCMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, float> PfoToFloatMap;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pMultiPfoList the list of pfos in multiple lar tpcs
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     */
    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap,
        PfoToFloatMap &stitchedPfosToX0Map);

protected:
    /**
     *  @brief  Run the stitching logic
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pStitchingOperations address of the calling algorithm's stitching operations implementation
     *  @param  pMultiPfoList the list of pfos in multiple lar tpcs
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     */
    virtual void RunStitching(const pandora::Algorithm *const pAlgorithm, const StitchingPfoOperations *const pStitchingOperations,
        const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_BASE_TOOL_H
