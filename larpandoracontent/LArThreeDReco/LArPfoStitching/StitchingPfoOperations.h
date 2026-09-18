/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h
 *
 *  @brief  Header file for the stitching pfo operations interface class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_PFO_OPERATIONS_H
#define LAR_STITCHING_PFO_OPERATIONS_H 1

#include "Pandora/AlgorithmHeaders.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  StitchingPfoOperations class
 */
class StitchingPfoOperations
{
public:
    typedef std::unordered_map<const pandora::ParticleFlowObject *, const pandora::LArTPC *> PfoToLArTPCMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, float> PfoToFloatMap;

    virtual ~StitchingPfoOperations() = default;

    /**
     *  @brief  Shift a Pfo hierarchy by a specified x0 value
     *
     *  @param  pParentPfo the address of the parent pfo
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  x0 the x0 correction relative to the input pfo
     */
    virtual void ShiftPfoHierarchy(const pandora::ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const = 0;

    /**
     *  @brief  Stitch together a pair of pfos
     *
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     */
    virtual void StitchPfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete,
        PfoToLArTPCMap &pfoToLArTPCMap) const = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_PFO_OPERATIONS_H
