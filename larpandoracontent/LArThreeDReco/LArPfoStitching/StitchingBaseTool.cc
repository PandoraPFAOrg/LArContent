/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.cc
 *
 *  @brief  Implementation file for the stitching tool base class
 *
 *  $Log: $
 */

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h"

using namespace pandora;

namespace lar_content
{

void StitchingBaseTool::Run(const Algorithm *const pAlgorithm, const PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap,
    PfoToFloatMap &stitchedPfosToX0Map)
{
    const StitchingPfoOperations *const pStitchingOperations{dynamic_cast<const StitchingPfoOperations *>(pAlgorithm)};
    if (!pStitchingOperations)
    {
        std::cout << "StitchingBaseTool: Calling algorithm must implement StitchingPfoOperations" << std::endl;
        PANDORA_THROW(STATUS_CODE_INVALID_PARAMETER);
    }

    this->RunStitching(pAlgorithm, pStitchingOperations, pMultiPfoList, pfoToLArTPCMap, stitchedPfosToX0Map);
}

} // namespace lar_content
