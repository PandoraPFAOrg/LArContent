/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverHitsTool.h
 *
 *  @brief  Header file for the plane solver-based hit creation tool.
 *
 *  $Log: $
 */
#ifndef LAR_PLANE_SOLVER_HITS_TOOL_H
#define LAR_PLANE_SOLVER_HITS_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar_content
{

/**
 *  @brief  PlaneSolverHitsTool class
 */
class PlaneSolverHitsTool : public HitCreationBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    PlaneSolverHitsTool();

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo the address of the pfo
     *  @param  inputTwoDHits the vector of input two dimensional hits
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_planeSolverContextName; ///< Name of the LArPlaneContextObject produced by PlaneSolverAlgorithm
};

} // namespace lar_content

#endif // #ifndef LAR_PLANE_SOLVER_HITS_TOOL_H
