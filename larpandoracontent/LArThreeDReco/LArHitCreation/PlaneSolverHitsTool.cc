/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverHitsTool.cc
 *
 *  @brief  Implementation of the plane solver-based hit creation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverHitsTool.h"

using namespace pandora;

namespace lar_content
{

PlaneSolverHitsTool::PlaneSolverHitsTool() :
    m_planeSolverContextName("PlaneContext")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverHitsTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
    const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector)
{

    const LArPlaneContextObject *pPlaneContextObject{
        dynamic_cast<const LArPlaneContextObject *>(PandoraContentApi::GetEventContextObject(*pAlgorithm, m_planeSolverContextName))};

    if (!pPlaneContextObject)
        return;

    for (auto const *pCaloHit2D : inputTwoDHits)
    {
        const LArPlaneContextObject::HitTriplet *pTriplet{pPlaneContextObject->GetHitTriplet(pCaloHit2D)};
        if (!pTriplet)
            continue;

        std::unordered_map<HitType, const CaloHit *> hitTypeToCaloHitMap;
        hitTypeToCaloHitMap[TPC_VIEW_U] = pTriplet->m_uHit;
        hitTypeToCaloHitMap[TPC_VIEW_V] = pTriplet->m_vHit;
        hitTypeToCaloHitMap[TPC_VIEW_W] = pTriplet->m_wHit;

        const HitType hitType{pCaloHit2D->GetHitType()};
        const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
        const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

        bool canUseHit1{false};
        if (hitTypeToCaloHitMap[hitType1] != nullptr)
            canUseHit1 = std::find(inputTwoDHits.begin(), inputTwoDHits.end(), hitTypeToCaloHitMap[hitType1]) != inputTwoDHits.end();

        bool canUseHit2{false};
        if (hitTypeToCaloHitMap[hitType2] != nullptr)
            canUseHit2 = std::find(inputTwoDHits.begin(), inputTwoDHits.end(), hitTypeToCaloHitMap[hitType2]) != inputTwoDHits.end();

        // If we don't have any hits then move on
        if (!canUseHit1 && !canUseHit2)
            continue;

        CartesianPointVector fitPositions1, fitPositions2;
        if (canUseHit1)
            fitPositions1.emplace_back(hitTypeToCaloHitMap[hitType1]->GetPositionVector());
        if (canUseHit2)
            fitPositions2.emplace_back(hitTypeToCaloHitMap[hitType2]->GetPositionVector());

        ProtoHit protoHit(pCaloHit2D);
        this->GetBestPosition3D(hitType1, hitType2, fitPositions1, fitPositions2, protoHit);

        if (protoHit.IsPositionSet() && (protoHit.GetChi2() < m_chiSquaredCut))
            protoHitVector.emplace_back(protoHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PlaneSolverHitsTool::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PlaneSolverContextName", m_planeSolverContextName));

    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
