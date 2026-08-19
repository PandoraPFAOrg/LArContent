/**
 *  @file   larpandoracontent/LArMonitoring/OpHitMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/OpHitMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

OpHitMonitoringAlgorithm::OpHitMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

OpHitMonitoringAlgorithm::~OpHitMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitList));

    const CaloHitList *pOpHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "OpticalHitList", pOpHitList));

    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &(*pCaloHitList), "Charge", GRAY));

    for (const CaloHit *const pCaloHit : *pOpHitList)
    {
        const LArOpHit *const pOpHit{dynamic_cast<const LArOpHit *>(pCaloHit)};
        if (!pOpHit)
            continue;
        const CartesianVector pos{pOpHit->GetPositionVector()};
        const float pe{pOpHit->GetInputEnergy() / 100.f};
        const int peRound{static_cast<int>(std::ceil(pe))};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "U true vertex", ORANGE, peRound));
        std::cout << "OpHit(" << pos.GetX() << "," << pos.GetY() << "," << pos.GetZ() << "): " << pOpHit->GetInputEnergy() << std::endl;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
