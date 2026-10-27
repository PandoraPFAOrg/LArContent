/**
 *  @file   larpandoracontent/LArMonitoring/OpHitMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_OP_HITMONITORING_ALGORITHM_H
#define LAR_OP_HITMONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  OpHitMonitoringAlgorithm class
 */
class OpHitMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    OpHitMonitoringAlgorithm();

    virtual ~OpHitMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // LAR_OP_HITMONITORING_ALGORITHM_H
