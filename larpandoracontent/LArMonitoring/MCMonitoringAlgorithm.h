/**
 *  @file   larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_MONITORING_ALGORITHM_H
#define LAR_MC_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MCMonitoringAlgorithm class
 */
class MCMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    MCMonitoringAlgorithm();

    virtual ~MCMonitoringAlgorithm();

private:
    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCHitsMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode BuildMCHitMap();
   
    std::string m_caloHitListName; ///< Name of the CaloHit list to be used
    std::string m_mcListName; ///< Name of the MC particle list to be used
    bool m_visualise{false}; ///< Flag to indicate whether to visualise the results
    bool m_colourByProcess{false}; ///< Flag to indicate whether to colour hits by process
    MCHitsMap m_mcHitsMap; ///< Map of MC particles to their associated hits
};

} // namespace lar_content

#endif // LAR_MC_MONITORING_ALGORITHM_H
