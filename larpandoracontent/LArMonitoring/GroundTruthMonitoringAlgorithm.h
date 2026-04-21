/**
 *  @file   larpandoracontent/LArMonitoring/GroundTruthMonitoringAlgorithm.h
 *
 *  @brief  An algorithm to visualise ground truth information.
 *
 *  $Log: $
 */
#ifndef LAR_GROUND_TRUTH_MONITORING_ALGORITHM_H
#define LAR_GROUND_TRUTH_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArUtility/RollUp.h"

namespace lar_content
{

/**
 *  @brief  GroundTruthMonitoringAlgorithm class
 */
class GroundTruthMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    GroundTruthMonitoringAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~GroundTruthMonitoringAlgorithm() = default;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Construct a map of MC particles to their hits
     *
     *  @param  caloHitList the list of hits to consider
     *  @param  mcMap the map to be filled with the MC particle to hit associations
     */
    void MakeMCToHitsMap(const pandora::CaloHitList &caloHitList, LArMCParticleHelper::MCContributionMap &mcMap) const;

    /**
     *  @brief  Partition the hits into their respective views
     *
     *  @param  hits2D the list of hits to partition
     *  @param  uHits the list to be filled with hits in the U view
     *  @param  vHits the list to be filled with hits in the V view
     *  @param  wHits the list to be filled with hits in the W view
     */
    void PartitionViews(const pandora::CaloHitList &hits2D, pandora::CaloHitList &uHits, pandora::CaloHitList &vHits, pandora::CaloHitList &wHits) const;

    /**
     *  @brief  Visualise the hits associated to each MC particle
     *
     *  @param  mcMap the map of MC particles to their associated hits
     */
    void Visualize(const LArMCParticleHelper::MCContributionMap &mcMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< Name of input list containing all 2D hits
    RollUpper m_rollUp;             ///< Roll-up object to determine the main contributing MC particle for each hit
};

} // namespace lar_content

#endif // LAR_GROUND_TRUTH_MONITORING_ALGORITHM_H
