/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/TutorialClusterMergingAlgorithm.h
 *
 *  @brief  A simple clustering algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_TUTORIAL_CLUSTER_MERGING_ALGORITHM_H
#define LAR_TUTORIAL_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TutorialClusterMergingAlgorithm class
 */
class TutorialClusterMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TutorialClusterMergingAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Determine if two clusters can be merged.
     *
     *  @param  pParentCluster the parent cluster that will be expanded if the clusters can be merged
     *  @param  pChildCluster the child cluster that will be merged into the parent cluster if the clusters can be merged
     *
     *  @return boolean indicating whether the clusters can be merged
     */
    bool AreClustersMergeable(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pChildCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputClusterListName;  ///< Name of input list containing the input 2D clusters
    std::string m_outputClusterListName;  ///< Name of the list to write the output clusters to
};

} // namespace lar_content

#endif // #ifndef LAR_TUTORIAL_CLUSTER_MERGING_ALGORITHM_H
