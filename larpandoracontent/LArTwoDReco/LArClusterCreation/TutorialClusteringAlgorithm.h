/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/TutorialClusteringAlgorithm.h
 *
 *  @brief  A simple clustering algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_TUTORIAL_CLUSTERING_ALGORITHM_H
#define LAR_TUTORIAL_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TutorialClusteringAlgorithm class
 */
class TutorialClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TutorialClusteringAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< Name of input list containing the 2D hits to cluster
    std::string m_outputClusterListName;  ///< Name of the list to write the output clusters to
};

} // namespace lar_content

#endif // #ifndef LAR_TUTORIAL_CLUSTERING_ALGORITHM_H
