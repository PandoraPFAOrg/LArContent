/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLThreeViewMergeAndCreateShowersTool.h
 *
 *  @brief  Header file for the three view merge and create showers tool class.
 *
 *  $Log: $
 */
#ifndef DL_THREE_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H
#define DL_THREE_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLThreeViewMergeAndCreateShowersTool class
 */
class DLThreeViewMergeAndCreateShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLThreeViewMergeAndCreateShowersTool();

    /**
     *  @brief  Create shower-like pfos from ambiguous L:M:N matched cluster groups
     *
     *  @param pAlgorithm address of the calling algorithm
     *  @param globalSimMatrix the output cluster->cluster->score mapping
     *
     *  @return whether shower-like pfos have been created
     */
    bool Run(DLMultiViewMatchingAlgorithm *const pAlgorithm,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create a shower-like particle from a L:M:N matched cluster group,
     *         where L,M,N >= 1 and at least two of L,M,N are >1
     *
     *  @param pAlgorithm the calling algorithm
     *  @param clusterGroup the input connected cluster group
     *  @param globalSimMatrix the input cluster->cluster->score mapping
     *
     *  @return whether a particle was created
     */  
    bool CreateShower(DLMultiViewMatchingAlgorithm *const pAlgorithm,
        const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

    /**
     *  @brief Within the ambiguous cluster group, pick out the 'best' cluster triplet
     *
     *  @param clusterGroup the input connected cluster group
     *  @param globalSimMatrix the input cluster->cluster->score mapping
     *  @param[out] the U-view seed cluster
     *  @param[out] the V-view seed cluster
     *  @param[out] the W-view seed cluster      
     *
     *  @return whether a seed cluster triplet was found i.e. in which all pairwise cluster pairs are 'similar'
     */      
    bool FindSeed(const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::Cluster *&pSeedU,
        const pandora::Cluster *&pSeedV, const pandora::Cluster *&pSeedW);

    /**
     *  @brief Grow the seed clusters by merging 'similar' same-view clusters within the connected group
     *
     *  @param pAlgorithm the calling algorithm
     *  @param clusterGroup the input connected cluster group
     *  @param globalSimMatrix the input cluster->cluster->score mapping
     *  @param the U-view seed cluster
     *  @param the V-view seed cluster
     *  @param the W-view seed cluster      
     */ 
    void MergeClusters(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::Cluster *const pSeedU,
        const pandora::Cluster *const pSeedV, const pandora::Cluster *const pSeedW);

    std::string m_clusterListNameU; ///< The input cluster list name for the U view
    std::string m_clusterListNameV; ///< The input cluster list name for the V view
    std::string m_clusterListNameW; ///< The input cluster list name for the W view
    float m_matchThreshold; ///< Threshold similarity score required for match
};

} // namespace lar_dl_content

#endif // #ifndef DL_THREE_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H

