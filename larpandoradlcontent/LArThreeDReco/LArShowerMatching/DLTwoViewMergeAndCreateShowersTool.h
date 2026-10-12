/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLTwoViewMergeAndCreateShowersTool.h
 *
 *  @brief  Header file for the two view merge and create showers tool class. 
 *
 *  $Log: $
 */
#ifndef DL_TWO_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H
#define DL_TWO_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLTwoViewMergeAndCreateShowersTool class
 */
class DLTwoViewMergeAndCreateShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLTwoViewMergeAndCreateShowersTool();

    /**
     *  @brief  Create shower-like pfos from ambiguous L:M:N matched cluster groups,
     *          where exactly one of L,M,N is zero and the other two are >=1
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
     *         where exactly one of L,M,N is zero and the other two are >=1
     *
     *  @param pAlgorithm the calling algorithm
     *  @param clusterGroup the input connected cluster group
     *  @param globalSimMatrix the input cluster->cluster->score mapping
     *  @param hitType1 the hit type of one of the two non-empty views
     *  @param hitType2 the other hit type of the two non-empty views
     *
     *  @return whether a particle was created
     */ 
    bool CreateShower(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix,
        const pandora::HitType hitType1, const pandora::HitType hitType2);

    /**
     *  @brief Within the ambiguous cluster group, pick out the 'best' cluster pair
     *
     *  @param clusterGroup the input connected cluster group
     *  @param globalSimMatrix the input cluster->cluster->score mapping
     *  @param hitType1 the hit type of one of the two non-empty views
     *  @param hitType2 the other hit type of the two non-empty views
     *  @param[out] the seed cluster in the first view
     *  @param[out] the seed cluster in the second view    
     *
     *  @return whether a seed cluster pair was found
     */      
    bool FindSeed(const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::HitType hitType1,
        const pandora::HitType hitType2, const pandora::Cluster *&pSeed1, const pandora::Cluster *&pSeed2);

    /**
     *  @brief Grow the seed clusters by merging 'similar' same-view clusters within the connected group
     *
     *  @param pAlgorithm the calling algorithm
     *  @param clusterGroup the input connected cluster group
     *  @param globalSimMatrix the input cluster->cluster->score mapping
     *  @param the seed cluster in the first view
     *  @param the seed cluster in the second view     
     */ 
    void MergeClusters(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::Cluster *const pSeed1,
        const pandora::Cluster *const pSeed2);

    std::string m_clusterListNameU; ///< The input cluster list name for the U view
    std::string m_clusterListNameV; ///< The input cluster list name for the V view
    std::string m_clusterListNameW; ///< The input cluster list name for the W view      
    float m_matchThreshold; ///< Threshold similarity score required for match
};

} // namespace lar_dl_content

#endif // #ifndef DL_TWO_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H

