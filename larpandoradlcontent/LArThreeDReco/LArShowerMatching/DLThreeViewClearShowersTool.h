/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLThreeViewClearShowersTool.h
 *
 *  @brief  Header file for the three view clear showers tool class.
 *
 *  $Log: $
 */
#ifndef DL_THREE_VIEW_CLEAR_SHOWERS_TOOL_H
#define DL_THREE_VIEW_CLEAR_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLThreeViewClearShowersTool class
 */
class DLThreeViewClearShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLThreeViewClearShowersTool();

    /**
     *  @brief  Create shower-like pfos from unambiguous 1:1:1 matched cluster groups
     *
     *  @param pAlgorithm address of the calling algorithm
     *  @param globalSimMatrix the output cluster->cluster->score mapping
     *
     *  @return whether shower-like pfos have been created
     */    
    bool Run(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create a shower-like particle from a 1:1:1 matched cluster group
     *
     *  @param pAlgorithm the calling algorithm
     *  @param clusterGroup the input connected cluster group
     */     
    void CreateShowers(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup);
};

} // namespace lar_dl_content

#endif // #ifndef DL_THREE_VIEW_CLEAR_SHOWERS_TOOL_H

