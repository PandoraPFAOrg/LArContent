/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLTwoViewClearShowersTool.h
 *
 *  @brief  Header file for the two view clear showers tool class.
 *
 *  $Log: $
 */
#ifndef DL_TWO_VIEW_CLEAR_SHOWERS_TOOL_H
#define DL_TWO_VIEW_CLEAR_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLTwoViewClearShowersTool class
 */
class DLTwoViewClearShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLTwoViewClearShowersTool();

    /**
     *  @brief  Create shower-like pfos from unambiguous 0:1:1/1:0:1/1:1:0 matched cluster groups
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
     *  @brief Create a shower-like particle from a 0:1:1/1:0:1/1:1:0 matched cluster group
     *
     *  @param pAlgorithm the calling algorithm
     *  @param clusterGroup the input connected cluster group
     */  
    void CreateShowers(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup);
};

} // namespace lar_dl_content

#endif // #ifndef DL_TWO_VIEW_CLEAR_SHOWERS_TOOL_H

