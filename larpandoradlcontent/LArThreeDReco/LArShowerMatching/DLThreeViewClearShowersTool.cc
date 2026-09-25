/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLThreeViewClearShowersTool.cc
 *
 *  @brief  Implementation of the three view clear showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLThreeViewClearShowersTool.h"

using namespace pandora;

namespace lar_dl_content
{

DLThreeViewClearShowersTool::DLThreeViewClearShowersTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLThreeViewClearShowersTool::Run(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    DLMultiViewMatchingAlgorithm::ClusterGroupVector clusterGroupVector;
    pAlgorithm->GetConnectedGroups(clusterGroupVector);    

    bool madeParticles(false);
    for (const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup : clusterGroupVector)
    {
        int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());                
        
        // Look for 1:1:1 unambiguous match
        if ((nU * nV * nW) != 1)
            continue;
            
        this->CreateShowers(pAlgorithm, clusterGroup);
        madeParticles = true;
    }

    return madeParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeViewClearShowersTool::CreateShowers(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup)
{
    pAlgorithm->CreatePfo({clusterGroup.m_clustersU.front(), clusterGroup.m_clustersV.front(), clusterGroup.m_clustersW.front()});
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeViewClearShowersTool::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content

