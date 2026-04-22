/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TutorialClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TutorialClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TutorialClusterMergingAlgorithm::TutorialClusterMergingAlgorithm() :
    m_inputClusterListName{""},
    m_outputClusterListName{""}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TutorialClusterMergingAlgorithm::Run()
{
    const ClusterList *pInputClusterList{nullptr};
    // Read the list of clusters from the XML configured list name into our local ClusterList
    // If the read fails, abort the event
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));
    // If there are no clusters, return and continue with the event (maybe we only have hits in two views)
    if (!pInputClusterList || pInputClusterList->empty())
        return STATUS_CODE_SUCCESS;
    // You could apply some pre-selection on the input clusters here, the key thing is, for the next stage, we want to work with a list we can
    // modify directly. Here, no selection is applied, so we'll just shallow copy the list
    ClusterVector selectedClusters(pInputClusterList->begin(), pInputClusterList->end());

    // We don't need any temporary list management here, because we're going to be updating existing clusters in an existing list
    ClusterList deletedClusters;

    // Loop over the input clusters and make (random) merge decisions
    for (const Cluster *const pParentCluster : selectedClusters)
    {
        if (std::find(deletedClusters.begin(), deletedClusters.end(), pParentCluster) != deletedClusters.end())
            continue;
        for (const Cluster *const pChildCluster : selectedClusters)
        {
            if (!pChildCluster || pParentCluster == pChildCluster)
                continue;
            if (std::find(deletedClusters.begin(), deletedClusters.end(), pChildCluster) != deletedClusters.end())
                continue;
            if (!this->AreClustersMergeable(pParentCluster, pChildCluster))
                continue;
            deletedClusters.push_back(pChildCluster);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pChildCluster));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TutorialClusterMergingAlgorithm::AreClustersMergeable([[maybe_unused]] const Cluster *const pParentCluster,
    [[maybe_unused]]const Cluster *const pChildCluster) const
{
    // This is where you would implement your merging decision. For this tutorial, we'll just merge clusters at random.
    return (rand() % 10) == 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TutorialClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read the string value from XML tag InputClusterListName into the member variable m_inputClusterListName
    // If the tag is missing, the PANDORA_RETURN_RESULT_IF macro will abort the event 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
