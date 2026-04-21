/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/TutorialClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/TutorialClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TutorialClusteringAlgorithm::TutorialClusteringAlgorithm() :
    m_caloHitListName{""}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TutorialClusteringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    // Read the list of calo hits from the XML configured list name into our local CaloHitList
    // If the read fails, abort the event
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    // If there are no hits, return and continue with the event (maybe we only have hits in two views)
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitList, "Hits", GRAY));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    // Create the temporary list to hold new clusters
    const ClusterList *pTemporaryList(nullptr);
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTemporaryList, temporaryListName));

    // Loop over the hits and add to the current cluster, creating a new cluster at random intervals
    const Cluster *pCluster(nullptr);
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;
        if (!pCluster || (rand() % 10) == 0)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.emplace_back(pCaloHit);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit));
        }
    }

    // Save the list of clusters
    if (!pTemporaryList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputClusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TutorialClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read the string value from XML tag CaloHitListName into the member variable m_caloHitListName
    // If the tag is missing, the PANDORA_RETURN_RESULT_IF macro will abort the event 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
