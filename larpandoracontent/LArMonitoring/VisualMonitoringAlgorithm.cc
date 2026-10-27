/**
 *  @file   larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the visual monitoring algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

VisualMonitoringAlgorithm::VisualMonitoringAlgorithm() :
    m_showCurrentMCParticles(false),
    m_showCurrentCaloHits(false),
    m_showCurrentTracks(false),
    m_showCurrentClusters(false),
    m_showCurrentPfos(false),
    m_showCurrentVertices(false),
    m_displayEvent(true),
    m_saveEventPath(""),
    m_showDetector(false),
    m_detectorView("xz"),
    m_showOnlyAvailable(false),
    m_showAssociatedTracks(false),
    m_hitColors("iterate"),
    m_thresholdEnergy(-1.f),
    m_transparencyThresholdE(-1.f),
    m_energyScaleThresholdE(1.f),
    m_scalingFactor(1.f),
    m_showPfoVertices(true),
    m_showPfoHierarchy(true),
    m_partitionCaloHitsByVolume(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), m_showDetector,
        (m_detectorView.find("xz") != std::string::npos)       ? DETECTOR_VIEW_XZ
            : (m_detectorView.find("xy") != std::string::npos) ? DETECTOR_VIEW_XY
                                                               : DETECTOR_VIEW_DEFAULT,
        m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));

    // Show current mc particles
    if (m_showCurrentMCParticles)
    {
        this->VisualizeMCParticleList(std::string());
    }

    // Show specified lists of mc particles
    for (StringVector::const_iterator iter = m_mcParticleListNames.begin(), iterEnd = m_mcParticleListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeMCParticleList(*iter);
    }

    // Show current calo hit list
    if (m_showCurrentCaloHits)
    {
        this->VisualizeCaloHitList(std::string());
    }

    // Show specified lists of calo hits
    for (StringVector::const_iterator iter = m_caloHitListNames.begin(), iterEnd = m_caloHitListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeCaloHitList(*iter);
    }

    // Show specified lists of op hits
    for (StringVector::const_iterator iter = m_opHitListNames.begin(), iterEnd = m_opHitListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeOpHitList(*iter);
    }

    // Show current cluster list
    if (m_showCurrentClusters)
    {
        this->VisualizeClusterList(std::string());
    }

    // Show specified lists of clusters
    for (StringVector::const_iterator iter = m_clusterListNames.begin(), iterEnd = m_clusterListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeClusterList(*iter);
    }

    // Show current track list
    if (m_showCurrentTracks)
    {
        this->VisualizeTrackList(std::string());
    }

    // Show specified lists of tracks
    for (StringVector::const_iterator iter = m_trackListNames.begin(), iterEnd = m_trackListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeTrackList(*iter);
    }

    // Show current particle flow objects
    if (m_showCurrentPfos)
    {
        this->VisualizeParticleFlowList(std::string());
    }

    // Show specified lists of pfo
    for (StringVector::const_iterator iter = m_pfoListNames.begin(), iterEnd = m_pfoListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeParticleFlowList(*iter);
    }

    // Show current vertex objects
    if (m_showCurrentVertices)
    {
        this->VisualizeVertexList(std::string());
    }

    // Show specified lists of vertices
    for (StringVector::const_iterator iter = m_vertexListNames.begin(), iterEnd = m_vertexListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeVertexList(*iter);
    }

    // Finally, display the event and pause application
    // Save the event displays of the event if given a path to save to.
    if (m_displayEvent && m_saveEventPath != "")
    {
        PANDORA_MONITORING_API(SaveAndViewEvent(this->GetPandora(), m_saveEventPath));
    }
    else if (m_displayEvent)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeMCParticleList(const std::string &listName) const
{
    const MCParticleList *pMCParticleList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pMCParticleList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: mc particle list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pMCParticleList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: mc particle list unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeMCParticles(
        this->GetPandora(), pMCParticleList, listName.empty() ? "CurrentMCParticles" : listName.c_str(), AUTO, &m_particleSuppressionMap));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeCaloHitList(const std::string &listName) const
{
    const CaloHitList *pCaloHitList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pCaloHitList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: current calo hit list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pCaloHitList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: calo hit list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    // Filter calo hit list
    CaloHitList caloHitList;

    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pCaloHit = *iter;

        if ((pCaloHit->GetElectromagneticEnergy() > m_thresholdEnergy) && (!m_showOnlyAvailable || PandoraContentApi::IsAvailable(*this, pCaloHit)))
        {
            caloHitList.push_back(pCaloHit);
        }
    }

    if (!m_partitionCaloHitsByVolume)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, listName.empty() ? "CurrentCaloHits" : listName.c_str(),
            (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY : GRAY)));
        return;
    }

    std::map<std::pair<unsigned int, unsigned int>, CaloHitList> volIdToCaloHits;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        PANDORA_THROW_IF(STATUS_CODE_FAILURE, !pLArCaloHit);
        const std::pair<unsigned int, unsigned int> volId{pLArCaloHit->GetLArTPCVolumeId(), pLArCaloHit->GetDaughterVolumeId()};
        volIdToCaloHits[volId].push_back(pCaloHit);
    }

#ifdef MONITORING
    Color colorIter{BLACK};
    for (const auto &[volId, caloHitListPartition] : volIdToCaloHits)
    {
        const std::string label{(listName.empty() ? "CurrentCaloHits" : listName) + ", volId (" + std::to_string(volId.first) + ", " +
            std::to_string(volId.second) + ")"};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListPartition, label.c_str(), colorIter));
        colorIter = static_cast<Color>(static_cast<int>(colorIter) + 1);
        if (colorIter >= AUTO)
            colorIter = BLACK;
    }
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeOpHitList(const std::string &listName) const
{
    const CaloHitList *pOpHitList{nullptr};

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pOpHitList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VisualMonitoringAlgorithm: optical hit list " << listName << " unavailable." << std::endl;
        return;
    }

    // Filter calo hit list
    CaloHitList opHitList;
    for (const CaloHit *const pOpHit : *pOpHitList)
        opHitList.push_back(pOpHit);

    PANDORA_MONITORING_API(VisualizeOpHits(this->GetPandora(), &opHitList, listName, -1.f, OPHIT_COLOR_SCALE_LOGARITHMIC));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeTrackList(const std::string &listName) const
{
    const TrackList *pTrackList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pTrackList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: current track list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pTrackList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: track list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    // Filter track list
    TrackList trackList;

    for (TrackList::const_iterator iter = pTrackList->begin(), iterEnd = pTrackList->end(); iter != iterEnd; ++iter)
    {
        const Track *const pTrack = *iter;

        if (!m_showOnlyAvailable || pTrack->IsAvailable())
            trackList.push_back(pTrack);
    }

    PANDORA_MONITORING_API(VisualizeTracks(this->GetPandora(), &trackList, listName.empty() ? "CurrentTracks" : listName.c_str(), GRAY));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeClusterList(const std::string &listName) const
{
    const ClusterList *pClusterList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pClusterList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: current cluster list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pClusterList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: cluster list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    // Filter cluster list
    ClusterList clusterList;

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!m_showOnlyAvailable || PandoraContentApi::IsAvailable(*this, pCluster))
            clusterList.push_back(pCluster);
    }

    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, listName.empty() ? "CurrentClusters" : listName.c_str(),
        (m_hitColors.find("particleid") != std::string::npos)    ? AUTOID
            : (m_hitColors.find("iterate") != std::string::npos) ? AUTOITER
            : (m_hitColors.find("energy") != std::string::npos)  ? AUTOENERGY
                                                                 : AUTO,
        m_showAssociatedTracks));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeParticleFlowList(const std::string &listName) const
{
    const PfoList *pPfoList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: current pfo list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: pfo list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, listName.empty() ? "CurrentPfos" : listName.c_str(),
        (m_hitColors.find("particleid") != std::string::npos)
            ? AUTOID
            : (m_hitColors.find("iterate") != std::string::npos ? AUTOITER : (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY : AUTO)),
        m_showPfoVertices, m_showPfoHierarchy));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeVertexList(const std::string &listName) const
{
    const VertexList *pVertexList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pVertexList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: current vertex list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pVertexList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VisualMonitoringAlgorithm: vertex list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    // Filter vertex list
    VertexList vertexList;

    for (VertexList::const_iterator iter = pVertexList->begin(), iterEnd = pVertexList->end(); iter != iterEnd; ++iter)
    {
        const Vertex *const pVertex = *iter;

        if (!m_showOnlyAvailable || pVertex->IsAvailable())
            vertexList.push_back(pVertex);
    }

    PANDORA_MONITORING_API(VisualizeVertices(this->GetPandora(), &vertexList, listName.empty() ? "CurrentVertices" : listName.c_str(), AUTO));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowCurrentMCParticles", m_showCurrentMCParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "MCParticleListNames", m_mcParticleListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowCurrentCaloHits", m_showCurrentCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "OpHitListNames", m_opHitListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowCurrentTracks", m_showCurrentTracks));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "TrackListNames", m_trackListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowCurrentClusters", m_showCurrentClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowCurrentPfos", m_showCurrentPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowCurrentVertices", m_showCurrentVertices));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "VertexListNames", m_vertexListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DisplayEvent", m_displayEvent));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SaveEventPath", m_saveEventPath));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowDetector", m_showDetector));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DetectorView", m_detectorView));
    std::transform(m_detectorView.begin(), m_detectorView.end(), m_detectorView.begin(), ::tolower);

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowOnlyAvailable", m_showOnlyAvailable));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowAssociatedTracks", m_showAssociatedTracks));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitColors", m_hitColors));
    std::transform(m_hitColors.begin(), m_hitColors.end(), m_hitColors.begin(), ::tolower);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ThresholdEnergy", m_thresholdEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowPfoVertices", m_showPfoVertices));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowPfoHierarchy", m_showPfoHierarchy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "PartitionCaloHitsByVolume", m_partitionCaloHitsByVolume));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "SuppressMCParticles", m_suppressMCParticles));

    for (StringVector::iterator iter = m_suppressMCParticles.begin(), iterEnd = m_suppressMCParticles.end(); iter != iterEnd; ++iter)
    {
        const std::string pdgEnergy(*iter);
        StringVector pdgEnergySeparated;
        const std::string delimiter = ":";
        XmlHelper::TokenizeString(pdgEnergy, pdgEnergySeparated, delimiter);

        if (pdgEnergySeparated.size() != 2)
            return STATUS_CODE_INVALID_PARAMETER;

        int pdgCode(0);
        float energy(0.f);

        if (!StringToType(pdgEnergySeparated.at(0), pdgCode) || !StringToType(pdgEnergySeparated.at(1), energy))
            return STATUS_CODE_INVALID_PARAMETER;

        m_particleSuppressionMap.insert(PdgCodeToEnergyMap::value_type(pdgCode, energy));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
