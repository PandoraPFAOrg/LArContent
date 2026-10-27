/**
 *  @file   larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.h
 *
 *  @brief  Header file for the visual monitoring algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_VISUAL_MONITORING_ALGORITHM_H
#define LAR_VISUAL_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief VisualMonitoringAlgorithm class
 */
class VisualMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VisualMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualize mc particle list
     *
     *  @param  listName the list name
     */
    void VisualizeMCParticleList(const std::string &listName) const;

    /**
     *  @brief  Visualize a specified calo hit list
     *
     *  @param  listName the list name
     */
    void VisualizeCaloHitList(const std::string &listName) const;

    /**
     *  @brief  Visualize a specified op hit list
     *
     *  @param  listName the list name
     */
    void VisualizeOpHitList(const std::string &listName) const;

    /**
     *  @brief  Visualize a specified track list
     *
     *  @param  listName the list name
     */
    void VisualizeTrackList(const std::string &listName) const;

    /**
     *  @brief  Visualize a specified calo hit list
     *
     *  @param  listName the list name
     */
    void VisualizeClusterList(const std::string &listName) const;

    /**
     *  @brief  Visualize a specified pfo list
     *
     *  @param  listName the list name
     */
    void VisualizeParticleFlowList(const std::string &listName) const;

    /**
     *  @brief  Visualize a specified vertex list
     *
     *  @param  listName the list name
     */
    void VisualizeVertexList(const std::string &listName) const;

    typedef std::map<int, float> PdgCodeToEnergyMap;

    bool m_showCurrentMCParticles;               ///< Whether to show current mc particles
    pandora::StringVector m_mcParticleListNames; ///< Names of mc particles lists to show

    bool m_showCurrentCaloHits;               ///< Whether to show current calohitlist
    pandora::StringVector m_caloHitListNames; ///< Names of calo hit lists to show

    pandora::StringVector m_opHitListNames;   ///< Names of optical hit lists to show

    bool m_showCurrentTracks;               ///< Whether to show current tracks
    pandora::StringVector m_trackListNames; ///< Names of track lists to show

    bool m_showCurrentClusters;               ///< Whether to show current clusters
    pandora::StringVector m_clusterListNames; ///< Names of cluster lists to show

    bool m_showCurrentPfos;               ///< Whether to show current particle flow object list
    pandora::StringVector m_pfoListNames; ///< Names of pfo lists to show

    bool m_showCurrentVertices;              ///< Whether to show current vertex list
    pandora::StringVector m_vertexListNames; ///< Names of vertex lists to show

    bool m_displayEvent;         ///< Whether to display the event
    std::string m_saveEventPath; ///< The path to save event displays to. m_displayEvent must also be set.
    bool m_showDetector;         ///< Whether to display the detector geometry
    std::string m_detectorView;  ///< The detector view, default, xy or xz

    bool m_showOnlyAvailable;       ///< Whether to show only available  (i.e. non-clustered) calohits and tracks
    bool m_showAssociatedTracks;    ///< Whether to display tracks associated to clusters when viewing cluster lists
    std::string m_hitColors;        ///< Define the hit coloring scheme (default: iterate, choices: iterate, energy, particleid)
    float m_thresholdEnergy;        ///< Cell energy threshold for display (em scale)
    float m_transparencyThresholdE; ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float m_energyScaleThresholdE;  ///< Cell energy for which color is at top end of continous color palette
    float m_scalingFactor;          ///< TEve works with [cm], Pandora usually works with [mm] (but LArContent went with cm too)

    bool m_showPfoVertices;           ///< Whether to display pfo vertices
    bool m_showPfoHierarchy;          ///< Whether to display daughter pfos only under parent pfo elements
    bool m_partitionCaloHitsByVolume; ///< Whether to partition calo hits into unique TPC + daughter volumes for display

    pandora::StringVector m_suppressMCParticles; ///< List of PDG numbers and energies for MC particles to be suppressed (e.g. " 22:0.1 2112:1.0 ")
    PdgCodeToEnergyMap m_particleSuppressionMap; ///< Map from pdg-codes to energy for suppression of particles types below specific energies
};

} // namespace lar_content

#endif // #ifndef LAR_VISUAL_MONITORING_ALGORITHM_H
