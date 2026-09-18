/**
 *  @file   larpandoracontent/LArPersistency/EventWritingAlgorithm.h
 *
 *  @brief  Header file for the event writing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_WRITING_ALGORITHM_H
#define LAR_EVENT_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Persistency/PandoraIO.h"

namespace pandora
{
class FileWriter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  EventWritingAlgorithm class
 */
class EventWritingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventWritingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventWritingAlgorithm();

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    /**
     *  @brief  Whether current event passes nuance code filter
     *
     *  @return boolean
     */
    bool PassNuanceCodeFilter() const;

    /**
     *  @brief  Whether current event passes mc particle constituent filter
     *
     *  @return boolean
     */
    bool PassMCParticleFilter() const;

    /**
     *  @brief  Whether current event passes neutrino vertex position filter (e.g. fiducial volume cut)
     *
     *  @return boolean
     */
    bool PassNeutrinoVertexFilter() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::FileType m_geometryFileType; ///< The geometry file type
    pandora::FileType m_eventFileType;    ///< The event file type

    pandora::FileWriter *m_pEventFileWriter;    ///< Address of the event file writer
    pandora::FileWriter *m_pGeometryFileWriter; ///< Address of the geometry file writer

    bool m_shouldWriteGeometry;     ///< Whether to write geometry to a specified file
    bool m_writtenGeometry;         ///< Whether geometry has been written
    std::string m_geometryFileName; ///< Name of the output geometry file

    bool m_shouldWriteEvents;    ///< Whether to write events to a specified file
    std::string m_eventFileName; ///< Name of the output event file

    bool m_writtenEventGlobalHeader; ///< Whether the global header has been written to the output event file

    bool m_shouldWriteMCRelationships;    ///< Whether to write mc relationship information to the events file
    bool m_shouldWriteTrackRelationships; ///< Whether to write track relationship information to the events file

    bool m_shouldOverwriteEventFile;    ///< Whether to overwrite existing event file with specified name, or append
    bool m_shouldOverwriteGeometryFile; ///< Whether to overwrite existing geometry file with specified name, or append

    bool m_useLArCaloHits;    ///< Whether to write lar calo hits, or standard pandora calo hits
    bool m_useLArMCParticles; ///< Whether to write lar mc particles, or standard pandora mc particles

    bool m_shouldFilterByNuanceCode; ///< Whether to filter output by nuance code
    int m_filterNuanceCode;          ///< The filter nuance code (required if specify filter by nuance code)

    bool m_shouldFilterByMCParticles;      ///< Whether to filter output by mc particle constituents
    bool m_neutrinoInducedOnly;            ///< Whether to consider only mc particles that were neutrino induced
    unsigned int m_matchingMinPrimaryHits; ///< The minimum number of mc primary hits used in matching scheme
    unsigned int m_nNonNeutrons;           ///< The requested number of mc primaries that are not neutrons
    unsigned int m_nMuons;                 ///< The requested number of mc primaries that are muons
    unsigned int m_nElectrons;             ///< The requested number of mc primaries that are electrons
    unsigned int m_nProtons;               ///< The requested number of mc primaries that are protons
    unsigned int m_nPhotons;               ///< The requested number of mc primaries that are photons
    unsigned int m_nChargedPions;          ///< The requested number of mc primaries that are charged pions

    bool m_shouldFilterByNeutrinoVertex; ///< Whether to filter output by neutrino vertex position (e.g. fiducial volume cut)
    float m_detectorHalfLengthX;         ///< Half length of detector in x dimension
    float m_detectorHalfLengthY;         ///< Half length of detector in y dimension
    float m_detectorHalfLengthZ;         ///< Half length of detector in z dimension
    float m_coordinateOffsetX;           ///< Origin offset (from detector corner) in x dimension
    float m_coordinateOffsetY;           ///< Origin offset (from detector corner) in y dimension
    float m_coordinateOffsetZ;           ///< Origin offset (from detector corner) in z dimension
    float m_selectedBorderX;             ///< Required distance from detector edge in x dimension
    float m_selectedBorderY;             ///< Required distance from detector edge in y dimension
    float m_selectedBorderZ;             ///< Required distance from detector edge in z dimension
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_WRITING_ALGORITHM_H
