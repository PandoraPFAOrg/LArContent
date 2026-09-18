/**
 *  @file   larpandoracontent/LArPersistency/EventReadingAlgorithm.h
 *
 *  @brief  Header file for the event reading algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_READING_ALGORITHM_H
#define LAR_EVENT_READING_ALGORITHM_H 1

#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "Pandora/PandoraInputTypes.h"

#include "Persistency/PandoraIO.h"

namespace pandora
{
class FileReader;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  EventReadingAlgorithm class
 */
class EventReadingAlgorithm : public pandora::ExternallyConfiguredAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventReadingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventReadingAlgorithm();

    /**
     *  @brief  External event reading parameters class
     */
    class ExternalEventReadingParameters : public pandora::ExternalParameters
    {
    public:
        std::string m_geometryFileName;   ///< Name of the file containing geometry information
        std::string m_eventFileNameList;  ///< Colon-separated list of file names to be processed
        pandora::InputUInt m_skipToEvent; ///< Index of first event to consider in input file
    };

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    /**
     *  @brief  Proceed to process next event file named in the input list
     */
    void MoveToNextEventFile();

    /**
     *  @brief  Replace the current event file reader with a new reader for the specified file
     *
     *  @param  fileName the file name
     */
    pandora::StatusCode ReplaceEventFileReader(const std::string &fileName);

    /**
     *  @brief  Analyze a provided file name to extract the file type/extension
     *
     *  @param  fileName the file name
     *
     *  @return the file type
     */
    pandora::FileType GetFileType(const std::string &fileName) const;

    /**
     *  @brief  Install the schema migrations this build knows about on a reader, before its global header is read.
     *
     *  @param  fileReader the reader to configure
     */
    void ConfigureReader(pandora::FileReader &fileReader) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_geometryFileName;              ///< Name of the file containing geometry information
    std::string m_eventFileName;                 ///< Name of the current file containing event information
    pandora::StringVector m_eventFileNameVector; ///< Vector of file names to be processed

    unsigned int m_skipToEvent;  ///< Index of first event to consider in first input file
    bool m_useLArCaloHits;       ///< Whether to read lar calo hits, or standard pandora calo hits
    bool m_useLArMCParticles;    ///< Whether to read lar mc particles, or standard pandora mc particles
    bool m_strictSchemaChecking; ///< Whether an unreadable schema version is an error rather than a warning

    pandora::FileReader *m_pEventFileReader; ///< Address of the event file reader
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_READING_ALGORITHM_H
