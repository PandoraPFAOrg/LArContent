/**
 *  @file   larpandoracontent/LArObjects/LArCaloHit.h
 *
 *  @brief  Header file for the lar calo hit class.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_HIT_H
#define LAR_CALO_HIT_H 1

#include "Objects/CaloHit.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/PandoraObjectFactories.h"

#include "Persistency/BinaryFileReader.h"
#include "Persistency/BinaryFileWriter.h"
#include "Persistency/XmlFileReader.h"
#include "Persistency/XmlFileWriter.h"

namespace lar_content
{

/**
 *  @brief  LAr calo hit parameters
 */
class LArHitParameters : public object_creation::CaloHit::Parameters
{
public:
    pandora::InputUInt m_larTPCVolumeId;    ///< The lar tpc volume id
    pandora::InputUInt m_daughterVolumeId;  ///< The daughter volume id
    pandora::InputFloat m_startTime;        ///< The start time of the optical hit
    pandora::InputFloat m_width;            ///< The width of the optical hit
    pandora::InputUInt  m_channel;          ///< The optical detector channel - ATTN: maybe daughter volume can be used instead
    pandora::FloatVector m_hitScores;       ///< Hit scores
    pandora::StringVector m_hitScoreLabels; ///< Labels for the hit scores
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LAr calo hit class
 */
class LArCaloHit : public object_creation::CaloHit::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar calo hit parameters
     */
    LArCaloHit(const LArHitParameters &parameters);

    /**
     *  @brief  Get the lar tpc volume id
     *
     *  @return the lar tpc volume id
     */
    unsigned int GetLArTPCVolumeId() const;

    /**
     *  @brief  Get the daughter volume id
     *
     *  @return the daughter volume id
     */
    unsigned int GetDaughterVolumeId() const;

    /**
     *  @brief  Fill the parameters associated with this calo hit
     *
     *  @param  parameters the output parameters
     */
    void FillParameters(LArHitParameters &parameters) const;

    /**
     *  @brief  Get scores for each hit
     *
     *  @return a vector of scores for each category
     */
    const pandora::FloatVector &GetHitScores() const;

    /**
     *  @brief  Get scores for each hit
     *
     *  @return a vector of score labels for each category
     */
    const pandora::StringVector &GetHitScoreLabels() const;

    /**
     *  @brief  Get the probability that the hit is track-like
     *
     *  @return the probability the hit is track-like
     */
    float GetTrackProbability() const;

    /**
     *  @brief  Get the probability that the hit is shower-like
     *
     *  @return the probability the hit is shower-like
     */
    float GetShowerProbability() const;

    /**
     *  @brief  Set the probability that the hit is track-like
     *
     *  @param  probability the probability the hit is track-like
     */
    void SetTrackProbability(const float probability);

    /**
     *  @brief  Set the probability that the hit is shower-like
     *
     *  @param  probability the probability the hit is shower-like
     */
    void SetShowerProbability(const float probability);

private:
    unsigned int m_larTPCVolumeId;          ///< The lar tpc volume id
    unsigned int m_daughterVolumeId;        ///< The daughter volume id
    pandora::FloatVector m_hitScores;       ///< Hit scores
    pandora::StringVector m_hitScoreLabels; ///< Labels for the hit scores
    pandora::InputFloat m_pTrack;           ///< The probability that the hit is track-like
    pandora::InputFloat m_pShower;          ///< The probability that the hit is shower-like
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LAr optical hit class
 */
class LArOpHit : public object_creation::CaloHit::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar op hit parameters
     */
    LArOpHit(const LArHitParameters &parameters);

    /**
     *  @brief  Get the start time of the optical hit
     *
     *  @return the start time of the optical hit
     */
    float GetStartTime() const;

    /**
     *  @brief  Get the width of the optical hit
     *
     *  @return the width of the optical hit
     */
    float GetWidth() const;

    /**
     *  @brief  Get the channel of the hit
     *
     *  @return the channel of the hit
     */
    unsigned int GetChannel() const;

    /**
     *  @brief  Fill the parameters associated with this calo hit
     *
     *  @param  parameters the output parameters
     */
    void FillParameters(LArHitParameters &parameters) const;

private:
    float m_startTime;      ///< The start time of the optical hit
    float m_width;          ///< The width of the optical hit
    unsigned int m_channel; ///< The channel of the hit
};


/**
 *  @brief  LArCaloHitFactory responsible for object creation
 */
class LArHitFactory : public pandora::ObjectFactory<object_creation::CaloHit::Parameters, object_creation::CaloHit::Object>
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  version the LAr Hit version
     */
    LArHitFactory(const unsigned int version = 2);

    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  fileReader the file reader, used to extract any additional parameters from file
     */
    pandora::StatusCode Read(Parameters &parameters, pandora::FileReader &fileReader) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fileWriter the file writer
     */
    pandora::StatusCode Write(const Object *const pObject, pandora::FileWriter &fileWriter) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const Parameters &parameters, const Object *&pObject) const;

private:
    unsigned int m_version; ///< The LAr hit version
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHit::LArCaloHit(const LArHitParameters &parameters) :
    object_creation::CaloHit::Object(parameters),
    m_larTPCVolumeId(parameters.m_larTPCVolumeId.Get()),
    m_daughterVolumeId(parameters.m_daughterVolumeId.IsInitialized() ? parameters.m_daughterVolumeId.Get() : 0),
    m_hitScores(parameters.m_hitScores),
    m_hitScoreLabels(parameters.m_hitScoreLabels)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetLArTPCVolumeId() const
{
    return m_larTPCVolumeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetDaughterVolumeId() const
{
    return m_daughterVolumeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::FloatVector &LArCaloHit::GetHitScores() const
{
    return m_hitScores;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::StringVector &LArCaloHit::GetHitScoreLabels() const
{
    return m_hitScoreLabels;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::FillParameters(LArHitParameters &parameters) const
{
    parameters.m_positionVector = this->GetPositionVector();
    parameters.m_expectedDirection = this->GetExpectedDirection();
    parameters.m_cellNormalVector = this->GetCellNormalVector();
    parameters.m_cellGeometry = this->GetCellGeometry();
    parameters.m_cellSize0 = this->GetCellSize0();
    parameters.m_cellSize1 = this->GetCellSize1();
    parameters.m_cellThickness = this->GetCellThickness();
    parameters.m_nCellRadiationLengths = this->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = this->GetNCellInteractionLengths();
    parameters.m_time = this->GetTime();
    parameters.m_inputEnergy = this->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = this->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = this->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = this->GetHadronicEnergy();
    parameters.m_isDigital = this->IsDigital();
    parameters.m_hitType = this->GetHitType();
    parameters.m_hitRegion = this->GetHitRegion();
    parameters.m_layer = this->GetLayer();
    parameters.m_isInOuterSamplingLayer = this->IsInOuterSamplingLayer();
    // ATTN Set the parent address to the original owner of the calo hit
    parameters.m_pParentAddress = static_cast<const void *>(this);
    parameters.m_larTPCVolumeId = this->GetLArTPCVolumeId();
    parameters.m_daughterVolumeId = this->GetDaughterVolumeId();
    parameters.m_hitScores = this->GetHitScores();
    parameters.m_hitScoreLabels = this->GetHitScoreLabels();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetTrackProbability() const
{
    return m_pTrack.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetShowerProbability() const
{
    return m_pShower.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::SetTrackProbability(const float probability)
{
    if (probability >= 0.f)
        m_pTrack = probability;
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::SetShowerProbability(const float probability)
{
    if (probability >= 0.f)
        m_pShower = probability;
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArOpHit::LArOpHit(const LArHitParameters &parameters) :
    object_creation::CaloHit::Object(parameters),
    m_startTime(parameters.m_startTime.Get()),
    m_width(parameters.m_width.Get()),
    m_channel(parameters.m_channel.Get())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArOpHit::GetStartTime() const
{
    return m_startTime;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArOpHit::GetWidth() const
{
    return m_width;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArOpHit::GetChannel() const
{
    return m_channel;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArOpHit::FillParameters(LArHitParameters &parameters) const
{
    parameters.m_positionVector = this->GetPositionVector();
    parameters.m_expectedDirection = this->GetExpectedDirection();
    parameters.m_cellNormalVector = this->GetCellNormalVector();
    parameters.m_cellGeometry = this->GetCellGeometry();
    parameters.m_cellSize0 = this->GetCellSize0();
    parameters.m_cellSize1 = this->GetCellSize1();
    parameters.m_cellThickness = this->GetCellThickness();
    parameters.m_nCellRadiationLengths = this->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = this->GetNCellInteractionLengths();
    parameters.m_time = this->GetTime();
    parameters.m_inputEnergy = this->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = this->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = this->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = this->GetHadronicEnergy();
    parameters.m_isDigital = this->IsDigital();
    parameters.m_hitType = this->GetHitType();
    parameters.m_hitRegion = this->GetHitRegion();
    parameters.m_layer = this->GetLayer();
    parameters.m_isInOuterSamplingLayer = this->IsInOuterSamplingLayer();
    // ATTN Set the parent address to the original owner of the calo hit
    parameters.m_pParentAddress = static_cast<const void *>(this);
    parameters.m_startTime = this->GetStartTime();
    parameters.m_width = this->GetWidth();
    parameters.m_channel = this->GetChannel();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArHitFactory::LArHitFactory(const unsigned int version) :
    m_version(version)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArHitFactory::Parameters *LArHitFactory::NewParameters() const
{
    return new LArHitParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArHitFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArHitParameters &larHitParameters(dynamic_cast<const LArHitParameters &>(parameters));

    switch (parameters.m_hitType.Get())
    {
        case pandora::OPTICAL_SIPM:
        case pandora::OPTICAL_TRAP:
        case pandora::OPTICAL_TPC:
        {
            pObject = new LArOpHit(larHitParameters);
            return pandora::STATUS_CODE_SUCCESS;
        }
        default:
        {
            pObject = new LArCaloHit(larHitParameters);
            return pandora::STATUS_CODE_SUCCESS;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArHitFactory::Read(Parameters &parameters, pandora::FileReader &fileReader) const
{
    switch (parameters.m_hitType.Get())
    {
        case pandora::OPTICAL_SIPM:
        case pandora::OPTICAL_TRAP:
        case pandora::OPTICAL_TPC:
        {
            float startTime(0.f);
            float width(0.f);
            unsigned int channel(0);

            if (pandora::BINARY == fileReader.GetFileType())
            {
                pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader &>(fileReader));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(startTime));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(width));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(channel));
            }
            else if (pandora::XML == fileReader.GetFileType())
            {
                pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader &>(fileReader));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("OpticalStartTime", startTime));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("OpticalWidth", width));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("OpticalChannel", channel));
            }
            else
            {
                return pandora::STATUS_CODE_INVALID_PARAMETER;
            }

            LArHitParameters &larOpHitParameters(dynamic_cast<LArHitParameters &>(parameters));
            larOpHitParameters.m_startTime = startTime;
            larOpHitParameters.m_width = width;
            larOpHitParameters.m_channel = channel;

            return pandora::STATUS_CODE_SUCCESS;
        }
        default:
        {
            unsigned int larTPCVolumeId(std::numeric_limits<unsigned int>::max());
            unsigned int daughterVolumeId(0);
            unsigned int nLabels(0);
            pandora::FloatVector hitScores;
            pandora::StringVector hitScoreLabels;

            if (pandora::BINARY == fileReader.GetFileType())
            {
                pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader &>(fileReader));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(larTPCVolumeId));
                if (m_version > 1)
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(daughterVolumeId));
                if (m_version > 2)
                {
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nLabels));
                    hitScores.reserve(nLabels);
                    hitScoreLabels.reserve(nLabels);
                    for (unsigned int i = 0; i < nLabels; ++i)
                    {
                        float score;
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(score));
                        hitScores.emplace_back(std::move(score));

                        std::string label;
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(label));
                        hitScoreLabels.emplace_back(std::move(label));
                    }
                }
            }
            else if (pandora::XML == fileReader.GetFileType())
            {
                pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader &>(fileReader));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("LArTPCVolumeId", larTPCVolumeId));
                if (m_version > 1)
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("DaughterVolumeId", daughterVolumeId));
                if (m_version > 2)
                {
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NHitLabels", nLabels));
                    hitScores.reserve(nLabels);
                    hitScoreLabels.reserve(nLabels);
                    for (unsigned int i = 0; i < nLabels; ++i)
                    {
                        float score;
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("HitScore", score));
                        hitScores.emplace_back(std::move(score));

                        std::string label;
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("HitScoreLabel", label));
                        hitScoreLabels.emplace_back(std::move(label));
                    }
                }
            }
            else
            {
                return pandora::STATUS_CODE_INVALID_PARAMETER;
            }

            LArHitParameters &larCaloHitParameters(dynamic_cast<LArHitParameters &>(parameters));
            larCaloHitParameters.m_larTPCVolumeId = larTPCVolumeId;
            larCaloHitParameters.m_daughterVolumeId = daughterVolumeId;
            larCaloHitParameters.m_hitScores = std::move(hitScores);
            larCaloHitParameters.m_hitScoreLabels = std::move(hitScoreLabels);

            return pandora::STATUS_CODE_SUCCESS;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArHitFactory::Write(const Object *const pObject, pandora::FileWriter &fileWriter) const
{
    // ATTN: To receive this call-back must have already set file writer mc particle factory to this factory
    const pandora::CaloHit *const pCaloHit(dynamic_cast<const pandora::CaloHit *>(pObject));
    if (!pCaloHit)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    switch (pCaloHit->GetHitType())
    {
        case pandora::OPTICAL_SIPM:
        case pandora::OPTICAL_TRAP:
        case pandora::OPTICAL_TPC:
        {
            const LArOpHit *const pLArOpHit(dynamic_cast<const LArOpHit *>(pObject));

            if (!pLArOpHit)
                return pandora::STATUS_CODE_INVALID_PARAMETER;

            if (pandora::BINARY == fileWriter.GetFileType())
            {
                pandora::BinaryFileWriter &binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter &>(fileWriter));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArOpHit->GetStartTime()));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArOpHit->GetWidth()));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArOpHit->GetChannel()));
            }
            else if (pandora::XML == fileWriter.GetFileType())
            {
                pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter &>(fileWriter));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("OpticalStartTime", pLArOpHit->GetWidth()));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("OpticalWidth", pLArOpHit->GetWidth()));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("OpticalChannel", pLArOpHit->GetChannel()));
            }
            else
            {
                return pandora::STATUS_CODE_INVALID_PARAMETER;
            }

            return pandora::STATUS_CODE_SUCCESS;
        }
        default:
        {
            const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pObject));

            if (!pLArCaloHit)
                return pandora::STATUS_CODE_INVALID_PARAMETER;

            if (pandora::BINARY == fileWriter.GetFileType())
            {
                pandora::BinaryFileWriter &binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter &>(fileWriter));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetLArTPCVolumeId()));
                if (m_version > 1)
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetDaughterVolumeId()));
                if (m_version > 2)
                {
                    const pandora::FloatVector &hitScores(pLArCaloHit->GetHitScores());
                    const pandora::StringVector &hitScoreLabels(pLArCaloHit->GetHitScoreLabels());
                    const unsigned int nLabels(hitScores.size());
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(nLabels));
                    for (unsigned int i = 0; i < nLabels; ++i)
                    {
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(hitScores.at(i)));
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(hitScoreLabels.at(i)));
                    }
                }
            }
            else if (pandora::XML == fileWriter.GetFileType())
            {
                pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter &>(fileWriter));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("LArTPCVolumeId", pLArCaloHit->GetLArTPCVolumeId()));
                if (m_version > 1)
                    PANDORA_RETURN_RESULT_IF(
                        pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("DaughterVolumeId", pLArCaloHit->GetDaughterVolumeId()));
                if (m_version > 2)
                {
                    const pandora::FloatVector &hitScores(pLArCaloHit->GetHitScores());
                    const pandora::StringVector &hitScoreLabels(pLArCaloHit->GetHitScoreLabels());
                    const unsigned int nLabels(hitScores.size());
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NHitLabels", nLabels));
                    for (unsigned int i = 0; i < nLabels; ++i)
                    {
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("HitScore", hitScores.at(i)));
                        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("HitScoreLabel", hitScoreLabels.at(i)));
                    }
                }
            }
            else
            {
                return pandora::STATUS_CODE_INVALID_PARAMETER;
            }

            return pandora::STATUS_CODE_SUCCESS;
        }
    }
}

} // namespace lar_content

#endif // #ifndef LAR_CALO_HIT_H
