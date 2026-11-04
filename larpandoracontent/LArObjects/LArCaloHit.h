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

#include "Persistency/FieldMap.h"

namespace lar_content
{

/**
 *  @brief  LAr calo hit parameters
 */
class LArCaloHitParameters : public object_creation::CaloHit::Parameters
{
public:
    pandora::InputUInt m_larTPCVolumeId;    ///< The lar tpc volume id
    pandora::InputUInt m_daughterVolumeId;  ///< The daughter volume id
    pandora::InputUInt m_channelId;         ///< The channel (wire) id
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
    LArCaloHit(const LArCaloHitParameters &parameters);

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
     *  @brief  Get the channel (wire) id
     *
     *  @return the channel id
     */
    unsigned int GetChannelId() const;

    /**
     *  @brief  Fill the parameters associated with this calo hit
     *
     *  @param  parameters the output parameters
     */
    void FillParameters(LArCaloHitParameters &parameters) const;

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
     *  @brief  Get the MIP score for the hit
     *
     *  @return the hit's MIP score
     */        
    float GetMIPScore() const;

    /**
     *  @brief  Get the HIP score for the hit
     *
     *  @return the hit's HIP score
     */           
    float GetHIPScore() const;

    /**
     *  @brief  Get the EM-track score for the hit
     *
     *  @return the hit's EM-track score
     */       
    float GetEMTrackScore() const;

    /**
     *  @brief  Get the EM-shower score for the hit
     *
     *  @return the hit's EM-shower score
     */    
    float GetEMShowerScore() const;

    /**
     *  @brief  Get the diffuse score for the hit
     *
     *  @return the hit's diffuse score
     */
    float GetDiffuseScore() const;        

    /**
     *  @brief  Get the track-like score for the hit
     *
     *  @return the hit's track-like score
     */
    float GetTrackScore() const;

    /**
     *  @brief  Get the shower-like score for the hit
     *
     *  @return the hit's shower-like score
     */    
    float GetShowerScore() const;    
    
    /**
     *  @brief  Add a hit score for a given label to the hit
     *
     *  @param  hitScore the score to add
     *  @param  scoreLabel the label for the score to add
     */
    void AddHitScore(const float &hitScore, const std::string &scoreLabel);

private:
    unsigned int m_larTPCVolumeId;          ///< The lar tpc volume id
    unsigned int m_daughterVolumeId;        ///< The daughter volume id
    unsigned int m_channelId;               ///< The channel (wire) id
    pandora::FloatVector m_hitScores;       ///< Hit scores
    pandora::StringVector m_hitScoreLabels; ///< Labels for the hit scores
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArCaloHitFactory responsible for object creation
 */
class LArCaloHitFactory : public pandora::ObjectFactory<object_creation::CaloHit::Parameters, object_creation::CaloHit::Object>
{
public:
    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read additional (LArCaloHit-specific) parameters from the FieldMap.
     *          Fields: larTPCVolumeId, daughterVolumeId, nHitScores,
     *                  hitScore_N, hitScoreLabel_N (indexed per score).
     *          All fields use GetOrDefault so files without them are handled gracefully.
     */
    pandora::StatusCode Read(Parameters &parameters, const pandora::FieldMap &fields) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fileWriter the file writer
     */
    pandora::StatusCode Write(const Object *const pObject, pandora::FieldMap &fields) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const Parameters &parameters, const Object *&pObject) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHit::LArCaloHit(const LArCaloHitParameters &parameters) :
    object_creation::CaloHit::Object(parameters),
    m_larTPCVolumeId(parameters.m_larTPCVolumeId.Get()),
    m_daughterVolumeId(parameters.m_daughterVolumeId.IsInitialized() ? parameters.m_daughterVolumeId.Get() : 0),
    m_channelId(parameters.m_channelId.IsInitialized() ? parameters.m_channelId.Get() : 0),
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

inline unsigned int LArCaloHit::GetChannelId() const
{
    return m_channelId;
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

inline float LArCaloHit::GetMIPScore() const
{
    const auto iter(std::find(m_hitScoreLabels.begin(), m_hitScoreLabels.end(), "MIPScore"));

    if (iter == m_hitScoreLabels.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    

    return (m_hitScores.at(iter - m_hitScoreLabels.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetHIPScore() const
{
    const auto iter(std::find(m_hitScoreLabels.begin(), m_hitScoreLabels.end(), "HIPScore"));

    if (iter == m_hitScoreLabels.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    

    return (m_hitScores.at(iter - m_hitScoreLabels.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetEMTrackScore() const
{
    const auto iter(std::find(m_hitScoreLabels.begin(), m_hitScoreLabels.end(), "EMTrackScore"));

    if (iter == m_hitScoreLabels.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    

    return (m_hitScores.at(iter - m_hitScoreLabels.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetEMShowerScore() const
{
    const auto iter(std::find(m_hitScoreLabels.begin(), m_hitScoreLabels.end(), "EMShowerScore"));

    if (iter == m_hitScoreLabels.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    

    return (m_hitScores.at(iter - m_hitScoreLabels.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetDiffuseScore() const
{
    const auto iter(std::find(m_hitScoreLabels.begin(), m_hitScoreLabels.end(), "DiffuseScore"));

    if (iter == m_hitScoreLabels.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    

    return (m_hitScores.at(iter - m_hitScoreLabels.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetTrackScore() const
{
    return (this->GetMIPScore() + this->GetHIPScore());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetShowerScore() const
{
    return (this->GetEMTrackScore() + this->GetEMShowerScore() + this->GetDiffuseScore());
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::FillParameters(LArCaloHitParameters &parameters) const
{
    parameters.m_cellGeometry = this->GetCellGeometry();
    parameters.m_positionVector = this->GetPositionVector();
    parameters.m_expectedDirection = this->GetExpectedDirection();
    parameters.m_cellNormalVector = this->GetCellNormalVector();
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
    parameters.m_channelId = this->GetChannelId();
    parameters.m_hitScores = this->GetHitScores();
    parameters.m_hitScoreLabels = this->GetHitScoreLabels();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::AddHitScore(const float &hitScore, const std::string &scoreLabel)
{
    if (std::find(m_hitScoreLabels.begin(), m_hitScoreLabels.end(), scoreLabel) != m_hitScoreLabels.end())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    m_hitScores.emplace_back(hitScore);
    m_hitScoreLabels.emplace_back(scoreLabel);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHitFactory::Parameters *LArCaloHitFactory::NewParameters() const
{
    return (new LArCaloHitParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArCaloHitParameters &larCaloHitParameters(dynamic_cast<const LArCaloHitParameters &>(parameters));
    pObject = new LArCaloHit(larCaloHitParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Read(Parameters &parameters, const pandora::FieldMap &fields) const
{
    LArCaloHitParameters &p(dynamic_cast<LArCaloHitParameters &>(parameters));

    p.m_larTPCVolumeId = fields.GetOrDefault<unsigned int>("larTPCVolumeId", 0u);
    p.m_daughterVolumeId = fields.GetOrDefault<unsigned int>("daughterVolumeId", 0u);
    p.m_channelId = fields.GetOrDefault<unsigned int>("channelId", 0u);

    const unsigned int nHitScores = fields.GetOrDefault<unsigned int>("nHitScores", 0u);
    p.m_hitScores.clear();
    p.m_hitScoreLabels.clear();
    p.m_hitScores.reserve(nHitScores);
    p.m_hitScoreLabels.reserve(nHitScores);

    for (unsigned int i = 0; i < nHitScores; ++i)
    {
        p.m_hitScores.push_back(fields.GetOrDefault<float>("hitScore_" + std::to_string(i), 0.f));
        p.m_hitScoreLabels.push_back(fields.GetOrDefault<std::string>("hitScoreLabel_" + std::to_string(i), std::string()));
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Write(const Object *const pObject, pandora::FieldMap &fields) const
{
    const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pObject));

    if (!pLArCaloHit)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    fields.Set("larTPCVolumeId", pLArCaloHit->GetLArTPCVolumeId());
    fields.Set("daughterVolumeId", pLArCaloHit->GetDaughterVolumeId());
    fields.Set("channelId", pLArCaloHit->GetChannelId());

    const pandora::FloatVector &hitScores(pLArCaloHit->GetHitScores());
    const pandora::StringVector &hitScoreLabels(pLArCaloHit->GetHitScoreLabels());
    const unsigned int nHitScores(static_cast<unsigned int>(hitScores.size()));

    fields.Set("nHitScores", nHitScores);

    for (unsigned int i = 0; i < nHitScores; ++i)
    {
        fields.Set("hitScore_" + std::to_string(i), hitScores.at(i));
        fields.Set("hitScoreLabel_" + std::to_string(i), hitScoreLabels.at(i));
    }

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_CALO_HIT_H
