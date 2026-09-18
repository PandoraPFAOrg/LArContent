/**
 *  @file   larpandoracontent/LArObjects/LArTrackPfo.h
 *
 *  @brief  Header file for the lar pfo class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_PFO_H
#define LAR_TRACK_PFO_H 1

#include "Objects/CartesianVector.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/TrackState.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/ObjectFactory.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

namespace pandora
{
class CaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  lar pfo parameters
 */
class LArTrackPfoParameters : public object_creation::ParticleFlowObject::Parameters
{
public:
    LArTrackStateVector m_trackStateVector;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object
 */
class LArTrackPfo : public pandora::ParticleFlowObject
{
public:
    /**
     *  @brief  Constructor
     */
    LArTrackPfo(const LArTrackPfoParameters &parameters);

    /**
     *  @brief  Get vertex position
     */
    const pandora::CartesianVector &GetVertexPosition() const;

    /**
     *  @brief  Get end position
     */
    const pandora::CartesianVector &GetEndPosition() const;

    /**
     *  @brief  Get vertex direction
     */
    const pandora::CartesianVector &GetVertexDirection() const;

    /**
     *  @brief  Get end direction
     */
    const pandora::CartesianVector &GetEndDirection() const;

public:
    const LArTrackStateVector m_trackStateVector; ///< The vector of track states

private:
    // OTHER MEMBER VARIABLES GO HERE
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object factory responsible for pfo creation
 */
class LArTrackPfoFactory : public pandora::ObjectFactory<object_creation::ParticleFlowObject::Parameters, pandora::ParticleFlowObject>
{
public:
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
     *  @param  fields the field map to read from
     */
    pandora::StatusCode Read(Parameters &parameters, const pandora::FieldMap &fields) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fields the field map to write to
     */
    pandora::StatusCode Write(const pandora::ParticleFlowObject *const pObject, pandora::FieldMap &fields) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const object_creation::ParticleFlowObject::Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTrackPfoFactory::Parameters *LArTrackPfoFactory::NewParameters() const
{
    return (new LArTrackPfoParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArTrackPfoFactory::Create(const Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const
{
    const LArTrackPfoParameters &larPfoParameters(dynamic_cast<const LArTrackPfoParameters &>(parameters));
    pObject = new LArTrackPfo(larPfoParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArTrackPfoFactory::Read(Parameters &, const pandora::FieldMap &) const
{
    // TODO: Provide this functionality when necessary

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArTrackPfoFactory::Write(const pandora::ParticleFlowObject *, pandora::FieldMap &) const
{
    // TODO: Provide this functionality when necessary

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_PFO_H
