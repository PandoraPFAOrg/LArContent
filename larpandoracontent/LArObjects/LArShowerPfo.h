/**
 *  @file   larpandoracontent/LArObjects/LArShowerPfo.h
 *
 *  @brief  Header file for the lar pfo class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_PFO_H
#define LAR_SHOWER_PFO_H 1

#include "Objects/ParticleFlowObject.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/ObjectFactory.h"

#include <string>

namespace lar_content
{

/**
 *  @brief  lar pfo parameters
 */
class LArShowerPfoParameters : public object_creation::ParticleFlowObject::Parameters
{
public:
    pandora::InputCartesianVector m_showerLength;          ///< Shower length and widths from 3d shower fit
    pandora::InputCartesianVector m_showerCentroid;        ///< Shower centroid from 3d shower fit
    pandora::InputFloat m_showerOpeningAngle;              ///< Shower opening angle
    pandora::InputCartesianVector m_showerDirection;       ///< Shower direction, also the primary eigen vector
    pandora::InputCartesianVector m_showerSecondaryVector; ///< Shower secondary eigen vector
    pandora::InputCartesianVector m_showerTertiaryVector;  ///< Shower teriary eigen vector
    pandora::InputCartesianVector m_showerEigenValues;     ///< Shower eigenvalues from 3d PCA
    pandora::InputCartesianVector m_showerVertex;          ///< Shower starting point
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object
 */
class LArShowerPfo : public object_creation::ParticleFlowObject::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar pfo parameters
     */
    LArShowerPfo(const LArShowerPfoParameters &parameters);

    /**
     *  @brief  Get the shower length and width from 3d shower fit
     *
     *  @return the shower length and width from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerLength() const;

    /**
     *  @brief  Get the shower centroid from the 3d shower fit
     *
     *  @return the shower centroid from the 3d shower fit
     */
    const pandora::CartesianVector &GetShowerCentroid() const;

    /**
     *  @brief  Get the shower opening angle from 3d shower fit
     *
     *  @return the shower opening angle from 3d shower fit
     */
    float GetShowerOpeningAngle() const;

    /**
     *  @brief  Get the shower direction, also the primary eigen vector from 3d shower fit
     *
     *  @return the shower direction, also the primary eigen vector from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerDirection() const;

    /**
     *  @brief  Get the shower secondary eigen vector from 3d shower fit
     *
     *  @return the shower secondary eigen vector from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerSecondaryVector() const;

    /**
     *  @brief  Get the shower tertiary eigen vector from 3d shower fit
     *
     *  @return the shower tertiary eigen vector from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerTertiaryVector() const;

    /**
     *  @brief  Get the shower eigen values from 3d PCA
     *
     *  @return the shower eigen values from 3d PCA
     */
    const pandora::CartesianVector &GetShowerEigenValues() const;

    /**
     *  @brief  Get the shower starting point from 3d shower fit
     *
     *  @return the shower starting point from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerVertex() const;

private:
    pandora::CartesianVector m_showerLength;          ///< Shower length and widths from 3d shower fit
    pandora::CartesianVector m_showerCentroid;        ///< Shower centroid from 3d shower fit
    float m_showerOpeningAngle;                       ///< Shower opening angle
    pandora::CartesianVector m_showerDirection;       ///< Shower direction, primary eigen vector
    pandora::CartesianVector m_showerSecondaryVector; ///< Shower secondary eigen vector
    pandora::CartesianVector m_showerTertiaryVector;  ///< Shower tertiary eigen vector
    pandora::CartesianVector m_showerEigenValues;     ///< Shower eigenvalues from 3d PCA
    pandora::CartesianVector m_showerVertex;          ///< Shower starting point
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object factory responsible for pfo creation
 */
class LArShowerPfoFactory : public pandora::ObjectFactory<object_creation::ParticleFlowObject::Parameters, object_creation::ParticleFlowObject::Object>
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

inline LArShowerPfo::LArShowerPfo(const LArShowerPfoParameters &parameters) :
    object_creation::ParticleFlowObject::Object(parameters),
    m_showerLength(parameters.m_showerLength.Get()),
    m_showerCentroid(parameters.m_showerCentroid.Get()),
    m_showerOpeningAngle(parameters.m_showerOpeningAngle.Get()),
    m_showerDirection(parameters.m_showerDirection.Get()),
    m_showerSecondaryVector(parameters.m_showerSecondaryVector.Get()),
    m_showerTertiaryVector(parameters.m_showerTertiaryVector.Get()),
    m_showerEigenValues(parameters.m_showerEigenValues.Get()),
    m_showerVertex(parameters.m_showerVertex.Get())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerLength() const
{
    return m_showerLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerCentroid() const
{
    return m_showerCentroid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArShowerPfo::GetShowerOpeningAngle() const
{
    return m_showerOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerDirection() const
{
    return m_showerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerSecondaryVector() const
{
    return m_showerSecondaryVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerTertiaryVector() const
{
    return m_showerTertiaryVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerEigenValues() const
{
    return m_showerEigenValues;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArShowerPfo::GetShowerVertex() const
{
    return m_showerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArShowerPfoFactory::Parameters *LArShowerPfoFactory::NewParameters() const
{
    return (new LArShowerPfoParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArShowerPfoFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArShowerPfoParameters &larPfoParameters(dynamic_cast<const LArShowerPfoParameters &>(parameters));
    pObject = new LArShowerPfo(larPfoParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArShowerPfoFactory::Read(Parameters &, const pandora::FieldMap &) const
{
    // TODO: Provide this functionality if necessary

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArShowerPfoFactory::Write(const pandora::ParticleFlowObject *, pandora::FieldMap &) const
{
    // TODO: Provide this functionality if necessary

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_PFO_H
