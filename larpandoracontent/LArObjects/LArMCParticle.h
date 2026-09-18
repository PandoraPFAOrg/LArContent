/**
 *  @file   larpandoracontent/LArObjects/LArMCParticle.h
 *
 *  @brief  Header file for the lar mc particle class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_PARTICLE_H
#define LAR_MC_PARTICLE_H 1

#include "Objects/MCParticle.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/PandoraObjectFactories.h"

#include "Persistency/FieldMap.h"

namespace lar_content
{

// Enumeration maps onto G4 process IDs from QGSP_BERT and EM standard physics lists, plus an ID for the incident neutrino
enum MCProcess
{
    MC_PROC_INCIDENT_NU = -1,
    MC_PROC_UNKNOWN,
    MC_PROC_PRIMARY,
    MC_PROC_COMPT,
    MC_PROC_PHOT,
    MC_PROC_ANNIHIL,
    MC_PROC_E_IONI,
    MC_PROC_E_BREM,
    MC_PROC_CONV,
    MC_PROC_MU_IONI,
    MC_PROC_MU_MINUS_CAPTURE_AT_REST,
    MC_PROC_NEUTRON_INELASTIC,
    MC_PROC_N_CAPTURE,
    MC_PROC_HAD_ELASTIC,
    MC_PROC_DECAY,
    MC_PROC_COULOMB_SCAT,
    MC_PROC_MU_BREM,
    MC_PROC_MU_PAIR_PROD,
    MC_PROC_PHOTON_INELASTIC,
    MC_PROC_HAD_IONI,
    MC_PROC_PROTON_INELASTIC,
    MC_PROC_PI_PLUS_INELASTIC,
    MC_PROC_CHIPS_NUCLEAR_CAPTURE_AT_REST,
    MC_PROC_PI_MINUS_INELASTIC,
    MC_PROC_TRANSPORTATION,
    MC_PROC_RAYLEIGH,
    MC_PROC_HAD_BREM,
    MC_PROC_HAD_PAIR_PROD,
    MC_PROC_ION_IONI,
    MC_PROC_NEUTRON_KILLER,
    MC_PROC_ION_INELASTIC,
    MC_PROC_HE3_INELASTIC,
    MC_PROC_ALPHA_INELASTIC,
    MC_PROC_ANTI_HE3_INELASTIC,
    MC_PROC_ANTI_ALPHA_INELASTIC,
    MC_PROC_HAD_FRITIOF_CAPTURE_AT_REST,
    MC_PROC_ANTI_DEUTERON_INELASTIC,
    MC_PROC_ANTI_NEUTRON_INELASTIC,
    MC_PROC_ANTI_PROTON_INELASTIC,
    MC_PROC_ANTI_TRITON_INELASTIC,
    MC_PROC_DEUTERON_INELASTIC,
    MC_PROC_ELECTRON_NUCLEAR,
    MC_PROC_PHOTON_NUCLEAR,
    MC_PROC_KAON_PLUS_INELASTIC,
    MC_PROC_KAON_MINUS_INELASTIC,
    MC_PROC_HAD_BERTINI_CAPTURE_AT_REST,
    MC_PROC_LAMBDA_INELASTIC,
    MC_PROC_MU_NUCLEAR,
    MC_PROC_TRITON_INELASTIC,
    MC_PROC_PRIMARY_BACKGROUND,
    MC_PROC_RADIOACTIVE_DECAY_BASE
};

/**
 *  @brief  LAr mc particle parameters
 */
class LArMCParticleParameters : public object_creation::MCParticle::Parameters
{
public:
    pandora::InputInt m_nuanceCode;               ///< The nuance code
    pandora::InputInt m_process;                  ///< The process creating the particle
    pandora::InputBool m_isCC;                    ///< Whether the neutrino interacts via a CC interaction
    pandora::InputFloat m_visibleEnergy;          ///< Energy 'seen' in the detector
    pandora::InputCartesianVector m_endDirection; ///< From momentum at penultimate trajectory point
    pandora::InputInt m_nTrajPoints;              ///< Number of trajectory points
    pandora::CartesianPointVector m_trajPoints;   ///< The MCParticle trajectory points
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LAr mc particle class
 */
class LArMCParticle : public object_creation::MCParticle::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar mc particle parameters
     */
    LArMCParticle(const LArMCParticleParameters &parameters);

    /**
     *  @brief  Get the nuance code
     *
     *  @return the nuance code
     */
    int GetNuanceCode() const;

    /**
     *  @brief  Fill the parameters associated with this MC particle
     *
     *  @param  parameters the output parameters
     */
    void FillParameters(LArMCParticleParameters &parameters) const;

    /**
     *  @brief  Get the process
     *
     *  @return the process
     */
    MCProcess GetProcess() const;

    /**
     *  @brief  Get Whether the neutrino interacts via a CC interaction
     *          ATTN: always false for non-neutrinos
     *
     *  @return whether the neutrino interacts via a CC interaction
     */
    bool GetIsCC() const;

    /**
     *  @brief  Get the visible energy
     *
     *  @return the visible energy
     */
    float GetVisibleEnergy() const;

    /**
     *  @brief  Get the end direction
     *
     *  @return the end direction
     */
    pandora::CartesianVector GetEndDirection() const;

    /**
     *  @brief  Get the number of trajectory points
     *
     *  @return the number of trajectory points
     */
    int GetNTrajPoints() const;

    /**
     *  @brief  Get the trajectory points of the MC trajectory
     *          note: trajectories are sparsified in the current sim workflow
     *
     *  @return the vector of trajectory points
     */
    const pandora::CartesianPointVector &GetTrajPoints() const;

private:
    int m_nuanceCode;                           ///< The nuance code
    int m_process;                              ///< The process creating the particle
    bool m_isCC;                                ///< Whether the neutrino interacts via a CC interaction
    float m_visibleEnergy;                      ///< Energy 'seen' in the detector
    pandora::CartesianVector m_endDirection;    ///< From momentum at penultimate trajectory point
    int m_nTrajPoints;                          ///< Number of trajectory points
    pandora::CartesianPointVector m_trajPoints; ///< The MCParticle trajectory points
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArMCParticleFactory responsible for object creation
 */
class LArMCParticleFactory : public pandora::ObjectFactory<object_creation::MCParticle::Parameters, object_creation::MCParticle::Object>
{
public:
    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read additional (LArMCParticle-specific) parameters from the FieldMap.
     *          Fields: nuanceCode, process, isCC, visibleEnergy, endDirection,
     *                  nTrajPoints, trajPoint_N (CartesianVector, indexed).
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

inline LArMCParticle::LArMCParticle(const LArMCParticleParameters &parameters) :
    object_creation::MCParticle::Object(parameters),
    m_nuanceCode(parameters.m_nuanceCode.Get()),
    m_process(parameters.m_process.Get()),
    m_isCC(parameters.m_isCC.Get()),
    m_visibleEnergy(parameters.m_visibleEnergy.Get()),
    m_endDirection(parameters.m_endDirection.Get()),
    m_nTrajPoints(parameters.m_nTrajPoints.Get()),
    m_trajPoints(parameters.m_trajPoints)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArMCParticle::GetNuanceCode() const
{
    return m_nuanceCode;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MCProcess LArMCParticle::GetProcess() const
{
    return MCProcess(m_process);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArMCParticle::GetIsCC() const
{
    return m_isCC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArMCParticle::GetVisibleEnergy() const
{
    return m_visibleEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArMCParticle::GetEndDirection() const
{
    return m_endDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArMCParticle::GetNTrajPoints() const
{
    return m_nTrajPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianPointVector &LArMCParticle::GetTrajPoints() const
{
    return m_trajPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArMCParticle::FillParameters(LArMCParticleParameters &parameters) const
{
    parameters.m_nuanceCode = this->GetNuanceCode();
    parameters.m_process = this->GetProcess();
    parameters.m_isCC = this->GetIsCC();
    parameters.m_visibleEnergy = this->GetVisibleEnergy();
    parameters.m_endDirection = this->GetEndDirection();
    parameters.m_nTrajPoints = this->GetNTrajPoints();
    parameters.m_trajPoints = this->GetTrajPoints();
    parameters.m_energy = this->GetEnergy();
    parameters.m_momentum = this->GetMomentum();
    parameters.m_vertex = this->GetVertex();
    parameters.m_endpoint = this->GetEndpoint();
    parameters.m_particleId = this->GetParticleId();
    parameters.m_mcParticleType = this->GetMCParticleType();
    // ATTN Set the parent address to the original owner of the mc particle
    parameters.m_pParentAddress = static_cast<const void *>(this);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticleFactory::Parameters *LArMCParticleFactory::NewParameters() const
{
    return (new LArMCParticleParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArMCParticleParameters &larMCParticleParameters(dynamic_cast<const LArMCParticleParameters &>(parameters));
    pObject = new LArMCParticle(larMCParticleParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Read(Parameters &parameters, const pandora::FieldMap &fields) const
{
    LArMCParticleParameters &p(dynamic_cast<LArMCParticleParameters &>(parameters));

    p.m_nuanceCode = fields.GetOrDefault<int>("nuanceCode", 0);
    p.m_process = fields.GetOrDefault<int>("process", static_cast<int>(MC_PROC_UNKNOWN));
    p.m_isCC = fields.GetOrDefault<bool>("isCC", false);
    p.m_visibleEnergy = fields.GetOrDefault<float>("visibleEnergy", 0.f);
    p.m_endDirection = fields.GetOrDefault<pandora::CartesianVector>("endDirection", pandora::CartesianVector(0.f, 0.f, 0.f));

    const int nTrajPoints = std::max(0, fields.GetOrDefault<int>("nTrajPoints", 0));
    p.m_nTrajPoints = nTrajPoints;
    p.m_trajPoints.clear();
    p.m_trajPoints.reserve(static_cast<std::size_t>(nTrajPoints));

    for (int i = 0; i < nTrajPoints; ++i)
    {
        p.m_trajPoints.push_back(
            fields.GetOrDefault<pandora::CartesianVector>("trajPoint_" + std::to_string(i), pandora::CartesianVector(0.f, 0.f, 0.f)));
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Write(const Object *const pObject, pandora::FieldMap &fields) const
{
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pObject));

    if (!pLArMCParticle)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    fields.Set("nuanceCode", pLArMCParticle->GetNuanceCode());
    fields.Set("process", static_cast<int>(pLArMCParticle->GetProcess()));
    fields.Set("isCC", pLArMCParticle->GetIsCC());
    fields.Set("visibleEnergy", pLArMCParticle->GetVisibleEnergy());
    fields.Set("endDirection", pLArMCParticle->GetEndDirection());

    const pandora::CartesianPointVector &trajPoints(pLArMCParticle->GetTrajPoints());
    const int nTrajPoints(pLArMCParticle->GetNTrajPoints());

    fields.Set("nTrajPoints", nTrajPoints);

    for (int i = 0; i < nTrajPoints; ++i)
        fields.Set("trajPoint_" + std::to_string(i), trajPoints.at(static_cast<std::size_t>(i)));

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_H
