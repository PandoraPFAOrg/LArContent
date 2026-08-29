/**
 *  @file   larpandoracontent/LArCheating/CheatingMvaPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cheated mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_MVA_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_CHEATING_MVA_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
 *  @brief  CheatingMvaPfoCharacterisationAlgorithm class
 */
class CheatingMvaPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingMvaPfoCharacterisationAlgorithm();

protected:
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusterCharacterisationFeatureTool::FeatureToolMap m_featureToolMap; ///< The feature tool map

    PfoCharacterisationFeatureTool::FeatureToolMap m_featureToolMapThreeD;       ///< FeatureToolMap as a map for 3D info
    PfoCharacterisationFeatureTool::FeatureToolMap m_featureToolMapNoChargeInfo; ///< FeatureToolMap as a map for missing W view

    pandora::StringVector m_algorithmToolNames; ///< Vector of strings saving feature tool order for use in feature calculation
    pandora::StringVector m_algorithmToolNamesNoChargeInfo; ///< Vector of strings saving feature tool order for use in feature calculation (missing W view)

    bool m_persistFeatures;             ///< Whether to write the features to the properties map
    bool m_useCaloHitsMatching;         ///< Wheter to use calorimetric hits for the truth matching (Andy C. suggestion)  
    float m_criticalEnergy;

    std::string m_caloHitListName;    ///< Name of input calo hit list
    std::string m_mcParticleListName; ///< Name of input MC particle list

    LArMCParticleHelper::PrimaryParameters m_primaryParameters; ///< The mc particle primary selection parameters

};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_MVA_PFO_CHARACTERISATION_ALGORITHM_H
