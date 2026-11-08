/**
 *  @file   larpandoracontent/LArCheating/CheatingMvaPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArCheating/CheatingMvaPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingMvaPfoCharacterisationAlgorithm::CheatingMvaPfoCharacterisationAlgorithm() :
    m_persistFeatures(false),
    m_useCaloHitsMatching(false),
    m_criticalEnergy(0.0305) // LAr Critical energy, standard value from https://lar.bnl.gov/properties/ (GeV)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingMvaPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingMvaPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    if (!LArPfoHelper::IsThreeD(pPfo))
    {
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["TrackScore"] = -1.f;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        
        return (pPfo->GetParticleId() == MU_MINUS);
    }
 
    // Charge related features are only calculated using hits in W view
    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    const PfoCharacterisationFeatureTool::FeatureToolMap &chosenFeatureToolMap(wClusterList.empty() ? m_featureToolMapNoChargeInfo : m_featureToolMapThreeD);
    const StringVector chosenFeatureToolOrder(wClusterList.empty() ? m_algorithmToolNamesNoChargeInfo : m_algorithmToolNames);
    StringVector featureOrder;
    const LArMvaHelper::MvaFeatureMap featureMap(
        LArMvaHelper::CalculateFeatures(chosenFeatureToolOrder, chosenFeatureToolMap, featureOrder, this, pPfo));

    for (auto const &[featureKey, featureValue] : featureMap)
    {
        (void)featureKey;

        if (!featureValue.IsInitialized())
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["TrackScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            return (pPfo->GetParticleId() == MU_MINUS);
        }
    }

    bool isTrueTrack(false);
    bool isMainMCParticleSet(false);
    bool isTrueTrackFromHits(false);

    try
    {
        // This define the truth trackScore value for the Pfo: 1 tracks, 0 showers
        const MCParticle *const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
        isTrueTrack = (
            (PHOTON != pMCParticle->GetParticleId()) && (
                E_MINUS != std::abs(pMCParticle->GetParticleId()) || pMCParticle->GetEnergy() <= m_criticalEnergy
            )
        );
        isMainMCParticleSet = (pMCParticle->GetParticleId() != 0);

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingMvaPfoCharacterisationAlgorithm: "
                      << "The MC particle found is with PID = " << pMCParticle->GetParticleId() 
                      << ", assigned TrackScore == " << static_cast<int>(isTrueTrack) << " ("
                      << (isTrueTrack ? "track" : "shower") << "), has Energy == " << pMCParticle->GetEnergy() << " GeV " << std::endl;

    }
    catch (const StatusCodeException &)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingMvaPfoCharacterisationAlgorithm: WARNING: I was unable to find the main MC particle for this Pfo" << std::endl;

    }

    CaloHitList pHitListAll;
    LArPfoHelper::GetAllCaloHits(pPfo, pHitListAll);
    
    const MCParticle *const pMCParticleFromHits(MCParticleHelper::GetMainMCParticle(&pHitListAll));
    isTrueTrackFromHits = (
        (PHOTON != pMCParticleFromHits->GetParticleId()) && (
            E_MINUS != std::abs(pMCParticleFromHits->GetParticleId()) || pMCParticleFromHits->GetEnergy() <= m_criticalEnergy
        )
    );
    
    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "CheatingMvaPfoCharacterisationAlgorithm: using calo hits matching: "
                  << " assigned TrackScore (from hits) == " << static_cast<int>(isTrueTrackFromHits) << " ("
                  << (isTrueTrackFromHits ? "track" : "shower") << ")" << std::endl;
    

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo() && !isMainMCParticleSet)
        std::cout << "CheatingMvaPfoCharacterisationAlgorithm: WARNING: No MC main particle found..." << std::endl;

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "CheatingMvaPfoCharacterisationAlgorithm: DEBUG: "
                  << " isTrueTrack == " << (isTrueTrack ? "true" : "false") 
                  << ", isTrueTrackFromHits == " << (isTrueTrackFromHits ? "true" : "false") 
                  << std::endl;

    if (isMainMCParticleSet && !m_useCaloHitsMatching)
    {
        // I found a MC particle (the main MC particle is set) so I can return the correct value
        // SIDE NOTE: if this association do not happen, the track score is not present in the pfp metadata,
        // so CAFAna will assign a default value of -5

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingMvaPfoCharacterisationAlgorithm: DEBUG: using pPfo truth matching" << std::endl;

        object_creation::ParticleFlowObject::Metadata metadata;
        const double score = static_cast<int>(isTrueTrack);
        metadata.m_propertiesToAdd["TrackScore"] = score;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingMvaPfoCharacterisationAlgorithm: assigned score == " << score << std::endl;

        if (m_persistFeatures)
        {
            for (auto const &[name, value] : featureMap)
            {
                metadata.m_propertiesToAdd[name] = value.Get();
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "CheatingMvaPfoCharacterisationAlgorithm: adding property : " << name << " with value " << value.Get() << std::endl;
            }
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        return isTrueTrack;
    }
    else
    {
        /* Here the LArMCParticleHelper::GetMainMCParticle somewhat failed, and so the reconstruction needs a more sofisticated procedure */

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingMvaPfoCharacterisationAlgorithm: DEBUG: using pHitListAll truth matching" << std::endl;
        
        isTrueTrack = isTrueTrackFromHits;

        object_creation::ParticleFlowObject::Metadata metadata;
        const double score = static_cast<int>(isTrueTrack);
        metadata.m_propertiesToAdd["TrackScore"] = score;

        if (m_persistFeatures)
        {
            for (auto const &[name, value] : featureMap)
            {
                metadata.m_propertiesToAdd[name] = value.Get();
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "CheatingMvaPfoCharacterisationAlgorithm: adding property : " << name << " with value " << value.Get() << std::endl;
            }
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        return isTrueTrack;
    }

    return isTrueTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingMvaPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Optional XML labels
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PersistFeatures", m_persistFeatures));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "UseCaloHitsMatching", m_useCaloHitsMatching));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CriticalEnergy", m_criticalEnergy));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));


    LArMvaHelper::AlgorithmToolMap algorithmToolMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        LArMvaHelper::ProcessAlgorithmToolListToMap(*this, xmlHandle, "FeatureTools", m_algorithmToolNames, algorithmToolMap));

    if (m_useThreeDInformation)
    {
        // and the map for NoChargeInfo
        LArMvaHelper::AlgorithmToolMap algorithmToolMapNoChargeInfo;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            LArMvaHelper::ProcessAlgorithmToolListToMap(
                *this, xmlHandle, "FeatureToolsNoChargeInfo", m_algorithmToolNamesNoChargeInfo, algorithmToolMapNoChargeInfo));

        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
            PANDORA_RETURN_RESULT_IF(
                STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMapThreeD));

        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMapNoChargeInfo)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMapNoChargeInfo));
    }
    else
    {
        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMap));
    }

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
