/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h
 *
 *  @brief  Header file for the shower hits base tool.
 *
 *  $Log: $
 */
#ifndef SHOWER_HITS_BASE_TOOL_H
#define SHOWER_HITS_BASE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ShowerHitsBaseTool class
 */
class ShowerHitsBaseTool : public HitCreationBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerHitsBaseTool();

    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector);

protected:
    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a list of candidate matched
     *          hits in the other two views
     *
     *  @param  caloHitVector1 the vector of candidate hits in view 1
     *  @param  caloHitVector2 the vector of candidate hits in view 2
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetShowerHit3D(const pandora::CaloHitVector &caloHitVector1, const pandora::CaloHitVector &caloHitVector2, ProtoHit &protoHit) const = 0;

    /**
     *  @brief  Create three dimensional hits, using a list of input two dimensional hits and the hits (contained in the same particle)
     *          from the other two views
     *
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  caloHitVector1 hits in the first alternate view
     *  @param  caloHitVector2 hits in the second alternate view
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    virtual void GetShowerHits3D(const pandora::CaloHitVector &inputTwoDHits, const pandora::CaloHitVector &caloHitVector1,
        const pandora::CaloHitVector &caloHitVector2, ProtoHitVector &protoHitVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief  Filter a list of calo hits to find those within a specified tolerance of a give x position
     *
     *  @param  x the x position
     *  @param  daughterVolumeId the daughter volume id of the hit in question
     *  @param  inputCaloHitVector the input calo hit vector
     *  @param  outputCaloHitVector to receive the output calo hit vector
     */
    void FilterCaloHits(const float x, const unsigned int daughterVolumeId, const pandora::CaloHitVector &inputCaloHitVector,
        pandora::CaloHitVector &outputCaloHitVector) const;

    float m_xTolerance;    ///< The x tolerance to use when looking for associated calo hits between views
    bool m_strictBounding; ///< Whether to apply daughter volume match requirements
};

} // namespace lar_content

#endif // #ifndef SHOWER_HITS_BASE_TOOL_H
