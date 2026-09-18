/**
 *  @file   larpandoracontent/LArHelpers/LArPfoHelper.h
 *
 *  @brief  Header file for the pfo helper class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_HELPER_H
#define LAR_PFO_HELPER_H 1

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Vertex.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

namespace lar_content
{

/**
 *  @brief  LArPfoHelper class
 */
class LArPfoHelper
{
public:
    /**
     *  @brief  Get a list of coordinates of a particular hit type from an input pfos
     *
     *  @param  pPfo the address of the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  coordinateVector the output list of coordinates
     */
    static void GetCoordinateVector(
        const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::CartesianPointVector &coordinateVector);

    /**
     *  @brief  Get a list of calo hits of a particular hit type from a list of pfos
     *
     *  @param  pfoList the input list of Pfos
     *  @param  hitType the cluster hit type
     *  @param  caloHitList the output list of calo hits
     */
    static void GetCaloHits(const pandora::PfoList &pfoList, const pandora::HitType &hitType, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of calo hits of a particular hit type from a given pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  caloHitList the output list of calo hits
     */
    static void GetCaloHits(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of isolated calo hits of a particular hit type from a list of pfos
     *
     *  @param  pfoList the input list of Pfos
     *  @param  hitType the cluster hit type
     *  @param  caloHitList the output list of calo hits
     */
    static void GetIsolatedCaloHits(const pandora::PfoList &pfoList, const pandora::HitType &hitType, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of isolated calo hits of a particular hit type from a given pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  caloHitList the output list of isolated calo hits
     */
    static void GetIsolatedCaloHits(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of all 2D calo hits (including isolated) of all types from a given pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  caloHitList the output list of calo hits
     */
    static void GetAllCaloHits2D(const pandora::ParticleFlowObject *pPfo, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of clusters of a particular hit type from a list of pfos
     *
     *  @param  pfoList the input list of Pfos
     *  @param  hitType the cluster hit type
     *  @param  clusterList the output list of clusters
     */
    static void GetClusters(const pandora::PfoList &pfoList, const pandora::HitType &hitType, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get a list of clusters of a particular hit type from a given pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  clusterList the output list of clusters
     */
    static void GetClusters(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::ClusterList &clusterList);

    /**
     * @brief Get the number of 2D hits of a PFO
     *
     * @param pPfo the pfo to check
     * @return int of number of 2D hits
     */
    static unsigned int GetNumberOfTwoDHits(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get the number of 3D hits owned by a pfo
     *
     *  @param  pPfo a pointer to the pfo
     *  @return The number of 3D hits
     */
    static unsigned int GetNumberOfThreeDHits(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief Get the list of 2D clusters from an input pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  clusterList the output list of clusters
     */
    static void GetTwoDClusterList(const pandora::ParticleFlowObject *const pPfo, pandora::ClusterList &clusterList);

    /**
     *  @brief Get the list of 3D clusters from an input pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  clusterList the output list of clusters
     */
    static void GetThreeDClusterList(const pandora::ParticleFlowObject *const pPfo, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get a flat list of all pfos, recursively including all daughters and parents associated with those pfos in an input list
     *
     *  @param  inputPfoList the input pfo list
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllConnectedPfos(const pandora::PfoList &inputPfoList, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get a flat list of all pfos, recursively including all daughters and parents associated with an input pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllConnectedPfos(const pandora::ParticleFlowObject *const pPfo, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get a flat list of all pfos, recursively, of all daughters associated with those pfos in an input list
     *
     *  @param  inputPfoList the input pfo list
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllDownstreamPfos(const pandora::PfoList &inputPfoList, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get a flat list of all pfos, recursively, of all daughters and parents associated with an input pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllDownstreamPfos(const pandora::ParticleFlowObject *const pPfo, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get flat lists of all downstream track pfos and also shower-like pfos.
     *          This method collects together all track-like particles downstream of the root particle, stopping at a leading shower and
     *          then storing that leading shower in a separate list.
     *
     *  @param  pPfo the input pfo
     *  @param  outputTrackPfoList the output list of descendent track-like particles
     *  @param  outputLeadingShowerParticles the output list of leading shower particles
     */
    static void GetAllDownstreamPfos(
        const pandora::ParticleFlowObject *const pPfo, pandora::PfoList &outputTrackPfoList, pandora::PfoList &outputLeadingShowerPfoList);

    /**
     *  @brief  Determine the position in the hierarchy for the MCParticle
     *
     *  @param  pPfo the input Pfo
     *
     *  @return integer
     */
    static int GetHierarchyTier(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Calculate length of Pfo using 2D clusters
     *
     *  @param  pPfo the input Pfo
     *
     *  @return  length variable
     */
    static float GetTwoDLengthSquared(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Calculate length of Pfo using 3D clusters
     *
     *  @param  pPfo the input Pfo
     *
     *  @return  length variable
     */
    static float GetThreeDLengthSquared(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get closest distance between Pfo and cluster
     *
     *  @param  pPfo the address of the input Pfo
     *  @param  pCluster the address of the input cluster
     */
    static float GetClosestDistance(const pandora::ParticleFlowObject *const pPfo, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get distance between two Pfos using 3D clusters
     *
     *  @param  pPfo the address of the first Pfo
     *  @param  pPfo the address of the second Pfo
     */
    static float GetThreeDSeparation(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2);

    /**
     *  @brief  Does Pfo contain 2D clusters
     *
     *  @param  pPfo the address of the Pfo
     */
    static bool IsTwoD(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Does Pfo contain 3D clusters
     *
     *  @param  pPfo the address of the Pfo
     */
    static bool IsThreeD(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Return track flag based on Pfo Particle ID
     *
     *  @param  pPfo the address of the Pfo
     */
    static bool IsTrack(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Return shower flag based on Pfo Particle ID
     *
     *  @param  pPfo the address of the Pfo
     */
    static bool IsShower(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get the track/shower score of the pfo
     *
     *  @param  pPfo the address pf the pfo
     *  @return the track score
     */
    static float GetTrackScore(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get primary neutrino or antineutrino
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return pdg code of neutrino (or zero, otherwise)
     */
    static int GetPrimaryNeutrino(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a primary parent particle
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsFinalState(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a final-state particle from a neutrino (or antineutrino) interaction
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsNeutrinoFinalState(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a neutrino or (antineutrino)
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsNeutrino(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a final-state particle from a test beam particle interaction
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsTestBeamFinalState(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a test beam particle
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsTestBeam(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get neutrino pfos from an input pfo list
     *
     *  @param  pPfoList the input pfo list
     *  @param  recoNeutrinos to receive the list of neutrino pfos
     */
    static void GetRecoNeutrinos(const pandora::PfoList *const pPfoList, pandora::PfoList &recoNeutrinos);

    /**
     *  @brief  Get the primary parent pfo
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return address of the primary parent pfo
     */
    static const pandora::ParticleFlowObject *GetParentPfo(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get primary neutrino or antineutrino
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return address of primary neutrino pfo
     */
    static const pandora::ParticleFlowObject *GetParentNeutrino(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get the pfo vertex
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return address of pfo vertex
     */
    static const pandora::Vertex *GetVertex(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get the pfo test beam interaction vertex
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return address of pfo vertex
     */
    static const pandora::Vertex *GetTestBeamInteractionVertex(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get the vertex with a specific vertex label in a given vertex list
     *
     *  @param  vertexList vertex list
     *  @param  vertexLabel target vertex label type
     *
     *  @return address of the desired vertex
     */
    static const pandora::Vertex *GetVertexWithLabel(const pandora::VertexList &vertexList, const pandora::VertexLabel vertexLabel);

    /**
     *  @brief  Apply 3D sliding fit to a set of 3D points and return track trajectory
     *
     *  @param  pointVector  the input list of 3D positions
     *  @param  vertexPosition  the input vertex position
     *  @param  layerWindow  size of half window for sliding linear fit
     *  @param  layerPitch  size of pitch for sliding linear fit
     *  @param  trackStateVector  the output track trajectory
     *  @param  pIndexVector lookup vector of spacepoint indices to store trajectory point sorting
     */
    static void GetSlidingFitTrajectory(const pandora::CartesianPointVector &pointVector, const pandora::CartesianVector &vertexPosition,
        const unsigned int layerWindow, const float layerPitch, LArTrackStateVector &trackStateVector, pandora::IntVector *const pIndexVector = nullptr);

    /**
     *  @brief  Apply 3D sliding fit to Pfo and return track trajectory
     *
     *  @param  pPfo  the address of the input Pfo
     *  @param  pVertex  the address of the input vertex
     *  @param  slidingFitHalfWindow  size of half window for sliding linear fit
     *  @param  layerPitch  size of pitch for sliding linear fit
     *  @param  trackStateVector  the output track trajectory
     */
    static void GetSlidingFitTrajectory(const pandora::ParticleFlowObject *const pPfo, const pandora::Vertex *const pVertex,
        const unsigned int slidingFitHalfWindow, const float layerPitch, LArTrackStateVector &trackStateVector);

    /**
     *  @brief  Apply 3D sliding fit to CaloHitList and return track trajectory
     *
     *  @param  pCaloHitList  the input list of calo hits
     *  @param  vertexPosition  the input vertex position
     *  @param  layerWindow  size of half window for sliding linear fit
     *  @param  layerPitch  size of pitch for sliding linear fit
     *  @param  trackStateVector  the output track trajectory
     *  @param  pIndexVector  lookup vector of spacepoint indices to store trajectory point sorting
     *  @param  return3DCaloHit  whether to tell TypeAdaptor to use the ParentAddress or use the 3D hit itself
     */
    static void GetSlidingFitTrajectory(const pandora::CaloHitList *const pCaloHitList, const pandora::CartesianVector &vertexPosition,
        const unsigned int layerWindow, const float layerPitch, LArTrackStateVector &trackStateVector,
        pandora::IntVector *const pIndexVector = nullptr, const bool return3DCaloHit = false);

    /**
     *  @brief  Perform PCA analysis on a set of 3D points and return results
     *
     *  @param  pointVector the input list of 3D positions
     *  @param  vertexPosition the input vertex position
     */
    static LArShowerPCA GetPrincipalComponents(const pandora::CartesianPointVector &pointVector, const pandora::CartesianVector &vertexPosition);

    /**
     *  @brief  Perform PCA analysis on Pfo and return results
     *
     *  @param  pPfo the address of the input Pfo
     *  @param  pVertex the address of the input vertex
     */
    static LArShowerPCA GetPrincipalComponents(const pandora::ParticleFlowObject *const pPfo, const pandora::Vertex *const pVertex);

    /**
     *  @brief  Sort pfos by number of constituent hits
     *
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByHitProjection(const LArTrackTrajectoryPoint &lhs, const LArTrackTrajectoryPoint &rhs);

    /**
     *  @brief  Sort pfos by number of constituent hits
     *
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    /**
     *  @brief  Retrieve a linearised representation of the PFO hierarchy in breadth first order. This iterates over the PFO hierarchy in a
     *          manor that sees primaries at the front of the list, with progressively deeper tiers later in the list. This is useful for
     *          some visualisation cases.
     *
     *  @param  pPfo a PFO in the hierarchy - can be any PFO
     *  @param  pfoList the output PFO list
     */
    static void GetBreadthFirstHierarchyRepresentation(const pandora::ParticleFlowObject *const pPfo, pandora::PfoList &pfoList);

private:
    /**
     *  @brief  Implementation of sliding fit trajectory extraction
     *
     *  @param  t the input information
     *  @param  pVertex the address of the input vertex
     *  @param  slidingFitHalfWindow  size of half window for sliding linear fit
     *  @param  layerPitch  size of pitch for sliding linear fit
     *  @param  trackStateVector the output track trajectory
     *  @param  pIndexVector lookup vector of spacepoint indices to store trajectory point sorting
     */
    template <typename T>
    static void SlidingFitTrajectoryImpl(const T *const pT, const pandora::CartesianVector &vertexPosition, const unsigned int layerWindow,
        const float layerPitch, LArTrackStateVector &trackStateVector, pandora::IntVector *const pIndexVector = nullptr,
        const bool return3DCaloHit = false);
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_HELPER_H
