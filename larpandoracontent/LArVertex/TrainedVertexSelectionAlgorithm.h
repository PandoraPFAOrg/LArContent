/**
 *  @file   larpandoracontent/LArVertex/TrainedVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the trained vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRAINED_VERTEX_SELECTION_ALGORITHM_H
#define LAR_TRAINED_VERTEX_SELECTION_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"
#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

#include <random>
#include <unordered_map>

namespace lar_content
{

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrainedVertexSelectionAlgorithm class
 */
class TrainedVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
{
public:
    /**
     *  @brief Vertex feature info class
     */
    class VertexFeatureInfo
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  beamDeweighting the beam deweighting feature
         *  @param  rPhiFeature the r/phi feature
         *  @param  energyKick the energy kick feature
         *  @param  localAsymmetry the local asymmetry feature
         *  @param  globalAsymmetry the global asymmetry feature
         *  @param  showerAsymmetry the shower asymmetry feature
         *  @param  dEdxAsymmetry the dE/dx asymmetry feature
         *  @param  vertexEnergy the vertex energy feature
         */
        VertexFeatureInfo(const float beamDeweighting, const float rPhiFeature, const float energyKick, const float localAsymmetry,
            const float globalAsymmetry, const float showerAsymmetry, const float dEdxAsymmetry, const float vertexEnergy);

        float m_beamDeweighting; ///< The beam deweighting feature
        float m_rPhiFeature;     ///< The r/phi feature
        float m_energyKick;      ///< The energy kick feature
        float m_localAsymmetry;  ///< The local asymmetry feature
        float m_globalAsymmetry; ///< The global asymmetry feature
        float m_showerAsymmetry; ///< The shower asymmetry feature
        float m_dEdxAsymmetry;   ///< The dE/dx asymmetry feature
        float m_vertexEnergy;    ///< The vertex energy feature
    };

    typedef std::map<const pandora::Vertex *const, VertexFeatureInfo> VertexFeatureInfoMap;

    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief Shared vertex feature info class
     */
    class VertexSharedFeatureInfo
    {
    public:
        /**
       *  @brief  Constructor
       *
       *  @param  m_separation the distance between the two vertices
       *  @param  m_axisHits the hit density along the axis between the two vertices
       */
        VertexSharedFeatureInfo(const float separation, const float axisHits);

        float m_separation; ///< The distance between the two vertices
        float m_axisHits;   ///< The hit density along the axis between the two vertices
    };

    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief Event feature info class
     */
    class EventFeatureInfo
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  eventShoweryness the event showeryness feature
         *  @param  eventEnergy the energy of the event
         *  @param  eventArea the area of the event
         *  @param  longitudinality the longitudinality of the event
         *  @param  nHits the number of hits in the event
         *  @param  nClusters the number of clusters in the event
         *  @param  nCandidates the total number of vertex candidates
         */
        EventFeatureInfo(const float eventShoweryness, const float eventEnergy, const float eventArea, const float longitudinality,
            const unsigned int nHits, const unsigned int nClusters, const unsigned int nCandidates);

        float m_eventShoweryness;   ///< The event showeryness feature
        float m_eventEnergy;        ///< The event energy
        float m_eventArea;          ///< The area of the event
        float m_longitudinality;    ///< The longitudinality of the event
        unsigned int m_nHits;       ///< The number of hits in the event
        unsigned int m_nClusters;   ///< The number of clusters in the event
        unsigned int m_nCandidates; ///< The total number of vertex candidates
    };

    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Default constructor
     */
    TrainedVertexSelectionAlgorithm();

protected:
    typedef std::pair<pandora::CartesianVector, pandora::CartesianVector> ClusterEndPoints;
    typedef std::map<const pandora::Cluster *const, ClusterEndPoints> ClusterEndPointsMap;
    typedef std::vector<LArMvaHelper::MvaFeatureVector> FeatureListVector;
    typedef std::vector<pandora::VertexVector> VectorOfVertexVectors;

    /**
     *  @brief  Get the vertex score list
     *
     *  @param  vertexVector the vector of vertices
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the hit kd tree for the U view
     *  @param  kdTreeV the hit kd tree for the V view
     *  @param  kdTreeW the hit kd tree for the W view
     *  @param  vertexScoreList the vertex score list to fill
     */
    virtual void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const = 0;

    /**
     *  @brief  Calculate the shower cluster map for a cluster list
     *
     *  @param  inputClusterList the input cluster list
     *  @param  showerClusterList the shower cluster list to populate
     */
    void CalculateShowerClusterList(const pandora::ClusterList &inputClusterList, ShowerClusterList &showerClusterList) const;

    /**
     *  @brief  Add the endpoints of any shower-like clusters to the map
     *
     *  @param  clusterList the list of clusters
     *  @param  showerLikeClusters the list of shower-like clusters to populate
     *  @param  clusterEndPointsMap the map of shower-like cluster endpoints to populate
     */
    void GetShowerLikeClusterEndPoints(
        const pandora::ClusterList &clusterList, pandora::ClusterList &showerLikeClusters, ClusterEndPointsMap &clusterEndPointsMap) const;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::unordered_map<const pandora::CaloHit *, const pandora::Cluster *> HitToClusterMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::ClusterSet> NearbyClusterCache;

    /**
     * @brief   Populate kd tree with information about hits in a provided list of clusters
     *
     * @param   clusterList the list of clusters
     * @param   kdTree to receive the populated kd tree
     * @param   hitToClusterMap to receive the populated hit to cluster map
     */
    void PopulateKdTree(const pandora::ClusterList &clusterList, HitKDTree2D &kdTree, HitToClusterMap &hitToClusterMap) const;

    /**
     *  @brief  Try to add an available cluster to a given shower cluster, using shower clustering approximation
     *
     *  @param  clusterEndPointsMap the map of shower-like cluster endpoints
     *  @param  availableShowerLikeClusters the list of shower-like clusters still available
     *  @param  pCluster the cluster in the shower cluster from which to consider distances
     *  @param  showerCluster the shower cluster
     *
     *  @return boolean
     */
    bool AddClusterToShower(const ClusterEndPointsMap &clusterEndPointsMap, pandora::ClusterList &availableShowerLikeClusters,
        const pandora::Cluster *const pCluster, pandora::ClusterList &showerCluster) const;

    /**
     *  @brief  Try to add an available cluster to a given shower cluster, using cluster hit positions cached in kd tree
     *
     *  @param  kdTree the kd tree, used purely for efficiency in events with large hit multiplicity
     *  @param  hitToClusterMap the hit to cluster map, used to interpret kd tree findings
     *  @param  nearbyClusterCache per-cluster cache of nearby clusters (the kd tree and hit map are static, so this never goes stale)
     *  @param  availableShowerLikeClusters the list of shower-like clusters still available
     *  @param  pCluster the cluster in the shower cluster from which to consider distances
     *  @param  showerCluster the shower cluster
     *
     *  @return boolean
     */
    bool AddClusterToShower(HitKDTree2D &kdTree, const HitToClusterMap &hitToClusterMap, NearbyClusterCache &nearbyClusterCache,
        pandora::ClusterList &availableShowerLikeClusters, const pandora::Cluster *const pCluster, pandora::ClusterList &showerCluster) const;

    /**
     *  @brief  Calculate the event parameters
     *
     *  @param  clusterListU the U-view cluster list
     *  @param  clusterListV the V-view cluster list
     *  @param  clusterListW the W-view cluster list
     *  @param  vertexVector the vector of vertex candidates
     *
     *  @return the EventFeatureInfo object
     */
    EventFeatureInfo CalculateEventFeatures(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW, const pandora::VertexVector &vertexVector) const;

    /**
     *  @brief  Increment the showery hit parameters for a cluster list
     *
     *  @param  clusterList the cluster list
     *  @param  nShoweryHits the number of showery hits
     *  @param  nHits the number of hits
     *  @param  eventEnergy the event energy
     */
    void IncrementShoweryParameters(const pandora::ClusterList &clusterList, unsigned int &nShoweryHits, unsigned int &nHits, float &eventEnergy) const;

    /**
     *  @brief  Find whether a cluster is shower-like
     *
     *  @param  pCluster the cluster
     *
     *  @return boolean
     */
    bool IsClusterShowerLike(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the event shape features
     *
     *  @param  clusterList the cluster list
     *  @param  eventVolume the event volume
     *  @param  longitudinality the event longitudinality
     */
    void GetLegacyEventShapeFeatures(const pandora::ClusterList &clusterList, float &eventVolume, float &longitudinality) const;

    /**
     *  @brief  Get the event shape features
     *
     *  @param  clusterList the map of cluster lists for each view
     *  @param  eventArea the event area
     *  @param  longitudinality the event longitudinality
     */
    void GetEventShapeFeatures(const ClusterListMap &clusterListMap, float &eventArea, float &longitudinality) const;

    /**
     *  @brief  Get the coordinate span in one view
     *
     *  @param  clusterList the cluster list in one view
     *  @param  xSpan the coordinate span in the drift direction
     *  @param  zSpan the coordinate span in the wire direction
     */
    void Get2DSpan(const pandora::ClusterList &clusterList, float &xSpan, float &zSpan) const;

    /**
     *  @brief  Update the min/max coordinate spans
     *
     *  @param  minPositionCoord the min position coordinate
     *  @param  maxPositionCoord the max position coordinate
     *  @param  minCoord the current min coordinate
     *  @param  maxCoord the current max coordinate
     */
    void UpdateSpanCoordinate(
        const float minPositionCoord, const float maxPositionCoord, pandora::InputFloat &minCoord, pandora::InputFloat &maxCoord) const;

    /**
     *  @brief  Get the coordinate span
     *
     *  @param  minCoord the min coordinate
     *  @param  maxCoord the max coordinate
     *
     *  @return the coordinate span
     */
    float GetCoordinateSpan(const pandora::InputFloat &minCoord, const pandora::InputFloat &maxCoord) const;

    /**
     *  @brief  Add the event features to a vector in the correct order
     *
     *  @param  eventFeatureInfo the event feature info
     *  @param  featureVector the vector of doubles to append
     */
    void AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo, LArMvaHelper::MvaFeatureVector &featureVector) const;

    /**
     *  @brief  Populate the vertex feature info map for a given vertex
     *
     *  @param  beamConstants the beam constants
     *  @param  clusterListMap the cluster list map
     *  @param  slidingFitDataListMap the sliding fit data list map
     *  @param  showerClusterListMap the shower cluster list map
     *  @param  kdTreeMap the kd tree map
     *  @param  pVertex the vertex
     *  @param  vertexFeatureInfoMap the map to populate
     */
    void PopulateVertexFeatureInfoMap(const BeamConstants &beamConstants, const ClusterListMap &clusterListMap,
        const SlidingFitDataListMap &slidingFitDataListMap, const ShowerClusterListMap &showerClusterListMap, const KDTreeMap &kdTreeMap,
        const pandora::Vertex *const pVertex, VertexFeatureInfoMap &vertexFeatureInfoMap) const;

    /**
     *  @brief  Populate the initial vertex score list for a given vertex
     *
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pVertex the vertex
     *  @param  initialScoreList the score list to populate
     */
    void PopulateInitialScoreList(VertexFeatureInfoMap &vertexFeatureInfoMap, const pandora::Vertex *const pVertex, VertexScoreList &initialScoreList) const;

    /**
     *  @brief  Get the list of top-N separated vertices
     *
     *  @param  initialScoreList the initial score list
     *  @param  bestRegionVertices the list of best region vertices populate
     */
    void GetBestRegionVertices(VertexScoreList &initialScoreList, pandora::VertexVector &bestRegionVertices) const;

    /**
     *  @brief  Produce the region and vertex training sets
     *
     *  @param  vertexVector the vector of all vertices
     *  @param  bestRegionVertices the best region vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  eventFeatureList the list of event features
     *  @param  kdTreeMap
     */
    void ProduceTrainingSets(const pandora::VertexVector &vertexVector, const pandora::VertexVector &bestRegionVertices,
        VertexFeatureInfoMap &vertexFeatureInfoMap, const LArMvaHelper::MvaFeatureVector &eventFeatureList, const KDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Calculate the r/phi scores for the vertices in a vector, possibly erasing those that fail the fast score test
     *
     *  @param  vertexVector the vector of vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  kdTreeMap the kd tree map
     */
    void CalculateRPhiScores(pandora::VertexVector &vertexVector, VertexFeatureInfoMap &vertexFeatureInfoMap, const KDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Get the interaction type string
     *
     *  @return the interaction type string
     */
    std::string GetInteractionType() const;

    /**
     *  @brief  Produce a set of training examples for a binary classifier
     *
     *  @param  vertexVector the vector of vertices to use
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  coinFlip a distribution for producing coin flips
     *  @param  generator the random number generator
     *  @param  interactionType the interaction type string
     *  @param  trainingOutputFile the training set output file
     *  @param  eventFeatureList the event feature list
     *  @param  kdTreeMap the map of 2D hit kd trees
     *  @param  maxRadius the maximum allowed radius for the 'best' vertex
     *  @param  useRPhi whether to include the r/phi feature
     *
     *  @return address of the vertex used as the 'best' vertex in the classifier
     */
    const pandora::Vertex *ProduceTrainingExamples(const pandora::VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        std::bernoulli_distribution &coinFlip, std::mt19937 &generator, const std::string &interactionType, const std::string &trainingOutputFile,
        const LArMvaHelper::MvaFeatureVector &eventFeatureList, const KDTreeMap &kdTreeMap, const float maxRadius, const bool useRPhi) const;

    /**
     *  @brief  Use the MC information to get the best vertex from a list
     *
     *  @param  vertexVector the vector of vertices
     *  @param  pBestVertex address of the best vertex
     *  @param  bestVertexDr dR of the best vertex
     */
    void GetBestVertex(const pandora::VertexVector &vertexVector, const pandora::Vertex *&pBestVertex, float &bestVertexDr) const;

    /**
     *  @brief  Calculates the shared features of a pair of vertex candidates
     *
     *  @param  pVertex1 the address of the first vertex
     *  @param  pVertex2 the address of the second vertex
     *  @param  kdTreeMap the map of 2D hit kd trees
     *  @param  separation the 3D distance between the two vertex candidates
     *  @param  axisHits the number of hits between the two candidates normalised to the distance between them
     */
    void GetSharedFeatures(const pandora::Vertex *const pVertex1, const pandora::Vertex *const pVertex2, const KDTreeMap &kdTreeMap,
        float &separation, float &axisHits) const;

    /**
     *  @brief  Increments the axis hits information for one view
     *
     *  @param  pos1 2D projected position of the first vertex
     *  @param  pos2 2D projected position of the second vertex
     *  @param  kdTree the kd tree of 2D hits
     *  @param  axisHits the number of hits between the two candidates
     */
    void IncrementSharedAxisValues(const pandora::CartesianVector pos1, const pandora::CartesianVector pos2, HitKDTree2D &kdTree, float &axisHits) const;

    /**
     *  @brief  Determines whether a hit lies within the box defined by four other positions
     *
     *  @param  hitPos the hit position
     *  @param  point1 the first corner of the box
     *  @param  point2 the second corner of the box
     *  @param  point3 the third corner of the box
     *  @param  point4 the fourth corner of the box
     *
     *  @return boolean
     */
    bool IsHitInBox(const pandora::CartesianVector &hitPos, const pandora::CartesianVector &point1, const pandora::CartesianVector &point2,
        const pandora::CartesianVector &point3, const pandora::CartesianVector &point4) const;

    /**
     *  @brief  Add the vertex features to a vector in the correct order
     *
     *  @param  vertexFeatureInfo the vertex feature info
     *  @param  featureVector the vector of floats to append
     *  @param  useRPhi whether to include the r/phi feature
     */
    void AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo, LArMvaHelper::MvaFeatureVector &featureVector, const bool useRPhi) const;

    /**
     *  @brief  Add the shared features to a vector in the correct order
     *
     *  @param  vertexSharedFeatureInfo the shared vertex feature info
     *  @param  featureVector the vector of floats to append
     */
    void AddSharedFeaturesToVector(const VertexSharedFeatureInfo &vertexSharedFeatureInfo, LArMvaHelper::MvaFeatureVector &featureVector) const;

    /**
     *  @brief  Populate the final vertex score list using the r/phi score to find the best vertex in the vicinity
     *
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pFavouriteVertex address of the favourite vertex
     *  @param  vertexVector the vector of all vertex candidates
     *  @param  finalVertexScoreList the final vertex score list to populate
     */
    void PopulateFinalVertexScoreList(const VertexFeatureInfoMap &vertexFeatureInfoMap, const pandora::Vertex *const pFavouriteVertex,
        const pandora::VertexVector &vertexVector, VertexScoreList &finalVertexScoreList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    VertexFeatureTool::FeatureToolVector m_featureToolVector; ///< The feature tool vector
    bool m_trainingSetMode;                                   ///< Whether to train
    bool m_allowClassifyDuringTraining;                       ///< Whether classification is allowed during training
    float m_mcVertexXCorrection;                              ///< The correction to the x-coordinate of the MC vertex position
    std::string m_trainingOutputFileRegion;                   ///< The training output file for the region mva
    std::string m_trainingOutputFileVertex;                   ///< The training output file for the vertex mva
    std::string m_mcParticleListName;                         ///< The MC particle list for creating training examples
    std::string m_caloHitListName;                            ///< The 2D CaloHit list name

    pandora::StringVector m_inputClusterListNames; ///< The list of cluster list names
    unsigned int m_minClusterCaloHits;             ///< The min number of hits parameter in the energy score
    unsigned int m_slidingFitWindow;               ///< The layer window for the sliding linear fits
    float m_minShowerSpineLength;                  ///< The minimum length at which all are considered to be tracks

    float m_beamDeweightingConstant; ///< The beam deweighting constant for the initial region score list
    float m_localAsymmetryConstant;  ///< The local asymmetry constant for the initial region score list
    float m_globalAsymmetryConstant; ///< The global asymmetry constant for the initial region score list
    float m_showerAsymmetryConstant; ///< The shower asymmetry constant for the initial region score list
    float m_energyKickConstant;      ///< The energy kick constant for the initial region score list

    float m_showerClusteringDistance;         ///< The shower clustering distance
    unsigned int m_minShowerClusterHits;      ///< The minimum number of shower cluster hits
    bool m_useShowerClusteringApproximation;  ///< Whether to use the shower clustering distance approximation
    float m_regionRadius;                     ///< The radius for a vertex region
    float m_rPhiFineTuningRadius;             ///< The maximum distance the r/phi tune can move a vertex
    float m_maxTrueVertexRadius;              ///< The maximum distance at which a vertex candidate can be considered the 'true' vertex
    bool m_useRPhiFeatureForRegion;           ///< Whether to use the r/phi feature for the region vertex
    bool m_dropFailedRPhiFastScoreCandidates; ///< Whether to drop candidates that fail the r/phi fast score test
    bool m_testBeamMode;                      ///< Test beam mode
    bool m_legacyEventShapes;                 ///< Whether to use the old event shapes calculation
    bool m_legacyVariables;                   ///< Whether to only use the old variables
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrainedVertexSelectionAlgorithm::VertexFeatureInfo::VertexFeatureInfo(const float beamDeweighting, const float rPhiFeature, const float energyKick,
    const float localAsymmetry, const float globalAsymmetry, const float showerAsymmetry, const float dEdxAsymmetry, const float vertexEnergy) :
    m_beamDeweighting(beamDeweighting),
    m_rPhiFeature(rPhiFeature),
    m_energyKick(energyKick),
    m_localAsymmetry(localAsymmetry),
    m_globalAsymmetry(globalAsymmetry),
    m_showerAsymmetry(showerAsymmetry),
    m_dEdxAsymmetry(dEdxAsymmetry),
    m_vertexEnergy(vertexEnergy)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrainedVertexSelectionAlgorithm::EventFeatureInfo::EventFeatureInfo(const float eventShoweryness, const float eventEnergy,
    const float eventArea, const float longitudinality, const unsigned int nHits, const unsigned int nClusters, const unsigned int nCandidates) :
    m_eventShoweryness(eventShoweryness),
    m_eventEnergy(eventEnergy),
    m_eventArea(eventArea),
    m_longitudinality(longitudinality),
    m_nHits(nHits),
    m_nClusters(nClusters),
    m_nCandidates(nCandidates)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrainedVertexSelectionAlgorithm::VertexSharedFeatureInfo::VertexSharedFeatureInfo(const float separation, const float axisHits) :
    m_separation(separation),
    m_axisHits(axisHits)
{
}

} // namespace lar_content

#endif // #ifndef LAR_TRAINED_VERTEX_SELECTION_ALGORITHM_H
