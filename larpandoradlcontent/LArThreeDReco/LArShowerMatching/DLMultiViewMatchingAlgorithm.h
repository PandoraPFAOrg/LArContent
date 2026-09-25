/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/DLMultiViewMatchingAlgorithm.h
 *
 *  @brief  Header file for the dl multi view matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MULTI_VIEW_MATCHING_ALGORITHM_H
#define LAR_MULTI_VIEW_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLShowerHelper.h"

namespace lar_dl_content
{

class DLShowerMatchingTool;    

/**
 *  @brief  DLMultiViewMatchingAlgorithm class
 */
class DLMultiViewMatchingAlgorithm : public pandora::Algorithm
{
public:
    struct ClusterGroup
    {
        pandora::ClusterList m_clustersU; ///< clusters in the U-view
        pandora::ClusterList m_clustersV; ///< clusters in the V-view
        pandora::ClusterList m_clustersW; ///< clusters in the W-view        
    };
    typedef std::vector<ClusterGroup> ClusterGroupVector;
    
    typedef std::map<const pandora::Cluster*, std::map<const pandora::Cluster*, float>> SimilarityMatrix;      
        
    /**
     *  @brief  Default constructor
     */
    DLMultiViewMatchingAlgorithm();

    /**
     *  @brief  Default destructor
     */    
    ~DLMultiViewMatchingAlgorithm();    

    /**
     *  @brief  Obtain groups of 'connected' clusters by following cluster navigations
     *
     *  @param[out] clusterGroupVector the output vector of ClusterGroups
     */    
    void GetConnectedGroups(ClusterGroupVector &clusterGroupVector);
    
    /**
     *  @brief  Create a pfo from a list of input clusters
     *
     *  @param clusters the input list of clusters
     */
    void CreatePfo(const pandora::ClusterList &clusters);

    /**
     *  @brief  Delete a cluster from all algorithm containers
     *
     *  @param pClusterToRemove a pointer to the cluster to remove
     */    
    void DeleteCluster(const pandora::Cluster *const pClusterToRemove);    

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const pandora::Cluster *, pandora::ClusterList> NavigationMap;
    typedef std::map<const pandora::Cluster *, std::pair<float, float>> ClusterExtentMap;
    typedef std::vector<DLShowerMatchingTool *> MatchingToolVector;     
    
    /**
     *  @brief  Fill an algorithm list
     *
     *  @param listName the name of the internal pandora list     
     *  @param[out] pList the algorithm list to fill
     *
     *  @return StatusCode whether the list has been filled
     */
    template <typename T>
    pandora::StatusCode GetList(const std::string listName, const T *&pList);

    /**
     *  @brief  Filter input clusters and prepare algorithm containers
     *
     *  @param pClusterListU the input list of U-view clusters
     *  @param pClusterListV the input list of V-view clusters
     *  @param pClusterListW the input list of W-view clusters     
     *  @param[out] clusterExtentMap the map of clusters->extremal drift coordinates
     */
    void PrepareClusters(const pandora::ClusterList *const pClusterListU, const pandora::ClusterList *const pClusterListV,
        const pandora::ClusterList *const pClusterListW, ClusterExtentMap &clusterExtentMap);

    /**
     *  @brief  Fill the algorithm navigation maps
     *
     *  @param clusterExtentMap the map of clusters->extremal drift coordinates     
     */
    void FillNavigationMaps(const ClusterExtentMap &clusterExtentMap);

    /**
     *  @brief  Whether the input cluster pair overlap (at all) in the wire and drift plane
     *
     *  @param pCluster1 the first cluster in the cluster pair match
     *  @param pCluster2 the second cluster in the cluster pair match
     *  @param clusterExtentMap the map of clusters->extremal drift coordinates     
     */
    bool DoClustersOverlap(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const ClusterExtentMap &clusterExtentMap);

    /**
     *  @brief  Whether the input cluster pair overlap (to a defined extent) in the wire plane
     *
     *  @param pCluster1 the first cluster in the cluster pair match
     *  @param pCluster2 the second cluster in the cluster pair match
     *  @param minX the lower bound of the drift overlap region
     *  @param maxX the upper bound of the drift overlap region     
     */    
    bool DoClustersOverlapInWire(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const float minX, const float maxX);

    /**
     *  @brief  Obtain group of 'connected' clusters by following the cluster navigations of an initial input cluster
     *
     *  @param pCluster the initial input cluster
     *  @param[out] clusterGroup the group of connected clusters
     *  @param[in,out] usedU the list of U-view clusters already assigned to a cluster group
     *  @param[in,out] usedV the list of V-view clusters already assigned to a cluster group
     *  @param[in,out] usedW the list of W-view clusters already assigned to a cluster group
     */    
    void GetConnectedGroup(const pandora::Cluster *const pCluster, ClusterGroup &clusterGroup,
        pandora::ClusterList &usedU, pandora::ClusterList &usedV, pandora::ClusterList &usedW);
    
    /**
     *  @brief  Apply the matching model to obtain pairwise similarity scores between all UVW clusters considered by the algorithm
     *
     *  @param clusterGroupVector the vector of ClusterGroups
     *  @param[out] globalSimMatrix the output cluster->cluster->score mapping
     */    
    void FillGlobalSimMatrix(const ClusterGroupVector &clusterGroupVector, SimilarityMatrix &globalSimMatrix);

    /**
     *  @brief  Apply the matching model to obtain pairwise similarity scores between UVW clusters within a common drift region
     *
     *  @param clusterListU the list of U clusters
     *  @param clusterListV the list of V clusters
     *  @param clusterListW the list of W clusters
     *  @param viewToVtxPos the hit type -> 2D nu vertex map
     *  @param detXGaps the container of drift-gaps within the detector
     *  @param[out] clusterSimMat the output cluster->cluster->score mapping
     */        
    pandora::StatusCode PredictClusterSimilarityMatrix(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW, const std::map<pandora::HitType, pandora::CartesianVector> &viewToVtxPos,
        const std::set<float> &detXGaps, SimilarityMatrix &clusterSimMat);

    /**
     *  @brief  Construct the tensor that encodes a cluster.
     *
     *  @param clusterFeatures vector of hit properties for the hits of one 2D cluster
     *  @param view the view of the 2D clusters
     *  @param nClusters the total number of 2D clusters in the view
     *  @param tensorCluster the tensor encoding the cluster as a sequence of hit feature vectors
     */    
    void MakeClusterTensor(const std::vector<LArDLShowerHelper::HitFeatures> &clusterFeatures, const pandora::HitType view,
        const int nClusters, torch::Tensor &tensorCluster) const;

    /**
     *  @brief  Populates a SimilarityMatrix from the Tensor output of model inference.
     *
     *  @param tensorSimMat the predicted similarity matrix tensor, shape [1, N, N] where N is the number of 2D clusters
     *  @param clusterList the list of 2D clusters
     *  @param[out] clusterSimMat the cluster similarity matrix object
     */    
    pandora::StatusCode PopulateClusterSimilarityMatrix(const torch::Tensor &tensorSimMat, const pandora::ClusterVector &clusterVector,
        SimilarityMatrix &clusterSimMat) const;

    /**
     *  @brief Update the algorithm's navigation maps, removing the links between clusters with poor similarity scores
     *
     *  @param globalSimMatrix the cluster->cluster->score mapping
     */    
    void UpdateNavigationMaps(const SimilarityMatrix &globalSimMatrix);

    /**
     *  @brief Reset algorithm containers
     */     
    void CleanUp();

    /**
     *  @brief Write, to a tree, the training dataset such that each entry in the tree corresponds to connected cluster group
     *
     *  @param clusterGroupVector the vector of ClusterGroups
     */       
    void PrepareTrainingSample(const ClusterGroupVector &clusterGroupVector);

    std::string m_nuVertexListName;   ///< The neutrino vertex list name
    std::string m_clusterListNameU;   ///< The input cluster list name for the U view
    std::string m_clusterListNameV;   ///< The input cluster list name for the V view
    std::string m_clusterListNameW;   ///< The input cluster list name for the W view
    std::string m_outputPfoListName;  ///< The name of the list in which to store the created pfos
    pandora::ClusterList m_filteredU; ///< The filtered list of U-view clusters to consider for matching
    pandora::ClusterList m_filteredV; ///< The filtered list of V-view clusters to consider for matching
    pandora::ClusterList m_filteredW; ///< The filtered list of W-view clusters to consider for matching
    NavigationMap m_navigationU; ///< The mapping between U-view clusters and 'matched' V/W clusters
    NavigationMap m_navigationV; ///< The mapping between V-view clusters and 'matched' U/W clusters
    NavigationMap m_navigationW; ///< The mapping between W-view clusters and 'matched' U/V clusters
    MatchingToolVector m_matchingToolVector; ///< The algorithm tool vector   
    bool m_trainingMode;               ///< Whether to run the algorithm in training mode
    std::string m_trainingFileName;    ///< The name of the produced training file
    std::string m_trainingTreeName;    ///< The name of the produced training tree
    unsigned int m_minNClusterHits;    ///< The threshold number of hits of a considered cluster
    unsigned int m_nMaxRepeats;        ///< The maximum number of times to repeat the matching tools
    float m_minWireOverlapFraction;    ///< The minimum required fraction of hits that exist on intersecting wires
    int m_hitFeatureDim;               ///< The number of hit features 
    float m_polarRScaleFactor;              ///< Scale factor for polar r coordinate input features
    float m_cartesianXScaleFactor;          ///< Scale factor for cartesian x coordinate input features
    float m_cartesianZScaleFactor;          ///< Scale factor for cartesian z coordinate input features
    float m_matchThreshold;                 ///< The threshold on the predicted similarity for a match to be made
    LArDLHelper::TorchModel m_modelEncoder; ///< TorchScript model for encoding hits in a cluster
    LArDLHelper::TorchModel m_modelAttn;    ///< TorchScript model for attention over encoded clusters in a view
    LArDLHelper::TorchModel m_modelSim;     ///< TorchScript model for pairwise similarities over clusters after attention    
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DLShowerTensorTool class
 */
class DLShowerMatchingTool : public pandora::AlgorithmTool
{
public:
   
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param pAlgorithm address of the calling algorithm
     *  @param globalSimMatrix the output cluster->cluster->score mapping
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix) = 0;
};    
    
} // namespace lar_dl_content

#endif // #ifndef LAR_MULTI_VIEW_MATCHING_ALGORITHM_H

