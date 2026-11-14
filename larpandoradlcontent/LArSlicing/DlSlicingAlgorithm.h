/**
 *  @file   larpandoradlcontent/LArSlicing/DLSlicingAlgorithm.h
 *
 *  @brief  Header file for the deep learning slicing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SLICING_ALGORITHM_H
#define LAR_DL_SLICING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"
#include "Objects/EventContext.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningSLicingAlgorithm class
 */
class DlSlicingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlSlicingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode Infer();

    /**
     *  @brief Get the node data (positions and features) from the input CaloHit list.
     *
     *  @param caloHits The input calo hits.
     *  @param pos The positions of the nodes.
     *  @param node_features The features for each node.
     */
    pandora::StatusCode GetNodeData(
        const pandora::CaloHitList &caloHits, std::vector<pandora::CartesianVector> &pos, std::vector<std::array<float, 1>> &node_features);

    /**
     *  @brief Process the given node data into the expected tensor format for the model, and insert into the input vector.
     *
     *  @param inputs The input vector to be filled with the node data.
     *  @param pos The positions of the nodes.
     *  @param node_features The features for each node.
     */
    pandora::StatusCode BuildInput(
        LArDLHelper::TorchInputVector &inputs, std::vector<pandora::CartesianVector> &pos, std::vector<std::array<float, 1>> &node_features);

    /**
     *  @brief Given the predicted slicing labels, perform topological
     *         splitting, to aide further post-processing and cluster classification.
     *         Clusters are split into either "anchor" clusters (large predicted
     *         clusters near to candidate vertices) or "debris" clusters (small
     *         predicted clusters far from candidate vertices).
     *
     *  @param positions The positions of the nodes.
     *  @param clusterLabels The predicted cluster labels for each node.
     *  @param candidateIndices The indices of the candidate vertices.
     *  @param newLabels The new cluster labels for each node, after topological splitting
     *  @param anchors The indices of the anchor nodes
     *  @param debris The indices of the debris nodes
     *  @param distanceThreshold The distance threshold to use when deciding whether to split a cluster.
     *  @param minAnchorSize The minimum number of nodes for an anchor cluster.
     */
    pandora::StatusCode SplitAndClassifyClusters(const std::vector<pandora::CartesianVector> &positions,
        const std::vector<int> &clusterLabels, const std::vector<int> &candidateIndices, std::vector<int> &newLabels,
        std::set<int> &anchors, std::set<int> &debris, const float distanceThreshold = 20.f, const int minAnchorSize = 20) const;

    /**
     *  @brief Starting from split anchors, perform a "flood fill" type
     *         algorithm to start attaching debris to the anchors, improving
     *         completeness whilst maintaining good purity.
     *         Bonuses are applied for attaching debris back to its original
     *         predicted cluster, similarly can gate connections with optional
     *         t0 information.
     *
     *  @param positions The positions of the nodes.
     *  @param t0s The t0 values of the nodes, if available.
     *  @param t0Valid The validity flags for each t0 value (i.e. missing / unset t0s will be marked as invalid).
     *  @param originalLabels The original predicted cluster labels for each node.
     *  @param finalLabels The new cluster labels for each node, after flood fill.
     *  @param anchors The indices of the anchor nodes.
     *  @param debris The indices of the debris nodes.
     *  @param baseGap The base distance gap to use when deciding whether to attach debris to an anchor.
     *  @param ogBonusGap The distance gap to use when applying a bonus for attaching debris back to its original predicted cluster.
     */
    pandora::StatusCode FloodFill(const std::vector<pandora::CartesianVector> &positions, const std::vector<float> &t0s,
        const std::vector<bool> &t0Valid, const std::vector<int> &originalLabels, std::vector<int> &finalLabels,
        const std::set<int> &anchors, const std::set<int> &debris, const float baseGap = 3.f, const float ogBonusGap = 15.f) const;

    /**
     *  @brief Final clean up algorithm to classify any remaining small
     *         clusters. Clusters are either: Restored back to their original,
     *         pre-split cluster, promoted up to a real cluster, or attached to the
     *         nearest large cluster.
     *
     *  @param positions The positions of the nodes.
     *  @param t0s The t0 values of the nodes, if available.
     *  @param t0Valid The validity flags for each t0 value (i.e. missing / unset t0s will be marked as invalid).
     *  @param originalLabels The original predicted cluster labels for each node.
     *  @param finalLabels The new cluster labels for each node, after clean up.
     *  @param minSize The min hit size for a cluster to be considered real, rather than debris.
     */
    pandora::StatusCode CleanSmallClusters(const std::vector<pandora::CartesianVector> &positions, const std::vector<float> &t0s,
        const std::vector<bool> &t0Valid, const std::vector<int> &originalLabels, std::vector<int> &finalLabels, const int minSize = 300) const;

    /**
     *  @brief Precompute the median t0 values for each cluster.
     *
     *  @param labels The cluster labels for each node.
     *  @param t0s The t0 values for each node.
     *  @param t0Valid The validity flags for each t0 value (i.e. missing / unset t0s will be marked as invalid).
     *
     *  @return A map from cluster IDs to their median t0 values.
     */
    std::unordered_map<int, float> PrecomputeT0Medians(
        const std::vector<int> &labels, const std::vector<float> &t0s, const std::vector<bool> &t0Valid) const;

    LArDLHelper::TorchModel m_modelFile; ///< The model to use.

    float m_scalingFactor;               ///< The scaling factor for the input.
    std::vector<float> m_thresholds;     ///< Distance Class Thresholds.
    int m_nDistanceClasses;              ///< The number of distance classes (derived from thresholds).
    std::string m_caloHitListName;       ///< The name of the input CaloHit list.
    std::string m_outputClusterListName; ///< The name of the output Cluster list to create with the predicted instances.
    std::string m_outputVertexListName;  ///< The name of the output Vertex list to create with the predicted candidate vertices.
    bool m_runPostProcessing;            ///< Whether to run the post-processing steps after inference.
};

} // namespace lar_dl_content

#endif // LAR_DL_SLICING_ALGORITHM_H
