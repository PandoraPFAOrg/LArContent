/**
 *  @file   larpandoradlcontent/LArSlicing/DLSlicingAlgorithm.h
 *
 *  @brief  Header file for the deep learning slicing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SLICING_ALGORITHM_H
#define LAR_DL_SLICING_ALGORITHM_H 1

#include "Objects/EventContext.h"
#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

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
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /**
     *  @brief Get the node data (positions and features) from the input CaloHit list.
     *
     *  @param caloHits The input calo hits.
     *  @param pos The positions of the nodes.
     *  @param nodeFeatures The features for each node.
     */
    pandora::StatusCode GetNodeData(
        const pandora::CaloHitList &caloHits, std::vector<pandora::CartesianVector> &pos, std::vector<std::array<float, 1>> &nodeFeatures);

    /**
     *  @brief Process the given node data into the expected tensor format for the model, and insert into the input vector.
     *
     *  @param inputs The input vector to be filled with the node data.
     *  @param pos The positions of the nodes.
     *  @param nodeFeatures The features for each node.
     *  @param xTensor The input tensor for the model.
     */
    pandora::StatusCode BuildInput(LArDLHelper::TorchInputVector &inputs, std::vector<pandora::CartesianVector> &pos,
        std::vector<std::array<float, 1>> &nodeFeatures, LArDLHelper::TorchInput &xTensor);

    /**
     *  @brief Extract candidate vertices from the model output, using NMS.
     *
     *  @param pos The positions of the nodes.
     *  @param semLogits The vertex logits for each node (shape [N, 1]).
     *  @param candidateHitIdx The indices of the hits chosen as candidate vertices, highest score first (cleared on entry).
     */
    pandora::StatusCode ExtractCandidateVertices(
        const std::vector<pandora::CartesianVector> &pos, const torch::Tensor &semLogits, std::vector<int> &candidateHitIdx) const;

    LArDLHelper::TorchModel m_modelFile; ///< The model to use.

    std::string m_caloHitListName;       ///< The name of the input CaloHit list.
    std::string m_outputClusterListName; ///< The name of the output Cluster list to create with the predicted instances.
    std::string m_outputVertexListName;  ///< The name of the output Vertex list to create with the predicted candidate vertices.
    float m_vertexScoreThreshold;        ///< The threshold for the vertex score to consider a candidate vertex as valid.
    float m_nmsRadius;                   ///< The radius for non-maximum suppression.

    bool m_trainingMode;            ///< Write the hits and their truth for training instead of running the model.
    std::string m_trainingFileName; ///< The CSV file the training examples are appended to.
};

} // namespace lar_dl_content

#endif // LAR_DL_SLICING_ALGORITHM_H
