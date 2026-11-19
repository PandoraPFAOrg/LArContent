/**
 *  @file   larpandoradlcontent/LArThreeDVertex/DLVertexingThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning 3D vertexing algorithm.
 *
 *  $Log: $
 */

#include <string>
#include <vector>

#include "Api/PandoraApi.h"
#include "Objects/CartesianVector.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArThreeDVertex/DlVertexingThreeDAlgorithm.h"
#include "larpandoradlcontent/LArThreeDVertex/HoughFinder.h"

#include <Eigen/Dense>
#include <c10/core/TensorOptions.h>
#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlThreeDVertexingAlgorithm::DlThreeDVertexingAlgorithm() :
    m_pass{1},
    m_scalingFactor{-1.0f},
    m_thresholds{},
    m_nDistanceClasses{-1}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::Run()
{
    if (m_pass == 1)
        return this->RunPass1();
    else if (m_pass == 2)
        return this->RunPass2();
    else
    {
        std::cout << "DlThreeDVertexingAlgorithm::Run - Invalid pass number: " << m_pass << ". Must be 1 or 2." << std::endl;
        return STATUS_CODE_FAILURE;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::RunPass1()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const unsigned int numHits = pCaloHitList->size();

    if (numHits == 0)
    {
        std::cout << "DlThreeDVertexingAlgorithm::RunPass1 - No calo hits in the input list, skipping." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    // ATTN: Load the Slicing Worker instance, to pull out the slice PFOs + their associated candidate vertices.
    //       We can use them as a start point for where to zoom in on.
    const auto parentInstance = MultiPandoraApi::GetPrimaryPandoraInstance(&(*this).GetPandora());
    const auto childInstanceList = MultiPandoraApi::GetDaughterPandoraInstanceList(parentInstance);
    const pandora::Pandora *pSlicingInstance{nullptr};

    for (const auto pInstance : childInstanceList)
    {
        if (pInstance->GetName() == "SlicingWorker")
        {
            pSlicingInstance = pInstance;
            break;
        }
    }

    if (pSlicingInstance != nullptr)
    {
        const PfoList *pSlicePfos(nullptr);
        const auto pfoResult = PandoraApi::GetCurrentPfoList(*pSlicingInstance, pSlicePfos);

        if (pfoResult != STATUS_CODE_SUCCESS || pSlicePfos == nullptr)
        {
            std::cout << "DlVertexingThreeDAlgorithm: Failed to get PFO list from slicing instance." << std::endl;
            return STATUS_CODE_FAILURE;
        }

        CartesianPointVector candidateVertices({});
        for (const auto pPfo : *pSlicePfos)
        {
            const auto vertices = pPfo->GetVertexList();
            for (const auto pVertex : vertices)
                candidateVertices.push_back(pVertex->GetPosition());
        }

        std::vector<CartesianVector> matchedVertices;
        const float matchEpsilonSq = 1e-4f;

        for (const CartesianVector &vertexPos : candidateVertices)
        {
            bool foundMatch = false;
            for (const auto pCaloHit : *pCaloHitList)
            {
                if (nullptr == pCaloHit)
                    continue;

                if ((pCaloHit->GetPositionVector() - vertexPos).GetMagnitudeSquared() < matchEpsilonSq)
                {
                    foundMatch = true;
                    break;
                }
            }

            if (foundMatch)
                matchedVertices.push_back(vertexPos);
        }

        if (!matchedVertices.empty())
        {
            if (matchedVertices.size() > 1)
            {
                std::cout << "DlVertexingThreeDAlgorithm: WARNING - Found " << matchedVertices.size()
                          << " candidate vertices in this slice. Proceeding with the first one." << std::endl;
            }

            const VertexList *pVertexList{nullptr};
            std::string temporaryListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

            // TODO: First vertex? All vertices? Can we pick the best or something?
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = matchedVertices.front();
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;

            const Vertex *pVertex(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));

            return STATUS_CODE_SUCCESS;
        }
    }

    // INFO: If we are here...then we failed to load or find a vertex for this slice. Run zoomed out inference.
    std::cout << "DlVertexingThreeDAlgorithm: No candidate vertex found from slicing. Inferring to produce candidate vertex." << std::endl;
    return this->RunModel(*pCaloHitList, nullptr);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::RunPass2()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const unsigned int numHits = pCaloHitList->size();

    if (numHits == 0)
    {
        std::cout << "DlThreeDVertexingAlgorithm::RunPass2 - No calo hits in the input list, skipping." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    const VertexList *pInputVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList));

    if (!pInputVertexList || pInputVertexList->empty())
    {
        std::cout << "DlThreeDVertexingAlgorithm::RunPass2 - No valid input vertices found, skipping." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    const auto cropCenter = (*pInputVertexList->begin())->GetPosition();
    return this->RunModel(*pCaloHitList, &cropCenter);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::RunModel(const pandora::CaloHitList &caloHits, const pandora::CartesianVector *cropVertex)
{
    std::vector<CartesianVector> nodes;
    std::vector<std::array<float, 1>> node_features;
    this->GetNodeData(caloHits, nodes, node_features, cropVertex);
    const int numHits(nodes.size());

    if (numHits == 0)
    {
        std::cout << "DLVertexingThreeDAlgorithm::RunModel - No calo hits to build graph with, skipping." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    LArDLHelper::TorchInputVector inputs;
    this->BuildInput(inputs, nodes, node_features);

    node_features.clear();
    node_features.shrink_to_fit();

    std::vector<CartesianVector> foundVertices;

    {
        LArDLHelper::TorchMultiOutput semanticOutput;
        LArDLHelper::Forward(m_modelFile, inputs, semanticOutput);

        const auto semanticLabels = semanticOutput.toTuple()->elements()[0].toTensor();
        semanticOutput = LArDLHelper::TorchMultiOutput();

        {
            const auto contiguousSemanticLabels = semanticLabels.contiguous();
            std::vector<float> semanticLabelsVec(
                contiguousSemanticLabels.data_ptr<float>(), contiguousSemanticLabels.data_ptr<float>() + (numHits * m_nDistanceClasses));

            ThreeDVertexingHoughTransform houghFinder(m_thresholds, m_scalingFactor, 0.5f, 3, 10.f, 0, true);
            foundVertices = houghFinder.Fit(nodes, semanticLabelsVec);
        }

        // TODO: We do also have that noise mask here if useful to be stored in the event context etc?

        if (foundVertices.size() == 0)
        {
            // INFO: Rather than fail, if we had a rough region...its better than literally nothing.
            if (cropVertex != nullptr)
            {
                std::cout << "DLVertexingThreeDAlgorithm::Infer - Hough returned empty. Reverting to slicing crop vertex." << std::endl;
                foundVertices.push_back(*cropVertex);
            }
            else
            {
                std::cout << "DLVertexingThreeDAlgorithm::Infer - no vertex candidates found and no crop vertex! Skipping writing list." << std::endl;
                return STATUS_CODE_SUCCESS;
            }
        }
    }

    const VertexList *pVertexList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const CartesianVector &position : foundVertices)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::GetNodeData(const CaloHitList &caloHits, std::vector<CartesianVector> &nodes,
    std::vector<std::array<float, 1>> &node_features, const pandora::CartesianVector *cropVertex)
{
    nodes.reserve(caloHits.size());
    node_features.reserve(caloHits.size());

    const float cropRadiusSq = 30.0f * 30.0f;

    for (const auto pCaloHit : caloHits)
    {
        if (nullptr == pCaloHit)
            continue;

        if (cropVertex == nullptr || ((pCaloHit->GetPositionVector() - *cropVertex).GetMagnitudeSquared() <= cropRadiusSq))
        {
            nodes.push_back(pCaloHit->GetPositionVector());
            node_features.push_back({pCaloHit->GetInputEnergy() / 10.f});
        }
    }

    nodes.shrink_to_fit();
    node_features.shrink_to_fit();

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::BuildInput(
    LArDLHelper::TorchInputVector &inputs, std::vector<CartesianVector> &nodes, std::vector<std::array<float, 1>> &node_features)
{
    const int numNodes{static_cast<int>(nodes.size())};
    const int numFeatures{static_cast<int>(node_features[0].size())};

    LArDLHelper::TorchInput posTensor, xTensor;
    posTensor = torch::empty({numNodes, 3}, torch::kFloat32);
    xTensor = torch::empty({numNodes, numFeatures}, torch::kFloat32);

    torch::Tensor batchTensor = torch::zeros(numNodes, torch::kLong);
    float *posTensorPtr = posTensor.data_ptr<float>();
    float *xTensorPtr = xTensor.data_ptr<float>();

    for (int i = 0; i < numNodes; ++i)
    {
        posTensorPtr[i * 3 + 0] = nodes[i].GetX();
        posTensorPtr[i * 3 + 1] = nodes[i].GetY();
        posTensorPtr[i * 3 + 2] = nodes[i].GetZ();
        xTensorPtr[i] = node_features[i][0];
    }

    inputs.insert(inputs.end(), {xTensor, posTensor, batchTensor});

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlThreeDVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    std::string modelName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
    modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
    LArDLHelper::LoadModel(modelName, m_modelFile);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Pass", m_pass));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    if (m_pass != 1 && m_pass != 2)
    {
        std::cout << "DlThreeDVertexingAlgorithm::ReadSettings - Invalid pass number: " << m_pass << ". Must be 1 or 2." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    if (m_pass == 2 && m_inputVertexListName.empty())
    {
        std::cout << "DlThreeDVertexingAlgorithm::ReadSettings - Pass 2 requires an input vertex list name, but it is missing." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));

    m_nDistanceClasses = m_thresholds.size() + 1; // We have one more class than thresholds, as the thresholds define the boundaries between classes.

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
