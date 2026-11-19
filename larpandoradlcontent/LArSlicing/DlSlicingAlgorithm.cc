/**
 *  @file   larpandoradlcontent/LArSlicing/DlSlicingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning slicing algorithm.
 *
 *  $Log: $
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "Api/PandoraApi.h"
#include "Objects/CartesianVector.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArSlicing/DlSlicingAlgorithm.h"

#include <Eigen/Dense>
#include <c10/core/TensorOptions.h>
#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSlicingAlgorithm::DlSlicingAlgorithm() :
    m_vertexScoreThreshold{0.5f},
    m_nmsRadius{8.0f},
    m_trainingMode{false}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::PrepareTrainingSample()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    const auto detectorBoundaries(LArGeometryHelper::GetDetectorBoundaries(this->GetPandora()));

    std::map<const MCParticle *, int> particleIds, neutrinoIds;
    const auto getId = [](std::map<const MCParticle *, int> &ids, const MCParticle *pMCParticle)
    { return ids.try_emplace(pMCParticle, static_cast<int>(ids.size())).first->second; };

    constexpr int missing{-999};
    constexpr int nValuesPerHit{14};
    const double nan{std::numeric_limits<double>::quiet_NaN()};

    const double run{static_cast<double>(this->GetPandora().GetRun())};
    const double subrun{static_cast<double>(this->GetPandora().GetSubrun())};
    const double event{static_cast<double>(this->GetPandora().GetEvent())};

    LArMvaHelper::MvaFeatureVector featureVector;
    featureVector.reserve(5 + nValuesPerHit * pCaloHitList->size());
    featureVector.emplace_back(run);
    featureVector.emplace_back(subrun);
    featureVector.emplace_back(event);
    featureVector.emplace_back(static_cast<double>(pCaloHitList->size()));
    featureVector.emplace_back(static_cast<double>(nValuesPerHit));

    unsigned int nHitsWithNeutrino{0};

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const CartesianVector &position{pCaloHit->GetPositionVector()};
        featureVector.emplace_back(static_cast<double>(position.GetX()));
        featureVector.emplace_back(static_cast<double>(position.GetY()));
        featureVector.emplace_back(static_cast<double>(position.GetZ()));
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetInputEnergy()));

        const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        featureVector.emplace_back(pLArCaloHit ? static_cast<double>(pLArCaloHit->GetLArTPCVolumeId()) : static_cast<double>(missing));

        const MCParticle *pDominant{nullptr};
        float dominantWeight{0.f}, totalWeight{0.f};
        if (pLArCaloHit)
        {
            for (const auto &[pMCParticle, weight] : pLArCaloHit->GetMCParticleWeightMap())
            {
                totalWeight += weight;
                if (weight > dominantWeight)
                {
                    dominantWeight = weight;
                    pDominant = pMCParticle;
                }
            }
        }

        featureVector.emplace_back(pDominant ? static_cast<double>(pDominant->GetParticleId()) : static_cast<double>(missing));
        featureVector.emplace_back(pDominant ? static_cast<double>(getId(particleIds, pDominant)) : static_cast<double>(missing));
        featureVector.emplace_back(totalWeight > 0.f ? static_cast<double>(dominantWeight / totalWeight) : 0.);

        const MCParticle *pNeutrino{nullptr};
        if (pDominant)
        {
            try
            {
                const MCParticle *const pParent{LArMCParticleHelper::GetParentMCParticle(pDominant)};
                if (LArMCParticleHelper::IsNeutrino(pParent))
                    pNeutrino = pParent;
            }
            catch (const StatusCodeException &)
            {
            }
        }

        if (pNeutrino)
        {
            const CartesianVector &vertex{pNeutrino->GetVertex()};
            featureVector.emplace_back(static_cast<double>(getId(neutrinoIds, pNeutrino)));
            featureVector.emplace_back(static_cast<double>(pNeutrino->GetParticleId()));
            featureVector.emplace_back(static_cast<double>(vertex.GetX()));
            featureVector.emplace_back(static_cast<double>(vertex.GetY()));
            featureVector.emplace_back(static_cast<double>(vertex.GetZ()));
            featureVector.emplace_back(LArGeometryHelper::IsInDetector(detectorBoundaries, vertex) ? 1. : 0.);
            ++nHitsWithNeutrino;
        }
        else
        {
            featureVector.emplace_back(static_cast<double>(missing));
            featureVector.emplace_back(static_cast<double>(missing));
            featureVector.emplace_back(nan);
            featureVector.emplace_back(nan);
            featureVector.emplace_back(nan);
            featureVector.emplace_back(0.);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_trainingFileName, true, featureVector));

    std::cout << "DlSlicingAlgorithm::PrepareTrainingSample - run " << run << ", subrun " << subrun << ", event " << event << ": wrote "
              << pCaloHitList->size() << " hits, " << nHitsWithNeutrino << " with a neutrino parent (" << neutrinoIds.size()
              << " neutrinos)" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Infer()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    std::vector<CartesianVector> nodes;
    std::vector<std::array<float, 1>> nodeFeatures;
    this->GetNodeData(*pCaloHitList, nodes, nodeFeatures);
    const unsigned int numHits = nodes.size();

    if (numHits == 0)
        return STATUS_CODE_SUCCESS;

    LArDLHelper::TorchInputVector inputs;
    LArDLHelper::TorchInput xTensor;
    this->BuildInput(inputs, nodes, nodeFeatures, xTensor);

    nodeFeatures.clear();
    nodeFeatures.shrink_to_fit();

    torch::Tensor instancePreds;

    std::vector<CartesianVector> foundVertices;
    std::vector<int> candidateIndices;

    {
        LArDLHelper::TorchMultiOutput semanticOutput;
        LArDLHelper::Forward(m_modelFile, inputs, semanticOutput);

        auto vertexLogits = semanticOutput.toTuple()->elements()[0].toTensor();
        auto hitEmbeddings = semanticOutput.toTuple()->elements()[1].toTensor();
        auto encodedPos = semanticOutput.toTuple()->elements()[2].toTensor();
        semanticOutput = LArDLHelper::TorchMultiOutput();

        if (this->ExtractCandidateVertices(nodes, vertexLogits, candidateIndices) != STATUS_CODE_SUCCESS)
            return STATUS_CODE_FAILURE;

        for (const int idx : candidateIndices)
            foundVertices.push_back(nodes[idx]);

        const int numCandidates = foundVertices.size();
        nodes.clear();
        nodes.shrink_to_fit();

        if (numCandidates == 0)
        {
            std::cout << "DLSlicingAlgorithm::Infer - no vertex candidates found, skipping instance segmentation step" << std::endl;
            return STATUS_CODE_SUCCESS;
        }

        constexpr int numCandidateColumns = 5;
        const auto asFloat = torch::TensorOptions().dtype(torch::kFloat32);
        LArDLHelper::TorchInput candidateTensor;
        LArDLHelper::InitialiseInput({numCandidates, numCandidateColumns}, candidateTensor, asFloat);
        float *candidateTensorData = candidateTensor.data_ptr<float>();

        const auto logits = vertexLogits.reshape({-1}).to(torch::kFloat32).contiguous();
        const float *logitData = logits.data_ptr<float>();

        for (int i = 0; i < numCandidates; ++i)
        {
            float *row = candidateTensorData + i * numCandidateColumns;
            row[0] = foundVertices[i].GetX();
            row[1] = foundVertices[i].GetY();
            row[2] = foundVertices[i].GetZ();

            row[3] = 1.f / (1.f + std::exp(-logitData[candidateIndices[i]]));

            float minDistance{std::numeric_limits<float>::max()};
            for (int j = 0; j < i; ++j)
                minDistance = std::min(minDistance, (foundVertices[i] - foundVertices[j]).GetMagnitude());
            row[4] = (i == 0) ? 3.f : std::min(minDistance / 100.f, 3.f);
        }

        // ATTN: Populate the final input tensor with per hit features, embeddings,
        // position encodings, vertex logits, and candidate vertex positions.
        inputs.clear();
        LArDLHelper::TorchInputVector fullInputTensor;
        fullInputTensor.push_back(std::move(xTensor));
        fullInputTensor.push_back(std::move(hitEmbeddings));
        fullInputTensor.push_back(std::move(encodedPos));
        fullInputTensor.push_back(std::move(vertexLogits));
        fullInputTensor.push_back(std::move(candidateTensor));

        LArDLHelper::TorchOutput instanceOutput;
        LArDLHelper::Forward(m_modelFile, fullInputTensor, instanceOutput, "predict_instances");
        fullInputTensor.clear();

        instancePreds = torch::argmax(instanceOutput, 1).contiguous();
        instanceOutput = at::Tensor();
    }

    std::vector<int> finalLabels(numHits, -1);
    const int64_t *predsData = instancePreds.data_ptr<int64_t>();

    for (unsigned int i = 0; i < numHits; ++i)
        finalLabels[i] = static_cast<int>(predsData[i]);

    std::map<int, std::list<const CaloHit *>> clusterHitsMap;
    int hitIdx = 0;
    for (const auto pCaloHit : *pCaloHitList)
    {
        if (nullptr == pCaloHit)
            continue;

        const auto clusterPrediction = finalLabels[hitIdx];
        if (clusterPrediction >= 0)
            clusterHitsMap[clusterPrediction].push_back(pCaloHit);

        ++hitIdx;
    }

    const ClusterList *pClusterList(nullptr);
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, temporaryListName));

    for (const auto &[clusterId, hits] : clusterHitsMap)
    {
        if (hits.empty())
            continue;

        PandoraContentApi::Cluster::Parameters clusterParameters;
        clusterParameters.m_caloHitList = hits;
        const Cluster *pCluster(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, clusterParameters, pCluster));
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputClusterListName));
    }

    const VertexList *pVertexList(nullptr);
    temporaryListName.clear();
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const auto &vertex : foundVertices)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = vertex;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::GetNodeData(const CaloHitList &caloHits, std::vector<CartesianVector> &pos, std::vector<std::array<float, 1>> &nodeFeatures)
{
    pos.reserve(caloHits.size());
    nodeFeatures.reserve(caloHits.size());

    for (const auto pCaloHit : caloHits)
    {
        if (nullptr == pCaloHit)
            continue;

        const CartesianVector &hitPos = pCaloHit->GetPositionVector();

        pos.push_back(hitPos);
        nodeFeatures.push_back({pCaloHit->GetInputEnergy() / 10.f});
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::BuildInput(LArDLHelper::TorchInputVector &inputs, std::vector<CartesianVector> &pos,
    std::vector<std::array<float, 1>> &nodeFeatures, LArDLHelper::TorchInput &xTensor)
{
    const int numNodes{static_cast<int>(pos.size())};
    const int numFeatures{static_cast<int>(nodeFeatures[0].size())};

    LArDLHelper::TorchInput posTensor;
    posTensor = torch::empty({numNodes, 3}, torch::kFloat32);
    xTensor = torch::empty({numNodes, numFeatures}, torch::kFloat32);

    // ATTN: Fake batch tensor (inference only has a single batch)
    torch::Tensor batchTensor = torch::zeros(numNodes, torch::kLong);

    float *posTensorPtr = posTensor.data_ptr<float>();
    float *xTensorPtr = xTensor.data_ptr<float>();

    for (int i = 0; i < numNodes; ++i)
    {
        posTensorPtr[i * 3 + 0] = pos[i].GetX();
        posTensorPtr[i * 3 + 1] = pos[i].GetY();
        posTensorPtr[i * 3 + 2] = pos[i].GetZ();
        xTensorPtr[i] = nodeFeatures[i][0];
    }

    inputs.insert(inputs.end(), {xTensor, posTensor, batchTensor});

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::ExtractCandidateVertices(
    const std::vector<CartesianVector> &pos, const torch::Tensor &vertexLogits, std::vector<int> &candidateHitIdx) const
{
    candidateHitIdx.clear();

    const auto logits = vertexLogits.reshape({-1}).to(torch::kFloat32).contiguous();
    if (static_cast<size_t>(logits.size(0)) != pos.size())
        return STATUS_CODE_INVALID_PARAMETER;

    const float *logitData = logits.data_ptr<float>();
    const float logitCut = std::log(m_vertexScoreThreshold / (1.f - m_vertexScoreThreshold));

    std::vector<int> valid;
    for (size_t i = 0; i < pos.size(); ++i)
    {
        if (logitData[i] > logitCut)
            valid.push_back(static_cast<int>(i));
    }

    std::stable_sort(valid.begin(), valid.end(), [&](const int a, const int b) { return logitData[a] > logitData[b]; });

    const float radiusSq = m_nmsRadius * m_nmsRadius;
    std::vector<bool> suppressed(valid.size(), false);

    for (size_t i = 0; i < valid.size(); ++i)
    {
        if (suppressed[i])
            continue;

        candidateHitIdx.push_back(valid[i]);

        for (size_t j = i + 1; j < valid.size(); ++j)
        {
            if (!suppressed[j] && (pos[valid[j]] - pos[valid[i]]).GetMagnitudeSquared() < radiusSq)
                suppressed[j] = true;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_caloHitListName));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
        return STATUS_CODE_SUCCESS;
    }

    std::string modelName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
    modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
    LArDLHelper::LoadModel(modelName, m_modelFile);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexScoreThreshold", m_vertexScoreThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NmsRadius", m_nmsRadius));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
