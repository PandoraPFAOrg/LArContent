/**
 *  @file   src/DLVertexingThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning 3D vertexing algorithm.
 *
 *  $Log: $
 */

#include <chrono>
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

#include "larpandoradlcontent/LArSlicing/DlVertexingThreeDAlgorithm.h"
#include "larpandoradlcontent/LArSlicing/HoughFinder.h"

#include <Eigen/Dense>
#include <c10/core/TensorOptions.h>
#include <torch/script.h>
#include <torch/torch.h>

#define DEBUG_MODE 0
#if DEBUG_MODE
#define HEP_EVD_PANDORA_HELPERS 1
#include "hep_evd.h"
#endif

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlThreeDVertexingAlgorithm::DlThreeDVertexingAlgorithm() :
    m_scalingFactor{-1.0f},
    m_thresholds{},
    m_nDistanceClasses{-1},
    m_pass{1}
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

    // Load the Slicing Worker instance, to pull out the slice PFOs + their associated candidate vertices.
    // We can use them as a start point for where to zoom in on.
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

    // We've found the slicing worker, and can now try and pull out the vertex for this slice.
    if (pSlicingInstance != nullptr)
    {
        const PfoList *pSlicePfos(nullptr);
        const auto slicingPfos = PandoraApi::GetCurrentPfoList(*pSlicingInstance, pSlicePfos);
        std::cout << "Slicing instance has " << pSlicePfos->size() << " PFOs." << std::endl;

        CartesianPointVector candidateVertices({});
        for (const auto pPfo : *pSlicePfos)
        {
            const auto vertices = pPfo->GetVertexList();
            for (const auto pVertex : vertices)
                candidateVertices.push_back(pVertex->GetPosition());
        }

        std::cout << "Found " << candidateVertices.size() << " candidate vertices from slicing." << std::endl;

        // Process the candidate vertices, and find the most appropriate for this slice.
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

            // Finally, write out the predicted vertices to the new output list.
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

    // INFO: If we are here...then we failed to load or find a vertex for this slice.
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
    std::cout << "DLVertexingThreeDAlgorithm::RunModel - Pass " << m_pass << " - Running model with " << caloHits.size() << " calo hits." << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<CartesianVector> nodes;
    std::vector<std::array<float, 1>> node_features;
    this->GetNodeData(caloHits, nodes, node_features, cropVertex);
    const int numHits(nodes.size());
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Getting node data took " << duration << " ms." << std::endl;

    if (numHits == 0)
    {
        std::cout << "DLVertexingThreeDAlgorithm::RunModel - No calo hits to build graph with, skipping." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    LArDLHelper::TorchInputVector inputs;
    t1 = std::chrono::high_resolution_clock::now();
    this->BuildInput(inputs, nodes, node_features);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Building input took " << duration << " ms." << std::endl;

    // Free the C++ intermediate containers: they are fully encoded in tensors now.
    node_features.clear();
    node_features.shrink_to_fit();

    // Store the final vertex candidates for saving to a Pandora list.
    std::vector<CartesianVector> foundVertices;

    {
        LArDLHelper::TorchMultiOutput semanticOutput;
        t1 = std::chrono::high_resolution_clock::now();
        LArDLHelper::Forward(m_modelFile, inputs, semanticOutput);
        t2 = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Inference took " << duration << " ms." << std::endl;

        // Pull out the semantic labels before clearing the container to free memory.
        const auto semanticLabels = semanticOutput.toTuple()->elements()[0].toTensor();
        semanticOutput = LArDLHelper::TorchMultiOutput();

        // Lets do some basic checks...
        std::cout << "Semantic Labels: " << semanticLabels.sizes() << ", " << semanticLabels.dtype() << std::endl;

#if DEBUG_MODE
        // DEBUG: Add visualization of the semantic labels to EVD.
        HepEVD::Hits hitsToVis;

        int evdHitIdx{0};
        for (const auto pCaloHit : caloHits)
        {
            if (nullptr == pCaloHit)
                continue;

            const auto x = pCaloHit->GetPositionVector().GetX();
            const auto y = pCaloHit->GetPositionVector().GetY();
            const auto z = pCaloHit->GetPositionVector().GetZ();
            const auto e = pCaloHit->GetInputEnergy();

            HepEVD::Hit *evdHit = new HepEVD::Hit({x, y, z}, e);

            hitsToVis.push_back(evdHit);

            evdHitIdx++;
        }

        HepEVD::setHepEVDGeometry(this->GetPandora().GetGeometry());
        HepEVD::getServer()->addHits(hitsToVis);
        HepEVD::saveState("RawHits");

        const auto argMaxLabels = torch::argmax(semanticLabels, 1);
        const torch::Tensor confidences = torch::softmax(semanticLabels, /*dim=*/1);
        evdHitIdx = 0;
        hitsToVis.clear();

        double lowestClassFound = std::numeric_limits<double>::max();
        double maxConfForLowestClass = -1.0;
        pandora::CartesianVector bestConfHitPos(0, 0, 0);

        for (const auto graphHit : nodes)
        {
            const double label = argMaxLabels[evdHitIdx].item<double>();
            const double confidence = confidences[evdHitIdx][label].item<double>();
            const double seedConfidence = confidences[evdHitIdx][0].item<double>();

            if (label < lowestClassFound)
            {
                lowestClassFound = label;
                maxConfForLowestClass = confidence;
                bestConfHitPos = pandora::CartesianVector(graphHit.GetX(), graphHit.GetY(), graphHit.GetZ());
            }
            else if (label == lowestClassFound && confidence > maxConfForLowestClass)
            {
                maxConfForLowestClass = confidence;
                bestConfHitPos = pandora::CartesianVector(graphHit.GetX(), graphHit.GetY(), graphHit.GetZ());
            }

            const auto x = graphHit.GetX();
            const auto y = graphHit.GetY();
            const auto z = graphHit.GetZ();
            const auto e = 1.0;

            HepEVD::Hit *evdHit = new HepEVD::Hit({x, y, z}, e);
            evdHit->addProperties({{"SemanticLabel", label}, {"SemanticConfidence", confidence}, {"SeedConfidence", seedConfidence}});

            if (label <= 0)
                evdHit->addProperties({{"SeedCandidate", 1}});

            hitsToVis.push_back(evdHit);

            evdHitIdx++;
        }

        HepEVD::getServer()->addHits(hitsToVis);
#endif

        // Next, process the semantic labels with the Hough Transform to find vertex
        // candidates.
        {
            const auto contiguousSemanticLabels = semanticLabels.contiguous();
            std::vector<float> semanticLabelsVec(
                contiguousSemanticLabels.data_ptr<float>(), contiguousSemanticLabels.data_ptr<float>() + (numHits * m_nDistanceClasses));

            // Setup and run the Hough Transform vertex finder.
            t1 = std::chrono::high_resolution_clock::now();
            ThreeDVertexingHoughTransform houghFinder(m_thresholds, m_scalingFactor, 0.5f, 3, 10.f, 0, true);
            foundVertices = houghFinder.Fit(nodes, semanticLabelsVec);
            t2 = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            std::cout << "Hough Transform vertex finding took " << duration << " ms." << std::endl;
            std::cout << "Found " << foundVertices.size() << " vertex candidates." << std::endl;
        }

        // TODO: We do also have that noise mask here if useful to be stored in the event context etc?

#if DEBUG_MODE
        // DEBUG: Add them to HepEVD.
        HepEVD::Markers pointsToVis;
        for (const auto &vertex : foundVertices)
        {
            HepEVD::Point *evdPoint = new HepEVD::Point({vertex.GetX(), vertex.GetY(), vertex.GetZ()});
            evdPoint->setColour("red");
            pointsToVis.push_back(*evdPoint);
        }

        if (cropVertex != nullptr)
        {
            HepEVD::Point *evdPoint = new HepEVD::Point({cropVertex->GetX(), cropVertex->GetY(), cropVertex->GetZ()});
            evdPoint->setColour("blue");
            pointsToVis.push_back(*evdPoint);
        }

        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
        LArMCParticleHelper::MCRelationMap mcPrimaryMap;
        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

        CaloHitList croppedHits;
        for (const auto &hit : caloHits)
        {
            for (const auto &pos : nodes)
            {
                if (hit->GetPositionVector() == pos)
                {
                    croppedHits.push_back(hit);
                    break;
                }
            }
        }

        LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
        LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
        LArMCParticleHelper::GetMCParticleToCaloHitMatches(&croppedHits, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);
        std::set<const MCParticle *> uniqueMCParticles;

        float bestCandidate = std::numeric_limits<float>::max();
        float bestNew = std::numeric_limits<float>::max();
        float bestHit = std::numeric_limits<float>::max();
        float bestCrop = std::numeric_limits<float>::max();
        float bestConfHitDist = std::numeric_limits<float>::max();

        for (const auto &caloHitMCPair : hitToMCMap)
        {
            const auto &pMCParticle = caloHitMCPair.second;
            const auto primary = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);

            if (uniqueMCParticles.find(primary) != uniqueMCParticles.end())
                continue;

            HepEVD::Point *evdPoint = new HepEVD::Point({primary->GetVertex().GetX(), primary->GetVertex().GetY(), primary->GetVertex().GetZ()});
            evdPoint->setColour("green");
            pointsToVis.push_back(*evdPoint);
            uniqueMCParticles.insert(primary);

            // Calculate how far the network's most confident hit (of the lowest available class) is from the true vertex
            const float squaredDistConf = bestConfHitPos.GetDistanceSquared(primary->GetVertex());
            if (squaredDistConf < bestConfHitDist)
                bestConfHitDist = squaredDistConf;

            if (cropVertex != nullptr)
            {
                const float squaredDist = cropVertex->GetDistanceSquared(primary->GetVertex());
                if (squaredDist < bestCandidate)
                    bestCandidate = squaredDist;
            }

            for (const auto &newVtx : foundVertices)
            {
                const float squaredDist = newVtx.GetDistanceSquared(primary->GetVertex());
                if (squaredDist < bestNew)
                    bestNew = squaredDist;
            }

            for (const auto &graphHit : nodes)
            {
                pandora::CartesianVector hitPos(graphHit.GetX(), graphHit.GetY(), graphHit.GetZ());
                const float squaredDist = hitPos.GetDistanceSquared(primary->GetVertex());
                if (squaredDist < bestCrop)
                    bestCrop = squaredDist;
            }

            for (const auto &caloHit : caloHits)
            {
                const float squaredDist = caloHit->GetPositionVector().GetDistanceSquared(primary->GetVertex());
                if (squaredDist < bestHit)
                    bestHit = squaredDist;
            }
        }

        std::cout << "[METRICS] Slice: " << std::sqrt(bestCandidate) << " | GNN Hough: " << std::sqrt(bestNew)
                  << " | Max Conf Hit: " << std::sqrt(bestConfHitDist) << " | Crop Limit: " << std::sqrt(bestCrop)
                  << " | Voxel Limit: " << std::sqrt(bestHit) << " cm" << std::endl;

        HepEVD::getServer()->addMarkers(pointsToVis);
        HepEVD::saveState("FoundVertices");
        HepEVD::startServer();
#endif

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

    // Finally, write out the predicted vertices to the new output list.
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

    // Perform the crop around the matched vertex.
    // TODO: Move to XML setting.
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

    // Shrink the vectors, in case cropping resulted in significantly fewer hits than the original list size.
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

    // Also create a batch tensor.
    // In python/training land...this is a tensor that tells the model how many graphs are in the batch, and which
    // nodes/edges belong to which graph.
    // In this case, we only have one graph, so we can just set it to 0 for all nodes and edges.
    torch::Tensor batchTensor = torch::zeros(numNodes, torch::kLong);

    // Use raw memory pointers to access the various tensors, to massively speed
    // up writing.
    float *posTensorPtr = posTensor.data_ptr<float>();
    float *xTensorPtr = xTensor.data_ptr<float>();

    // Fill in the position and node feature tensors...
    for (int i = 0; i < numNodes; ++i)
    {
        posTensorPtr[i * 3 + 0] = nodes[i].GetX();
        posTensorPtr[i * 3 + 1] = nodes[i].GetY();
        posTensorPtr[i * 3 + 2] = nodes[i].GetZ();
        xTensorPtr[i] = node_features[i][0];
    }

    // Finally, stick them together into the input vector.
    inputs.insert(inputs.end(), {xTensor, posTensor, batchTensor});

    // Print some debug information
    std::cout << "Nodes: " << posTensor.sizes() << ", " << posTensor.dtype() << std::endl;
    std::cout << "Features: " << xTensor.sizes() << ", " << xTensor.dtype() << std::endl;

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

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
