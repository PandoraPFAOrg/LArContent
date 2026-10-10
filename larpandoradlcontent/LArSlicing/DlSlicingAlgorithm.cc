/**
 *  @file   larpandoradlcontent/LArSlicing/DlSlicingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning slicing algorithm.
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

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArSlicing/DlSlicingAlgorithm.h"
#include "larpandoradlcontent/LArSlicing/HoughFinder.h"
#include "larpandoradlcontent/LArSlicing/KnnKDTree.h"

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

DlSlicingAlgorithm::DlSlicingAlgorithm() :
    m_scalingFactor{-1.0f},
    m_thresholds{},
    m_nDistanceClasses{-1},
    m_runPostProcessing{false}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Run()
{
    std::cout << "Starting DL Slicing Algorithm..." << std::endl;
    return this->Infer();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Infer()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const unsigned int numHits = pCaloHitList->size();

    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<CartesianVector> nodes;
    std::vector<std::array<float, 1>> node_features;
    this->GetNodeData(*pCaloHitList, nodes, node_features);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Getting node data took " << duration << " ms." << std::endl;

    LArDLHelper::TorchInputVector inputs;
    t1 = std::chrono::high_resolution_clock::now();
    this->BuildInput(inputs, nodes, node_features);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Building input took " << duration << " ms." << std::endl;

    // Free the C++ intermediate containers: they are fully encoded in tensors now.
    node_features.clear();
    node_features.shrink_to_fit();

    // Tensors for inference results...
    torch::Tensor instancePreds;
    torch::Tensor fullInstancePreds;

    // Data structures for post-processing...
    std::vector<CartesianVector> foundVertices;
    std::vector<int> candidateIndices;
    std::vector<int> noiseMask(numHits, 0);

    {
        LArDLHelper::TorchMultiOutput semanticOutput;
        t1 = std::chrono::high_resolution_clock::now();
        LArDLHelper::Forward(m_modelFile, inputs, semanticOutput);
        t2 = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Inference took " << duration << " ms." << std::endl;

        // Now, we can process the outputs.
        // We should have 3 things:
        //
        // 1) Semantic distance labels for each node (i.e. hit) which is a class up
        //  to m_numSemanticClasses. This represents the distance of each hit from
        //  its primary neutrino vertex.
        //
        // 2) The raw embeddings for each node. We don't need to do anything with
        // these, except use them later on.
        //
        // 3) The encoded positions, same as above. Just saves re-computing them.
        const auto semanticLabels = semanticOutput.toTuple()->elements()[0].toTensor();
        const auto hitEmbeddings = semanticOutput.toTuple()->elements()[1].toTensor();
        const auto encodedPos = semanticOutput.toTuple()->elements()[2].toTensor();
        semanticOutput = LArDLHelper::TorchMultiOutput();

        // Lets do some basic checks...
        std::cout << "Semantic Labels: " << semanticLabels.sizes() << ", " << semanticLabels.dtype() << std::endl;
        std::cout << "Raw Embeddings: " << hitEmbeddings.sizes() << ", " << hitEmbeddings.dtype() << std::endl;
        std::cout << "Pos Embeddings: " << encodedPos.sizes() << ", " << encodedPos.dtype() << std::endl;

#if DEBUG_MODE
        // DEBUG: Add visualization of the semantic labels to EVD, to check they
        // look sensible before we try to do any more complicated processing.
        const auto argMaxLabels = torch::argmax(semanticLabels, 1);
        HepEVD::Hits hitsToVis;

        int evdHitIdx{0};
        for (const auto pCaloHit : *pCaloHitList)
        {
            if (nullptr == pCaloHit)
                continue;

            const double label = argMaxLabels[evdHitIdx].item<double>();

            const auto x = pCaloHit->GetPositionVector().GetX();
            const auto y = pCaloHit->GetPositionVector().GetY();
            const auto z = pCaloHit->GetPositionVector().GetZ();
            const auto e = pCaloHit->GetInputEnergy();

            HepEVD::Hit *evdHit = new HepEVD::Hit({x, y, z}, e);
            evdHit->addProperties({{"SemanticLabel", label}});

            if (label <= 2)
                evdHit->addProperties({{"SeedCandidate", 1}});

            hitsToVis.push_back(evdHit);

            evdHitIdx++;
        }

        HepEVD::setHepEVDGeometry(this->GetPandora().GetGeometry());
        HepEVD::getServer()->addHits(hitsToVis);
#endif

        // Next, process the semantic labels with the Hough Transform to find vertex
        // candidates. Scope the working buffers so they're freed before instance seg.
        {
            const auto contiguousSemanticLabels = semanticLabels.contiguous();
            std::vector<float> semanticLabelsVec(
                contiguousSemanticLabels.data_ptr<float>(), contiguousSemanticLabels.data_ptr<float>() + (numHits * m_nDistanceClasses));

            // Setup and run the Hough Transform vertex finder.
            t1 = std::chrono::high_resolution_clock::now();
            ThreeDVertexingHoughTransform houghFinder(m_thresholds, m_scalingFactor);
            foundVertices = houghFinder.Fit(nodes, semanticLabelsVec);
            t2 = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            std::cout << "Hough Transform vertex finding took " << duration << " ms." << std::endl;
            std::cout << "Found " << foundVertices.size() << " vertex candidates." << std::endl;
        }

        // Store some info about the candidates for later post processing...
        for (int i = 0; i < numHits; ++i)
        {
            for (const auto &vtx : foundVertices)
            {
                if ((nodes[i] - vtx).GetMagnitudeSquared() < 100.f)
                {
                    candidateIndices.push_back(i);
                    break;
                }
            }
        }

        // Similarly, pull out a noise mask before we clear that.
        const int noiseClass = m_nDistanceClasses - 1;
        auto semanticPreds = torch::argmax(semanticLabels, 1);
        for (int i = 0; i < numHits; ++i)
        {
            if (semanticPreds[i].item<int>() == noiseClass)
            {
                noiseMask[i] = 1;
            }
        }

        // nodes is no longer needed after the Hough finder.
        nodes.clear();
        nodes.shrink_to_fit();

#if DEBUG_MODE
        // DEBUG: Add them to HepEVD.
        HepEVD::Markers pointsToVis;
        for (const auto &vertex : foundVertices)
        {
            HepEVD::Point *evdPoint = new HepEVD::Point({vertex.GetX(), vertex.GetY(), vertex.GetZ()});
            pointsToVis.push_back(*evdPoint);
        }

        HepEVD::getServer()->addMarkers(pointsToVis);
        HepEVD::saveState("FoundVertices");
#endif

        // Start setting up the inputs for instance segmentation.
        const int numCandidates = foundVertices.size();

        if (numCandidates == 0)
        {
            std::cout << "DLSlicingAlgorithm::Infer - no vertex candidates found, skipping instance segmentation step" << std::endl;
            // TODO: What to do here?
            return STATUS_CODE_SUCCESS;
        }

        const auto asFloat = torch::TensorOptions().dtype(torch::kFloat32);
        LArDLHelper::TorchInput candidateTensor;
        LArDLHelper::InitialiseInput({numCandidates, 3}, candidateTensor, asFloat);
        float *candidateTensorData = candidateTensor.data_ptr<float>();

        // Populate the candidate tensor with the positions of the found vertices.
        for (unsigned int i = 0; i < numCandidates; ++i)
        {
            candidateTensorData[i * 3 + 0] = foundVertices[i].GetX();
            candidateTensorData[i * 3 + 1] = foundVertices[i].GetY();
            candidateTensorData[i * 3 + 2] = foundVertices[i].GetZ();
        }

        // Now, populate the full input tensor with all the required data:
        // 1) The position embeddings for each hit.
        // 2) The raw embeddings for each hit.
        // 3) The semantic distance logits for each hit.
        // 4) The candidate vertex positions.
        inputs.clear();
        LArDLHelper::TorchInputVector fullInputTensor;
        fullInputTensor.push_back(std::move(encodedPos));
        fullInputTensor.push_back(std::move(hitEmbeddings));
        fullInputTensor.push_back(std::move(semanticLabels));
        fullInputTensor.push_back(std::move(candidateTensor));

        // Okay, we are good to go!
        t1 = std::chrono::high_resolution_clock::now();
        LArDLHelper::TorchOutput instanceOutput;
        LArDLHelper::Forward(m_modelFile, fullInputTensor, instanceOutput, "predict_instances");
        fullInputTensor.clear();
        t2 = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Instance segmentation inference took " << duration << " ms." << std::endl;

        std::cout << "DLSlicingAlgorithm::Infer - instance segmentation output: " << instanceOutput.sizes() << ", "
                  << instanceOutput.dtype() << std::endl;
        instancePreds = std::get<1>(torch::max(torch::sigmoid(instanceOutput), 1));
        fullInstancePreds = torch::sigmoid(instanceOutput);
        instanceOutput = at::Tensor();
    }

#if DEBUG_MODE
    // DEBUG: Visualise the pre-post-processing clusters.
    std::map<int, HepEVD::Hits> instanceHitsMap;
    int evdHitIdx{0};

    for (const auto pCaloHit : *pCaloHitList)
    {
        if (nullptr == pCaloHit)
            continue;

        const auto x = pCaloHit->GetPositionVector().GetX();
        const auto y = pCaloHit->GetPositionVector().GetY();
        const auto z = pCaloHit->GetPositionVector().GetZ();
        const auto e = pCaloHit->GetInputEnergy();

        const auto clusterPrediction = instancePreds[evdHitIdx].item<int>();

        HepEVD::Hit *evdHit = new HepEVD::Hit({x, y, z}, e);

        for (int i = 0; i < fullInstancePreds.size(1); ++i)
        {
            const auto confidence = fullInstancePreds[evdHitIdx][i].item<float>();
            evdHit->addProperties({{"InstanceConfidence_" + std::to_string(i), confidence}});
        }

        instanceHitsMap[clusterPrediction].push_back(evdHit);

        ++evdHitIdx;
    }

    // DEBUG: Flatten to particles.
    HepEVD::Particles particlesToVis;
    unsigned int clusterId = 0;
    for (const auto &[clusterId, hits] : instanceHitsMap)
    {
        std::cout << "Cluster " << clusterId << ": " << hits.size() << " hits" << std::endl;
        HepEVD::Particle *evdParticle = new HepEVD::Particle(hits, std::to_string(clusterId));
        particlesToVis.push_back(evdParticle);
    }

    HepEVD::getServer()->addParticles(particlesToVis);
    HepEVD::saveState("Slicing Result");
#endif

    // Final cluster labels, to be updated by the post-processing if ran.
    std::vector<int> finalLabels(instancePreds.size(0), -1);

    if (m_runPostProcessing)
    {

        std::cout << "Starting Post-Processing..." << std::endl;
        t1 = std::chrono::high_resolution_clock::now();

        // Rebuild the positions array from the native Pandora hits, as it is needed for the post-processing steps.
        std::vector<pandora::CartesianVector> positions;
        positions.reserve(pCaloHitList->size());
        for (const auto pCaloHit : *pCaloHitList)
        {
            if (pCaloHit)
                positions.push_back(pCaloHit->GetPositionVector());
        }

        std::vector<int> originalLabels(numHits);
        std::vector<int> cleanLabels(numHits);

        // TODO: Populate the t0 information once a sensible input format is decided on.
        std::vector<float> t0s;
        std::vector<bool> t0Valid;

        // Unpack DL output and apply the noise mask we saved earlier
        for (int i = 0; i < numHits; ++i)
        {
            int label = instancePreds[i].item<int>();
            originalLabels[i] = label;
            cleanLabels[i] = (noiseMask[i] == 1) ? -1 : label;
        }

        // Perform the stages of actual post processing...
        // First, split up the large clusters, to highly pure anchor and debris clusters...
        t2 = std::chrono::high_resolution_clock::now();
        auto postProcessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Pre-processing for post-processing took " << postProcessingDuration << " ms." << std::endl;
        t1 = std::chrono::high_resolution_clock::now();

        std::vector<int> splitLabels;
        std::set<int> anchors, debris;
        this->SplitAndClassifyClusters(positions, cleanLabels, candidateIndices, splitLabels, anchors, debris, 10.0f, 50);
        t2 = std::chrono::high_resolution_clock::now();
        postProcessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Split and classify clusters took " << postProcessingDuration << " ms." << std::endl;

        // Then, start to attach the split clusters back together...
        std::vector<int> floodLabels = splitLabels;
        t1 = std::chrono::high_resolution_clock::now();
        this->FloodFill(positions, t0s, t0Valid, cleanLabels, floodLabels, anchors, debris, 20.f, 25.0f);
        t2 = std::chrono::high_resolution_clock::now();
        postProcessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Flood fill took " << postProcessingDuration << " ms." << std::endl;

        // And finally, perform clean up, ensuring every hit ends up in a cluster.
        finalLabels = floodLabels;
        t1 = std::chrono::high_resolution_clock::now();
        this->CleanSmallClusters(positions, t0s, t0Valid, originalLabels, finalLabels, 450);
        t2 = std::chrono::high_resolution_clock::now();
        postProcessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Clean small clusters took " << postProcessingDuration << " ms." << std::endl;
    }
    else
    {
        for (int i = 0; i < numHits; ++i)
        {
            int label = instancePreds[i].item<int>();
            finalLabels[i] = (noiseMask[i] == 1) ? -1 : label;
        }
    }

    // Build the final cluster -> hit map, so we can save it.
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

#if DEBUG_MODE
    // DEBUG: Visualise the post-post-processing clusters.
    instanceHitsMap.clear();
    evdHitIdx = 0;

    for (const auto pCaloHit : *pCaloHitList)
    {
        if (nullptr == pCaloHit)
            continue;

        const auto x = pCaloHit->GetPositionVector().GetX();
        const auto y = pCaloHit->GetPositionVector().GetY();
        const auto z = pCaloHit->GetPositionVector().GetZ();
        const auto e = pCaloHit->GetInputEnergy();

        const auto clusterPrediction = finalLabels[evdHitIdx];

        HepEVD::Hit *evdHit = new HepEVD::Hit({x, y, z}, e);
        instanceHitsMap[clusterPrediction].push_back(evdHit);

        ++evdHitIdx;
    }

    // DEBUG: Flatten to particles.
    particlesToVis.clear();
    clusterId = 0;

    for (const auto &[clusterId, hits] : instanceHitsMap)
    {
        std::cout << "Cluster " << clusterId << ": " << hits.size() << " hits" << std::endl;
        HepEVD::Particle *evdParticle = new HepEVD::Particle(hits, std::to_string(clusterId));
        particlesToVis.push_back(evdParticle);
    }

    HepEVD::getServer()->addParticles(particlesToVis);
    HepEVD::saveState("Post Processing Result");
    HepEVD::startServer();
#endif

    // For now...lets just write out a new Cluster list, with one cluster per
    // predicted instance.
    // This will likely eventually need to be made into a LArRecoND algorithm,
    // that is based on EventSlicingThreeD, but this will work for now
    // and we can have a basic LArRecoND tool that just loads this cluster list.
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

    // We should also write out the predicted vertices to a new output list, as
    // they may be useful for seeding later algorithms, even if the instance
    // segmentation fails.
    const VertexList *pVertexList(nullptr);
    temporaryListName.clear();
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));
    std::cout << "Writing out " << foundVertices.size() << " vertex candidates to " << m_outputVertexListName << std::endl;

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

StatusCode DlSlicingAlgorithm::GetNodeData(const CaloHitList &caloHits, std::vector<CartesianVector> &pos, std::vector<std::array<float, 1>> &node_features)
{
    pos.reserve(caloHits.size());
    node_features.reserve(caloHits.size());

    // Populate the positional and node feature vectors from the CaloHits.
    int hitIdx{0};
    for (const auto pCaloHit : caloHits)
    {
        if (nullptr == pCaloHit)
            continue;

        const CartesianVector &hitPos = pCaloHit->GetPositionVector();

        pos.push_back(hitPos);
        node_features.push_back({pCaloHit->GetInputEnergy() / 10.f});
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::BuildInput(
    LArDLHelper::TorchInputVector &inputs, std::vector<CartesianVector> &pos, std::vector<std::array<float, 1>> &node_features)
{
    const int numNodes{static_cast<int>(pos.size())};
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
        posTensorPtr[i * 3 + 0] = pos[i].GetX();
        posTensorPtr[i * 3 + 1] = pos[i].GetY();
        posTensorPtr[i * 3 + 2] = pos[i].GetZ();
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

StatusCode DlSlicingAlgorithm::SplitAndClassifyClusters(const std::vector<pandora::CartesianVector> &positions,
    const std::vector<int> &clusterLabels, const std::vector<int> &candidateIndices, std::vector<int> &newLabels, std::set<int> &anchors,
    std::set<int> &debris, const float distanceThreshold, const int minAnchorSize) const
{
    const int nHits = positions.size();
    newLabels.assign(nHits, -1);
    std::vector<bool> visited(nHits, false);
    int currentLabel = 0;
    const float thresholdSq = distanceThreshold * distanceThreshold;

    // Identify unique clusters, ignoring noise hits (label < 0)
    std::set<int> uniqueClusters;
    for (int label : clusterLabels)
    {
        if (label >= 0)
            uniqueClusters.insert(label);
    }

    // Split the clusters using BFS.
    //
    // Essentially, for each unvisted hit in a cluster, do a BFS to find all
    // hits in the same cluster within the distance threshold, and label them
    // with the same new cluster ID.
    for (int origCluster : uniqueClusters)
    {
        for (int i = 0; i < nHits; ++i)
        {
            if (clusterLabels[i] != origCluster || visited[i])
                continue;

            std::queue<int> queue;
            queue.push(i);
            visited[i] = true;
            newLabels[i] = currentLabel;

            while (!queue.empty())
            {
                int currIdx = queue.front();
                queue.pop();
                const pandora::CartesianVector &currPos = positions[currIdx];

                for (int j = 0; j < nHits; ++j)
                {
                    if (clusterLabels[j] != origCluster || visited[j])
                        continue;

                    if ((positions[j] - currPos).GetMagnitudeSquared() < thresholdSq)
                    {
                        visited[j] = true;
                        newLabels[j] = currentLabel;
                        queue.push(j);
                    }
                }
            }
            currentLabel++;
        }
    }

    // Add the noise hits back in with their original labels.
    for (int i = 0; i < nHits; ++i)
    {
        if (clusterLabels[i] < 0)
            newLabels[i] = clusterLabels[i];
    }

    // Classify the clusters as either anchors or debris, based on their proximity to candidate vertices and their size.
    // There should only be a single anchor per original cluster.
    // Later algorithms can then use the anchors as a seed for growing the clusters back out.
    for (int origId : uniqueClusters)
    {
        std::map<int, std::vector<int>> subClusters;
        for (int i = 0; i < nHits; ++i)
        {
            if (clusterLabels[i] == origId && newLabels[i] >= 0)
            {
                subClusters[newLabels[i]].push_back(i);
            }
        }

        if (subClusters.empty())
            continue;

        int anchorId = -1;
        int maxVertexHits = 0;
        int maxTotalHits = 0;
        int fallbackAnchorId = -1;

        for (const auto &[scId, scHits] : subClusters)
        {
            int vertexHitCount = 0;
            for (int hitIdx : scHits)
            {
                if (std::find(candidateIndices.begin(), candidateIndices.end(), hitIdx) != candidateIndices.end())
                    vertexHitCount++;
            }

            if (vertexHitCount > maxVertexHits && scHits.size() >= minAnchorSize)
            {
                maxVertexHits = vertexHitCount;
                anchorId = scId;
            }

            if (scHits.size() > maxTotalHits && scHits.size() >= minAnchorSize)
            {
                maxTotalHits = scHits.size();
                fallbackAnchorId = scId;
            }
        }

        if (anchorId == -1)
            anchorId = fallbackAnchorId;

        if (anchorId != -1)
            anchors.insert(anchorId);

        for (const auto &[scId, scHits] : subClusters)
        {
            if (scId != anchorId)
                debris.insert(scId);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float GetMedianT0(std::vector<float> &t0s)
{
    if (t0s.empty())
        return -999.f;
    size_t n = t0s.size() / 2;
    std::nth_element(t0s.begin(), t0s.begin() + n, t0s.end());

    if (t0s.size() % 2 == 0)
    {
        auto max_it = std::max_element(t0s.begin(), t0s.begin() + n);
        return (*max_it + t0s[n]) / 2.0f;
    }

    return t0s[n];
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::unordered_map<int, float> DlSlicingAlgorithm::PrecomputeT0Medians(
    const std::vector<int> &labels, const std::vector<float> &t0s, const std::vector<bool> &t0Valid) const
{
    std::unordered_map<int, float> medians;
    if (t0s.empty() || t0Valid.empty())
        return medians;

    std::unordered_map<int, std::vector<float>> t0Map;
    for (size_t i = 0; i < labels.size(); ++i)
    {
        if (labels[i] >= 0 && t0Valid[i])
        {
            t0Map[labels[i]].push_back(t0s[i]);
        }
    }
    for (auto &[cid, t0List] : t0Map)
    {
        medians[cid] = GetMedianT0(t0List);
    }
    return medians;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::FloodFill(const std::vector<pandora::CartesianVector> &positions, const std::vector<float> &t0s,
    const std::vector<bool> &t0Valid, const std::vector<int> &originalLabels, std::vector<int> &finalLabels, const std::set<int> &anchors,
    const std::set<int> &debris, const float baseGap, const float ogBonusGap) const
{
    const bool useT0 = !t0s.empty();
    auto clusterMedians = this->PrecomputeT0Medians(finalLabels, t0s, t0Valid);

    // Pre-group hits by their current cluster ID to speed up later usage.
    std::unordered_map<int, std::vector<int>> clusterToHits;
    for (size_t i = 0; i < finalLabels.size(); ++i)
    {
        if (finalLabels[i] >= 0)
            clusterToHits[finalLabels[i]].push_back(i);
    }

    std::set<int> activeDebris = debris;

    for (int anchorId : anchors)
    {
        if (!clusterToHits.count(anchorId))
            continue;
        std::set<int> activeClusters = {anchorId};

        // Find OG ID for gap bonus.
        std::unordered_map<int, int> origCounts;
        int bestOrigId = -1, maxCount = 0;
        for (int idx : clusterToHits[anchorId])
        {
            int oId = originalLabels[idx];
            if (oId >= 0 && ++origCounts[oId] > maxCount)
            {
                maxCount = origCounts[oId];
                bestOrigId = oId;
            }
        }

        std::unordered_set<int> ogDebris;
        for (size_t i = 0; i < finalLabels.size(); ++i)
        {
            if (originalLabels[i] == bestOrigId && debris.count(finalLabels[i]))
            {
                ogDebris.insert(finalLabels[i]);
            }
        }

        bool addedThisRound = true;
        while (addedThisRound)
        {
            addedThisRound = false;

            // Build Tree of active clusters
            std::vector<KnnKdTree::KnnNode> activeNodes;
            for (int cId : activeClusters)
            {
                for (int idx : clusterToHits[cId])
                {
                    KnnKdTree::KnnNode node;
                    node.coords[0] = positions[idx].GetX();
                    node.coords[1] = positions[idx].GetY();
                    node.coords[2] = positions[idx].GetZ();
                    node.original_id = idx;
                    activeNodes.push_back(node);
                }
            }

            if (activeNodes.empty())
                break;
            KnnKdTree tree(activeNodes);

            // Calculate current anchor median
            float currentAnchorMedian = -999.f;
            if (useT0)
            {
                std::vector<float> anchorT0s;
                for (int cId : activeClusters)
                {
                    for (int idx : clusterToHits[cId])
                    {
                        if (t0Valid[idx])
                            anchorT0s.push_back(t0s[idx]);
                    }
                }
                currentAnchorMedian = GetMedianT0(anchorT0s);
            }

            std::vector<int> toRemove;
            for (int candId : activeDebris)
            {
                if (!clusterToHits.count(candId))
                    continue;

                float minDist = 99999.f;
                for (int idx : clusterToHits[candId])
                {
                    KnnKdTree::KnnNode qNode;
                    qNode.coords[0] = positions[idx].GetX();
                    qNode.coords[1] = positions[idx].GetY();
                    qNode.coords[2] = positions[idx].GetZ();

                    auto nn = tree.FindNearestNeighbours(qNode, 1);
                    if (!nn.empty())
                    {
                        float dist = (positions[nn[0]] - positions[idx]).GetMagnitude();
                        if (dist < minDist)
                            minDist = dist;
                    }
                }

                // Is the closest approach to the active clusters within the allowed gap? If not, skip it.
                float allowedGap = ogDebris.count(candId) ? ogBonusGap : baseGap;
                if (minDist >= allowedGap)
                    continue;

                // And then check the t0, if in use.
                if (useT0 && currentAnchorMedian > -900.f && clusterMedians.count(candId))
                {
                    float candMedian = clusterMedians[candId];
                    if (candMedian > -900.f && std::abs(currentAnchorMedian - candMedian) > 1.0f)
                        continue;
                }

                // Survived all checks! Merge it.
                activeClusters.insert(candId);
                toRemove.push_back(candId);
                addedThisRound = true;
            }

            for (int id : toRemove)
                activeDebris.erase(id);
        }

        // Apply final merges for this anchor
        for (int sc : activeClusters)
        {
            if (sc == anchorId)
                continue;

            for (int idx : clusterToHits[sc])
                finalLabels[idx] = anchorId;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::CleanSmallClusters(const std::vector<pandora::CartesianVector> &positions, const std::vector<float> &t0s,
    const std::vector<bool> &t0Valid, const std::vector<int> &originalLabels, std::vector<int> &finalLabels, const int minSize) const
{
    // Get t0 medians, if using.
    auto clusterMedians = this->PrecomputeT0Medians(finalLabels, t0s, t0Valid);

    // Pre-group hits by their current and original clusters to speed up later usage.
    std::unordered_map<int, std::vector<int>> currentClusterHits;
    std::unordered_map<int, std::unordered_map<int, int>> origToCurrCounts;

    for (size_t i = 0; i < finalLabels.size(); ++i)
    {
        if (finalLabels[i] >= 0)
            currentClusterHits[finalLabels[i]].push_back(i);

        if (originalLabels[i] >= 0 && finalLabels[i] >= 0)
            origToCurrCounts[originalLabels[i]][finalLabels[i]]++;
    }

    // Find dominant current ID for each original ID.
    std::unordered_map<int, int> origToDominant;
    for (const auto &[origId, currCounts] : origToCurrCounts)
    {
        int bestCurr = -1, maxC = 0;
        for (const auto &[cId, count] : currCounts)
        {
            if (count > maxC)
            {
                maxC = count;
                bestCurr = cId;
            }
        }
        if (bestCurr >= 0)
            origToDominant[origId] = bestCurr;
    }

    // Process small clusters
    for (const auto &[currId, maskIndices] : currentClusterHits)
    {
        if (maskIndices.size() >= minSize)
            continue;

        // Find majority original ID
        std::unordered_map<int, int> origCounts;
        int bestOrig = -1, maxC = 0;
        for (int idx : maskIndices)
        {
            int oId = originalLabels[idx];
            if (oId >= 0 && ++origCounts[oId] > maxC)
            {
                maxC = origCounts[oId];
                bestOrig = oId;
            }
        }

        if (bestOrig < 0 || !origToDominant.count(bestOrig))
            continue;

        int targetLabel = origToDominant[bestOrig];

        // Apply t0 veto, if in use, to avoid merging things with very different predicted t0s.
        if (!t0s.empty())
        {
            float currMedian = clusterMedians.count(currId) ? clusterMedians[currId] : -999.f;
            float targetMedian = clusterMedians.count(targetLabel) ? clusterMedians[targetLabel] : -999.f;

            if (currMedian > -900.f && targetMedian > -900.f && std::abs(currMedian - targetMedian) > 2.0f)
                continue;
        }

        // Survived all checks, apply merge!
        for (int idx : maskIndices)
            finalLabels[idx] = targetLabel;
    }

    // Final noise sweep...every hit needs a home!
    std::vector<KnnKdTree::KnnNode> signalNodes;
    std::vector<int> noiseIndices;

    for (size_t i = 0; i < finalLabels.size(); ++i)
    {
        if (finalLabels[i] >= 0)
        {
            KnnKdTree::KnnNode node;
            node.coords[0] = positions[i].GetX();
            node.coords[1] = positions[i].GetY();
            node.coords[2] = positions[i].GetZ();
            node.original_id = i;
            signalNodes.push_back(node);
        }
        else
            noiseIndices.push_back(i);
    }

    if (!signalNodes.empty() && !noiseIndices.empty())
    {
        KnnKdTree tree(signalNodes);

        for (int noiseIdx : noiseIndices)
        {
            KnnKdTree::KnnNode qNode;
            qNode.coords[0] = positions[noiseIdx].GetX();
            qNode.coords[1] = positions[noiseIdx].GetY();
            qNode.coords[2] = positions[noiseIdx].GetZ();

            auto nn = tree.FindNearestNeighbours(qNode, 1);
            if (!nn.empty())
                finalLabels[noiseIdx] = finalLabels[nn[0]];
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    std::string modelName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
    modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
    LArDLHelper::LoadModel(modelName, m_modelFile);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));

    m_nDistanceClasses = m_thresholds.size() + 1; // We have one more class than thresholds, as the thresholds define the boundaries between classes.

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
