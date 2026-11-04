/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <unistd.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlHitTrackShowerIdAlgorithm::DlHitTrackShowerIdAlgorithm() :
    m_caloHitListName(""),
    m_trainingMode(false),
    m_maxSeqLen(4096)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlHitTrackShowerIdAlgorithm::~DlHitTrackShowerIdAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlHitTrackShowerIdAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::PrepareTrainingSample()
{
    // Get hits
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Map truth information for hits
    // mcToHitsMap (to OG hits) leadingToHitsMap (fold up showers) instanceToHitsMap (reassign collinear hits)
    // Simultaneously work out how close hits are to vertex
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, false);
    MCFoldingMap mcFoldingMap;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const MCParticle *const pLeading{LArMCParticleHelper::GetLeadingEMParticle(pMC)};
        mcFoldingMap[pMC] = pLeading ? pLeading : pMC;
    }

    LArMCParticleHelper::MCContributionMap leadingToHitsMap;
    this->FoldToLeading(mcToHitsMap, mcFoldingMap, leadingToHitsMap);
    LArMCParticleHelper::MCContributionMap instanceToHitsMap;
    this->ConsolidateInstances(leadingToHitsMap, instanceToHitsMap);

    // Now assign each MCParticle and instance label (ID)
    int inst{1};
    std::unordered_map<const MCParticle *, int> mcToInstanceMap;
    for (const auto &[pMC, dummy] : instanceToHitsMap)
    {
        mcToInstanceMap[pMC] = inst;
        ++inst;
    }

    // Now assign mcparticles/hits to categories
    std::unordered_map<const CaloHit *, const MCParticle *> hitToMCParticle;
    std::unordered_map<const MCParticle *, Category> mcToCategoryMap;
     
    for (const auto &[pMC, caloHitList] : instanceToHitsMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};

        if (!pLArMC)
            continue;

        for (const CaloHit *const pCaloHit : caloHitList)
            hitToMCParticle[pCaloHit] = pMC;
                
        if (this->IsHip(pLArMC))
        {
            mcToCategoryMap[pMC] = HIP;
        }
        else if (this->IsMip(pLArMC))
        {
           mcToCategoryMap[pMC] = MIP;
        }
        else if (caloHitList.size() == 1 || this->IsDiffuse(pLArMC, mcToHitsMap))
        {
            mcToCategoryMap[pMC] = DIFFUSE;
        }
        else if (this->IsShowerLikeEM(pLArMC, mcToHitsMap))
        {
            mcToCategoryMap[pMC] = EM_SHOWER;
        }
        else if (this->IsTrackLikeEM(pLArMC))
        {
            mcToCategoryMap[pMC] = EM_TRACK;
        }
        else
        {
            std::cout << "Warning: Found unclassifiable particle, default to MIP" << std::endl;
            std::cout << "  " << pLArMC->GetParticleId() << " " << pLArMC->GetProcess() << ": " << caloHitList.size() << " hits" << std::endl;
            std::cout << "  Parent PDG code: " << (pLArMC->GetParentList().empty() ? 0 : pLArMC->GetParentList().front()->GetParticleId()) << std::endl;
            mcToCategoryMap[pMC] = MIP;
        }
    }

    if (hitToMCParticle.empty())
        return STATUS_CODE_SUCCESS;
    
    // Prepare input feature vectors
    CaloHitList trainingHitList;    
    for (const auto &[pCaloHit, dummy] : hitToMCParticle)
        trainingHitList.emplace_back(pCaloHit);

    // Get vertex2D
    float vx{0.f}, vz{0.f};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, != , this->GetVertexPlanarCoordinates(trainingHitList.front()->GetHitType(), vx, vz));    

    // Get bounds
    Bounds bounds;
    this->GetCoordinateExtrema(trainingHitList, vx, vz, bounds);

    // Get input vectors
    CaloHitVector orderedTrainingHitVector;
    FloatVector x_rel, z_rel, x_abs, z_abs, rr, cosTheta, sinTheta, wirePitch, wireAngle, adc, width;
    this->PopulateInputVectors(trainingHitList, bounds, vx, vz, x_rel, z_rel, x_abs, z_abs, rr, cosTheta, sinTheta, wirePitch, wireAngle, adc, width, orderedTrainingHitVector);
    
    // Get training
    // Collect hits and respective categories into feature vectors for training
    size_t nHits(trainingHitList.size());
    IntVector semanticLabel(nHits, UNINITIALISED);
    IntVector instanceLabel(nHits, UNINITIALISED);

    for (size_t i = 0; i < nHits; ++i)
    {
        const CaloHit *const pHit = orderedTrainingHitVector[i];
        const MCParticle *const pMC = hitToMCParticle.at(pHit);

        semanticLabel[i] = mcToCategoryMap.at(pMC);
        instanceLabel[i] = mcToInstanceMap.at(pMC);    
    }
    
    // Now save...
    const int event(this->GetPandora().GetEvent()), subrun(this->GetPandora().GetSubrun()), run(this->GetPandora().GetRun());       
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "run", run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "subrun", subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "event", event));

    for (size_t s = 0; s < nHits; s += m_maxSeqLen)
    {
        size_t end{std::min(s + m_maxSeqLen, nHits)};
        FloatVector xRelChunk(x_rel.begin() + s, x_rel.begin() + end);
        FloatVector zRelChunk(z_rel.begin() + s, z_rel.begin() + end);
        FloatVector xAbsChunk(x_abs.begin() + s, x_abs.begin() + end);
        FloatVector zAbsChunk(z_abs.begin() + s, z_abs.begin() + end);
        FloatVector wpChunk(wirePitch.begin() + s, wirePitch.begin() + end);
        FloatVector waChunk(wireAngle.begin() + s, wireAngle.begin() + end);
        FloatVector rrChunk(rr.begin() + s, rr.begin() + end);
        FloatVector cosThetaChunk(cosTheta.begin() + s, cosTheta.begin() + end);
        FloatVector sinThetaChunk(sinTheta.begin() + s, sinTheta.begin() + end);
        FloatVector adcChunk(adc.begin() + s, adc.begin() + end);
        FloatVector widthChunk(width.begin() + s, width.begin() + end);
        IntVector semanticLabelChunk(semanticLabel.begin() + s, semanticLabel.begin() + end);
        IntVector instanceLabelChunk(instanceLabel.begin() + s, instanceLabel.begin() + end);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "wire_pitch", &wpChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "wire_angle", &waChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "x_rel", &xRelChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "z_rel", &zRelChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "x_abs", &xAbsChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "z_abs", &zAbsChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "r", &rrChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "cos_theta", &cosThetaChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "sin_theta", &sinThetaChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "adc", &adcChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "width", &widthChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "semantic_label", &semanticLabelChunk));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "instance_label", &instanceLabelChunk));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Infer()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Get projected vertex coordinates if needed
    float vx{0.f}, vz{0.f};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, != , this->GetVertexPlanarCoordinates(pCaloHitList->front()->GetHitType(), vx, vz));

    // Get coordinate minima and maxima
    Bounds bounds;
    this->GetCoordinateExtrema(*pCaloHitList, vx, vz, bounds);

    // Prepare input feature vectors
    CaloHitVector caloHitVector;
    FloatVector x_rel, z_rel, x_abs, z_abs, rr, cosTheta, sinTheta, wirePitch, wireAngle, adc, width;
    this->PopulateInputVectors(*pCaloHitList, bounds, vx, vz, x_rel, z_rel, x_abs, z_abs, rr, cosTheta, sinTheta, wirePitch, wireAngle, adc, width, caloHitVector);

    // Ensure each run conforms to the maximum sequence length
    size_t totalSize{rr.size()};
    for (size_t s = 0; s < totalSize; s += m_maxSeqLen)
    {
        FloatVector xRelBatch(m_maxSeqLen, -999.f), zRelBatch(m_maxSeqLen, -999.f), xAbsBatch(m_maxSeqLen, -999.f);
        FloatVector zAbsBatch(m_maxSeqLen, -999.f), wpBatch(m_maxSeqLen, -999.f), waBatch(m_maxSeqLen, -999.f);
        FloatVector rrBatch(m_maxSeqLen, -999.f), cosBatch(m_maxSeqLen, -999.f), sinBatch(m_maxSeqLen, -999.f);
        FloatVector adcBatch(m_maxSeqLen, -999.f), widthBatch(m_maxSeqLen, -999.f);
        IntVector maskBatch(m_maxSeqLen, 0);    

        for (int i = 0; i < m_maxSeqLen; ++i)
        {
            size_t iInput(s + i);

            if (iInput < totalSize)
            {
                xRelBatch[i] = x_rel[iInput];
                zRelBatch[i] = z_rel[iInput];
                xAbsBatch[i] = x_abs[iInput];
                zAbsBatch[i] = z_abs[iInput];
                wpBatch[i] = wirePitch[iInput];
                waBatch[i] = wireAngle[iInput];
                rrBatch[i] = rr[iInput];
                cosBatch[i] = cosTheta[iInput];
                sinBatch[i] = sinTheta[iInput];
                adcBatch[i] = adc[iInput];
                widthBatch[i] = width[iInput];
                maskBatch[i] = 1;
            }
        }

        std::vector<LArDLHelper::TorchInput> tensors{
            torch::from_blob(xRelBatch.data(), {m_maxSeqLen}),
            torch::from_blob(zRelBatch.data(), {m_maxSeqLen}),
            torch::from_blob(xAbsBatch.data(), {m_maxSeqLen}),
            torch::from_blob(zAbsBatch.data(), {m_maxSeqLen}),
            torch::from_blob(widthBatch.data(), {m_maxSeqLen}),
            torch::from_blob(adcBatch.data(), {m_maxSeqLen}),
            torch::from_blob(rrBatch.data(), {m_maxSeqLen}),
            torch::from_blob(cosBatch.data(), {m_maxSeqLen}),
            torch::from_blob(sinBatch.data(), {m_maxSeqLen}),            
            torch::from_blob(wpBatch.data(), {m_maxSeqLen}),
            torch::from_blob(waBatch.data(), {m_maxSeqLen})        
        };
        
        // Collate input
        LArDLHelper::TorchInputVector inputs;
        LArDLHelper::TorchInput input{torch::stack(tensors, 1)};
        input = input.unsqueeze(0); // Add batch dimension
        inputs.emplace_back(input);
        
        // Form input mask
        LArDLHelper::TorchInput inputMask = torch::from_blob(maskBatch.data(), {m_maxSeqLen}).to(torch::kBool);
        inputMask = inputMask.unsqueeze(0); // Add batch dimension
        inputs.emplace_back(inputMask);

        // Define output
        LArDLHelper::TorchOutput output;
        
        // Run the input through the trained model
        LArDLHelper::Forward(m_model, inputs, output);
        
        // Get the outputs as tensors - we're assuming batch size of 1, hence squeeze(0) for convenience
        LArDLHelper::TorchOutput classLogits{output.squeeze(0)}; // [N, num_classes]
        LArDLHelper::TorchOutput classProbs{torch::softmax(classLogits, 1)};        
        
        // Save scores
        std::vector<std::string> classLabels({"MIPScore", "HIPScore", "EMTrackScore", "EMShowerScore", "DiffuseScore"});
        for (int i = 0; i < m_maxSeqLen; ++i)
        {
            size_t iInput(s + i);

            if (iInput >= totalSize)
                continue;
            
            const CaloHit *const pCaloHit{caloHitVector[iInput]};
            LArCaloHit *const pLArCaloHit{const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};

            for (unsigned int classIndex = 0; classIndex < classLabels.size(); ++classIndex)
            {
                const float classScore(classProbs[i][classIndex].item<float>());
                const std::string classLabel(classLabels.at(classIndex));
                pLArCaloHit->AddHitScore(classScore, classLabel);
            }            
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::FoldToLeading(const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, const MCFoldingMap &mcFoldingMap,
    lar_content::LArMCParticleHelper::MCContributionMap &leadingHitMap) const
{
    leadingHitMap.clear();
    for (const auto &[pMC, caloHitList] : mcHitMap)
    {
        const MCParticle *const pLeading{mcFoldingMap.at(pMC)};
        leadingHitMap[pLeading].insert(leadingHitMap[pLeading].end(), caloHitList.begin(), caloHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::ConsolidateInstances(const lar_content::LArMCParticleHelper::MCContributionMap &leadingHitMap, 
    lar_content::LArMCParticleHelper::MCContributionMap &instanceHitMap) const
{
    instanceHitMap.clear();
    for (const auto &[pMC, caloHitList] : leadingHitMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC)
        {
            // Ensure delta ray hits along (colinear) muon/pion tracks are assigned to the parent particle
            const int pdg{std::abs(pLArMC->GetParticleId())};
            if (pdg == E_MINUS || pdg == PHOTON)
            {
                CaloHitList particleOwnedHits, parentOwnedHits;
                const LArMCParticle *const pParent{dynamic_cast<const LArMCParticle *>(this->AllocateHitOwner(pLArMC, leadingHitMap,
                    particleOwnedHits, parentOwnedHits))};
                if (!parentOwnedHits.empty())
                    instanceHitMap[pParent].insert(instanceHitMap[pParent].end(), parentOwnedHits.begin(), parentOwnedHits.end());
                if (!particleOwnedHits.empty())
                    instanceHitMap[pLArMC].insert(instanceHitMap[pLArMC].end(), particleOwnedHits.begin(), particleOwnedHits.end());
            }
            else
            {
                if (!caloHitList.empty())
                    instanceHitMap[pLArMC].insert(instanceHitMap[pLArMC].end(), caloHitList.begin(), caloHitList.end());
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsDiffuse(const LArMCParticle *const pLArMC, const LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    if (!pLArMC)
        return false;

    return this->AllHitsCompton(pLArMC, mcToHitsMap, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsHip(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const int pdg{std::abs(pLArMC->GetParticleId())};
    return (pdg > 1e9 || pdg == PROTON || pdg == K_PLUS || pdg == SIGMA_PLUS || pdg == SIGMA_MINUS || pdg == HYPERON_MINUS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsMip(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const int pdg{std::abs(pLArMC->GetParticleId())};
    return (pdg == MU_MINUS || pdg == PI_PLUS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsShowerLikeEM(const LArMCParticle *const pLArMC, const LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    if (!pLArMC)
        return false;

    const int pdg{std::abs(pLArMC->GetParticleId())};

    if (!(pdg == PHOTON || pdg == E_MINUS))
        return false;

    // If it Bremms, it's a shower
    return this->UndergoesBremmstrahlung(pLArMC, mcToHitsMap, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsTrackLikeEM(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const int pdg{std::abs(pLArMC->GetParticleId())};

    // Assumes we've already ruled out diffuse and shower-like EM particles
    return pdg == PHOTON || pdg == E_MINUS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::UndergoesBremmstrahlung(const LArMCParticle *const pRoot, const LArMCParticleHelper::MCContributionMap &mcToHitsMap,
    bool onBremmBranch) const
{
    if (!pRoot)
        return false;
    MCProcess process{pRoot->GetProcess()};
    const CaloHitList &hits{mcToHitsMap.find(pRoot) != mcToHitsMap.end() ? mcToHitsMap.at(pRoot) : CaloHitList()};
    if ((!hits.empty() && process == MC_PROC_E_BREM) || (onBremmBranch && !hits.empty()))
        return true;

    const MCParticleList &childList{pRoot->GetDaughterList()};
    for (const MCParticle *const pChild : childList)
    {
        if (this->UndergoesBremmstrahlung(dynamic_cast<const LArMCParticle *>(pChild), mcToHitsMap, onBremmBranch || process == MC_PROC_E_BREM))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::AllHitsCompton(const LArMCParticle *const pRoot, const LArMCParticleHelper::MCContributionMap &mcToHitsMap,
    bool onBremmBranch) const
{
    if (!pRoot)
        return false;
    MCProcess process{pRoot->GetProcess()};
    const CaloHitList &hits{mcToHitsMap.find(pRoot) != mcToHitsMap.end() ? mcToHitsMap.at(pRoot) : CaloHitList()};
    if ((!hits.empty() && process != MC_PROC_COMPT) || (onBremmBranch && !hits.empty()))
        return false;

    const MCParticleList &childList{pRoot->GetDaughterList()};
    for (const MCParticle *const pChild : childList)
    {
        if (!this->AllHitsCompton(dynamic_cast<const LArMCParticle *>(pChild), mcToHitsMap, onBremmBranch || process == MC_PROC_E_BREM))
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArMCParticle *DlHitTrackShowerIdAlgorithm::AllocateHitOwner(const pandora::MCParticle *const pParticle, const MCFoldingMap &mcFoldingMap,
    const LArMCParticleHelper::MCContributionMap &mcHitMap, CaloHitList &particleOwnedHits, CaloHitList &parentOwnedHits) const
{
    const MCParticleList &parentList{mcFoldingMap.at(pParticle)->GetParentList()};
    if (parentList.empty() || mcFoldingMap.find(parentList.front()) == mcFoldingMap.end())
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    const LArMCParticle *pParent{parentList.empty() ? nullptr : dynamic_cast<const LArMCParticle *>(mcFoldingMap.at(parentList.front()))};
    if (pParent == pParticle)
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    while (mcHitMap.find(pParent) == mcHitMap.end())
    {
        pParent = pParent->GetParentList().empty() ? nullptr : dynamic_cast<const LArMCParticle *>(mcFoldingMap.at(pParent->GetParentList().front()));
        if (!pParent)
            return nullptr;
    }
    if (!(this->IsHip(dynamic_cast<const LArMCParticle *>(pParent)) || this->IsMip(dynamic_cast<const LArMCParticle *>(pParent))))
        return nullptr;
    const CaloHitList &particleHits{mcHitMap.at(pParticle)};
    CaloHitList parentHits(mcHitMap.at(pParent).begin(), mcHitMap.at(pParent).end());
    CartesianVector centroid(0, 0, 0);
    LArPcaHelper::EigenValues eigenValues(0, 0, 0);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(parentHits, centroid, eigenValues, eigenVectors);
    FloatVector projections;
    for (const CaloHit *const pCaloHit : parentHits)
    {
        const CartesianVector dir{pCaloHit->GetPositionVector() - centroid};
        projections.emplace_back(dir.GetDotProduct(eigenVectors[0]));
    }
    std::vector<size_t> indices(projections.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&projections](size_t i, size_t j) { return projections[i] < projections[j]; });

    CaloHitVector parentHitsVector(parentHits.begin(), parentHits.end());
    std::vector<std::pair<size_t, const CaloHit *>> paired;
    for (size_t i = 0; i < indices.size(); ++i)
        paired.emplace_back(indices[i], parentHitsVector[i]);
    std::sort(paired.begin(), paired.end(), [](const std::pair<size_t, const CaloHit *> &a, const std::pair<size_t, const CaloHit *> &b)
    {
        return a.first > b.first;
    });

    for (const CaloHit *const pCaloHit : particleHits)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        const double xmin{pos.GetX() - 0.5 * pCaloHit->GetCellSize1()}, xmax{pos.GetX() + 0.5 * pCaloHit->GetCellSize1()};
        const double zmin{pos.GetZ() - 0.5 * pCaloHit->GetCellSize0()}, zmax{pos.GetZ() + 0.5 * pCaloHit->GetCellSize0()};
        bool hasInterveningHit{false};
        for (size_t i = 0; i < paired.size() - 1; ++i)
        {
            const CaloHit *const pParentHit1{paired[i].second};
            const CartesianVector &pos1{pParentHit1->GetPositionVector()};
            const CaloHit *const pParentHit2{paired[i + 1].second};
            const CartesianVector &pos2{pParentHit2->GetPositionVector()};
            if (pCaloHit == pParentHit1 || pCaloHit == pParentHit2)
                continue;

            double entry{0.}, exit{1.};
            auto check_axis = [&](double p1, double p2, double minB, double maxB)
            {
                double t1{(minB - p1) / (p2 - p1)};
                double t2{(maxB - p1) / (p2 - p1)};

                if (t1 > t2)
                    std::swap(t1, t2);

                entry = std::max(entry, t1);
                exit = std::min(exit, t2);

                return entry <= exit;
            };

            if (check_axis(pos1.GetX(), pos2.GetX(), xmin, xmax) && check_axis(pos1.GetZ(), pos2.GetZ(), zmin, zmax))
            {
                hasInterveningHit = true;
                break;
            }
        }
        if (!hasInterveningHit)
        {
            particleOwnedHits.emplace_back(pCaloHit);
        }
        else
        {
            parentOwnedHits.emplace_back(pCaloHit);
        }
    }

    return pParent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArMCParticle *DlHitTrackShowerIdAlgorithm::AllocateHitOwner(const pandora::MCParticle *const pParticle,
    const LArMCParticleHelper::MCContributionMap &mcHitMap, CaloHitList &particleOwnedHits, CaloHitList &parentOwnedHits) const
{
    const MCParticleList &parentList{pParticle->GetParentList()};
    if (parentList.empty())
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    const LArMCParticle *pParent{dynamic_cast<const LArMCParticle *>(parentList.front())};
    if (!(this->IsHip(pParent) || this->IsMip(pParent)))
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    const CaloHitList &particleHits{mcHitMap.at(pParticle)};
    // Not guaranteed that a parent particle has hits, and so may not be in the hit map
    if (mcHitMap.find(pParent) == mcHitMap.end())
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    CaloHitList parentHits(mcHitMap.at(pParent).begin(), mcHitMap.at(pParent).end());
    CartesianVector centroid(0, 0, 0);
    LArPcaHelper::EigenValues eigenValues(0, 0, 0);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(parentHits, centroid, eigenValues, eigenVectors);
    FloatVector projections;
    for (const CaloHit *const pCaloHit : parentHits)
    {
        const CartesianVector dir{pCaloHit->GetPositionVector() - centroid};
        projections.emplace_back(dir.GetDotProduct(eigenVectors[0]));
    }
    std::vector<size_t> indices(projections.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&projections](size_t i, size_t j) { return projections[i] < projections[j]; });

    CaloHitVector parentHitsVector(parentHits.begin(), parentHits.end());
    std::vector<std::pair<size_t, const CaloHit *>> paired;
    for (size_t i = 0; i < indices.size(); ++i)
        paired.emplace_back(indices[i], parentHitsVector[i]);
    std::sort(paired.begin(), paired.end(), [](const std::pair<size_t, const CaloHit *> &a, const std::pair<size_t, const CaloHit *> &b)
    {
        return a.first > b.first;
    });

    for (const CaloHit *const pCaloHit : particleHits)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        const double xmin{pos.GetX() - 0.5 * pCaloHit->GetCellSize1()}, xmax{pos.GetX() + 0.5 * pCaloHit->GetCellSize1()};
        const double zmin{pos.GetZ() - 0.5 * pCaloHit->GetCellSize0()}, zmax{pos.GetZ() + 0.5 * pCaloHit->GetCellSize0()};
        bool hasInterveningHit{false};
        for (size_t i = 0; i < paired.size() - 1; ++i)
        {
            const CaloHit *const pParentHit1{paired[i].second};
            const CartesianVector &pos1{pParentHit1->GetPositionVector()};
            const CaloHit *const pParentHit2{paired[i + 1].second};
            const CartesianVector &pos2{pParentHit2->GetPositionVector()};
            if (pCaloHit == pParentHit1 || pCaloHit == pParentHit2)
                continue;

            double entry{0.}, exit{1.};
            auto check_axis = [&](double p1, double p2, double minB, double maxB)
            {
                double t1{(minB - p1) / (p2 - p1)};
                double t2{(maxB - p1) / (p2 - p1)};

                if (t1 > t2)
                    std::swap(t1, t2);

                entry = std::max(entry, t1);
                exit = std::min(exit, t2);

                return entry <= exit;
            };

            if (check_axis(pos1.GetX(), pos2.GetX(), xmin, xmax) && check_axis(pos1.GetZ(), pos2.GetZ(), zmin, zmax))
            {
                hasInterveningHit = true;
                break;
            }
        }
        if (!hasInterveningHit)
        {
            particleOwnedHits.emplace_back(pCaloHit);
        }
        else
        {
            parentOwnedHits.emplace_back(pCaloHit);
        }
    }

    return pParent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::GetVertexPlanarCoordinates(const HitType view, float &vx, float &vz) const
{
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    if (pVertexList->empty())
    {
        std::cout << "DlHitTrackShowerIdAlgorithm: No vertices found in the vertex list" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }
    const CartesianVector &vertex{pVertexList->front()->GetPosition()};
    vx = vertex.GetX();
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    switch (view)
    {
        case TPC_VIEW_U:
            vz = transform->YZtoU(vertex.GetY(), vertex.GetZ());
            break;
        case TPC_VIEW_V:
            vz = transform->YZtoV(vertex.GetY(), vertex.GetZ());
            break;
        case TPC_VIEW_W:
            vz = transform->YZtoW(vertex.GetY(), vertex.GetZ());
            break;
        default:
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::GetCoordinateExtrema(const CaloHitList &caloHitList, const float vx, const float vz, Bounds &bounds) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        bounds.xMin = std::min(pos.GetX(), bounds.xMin);
        bounds.xMax = std::max(pos.GetX(), bounds.xMax);
        bounds.zMin = std::min(pos.GetZ(), bounds.zMin);
        bounds.zMax = std::max(pos.GetZ(), bounds.zMax);
        const float x{pos.GetX() - vx};
        const float z{pos.GetZ() - vz};
        const float r{std::sqrt(x * x + z * z)};
        bounds.rMax = std::max(r, bounds.rMax);
    }
    bounds.xRange = bounds.xMax - bounds.xMin > 0 ? bounds.xMax - bounds.xMin : 1.f;
    bounds.zRange = bounds.zMax - bounds.zMin > 0 ? bounds.zMax - bounds.zMin : 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::PopulateInputVectors(const CaloHitList &caloHitList, const Bounds &bounds, const float vx, const float vz,
    FloatVector &x_rel, FloatVector &z_rel, FloatVector &x_abs, FloatVector &z_abs, FloatVector &rr, FloatVector &cosTheta, FloatVector &sinTheta,
    FloatVector &wirePitch, FloatVector &wireAngle, FloatVector &adc, FloatVector &width, CaloHitVector &sortedCaloHitList) const
{
    if (caloHitList.empty())
        return;

    LArGeometryHelper::DetectorBoundaries detBounds{LArGeometryHelper::GetDetectorBoundaries(this->GetPandora())};
    const float detSpanX{detBounds.m_xBoundaries.second - detBounds.m_xBoundaries.first};
    const float detSpanZ{detBounds.m_zBoundaries.second - detBounds.m_zBoundaries.first};
    const LArTPC *const pLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);

    const size_t nHits{caloHitList.size()};

    x_rel.resize(nHits, 0);
    z_rel.resize(nHits, 0);
    x_abs.resize(nHits, 0);
    z_abs.resize(nHits, 0);
    rr.resize(nHits, 0);
    cosTheta.resize(nHits, 0);
    sinTheta.resize(nHits, 0);
    wirePitch.resize(nHits, 0);
    wireAngle.resize(nHits, 0);
    adc.resize(nHits, 0);
    width.resize(nHits, 0);

    size_t i{0};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        const HitType view{pCaloHit->GetHitType()};
        switch (view)
        {
            case TPC_VIEW_U:
                wirePitch[i] = pLArTPC->GetWirePitchU();
                wireAngle[i] = pLArTPC->GetWireAngleU();
                break;
            case TPC_VIEW_V:
                wirePitch[i] = pLArTPC->GetWirePitchV();
                wireAngle[i] = pLArTPC->GetWireAngleV();
                break;
            case TPC_VIEW_W:
                wirePitch[i] = pLArTPC->GetWirePitchW();
                wireAngle[i] = pLArTPC->GetWireAngleW();
                break;
            default:
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }
        const float dx{x - vx};
        const float dz{z - vz};
        x_rel[i] = dx / bounds.rMax;
        z_rel[i] = dz / bounds.rMax;
        x_abs[i] = x / detSpanX;
        z_abs[i] = z / detSpanZ;

        const float r{std::sqrt(dx * dx + dz * dz)};
        if (r > 0.f)
        {
            rr[i] = r / bounds.rMax;
            // Scale to [0,1]
            cosTheta[i] = dx / (r + std::numeric_limits<float>::epsilon());
            sinTheta[i] = dz / (r + std::numeric_limits<float>::epsilon());
        }
        else
        {
            rr[i] = 0.f;
            cosTheta[i] = 0;
            sinTheta[i] = 0;
        }

        adc[i] = std::log(1.f + pCaloHit->GetInputEnergy()) * 0.1; // Log scale and nomralization to keep expected values in [0,1]
        width[i] = pCaloHit->GetCellSize1() / detSpanX;
        ++i;
    }

    // Sort the vectors by r
    IndexVector indices(adc.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&rr](size_t a, size_t b) { return rr[a] < rr[b]; });

    auto reorder = [&](auto &vec)
    {
        using VecType = std::decay_t<decltype(vec)>;
        VecType sorted(vec.size());
        for (size_t a = 0; a < indices.size(); ++a)
            sorted[a] = vec[indices[a]];
        vec.swap(sorted);
    };

    reorder(x_rel);
    reorder(z_rel);
    reorder(x_abs);
    reorder(z_abs);
    reorder(wirePitch);
    reorder(wireAngle);
    reorder(rr);
    reorder(cosTheta);
    reorder(sinTheta);
    reorder(adc);
    reorder(width);

    sortedCaloHitList.clear();
    sortedCaloHitList.insert(sortedCaloHitList.end(), caloHitList.begin(), caloHitList.end());
    reorder(sortedCaloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::TraverseChildren(const MCParticle *const pRoot, const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const std::string &tab) const
{
    const MCParticleList &childList{pRoot->GetDaughterList()};
    for (const MCParticle *const pChild : childList)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pChild)};
        if (!pLArMC)
            continue;
        MCProcess process{pLArMC->GetProcess()};
        const int nHits{mcToHitsMap.find(pChild) != mcToHitsMap.end() ? static_cast<int>(mcToHitsMap.at(pChild).size()) : 0};
        std::cout << (tab + "  ") << "PDG: " << pLArMC->GetParticleId() << " Proc: " << process << " Hits: " << nHits << std::endl;

        this->TraverseChildren(pChild, mcToHitsMap, tab + "  ");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeqLen", m_maxSeqLen));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_model);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
