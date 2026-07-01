/**
 *  @file   larpandoracontent/LArMetrics/ShowerValidationTool.cc
 *
 *  @brief  Implementation of the shower validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMetrics/ShowerValidationTool.h"

using namespace pandora;

namespace lar_content
{

ShowerValidationTool::ShowerValidationTool() :
    m_trueLengthEnergyFrac(0.90),
    m_initialRegion3D(14.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerValidationTool::Run(const Algorithm *const pAlgorithm, [[maybe_unused]] const MCParticle *const pMCNu,
    const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ShowerTreeVars showerTreeVars;
    showerTreeVars.m_run = this->GetPandora().GetRun();
    showerTreeVars.m_subrun = this->GetPandora().GetSubrun();
    showerTreeVars.m_event = this->GetPandora().GetEvent();

    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMC(targetMC.at(i));
        bool doesMCHaveDir(pMC->GetMomentum().GetMagnitudeSquared() > std::numeric_limits<float>::epsilon());
        const Pfo *const pPfo(bestRecoMatch.at(i));

        // Fill MC-related vars
        if (doesMCHaveDir)
        {
            this->GetTrueLength(mcMatchesVec, pMC, showerTreeVars);
            this->GetInitialRegionVars(mcMatchesVec, pMC, pPfo, showerTreeVars);
        }
        else
        {
            this->FillForNullMCDir(showerTreeVars);
        }

        // Perform shower fit
        float recoShrLength(m_invalidSmallFloat);
        CartesianVector recoShrVtx(m_invalidLargeFloat, m_invalidLargeFloat, m_invalidLargeFloat);
        CartesianVector recoShrDir(m_invalidLargeFloat, m_invalidLargeFloat, m_invalidLargeFloat);
        bool fitSuccess(pPfo && (this->FitShower(pPfo, recoShrVtx, recoShrDir, recoShrLength)));

        if (fitSuccess)
        {
            this->GetRecoVertexInfo(recoShrVtx, recoShrDir, recoShrLength, (doesMCHaveDir ? pMC : nullptr), showerTreeVars);
            this->GetMoliere(pPfo, recoShrVtx, recoShrDir, showerTreeVars);
        }
        else
        {
            this->FillForFailedPfo(showerTreeVars);
        }
    }

    this->FillTree(showerTreeVars);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetTrueLength(
    const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticle *const pMCParticle, ShowerTreeVars &showerTreeVars)
{
    // Get 3D positions/directions
    const CartesianVector trueShrVtx((pMCParticle->GetParticleId() == PHOTON) ? pMCParticle->GetEndpoint() : pMCParticle->GetVertex());
    const CartesianVector trueShrDir(pMCParticle->GetMomentum().GetUnitVector());
    const CartesianVector trueShrDirSeed(trueShrVtx + (trueShrDir * 10.0f));

    // Need to loop over matches to get MCNode
    CaloHitList mcHits;
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() == pMCParticle)
            mcHits = mcMatches.GetMC()->GetCaloHits();
    }

    // Find X% containment in longitudinal
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        float totalEnergy(0.f); // need for function call
        CaloHitVector viewMCHits;
        this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);

        const CartesianVector trueShrVtx2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrVtx, hitType));
        const CartesianVector trueShrDirSeed2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrDirSeed, hitType));
        const CartesianVector trueShrDir2D((trueShrDirSeed2D - trueShrVtx2D).GetUnitVector());

        // Order hits wrt l
        std::sort(viewMCHits.begin(), viewMCHits.end(),
            [&trueShrDir2D, &trueShrVtx2D](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
            {
                const CartesianVector positionA(pCaloHitA->GetPositionVector() - trueShrVtx2D);
                const CartesianVector positionB(pCaloHitB->GetPositionVector() - trueShrVtx2D);

                const float lA(trueShrDir2D.GetDotProduct(positionA));
                const float lB(trueShrDir2D.GetDotProduct(positionB));

                return lA < lB;
            });

        // Calculate 'true length'
        FloatVector &trueLength(hitType == TPC_VIEW_U ? showerTreeVars.m_coreTrueLengthFromU
                : hitType == TPC_VIEW_V               ? showerTreeVars.m_coreTrueLengthFromV
                                                      : showerTreeVars.m_coreTrueLengthFromW);

        float runningEnergySum(0.f), endpointL(m_invalidLargeFloat);
        for (const CaloHit *const pHit2D : viewMCHits)
        {
            const float hitEnergy(std::fabs(pHit2D->GetElectromagneticEnergy()));
            runningEnergySum += hitEnergy;

            if ((totalEnergy > std::numeric_limits<float>::epsilon()) && ((runningEnergySum / totalEnergy) > m_trueLengthEnergyFrac))
            {
                const CartesianVector displacement(pHit2D->GetPositionVector() - trueShrVtx2D);
                endpointL = trueShrDir2D.GetDotProduct(displacement);
                break;
            }
        }

        if (fabs(endpointL - m_invalidLargeFloat) > std::numeric_limits<float>::epsilon())
        {
            const CartesianVector showerEndpoint2D(trueShrVtx2D + (trueShrDir2D * endpointL));
            const float scale3D(trueShrDir.GetX() > std::numeric_limits<float>::epsilon()
                    ? (showerEndpoint2D.GetX() - trueShrVtx2D.GetX()) / trueShrDir.GetX()
                    : trueShrDir.GetY() > std::numeric_limits<float>::epsilon()
                    ? (showerEndpoint2D.GetY() - trueShrVtx2D.GetY()) / trueShrDir.GetY()
                    : (showerEndpoint2D.GetZ() - trueShrVtx2D.GetZ()));
            trueLength.push_back(scale3D);
        }
        else
        {
            trueLength.push_back(m_invalidSmallFloat);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetInitialRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec,
    const MCParticle *const pMCParticle, const Pfo *const pPfo, ShowerTreeVars &showerTreeVars)
{
    // Loop over matches to get MCNode
    CaloHitList mcHits;
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() == pMCParticle)
        {
            mcHits = mcMatches.GetMC()->GetCaloHits();
            break;
        }
    }

    // Get 3D positions/directions
    const CartesianVector trueShrVtx((pMCParticle->GetParticleId() == PHOTON) ? pMCParticle->GetEndpoint() : pMCParticle->GetVertex());
    const CartesianVector trueShrDir(pMCParticle->GetMomentum().GetUnitVector());
    const CartesianVector trueShrDirSeed(trueShrVtx + (trueShrDir * m_initialRegion3D));

    int nTotalInitialMCHits(0), nTotalInitialPfoHits(0), nTotalInitialSharedHits(0);
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const CartesianVector trueShrVtx2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrVtx, hitType));
        const CartesianVector trueShrDirSeed2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrDirSeed, hitType));
        const CartesianVector trueShrDir2D((trueShrDirSeed2D - trueShrVtx2D).GetUnitVector());
        const float initialRegionL(trueShrDir2D.GetDotProduct(trueShrDirSeed2D - trueShrVtx2D));

        // Find initial region MC hits
        float totalEnergy(0.f); // need for function call
        CaloHitVector viewMCHits, initialMCHits;
        this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);
        for (const CaloHit *const pCaloHit : viewMCHits)
        {
            const float l(trueShrDir2D.GetDotProduct(pCaloHit->GetPositionVector() - trueShrVtx2D));

            if (l < initialRegionL)
                initialMCHits.push_back(pCaloHit);
        }
        nTotalInitialMCHits += initialMCHits.size();

        int nInitialPfoHits(0), nSharedHits(0);
        if (pPfo)
        {
            // Find initial region pfo hits
            CaloHitList viewPfoHitList;
            LArPfoHelper::GetCaloHits(pPfo, hitType, viewPfoHitList);
            CaloHitVector viewPfoHits(viewPfoHitList.begin(), viewPfoHitList.end()), initialPfoHits;
            for (const CaloHit *const pCaloHit : viewPfoHits)
            {
                const float l(trueShrDir2D.GetDotProduct(pCaloHit->GetPositionVector() - trueShrVtx2D));

                if (l < initialRegionL)
                    initialPfoHits.push_back(pCaloHit);
            }
            nInitialPfoHits = initialPfoHits.size();
            nTotalInitialPfoHits += nInitialPfoHits;

            // Get shared hits
            for (const CaloHit *const pCaloHit : initialMCHits)
            {
                if (std::find(initialPfoHits.begin(), initialPfoHits.end(), pCaloHit) != initialPfoHits.end())
                    ++nSharedHits;
            }
            nTotalInitialSharedHits += nSharedHits;
        }

        IntVector &viewNInitialMCHits(hitType == TPC_VIEW_U ? showerTreeVars.m_nInitialMCHitsU
                : hitType == TPC_VIEW_V                     ? showerTreeVars.m_nInitialMCHitsV
                                                            : showerTreeVars.m_nInitialMCHitsW);
        IntVector &viewNInitialPfoHits(hitType == TPC_VIEW_U ? showerTreeVars.m_nInitialPfoHitsU
                : hitType == TPC_VIEW_V                      ? showerTreeVars.m_nInitialPfoHitsV
                                                             : showerTreeVars.m_nInitialPfoHitsW);
        FloatVector &viewCompleteness(hitType == TPC_VIEW_U ? showerTreeVars.m_initialCompletenessU
                : hitType == TPC_VIEW_V                     ? showerTreeVars.m_initialCompletenessV
                                                            : showerTreeVars.m_initialCompletenessW);
        FloatVector &viewPurity(hitType == TPC_VIEW_U ? showerTreeVars.m_initialPurityU
                : hitType == TPC_VIEW_V               ? showerTreeVars.m_initialPurityV
                                                      : showerTreeVars.m_initialPurityW);

        const float thisCompleteness(initialMCHits.size() == 0 ? m_invalidSmallFloat : static_cast<float>(nSharedHits) / initialMCHits.size());
        const float thisPurity(nInitialPfoHits == 0 ? m_invalidSmallFloat : static_cast<float>(nSharedHits) / nInitialPfoHits);
        viewNInitialMCHits.push_back(initialMCHits.size());
        viewNInitialPfoHits.push_back(nInitialPfoHits);
        viewCompleteness.push_back(thisCompleteness);
        viewPurity.push_back(thisPurity);
    }

    const float completeness(nTotalInitialMCHits == 0 ? 0.f : static_cast<float>(nTotalInitialSharedHits) / nTotalInitialMCHits);
    const float purity(nTotalInitialPfoHits == 0 ? 0.f : static_cast<float>(nTotalInitialSharedHits) / nTotalInitialPfoHits);

    showerTreeVars.m_nInitialMCHits.push_back(nTotalInitialMCHits);
    showerTreeVars.m_nInitialPfoHits.push_back(nTotalInitialPfoHits);
    showerTreeVars.m_initialCompleteness.push_back(completeness);
    showerTreeVars.m_initialPurity.push_back(purity);
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Replicate PandoraShower fitting
bool ShowerValidationTool::FitShower(const Pfo *const pPfo, CartesianVector &showerVertex, CartesianVector &showerDirection, float &showerLength)
{
    CartesianPointVector positions3D;
    LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, positions3D);
    if (positions3D.empty())
    {
        return false;
    }

    try
    {
        const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(pPfo));
        const CartesianVector vertexPosition(pRecoVertex->GetPosition());
        const LArShowerPCA initialLArShowerPCA(LArPfoHelper::GetPrincipalComponents(positions3D, vertexPosition));
        const CartesianVector &centroid(initialLArShowerPCA.GetCentroid());
        const CartesianVector &primaryAxis(initialLArShowerPCA.GetPrimaryAxis());
        const CartesianVector &secondaryAxis(initialLArShowerPCA.GetSecondaryAxis());
        const CartesianVector &tertiaryAxis(initialLArShowerPCA.GetTertiaryAxis());
        const CartesianVector &eigenvalues(initialLArShowerPCA.GetEigenValues());

        // Project the PFParticle vertex onto the PCA axis
        const CartesianVector projectedVertexPosition(centroid - primaryAxis.GetUnitVector() * (centroid - vertexPosition).GetDotProduct(primaryAxis));

        // By convention, principal axis should always point away from vertex
        const float testProjection(primaryAxis.GetDotProduct(projectedVertexPosition - centroid));
        const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

        const LArShowerPCA larShowerPCA(centroid, primaryAxis * directionScaleFactor, secondaryAxis * directionScaleFactor,
            tertiaryAxis * directionScaleFactor, eigenvalues);

        showerVertex = projectedVertexPosition;
        showerDirection = larShowerPCA.GetPrimaryAxis();
        showerLength = larShowerPCA.GetAxisLengths().GetX();
    }
    catch (...)
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetRecoVertexInfo(const CartesianVector &recoShrVtx, const CartesianVector &recoShrDir,
    const float recoShrLength, const MCParticle *const pMC, ShowerTreeVars &showerTreeVars)
{
    showerTreeVars.m_recoShrVtxX.push_back(recoShrVtx.GetX());
    showerTreeVars.m_recoShrVtxY.push_back(recoShrVtx.GetY());
    showerTreeVars.m_recoShrVtxZ.push_back(recoShrVtx.GetZ());
    showerTreeVars.m_recoShrDirX.push_back(recoShrDir.GetX());
    showerTreeVars.m_recoShrDirY.push_back(recoShrDir.GetY());
    showerTreeVars.m_recoShrDirZ.push_back(recoShrDir.GetZ());
    showerTreeVars.m_recoShrLength.push_back(recoShrLength);

    if (pMC)
    {
        const CartesianVector trueDir(pMC->GetMomentum().GetUnitVector());
        float recoShrDirAcc(m_invalidAngle);
        try { recoShrDirAcc = trueDir.GetOpeningAngle(recoShrDir); } 
        catch (StatusCodeException &) { recoShrDirAcc = m_invalidLargeFloat; }
        showerTreeVars.m_recoShrDirAcc.push_back(recoShrDirAcc);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetMoliere(
    const Pfo *const pPfo, const CartesianVector &showerVertex, const CartesianVector &showerDirection, ShowerTreeVars &showerTreeVars)
{
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHits3D);

    // Get 3D hits created from collection-view hits
    float totalEnergy(0.f);
    CaloHitVector caloHits3DFromW;

    for (const CaloHit *const pHit3D : caloHits3D)
    {
        const CaloHit *const pHit2D{static_cast<const CaloHit *>(pHit3D->GetParentAddress())};

        if (pHit2D->GetHitType() == TPC_VIEW_W)
        {
            totalEnergy += pHit2D->GetElectromagneticEnergy();
            caloHits3DFromW.push_back(pHit3D);
        }
    }

    // Now do Molliere
    float runningEnergySum(0.f), moliereRadius(m_invalidSmallFloat);

    std::sort(caloHits3DFromW.begin(), caloHits3DFromW.end(),
        [&showerDirection, &showerVertex](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
        {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - showerVertex);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - showerVertex);

            const float tA(showerDirection.GetCrossProduct(positionA).GetMagnitude());
            const float tB(showerDirection.GetCrossProduct(positionB).GetMagnitude());

            return tA < tB;
        });

    for (const CaloHit *const pHit3D : caloHits3DFromW)
    {
        const CaloHit *const pHit2D{static_cast<const CaloHit *>(pHit3D->GetParentAddress())};
        const float hitEnergy(std::fabs(pHit2D->GetElectromagneticEnergy()));
        runningEnergySum += hitEnergy;

        if ((totalEnergy > std::numeric_limits<float>::epsilon()) && ((runningEnergySum / totalEnergy) > 0.9f))
        {
            const CartesianVector displacement(pHit3D->GetPositionVector() - showerVertex);
            moliereRadius = showerDirection.GetCrossProduct(displacement).GetMagnitude();
            break;
        }
    }

    showerTreeVars.m_moliereRadius.push_back(moliereRadius);

    // Now do core reco length
    runningEnergySum = 0.f;
    float coreShowerLength = 0.f;

    std::sort(caloHits3DFromW.begin(), caloHits3DFromW.end(),
        [&showerDirection, &showerVertex](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
        {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - showerVertex);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - showerVertex);

            const float lA(showerDirection.GetDotProduct(positionA));
            const float lB(showerDirection.GetDotProduct(positionB));

            return lA < lB;
        });

    for (const CaloHit *const pHit3D : caloHits3DFromW)
    {
        const CaloHit *const pHit2D{static_cast<const CaloHit *>(pHit3D->GetParentAddress())};
        const float hitEnergy(std::fabs(pHit2D->GetElectromagneticEnergy()));
        runningEnergySum += hitEnergy;

        if ((totalEnergy > std::numeric_limits<float>::epsilon()) && ((runningEnergySum / totalEnergy) > m_trueLengthEnergyFrac))
        {
            const CartesianVector displacement(pHit3D->GetPositionVector() - showerVertex);
            coreShowerLength = showerDirection.GetDotProduct(displacement);
            break;
        }
    }

    showerTreeVars.m_coreRecoLength.push_back(coreShowerLength);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::FillForNullMCDir(ShowerTreeVars &showerTreeVars)
{
    showerTreeVars.m_coreTrueLengthFromU.push_back(m_invalidSmallFloat);
    showerTreeVars.m_coreTrueLengthFromV.push_back(m_invalidSmallFloat);
    showerTreeVars.m_coreTrueLengthFromW.push_back(m_invalidSmallFloat);
    showerTreeVars.m_nInitialMCHits.push_back(m_invalidInt);
    showerTreeVars.m_nInitialMCHitsU.push_back(m_invalidInt);
    showerTreeVars.m_nInitialMCHitsV.push_back(m_invalidInt);
    showerTreeVars.m_nInitialMCHitsW.push_back(m_invalidInt);
    showerTreeVars.m_nInitialPfoHits.push_back(m_invalidInt);
    showerTreeVars.m_nInitialPfoHitsU.push_back(m_invalidInt);
    showerTreeVars.m_nInitialPfoHitsV.push_back(m_invalidInt);
    showerTreeVars.m_nInitialPfoHitsW.push_back(m_invalidInt);
    showerTreeVars.m_initialCompleteness.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialCompletenessU.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialCompletenessV.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialCompletenessW.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialPurity.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialPurityU.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialPurityV.push_back(m_invalidSmallFloat);
    showerTreeVars.m_initialPurityW.push_back(m_invalidSmallFloat);
    showerTreeVars.m_recoShrDirAcc.push_back(m_invalidAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::FillForFailedPfo(ShowerTreeVars &showerTreeVars)
{
    showerTreeVars.m_recoShrVtxX.push_back(m_invalidLargeFloat);
    showerTreeVars.m_recoShrVtxY.push_back(m_invalidLargeFloat);
    showerTreeVars.m_recoShrVtxZ.push_back(m_invalidLargeFloat);
    showerTreeVars.m_recoShrDirX.push_back(m_invalidLargeFloat);
    showerTreeVars.m_recoShrDirY.push_back(m_invalidLargeFloat);
    showerTreeVars.m_recoShrDirZ.push_back(m_invalidLargeFloat);
    showerTreeVars.m_recoShrLength.push_back(m_invalidSmallFloat);
    showerTreeVars.m_recoShrDirAcc.push_back(m_invalidAngle);
    showerTreeVars.m_moliereRadius.push_back(m_invalidSmallFloat);
    showerTreeVars.m_coreRecoLength.push_back(m_invalidSmallFloat);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::FillTree(ShowerTreeVars &showerTreeVars)
{
    FloatVector &coreTrueLengthFromU = showerTreeVars.m_coreTrueLengthFromU;
    FloatVector &coreTrueLengthFromV = showerTreeVars.m_coreTrueLengthFromV;
    FloatVector &coreTrueLengthFromW = showerTreeVars.m_coreTrueLengthFromW;
    FloatVector &recoShrVtxX = showerTreeVars.m_recoShrVtxX;
    FloatVector &recoShrVtxY = showerTreeVars.m_recoShrVtxY;
    FloatVector &recoShrVtxZ = showerTreeVars.m_recoShrVtxZ;
    FloatVector &recoShrDirX = showerTreeVars.m_recoShrDirX;
    FloatVector &recoShrDirY = showerTreeVars.m_recoShrDirY;
    FloatVector &recoShrDirZ = showerTreeVars.m_recoShrDirZ;
    FloatVector &coreRecoLength = showerTreeVars.m_coreRecoLength;
    FloatVector &recoShrLength = showerTreeVars.m_recoShrLength;
    FloatVector &recoShrDirAcc = showerTreeVars.m_recoShrDirAcc;
    FloatVector &moliereRadius = showerTreeVars.m_moliereRadius;
    IntVector &nInitialMCHits = showerTreeVars.m_nInitialMCHits;
    IntVector &nInitialMCHitsU = showerTreeVars.m_nInitialMCHitsU;
    IntVector &nInitialMCHitsV = showerTreeVars.m_nInitialMCHitsV;
    IntVector &nInitialMCHitsW = showerTreeVars.m_nInitialMCHitsW;
    IntVector &nInitialPfoHits = showerTreeVars.m_nInitialPfoHits;
    IntVector &nInitialPfoHitsU = showerTreeVars.m_nInitialPfoHitsU;
    IntVector &nInitialPfoHitsV = showerTreeVars.m_nInitialPfoHitsV;
    IntVector &nInitialPfoHitsW = showerTreeVars.m_nInitialPfoHitsW;
    FloatVector &initialCompleteness = showerTreeVars.m_initialCompleteness;
    FloatVector &initialCompletenessU = showerTreeVars.m_initialCompletenessU;
    FloatVector &initialCompletenessV = showerTreeVars.m_initialCompletenessV;
    FloatVector &initialCompletenessW = showerTreeVars.m_initialCompletenessW;
    FloatVector &initialPurity = showerTreeVars.m_initialPurity;
    FloatVector &initialPurityU = showerTreeVars.m_initialPurityU;
    FloatVector &initialPurityV = showerTreeVars.m_initialPurityV;
    FloatVector &initialPurityW = showerTreeVars.m_initialPurityW;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Run", showerTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Subrun", showerTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Event", showerTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_TrueCoreLengthFromU", &coreTrueLengthFromU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_TrueCoreLengthFromV", &coreTrueLengthFromV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_TrueCoreLengthFromW", &coreTrueLengthFromW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoVtxX", &recoShrVtxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoVtxY", &recoShrVtxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoVtxZ", &recoShrVtxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoDirX", &recoShrDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoDirY", &recoShrDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoDirZ", &recoShrDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoCoreLength", &coreRecoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_RecoLength", &recoShrLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_DirAcc", &recoShrDirAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_MoliereRadius", &moliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_InitialMCHits", &nInitialMCHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_InitialMCHitsU", &nInitialMCHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_InitialMCHitsV", &nInitialMCHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "MCP_InitialMCHitsW", &nInitialMCHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPfoHits", &nInitialPfoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPfoHitsU", &nInitialPfoHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPfoHitsV", &nInitialPfoHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPfoHitsW", &nInitialPfoHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialCompleteness", &initialCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialCompletenessU", &initialCompletenessU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialCompletenessV", &initialCompletenessV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialCompletenessW", &initialCompletenessW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPurity", &initialPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPurityU", &initialPurityU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPurityV", &initialPurityV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "BM_InitialPurityW", &initialPurityW));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ShowerTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InitialRegion", m_initialRegion3D));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrueLengthEnergyFrac", m_trueLengthEnergyFrac));

    return BaseValidationTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
