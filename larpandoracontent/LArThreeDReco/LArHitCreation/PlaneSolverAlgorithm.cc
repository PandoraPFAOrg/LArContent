/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverAlgorithm.h"

#include <ranges>

using namespace pandora;

namespace lar_content
{

PlaneSolverAlgorithm::PlaneSolverAlgorithm() :
    m_chi2Threshold(6.f),
    m_caloHitListName{"CaloHitList2D"},
    m_eventContextName{"PlaneContext"},
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PlaneSolverAlgorithm::Run()
{
    m_planeToReadoutMap.clear();
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    this->FillHitMap(*pCaloHitList);
    this->Solve();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::FillHitMap(const CaloHitList &caloHitList)
{
    for (const auto *pCaloHit : caloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));
        if (!pLArCaloHit)
        {
            std::cout << "PlaneSolverAlgorithm - non LArCaloHit found in calo hit list: " << m_caloHitListName << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
        const unsigned int tpc{pLArCaloHit->GetLArTPCVolumeId()};
        const unsigned int childVolume{pLArCaloHit->GetDaughterVolumeId()};
        const Volume volume{tpc, childVolume};

        const HitType view(pLArCaloHit->GetHitType());
        m_planeToReadoutMap[volume][view].emplace_back(pLArCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::BinHitsByDrift(const PlaneToHitsMap &hitsByPlane, DriftBinToHitsMap &driftBinToHitsMap) const
{
    driftBinToHitsMap.clear();
    float xMin{std::numeric_limits<float>::max()}, xMax{std::numeric_limits<float>::lowest()};
    for (const auto &[view, hits] : hitsByPlane)
    {
        for (const auto *pHit : hits)
        {
            const float x{pHit->GetPositionVector().GetX()};
            xMin = std::min(xMin, x);
            xMax = std::max(xMax, x);
        }
    }

    // Adjust the minimum and maximum to sit on 0.5 cm boundaries, and ensure that the bin region extends beyond the range of hits
    float xLowBin{std::floor(xMin * 2.f) / 2.f};
    float xHighBin{std::ceil(xMax * 2.f) / 2.f};
    xHighBin = xHighBin == xMax ? xHighBin + 0.5f : xHighBin;
    const int N{static_cast<int>((xHighBin - xLowBin) * 2.f)};

    for (int i = 0, b = 0; i < N; i += 3, ++b)
    {
        const float xBinLow{xLowBin + 0.5f * i};
        const float xBinHigh{xLowBin + 0.5f * (i + 3)};
        for (const auto &[view, hits] : hitsByPlane)
        {
            for (const auto *pHit : hits)
            {
                const float x{pHit->GetPositionVector().GetX()};
                if (x >= xBinLow && x < xBinHigh)
                    driftBinToHitsMap[b][view].emplace_back(pHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::Solve() const
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }
    LArPlaneContextObject *pPlaneContextObject{new LArPlaneContextObject()};
    // Loop over each read out volume (e.g. an APA in a horizontal drift detector) and solve for the optimal set of triplet and doublet
    // relationships between the 2D hits in each unit.
    for (const auto &[_, readout] : m_planeToReadoutMap)
    {
        if (readout.size() < 3)
            continue;

        DriftBinToHitsMap driftBinToHitsMap;
        this->BinHitsByDrift(readout, driftBinToHitsMap);
        for (const auto &[_, driftBinHits] : driftBinToHitsMap)
        {
            if (driftBinHits.size() < 3)
                continue;

            // Set the used view to an "invalid" type for the first pass
            HitType usedView{HIT_CUSTOM};
            CaloHitSet usedHits;
            TripletVector fallbackTriplets;
            for (int pass = 1; pass <= 2; ++pass)
            {
                int nHitsU{static_cast<int>(usedView != TPC_VIEW_U ? driftBinHits.at(TPC_VIEW_U).size() : 0)};
                int nHitsV{static_cast<int>(usedView != TPC_VIEW_V ? driftBinHits.at(TPC_VIEW_V).size() : 0)};
                int nHitsW{static_cast<int>(usedView != TPC_VIEW_W ? driftBinHits.at(TPC_VIEW_W).size() : 0)};
                if (pass == 2 && nHitsU != 0)
                {
                    nHitsU = 0;
                    for (const auto *pHit : driftBinHits.at(TPC_VIEW_U))
                        nHitsU += !usedHits.count(pHit);
                }
                if (pass == 2 && nHitsV != 0)
                {
                    nHitsV = 0;
                    for (const auto *pHit : driftBinHits.at(TPC_VIEW_V))
                        nHitsV += !usedHits.count(pHit);
                }
                if (pass == 2 && nHitsW != 0)
                {
                    nHitsW = 0;
                    for (const auto *pHit : driftBinHits.at(TPC_VIEW_W))
                        nHitsW += !usedHits.count(pHit);
                }

                HitType constraintView{HIT_CUSTOM};
                if (pass == 1)
                {
                    constraintView = nHitsW >= nHitsV && nHitsW >= nHitsU ? TPC_VIEW_W : nHitsV >= nHitsU ? TPC_VIEW_V : TPC_VIEW_U;
                }
                else
                {
                    if (usedView == TPC_VIEW_U)
                        constraintView = nHitsW <= nHitsV ? TPC_VIEW_W : TPC_VIEW_V;
                    else if (usedView == TPC_VIEW_V)
                        constraintView = nHitsW <= nHitsU ? TPC_VIEW_W : TPC_VIEW_U;
                    else if (usedView == TPC_VIEW_W)
                        constraintView = nHitsV <= nHitsU ? TPC_VIEW_V : TPC_VIEW_U;
                }

                CostMatrix costMatrix{this->ComputeCostMatrix(driftBinHits, 100.f, constraintView, usedHits)};
                const IntVector assignment{this->KuhnMunkres(costMatrix)};
                HitType viewA, viewB;
                this->SelectViewPair(constraintView, viewA, viewB);
                int nHitsA{static_cast<int>(driftBinHits.at(viewA).size())};
                int nHitsB{static_cast<int>(driftBinHits.at(viewB).size())};
                if (nHitsA == 0 || nHitsB == 0)
                    continue;
                const PairVector pairs{this->BuildPairs(assignment, nHitsA, nHitsB, costMatrix, m_chi2Threshold)};
                CostMatrix tripletCostMatrix{this->ComputeTripletCostMatrix(pairs, driftBinHits, 100.f, constraintView, usedHits)};
                const IntVector tripletAssignment{this->KuhnMunkres(tripletCostMatrix)};
                int nHitsC{static_cast<int>(driftBinHits.at(constraintView).size())};
                const TripletVector triplets{this->BuildTriplets(pairs, tripletAssignment, nHitsC, tripletCostMatrix, constraintView, m_chi2Threshold)};
                for (size_t i = 0; i < triplets.size(); ++i)
                {
                    if (triplets[i].m_uIndex >= 0 && triplets[i].m_vIndex >= 0 && triplets[i].m_wIndex >= 0)
                    {
                        const CaloHit *pHitU{driftBinHits.at(TPC_VIEW_U)[triplets[i].m_uIndex]};
                        const CaloHit *pHitV{driftBinHits.at(TPC_VIEW_V)[triplets[i].m_vIndex]};
                        const CaloHit *pHitW{driftBinHits.at(TPC_VIEW_W)[triplets[i].m_wIndex]};

                        if (pass == 1 && (usedHits.count(pHitU) || usedHits.count(pHitV) || usedHits.count(pHitW)))
                            continue;
                        if (pass == 2)
                        {
                            // If the matched constraint hit is already used, we can still use the triplet if the other two hits are unused
                            if (constraintView == TPC_VIEW_U && usedHits.count(pHitU))
                                pHitU = nullptr;
                            else if (constraintView == TPC_VIEW_V && usedHits.count(pHitV))
                                pHitV = nullptr;
                            else if (constraintView == TPC_VIEW_W && usedHits.count(pHitW))
                                pHitW = nullptr;

                            if ((pHitU && usedHits.count(pHitU)) || (pHitV && usedHits.count(pHitV)) || (pHitW && usedHits.count(pHitW)))
                                continue;
                        }
                        if (pass == 1 && pPlaneContextObject->AddHitTriplet(pHitU, pHitV, pHitW))
                        {
                            usedHits.insert(pHitU);
                            usedHits.insert(pHitV);
                            usedHits.insert(pHitW);
                        }
                    }
                    else
                    {
                        const CaloHit *pHitU{(triplets[i].m_uIndex >= 0) ? driftBinHits.at(TPC_VIEW_U)[triplets[i].m_uIndex] : nullptr};
                        const CaloHit *pHitV{(triplets[i].m_vIndex >= 0) ? driftBinHits.at(TPC_VIEW_V)[triplets[i].m_vIndex] : nullptr};
                        const CaloHit *pHitW{(triplets[i].m_wIndex >= 0) ? driftBinHits.at(TPC_VIEW_W)[triplets[i].m_wIndex] : nullptr};
                        if ((pHitU && usedHits.count(pHitU)) || (pHitV && usedHits.count(pHitV)) || (pHitW && usedHits.count(pHitW)))
                            continue;
                        if (pass == 1)
                            fallbackTriplets.push_back({triplets[i].m_uIndex, triplets[i].m_vIndex, triplets[i].m_wIndex, triplets[i].m_cost});
                        if (pass == 2 && pPlaneContextObject->AddHitTriplet(pHitU, pHitV, pHitW))
                        {
                            usedHits.insert(pHitU);
                            usedHits.insert(pHitV);
                            usedHits.insert(pHitW);
                        }
                    }
                }
                usedView = constraintView;
            }
            // Add in the fallback triplets that were not picked up in the second pass
            for (const auto &triplet : fallbackTriplets)
            {
                const CaloHit *pHitU{(triplet.m_uIndex >= 0) ? driftBinHits.at(TPC_VIEW_U)[triplet.m_uIndex] : nullptr};
                const CaloHit *pHitV{(triplet.m_vIndex >= 0) ? driftBinHits.at(TPC_VIEW_V)[triplet.m_vIndex] : nullptr};
                const CaloHit *pHitW{(triplet.m_wIndex >= 0) ? driftBinHits.at(TPC_VIEW_W)[triplet.m_wIndex] : nullptr};
                if ((pHitU && usedHits.count(pHitU)) || (pHitV && usedHits.count(pHitV)) || (pHitW && usedHits.count(pHitW)))
                    continue;

                pPlaneContextObject->AddHitTriplet(pHitU, pHitV, pHitW);
            }
        }
    }
    PandoraContentApi::AddEventContextObject(*this, m_eventContextName, pPlaneContextObject);

    if (m_visualize)
    {
        CaloHitList hits3D;
        for (const LArPlaneContextObject::HitTriplet *pHitTriplet : *pPlaneContextObject)
        {
            const CaloHit *pHitU{pHitTriplet->m_uHit};
            const CaloHit *pHitV{pHitTriplet->m_vHit};
            const CaloHit *pHitW{pHitTriplet->m_wHit};
            const CaloHit *pCaloHit3D{nullptr};
            this->CreateThreeDHit(pHitU, pHitV, pHitW, pCaloHit3D);
            if (pCaloHit3D)
                hits3D.push_back(pCaloHit3D);
        }

        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits3D, "3D", BLACK));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::CostMatrix PlaneSolverAlgorithm::ComputeCostMatrix(
    const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost, const HitType constraintView, const CaloHitSet &usedHits) const
{
    HitType viewA, viewB;
    this->SelectViewPair(constraintView, viewA, viewB);
    const CaloHitVector &aHits(planeToHitsMap.at(viewA));
    const CaloHitVector &bHits(planeToHitsMap.at(viewB));
    const CaloHitVector &cHits(planeToHitsMap.at(constraintView));
    const int nA{static_cast<int>(aHits.size())};
    const int nB{static_cast<int>(bHits.size())};
    const int nC{static_cast<int>(cHits.size())};
    const int N{std::max(nA, nB)};

    // Initialize the cost matrix with the cost for unmatched hits.
    CostMatrix C(N, FloatVector(N, unmatchedCost));

    // Compute all chi-squared values
    for (int i = 0; i < nA; ++i)
    {
        if (usedHits.count(aHits[i]))
            continue;

        for (int j = 0; j < nB; ++j)
        {
            if (usedHits.count(bHits[j]))
                continue;

            float bestChi2 = std::numeric_limits<float>::max();
            int bestK{-1};

            for (int k = 0; k < nC; ++k)
            {
                // We allow the constraint hit to be used for cost calculation even if it's already used (it can validate a doublet match)
                float chi2{std::numeric_limits<float>::max()};
                switch (constraintView)
                {
                    case TPC_VIEW_W:
                        chi2 = LArGeometryHelper::CalculateChiSquared(this->GetPandora(), aHits[i], bHits[j], cHits[k]);
                        break;
                    case TPC_VIEW_U:
                        chi2 = LArGeometryHelper::CalculateChiSquared(this->GetPandora(), cHits[k], aHits[i], bHits[j]);
                        break;
                    case TPC_VIEW_V:
                        chi2 = LArGeometryHelper::CalculateChiSquared(this->GetPandora(), aHits[i], cHits[k], bHits[j]);
                        break;
                    default:
                        PANDORA_THROW(STATUS_CODE_INVALID_PARAMETER);
                }
                if (chi2 < bestChi2)
                {
                    bestChi2 = chi2;
                    bestK = k;
                }
            }

            if (bestK >= 0)
                C[i][j] = bestChi2;
        }
    }

    return C;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::CostMatrix PlaneSolverAlgorithm::ComputeTripletCostMatrix(const PairVector &pairs,
    const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost, const HitType constraintView, const CaloHitSet &usedHits) const
{
    HitType viewA, viewB;
    this->SelectViewPair(constraintView, viewA, viewB);
    const CaloHitVector &aHits{planeToHitsMap.at(viewA)};
    const CaloHitVector &bHits{planeToHitsMap.at(viewB)};
    const CaloHitVector &cHits{planeToHitsMap.at(constraintView)};
    int nPairs{static_cast<int>(pairs.size())};
    int nC{static_cast<int>(cHits.size())};
    int N{std::max(nPairs, nC)};

    CostMatrix C(N, FloatVector(N, unmatchedCost));
    for (int p = 0; p < nPairs; ++p)
    {
        int i{pairs[p].m_aIndex};
        int j{pairs[p].m_bIndex};
        if (usedHits.count(aHits[i]) || usedHits.count(bHits[j]))
            continue;

        for (int k = 0; k < nC; ++k)
        {
            float chi2{std::numeric_limits<float>::max()};
            switch (constraintView)
            {
                case TPC_VIEW_W:
                    chi2 = LArGeometryHelper::CalculateChiSquared(this->GetPandora(), aHits[i], bHits[j], cHits[k]);
                    break;
                case TPC_VIEW_U:
                    chi2 = LArGeometryHelper::CalculateChiSquared(this->GetPandora(), cHits[k], aHits[i], bHits[j]);
                    break;
                case TPC_VIEW_V:
                    chi2 = LArGeometryHelper::CalculateChiSquared(this->GetPandora(), aHits[i], cHits[k], bHits[j]);
                    break;
                default:
                    PANDORA_THROW(STATUS_CODE_INVALID_PARAMETER);
            }

            C[p][k] = chi2;
        }
    }

    return C;
}

//------------------------------------------------------------------------------------------------------------------------------------------

IntVector PlaneSolverAlgorithm::KuhnMunkres(const CostMatrix &cost) const
{
    const int N{static_cast<int>(cost.size())};
    // Define the dual variables for rows (u) and columns (v), enforcing u[i] + v[j] <= cost[i][j]
    // We are essentially trying to maximise the sum of the duals over all rows/columns, while enforcing the constraint above
    FloatVector u(N + 1), v(N + 1);
    // Define the matching from a row to a column (p[col] = row) and the way matrix to reconstruct the augmenting path - indexing is
    // 1-based, zero is a dummy node
    IntVector p(N + 1), way(N + 1);

    // Iterate over the rows
    for (int i = 1; i <= N; ++i)
    {
        // Start from the dummy column zero and initialise lowest cost and element usage
        int j0 = 0;
        p[j0] = i;
        FloatVector minReducedCost(N + 1, std::numeric_limits<float>::max());
        IntVector used(N + 1, 0);

        // Build the augmenting path (Dijkstra-like)
        do
        {
            // Mark the current column (j0) and the current row (i0)
            used[j0] = 1;
            int i0{p[j0]}, j1{0};
            float delta{std::numeric_limits<float>::max()};

            // Evaluate all possible assignments - trying to find a path from the current row to the best column with the minimum reduced cost
            for (int j = 1; j <= N; ++j)
            {
                if (used[j])
                    continue;
                // Key step, computing the reduced cost to maintain feasibility
                // When an optimal edge between hits is found, the reduced cost will be zero
                float reducedCost{cost[i0 - 1][j - 1] - u[i0] - v[j]};
                if (reducedCost < minReducedCost[j])
                {
                    // Keep track of the best cost for reaching column j so far
                    minReducedCost[j] = reducedCost;
                    // ... and how we got there. We reach column j via column j0
                    way[j] = j0;
                }
                if (minReducedCost[j] < delta)
                {
                    // Keep track of the best column to expand to
                    delta = minReducedCost[j];
                    j1 = j;
                }
            }

            // Update the potentials (ensures reduced cost of zero for at least one edge)
            // Opposing deltas also ensure that u_i + v_j = (u_i + delta) + (v_j - delta), so the constraint still holds
            for (int j = 0; j <= N; ++j)
            {
                if (used[j])
                {
                    u[p[j]] += delta;
                    v[j] -= delta;
                }
                else
                {
                    // Reduce the slack for unmatched columns
                    minReducedCost[j] -= delta;
                }
            }

            // Step j0 forward - continue the search
            j0 = j1;
            // Stop when we get to an unmatched column, implying we've found an augmenting path. Otherwise, we have a clash and should look
            // to reassign the previous row (which happens because we've now set j0 to j1 for the next loop iteration)
        } while (p[j0] != 0);

        // Walk back through way to update the matching along the augmenting path
        do
        {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    // p[j] = i means row i is matched to column j, so we need to invert this to get the assignment of columns to rows.
    // Unmatched columns will have p[j] = 0, so we assign -1 in the output.
    IntVector assignment(N, -1);
    for (int j = 1; j <= N; ++j)
        if (p[j] != 0)
            assignment[p[j] - 1] = j - 1;

    return assignment;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::PairVector PlaneSolverAlgorithm::BuildPairs(
    const IntVector &assignment, int nA, int nB, const CostMatrix &costMatrix, float chi2Threshold) const
{
    PairVector pairs;

    for (int i = 0; i < nA; ++i)
    {
        const int j{assignment[i]};

        if (j >= 0 && j < nB && costMatrix[i][j] < chi2Threshold)
            pairs.push_back({i, j, costMatrix[i][j]});
    }

    return pairs;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::TripletVector PlaneSolverAlgorithm::BuildTriplets(const PairVector &pairs, const IntVector &assignment, int nC,
    const CostMatrix &costMatrix, const HitType constraintView, float chi2Threshold) const
{
    TripletVector result;
    const int nPairs{static_cast<int>(pairs.size())};

    for (int p = 0; p < nPairs; ++p)
    {
        const int k{assignment[p]};

        if (k >= 0 && k < nC && costMatrix[p][k] < chi2Threshold)
        {
            switch (constraintView)
            {
                case TPC_VIEW_U:
                    result.push_back({k, pairs[p].m_aIndex, pairs[p].m_bIndex, costMatrix[p][k]});
                    break;
                case TPC_VIEW_V:
                    result.push_back({pairs[p].m_aIndex, k, pairs[p].m_bIndex, costMatrix[p][k]});
                    break;
                case TPC_VIEW_W:
                    result.push_back({pairs[p].m_aIndex, pairs[p].m_bIndex, k, costMatrix[p][k]});
                    break;
                default:
                    PANDORA_THROW(STATUS_CODE_INVALID_PARAMETER);
            }
        }
        else if (pairs[p].m_cost < chi2Threshold)
        {
            switch (constraintView)
            {
                case TPC_VIEW_U:
                    result.push_back({-1, pairs[p].m_aIndex, pairs[p].m_bIndex, pairs[p].m_cost});
                    break;
                case TPC_VIEW_V:
                    result.push_back({pairs[p].m_aIndex, -1, pairs[p].m_bIndex, pairs[p].m_cost});
                    break;
                case TPC_VIEW_W:
                    result.push_back({pairs[p].m_aIndex, pairs[p].m_bIndex, -1, pairs[p].m_cost});
                    break;
                default:
                    PANDORA_THROW(STATUS_CODE_INVALID_PARAMETER);
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::CreateThreeDHit(
    const CaloHit *const pCaloHitU, const CaloHit *const pCaloHitV, const CaloHit *const pCaloHitW, const CaloHit *&pCaloHit3D) const
{
    const int nValidHits = (pCaloHitU ? 1 : 0) + (pCaloHitV ? 1 : 0) + (pCaloHitW ? 1 : 0);
    if (nValidHits < 2)
    {
        pCaloHit3D = nullptr;
        return;
    }
    float chi2{0.f};
    CartesianVector pos3D(0, 0, 0);
    PandoraContentApi::CaloHit::Parameters parameters;
    if (pCaloHitU && pCaloHitV && pCaloHitW)
    {
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, pCaloHitU->GetPositionVector(),
            pCaloHitV->GetPositionVector(), pCaloHitW->GetPositionVector(), pos3D, chi2);
    }
    else
    {
        const CaloHit *const pHitA(pCaloHitU ? pCaloHitU : pCaloHitV);
        const CaloHit *const pHitB(pCaloHitW ? pCaloHitW : pCaloHitV);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), pHitA->GetHitType(), pHitB->GetHitType(), pHitA->GetPositionVector(),
            pHitB->GetPositionVector(), pos3D, chi2);
    }

    parameters.m_positionVector = pos3D;
    parameters.m_hitType = TPC_3D;

    const CaloHit *const pCaloHit2D(pCaloHitW ? pCaloHitW : pCaloHitU ? pCaloHitU : pCaloHitV);
    parameters.m_pParentAddress = static_cast<const void *>(pCaloHit2D);
    parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
    parameters.m_cellGeometry = RECTANGULAR;
    parameters.m_cellSize0 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellSize1 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellNormalVector = pCaloHit2D->GetCellNormalVector();
    parameters.m_expectedDirection = pCaloHit2D->GetExpectedDirection();
    parameters.m_nCellRadiationLengths = pCaloHit2D->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit2D->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit2D->GetTime();
    parameters.m_inputEnergy = pCaloHit2D->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit2D->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit2D->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit2D->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit2D->IsDigital();
    parameters.m_hitRegion = pCaloHit2D->GetHitRegion();
    parameters.m_layer = pCaloHit2D->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::SelectViewPair(const HitType constraintView, HitType &viewA, HitType &viewB) const
{
    switch (constraintView)
    {
        case TPC_VIEW_U:
            viewA = TPC_VIEW_V;
            viewB = TPC_VIEW_W;
            break;
        case TPC_VIEW_V:
            viewA = TPC_VIEW_U;
            viewB = TPC_VIEW_W;
            break;
        case TPC_VIEW_W:
            viewA = TPC_VIEW_U;
            viewB = TPC_VIEW_V;
            break;
        default:
            PANDORA_THROW(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PlaneSolverAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Chi2Threshold", m_chi2Threshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventContextName", m_eventContextName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
