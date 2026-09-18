/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverAlgorithm.h
 *
 *  @brief  The Plane Solver algorithm takes the collection of 2D hits across all planes and constructs an optimal set of triplet and doublet
 *          relationships between them. These relationships are persisted via the EventContext object so that they are available to subsequent
 *          algorithms.
 *
 *  $Log: $
 */
#ifndef LAR_PLANE_SOLVER_ALGORITHM_H
#define LAR_PLANE_SOLVER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  An algorithm to construct hit triplets and doublets representing candidate 3D hits
 */
class PlaneSolverAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PlaneSolverAlgorithm();

private:
    struct Volume
    {
        unsigned int m_tpc;
        unsigned int m_childVolume;

        bool operator==(const Volume &other) const
        {
            return m_tpc == other.m_tpc && m_childVolume == other.m_childVolume;
        }
    };

    struct VolumeHash
    {
        std::size_t operator()(const Volume &v) const noexcept
        {
            return std::hash<int>{}(v.m_tpc) ^ (std::hash<int>{}(v.m_childVolume) << 1);
        }
    };

    struct Pair
    {
        int m_aIndex;
        int m_bIndex;
        float m_cost;
    };

    struct Triplet
    {
        int m_uIndex;
        int m_vIndex;
        int m_wIndex;
        float m_cost;
    };

    typedef std::unordered_map<pandora::HitType, pandora::CaloHitVector> PlaneToHitsMap;
    typedef std::unordered_map<Volume, PlaneToHitsMap, VolumeHash> VolumeToReadoutMap;
    typedef std::unordered_map<int, PlaneToHitsMap> DriftBinToHitsMap;
    typedef std::vector<std::vector<int>> IndexMatrix;
    typedef std::vector<std::vector<float>> CostMatrix;
    typedef std::vector<Pair> PairVector;
    typedef std::vector<Triplet> TripletVector;

    pandora::StatusCode Run();

    /**
     *  @brief  Fills the map of the 2D hits, keyed by volume
     *
     *  @param  caloHitList the list of 2D hits to be mapped
     */
    void FillHitMap(const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Bins the hits by their drift time, to allow for a more efficient matching of hits in the same drift bin.
     *
     *  @param  hitsByPlane the map of 2D hits to be binned, by plane
     *  @param[out] driftBinToHitsMap the map of drift bins to the 2D hits in that bin
     */
    void BinHitsByDrift(const PlaneToHitsMap &hitsByPlane, DriftBinToHitsMap &driftBinToHitsMap) const;

    /**
     *  @brief  Solve for the optimal set of triplet and doublet relationships between the 2D hits
     */
    void Solve() const;

    /**
     *  @brief  Compute the cost matrix for the hits, where the cost is based on the chi-squared value for hit triplets.
     *          This function computes the cost of matching a pair, using the constraint hits to provide a constraint, but allows for the
     *          constraint hits to be used more than once.
     *
     *  @param  planetoHitsMap the map of the 2D hits, keyed by view
     *  @param  unmatchedCost the cost to be assigned to unmatched hits
     *  @param  constraintView the view to be used as the constraint in the chi-squared calculation
     *  @param  usedHits the set of hits that have already been used in a match
     *
     *  @return The cost matrix for the hits in the slice
     */
    CostMatrix ComputeCostMatrix(const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost, const pandora::HitType constraintView,
        const pandora::CaloHitSet &usedHits) const;

    /**
     *  @brief  Compute the cost matrix for the hits, where the cost is based on the chi-squared value for hit triplets.
     *          This function starts with the pairs from the first pass of the Kuhn-Munkres algorithm, and computes the cost of matching
     *          each pair to a unique constraint hit, constructing triplet matches where possible, but retaining the doublet matches otherwise.
     *
     *  @param  pairs the list of pairs to be considered for triplet matching
     *  @param  planetoHitsMap the map of the 2D hits, keyed by view
     *  @param  unmatchedCost the cost to be assigned to unmatched hits
     *  @param  constraintView the view to be used as the constraint in the chi-squared calculation
     *  @param  usedHits the set of hits that have already been used in a match
     *
     *  @return The cost matrix for the hits in the slice
     */
    CostMatrix ComputeTripletCostMatrix(const PairVector &pairs, const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost,
        const pandora::HitType constraintView, const pandora::CaloHitSet &usedHits) const;

    /**
     *  @brief  Implementation of the Kunhne-Munkres (aka Hungarian) algorithm to solve the optimal matching between UV pairs and W hits.
     *          For each row, attempt to assign it and if there is a conflict, reroute via an augmenting path.
     *
     *          Based on Andrey Lopatin's implementation as described at https://cp-algorithms.com/graph/hungarian-algorithm.html
     *
     *          H. W. Kuhn, “The Hungarian Method for the Assignment Problem,” Naval Research Logistics Quarterly, 1955.
     *          J. Munkres, “Algorithms for the Assignment and Transportation Problems,” 1957.
     *
     *  @param  cost the cost matrix for the hits in a volume
     *
     *  @return The optimal matching between UV pairs and W hits, where the vector elements align with the U hits and indicate matched index
     *          in the V hits (or -1 if unmatched)
     */
    pandora::IntVector KuhnMunkres(const CostMatrix &cost) const;

    /**
     *  @brief  Build the list of AB pairs based on the assignment of B hits to A hits from the Kuhn-Munkres algorithm
     *
     *  @param  assignment the output of the Kuhn-Munkres algorithm, indicating the index of the assigned B hit for each A hit (or -1 if unmatched)
     *  @param  nA the number of A hits
     *  @param  nB the number of B hits
     *  @param  costMatrix the cost matrix that was used as input to the Kuhn-Munkres algorithm
     *  @param  chi2Threshold the chi-squared threshold below which a pair is considered a valid match
     *
     *  @return The list of UV pairs
     */
    PairVector BuildPairs(const pandora::IntVector &assignment, int nA, int nB, const CostMatrix &costMatrix, float chi2Threshold) const;

    /**
     *  @brief  Build the list of triplets based on the assignment of constraint hits to AB pairs from the second pass of the Kuhn-Munkres algorithm
     *
     *  @param  pairs the list of AB pairs to be considered for triplet matching
     *  @param  assignment the output of the Kuhn-Munkres algorithm, indicating the index of the assigned constraint hit for each AB pair
     *          (or -1 if unmatched)
     *  @param  nC the number of constraint hits
     *  @param  costMatrix the cost matrix that was used as input to the Kuhn-Munkres algorithm
     *  @param  constraintView the view that was used as the constraint in the chi-squared calculation
     *  @param  chi2Threshold the chi-squared threshold below which a triplet is considered a valid match
     *
     *  @return The list of triplets, where unmatched AB pairs are indicated with a constraint index of -1
     */
    TripletVector BuildTriplets(const PairVector &pairs, const pandora::IntVector &assignment, int nC, const CostMatrix &costMatrix,
        const pandora::HitType constraintView, float chi2Threshold) const;

    /**
     *  @brief  Create a 3D hit from a triplet of U, V and W hits
     *
     *  @param  pCaloHitU the U view hit
     *  @param  pCaloHitV the V view hit
     *  @param  pCaloHitW the W view hit
     *  @param[out] pCaloHit3D the 3D hit to be created
     */
    void CreateThreeDHit(const pandora::CaloHit *const pCaloHitU, const pandora::CaloHit *const pCaloHitV,
        const pandora::CaloHit *const pCaloHitW, const pandora::CaloHit *&pCaloHit3D) const;

    /**
     *  @brief  Select the pairing views based on the specified constraint view
     *
     *  @param  constraintView the view to be used as the constraint in the chi-squared calculation
     *  @param[out] viewA the first view to be used for pairing
     *  @param[out] viewB the second view to be used for pairing
     */
    void SelectViewPair(const pandora::HitType constraintView, pandora::HitType &viewA, pandora::HitType &viewB) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_chi2Threshold;          ///< The chi-squared threshold below which a triplet is considered a valid match
    std::string m_caloHitListName;  ///< The name of the calo hit list containing all of the 2D hits
    std::string m_eventContextName; ///< The name of the EventContext object in which to persist the relationships between the 2D hits
    VolumeToReadoutMap m_planeToReadoutMap; ///< A map from volume, to readout plane, to the 2D hits in that plane
    bool m_visualize;                       ///< Whether to visualize the hit matching process
};

} // namespace lar_content

#endif // #ifndef LAR_PLANE_SOLVER_ALGORITHM_H
