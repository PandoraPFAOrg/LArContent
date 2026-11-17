/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_CREATION_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerToolsT.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  SimpleClusterCreationAlgorithm class
 */
class SimpleClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleClusterCreationAlgorithm();

protected:
    typedef std::unordered_map<const pandora::CaloHit *, pandora::CaloHitVector> HitAssociationMap;

    /**
     *  @brief  Find the hits within the clustering window of each hit
     *
     *  @param  caloHitVector the calo hits
     *  @param  hitAssociationMap to receive, for each hit with any, its nearby hits sorted by position
     */
    virtual void BuildAssociationMap(const pandora::CaloHitVector &caloHitVector, HitAssociationMap &hitAssociationMap) const;

    /**
     *  @brief  Implementation of BuildAssociationMap, using a kd tree to find the hits within the clustering window of each hit
     *
     *  @tparam DIM the number of coordinates to search in: 2 for (x, z), or 3 for (x, y, z)
     *
     *  @param  caloHitVector the calo hits
     *  @param  hitAssociationMap to receive, for each hit with any, its nearby hits sorted by position
     */
    template <unsigned int DIM>
    void BuildAssociationMapKDTree(const pandora::CaloHitVector &caloHitVector, HitAssociationMap &hitAssociationMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_clusteringWindowSquared;     ///< Maximum distance (squared) for two hits to be joined
    std::string m_inputCaloHitListName;  ///< Name of the input calo hit list (empty: use the current list)
    std::string m_outputClusterListName; ///< Name to save the new clusters under (empty: leave in the current cluster list)

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Create a cluster from each group of hits joined together by the associations
     *
     *  @param  caloHitVector the calo hits, sorted by position
     *  @param  hitAssociationMap the associations between hits
     */
    void CreateClusters(const pandora::CaloHitVector &caloHitVector, const HitAssociationMap &hitAssociationMap) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <unsigned int DIM>
void SimpleClusterCreationAlgorithm::BuildAssociationMapKDTree(const pandora::CaloHitVector &caloHitVector, HitAssociationMap &hitAssociationMap) const
{
    static_assert(DIM == 2 || DIM == 3, "Only 2D (x, z) and 3D (x, y, z) searches are supported");

    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, DIM> HitKDTree;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, DIM> HitKDNode;

    const pandora::CaloHitList caloHitList(caloHitVector.begin(), caloHitVector.end());
    std::vector<HitKDNode> hitKDNodes;
    HitKDTree kdTree;

    if constexpr (DIM == 2)
        kdTree.build(hitKDNodes, fill_and_bound_2d_kd_tree(caloHitList, hitKDNodes));
    else
        kdTree.build(hitKDNodes, fill_and_bound_3d_kd_tree(caloHitList, hitKDNodes));

    const float window(std::sqrt(m_clusteringWindowSquared));

    for (const pandora::CaloHit *const pCaloHitI : caloHitVector)
    {
        std::vector<HitKDNode> found;

        if constexpr (DIM == 2)
            kdTree.search(build_2d_kd_search_region(pCaloHitI, window, window), found);
        else
            kdTree.search(build_3d_kd_search_region(pCaloHitI, window, window, window), found);

        // The kd tree search region is a box, so apply the exact distance cut
        pandora::CaloHitVector nearbyHits;

        for (const HitKDNode &node : found)
        {
            const pandora::CaloHit *const pCaloHitJ(node.data);

            if (pCaloHitJ != pCaloHitI && (pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_clusteringWindowSquared)
                nearbyHits.push_back(pCaloHitJ);
        }

        if (!nearbyHits.empty())
        {
            std::sort(nearbyHits.begin(), nearbyHits.end(), LArClusterHelper::SortHitsByPosition);
            hitAssociationMap.emplace(pCaloHitI, std::move(nearbyHits));
        }
    }
}

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_CREATION_ALGORITHM_H
