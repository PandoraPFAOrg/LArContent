/**
 *  @file   larpandoradlcontent/LArSlicing/KnnKDTree.h
 *
 *  @brief  Header file for a KD-Tree
 *
 *  $Log: $
 */
#ifndef LAR_KNN_KD_TREE_H
#define LAR_KNN_KD_TREE_H 1

#include <queue>
#include <utility>
#include <vector>

namespace lar_content
{

/**
 * @brief  KnnKdTree class
 */
class KnnKdTree
{
public:
    /**
     * @brief Simple struct to hold node information for the KD-Tree.
     */
    struct KnnNode
    {
        float coords[3];
        int original_id;
    };

    /**
     * @brief  Constructor builds the KD-Tree immediately from the provided nodes
     *
     * @param  inputNodes The flat list of nodes to build the tree from
     */
    KnnKdTree(const std::vector<KnnNode> &inputNodes);

    /**
     * @brief  Search for the k nearest neighbors of a query point
     *
     * @param  query The node to search around
     * @param  k The number of neighbors to find
     *
     * @return A vector of the original_ids of the nearest neighbors
     */
    std::vector<int> FindNearestNeighbours(const KnnNode &query, const int k) const;

    /**
     * @brief  Search for all neighbors within a certain radius of a query point.
     *
     * @param  query The node to search around
     * @param  radius The radius within which to search for neighbors
     *
     * @return A vector of the original_ids of the neighbors within the radius
     */
    std::vector<int> FindWithinRadius(const KnnNode &query, const float radius) const;

    /**
     * @brief  Search for all neighbors within a certain radius of a query point.
     *
     * @param  query The node to search around
     * @param  radiusSqd The square of the radius within which to search
     *
     * @return A vector of the original_ids of the neighbors within the radius
     */
    std::vector<int> FindWithinRadiusSqd(const KnnNode &query, const float radiusSqd) const;

private:
    /**
     * @brief Recursive tree builder
     *
     * @param  start The starting index of the current subtree in m_nodes
     * @param  end The ending index of the current subtree in m_nodes
     * @param  depth The current depth in the tree, used to determine splitting dimension
     */
    void BuildTree(int start, int end, int depth);

    /**
     * @brief Recursive tree searcher
     *
     * @param  start The starting index of the current subtree in m_nodes
     * @param  end The ending index of the current subtree in m_nodes
     * @param  depth The current depth in the tree, used to determine splitting dimension
     * @param  query The node to search around
     * @param  k The number of neighbors to find
     * @param  pq The priority queue to store the nearest neighbors
     */
    void SearchTree(int start, int end, int depth, const KnnNode &query, const int k, std::priority_queue<std::pair<float, int>> &pq) const;

    /**
     * @brief Recursive radius searcher
     *
     * @param  start The starting index of the current subtree in m_nodes
     * @param  end The ending index of the current subtree in m_nodes
     * @param  depth The current depth in the tree, used to determine splitting dimension
     * @param  query The node to search around
     * @param  radiusSq The square of the radius within which to search for neighbors
     * @param  neighbors The vector to store the original_ids of the neighbors within the radius
     */
    void SearchRadius(int start, int end, int depth, const KnnNode &query, const float radiusSq, std::vector<int> &neighbors) const;

    std::vector<KnnNode> m_nodes;
};

} // namespace lar_content

#endif // LAR_KNN_KD_TREE_H
