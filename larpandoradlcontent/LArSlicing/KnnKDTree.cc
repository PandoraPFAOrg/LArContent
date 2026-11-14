/**
 *  @file   larpandoradlcontent/LArSlicing/KnnKDTree.cc
 *
 *  @brief  Implementation of k-NN KD-Tree.
 *
 *  $Log: $
 */

#include <algorithm>

#include "larpandoradlcontent/LArSlicing/KnnKDTree.h"

namespace lar_content
{

KnnKdTree::KnnKdTree(const std::vector<KnnNode> &inputNodes) :
    m_nodes(std::move(inputNodes))
{
    this->BuildTree(0, m_nodes.size(), 0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void KnnKdTree::BuildTree(int start, int end, int depth)
{
    if (start >= end)
        return;

    const int axis = depth % 3;
    const int mid = start + (end - start) / 2;

    std::nth_element(m_nodes.begin() + start, m_nodes.begin() + mid, m_nodes.begin() + end,
        [axis](const KnnNode &a, const KnnNode &b) { return a.coords[axis] < b.coords[axis]; });

    this->BuildTree(start, mid, depth + 1);
    this->BuildTree(mid + 1, end, depth + 1);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> KnnKdTree::FindNearestNeighbours(const KnnNode &query, const int k) const
{
    std::priority_queue<std::pair<float, int>> pq;
    this->SearchTree(0, m_nodes.size(), 0, query, k, pq);

    std::vector<int> neighbors;
    neighbors.reserve(pq.size());
    while (!pq.empty())
    {
        // Extract the original_id from the stored index
        neighbors.push_back(m_nodes[pq.top().second].original_id);
        pq.pop();
    }

    return neighbors;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> KnnKdTree::FindWithinRadius(const KnnNode &query, const float radius) const
{
    std::vector<int> neighbors;
    this->SearchRadius(0, m_nodes.size(), 0, query, radius * radius, neighbors);
    return neighbors;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> KnnKdTree::FindWithinRadiusSqd(const KnnNode &query, const float radiusSqd) const
{
    std::vector<int> neighbors;
    this->SearchRadius(0, m_nodes.size(), 0, query, radiusSqd, neighbors);
    return neighbors;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void KnnKdTree::SearchTree(int start, int end, int depth, const KnnNode &query, const int k, std::priority_queue<std::pair<float, int>> &pq) const
{
    if (start >= end)
        return;

    const int axis = depth % 3;
    const int mid = start + (end - start) / 2;
    const KnnNode &currentNode = m_nodes[mid];

    float dist_sq = 0.f;
    for (int i = 0; i < 3; ++i)
    {
        const float diff = query.coords[i] - currentNode.coords[i];
        dist_sq += diff * diff;
    }

    // Don't self-match
    if (query.original_id != currentNode.original_id)
    {
        if (pq.size() < static_cast<size_t>(k))
        {
            pq.emplace(dist_sq, mid);
        }
        else if (dist_sq < pq.top().first)
        {
            pq.pop();
            pq.emplace(dist_sq, mid);
        }
    }

    const float axis_diff = query.coords[axis] - currentNode.coords[axis];
    int first_half_start = start, first_half_end = mid;
    int second_half_start = mid + 1, second_half_end = end;

    if (axis_diff > 0.f)
    {
        std::swap(first_half_start, second_half_start);
        std::swap(first_half_end, second_half_end);
    }

    this->SearchTree(first_half_start, first_half_end, depth + 1, query, k, pq);

    if (pq.size() < static_cast<size_t>(k) || (axis_diff * axis_diff) < pq.top().first)
        this->SearchTree(second_half_start, second_half_end, depth + 1, query, k, pq);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void KnnKdTree::SearchRadius(int start, int end, int depth, const KnnNode &query, const float radiusSq, std::vector<int> &neighbors) const
{
    if (start >= end)
        return;

    const int axis = depth % 3;
    const int mid = start + (end - start) / 2;
    const KnnNode &currentNode = m_nodes[mid];

    float dist_sq = 0.f;
    for (int i = 0; i < 3; ++i)
    {
        const float diff = query.coords[i] - currentNode.coords[i];
        dist_sq += diff * diff;
    }

    if (dist_sq <= radiusSq)
        neighbors.push_back(currentNode.original_id);

    const float axis_diff = query.coords[axis] - currentNode.coords[axis];

    // Query right or left first, depending on which side of the splitting plane it's on
    if (axis_diff > 0.f)
    {
        this->SearchRadius(mid + 1, end, depth + 1, query, radiusSq, neighbors);
        if ((axis_diff * axis_diff) <= radiusSq)
            this->SearchRadius(start, mid, depth + 1, query, radiusSq, neighbors);
    }
    else
    {
        this->SearchRadius(start, mid, depth + 1, query, radiusSq, neighbors);
        if ((axis_diff * axis_diff) <= radiusSq)
            this->SearchRadius(mid + 1, end, depth + 1, query, radiusSq, neighbors);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
