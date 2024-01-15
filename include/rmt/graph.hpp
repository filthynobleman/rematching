/**
 * @file        graph.hpp
 * 
 * @brief       Declaration of a graph like data structure.
 * 
 * @details     This file contains the declaration of a class representing a graph
 *              embedded in 3D space.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-17
 */
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <set>

namespace rmt
{

/**
 * @brief       Weighted edge of a graph.
 * 
 * @details     A weighted edge of a graph is defined as a destintation and a cost.
 */
typedef std::pair<int, double> WEdge;

/**
 * @brief       Path inside a weighted graph.
 * 
 * @details     A path inside a weighted graph is defined as sequence of weighted edges
 *              and the sum of the total weights. The starting node has cost zero.
 */
typedef std::pair<double, std::vector<rmt::WEdge>> Path;

/**
 * @brief       A graph-like data structure.
 * 
 * @details     This class represents a graph embedded in 3D space.\n 
 *              The embedding of the graph determines the weights of the edges, since
 *              the weight of each edge is defined as its Euclidean length.
 */
class Graph
{
private:
    // std::vector<Eigen::Vector3d> m_Verts;
    std::vector<int> m_Idxs;
    std::vector<WEdge> m_Adjs;

public:
    Graph(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
    Graph(const Eigen::MatrixXd& V, const std::vector<std::pair<int, int>>& E);
    Graph(const Eigen::MatrixXd& V, const std::set<std::pair<int, int>>& E);
    Graph(const rmt::Graph& G);
    Graph& operator=(const rmt::Graph& G);
    Graph(rmt::Graph&& G);
    Graph& operator=(rmt::Graph&& G);
    ~Graph();

    int NumVertices() const;
    int NumAdjacents(int i) const;
    int NumEdges() const;

    // const Eigen::Vector3d& GetVertex(int i) const;
    // const std::vector<Eigen::Vector3d>& GetVertices() const;

    const WEdge& GetAdjacent(int node_i, int adj_i) const;

    rmt::Path DijkstraPath(int src, int dst) const;
    Eigen::VectorXd DijkstraDistance(int src) const;
    void DijkstraDistance(int src, Eigen::VectorXd& Distances) const;
    int FarthestFiltered(int src, const std::vector<int>& Tag, int Filter) const;
    int FarthestAtBoundary(int src, const std::vector<int>& Tag, int Region, int Neighbor) const;
    std::vector<int> ConnectedComponents() const;
};

} // namespace rmt
