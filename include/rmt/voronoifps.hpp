/**
 * @file        voronoifps.hpp
 * 
 * @brief       Declaration of class rmt::VoronoiPartitioning.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-15
 */
#pragma once

#include <rmt/graph.hpp>
#include <rmt/mesh.hpp>
#include <cut/cut.hpp>


namespace rmt
{
    
class VoronoiPartitioning
{
private:
    rmt::Graph m_G;
    std::vector<int> m_Samples;
    Eigen::VectorXi m_Partitions;
    Eigen::VectorXd m_Distances;
    cut::MinHeap* m_HDists;


public:
    VoronoiPartitioning(const rmt::Mesh& M);
    VoronoiPartitioning(rmt::VoronoiPartitioning&& VP);
    rmt::VoronoiPartitioning& operator=(rmt::VoronoiPartitioning&& VP);
    ~VoronoiPartitioning();

    double GetDistance(int i) const;
    const Eigen::VectorXd& GetDistances() const;
    int GetPartition(int i) const;
    const Eigen::VectorXi& GetPartitions() const;
    
    int NumSamples() const;
    int GetSample(int i) const;
    const std::vector<int>& GetSamples() const;

    int FarthestVertex() const;
    void AddSample(int NewSample);
};



} // namespace rmt
