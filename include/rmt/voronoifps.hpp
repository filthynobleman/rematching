/**
 * @file        voronoifps.hpp
 * 
 * @brief       Computation of fast geodesic Voronoi FPS.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-17
 */
#pragma once

#include <cut/cut.hpp>
#include <rmt/graph.hpp>


namespace rmt
{

struct VoronoiPartitioning
{
    std::vector<int> Samples;
    std::vector<int> Partition;
    std::vector<double> Distances;
};
    
std::vector<double> Distances(const Graph& G, int N);
rmt::VoronoiPartitioning VoronoiFPS(const rmt::Graph& G, 
                                     int NSamples, 
                                     int Seed = 0);
void UpdateVoronoi(const rmt::Graph& G,
                   std::vector<double>& Dists,
                   std::vector<int>& Partition,
                   cut::MinHeap& HDists,
                   int p);

} // namespace rmt
