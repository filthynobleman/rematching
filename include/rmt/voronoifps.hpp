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

#include <rmt/graph.hpp>


namespace rmt
{
    
std::vector<double> Distances(const Graph& G, int N);
std::pair<std::vector<int>, std::vector<int>> VoronoiFPS(const rmt::Graph& G, 
                                                         int NSamples, 
                                                         int Seed = 0);

} // namespace rmt
