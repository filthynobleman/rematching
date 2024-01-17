/**
 * @file        flatunion.hpp
 * 
 * @brief       Declaration of the class rmt::FlatUnion, which verifies and enforce
 *              the flat union property.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-15
 */
#pragma once

#include <rmt/mesh.hpp>
#include <rmt/voronoifps.hpp>
#include <rmt/region.hpp>
#include <rmt/utils.hpp>

namespace rmt
{

class FlatUnion
{
private:
    const rmt::Mesh& m_Mesh;
    rmt::VoronoiPartitioning& m_VPart;
    rmt::RegionDictionary m_RDict;

    std::unordered_map<std::tuple<int, int, int>,
                       std::vector<int>,
                       rmt::TripleHash<int>> m_Midpoints;

    std::unordered_map<std::pair<int, int>,
                       std::pair<int, double>,
                       rmt::PairHash<int>> m_BoundBreak;

    std::vector<std::pair<int, double>> m_Farthests;

public:
    FlatUnion(const rmt::Mesh& M,
              rmt::VoronoiPartitioning& VPart);
    ~FlatUnion();

    void DetermineRegions();
    void ComputeTopologies();
    bool FixIssues();
};
    
} // namespace rmt
