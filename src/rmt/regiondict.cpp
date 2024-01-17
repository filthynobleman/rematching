/**
 * @file        regiondict.cpp
 * 
 * @brief       Implementation of class rmt::RegionDictionary.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#include <rmt/region.hpp>
#include <cut/cut.hpp>
#include <iostream>

rmt::RegionDictionary::RegionDictionary() { }

rmt::RegionDictionary::RegionDictionary(size_t NumSamples)
{
    Clear(NumSamples);
}


rmt::RegionDictionary::~RegionDictionary() { }

void rmt::RegionDictionary::Clear()
{
    size_t NumSamples = m_VMap.size();
    Clear(NumSamples);
}

void rmt::RegionDictionary::Clear(size_t NumSamples)
{
    m_VMap.clear();
    m_Regions.clear();
    m_RegSamples.clear();

    m_ESet.clear();
    m_ESet.reserve(3 * NumSamples);
    m_TSet.clear();
    m_TSet.reserve(3 * NumSamples);

    m_Regions.reserve(3 * NumSamples);
    m_RegSamples.reserve(3 * NumSamples);
    m_VMap.resize(NumSamples);
    for (auto& RegList : m_VMap)
        RegList.reserve(6);
}


void rmt::RegionDictionary::OrderIndices(int& pi, int& pj)
{
    if (pi > pj)
        std::swap(pi, pj);
}

void rmt::RegionDictionary::OrderIndices(int& pi, int& pj, int& pk)
{
    if (pi > pj)
        std::swap(pi, pj);
    // now pi <= pj
    if (pj > pk)
    {
        std::swap(pj, pk);
        // now pi <= pk and pj <= pk
        if (pi > pj)
            std::swap(pi, pj);
        // now pi <= pj <= pk
    }
}


void rmt::RegionDictionary::AddRegion(int pi, int pj, int pk)
{
    if (pi == pj || pj == pk || pk == pi)
        return;

    OrderIndices(pi, pj, pk);
    if (HasRegion(pi, pj, pk))
        return;


    m_Regions.emplace_back(pi, pj, pk);
    auto& Reg = m_Regions[(int)m_Regions.size() - 1];
    m_RegSamples.emplace(Reg.GetSamples());

    m_VMap[pi].emplace_back((int)m_Regions.size() - 1);
    m_VMap[pj].emplace_back((int)m_Regions.size() - 1);
    m_VMap[pk].emplace_back((int)m_Regions.size() - 1);

    m_ESet.emplace(pi, pj);
    m_ESet.emplace(pj, pk);
    m_ESet.emplace(pi, pk);

    m_TSet.emplace(pi, pj, pk);
}

void rmt::RegionDictionary::AddRegion(int pi, int pj)
{
    if (pi == pj)
        return;

    OrderIndices(pi, pj);
    if (HasRegion(pi, pj))
        return;

    m_Regions.emplace_back(pi, pj);
    auto& Reg = m_Regions[(int)m_Regions.size() - 1];
    m_RegSamples.emplace(Reg.GetSamples());

    m_VMap[pi].emplace_back((int)m_Regions.size() - 1);
    m_VMap[pj].emplace_back((int)m_Regions.size() - 1);

    m_ESet.emplace(pi, pj);
}

void rmt::RegionDictionary::AddRegion(int pi)
{
    if (HasRegion(pi))
        return;

    m_Regions.emplace_back(pi);
    auto& Reg = m_Regions[(int)m_Regions.size() - 1];
    m_RegSamples.emplace(Reg.GetSamples());
    
    m_VMap[pi].emplace_back((int)m_Regions.size() - 1);
    
    // m_TMap.emplace(Reg.GetSamples(), (int)m_Regions.size() - 1);
}

bool rmt::RegionDictionary::HasRegion(int pi, int pj, int pk) const
{
    OrderIndices(pi, pj, pk);
    return m_RegSamples.find({ pi, pj, pk }) != m_RegSamples.end();
}

bool rmt::RegionDictionary::HasRegion(int pi, int pj) const
{
    return HasRegion(pi, pj, std::numeric_limits<int>::max());
}

bool rmt::RegionDictionary::HasRegion(int pi) const
{
    return HasRegion(pi, std::numeric_limits<int>::max() - 1);
}



void rmt::RegionDictionary::AddVertex(int pi)
{
    // Adding a primal vertex means adding a face
    for (int RID : m_VMap[pi])
        m_Regions[RID].AddFace();
}

// void rmt::RegionDictionary::AddBoundaryVertex(int pi)
// {
//     // Adding a primal vertex means adding a face
//     for (int RID : m_VMap[pi])
//         m_Regions[RID].AddBoundaryFace();
// }

void rmt::RegionDictionary::AddEdge(int pi, int pj)
{
    // Adding a primal edge means adding a dual edge
    // If they are the same sample, we add an edge to all the regions that contain that sample
    if (pi == pj)
    {
        for (int RID : m_VMap[pi])
            m_Regions[RID].AddEdge();
        return;
    }

    // Otherwise add an edge to the regions that contain at leas one of the samples
    OrderIndices(pi, pj);
    std::pair<int, int> e{pi, pj};
    for (int i : m_EMap[e])
        m_Regions[i].AddEdge();
}

// void rmt::RegionDictionary::AddBoundaryEdge(int pi, int pj)
// {
//     // Adding a primal edge means adding a dual edge
//     // If they are the same sample, we add an edge to all the regions that contain that sample
//     if (pi == pj)
//     {
//         for (int RID : m_VMap[pi])
//             m_Regions[RID].AddBoundaryEdge();
//         return;
//     }

//     // Otherwise add an edge to the regions that contain at leas one of the samples
//     OrderIndices(pi, pj);
//     std::pair<int, int> e{pi, pj};
//     for (int i : m_EMap[e])
//         m_Regions[i].AddBoundaryEdge();
// }

void rmt::RegionDictionary::AddTriangle(int pi, int pj, int pk)
{
    // Adding a primal triangle means adding a dual vertex
    // If all vertices belongs to the same partition, we add a vertex to all
    // the regions containing that partition
    if (pi == pj && pj == pk)
    {
        for (int RID : m_VMap[pi])
            m_Regions[RID].AddVertex();
        return;
    }

    // Ease the identification of different partitions
    OrderIndices(pi, pj, pk);
    // If two vertices belongs to the same partition, we add a vertex to all
    // the regions containing at least one between that partition and the other
    if (pi == pj || pj == pk)
    {
        std::pair<int, int> e{pi, pk};
        for (int i : m_EMap[e])
            m_Regions[i].AddVertex();
        return;
    }

    // If all vertices belongs to different partitions, we add a vertex to all
    // the regions containing at least one of those partitions
    std::tuple<int, int, int> t{pi, pj, pk};
    for (int i : m_TMap[t])
        m_Regions[i].AddVertex();
}



size_t rmt::RegionDictionary::NumRegions() const { return m_Regions.size(); }

const rmt::SurfaceRegion& rmt::RegionDictionary::GetRegion(size_t i) const
{
    CUTCheckLess(i, NumRegions());
    return m_Regions[i];
}

bool rmt::RegionDictionary::IsClosed2Ball(size_t i) const
{
    CUTCheckLess(i, NumRegions());
    return m_Regions[i].EulerCharacteristic() == 1;
}


void rmt::RegionDictionary::BuildRegionMaps()
{
    m_EMap.clear();
    m_EMap.reserve(m_ESet.size());
    for (auto e : m_ESet)
    {
        m_EMap.emplace(e, std::vector<int>{});
        std::vector<int>& Pe = m_EMap[e];
        Pe.insert(Pe.end(), m_VMap[e.first].begin(), m_VMap[e.first].end());
        Pe.insert(Pe.end(), m_VMap[e.second].begin(), m_VMap[e.second].end());
        std::sort(Pe.begin(), Pe.end());
        auto PeEnd = std::unique(Pe.begin(), Pe.end());
        Pe.erase(PeEnd, Pe.end());
    }
    
    m_TMap.clear();
    m_TMap.reserve(m_TSet.size());
    for (auto t : m_TSet)
    {
        m_TMap.emplace(t, std::vector<int>{});
        std::vector<int>& Pt = m_TMap[t];
        Pt.insert(Pt.end(), m_VMap[std::get<0>(t)].begin(), m_VMap[std::get<0>(t)].end());
        Pt.insert(Pt.end(), m_VMap[std::get<1>(t)].begin(), m_VMap[std::get<1>(t)].end());
        Pt.insert(Pt.end(), m_VMap[std::get<2>(t)].begin(), m_VMap[std::get<2>(t)].end());
        std::sort(Pt.begin(), Pt.end());
        auto PtEnd = std::unique(Pt.begin(), Pt.end());
        Pt.erase(PtEnd, Pt.end());
    }
}