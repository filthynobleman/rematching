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


rmt::RegionDictionary::RegionDictionary() { }

rmt::RegionDictionary::RegionDictionary(size_t NumSamples)
{
    m_Regions.reserve(3 * NumSamples);
    m_VMap.resize(NumSamples);
    for (auto& RegList : m_VMap)
        RegList.reserve(6);
    m_EMap.reserve(3 * NumSamples);
    m_TMap.reserve(3 * NumSamples);
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
    m_EMap.clear();
    m_TMap.clear();
    m_Regions.clear();


    m_Regions.reserve(3 * NumSamples);
    m_VMap.resize(NumSamples);
    for (auto& RegList : m_VMap)
        RegList.reserve(6);
    m_EMap.reserve(3 * NumSamples);
    m_TMap.reserve(3 * NumSamples);
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
    auto Reg = &m_Regions[m_Regions.size() - 1];

    m_VMap[pi].emplace_back(Reg);
    m_VMap[pj].emplace_back(Reg);
    m_VMap[pk].emplace_back(Reg);

    m_EMap[{ pi, pj }].emplace_back(Reg);
    m_EMap[{ pj, pk }].emplace_back(Reg);
    m_EMap[{ pi, pk }].emplace_back(Reg);

    std::tuple<int, int, int> T{pi, pj, pk};
    m_TMap.emplace(T, Reg);
}

bool rmt::RegionDictionary::HasRegion(int pi, int pj, int pk) const
{
    OrderIndices(pi, pj, pk);
    return m_TMap.find({ pi, pj, pk }) != m_TMap.end();
}



void rmt::RegionDictionary::AddVertex(int pi)
{
    // Adding a primal vertex means adding a face
    for (auto& Reg : m_VMap[pi])
        Reg->AddFace();
}

void rmt::RegionDictionary::AddBoundaryVertex(int pi)
{
    // Adding a primal vertex means adding a face
    for (auto& Reg : m_VMap[pi])
        Reg->AddBoundaryFace();
}

void rmt::RegionDictionary::AddEdge(int pi, int pj)
{
    // Adding a primal edge means adding a dual edge
    // If they are the same sample, we add an edge to all the regions that contain that sample
    if (pi == pj)
    {
        for (auto& Reg : m_VMap[pi])
            Reg->AddEdge();
        return;
    }

    // Otherwise add an edge to the regions that contain at leas one of the samples
    std::vector<rmt::SurfaceRegion*> Pij = m_VMap[pi];
    Pij.insert(Pij.end(), m_VMap[pj].begin(), m_VMap[pj].end());
    std::sort(Pij.begin(), Pij.end());
    auto PijEnd = std::unique(Pij.begin(), Pij.end());
    for (auto it = Pij.begin(); it != PijEnd; ++it)
        (*it)->AddEdge();
}

void rmt::RegionDictionary::AddBoundaryEdge(int pi, int pj)
{
    // Adding a primal edge means adding a dual edge
    // If they are the same sample, we add an edge to all the regions that contain that sample
    if (pi == pj)
    {
        for (auto& Reg : m_VMap[pi])
            Reg->AddBoundaryEdge();
        return;
    }

    // Otherwise add an edge to the regions that contain at leas one of the samples
    std::vector<rmt::SurfaceRegion*> Pij = m_VMap[pi];
    Pij.insert(Pij.end(), m_VMap[pj].begin(), m_VMap[pj].end());
    std::sort(Pij.begin(), Pij.end());
    auto PijEnd = std::unique(Pij.begin(), Pij.end());
    for (auto it = Pij.begin(); it != PijEnd; ++it)
        (*it)->AddBoundaryEdge();
}

void rmt::RegionDictionary::AddTriangle(int pi, int pj, int pk)
{
    // Adding a primal triangle means adding a dual vertex
    // If all vertices belongs to the same partition, we add a vertex to all
    // the regions containing that partition
    if (pi == pj && pj == pk)
    {
        for (auto& Reg : m_VMap[pi])
            Reg->AddVertex();
        return;
    }

    // Ease the identification of different partitions
    OrderIndices(pi, pj, pk);
    // If two vertices belongs to the same partition, we add a vertex to all
    // the regions containing at least one between that partition and the other
    if (pi == pj || pj == pk)
    {
        std::vector<rmt::SurfaceRegion*> Pik = m_VMap[pi];
        Pik.insert(Pik.end(), m_VMap[pk].begin(), m_VMap[pk].end());
        std::sort(Pik.begin(), Pik.end());
        auto PikEnd = std::unique(Pik.begin(), Pik.end());
        for (auto it = Pik.begin(); it != PikEnd; ++it)
            (*it)->AddVertex();
        return;
    }

    // If all vertices belongs to different partitions, we add a vertex to all
    // the regions containing at least one of those partitions
    std::vector<rmt::SurfaceRegion*> Pijk = m_VMap[pi];
    Pijk.insert(Pijk.end(), m_VMap[pj].begin(), m_VMap[pj].end());
    Pijk.insert(Pijk.end(), m_VMap[pk].begin(), m_VMap[pk].end());
    std::sort(Pijk.begin(), Pijk.end());
    auto PijkEnd = std::unique(Pijk.begin(), Pijk.end());
    for (auto it = Pijk.begin(); it != PijkEnd; ++it)
        (*it)->AddVertex();
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