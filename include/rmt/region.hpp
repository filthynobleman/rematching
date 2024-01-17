/**
 * @file        region.hpp
 * 
 * @brief       A class representing a region over a triangulated surface.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <rmt/utils.hpp>


namespace rmt
{
    
class SurfaceRegion
{
private:
    // // The number of vertices in the region
    // int m_NVerts;
    // // **TWICE** the number of edges in the region
    // int m_NEdges;
    // // The number of faces in the region
    // int m_NFaces;
    int m_chi;

    // Samples generating the regions
    std::tuple<int, int, int> m_Samples;

public:
    SurfaceRegion(int pi);
    SurfaceRegion(int pi, int pj);
    SurfaceRegion(int pi, int pj, int pk);
    SurfaceRegion(const rmt::SurfaceRegion& SR);
    rmt::SurfaceRegion& operator=(const rmt::SurfaceRegion& SR);
    ~SurfaceRegion();

    std::tuple<int, int, int> GetSamples() const;
    // int NumVertices() const;
    // int NumEdges() const;
    // int NumFaces() const;
    int EulerCharacteristic() const;

    void AddFace();
    // void AddBoundaryFace();
    void AddEdge();
    // void AddBoundaryEdge();
    void AddVertex();

    friend bool operator==(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
    friend bool operator!=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
    friend bool operator<(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
    friend bool operator<=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
    friend bool operator>(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
    friend bool operator>=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
};

bool operator==(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
bool operator!=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
bool operator<(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
bool operator<=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
bool operator>(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);
bool operator>=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR3);


class RegionDictionary
{
private:
    // Vector of regions as unions of three Voronoi regions
    std::vector<rmt::SurfaceRegion> m_Regions;
    std::unordered_set<std::tuple<int, int, int>, rmt::TripleHash<int>> m_RegSamples;

    // Map samples into regions containing them
    std::vector<std::vector<int>> m_VMap;

    // Map couples of samples into regions containing at least one
    std::unordered_map<std::pair<int, int>, 
                       std::vector<int>,
                       rmt::PairHash<int>> m_EMap;
    std::unordered_set<std::pair<int, int>, rmt::PairHash<int>> m_ESet;

    // Map triples of samples into regions containing at least one
    std::unordered_map<std::tuple<int, int, int>,
                       std::vector<int>,
                       rmt::TripleHash<int>> m_TMap;
    std::unordered_set<std::tuple<int, int, int>, rmt::TripleHash<int>> m_TSet;

public:
    RegionDictionary();
    RegionDictionary(size_t NumSamples);
    ~RegionDictionary();

    void Clear();
    void Clear(size_t NumSamples);
    void AddRegion(int pi);
    void AddRegion(int pi, int pj);
    void AddRegion(int pi, int pj, int pk);

    bool HasRegion(int pi) const;
    bool HasRegion(int pi, int pj) const;
    bool HasRegion(int pi, int pj, int pk) const;

    void AddVertex(int pi);
    // void AddBoundaryVertex(int pi);
    void AddEdge(int pi, int pj);
    // void AddBoundaryEdge(int pi, int pj);
    void AddTriangle(int pi, int pj, int pk);

    void BuildRegionMaps();
    size_t NumRegions() const;
    const rmt::SurfaceRegion& GetRegion(size_t i) const;
    bool IsClosed2Ball(size_t i) const;


    static void OrderIndices(int& pi, int& pj);
    static void OrderIndices(int& pi, int& pj, int& pk);
};



} // namespace rmt
