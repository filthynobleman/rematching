/**
 * @file        flatunion.cpp
 * 
 * @brief       Implementation of rmt::FlatUnion.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-15
 */
#include <rmt/flatunion.hpp>


rmt::FlatUnion::FlatUnion(const rmt::Mesh& M, rmt::VoronoiPartitioning& VPart)
    : m_Mesh(M), m_VPart(VPart) { }

rmt::FlatUnion::~FlatUnion() { }


void rmt::FlatUnion::DetermineRegions()
{
    int NSamples = m_VPart.NumSamples();
    m_Midpoints.clear();
    m_Midpoints.reserve(3 * NSamples);
    m_BoundBreak.clear();
    m_BoundBreak.reserve(3 * NSamples);
    m_RDict.Clear(NSamples);

    // Add a region for each sample point
    for (int i = 0; i < NSamples; ++i)
        m_RDict.AddRegion(i);

    m_Farthests.clear();
    m_Farthests.resize(NSamples);
    for (int i = 0; i < NSamples; ++i)
        m_Farthests[i] = { m_VPart.GetSample(i), 0.0 };

    int NFaces = m_Mesh.NumTriangles();
    const Eigen::MatrixXi& F = m_Mesh.GetTriangles();
    const Eigen::VectorXi& P = m_VPart.GetPartitions();
    const Eigen::VectorXd& D = m_VPart.GetDistances();
    for (int i = 0; i < NFaces; ++i)
    {
        // Update farthests inside regions
        for (int j = 0; j < 3; ++j)
        {
            int v = F(i, j);
            int p = P[v];
            if (D[v] > m_Farthests[p].second)
                m_Farthests[p] = { v, D[v] };
        }

        int p0 = P[F(i, 0)];
        int p1 = P[F(i, 1)];
        int p2 = P[F(i, 2)];

        // If triangle crosses two regions, add the union region 
        // and update the eventual breakpoint
        if (p0 != p1)
        {
            std::pair<int, int> p01{p0, p1};
            if (p1 < p0)
                std::swap(p01.first, p01.second);
            std::pair<int, double> v01{F(i, 0), D[F(i, 0)]};
            if (D[F(i, 1)] > v01.second)
                v01 = { F(i, 1), D[F(i, 1)] };
            if (m_BoundBreak.find(p01) == m_BoundBreak.end())
            {
                m_RDict.AddRegion(p0, p1);
                m_BoundBreak.emplace(p01, v01);
            }
            if (v01.second > m_BoundBreak[p01].second)
                m_BoundBreak[p01] = v01;
        }
        if (p1 != p2)
        {
            std::pair<int, int> p12{p1, p2};
            if (p2 < p1)
                std::swap(p12.first, p12.second);
            std::pair<int, double> v12{F(i, 1), D[F(i, 1)]};
            if (D[F(i, 2)] > v12.second)
                v12 = { F(i, 2), D[F(i, 2)] };
            if (m_BoundBreak.find(p12) == m_BoundBreak.end())
            {
                m_RDict.AddRegion(p1, p2);
                m_BoundBreak.emplace(p12, v12);
            }
            if (v12.second > m_BoundBreak[p12].second)
                m_BoundBreak[p12] = v12;
        }
        if (p2 != p0)
        {
            std::pair<int, int> p20{p2, p0};
            if (p0 < p2)
                std::swap(p20.first, p20.second);
            std::pair<int, double> v20{F(i, 2), D[F(i, 2)]};
            if (D[F(i, 0)] > v20.second)
                v20 = { F(i, 0), D[F(i, 0)] };
            if (m_BoundBreak.find(p20) == m_BoundBreak.end())
            {
                m_RDict.AddRegion(p2, p0);
                m_BoundBreak.emplace(p20, v20);
            }
            if (v20.second > m_BoundBreak[p20].second)
                m_BoundBreak[p20] = v20;
        }

        // If triangle not crossing three regions, skip
        if (p0 == p1)
            continue;
        if (p1 == p2)
            continue;
        if (p2 == p0)
            continue;

        // Add the midpoint
        rmt::RegionDictionary::OrderIndices(p0, p1, p2);
        std::tuple<int, int, int> T(p0, p1, p2);
        if (m_Midpoints.find(T) == m_Midpoints.end())
            m_Midpoints.emplace(T, std::vector<int>{});
        // Order the indices by distance to the sample set
        int v0 = F(i, 0);
        int v1 = F(i, 1);
        int v2 = F(i, 2);
        if (D[v0] < D[v1])
            std::swap(v0, v1);
        if (D[v1] < D[v2])
        {
            std::swap(v1, v2);
            if (D[v0] < D[v1])
                std::swap(v0, v1);
        }
        // Additional check: if all the vertices are samples, this is going to be
        // a triangle in the final mesh, and we are not going to add the midpoint
        // to the list of new samples candidates
        if (v0 != m_VPart.GetSample(P[v0]))
            m_Midpoints[T].emplace_back(v0);
        else if (v1 != m_VPart.GetSample(P[v1]))
            m_Midpoints[T].emplace_back(v1);
        else if (v2 != m_VPart.GetSample(P[v2]))
            m_Midpoints[T].emplace_back(v2);

        // Add the region
        m_RDict.AddRegion(p0, p1, p2);
    }

    m_RDict.BuildRegionMaps();
}




void rmt::FlatUnion::ComputeTopologies()
{
    const Eigen::VectorXi& P = m_VPart.GetPartitions();


    int NVerts = m_Mesh.NumVertices();
    const Eigen::VectorXi& BV = m_Mesh.GetBoundaryVertices();
    for (int i = 0; i < NVerts; ++i)
    {
        CUTAssert(P[i] < m_VPart.NumSamples());
        CUTAssert(P[i] >= 0);
        if (!BV[i])
            m_RDict.AddVertex(P[i]);
        // Boundary vertices do not contribute to topological changes
        // else
        //     m_RDict.AddBoundaryVertex(P[i]);
    }

    const Eigen::MatrixXi& E = m_Mesh.GetEdges();
    const Eigen::VectorXi& BE = m_Mesh.GetBoundaryEdges();
    int NEdges = E.rows();
    for (int i = 0; i < NEdges; ++i)
    {
        if (!BE[i])
            m_RDict.AddEdge(P[E(i, 0)], P[E(i, 1)]);
        // Boundary edges do not contribute to topological changes
        // else
        //     m_RDict.AddBoundaryEdge(P[E(i, 0)], P[E(i, 1)]);
    }

    const Eigen::MatrixXi& F = m_Mesh.GetTriangles();
    int NTris = F.rows();
    for (int i = 0; i < NTris; ++i)
        m_RDict.AddTriangle(P[F(i, 0)], P[F(i, 1)], P[F(i, 2)]);
}


bool rmt::FlatUnion::FixIssues()
{
    std::set<int> NewSamples;

    // For each region that is not a closed 2-ball, fix it
    size_t NRegions = m_RDict.NumRegions();
    for (size_t i = 0; i < NRegions; ++i)
    {
        if (m_RDict.IsClosed2Ball(i))
            continue;

        // Get the samples
        std::tuple<int, int, int> T = m_RDict.GetRegion(i).GetSamples();
        
        // If is a Voronoi texel, add a sample from its boundary
        if (std::get<1>(T) >= m_VPart.NumSamples())
        {
            int p = std::get<0>(T);
            NewSamples.insert(m_Farthests[p].first);
            continue;
        }
        // If is a union of two texels, add a point at the boundary
        if (std::get<2>(T) >= m_VPart.NumSamples())
        {
            int p0 = std::get<0>(T);
            int p1 = std::get<1>(T);
            int v = m_BoundBreak[{ p0, p1 }].first;
            // A sample cannot insert itself
            if (m_VPart.GetSample(p0) == v)
                continue;
            if (m_VPart.GetSample(p1) == v)
                continue;
            NewSamples.insert(v);
            continue;
        }
        // If is a union of three texels, add all the midpoints
        NewSamples.insert(m_Midpoints[T].begin(), m_Midpoints[T].end());
    }

    // Add the samples
    for (int v : NewSamples)
        m_VPart.AddSample(v);

    return NewSamples.size() == 0;
}