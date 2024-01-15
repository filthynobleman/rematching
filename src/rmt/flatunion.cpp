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
    m_RDict.Clear(NSamples);

    m_Neighs.resize(NSamples);
    m_Farthests.clear();
    m_Farthests.resize(NSamples);
    m_FarDists.clear();
    m_FarDists.resize(NSamples);
    for (int i = 0; i < NSamples; ++i)
    {
        m_Farthests[i] = { m_VPart.GetSample(i), m_VPart.GetSample(i), m_VPart.GetSample(i) };
        m_FarDists[i] = { 0.0, 0.0, 0.0 };
    }

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
            if (D[v] > std::get<0>(m_FarDists[p]))
            {
                std::get<2>(m_Farthests[p]) = std::get<1>(m_Farthests[p]);
                std::get<2>(m_FarDists[p]) = std::get<1>(m_FarDists[p]);
                std::get<1>(m_Farthests[p]) = std::get<0>(m_Farthests[p]);
                std::get<1>(m_FarDists[p]) = std::get<0>(m_FarDists[p]);
                std::get<0>(m_FarDists[p]) = D[v];
                std::get<0>(m_Farthests[p]) = v;
                continue;
            }
            if (D[v] > std::get<1>(m_FarDists[p]))
            {
                std::get<2>(m_Farthests[p]) = std::get<1>(m_Farthests[p]);
                std::get<2>(m_FarDists[p]) = std::get<1>(m_FarDists[p]);
                std::get<1>(m_FarDists[p]) = D[v];
                std::get<1>(m_Farthests[p]) = v;
                continue;
            }
            if (D[v] > std::get<2>(m_FarDists[p]))
            {
                std::get<2>(m_FarDists[p]) = D[v];
                std::get<2>(m_Farthests[p]) = v;
                continue;
            }
        }

        int p0 = P[F(i, 0)];
        int p1 = P[F(i, 1)];
        int p2 = P[F(i, 2)];

        // If triangle crosses two regions, add the edge
        if (p0 != p1)
        {
            m_Neighs[p0].insert(p1);
            m_Neighs[p1].insert(p0);
        }
        if (p1 != p2)
        {
            m_Neighs[p1].insert(p2);
            m_Neighs[p2].insert(p1);
        }
        if (p2 != p0)
        {
            m_Neighs[p2].insert(p0);
            m_Neighs[p0].insert(p2);
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
}



void rmt::FlatUnion::ComputeTopologies()
{
    const Eigen::VectorXi& P = m_VPart.GetPartitions();

    int NVerts = m_Mesh.NumVertices();
    const Eigen::VectorXi& BV = m_Mesh.GetBoundaryVertices();
    for (int i = 0; i < NVerts; ++i)
    {
        if (BV[i])
            m_RDict.AddBoundaryVertex(P[i]);
        else
            m_RDict.AddVertex(P[i]);
    }

    const Eigen::MatrixXi& E = m_Mesh.GetEdges();
    const Eigen::VectorXi& BE = m_Mesh.GetBoundaryEdges();
    int NEdges = E.rows();
    for (int i = 0; i < NEdges; ++i)
    {
        if (BE[i])
            m_RDict.AddBoundaryEdge(P[E(i, 0)], P[E(i, 1)]);
        else
            m_RDict.AddEdge(P[E(i, 0)], P[E(i, 1)]);
    }

    const Eigen::MatrixXi& F = m_Mesh.GetTriangles();
    int NTris = F.rows();
    for (int i = 0; i < NTris; ++i)
        m_RDict.AddTriangle(P[F(i, 0)], P[F(i, 1)], P[F(i, 2)]);
}


bool rmt::FlatUnion::FixIssues()
{
    std::set<int> NewSamples;

    // For each sample with less than three adjacents, add samples to the boundary
    // until three adjacents are reached
    int NSamples = m_VPart.NumSamples();
    for (int i = 0; i < NSamples; ++i)
    {
        if (m_Neighs[i].size() >= 3)
            continue;
        if (std::get<0>(m_Farthests[i]) != m_VPart.GetSample(i))
            NewSamples.insert(std::get<0>(m_Farthests[i]));
        if (m_Neighs[i].size() < 2 && std::get<1>(m_Farthests[i]) != m_VPart.GetSample(i))
            NewSamples.insert(std::get<1>(m_Farthests[i]));
        if (m_Neighs[i].size() < 1 && std::get<2>(m_Farthests[i]) != m_VPart.GetSample(i))
            NewSamples.insert(std::get<2>(m_Farthests[i]));
    }


    // For each region that is not a closed 2-ball, break the connection
    size_t NRegions = m_RDict.NumRegions();
    for (size_t i = 0; i < NRegions; ++i)
    {
        if (m_RDict.IsClosed2Ball(i))
            continue;
        std::tuple<int, int, int> T = m_RDict.GetRegion(i).GetSamples();
        NewSamples.insert(m_Midpoints[T].begin(), m_Midpoints[T].end());
    }

    // Add the samples
    for (int v : NewSamples)
        m_VPart.AddSample(v);

    return NewSamples.size() == 0;
}