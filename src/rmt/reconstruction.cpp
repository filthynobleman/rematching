/**
 * @file        reconstruction.cpp
 * 
 * @brief       Implements MeshFromVoronoi().
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-17
 */
#define NOMINMAX
#include <rmt/reconstruction.hpp>
#include <rmt/voronoifps.hpp>
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <unordered_map>
#include <unordered_set>
#include <set>

using namespace rmt;

typedef std::pair<int, std::pair<int, int>> Tri;

struct MyPairHash
{
    std::size_t operator()(const std::pair<int, int> v) const noexcept
    {
        return (size_t)v.first ^ ((size_t)v.second << 1); // or use boost::hash_combine
    }
};

void rmt::MeshFromVoronoi(const Graph & G, 
                          rmt::VoronoiPartitioning& Parts,
                          Eigen::MatrixXd& V, 
                          Eigen::MatrixXi& F)
{
    auto& Samples = Parts.Samples;
    auto& Partition = Parts.Partition;
    auto& Dists = Parts.Distances;

    // Retrieve vertices and remap indices
    int nSamples = Samples.size();
    V.resize(nSamples, 3);
    std::unordered_map<int, int> VMap;
    VMap.reserve(nSamples);
    for (int i = 0; i < nSamples; ++i)
    {
        V.row(i) = G.GetVertex(Samples[i]);
        VMap.emplace(Samples[i], i);
    }

    // Determine edges
    std::set<std::pair<int, int>> Edges;
    int nVerts = G.NumVertices();
    for (int i = 0; i < nVerts; ++i)
    {
        int nAdjs = G.NumAdjacents(i);
        for (int jj = 0; jj < nAdjs; ++jj)
        {
            int j = G.GetAdjacent(i, jj).first;
            if (Partition[i] != Partition[j])
                Edges.emplace(VMap[Partition[i]], VMap[Partition[j]]);
        }
    }

    // Create geodesic Delaunay graph
    Graph DG(V, Edges);

    // Find all 3-cycles in graph
    std::set<Tri> Tris;
    for (int i = 0; i < nSamples; ++i)
    {
        int nAdjsI = DG.NumAdjacents(i);
        for (int jj = 0; jj < nAdjsI; ++jj)
        {
            int j = DG.GetAdjacent(i, jj).first;
            // If j < i, we already covered the case when iterating over j
            if (j < i)
                continue;
            int nAdjsJ = DG.NumAdjacents(j);
            for (int kk = jj + 1; kk < nAdjsI; ++kk)
            {
                int k = DG.GetAdjacent(i, kk).first;
                for (int hh = 0; hh < nAdjsJ; ++hh)
                {
                    int h = DG.GetAdjacent(j, hh).first;
                    if (h == k)
                    {
                        // i, j, k is a triangle
                        Tri t(i, { j, k });
                        Tris.insert(t);
                    }
                }
            }
        }
    }

    // Find all non-manifold edges
    std::unordered_map<std::pair<int, int>, int, MyPairHash> ECount;
    ECount.reserve(Samples.size() * 3);
    std::unordered_set<std::pair<int, int>, MyPairHash> NMEdges;
    NMEdges.reserve(Samples.size() * 3);
    std::pair<int, int> ETmp[3];
    for (auto t : Tris)
    {
        ETmp[0] = { std::min(t.first, t.second.first), std::max(t.first, t.second.first) };
        ETmp[1] = { std::min(t.second.first, t.second.second), std::max(t.second.first, t.second.second) };
        ETmp[2] = { std::min(t.first, t.second.second), std::max(t.first, t.second.second) };

        for (int j = 0; j < 3; ++j)
        {
            if (ECount.find(ETmp[j]) == ECount.end())
                ECount.emplace(ETmp[j], 0);
            ECount[ETmp[j]] += 1;
            if (ECount[ETmp[j]] > 2)
                NMEdges.insert(ETmp[j]);
        }
    }

    // Compute the manifold triangles
    if (NMEdges.size() > 0)
    {
        std::set<Tri> ManTris;
        for (auto t : Tris)
        {
            ETmp[0] = { std::min(t.first, t.second.first), std::max(t.first, t.second.first) };
            ETmp[1] = { std::min(t.second.first, t.second.second), std::max(t.second.first, t.second.second) };
            ETmp[2] = { std::min(t.first, t.second.second), std::max(t.first, t.second.second) };

            // Count incident non-manifold edges
            int NMCount = 0;
            for (int j = 0; j < 3; ++j)
            {
                if (NMEdges.find(ETmp[j]) != NMEdges.end())
                    NMCount += 1;
            }

            // If all its edges are non-manifold, it is a repeated patch and must be removed
            if (NMCount == 3)
                continue;

            // Otherwise, we keep it for now
            ManTris.insert(t);
        }
        Tris = std::move(ManTris);
    }


    F.resize(Tris.size(), 3);
    int i = 0;
    for (auto it = Tris.begin(); it != Tris.end(); it++)
        F.row(i++) = Eigen::Vector3i(it->first, it->second.first, it->second.second);
}

void rmt::Refine(const Eigen::MatrixXd& VOld,
                 const Eigen::MatrixXi& FOld,
                 const rmt::Graph& G,
                 rmt::VoronoiPartitioning& Parts,
                 Eigen::MatrixXd& V,
                 Eigen::MatrixXi& F)
{
    do { } 
    while (rmt::_RefineStep(VOld, FOld, G, Parts, V, F));
}


bool rmt::_RefineStep(const Eigen::MatrixXd& VOld,
                      const Eigen::MatrixXi& FOld,
                      const rmt::Graph& G,
                      rmt::VoronoiPartitioning& Parts,
                      Eigen::MatrixXd& V,
                      Eigen::MatrixXi& F)
{
    auto& Samples = Parts.Samples;
    auto& Partition = Parts.Partition;
    auto& Dists = Parts.Distances;

    // Reorient faces
    rmt::ReorientFaces(Samples, VOld, FOld, V, F);

    // Check non-manifoldness
    bool IsManifold = true;
    // Check for non-manifold edges
    Eigen::MatrixXi E, BF;
    Eigen::VectorXi EMAP, EM;
    if (!igl::is_edge_manifold(F, BF, E, EMAP, EM))
    {
        IsManifold = false;
        std::cout << "Non manifold edges: " << (1 - EM.array()).sum() << std::endl;
        cut::MinHeap HM(Dists, true);
        auto OldP = Partition;
        for (int i = 0; i < EM.rows(); ++i)
        {
            if (EM[i]) 
                continue;

            auto PE = G.DijkstraPath(Samples[E(i, 0)], Samples[E(i, 1)]);
            double l = 0;
            int newp = -1;
            for (auto we : PE.second)
            {
                l += we.second;
                if (2 * l >= PE.first)
                {
                    newp = we.first;
                    break;
                }
            }
            CUTAssert(newp != -1);
            CUTAssert(newp != Samples[E(i, 0)]);
            CUTAssert(newp != Samples[E(i, 1)]);
            Samples.emplace_back(newp);
            rmt::UpdateVoronoi(G, Dists, Partition, HM, newp);
        }
    }

    // Check for non-manifold vertices
    Eigen::VectorXi VM;
    if (!igl::is_vertex_manifold(F, VM))
    {
        IsManifold = false;
        std::cout << "Non manifold vertices: " << (1 - VM.array()).sum() << std::endl;
        cut::MinHeap HM(Dists, true);
        for (int i = 0; i < VM.rows(); ++i)
        {
            if (VM[i])
                continue;
            
            int newp = G.FarthestFiltered(Samples[i], Partition, Samples[i]);
            CUTAssert(newp != Samples[i]);
            Samples.emplace_back(newp);
            rmt::UpdateVoronoi(G, Dists, Partition, HM, newp);
        }
    }

    // Try reconstructing again
    if (!IsManifold)
        rmt::MeshFromVoronoi(G, Parts, V, F);
    // If input is manifold, no changes are required
    return !IsManifold;
}



void rmt::ReorientFaces(const std::vector<int>& Samples,
                               const Eigen::MatrixXd& VOld,
                               const Eigen::MatrixXi& FOld,
                               const Eigen::MatrixXd& VNew,
                               Eigen::MatrixXi& FNew)
{
    // Compute normals of original mesh
    int noFaces = FOld.rows();
    Eigen::MatrixXd NOld;
    // igl::per_face_normals_stable(VOld, FOld, NOld);
    igl::per_vertex_normals(VOld, FOld, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, NOld);
    // NOld.resize(VOld.rows(), 3);
    // NOld.setZero();
    // Eigen::Vector3d Tmp1, Tmp2;
    // for (int i = 0; i < noFaces; ++i)
    // {
    //     Eigen::Vector3d n(0.0, 0.0, 0.0);
    //     Tmp1 = (VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).transpose();
    //     Tmp2 = (VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).transpose();
    //     n += Tmp1.cross(Tmp2).normalized();

    //     Tmp1 = (VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).transpose();
    //     Tmp2 = (VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).transpose();
    //     n += Tmp1.cross(Tmp2).normalized();

    //     Tmp1 = (VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).transpose();
    //     Tmp2 = (VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).transpose();
    //     n += Tmp1.cross(Tmp2).normalized();


    //     // n += (VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).cross(VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).transpose().normalized();
    //     // n += (VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).cross(VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).transpose().normalized();
    //     // n += (VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).cross(VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).transpose().normalized();
    //     n /= 3.0;
    //     NOld.row(FOld(i, 0)) += n;
    //     NOld.row(FOld(i, 1)) += n;
    //     NOld.row(FOld(i, 2)) += n;
    // }


    // Compute normals of new faces and reorient if needed
    int nFaces = FNew.rows();
    Eigen::MatrixXd N;
    igl::per_face_normals_stable(VNew, FNew, N);
    // N.resize(nFaces, 3);
    // N.setZero();
    for (int i = 0; i < nFaces; ++i)
    {
        // Eigen::Vector3d v[3];
        // for (int j = 0; j < 3; ++j)
        //     v[j] = VNew.row(FNew(i, j)).segment<3>(0);
        // Eigen::Vector3d e01, e12, e20;
        // e01 = v[1] - v[0];
        // e12 = v[2] - v[1];
        // e20 = v[0] - v[2];
        Eigen::Vector3d n(0.0, 0.0, 0.0);
        // n += e01.cross(e12).normalized();
        // n += e12.cross(e20).normalized();
        // n += e20.cross(e01).normalized();
        // n /= 3.0;
        // N.row(i) = n.transpose();

        n = NOld.row(Samples[FNew(i, 0)]) + NOld.row(Samples[FNew(i, 1)]) + NOld.row(Samples[FNew(i, 2)]);
        n /= 3.0;

        if (n.dot(N.row(i)) < 0.0)
        {
            FNew.row(i).reverseInPlace();
            N.row(i) = -N.row(i);
        }
    }
}