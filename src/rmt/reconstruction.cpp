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
#include <rmt/reconstruction.hpp>
#include <unordered_map>
#include <set>

using namespace rmt;

typedef std::pair<int, std::pair<int, int>> Tri;

void rmt::MeshFromVoronoi(const Graph & G, 
                          const std::vector<int>& Samples, 
                          const std::vector<int>& Partition, 
                          Eigen::MatrixXd& V, 
                          Eigen::MatrixXi& F)
{
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

    F.resize(Tris.size(), 3);
    int i = 0;
    for (auto it = Tris.begin(); it != Tris.end(); it++)
        F.row(i++) = Eigen::Vector3i(it->first, it->second.first, it->second.second);
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
    NOld.resize(VOld.rows(), 3);
    NOld.setZero();
    Eigen::Vector3d Tmp1, Tmp2;
    for (int i = 0; i < noFaces; ++i)
    {
        Eigen::Vector3d n(0.0, 0.0, 0.0);
        Tmp1 = (VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).transpose();
        Tmp2 = (VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).transpose();
        n += Tmp1.cross(Tmp2).normalized();

        Tmp1 = (VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).transpose();
        Tmp2 = (VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).transpose();
        n += Tmp1.cross(Tmp2).normalized();

        Tmp1 = (VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).transpose();
        Tmp2 = (VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).transpose();
        n += Tmp1.cross(Tmp2).normalized();


        // n += (VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).cross(VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).transpose().normalized();
        // n += (VOld.row(FOld(i, 2)) - VOld.row(FOld(i, 1))).segment<3>(0, 3).cross(VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).transpose().normalized();
        // n += (VOld.row(FOld(i, 0)) - VOld.row(FOld(i, 2))).segment<3>(0, 3).cross(VOld.row(FOld(i, 1)) - VOld.row(FOld(i, 0))).segment<3>(0, 3).transpose().normalized();
        n /= 3.0;
        NOld.row(FOld(i, 0)) += n;
        NOld.row(FOld(i, 1)) += n;
        NOld.row(FOld(i, 2)) += n;
    }


    // Compute normals of new faces and reorient if needed
    int nFaces = FNew.rows();
    Eigen::MatrixXd N;
    N.resize(nFaces, 3);
    N.setZero();
    for (int i = 0; i < nFaces; ++i)
    {
        Eigen::Vector3d v[3];
        for (int j = 0; j < 3; ++j)
            v[j] = VNew.row(FNew(i, j)).segment<3>(0);
        Eigen::Vector3d e01, e12, e20;
        e01 = v[1] - v[0];
        e12 = v[2] - v[1];
        e20 = v[0] - v[2];
        Eigen::Vector3d n(0.0, 0.0, 0.0);
        n += e01.cross(e12).normalized();
        n += e12.cross(e20).normalized();
        n += e20.cross(e01).normalized();
        n /= 3.0;
        N.row(i) = n.transpose();

        n = NOld.row(Samples[FNew(i, 0)]) + NOld.row(Samples[FNew(i, 1)]) + NOld.row(Samples[FNew(i, 2)]);
        n /= 3.0;

        if (n.dot(N.row(i)) < 0.0)
        {
            FNew.row(i).reverseInPlace();
            N.row(i) = -N.row(i);
        }
    }
}