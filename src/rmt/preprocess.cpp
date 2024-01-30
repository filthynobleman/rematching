/**
 * @file        preprocess.cpp
 * 
 * @brief       Implements MaxEdgeLength() and ResampleMesh().
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-26
 */
#include <rmt/preprocess.hpp>
#include <rmt/utils.hpp>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <queue>
#define NOMINMAX
#include <igl/doublearea.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/facet_components.h>

struct MyPairHash
{
    std::size_t operator()(Eigen::Vector2i const& v) const noexcept
    {
        return (size_t)v[0] ^ ((size_t)v[1] << 1); // or use boost::hash_combine
    }
};

void rmt::RescaleInsideUnitBox(Eigen::MatrixXd& V)
{
    int nVerts = V.rows();
    double L = V.cwiseAbs().maxCoeff();
    V /= L;
}


double rmt::MaxEdgeLength(const Eigen::MatrixXd& V,
                                 const Eigen::MatrixXi& F,
                                 int OutputSize)
{
    Eigen::VectorXd TriAreas;
    igl::doublearea(V, F, TriAreas);
    // We expect output mesh to have about F = 2 * V
    // So, average triangle area is A / (2 * V)
    // We consider an extra factor of two because we computed double areas
    double AvgArea = TriAreas.sum() / (4 * OutputSize);
    // We expect/aim to equilateral triangles, whose area is
    // l^2 * sqrt(3) / 4. So edge length is 2 * sqrt(A / sqrt(3))
    return 2 * std::sqrt(AvgArea / std::sqrt(3));
}

void rmt::ResampleMesh(Eigen::MatrixXd &V, 
                              Eigen::MatrixXi &F, 
                              double MaxEdgeLen)
{
    // Find edges to split
    std::unordered_map<Eigen::Vector2i, int, MyPairHash> EdgeMap;
    int nFaces = F.rows();
    int nVerts = V.rows();
    std::vector<Eigen::Vector3d> NewV;
    EdgeMap.reserve(3 * nFaces);
    NewV.reserve(3 * nFaces);
    for (int i = 0; i < nFaces; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int j0 = j;
            int j1 = (j + 1) % 3;
            const Eigen::RowVector3d& v0 = V.row(F(i, j0));
            const Eigen::RowVector3d& v1 = V.row(F(i, j1));
            Eigen::Vector2i E(std::min(F(i, j0), F(i, j1)), std::max(F(i, j0), F(i, j1)));
            if (EdgeMap.find(E) != EdgeMap.end())
                continue;
            double ELen = (v1 - v0).norm();
            if (ELen > MaxEdgeLen)
            {
                EdgeMap.emplace(E, (int)NewV.size());
                NewV.emplace_back(0.5 * (v0 + v1));
            }
            else
            {
                EdgeMap.emplace(E, -1);
            }
        }
    }

    if (NewV.size() == 0)
        return;

    // Add the vertices
    V.conservativeResize(nVerts + NewV.size(), 3);
    int NewVLen = NewV.size();
    for (int i = 0; i < NewVLen; ++i)
        V.row(nVerts + i) = NewV[i];

    // Resample faces
    std::vector<Eigen::RowVector3i> NewF;
    NewF.reserve(3 * nFaces);
    int NRes = 0;
    for (int i = 0; i < nFaces; ++i)
    {
        // Which edges to split in the triangle
        int NewVerts[3] = { -1, -1, -1 };
        for (int j = 0; j < 3; ++j)
        {
            int j0 = j;
            int j1 = (j + 1) % 3;
            int j2 = (j + 2) % 3;
            Eigen::Vector2i E(std::min(F(i, j2), F(i, j1)), std::max(F(i, j2), F(i, j1)));
            NewVerts[j] = EdgeMap[E];
        }

        int NumSplits = 0;
        for (int j = 0; j < 3; ++j)
        {
            if (NewVerts[j] == -1)
                continue;
            NewVerts[j] += nVerts;
            NumSplits += 1;
        }

        if (NumSplits == 0)
        {
            // No splits, triangle is good as is
            continue;
        }

        Eigen::Vector3i f = F.row(i);
        if (NumSplits == 3)
        {
            // All splits, divide into 4 triangles
            F.row(i) = Eigen::RowVector3i(NewVerts[0], NewVerts[1], NewVerts[2]);
            NewF.emplace_back(NewVerts[2], NewVerts[1], f[0]);
            NewF.emplace_back(NewVerts[0], NewVerts[2], f[1]);
            NewF.emplace_back(NewVerts[1], NewVerts[0], f[2]);
        }
        else if (NumSplits == 1)
        {
            // One split, connect to opposite vertex
            for (int j = 0; j < 3; ++j)
            {
                if (NewVerts[j] == -1)
                    continue;
                
                F.row(i) = Eigen::RowVector3i(f[j], f[(j + 1) % 3], NewVerts[j]);
                NewF.emplace_back(f[(j + 2) % 3], f[j], NewVerts[j]);
            }
        }
        else if (NumSplits == 2)
        {
            // Two long edges, connect them and cut the quad along a diagonal
            int j0 = 0;
            int j1 = 1;
            int j2 = 2;
            if (NewVerts[0] == -1)
            {
                j0 = 1;
                j1 = 2;
                j2 = 0;
            }
            else if (NewVerts[1] == -1)
            {
                j0 = 2;
                j1 = 0;
                j2 = 1;
            }

            // Connect the long edges to form a new triangle
            F.row(i) = Eigen::RowVector3i(NewVerts[j1], NewVerts[j0], f[j2]);
            double e0 = (V.row(f[j0]) - V.row(NewVerts[j0])).norm();
            double e1 = (V.row(f[j1]) - V.row(NewVerts[j1])).norm();
            if (e1 < e0)
            {
                NewF.emplace_back(f[j0], f[j1], NewVerts[j1]);
                NewF.emplace_back(f[j1], NewVerts[j0], NewVerts[j1]);
            }
            else
            {
                NewF.emplace_back(f[j0], f[j1], NewVerts[j0]);
                NewF.emplace_back(f[j0], NewVerts[j0], NewVerts[j1]);
            }
        }
    }

    // Add all the new faces to the mesh
    F.conservativeResize(nFaces + NewF.size(), 3);
    int NewFLen = NewF.size();
    for (int i = 0; i < NewFLen; ++i)
        F.row(nFaces + i) = NewF[i];

    // Call recursively
    ResampleMesh(V, F, MaxEdgeLen);
}