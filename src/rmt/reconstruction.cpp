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
#include <igl/boundary_loop.h>
#include <igl/euler_characteristic.h>
#include <igl/edge_flaps.h>
#include <igl/unique_edge_map.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <unordered_map>
#include <unordered_set>
#include <set>


using namespace rmt;

typedef std::tuple<int, int, int> Tri;

struct MyPairHash
{
    std::size_t operator()(const std::pair<int, int>& v) const noexcept
    {
        return std::hash<int>{}(v.first) ^ (std::hash<int>{}(v.second) << 1); // or use boost::hash_combine
    }
};

struct MyTriHash
{
    std::size_t operator()(const Tri& t) const noexcept
    {
        return std::hash<int>{}(std::get<0>(t)) ^ ((std::get<1>(t) ^ (std::get<2>(t) << 1)) << 1);
    }
};

void rmt::MeshFromVoronoi(const Eigen::MatrixXd& VOld,
                          const Eigen::MatrixXi& FOld,
                          rmt::VoronoiPartitioning& Parts,
                          Eigen::MatrixXd& V,
                          Eigen::MatrixXi& F)
{
    auto& Samples = Parts.GetSamples();
    auto& Partition = Parts.GetPartitions();
    auto& Dists = Parts.GetDistances();

    // Retrieve vertices and remap indices
    int nSamples = Samples.size();
    V.resize(nSamples, 3);
    for (int i = 0; i < nSamples; ++i)
        V.row(i) = VOld.row(Samples[i]);

    // Iterate over original faces and compute triangles incident on three partitions
    std::unordered_set<Tri, MyTriHash> Tris;
    Tris.reserve(3 * nSamples);
    for (int i = 0; i < FOld.rows(); ++i)
    {
        if (Partition[FOld(i, 0)] == Partition[FOld(i, 1)])
            continue;
        if (Partition[FOld(i, 1)] == Partition[FOld(i, 2)])
            continue;
        if (Partition[FOld(i, 2)] == Partition[FOld(i, 0)])
            continue;

        Tri t = std::tie(Partition[FOld(i, 0)], Partition[FOld(i, 1)], Partition[FOld(i, 2)]);
        if (std::get<1>(t) < std::get<0>(t) && std::get<1>(t) < std::get<2>(t))
            t = Tri(std::get<1>(t), std::get<2>(t), std::get<0>(t));
        else if (std::get<2>(t) < std::get<0>(t) && std::get<2>(t) < std::get<1>(t))
            t = Tri(std::get<2>(t), std::get<0>(t), std::get<1>(t));

        Tris.insert(t);
    }

    // Create the new faces
    F.resize(Tris.size(), 3);
    int i = 0;
    for (const Tri& t : Tris)
        F.row(i++) = Eigen::RowVector3i(std::get<0>(t), std::get<1>(t), std::get<2>(t));


    
    // Eigen::MatrixXi E;
    // Eigen::MatrixXi uE;
    // Eigen::VectorXi EMAP;
    // Eigen::MatrixXi EF;
    // Eigen::MatrixXi EI;
    // std::vector<std::vector<int>> uE2E;
    // igl::unique_edge_map(F, E, uE, EMAP, uE2E);
    // // igl::edge_flaps(F, (const Eigen::MatrixXi)uE, (const Eigen::VectorXi)EMAP, EF, EI);
    // Eigen::MatrixXd l;
    // igl::edge_lengths(V, F, l);
    // Eigen::MatrixXd l_intrinsic;
    // Eigen::MatrixXi F_intrinsic;
    // igl::intrinsic_delaunay_triangulation(l, F, l_intrinsic, F_intrinsic, E, uE, EMAP, uE2E);

    // Eigen::MatrixXd SV = V;
    // Eigen::VectorXi SVI, SVJ;
    // igl::remove_duplicate_vertices(SV, F_intrinsic, 1e-2 * l.mean(), V, SVI, SVJ, F);
    // SV = V;
    // F_intrinsic.resize(F.rows(), 3);
    // int LastFace = 0;
    // for (int i = 0; i < F.rows(); ++i)
    // {
    //     if (F(i, 0) == F(i, 1) || F(i, 1) == F(i, 2) || F(i, 2) == F(i, 0))
    //         continue;
    //     F_intrinsic.row(LastFace++) = F.row(i);
    // }
    // igl::remove_unreferenced(SV, F_intrinsic, V, F, SVI, SVJ);
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