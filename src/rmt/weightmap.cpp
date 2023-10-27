/**
 * @file        weightmap.cpp
 * 
 * @brief       Implements rmt::WeightMap().
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-25
 */
#include <rmt/weightmap.hpp>
#include <rmt/voronoifps.hpp>
#include <set>
#include <unordered_map>
#include <queue>
#include <assert.h>

#include <igl/barycentric_coordinates.h>
#include <igl/point_mesh_squared_distance.h>


Eigen::SparseMatrix<double> rmt::WeightMap(const Eigen::MatrixXd& P,
                                           const Eigen::MatrixXd& V,
                                           const Eigen::MatrixXi& F,
                                           int nVerts)
{
    if (nVerts < 0)
        nVerts = P.rows();

    Eigen::SparseMatrix<double> W;
    W.resize(nVerts, V.rows());

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

    std::vector<Eigen::Triplet<double>> Triplets;
    Triplets.reserve(3 * nVerts);
    for (int i = 0; i < nVerts; ++i)
    {
        Eigen::RowVector3d L;
        igl::barycentric_coordinates(C.row(i),
                                     V.row(F(I[i], 0)),
                                     V.row(F(I[i], 1)),
                                     V.row(F(I[i], 2)),
                                     L);
        for (int j = 0; j < 3; ++j)
            Triplets.emplace_back(i, F(I[i], j), L[j]);
    }

    W.setFromTriplets(Triplets.begin(), Triplets.end());
    return W;
}