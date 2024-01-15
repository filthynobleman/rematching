/**
 * @file        eval.cpp
 * 
 * @brief       Implements evaluation metrics.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-27
 */
#define NOMINMAX
#include <rmt/eval.hpp>
#include <assert.h>
#include <igl/hausdorff.h>
#include <igl/point_mesh_squared_distance.h>

using namespace rmt;


double rmt::Hausdorff(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXd& VRem,
                             const Eigen::MatrixXi& FRem)
{
    double HD = std::numeric_limits<double>::infinity();
    igl::hausdorff(V, F, VRem, FRem, HD);

    return HD;
}

double rmt::Chamfer(const Eigen::MatrixXd& V,
                           const Eigen::MatrixXi& F,
                           const Eigen::MatrixXd& VRem,
                           const Eigen::MatrixXi& FRem)
{
    double CD = 0.0;

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    igl::point_mesh_squared_distance(V, VRem, FRem, sqrD, I, C);
    CD += sqrD.sum() / V.rows();

    igl::point_mesh_squared_distance(VRem, V, F, sqrD, I, C);
    CD += sqrD.sum() / VRem.rows();

    return CD;
}


EvaluationMetrics rmt::Evaluate(const Eigen::MatrixXd& V,
                                       const Eigen::MatrixXi& F,
                                       const Eigen::MatrixXd& VRem,
                                       const Eigen::MatrixXi& FRem,
                                       int NumOrigVerts)
{
    static const double one_sixth = 1 / 6.0;
    static const double two_sqrt3 = 2.0 * std::sqrt(3.0);

    EvaluationMetrics Metrics;
    auto VV = V.block(0, 0, NumOrigVerts, V.cols());

    Metrics.Hausdorff = Hausdorff(VV, F, VRem, FRem);
    Metrics.Chamfer = Chamfer(VV, F, VRem, FRem);

    Metrics.MinArea = std::numeric_limits<double>::infinity();
    Metrics.MaxArea = - std::numeric_limits<double>::infinity();
    Metrics.AvgArea = 0.0;
    Metrics.StdArea = 0.0;
    Metrics.MinQuality = std::numeric_limits<double>::infinity();
    Metrics.MaxQuality = - std::numeric_limits<double>::infinity();
    Metrics.AvgQuality = 0.0;
    Metrics.StdQuality = 0.0;

    int nFaces = FRem.rows();
    std::vector<double> Areas;
    std::vector<double> Qualities;
    Areas.resize(nFaces);
    Qualities.resize(nFaces);
    for (int i = 0; i < nFaces; ++i)
    {
        Eigen::Vector3d e01 = VRem.row(FRem(i, 1)) - VRem.row(FRem(i, 0));
        Eigen::Vector3d e12 = VRem.row(FRem(i, 2)) - VRem.row(FRem(i, 1));
        Eigen::Vector3d e20 = VRem.row(FRem(i, 0)) - VRem.row(FRem(i, 2));

        double A = 0.0;
        A += e01.cross(e12).norm();
        A += e12.cross(e20).norm();
        A += e20.cross(e01).norm();
        A *= one_sixth;
        Areas[i] = A;

        double l01 = e01.norm();
        double l12 = e12.norm();
        double l20 = e20.norm();
        double s = 0.5 * (l01 + l12 + l20);
        double lm = std::max(l01, std::max(l12, l20));
        double Q = two_sqrt3 * A / (s * lm);
        Qualities[i] = Q;

        Metrics.MinArea = std::min(Metrics.MinArea, A);
        Metrics.MaxArea = std::max(Metrics.MaxArea, A);
        Metrics.AvgArea += A;

        Metrics.MinQuality = std::min(Metrics.MinQuality, Q);
        Metrics.MaxQuality = std::max(Metrics.MaxQuality, Q);
        Metrics.AvgQuality += Q;
    }

    Metrics.AvgArea /= nFaces;
    Metrics.AvgQuality /= nFaces;

    for (int i = 0; i < nFaces; ++i)
    {
        double d = Metrics.AvgArea - Areas[i];
        Metrics.StdArea += d * d;

        d = Metrics.AvgQuality - Qualities[i];
        Metrics.StdQuality += d * d;
    }
    Metrics.StdArea = std::sqrt(Metrics.StdArea / nFaces);
    Metrics.StdQuality = std::sqrt(Metrics.StdQuality / nFaces);

    return Metrics;
}

