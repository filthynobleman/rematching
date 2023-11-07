/**
 * @file        eval.hpp
 * 
 * @brief       Functions for quality evaluation.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-27
 */
#pragma once


#include <Eigen/Dense>


namespace rmt
{

struct EvaluationMetrics
{
    double Hausdorff;
    double Chamfer;
    double MinArea;
    double MaxArea;
    double AvgArea;
    double StdArea;
    double MinQuality;
    double MaxQuality;
    double AvgQuality;
    double StdQuality;
};

    
double Hausdorff(const Eigen::MatrixXd& V,
                 const Eigen::MatrixXi& F,
                 const Eigen::MatrixXd& VRem,
                 const Eigen::MatrixXi& FRem);
double Chamfer(const Eigen::MatrixXd& V,
               const Eigen::MatrixXi& F,
               const Eigen::MatrixXd& VRem,
               const Eigen::MatrixXi& FRem);

EvaluationMetrics Evaluate(const Eigen::MatrixXd& V,
                           const Eigen::MatrixXi& F,
                           const Eigen::MatrixXd& VRem,
                           const Eigen::MatrixXi& FRem,
                           int NumOrigVerts);

} // namespace rmt
