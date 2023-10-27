/**
 * @file        preprocess.hpp
 * 
 * @brief       Declaration of preprocessing functions.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-26
 */
#pragma once

#include <Eigen/Dense>


namespace rmt
{

double MaxEdgeLength(const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     int OutputSize);

void ResampleMesh(Eigen::MatrixXd& V,
                  Eigen::MatrixXi& F,
                  double MaxEdgeLen);

void RescaleInsideUnitBox(Eigen::MatrixXd& V);

} // namespace rmt
