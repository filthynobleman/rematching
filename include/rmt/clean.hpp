/**
 * @file        clean.hpp
 * 
 * @brief       Functions for proicessing and cleaning meshes.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-30
 */
#pragma once

#include <Eigen/Dense>

namespace rmt
{

void MakeManifold(Eigen::MatrixXd& V,
                  Eigen::MatrixXi& F);

void RemoveSmallComponents(Eigen::MatrixXd& V,
                           Eigen::MatrixXi& F,
                           double AreaFraction = 1e-2);

void RemoveDegeneracies(Eigen::MatrixXd& V,
                        Eigen::MatrixXi& F,
                        double DistThreshold = 1e-4);

void CleanUp(Eigen::MatrixXd& V,
             Eigen::MatrixXi& F,
             double AreaFraction = 1e-2,
             double DistanceThreshold = 1e-4);
    
} // namespace rmt
