/**
 * @file        io.hpp
 * 
 * @brief       Functions for loading and exporting data.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace rmt
{

bool LoadMesh(const std::string& Filename,
              Eigen::MatrixXd& V,
              Eigen::MatrixXi& F);

bool ExportMesh(const std::string& Filename,
                const Eigen::MatrixXd& V,
                const Eigen::MatrixXi& F);


bool ExportWeightmap(const std::string& Filename,
                     const Eigen::SparseMatrix<double>& WM);
    
} // namespace rmt
