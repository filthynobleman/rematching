/**
 * @file        weightmap.hpp
 * 
 * @brief       Declaration of the function that builds the mapping from old to new mesh.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-25
 */
#pragma once

#include <rmt/graph.hpp>
#include <Eigen/Sparse>

namespace rmt
{
    
Eigen::SparseMatrix<double> WeightMap(const Eigen::MatrixXd& P,
                                      const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      int nVerts = -1);

} // namespace rmt
