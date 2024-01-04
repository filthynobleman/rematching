/**
 * @file        reconstruction.hpp
 * 
 * @brief       Building mesh from VoronoiFPS.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-17
 */
#pragma once

#include <rmt/graph.hpp>
#include <rmt/voronoifps.hpp>
#include <vector>


namespace rmt
{
    
void MeshFromVoronoi(const Graph& G,
                     rmt::VoronoiPartitioning& Parts,
                     Eigen::MatrixXd& V,
                     Eigen::MatrixXi& F);

void Refine(const Eigen::MatrixXd& VOld,
            const Eigen::MatrixXi& FOld,
            const rmt::Graph& G,
            rmt::VoronoiPartitioning& Parts,
            Eigen::MatrixXd& V,
            Eigen::MatrixXi& F);

bool _RefineStep(const Eigen::MatrixXd& VOld,
                 const Eigen::MatrixXi& FOld,
                 const rmt::Graph& G,
                 rmt::VoronoiPartitioning& Parts,
                 Eigen::MatrixXd& V,
                 Eigen::MatrixXi& F);

void ReorientFaces(const std::vector<int>& Samples,
                   const Eigen::MatrixXd& VOld,
                   const Eigen::MatrixXi& FOld,
                   const Eigen::MatrixXd& VNew,
                   Eigen::MatrixXi& FNew);

} // namespace rmt
