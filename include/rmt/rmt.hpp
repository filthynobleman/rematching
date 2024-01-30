/**
 * @file        rmt.hpp
 * 
 * @brief       Include whole library.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-31
 */
#pragma once

#include <rmt/graph.hpp>
#include <rmt/mesh.hpp>
#include <rmt/preprocess.hpp>
#include <rmt/voronoifps.hpp>
#include <rmt/flatunion.hpp>
#include <rmt/reconstruction.hpp>
#include <rmt/weightmap.hpp>
#include <rmt/eval.hpp>
#include <rmt/io.hpp>
#include <rmt/clean.hpp>

#include <cassert>


namespace rmt
{
    
void Remesh(const Eigen::MatrixXd& Vin,
            const Eigen::MatrixXi& Fin,
            int NSamples,
            Eigen::MatrixXd& Vout,
            Eigen::MatrixXi& Fout);

void Remesh(const Eigen::MatrixXd& Vin,
            const Eigen::MatrixXi& Fin,
            int NSamples,
            Eigen::MatrixXd& Vout,
            Eigen::MatrixXi& Fout,
            Eigen::VectorXi& Vidx);

} // namespace rmt
