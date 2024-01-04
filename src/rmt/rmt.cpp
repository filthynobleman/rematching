/**
 * @file        rmt.cpp
 * 
 * @brief       ReMatching exposed function.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#include <rmt/rmt.hpp>


void rmt::Remesh(const Eigen::MatrixXd & Vin, 
                 const Eigen::MatrixXi & Fin, 
                 int NSamples, 
                 Eigen::MatrixXd & Vout, 
                 Eigen::MatrixXi & Fout)
{
    rmt::Graph Graph(Vin, Fin);
    auto VFPS = rmt::VoronoiFPS(Graph, NSamples);
    rmt::MeshFromVoronoi(Graph, VFPS, Vout, Fout);
    rmt::Refine(Vin, Fin, Graph, VFPS, Vout, Fout);
    // rmt::ReorientFaces(VFPS.Samples, Vin, Fin, Vout, Fout);
}