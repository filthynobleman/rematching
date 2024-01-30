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
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>


void rmt::Remesh(const Eigen::MatrixXd & Vin, 
                 const Eigen::MatrixXi & Fin, 
                 int NSamples, 
                 Eigen::MatrixXd & Vout, 
                 Eigen::MatrixXi & Fout)
{
    rmt::Mesh M(Vin, Fin);
    M.ComputeEdgesAndBoundaries();
    rmt::VoronoiPartitioning VPart(M);
    while (VPart.NumSamples() < NSamples)
        VPart.AddSample(VPart.FarthestVertex());
    if (igl::is_edge_manifold(Fin) && igl::is_vertex_manifold(Fin))
    {
        rmt::FlatUnion FU(M, VPart);
        do
        {
            FU.DetermineRegions();
            FU.ComputeTopologies();
        } while (!FU.FixIssues());
    }
    
    rmt::MeshFromVoronoi(Vin, Fin, VPart, Vout, Fout);
}

void rmt::Remesh(const Eigen::MatrixXd & Vin, 
                 const Eigen::MatrixXi & Fin, 
                 int NSamples, 
                 Eigen::MatrixXd & Vout, 
                 Eigen::MatrixXi & Fout,
                 Eigen::VectorXi & Vidx)
{
    rmt::Mesh M(Vin, Fin);
    M.ComputeEdgesAndBoundaries();
    rmt::VoronoiPartitioning VPart(M);
    while (VPart.NumSamples() < NSamples)
        VPart.AddSample(VPart.FarthestVertex());
    if (igl::is_edge_manifold(Fin) && igl::is_vertex_manifold(Fin))
    {
        rmt::FlatUnion FU(M, VPart);
        do
        {
            FU.DetermineRegions();
            FU.ComputeTopologies();
        } while (!FU.FixIssues());
    }
    
    rmt::MeshFromVoronoi(Vin, Fin, VPart, Vout, Fout);
    Vidx.setZero(VPart.NumSamples());
    for (int i = 0; i < Vidx.rows(); ++i)
        Vidx[i] = VPart.GetSample(i);
}