/**
 * @file        rmt.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#include <rmt/mex/rmt.hpp>
#include <rmt/rmt.hpp>

void rmt::Remesh(const double* Vin, 
                 int NVerts,
                 const int* Fin, 
                 int NFaces,
                 int NSamples, 
                 double** Vout, 
                 size_t* NVout,
                 int** Fout,
                 size_t* NFout)
{
    Eigen::MatrixXd V(Eigen::MatrixXd::Map(Vin, NVerts, 3));
    Eigen::MatrixXi F(Eigen::MatrixXi::Map(Fin, NFaces, 3));
    Eigen::MatrixXd VV;
    Eigen::MatrixXi FF;
    F = F.array() - 1;
    rmt::Remesh(V, F, NSamples, VV, FF);
    FF = FF.array() + 1;

    *Vout = (double*)std::malloc(VV.size() * sizeof(double));
    std::memcpy(*Vout, VV.data(), VV.size() * sizeof(double));
    *NVout = VV.rows();
    *Fout = (int*)std::malloc(FF.size() * sizeof(int));
    std::memcpy(*Fout, FF.data(), FF.size() * sizeof(int));
    *NFout = FF.rows();
}


void rmt::Resample(const double* Vin, 
                   int NVerts, 
                   const int* Fin, 
                   int NFaces, 
                   int NSamples, 
                   double** Vout, 
                   size_t* NVout, 
                   int** Fout, 
                   size_t* NFout)
{
    Eigen::MatrixXd V(Eigen::MatrixXd::Map(Vin, NVerts, 3));
    Eigen::MatrixXi F(Eigen::MatrixXi::Map(Fin, NFaces, 3));
    F = F.array() - 1;
    double MaxEdge = rmt::MaxEdgeLength(V, F, NSamples);
    rmt::ResampleMesh(V, F, MaxEdge);
    F = F.array() + 1;

    *Vout = (double*)std::malloc(V.size() * sizeof(double));
    std::memcpy(*Vout, V.data(), V.size() * sizeof(double));
    *NVout = V.rows();
    *Fout = (int*)std::malloc(F.size() * sizeof(int));
    std::memcpy(*Fout, F.data(), F.size() * sizeof(int));
    *NFout = F.rows();
}