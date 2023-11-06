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
#include <assert.h>

#include <igl/barycentric_coordinates.h>
#include <igl/point_mesh_squared_distance.h>

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


void rmt::WeightMap(const double* VSrc,
                    int NVerts,
                    const int* FSrc,
                    int NFaces,
                    const double* VTrg,
                    int NTVerts,
                    int** UI,
                    int** UJ,
                    double** UV,
                    size_t* NNZ)
{
    Eigen::MatrixXd V(Eigen::MatrixXd::Map(VSrc, NVerts, 3));
    Eigen::MatrixXd P(Eigen::MatrixXd::Map(VTrg, NTVerts, 3));
    Eigen::MatrixXi F(Eigen::MatrixXi::Map(FSrc, NFaces, 3));
    F = F.array() - 1;

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

    Eigen::MatrixXi F2 = F(I, Eigen::placeholders::all);
    Eigen::MatrixXd Va = V(F2.col(0), Eigen::placeholders::all);
    Eigen::MatrixXd Vb = V(F2.col(1), Eigen::placeholders::all);
    Eigen::MatrixXd Vc = V(F2.col(2), Eigen::placeholders::all);
    Eigen::MatrixXd L;
    igl::barycentric_coordinates(C, Va, Vb, Vc, L);

    *UI = (int*)malloc(3 * NTVerts * sizeof(int));
    *UJ = (int*)malloc(3 * NTVerts * sizeof(int));
    *UV = (double*)malloc(3 * NTVerts * sizeof(double));
    *NNZ = 3 * NTVerts;

    for (int i = 0; i < NTVerts; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            (*UI)[3 * i + j] = i + 1;
            (*UJ)[3 * i + j] = F2(i, j) + 1;
            (*UV)[3 * i + j] = L(i, j);
        }
    }

}