/**
 * @file        rmt.hpp
 * 
 * @brief       ReMatching interface to MATLAB Mex.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#pragma once

#include <cstddef>

namespace rmt
{
    

void Remesh(const double* Vin,
            int NVerts,
            const int* Fin,
            int NFaces,
            int NSamples,
            double** Vout,
            size_t* NVout,
            int** Fout,
            size_t* NFout);

void Resample(const double* Vin,
              int NVerts,
              const int* Fin,
              int NFaces,
              int NSamples,
              double** Vout,
              size_t* NVout,
              int** Fout,
              size_t* NFout);

void WeightMap(const double* VSrc,
               int NVerts,
               const int* FSrc,
               int NFaces,
               const double* VTrg,
               int NTVerts,
               int** UI,
               int** UJ,
               double** UV,
               size_t* NNZ);


} // namespace rmt
