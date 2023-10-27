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


} // namespace rmt
