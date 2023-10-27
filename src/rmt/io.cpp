/**
 * @file        io.cpp
 * 
 * @brief       Implements I/O functionalities.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#include <rmt/io.hpp>

#include <filesystem>
#include <algorithm>
#include <cctype>

#include <unsupported/Eigen/SparseExtra>

#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>

#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writePLY.h>

bool rmt::ExportWeightmap(const std::string & Filename, 
                          const Eigen::SparseMatrix<double>& WM)
{
    return Eigen::saveMarket(WM, Filename);
}


bool rmt::LoadMesh(const std::string& Filename,
                   Eigen::MatrixXd& V,
                   Eigen::MatrixXi& F)
{
    std::string Ext;
    Ext = std::filesystem::path(Filename).extension().string();
    std::transform(Ext.begin(), Ext.end(), Ext.begin(), [](int c) { return std::tolower(c); });
    
    if (Ext == ".obj")
        return igl::readOBJ(Filename, V, F);
    else if (Ext == ".off")
        return igl::readOFF(Filename, V, F);
    else if (Ext == ".ply")
        return igl::readPLY(Filename, V, F);

    return false;
}


bool rmt::ExportMesh(const std::string & Filename, 
                     const Eigen::MatrixXd & V, 
                     const Eigen::MatrixXi & F)
{
    std::string Ext;
    Ext = std::filesystem::path(Filename).extension().string();
    std::transform(Ext.begin(), Ext.end(), Ext.begin(), [](int c) { return std::tolower(c); });

    if (Ext == ".obj")
        return igl::writeOBJ(Filename, V, F);
    else if (Ext == ".off")
        return igl::writeOFF(Filename, V, F);
    else if (Ext == ".ply")
        return igl::writePLY(Filename, V, F);
    
    return false;
}