/**
 * @file        mesh.hpp
 * 
 * @brief       A data structure representing a triangular mesh.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#pragma once

#include <Eigen/Dense>
#include <string>

namespace rmt
{
    
class Mesh
{
private:
    // Mesh vertices #V x 3
    Eigen::MatrixXd m_V;

    // Mesh triangles #F x 3
    Eigen::MatrixXi m_F;

    // Mesh (unique) edges #E x 2
    Eigen::MatrixXi m_E;

    // Identifies boundary (unique) edges #E x 1
    Eigen::VectorXi m_BE;


public:
    Mesh(const Eigen::MatrixXd& V,
         const Eigen::MatrixXi& F);
    Mesh(const std::string& Filename);
    Mesh(const rmt::Mesh& M);
    Mesh(rmt::Mesh&& M);
    rmt::Mesh& operator=(const rmt::Mesh& M);
    rmt::Mesh& operator=(rmt::Mesh&& M);
    ~Mesh();

    const Eigen::MatrixXd& GetVertices() const;
    const Eigen::MatrixXi& GetTriangles() const;
    const Eigen::MatrixXi& GetEdges() const;
    const Eigen::VectorXi& GetBoundaryEdges() const;

    void ComputeEdges();
};


} // namespace rmt
