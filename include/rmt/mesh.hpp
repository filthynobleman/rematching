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

    // Identifies boundary vertices #V x 1
    Eigen::VectorXi m_BV;


public:
    Mesh(const Eigen::MatrixXd& V,
         const Eigen::MatrixXi& F);
    Mesh(const std::string& Filename);
    Mesh(const rmt::Mesh& M);
    Mesh(rmt::Mesh&& M);
    rmt::Mesh& operator=(const rmt::Mesh& M);
    rmt::Mesh& operator=(rmt::Mesh&& M);
    ~Mesh();

    int NumVertices() const;
    int NumBoundaryVertices() const;
    int NumEdges() const;
    int NumBoundaryEdges() const;
    int NumTriangles() const;

    const Eigen::MatrixXd& GetVertices() const;
    const Eigen::MatrixXi& GetTriangles() const;
    const Eigen::MatrixXi& GetEdges() const;
    const Eigen::VectorXi& GetBoundaryEdges() const;
    const Eigen::VectorXi& GetBoundaryVertices() const;

    void ComputeEdgesAndBoundaries();

    void Scale(double Alpha);
    void Translate(const Eigen::Vector3d& Movement);
    void Translate(const Eigen::VectorXd& Movement);

    void CenterAtOrigin();
    void RescaleInsideUnitBox();
    void RescaleInsideUnitSphere();

    void Resample(int OutputSize);
};


} // namespace rmt
