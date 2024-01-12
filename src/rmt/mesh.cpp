/**
 * @file        mesh.cpp
 * 
 * @brief       Implementation of class rmt::Mesh.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#include <rmt/mesh.hpp>
#include <rmt/io.hpp>
#include <cut/cut.hpp>
#define NOMINMAX
#include <igl/unique_edge_map.h>

rmt::Mesh::Mesh(const Eigen::MatrixXd& V,
                const Eigen::MatrixXi& F)
    : m_V(V), m_F(F)
{
    m_E.resize(0);
    m_BE.resize(0);
}

rmt::Mesh::Mesh(const std::string& Filename)
{
    CUTAssert(rmt::LoadMesh(Filename, m_V, m_F));
    m_E.resize(0);
    m_BE.resize(0);
}

rmt::Mesh::Mesh(const rmt::Mesh& M)
    : m_V(M.m_V), m_F(M.m_F), m_E(M.m_E), m_BE(M.m_BE)
{ }

rmt::Mesh& rmt::Mesh::operator=(const rmt::Mesh& M)
{
    m_V = M.m_V;
    m_F = M.m_F;
    m_E = M.m_E;
    m_BE = M.m_BE;
    return *this;
}

rmt::Mesh::Mesh(rmt::Mesh&& M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_E = std::move(M.m_E);
    m_BE = std::move(M.m_BE);
}

rmt::Mesh& rmt::Mesh::operator=(rmt::Mesh&& M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_E = std::move(M.m_E);
    m_BE = std::move(M.m_BE);
    return *this;
}

rmt::Mesh::~Mesh() { }


const Eigen::MatrixXd& rmt::Mesh::GetVertices() const       { return m_V; }
const Eigen::MatrixXi& rmt::Mesh::GetTriangles() const      { return m_F; }
const Eigen::MatrixXi& rmt::Mesh::GetEdges() const
{ 
    CUTAssert(m_E.size() > 0);
    return m_E;
}
const Eigen::VectorXi& rmt::Mesh::GetBoundaryEdges() const
{
    CUTAssert(m_BE.size() > 0);
    return m_BE;
}


void rmt::Mesh::ComputeEdges()
{
    // Compute all the unique edges and the number of occurrences of each edge
    Eigen::MatrixXi AllE;
    Eigen::VectorXi EMap, uEC, uEE;
    igl::unique_edge_map(m_F, AllE, m_E, EMap, uEC, uEE);

    // Boundary edges are edges which occur only once
    m_BE = uEC(Eigen::seq(1, uEC.rows() - 1)) - uEC(Eigen::seq(0, uEC.rows() - 2));
    for (int i = 0; i < m_BE.rows(); ++i)
        m_BE[i] = m_BE[i] == 1 ? 1 : 0;
}