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
#include <rmt/preprocess.hpp>
#include <cut/cut.hpp>
#define NOMINMAX
#include <igl/unique_edge_map.h>
#include <igl/barycenter.h>
#include <igl/doublearea.h>

rmt::Mesh::Mesh(const Eigen::MatrixXd& V,
                const Eigen::MatrixXi& F)
    : m_V(V), m_F(F)
{
    m_E.resize(0, 2);
    m_BE.resize(0);
    m_BV.resize(0);
}

rmt::Mesh::Mesh(const std::string& Filename)
{
    CUTAssert(rmt::LoadMesh(Filename, m_V, m_F));
    m_E.resize(0, 2);
    m_BE.resize(0);
    m_BV.resize(0);
}

rmt::Mesh::Mesh(const rmt::Mesh& M)
    : m_V(M.m_V), m_F(M.m_F), m_E(M.m_E), m_BE(M.m_BE), m_BV(M.m_BV)
{ }

rmt::Mesh& rmt::Mesh::operator=(const rmt::Mesh& M)
{
    m_V = M.m_V;
    m_F = M.m_F;
    m_E = M.m_E;
    m_BE = M.m_BE;
    m_BV = M.m_BV;
    return *this;
}

rmt::Mesh::Mesh(rmt::Mesh&& M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_E = std::move(M.m_E);
    m_BE = std::move(M.m_BE);
    m_BV = std::move(M.m_BV);
}

rmt::Mesh& rmt::Mesh::operator=(rmt::Mesh&& M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_E = std::move(M.m_E);
    m_BE = std::move(M.m_BE);
    m_BV = std::move(M.m_BV);
    return *this;
}

rmt::Mesh::~Mesh() { }


int rmt::Mesh::NumVertices() const          { return m_V.rows(); }
int rmt::Mesh::NumEdges() const             { return m_E.rows(); }
int rmt::Mesh::NumTriangles() const         { return m_F.rows(); }
int rmt::Mesh::NumBoundaryEdges() const     { return m_BE.sum(); }
int rmt::Mesh::NumBoundaryVertices() const  { return m_BV.sum(); }


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
const Eigen::VectorXi& rmt::Mesh::GetBoundaryVertices() const
{
    CUTAssert(m_BV.size() > 0);
    return m_BV;
}


void rmt::Mesh::ComputeEdgesAndBoundaries()
{
    // Compute all the unique edges and the number of occurrences of each edge
    Eigen::MatrixXi AllE;
    Eigen::VectorXi EMap, uEC, uEE;
    igl::unique_edge_map(m_F, AllE, m_E, EMap, uEC, uEE);

    // Boundary edges are edges which occur only once
    m_BE = uEC(Eigen::seq(1, uEC.rows() - 1)) - uEC(Eigen::seq(0, uEC.rows() - 2));
    m_BV.setZero(NumVertices());
    for (int i = 0; i < m_BE.rows(); ++i)
    {
        m_BE[i] = m_BE[i] == 1 ? 1 : 0;
        if (m_BE[i])
        {
            m_BV[m_E(i, 0)] = 1;
            m_BV[m_E(i, 1)] = 1;
        }
    }
}


void rmt::Mesh::Scale(double Alpha)
{
    m_V *= Alpha;
}


void rmt::Mesh::Translate(const Eigen::VectorXd& Movement)
{
    Eigen::Vector3d Mov3 = Movement.segment<3>(0);
    Translate(Mov3);
}

void rmt::Mesh::Translate(const Eigen::Vector3d& Movement)
{
    m_V.rowwise() += Movement.transpose();
}



void rmt::Mesh::CenterAtOrigin()
{
    Eigen::MatrixXd BC;
    igl::barycenter(m_V, m_F, BC);
    Eigen::VectorXd Center = BC.colwise().mean();
    Translate(Center);
}


void rmt::Mesh::RescaleInsideUnitBox()
{
    double L = m_V.cwiseAbs().maxCoeff();
    m_V /= L;
}

void rmt::Mesh::RescaleInsideUnitSphere()
{
    double L = std::sqrt(m_V.rowwise().squaredNorm().maxCoeff());
    m_V /= L;
}


void rmt::Mesh::Resample(int OutputSize)
{
    double MEL = rmt::MaxEdgeLength(m_V, m_F, OutputSize);
    rmt::ResampleMesh(m_V, m_F, MEL);
}