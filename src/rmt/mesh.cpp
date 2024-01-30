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
#include <rmt/clean.hpp>
#include <rmt/utils.hpp>
#include <cut/cut.hpp>
#include <cassert>
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
    // Compute all the edges and their number
    std::unordered_map<std::pair<int, int>, std::pair<int, int>, rmt::PairHash<int>> m_Ecount;
    m_Ecount.reserve(2 * NumTriangles());
    for (int i = 0; i < m_F.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int j1 = j + 1;
            if (j1 > 2)
                j1 = 0;
            std::pair<int, int> e{m_F(i, j), m_F(i, j1)};
            if (e.first > e.second)
                std::swap(e.first, e.second);
            if (m_Ecount.find(e) == m_Ecount.end())
            {
                std::pair<int, int> cnt{m_Ecount.size(), 0};
                m_Ecount.emplace(e, cnt);
            }
            m_Ecount[e].second += 1;
        }
    }

    m_E.resize(m_Ecount.size(), 2);
    m_BE.setZero(m_Ecount.size());
    m_BV.setZero(NumVertices());
    for (auto ec : m_Ecount)
    {
        m_E(ec.second.first, 0) = ec.first.first;
        m_E(ec.second.first, 1) = ec.first.second;
        if (ec.second.second == 1)
        {
            m_BE[ec.second.first] = 1;
            m_BV[m_E(ec.second.first, 0)] = 1;
            m_BV[m_E(ec.second.first, 1)] = 1;
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


void rmt::Mesh::MakeManifold()
{
    rmt::MakeManifold(m_V, m_F);
}

void rmt::Mesh::RemoveSmallComponents(double AreaFraction)
{
    rmt::RemoveSmallComponents(m_V, m_F, AreaFraction);
}

void rmt::Mesh::RemoveDegenaracies(double DistanceThreshold)
{
    rmt::RemoveDegeneracies(m_V, m_F, DistanceThreshold);
}

void rmt::Mesh::CleanUp(double AreaFraction, double DistanceThreshold)
{
    rmt::CleanUp(m_V, m_F, AreaFraction, DistanceThreshold);
}