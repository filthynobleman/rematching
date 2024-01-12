/**
 * @file        region.cpp
 * 
 * @brief       Implementation of class rmt::SurfaceRegion.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#include <rmt/region.hpp>


rmt::SurfaceRegion::SurfaceRegion(int pi, int pj, int pk)
{
    m_NVerts = 0;
    m_NEdges = 0;
    m_NFaces = 0;

    rmt::RegionDictionary::OrderIndices(pi, pj, pk);
    m_Samples = { pi, pj, pk };
}

rmt::SurfaceRegion::SurfaceRegion(const rmt::SurfaceRegion& SR)
{
    m_NVerts = SR.m_NVerts;
    m_NEdges = SR.m_NEdges;
    m_NFaces = SR.m_NFaces;
    m_Samples = SR.m_Samples;
}

rmt::SurfaceRegion& rmt::SurfaceRegion::operator=(const rmt::SurfaceRegion& SR)
{
    m_NVerts = SR.m_NVerts;
    m_NEdges = SR.m_NEdges;
    m_NFaces = SR.m_NFaces;
    m_Samples = SR.m_Samples;

    return *this;
}

rmt::SurfaceRegion::~SurfaceRegion() { }


std::tuple<int, int, int> rmt::SurfaceRegion::GetSamples() const { return m_Samples; }

int rmt::SurfaceRegion::NumVertices() const { return m_NVerts; }
int rmt::SurfaceRegion::NumFaces() const    { return m_NFaces; }
int rmt::SurfaceRegion::NumEdges() const    { return m_NEdges >> 1; } // divide by two
int rmt::SurfaceRegion::EulerCharacteristic() const
{
    return NumVertices() + NumFaces() - NumEdges();
}

void rmt::SurfaceRegion::AddFace()          { m_NFaces += 1; }
void rmt::SurfaceRegion::AddVertex()        { m_NVerts += 1; }
void rmt::SurfaceRegion::AddEdge()          { m_NEdges += 1; }
void rmt::SurfaceRegion::AddBoundaryEdge()  { m_NEdges += 2; }


bool rmt::operator==(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    return SR1.m_Samples == SR2.m_Samples;
}
bool rmt::operator!=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    return SR1.m_Samples != SR2.m_Samples;
}
bool rmt::operator<=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    return SR1.m_Samples <= SR2.m_Samples;
}
bool rmt::operator<(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    return SR1.m_Samples < SR2.m_Samples;
}
bool rmt::operator>=(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    return SR1.m_Samples >= SR2.m_Samples;
}
bool rmt::operator>(const rmt::SurfaceRegion& SR1, const rmt::SurfaceRegion& SR2)
{
    return SR1.m_Samples > SR2.m_Samples;
}