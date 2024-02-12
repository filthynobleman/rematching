/**
 * @file        clean.cpp
 * 
 * @brief       Implements cleaning operations.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-30
 */
#include <rmt/clean.hpp>
#include <rmt/utils.hpp>

#include <unordered_map>
#include <set>
#include <queue>

#define NOMINMAX
#include <igl/doublearea.h>
#include <igl/facet_components.h>
namespace Eigen
{
    Eigen::internal::all_t all = Eigen::placeholders::all;
};
#include <igl/remove_unreferenced.h>
#include <igl/remove_duplicate_vertices.h>


typedef std::pair<int, int> IntPair;
typedef std::tuple<int, int, int> IntTriple;
typedef std::vector<int> IntVec;

typedef std::unordered_map<IntPair, IntVec, rmt::PairHash<int>> E2T_t;
typedef std::vector<IntVec> V2T_t;
typedef std::vector<IntTriple> T2T_t;


void DeleteTriangles(Eigen::MatrixXi& F,
                     const std::set<int>& TriDelete)
{
    int EndPtr = F.rows() - 1;
    for (int t : TriDelete)
    {
        // EndPtr is the last manifold triangle of the mesh
        while (TriDelete.find(EndPtr) != TriDelete.end())
            EndPtr -= 1;
        
        // If triangle is already in the to-remove region, good
        if (t > EndPtr)
            continue;
        
        // Swap the triangle to remove with the last manifold triangle
        F.row(t) = F.row(EndPtr--);
    }

    // Resize F to remove the non-manifold triangles
    F.conservativeResize(EndPtr + 1, 3);
}

void DeleteVertices(Eigen::MatrixXd& V,
                    Eigen::MatrixXi& F,
                    const std::set<int>& VertDelete)
{
    // Mapping from old to new vertex indices
    Eigen::VectorXi VMap = Eigen::VectorXi::LinSpaced(V.rows(), 0, V.rows() - 1);
    for (int v : VertDelete)
        VMap[v] = -1;


    int EndPtr = V.rows() - 1;
    for (int v : VertDelete)
    {
        // EndPtr is the last manifold vertex of the mesh
        while (VertDelete.find(EndPtr) != VertDelete.end())
            EndPtr -= 1;
        
        // If vertex is already in the to-remove region, good
        if (v > EndPtr)
            continue;
        
        // Swap the vertex to remove with the last manifold vertex
        // Change the indices, so that the vertex that previously was at EndPtr
        // goes to v, and v goes to -1 (so that we can also remove faces pointing
        // to deleted vertices)
        VMap[EndPtr] = v;
        VMap[v] = -1;
        // Swap the vertex content
        V.row(v) = V.row(EndPtr--);
    }

    // Resize V to remove the non-manifold vertices
    V.conservativeResize(EndPtr + 1, 3);

    // Remap all the indices in F
    std::set<int> TriDelete;
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            F(i, j) = VMap[F(i,j)];
            if (F(i, j) == -1)
            {
                TriDelete.emplace(i);
                break;
            }
        }
    }

    // Eventually, remove triangles referring to deleted vertices
    DeleteTriangles(F, TriDelete);
}




void RepairNonManifoldEdges(Eigen::MatrixXd& V, 
                            Eigen::MatrixXi& F)
{
    bool NonManifold = false;

    // Map edges to triangles
    E2T_t E2T;
    E2T.reserve((3 * F.rows()) / 2);
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j0 = 0; j0 < 3; ++j0)
        {
            int j1 = j0 + 1;
            if (j1 == 3)
                j1 = 0;
            
            std::pair<int, int> e{F(i, j0), F(i, j1)};
            if (e.second < e.first)
                std::swap(e.first, e.second);
            
            if (E2T.find(e) == E2T.end())
            {
                E2T.emplace(e, std::vector<int>{});
                E2T[e].reserve(3);
            }
            E2T[e].emplace_back(i);

            if (E2T[e].size() > 2)
                NonManifold = true;
        }
    }

    if (!NonManifold)
        return;

    // Compute triangle areas
    Eigen::VectorXd Areas;
    igl::doublearea(V, F, Areas);

    // For each non manifold edge, keep only the two triangles with larger areas
    std::vector<std::pair<double, int>> area_loc;
    area_loc.reserve(10);
    std::set<int> TriDelete;
    for (auto p : E2T)
    {
        const std::pair<int, int>& e = p.first;
        const std::vector<int>& tris = p.second;
        if (tris.size() < 2)
            continue;
        area_loc.clear();
        area_loc.reserve(tris.size());
        for (int t : tris)
            area_loc.emplace_back(Areas[t], t);
        std::sort(area_loc.begin(), area_loc.end(), std::greater<>{});
        for (int i = 2; i < area_loc.size(); ++i)
            TriDelete.emplace(area_loc[i].second);
    }

    // Effectively remove the triangles from the list
    DeleteTriangles(F, TriDelete);
}


void RepairNonManifoldVertices(Eigen::MatrixXd& V,
                               Eigen::MatrixXi& F)
{
    // Map edges to triangles, and we assume edges are all manifold
    // We also map vertices to triangles
    E2T_t E2T;
    E2T.reserve((3 * F.rows()) / 2);
    std::vector<std::vector<int>> V2T;
    V2T.resize(V.rows());
    for (auto& TList : V2T)
        TList.reserve(10);
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j0 = 0; j0 < 3; ++j0)
        {
            V2T[F(i, j0)].emplace_back(i);

            int j1 = j0 + 1;
            if (j1 == 3)
                j1 = 0;
            
            std::pair<int, int> e{F(i, j0), F(i, j1)};
            if (e.second < e.first)
                std::swap(e.first, e.second);
            
            if (E2T.find(e) == E2T.end())
            {
                IntVec ts = { i, -1 };
                E2T.emplace(e, ts);
            }
            else
                E2T[e][1] = i;
        }
    }

    // Build triangle-triangle adjacency
    T2T_t T2T;
    T2T.resize(F.rows(), { -1, -1, -1 });
    for (auto p : E2T)
    {
        const std::pair<int, int>& e = p.first;
        const auto& tris = p.second;

        // Boundary edges only have one incident triangle
        if (tris[1] == -1)
            continue;

        std::tuple<int, int, int>& t0 = T2T[tris[0]];
        if (std::get<0>(t0) == -1)
            std::get<0>(t0) = tris[1];
        else if (std::get<1>(t0) == -1)
            std::get<1>(t0) = tris[1];
        else
            std::get<2>(t0) = tris[1];
        
        std::tuple<int, int, int>& t1 = T2T[tris[1]];
        if (std::get<0>(t1) == -1)
            std::get<0>(t1) = tris[0];
        else if (std::get<1>(t1) == -1)
            std::get<1>(t1) = tris[0];
        else
            std::get<2>(t1) = tris[0];
    }

    // Check if vertex manifold
    std::map<int, std::vector<int>> NewVs;
    int LastV = V.rows();
    for (int v = 0; v < V.rows(); ++v)
    {
        // A vertex is manifold if all the triangles incident on it form a 
        // connected component via edge connection
        std::set<int> Visited;
        Visited.clear();
        for (int t : V2T[v])
        {
            if (Visited.find(t) != Visited.end())
                continue;

            std::set<int> VisitedFan;
            std::queue<int> Q;
            Q.push(t);
            while (!Q.empty())
            {
                int ct = Q.front();
                Q.pop();
                if (VisitedFan.find(ct) != VisitedFan.end())
                    continue;
                VisitedFan.emplace(ct);
                
                int tadj;
                tadj = std::get<0>(T2T[ct]);
                if (std::find(V2T[v].begin(), V2T[v].end(), tadj) != V2T[v].end())
                    Q.push(tadj);
                tadj = std::get<1>(T2T[ct]);
                if (std::find(V2T[v].begin(), V2T[v].end(), tadj) != V2T[v].end())
                    Q.push(tadj);
                tadj = std::get<2>(T2T[ct]);
                if (std::find(V2T[v].begin(), V2T[v].end(), tadj) != V2T[v].end())
                    Q.push(tadj);
            }

            // Update the visited triangles
            Visited.insert(VisitedFan.begin(), VisitedFan.end());

            // If all triangles have been visited, the vertex is manifold or has become manifold
            // We don't have to chrck further
            if (Visited.size() == V2T[v].size())
                break;

            // Otherwise, create a new vertex and disconnect the visited fan
            if (NewVs.find(v) == NewVs.end())
                NewVs.emplace(v, std::vector<int>{});
            for (int t : VisitedFan)
            {
                for (int j = 0; j < 3; ++j)
                {
                    if (F(t, j) == v)
                        F(t, j) = LastV;
                }
            }
            NewVs[v].emplace_back(LastV++);
        }
    }

    // If we created some new vertex, we duplicate the coordinates
    if (NewVs.size() == 0)
        return;

    V.conservativeResize(LastV, 3);
    for (auto p : NewVs)
    {
        for (int v : p.second)
            V.row(v) = V.row(p.first);
    }
}



void rmt::MakeManifold(Eigen::MatrixXd& V,
                       Eigen::MatrixXi& F)
{
    RepairNonManifoldEdges(V, F);
    RepairNonManifoldVertices(V, F);
}


void rmt::RemoveSmallComponents(Eigen::MatrixXd& V,
                                Eigen::MatrixXi& F,
                                double AreaFraction)
{
    Eigen::VectorXd dblA;
    igl::doublearea(V, F, dblA);
    double TotArea = dblA.sum();
    Eigen::VectorXi C;
    std::vector<int> CCount;
    std::vector<double> CArea;
    CArea.resize(igl::facet_components(F, C));
    CCount.resize(CArea.size());
    for (int i = 0; i < C.rows(); ++i)
    {
        CArea[C[i]] += dblA[i];
        CCount[C[i]] += 1;
    }
    std::set<int> CDelete;
    for (int i = 0; i < CCount.size(); ++i)
    {
        if (CArea[i] <= AreaFraction * TotArea || CCount[i] <= 3)
            CDelete.emplace(i);
    }
    std::set<int> TriDelete;
    for (int i = 0; i < C.rows(); ++i)
    {
        if (CDelete.find(C[i]) != CDelete.end())
            TriDelete.emplace(i);
    }
    std::set<int> VertDelete;
    for (int t : TriDelete)
    {
        for (int j = 0; j < 3; ++j)
            VertDelete.emplace(F(t, j));
    }
    DeleteTriangles(F, TriDelete);
    Eigen::MatrixXd VTmp = V;
    Eigen::MatrixXi FTmp = F;
    Eigen::VectorXi I;
    igl::remove_unreferenced(VTmp, FTmp, V, F, I);
}


void rmt::RemoveDegeneracies(Eigen::MatrixXd& V,
                             Eigen::MatrixXi& F,
                             double DistThreshold)
{
    Eigen::MatrixXd VTmp;
    Eigen::MatrixXi FTmp;
    Eigen::VectorXi I, J;
    igl::remove_unreferenced(V, F, VTmp, FTmp, I);
    igl::remove_duplicate_vertices(VTmp, FTmp, 1e-4, V, I, J, F);
    FTmp.resize(F.rows(), 3);
    int LastFace = 0;
    for (int i = 0; i < F.rows(); ++i)
    {
        if (F(i, 0) == F(i, 1) || F(i, 1) == F(i, 2) || F(i, 2) == F(i, 0))
            continue;
        FTmp.row(LastFace++) = F.row(i);
    }
    F = FTmp.block(0, 0, LastFace, 3);
}


void rmt::CleanUp(Eigen::MatrixXd& V,
                  Eigen::MatrixXi& F,
                  double AreaFraction,
                  double DistanceThreshold)
{
    RemoveDegeneracies(V, F, DistanceThreshold);
    MakeManifold(V, F);
    RemoveSmallComponents(V, F, AreaFraction);
}