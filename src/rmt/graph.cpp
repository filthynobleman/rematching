/**
 * @file        graph.cpp
 * 
 * @brief       Impllements rmt::Graph.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-17
 */
#include <rmt/graph.hpp>
#include <cut/cut.hpp>
#include <set>
#include <queue>
#include <algorithm>


using namespace rmt;

typedef std::pair<int, int> Edge;  // Mesh edges


Graph::Graph(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    int nVerts = V.rows();
    // m_Verts.resize(nVerts);
    // for (int i = 0; i < nVerts; ++i)
    //     m_Verts[i] = V.row(i).segment<3>(0);

    // std::set<Edge> Edges;
    std::vector<Edge> Edges;
    Edges.reserve(6 * F.rows());
    int nTris = F.rows();
    for (int i = 0; i < nTris; ++i)
    {
        Edges.emplace_back(F(i, 0), F(i, 1));
        Edges.emplace_back(F(i, 1), F(i, 2));
        Edges.emplace_back(F(i, 2), F(i, 0));

        Edges.emplace_back(F(i, 1), F(i, 0));
        Edges.emplace_back(F(i, 2), F(i, 1));
        Edges.emplace_back(F(i, 0), F(i, 2));
    }
    std::sort(Edges.begin(), Edges.end());
    auto EEnd = std::unique(Edges.begin(), Edges.end());

    m_Idxs.resize(nVerts + 1);
    // m_Adjs.reserve(Edges.size());
    m_Adjs.reserve(std::distance(Edges.begin(), EEnd) + 1);
    int CurNode = 0;
    Eigen::Vector3d CurVert = V.row(0);// m_Verts[0];
    for (auto it = Edges.begin(); it != EEnd; it++)
    {
        if (it->first != CurNode)
        {
            CurNode++;
            m_Idxs[CurNode] = m_Adjs.size();
            CurVert = V.row(CurNode);// m_Verts[CurNode];
        }

        int ONode = it->second;
        Eigen::Vector3d OVert = V.row(ONode);// m_Verts[ONode];
        m_Adjs.emplace_back(ONode, (CurVert - OVert).norm());
    }
    m_Idxs[CurNode + 1] = m_Adjs.size();
}

Graph::Graph(const Eigen::MatrixXd& V, const std::vector<std::pair<int, int>>& E)
{
    int nVerts = V.rows();
    // m_Verts.resize(nVerts);
    // for (int i = 0; i < nVerts; ++i)
    //     m_Verts[i] = V.row(i).segment<3>(0);

    std::set<Edge> Edges;
    int nEdges = E.size();
    for (int i = 0; i < nEdges; ++i)
    {
        Edges.emplace(std::min(E[i].first, E[i].second), std::max(E[i].first, E[i].second));
        Edges.emplace(std::max(E[i].first, E[i].second), std::min(E[i].first, E[i].second));
    }

    m_Idxs.resize(nVerts + 1);
    m_Adjs.reserve(Edges.size());
    int CurNode = 0;
    Eigen::Vector3d CurVert = V.row(0);// m_Verts[0];
    for (auto it = Edges.begin(); it != Edges.end(); it++)
    {
        if (it->first != CurNode)
        {
            CurNode++;
            m_Idxs[CurNode] = m_Adjs.size();
            CurVert = V.row(CurNode);// m_Verts[CurNode];
        }

        int ONode = it->second;
        Eigen::Vector3d OVert = V.row(ONode);// m_Verts[ONode];
        m_Adjs.emplace_back(ONode, (CurVert - OVert).norm());
    }
    m_Idxs[CurNode + 1] = m_Adjs.size();
}

Graph::Graph(const Eigen::MatrixXd& V, const std::set<std::pair<int, int>>& E)
{
    int nVerts = V.rows();
    // m_Verts.resize(nVerts);
    // for (int i = 0; i < nVerts; ++i)
    //     m_Verts[i] = V.row(i).segment<3>(0);

    std::set<Edge> Edges;
    int nEdges = E.size();
    for (auto it = E.begin(); it != E.end(); it++)
    {
        Edges.emplace(std::min(it->first, it->second), std::max(it->first, it->second));
        Edges.emplace(std::max(it->first, it->second), std::min(it->first, it->second));
    }

    m_Idxs.resize(nVerts + 1);
    m_Adjs.reserve(Edges.size());
    int CurNode = 0;
    Eigen::Vector3d CurVert = V.row(0);// m_Verts[0];
    for (auto it = Edges.begin(); it != Edges.end(); it++)
    {
        if (it->first != CurNode)
        {
            CurNode++;
            m_Idxs[CurNode] = m_Adjs.size();
            CurVert = V.row(CurNode);// m_Verts[CurNode];
        }

        int ONode = it->second;
        Eigen::Vector3d OVert = V.row(ONode);// m_Verts[ONode];
        m_Adjs.emplace_back(ONode, (CurVert - OVert).norm());
    }
    m_Idxs[CurNode + 1] = m_Adjs.size();
}


Graph::Graph(const Graph& G)
{
    // m_Verts = G.m_Verts;
    m_Idxs = G.m_Idxs;
    m_Adjs = G.m_Adjs;
}

Graph& Graph::operator=(const Graph& G)
{
    // m_Verts = G.m_Verts;
    m_Idxs = G.m_Idxs;
    m_Adjs = G.m_Adjs;

    return *this;
}

Graph::Graph(Graph&& G)
{
    // m_Verts = std::move(G.m_Verts);
    m_Idxs = std::move(G.m_Idxs);
    m_Adjs = std::move(G.m_Adjs);
}

Graph& Graph::operator=(Graph&& G)
{
    // m_Verts = std::move(G.m_Verts);
    m_Idxs = std::move(G.m_Idxs);
    m_Adjs = std::move(G.m_Adjs);

    return *this;
}

Graph::~Graph() { }



// int Graph::NumVertices() const { return m_Verts.size(); }
int Graph::NumVertices() const { return m_Idxs.size() - 1; }
int Graph::NumEdges() const { return m_Adjs.size() / 2; }
int Graph::NumAdjacents(int i) const { return m_Idxs[i + 1] - m_Idxs[i]; }

// const Eigen::Vector3d& Graph::GetVertex(int i) const { return m_Verts[i]; }
// const std::vector<Eigen::Vector3d>& Graph::GetVertices() const { return m_Verts; }

const WEdge& Graph::GetAdjacent(int node_i, int adj_i) const
{
    return m_Adjs[m_Idxs[node_i] + adj_i];
}


rmt::Path Graph::DijkstraPath(int src, int dst) const
{
    Eigen::VectorXd Dists;
    Dists.setConstant(NumVertices(), std::numeric_limits<double>::infinity());

    Eigen::VectorXi Parent;
    Parent.setConstant(NumVertices(), -1);

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i;
        double wi;
        std::tie(wi, i) = Q.top();
        Q.pop();

        if (i == dst)
            break;
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Parent[j] = i;
            Q.emplace(Dists[j], j);
        }
    }

    std::vector<rmt::WEdge> Path;
    double Length = Dists[dst];
    while (Parent[dst] != -1)
    {
        Path.emplace_back(dst, Dists[dst] - Dists[Parent[dst]]);
        dst = Parent[dst];
    }
    Path.emplace_back(dst, Dists[dst]);
    std::reverse(Path.begin(), Path.end());
    CUTAssert(Path[0].first == src);
    return { Length, Path };
}

Eigen::VectorXd rmt::Graph::DijkstraDistance(int src) const
{
    Eigen::VectorXd D;
    DijkstraDistance(src, D);
    return D;
}

void rmt::Graph::DijkstraDistance(int src, Eigen::VectorXd& Dists) const
{
    Dists.setConstant(NumVertices(), std::numeric_limits<double>::infinity());

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i;
        double wi;
        std::tie(wi, i) = Q.top();
        Q.pop();
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Q.emplace(Dists[j], j);
        }
    }
}


int rmt::Graph::FarthestFiltered(int src, const std::vector<int>& Tag, int Filter) const
{
    Eigen::VectorXd Dists;
    Dists.setConstant(NumVertices(), std::numeric_limits<double>::infinity());

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    int Farthest = src;
    double MaxDist = 0.0;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i;
        double wi;
        std::tie(wi, i) = Q.top();
        Q.pop();
        if (wi > MaxDist)
        {
            MaxDist = wi;
            Farthest = i;
        }
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Tag[j] != Filter)
                continue;
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Q.emplace(Dists[j], j);
        }
    }

    return Farthest;
}

int rmt::Graph::FarthestAtBoundary(int src, const std::vector<int>& Tag, int Region, int Neighbor) const
{
    Eigen::VectorXd Dists;
    Dists.setConstant(NumVertices(), std::numeric_limits<double>::infinity());

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    int Farthest = src;
    double MaxDist = 0.0;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i;
        double wi;
        std::tie(wi, i) = Q.top();
        Q.pop();
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Tag[j] == Neighbor && wi > MaxDist)
            {
                MaxDist = wi;
                Farthest = i;
            }
            if (Tag[j] != Region)
                continue;
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Q.emplace(Dists[j], j);
        }
    }

    return Farthest;
}


std::vector<int> Graph::ConnectedComponents() const
{
    std::vector<int> CC;
    CC.resize(NumVertices(), 0);
    int CurCC = 0;

    Eigen::VectorXi Visited;
    Visited.resize(NumVertices());
    Visited.setZero();

    while (Visited.sum() < NumVertices())
    {
        int Root = 0;
        for (int i = 0; i < NumVertices(); ++i)
        {
            if (Visited[i] == 0)
            {
                Root = i;
                break;
            }
        }

        std::queue<int> Q;
        Q.emplace(Root);
        while (!Q.empty())
        {
            int N = Q.front();
            Q.pop();
            
            if (Visited[N] != 0)
                continue;

            CC[N] = CurCC;
            Visited[N] = 1;
            int DegN = NumAdjacents(N);
            for (int jj = 0; jj < DegN; ++jj)
            {
                int j = GetAdjacent(N, jj).first;
                Q.emplace(j);
            }
        }

        CurCC += 1;
    }

    return CC;
}