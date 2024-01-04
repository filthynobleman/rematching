/**
 * @file        voronoifps.cpp
 * 
 * @brief       Implements VoronoiFPS().
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-17
 */
#include <rmt/voronoifps.hpp>
#include <random>
#include <queue>
#include <igl/exact_geodesic.h>
#include <cut/algo/minheap.hpp>

using namespace rmt;


std::vector<double> rmt::Distances(const Graph& G, int N)
{
    std::vector<double> Dists;
    Dists.resize(G.NumVertices(), std::numeric_limits<double>::infinity());
    Dists[N] = 0;

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    Q.emplace(0, N);
    while (!Q.empty())
    {
        std::pair<double, int> Next = Q.top();
        Q.pop();

        int Cur = Next.second;
        double W = Next.first;
        int Deg = G.NumAdjacents(Cur);
        for (int j = 0; j < Deg; ++j)
        {
            WEdge Neig = G.GetAdjacent(Cur, j);
            if (Dists[Neig.first] <= W + Neig.second)
                continue;
            Dists[Neig.first] = W + Neig.second;
            Q.emplace(Dists[Neig.first], Neig.first);
        }
    }

    return Dists;
}


rmt::VoronoiPartitioning rmt::VoronoiFPS(const Graph& G, 
                                          int NSamples, 
                                          int Seed)
{
    if (NSamples >= G.NumVertices())
    {
        NSamples = G.NumVertices();
        std::vector<int> Partition;
        std::vector<int> Samples;
        Partition.resize(NSamples);
        Samples.resize(NSamples);
        for (int i = 0; i < NSamples; ++i)
        {
            Samples[i] = i;
            Partition[i] = i;
        }
        return { Samples, Partition };
    }

    std::mt19937 Eng(Seed);
    std::uniform_int_distribution<int> Distr(0, G.NumVertices() - 1);

    int NVerts = G.NumVertices();
    
    int p = Distr(Eng);
    std::vector<double> Dists = Distances(G, p);
    std::vector<int> Partition;
    Partition.resize(NVerts, p);
    std::vector<int> Samples;
    Samples.resize(NSamples, p);

    cut::MinHeap HDists(Dists, true);

    for (int h = 1; h < NSamples; ++h)
    {
        // p = argmax_j Dists[j]
        p = HDists.FindMin().second;
        HDists.SetKey(p, 0);
        Samples[h] = p;

        rmt::UpdateVoronoi(G, Dists, Partition, HDists, p);
    }

    return { Samples, Partition, Dists };
}

void rmt::UpdateVoronoi(const rmt::Graph& G,
                        std::vector<double>& Dists,
                        std::vector<int>& Partition,
                        cut::MinHeap& HDists,
                        int p)
{
    std::priority_queue<std::pair<double, int>,
                            std::vector<std::pair<double, int>>,
                            std::greater<std::pair<double, int>>> Q;
    Dists[p]= 0;
    Partition[p] = p;
    Q.emplace(0, p);
    while (!Q.empty())
    {
        std::pair<double, int> Next = Q.top();
        Q.pop();

        int Cur = Next.second;
        double W = Next.first;

        if (W < HDists.GetKey(Cur))
            HDists.SetKey(Cur, W);

        int Deg = G.NumAdjacents(Cur);
        for (int j = 0; j < Deg; ++j)
        {
            WEdge Neig = G.GetAdjacent(Cur, j);
            if (Dists[Neig.first] <= W + Neig.second)
                continue;
            Dists[Neig.first] = W + Neig.second;
            Q.emplace(Dists[Neig.first], Neig.first);
            Partition[Neig.first] = p;
        }
    }
}