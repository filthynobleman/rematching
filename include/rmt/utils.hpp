/**
 * @file        utils.hpp
 * 
 * @brief       Defines some utilities for this project.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-12
 */
#pragma once

#include <utility>
#include <functional>


namespace rmt
{

template<typename T>
struct PairHash
{
    std::size_t operator()(const std::pair<T, T>& X) const noexcept
    {
        std::hash<T> h{};
        return h(X.first) ^ h(X.second);
    }
};

template<typename T>
struct TripleHash
{
    std::size_t operator()(const std::tuple<T, T, T>& X) const noexcept
    {
        std::hash<T> h{};
        return h(std::get<0>(X)) ^ (h(std::get<1>(X) ^ (h(std::get<2>(X)) << 1)) << 1);
    }
};




} // namespace rmt
