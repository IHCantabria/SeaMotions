/*
 * Copyright (c) 2025 Sergio Fern\u00e1ndez Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "../config.hpp"
#include "../math/gauss_t.hpp"

#if __has_include("../math/cusvector.hpp")
#include "../math/cusvector.hpp"
#else
template<typename T>
using CusVector = std::vector<T>;
#endif

/**
 * @brief Angular quadrature rule for annulus sampling.
 */
enum class AngularRule
{
    Trapezoidal,
    GaussChebyshev
};

/**
 * @brief Node container for annulus integration sampling.
 *
 * Stores polar and Cartesian coordinates plus separable and combined weights
 * for numerical integration over an annulus [r0, r1] x [0, 2*pi].
 */
struct AnnulusNodeCloud
{
    int nodes_np = 0;
    int radial_np = 0;
    int theta_np = 0;

    CusVector<cusfloat> r;
    CusVector<cusfloat> theta;
    CusVector<cusfloat> x;
    CusVector<cusfloat> y;
    CusVector<cusfloat> z;
    CusVector<cusfloat> w_r;
    CusVector<cusfloat> w_theta;
    CusVector<cusfloat> w;
};

/**
 * @brief Build annulus sampling nodes for quadrature.
 *
 * Radial points use Gauss-Legendre of order @p RadialOrder.
 * Angular points use either a trapezoidal rule or Gauss-Chebyshev of the
 * first kind mapped to [0, 2*pi]. The final weight is:
 *
 *   w = w_r[i] * w_theta[j] * r_i
 *
 * which includes the Jacobian for polar coordinates.
 *
 * @tparam RadialOrder Gauss-Legendre order for the radial coordinate (1..19).
 * @param r0 Inner radius.
 * @param r1 Outer radius.
 * @param theta_np Number of angular points.
 * @param rule Angular quadrature rule.
 * @return AnnulusNodeCloud with nodes and weights.
 */
template<std::size_t RadialOrder>
inline AnnulusNodeCloud create_annulus_nodes(
    cusfloat r0,
    cusfloat r1,
    int theta_np,
    AngularRule rule = AngularRule::Trapezoidal)
{
    if (r1 <= r0)
    {
        throw std::runtime_error("create_annulus_nodes: r1 must be greater than r0.");
    }
    if (theta_np <= 0)
    {
        throw std::runtime_error("create_annulus_nodes: theta_np must be positive.");
    }

    constexpr cusfloat kPi = static_cast<cusfloat>(3.14159265358979323846);

    AnnulusNodeCloud cloud;
    cloud.radial_np = static_cast<int>(GaussPointsT<1, RadialOrder>::N);
    cloud.theta_np = theta_np;
    cloud.nodes_np = cloud.radial_np * cloud.theta_np;

    cloud.r.reserve(static_cast<size_t>(cloud.radial_np));
    cloud.theta.reserve(static_cast<size_t>(cloud.theta_np));
    cloud.w_r.reserve(static_cast<size_t>(cloud.radial_np));
    cloud.w_theta.reserve(static_cast<size_t>(cloud.theta_np));
    cloud.x.reserve(static_cast<size_t>(cloud.nodes_np));
    cloud.y.reserve(static_cast<size_t>(cloud.nodes_np));
    cloud.z.reserve(static_cast<size_t>(cloud.nodes_np));
    cloud.w.reserve(static_cast<size_t>(cloud.nodes_np));

    const cusfloat half = static_cast<cusfloat>(0.5) * (r1 - r0);
    const cusfloat mid = static_cast<cusfloat>(0.5) * (r1 + r0);

    for (int i = 0; i < cloud.radial_np; ++i)
    {
        const cusfloat xi = GaussPointsT<1, RadialOrder>::roots_x[i];
        const cusfloat wi = GaussPointsT<1, RadialOrder>::weights_x[i];
        const cusfloat ri = half * xi + mid;
        const cusfloat wri = half * wi;
        cloud.r.push_back(ri);
        cloud.w_r.push_back(wri);
    }

    if (rule == AngularRule::Trapezoidal)
    {
        const cusfloat wtheta = static_cast<cusfloat>(2.0) * kPi / static_cast<cusfloat>(theta_np);
        for (int j = 0; j < theta_np; ++j)
        {
            const cusfloat theta = static_cast<cusfloat>(2.0) * kPi * static_cast<cusfloat>(j) / static_cast<cusfloat>(theta_np);
            cloud.theta.push_back(theta);
            cloud.w_theta.push_back(wtheta);
        }
    }
    else
    {
        // Gauss-Chebyshev (first kind) mapped to [0, 2*pi] with weights adjusted to
        // approximate integral of f(theta) over theta without Chebyshev weight.
        // theta = pi * (1 + x), dtheta = pi dx, x in [-1, 1].
        // Integral f(theta) dtheta = pi * integral_{-1}^1 f(pi(1+x)) dx.
        // We approximate integral_{-1}^1 g(x) dx by Gauss-Chebyshev on g(x)*sqrt(1-x^2).
        const cusfloat base_w = kPi * kPi / static_cast<cusfloat>(theta_np);
        for (int j = 0; j < theta_np; ++j)
        {
            const cusfloat x = std::cos((static_cast<cusfloat>(2 * j + 1) * kPi) / (static_cast<cusfloat>(2 * theta_np)));
            const cusfloat theta = kPi * (static_cast<cusfloat>(1.0) + x);
            const cusfloat wtheta = base_w * std::sqrt(std::max(static_cast<cusfloat>(0.0), static_cast<cusfloat>(1.0) - x * x));
            cloud.theta.push_back(theta);
            cloud.w_theta.push_back(wtheta);
        }
    }

    for (int i = 0; i < cloud.radial_np; ++i)
    {
        const cusfloat ri = cloud.r[static_cast<size_t>(i)];
        const cusfloat wri = cloud.w_r[static_cast<size_t>(i)];
        for (int j = 0; j < cloud.theta_np; ++j)
        {
            const cusfloat theta = cloud.theta[static_cast<size_t>(j)];
            const cusfloat wtheta = cloud.w_theta[static_cast<size_t>(j)];
            cloud.x.push_back(ri * std::cos(theta));
            cloud.y.push_back(ri * std::sin(theta));
            cloud.z.push_back(static_cast<cusfloat>(0.0));
            cloud.w.push_back(wri * wtheta * ri);
        }
    }

    return cloud;
}

/**
 * @brief Runtime-dispatched annulus node generator.
 *
 * Selects the Gauss-Legendre order at runtime and forwards to the
 * corresponding create_annulus_nodes instantiation.
 */
inline AnnulusNodeCloud create_annulus_nodes_runtime(
    cusfloat r0,
    cusfloat r1,
    int theta_np,
    int radial_order,
    AngularRule rule = AngularRule::Trapezoidal)
{
    switch (radial_order)
    {
        case 1:  return create_annulus_nodes<1>(r0, r1, theta_np, rule);
        case 2:  return create_annulus_nodes<2>(r0, r1, theta_np, rule);
        case 3:  return create_annulus_nodes<3>(r0, r1, theta_np, rule);
        case 4:  return create_annulus_nodes<4>(r0, r1, theta_np, rule);
        case 5:  return create_annulus_nodes<5>(r0, r1, theta_np, rule);
        case 6:  return create_annulus_nodes<6>(r0, r1, theta_np, rule);
        case 7:  return create_annulus_nodes<7>(r0, r1, theta_np, rule);
        case 8:  return create_annulus_nodes<8>(r0, r1, theta_np, rule);
        case 9:  return create_annulus_nodes<9>(r0, r1, theta_np, rule);
        case 10: return create_annulus_nodes<10>(r0, r1, theta_np, rule);
        case 11: return create_annulus_nodes<11>(r0, r1, theta_np, rule);
        case 12: return create_annulus_nodes<12>(r0, r1, theta_np, rule);
        case 13: return create_annulus_nodes<13>(r0, r1, theta_np, rule);
        case 14: return create_annulus_nodes<14>(r0, r1, theta_np, rule);
        case 15: return create_annulus_nodes<15>(r0, r1, theta_np, rule);
        case 16: return create_annulus_nodes<16>(r0, r1, theta_np, rule);
        case 17: return create_annulus_nodes<17>(r0, r1, theta_np, rule);
        case 18: return create_annulus_nodes<18>(r0, r1, theta_np, rule);
        case 19: return create_annulus_nodes<19>(r0, r1, theta_np, rule);
        default:
            return create_annulus_nodes<6>(r0, r1, theta_np, rule);
    }
}
