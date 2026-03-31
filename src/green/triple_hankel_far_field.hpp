
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
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


/**
 * @file triple_hankel_far_field.hpp
 * @brief Numerical integration utilities for far-field triple Hankel products.
 *
 * This module provides numerical evaluators for integrals of the form
 *
 *   I = \int_{r_min}^{\infty} r H_{n1}^{(k0)}(a r) H_{n2}^{(k1)}(b r) H_{n3}^{(2)}(c r) dr,
 *
 * using two strategies:
 * - Direct segmented integration on the real axis.
 * - A mixed method that integrates the near field on the real axis and the
 *   tail with a contour-rotated asymptotic representation.
 *
 * The API is designed for validation and database generation workflows where
 * repeated evaluations are required under configurable tolerance and segmentation
 * controls.
 */
// Numerical integration for triple Hankel products over [1, inf).
#pragma once

// Include general usage libraries
#include <complex>
#include <map>
#include <tuple>
#include <vector>

// Include local modules
#include "../config.hpp"

/**
 * @brief Configuration options for triple Hankel integrations.
 *
 * These options control absolute/relative stopping criteria, integration
 * segmentation, Romberg sub-integration depth, Hankel kinds, and diagnostics.
 *
 * @note Default values are tuned for robust reference computations rather than
 * maximum speed.
 */
struct TripleHankelIO
{
        /** @brief Lower integration bound. Values below 1.0 are accepted by API but
            * practical workflows typically use r_min = 1.0 to avoid near-origin issues.
            */
    cusfloat    r_min = 1.0;

        /** @brief Relative tolerance for segment convergence checks. */
    cusfloat    rtol = 1e-6;

        /** @brief Absolute tolerance for segment convergence checks. */
    cusfloat    atol = 1e-10;

        /** @brief Maximum number of integration segments for finite or tail parts. */
    int         max_segments = 10000;

        /** @brief Minimum number of segments evaluated before early-stop tests. */
    int         min_segments = 4;

        /** @brief Number of dominant oscillation cycles per segment in direct mode. */
    cusfloat    segment_cycles = 2.0;

        /** @brief Optional exponential damping factor on the real-axis integrand. */
    cusfloat    damping = 0.0;

        /** @brief Maximum Romberg extrapolation level used in sub-integrations. */
    int         romberg_max_level = 8;

        /** @brief Hankel kind for the first factor: 1 => H^(1), 2 => H^(2). */
    int         hkind0 = 1;

        /** @brief Hankel kind for the second factor: 1 => H^(1), 2 => H^(2). */
    int         hkind1 = 1;

        /** @brief Tail extent control (in effective oscillatory cycles) for rotated contour setup. */
    cusfloat    rotated_tail_cycles = 30.0;

        /** @brief Enables per-segment diagnostic output when true. */
    bool        verbose = false;
};

/**
 * @brief Integrate the triple Hankel product on the real axis only.
 *
 * The integral is evaluated by splitting [r_min, +inf) into fixed-length
 * segments determined from dominant oscillation scales. Each segment is
 * integrated with Romberg quadrature and accumulated until convergence criteria
 * based on @p atol and @p rtol are met.
 *
 * @param n1 Order of the first Hankel factor.
 * @param n2 Order of the second Hankel factor.
 * @param n3 Order of the third Hankel factor.
 * @param a Wavenumber scaling of the first Hankel factor.
 * @param b Wavenumber scaling of the second Hankel factor.
 * @param c Wavenumber scaling of the third Hankel factor.
 * @param options Numerical integration controls (tolerances, segmentation, kinds).
 *
 * @return Complex integral value over [r_min, +inf).
 */
cuscomplex  integrate_triple_hankel(
                                            int                     n1,
                                            int                     n2,
                                            int                     n3,
                                            cusfloat                a,
                                            cusfloat                b,
                                            cusfloat                c,
                                            const TripleHankelIO&   options = { }
                                    );

                        /**
                         * @brief Integrate the triple Hankel product with a rotated-contour tail treatment.
                         *
                         * This method computes:
                         * - A finite near-field contribution on the real axis using full Hankel evaluations.
                         * - A semi-infinite far-field contribution using asymptotic/contour-rotated form.
                         *
                         * The final value is the sum of both parts including the asymptotic scaling phase
                         * factor used to map the rotated tail back to the target integral.
                         *
                         * @param n1 Order of the first Hankel factor.
                         * @param n2 Order of the second Hankel factor.
                         * @param n3 Order of the third Hankel factor.
                         * @param a Wavenumber scaling of the first Hankel factor.
                         * @param b Wavenumber scaling of the second Hankel factor.
                         * @param c Wavenumber scaling of the third Hankel factor.
                         * @param options Numerical integration controls, including rotated-tail settings.
                         *
                         * @return Complex integral value over [r_min, +inf).
                         */
cuscomplex  integrate_triple_hankel_mod(
                                            int                    n1,
                                            int                    n2,
                                            int                    n3,
                                            cusfloat               a,
                                            cusfloat               b,
                                            cusfloat               c,
                                            const TripleHankelIO&  options = {}
                                        );

/**
 * @brief Evaluate triple Hankel integrals for all combinations of selected orders.
 *
 * For each tuple (n1, n2, n3) in @p orders x @p orders x @p orders, this function
 * calls integrate_triple_hankel and stores the result in a map indexed by the
 * corresponding order tuple.
 *
 * @param orders Set of Hankel orders used for all three factors.
 * @param a Wavenumber scaling of the first Hankel factor.
 * @param b Wavenumber scaling of the second Hankel factor.
 * @param c Wavenumber scaling of the third Hankel factor.
 * @param options Numerical integration controls passed through to each evaluation.
 *
 * @return Map from (n1, n2, n3) to complex integral value.
 */
std::map<std::tuple<int, int, int>, cuscomplex> integrate_triple_hankel_orders(
                                                                                            const std::vector<int>& orders,
                                                                                            cusfloat                a,
                                                                                            cusfloat                b,
                                                                                            cusfloat                c,
                                                                                            const TripleHankelIO&   options = {}
                                                                                        );
