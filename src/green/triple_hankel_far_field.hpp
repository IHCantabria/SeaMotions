
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


// Numerical integration for triple Hankel products over [1, inf).
#pragma once

// Include general usage libraries
#include <complex>
#include <map>
#include <tuple>
#include <vector>

// Include local modules
#include "../config.hpp"

struct TripleHankelIO 
{
    cusfloat    r_min = 1.0;
    cusfloat    rtol = 1e-6;
    cusfloat    atol = 1e-10;
    int         max_segments = 10000;
    int         min_segments = 4;
    cusfloat    segment_cycles = 2.0;
    cusfloat    damping = 0.0;
    int         romberg_max_level = 8;
    int         hkind0 = 1;
    int         hkind1 = 1;
    cusfloat    rotated_tail_cycles = 30.0;
    bool        verbose = false;
};

cuscomplex  integrate_triple_hankel(
                                            int                     n1,
                                            int                     n2,
                                            int                     n3,
                                            cusfloat                a,
                                            cusfloat                b,
                                            cusfloat                c,
                                            const TripleHankelIO&   options = { }
                                    );

cuscomplex  integrate_triple_hankel_mod(
                                            int                    n1,
                                            int                    n2,
                                            int                    n3,
                                            cusfloat               a,
                                            cusfloat               b,
                                            cusfloat               c,
                                            const TripleHankelIO&  options = {}
                                        );

std::map<std::tuple<int, int, int>, cuscomplex> integrate_triple_hankel_orders(
                                                                                            const std::vector<int>& orders,
                                                                                            cusfloat                a,
                                                                                            cusfloat                b,
                                                                                            cusfloat                c,
                                                                                            const TripleHankelIO&   options = {}
                                                                                        );
