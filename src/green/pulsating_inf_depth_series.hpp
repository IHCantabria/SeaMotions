
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

#pragma once

// Include local modules
#include "../../src/config.hpp"


void domain_inf_fit(cusfloat x, cusfloat y, cusfloat &xl, cusfloat &yl, cusfloat &jac);
cusfloat eval_chebyshev_fit(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y);
cusfloat eval_chebyshev_fit_dxndim(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y);
void get_inf_domain_bounds(cusfloat x, cusfloat y, cusfloat &x0, cusfloat &x1, cusfloat &y0, cusfloat &y1);
cusfloat expint_inf_depth_num(cusfloat X, cusfloat Y);
cusfloat expint_inf_depth_num_dxndim(cusfloat X, cusfloat Y);
cusfloat expint_inf_depth_num_dxtndim(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_series(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_dxndim_series(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_dyndim_series(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num_dxndim(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num_dyndim(cusfloat X, cusfloat Y);
