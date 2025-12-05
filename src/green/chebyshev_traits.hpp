
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
#include "chebyshev_traits_macros.hpp"
#include "./fin_depth_coeffs/L1.hpp"
#include "./fin_depth_coeffs/L1_dA.hpp"
#include "./fin_depth_coeffs/L1_dB.hpp"
#include "./fin_depth_coeffs/L2.hpp"
#include "./fin_depth_coeffs/L3.hpp"
#include "./fin_depth_coeffs/L3_dA.hpp"
#include "./fin_depth_coeffs/L3_dB.hpp"
#include "./fin_depth_coeffs/M1.hpp"
#include "./fin_depth_coeffs/M1_dA.hpp"
#include "./fin_depth_coeffs/M1_dB.hpp"
#include "./fin_depth_coeffs/M2.hpp"
#include "./fin_depth_coeffs/M3.hpp"
#include "./fin_depth_coeffs/M3_dA.hpp"
#include "./fin_depth_coeffs/M3_dB.hpp"
#include "./inf_depth_coeffs/R11.hpp"
#include "./inf_depth_coeffs/R00_dX.hpp"
#include "./inf_depth_coeffs/R11_dX.hpp"


// Generic trait declaration
template<typename NS>
struct ChebyshevTraits;

// Expand macros for infinite water depth traits
CHEBYSHEV_2D_TRAITS( R11C );
CHEBYSHEV_2D_TRAITS( R00_dXC );
CHEBYSHEV_2D_TRAITS( R11_dXC );

// Expand macros for finit water depth traits
CHEBYSHEV_3DF_TRAITS( L1C )
CHEBYSHEV_3DF_TRAITS( L1_dAC )
CHEBYSHEV_3DF_TRAITS( L1_dBC )
CHEBYSHEV_1DF_TRAITS( L2C )
CHEBYSHEV_3DF_TRAITS( L3C )
CHEBYSHEV_3DF_TRAITS( L3_dAC )
CHEBYSHEV_3DF_TRAITS( L3_dBC )

CHEBYSHEV_3DF_TRAITS( M1C )
CHEBYSHEV_3DF_TRAITS( M1_dAC )
CHEBYSHEV_3DF_TRAITS( M1_dBC )
CHEBYSHEV_1DF_TRAITS( M2C )
CHEBYSHEV_3DF_TRAITS( M3C )
CHEBYSHEV_3DF_TRAITS( M3_dAC )
CHEBYSHEV_3DF_TRAITS( M3_dBC )
