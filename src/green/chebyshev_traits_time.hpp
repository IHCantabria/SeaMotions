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
#include "chebyshev_traits_macros_time.hpp"
#include "td_coeffs_files.hpp"

// Generic trait declaration (re-declaration is harmless; defined in chebyshev_traits.hpp)
template<typename NS>
struct ChebyshevTraits;

/***************************************************************/
/******* Trait specializations: 2D foldable time types *********/
/*** (dGdt, dGdtt, dGdtx, dGdtxx, dGdttx, dGdttxx)          ***/
/*** Uses CHEBYSHEV_TIME_2DF_TRAITS: exposes full 2D read-only***/
/*** metadata, G0 near-field correction data, and mutable     ***/
/*** inline-static storage for fold_time_residual_x results. ***/
/***************************************************************/
CHEBYSHEV_TIME_2DF_TRAITS( dGdtC   )
CHEBYSHEV_TIME_2DF_TRAITS( dGdttC  )
CHEBYSHEV_TIME_2DF_TRAITS( dGdtxC  )
CHEBYSHEV_TIME_2DF_TRAITS( dGdtxxC )
CHEBYSHEV_TIME_2DF_TRAITS( dGdttxC )
CHEBYSHEV_TIME_2DF_TRAITS( dGdttxxC)

/***************************************************************/
/******* Trait specializations: 1D time types             ******/
/*** (dGdtA0, dGdttA0/A1/A2, dGdtxA0/A1/A2,              ***/
/***  dGdtxxA0/A1/A2, dGdttxA0/A1/A2, dGdttxxA0/A1/A2)  ***/
/*** Uses CHEBYSHEV_1D_TRAITS: exposes 1D region grid,    ***/
/*** rise/fall block arrays, ncx.                        ***/
/***************************************************************/
CHEBYSHEV_1D_TRAITS( dGdtA0C   )

CHEBYSHEV_1D_TRAITS( dGdttA0C  )
CHEBYSHEV_1D_TRAITS( dGdttA1C  )
CHEBYSHEV_1D_TRAITS( dGdttA2C  )

CHEBYSHEV_1D_TRAITS( dGdtxA0C  )
CHEBYSHEV_1D_TRAITS( dGdtxA1C  )
CHEBYSHEV_1D_TRAITS( dGdtxA2C  )

CHEBYSHEV_1D_TRAITS( dGdtxxA0C )
CHEBYSHEV_1D_TRAITS( dGdtxxA1C )
CHEBYSHEV_1D_TRAITS( dGdtxxA2C )

CHEBYSHEV_1D_TRAITS( dGdttxA0C )
CHEBYSHEV_1D_TRAITS( dGdttxA1C )
CHEBYSHEV_1D_TRAITS( dGdttxA2C )

CHEBYSHEV_1D_TRAITS( dGdttxxA0C)
CHEBYSHEV_1D_TRAITS( dGdttxxA1C)
CHEBYSHEV_1D_TRAITS( dGdttxxA2C)
