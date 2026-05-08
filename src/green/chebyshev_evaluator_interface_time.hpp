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
#include "chebyshev_evaluator_base.hpp"
#include "chebyshev_traits_time.hpp"

/***************************************************************/
/******* Alias Vector Classes (2D foldable time types) *********/
/*** Instantiate ChebyshevEvaluatorBaseVector for each       ***/
/*** 2D coefficient set.  Template parameters:              ***/
/***   N         – number of evaluation points per call     ***/
/***   mode_loop – loop-unrolling / SIMD mode selector      ***/
/***************************************************************/

template<std::size_t N, int mode_loop>
using dGdtCEV       = ChebyshevEvaluatorBaseVector<ChebyshevTraits<dGdtC>,    N, mode_loop>;

template<std::size_t N, int mode_loop>
using dGdttCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<dGdttC>,   N, mode_loop>;

template<std::size_t N, int mode_loop>
using dGdtxCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<dGdtxC>,   N, mode_loop>;

template<std::size_t N, int mode_loop>
using dGdtxxCEV     = ChebyshevEvaluatorBaseVector<ChebyshevTraits<dGdtxxC>,  N, mode_loop>;

template<std::size_t N, int mode_loop>
using dGdttxCEV     = ChebyshevEvaluatorBaseVector<ChebyshevTraits<dGdttxC>,  N, mode_loop>;

template<std::size_t N, int mode_loop>
using dGdttxxCEV    = ChebyshevEvaluatorBaseVector<ChebyshevTraits<dGdttxxC>, N, mode_loop>;

/***************************************************************/
/******* Alias Scalar Classes (1D time types)            *******/
/*** Instantiate ChebyshevEvaluatorBase for each 1D     ***/
/*** coefficient set (single-point evaluation).         ***/
/***************************************************************/

using dGdtA0CE      = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtA0C>>;

using dGdttA0CE     = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttA0C>>;
using dGdttA1CE     = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttA1C>>;
using dGdttA2CE     = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttA2C>>;

using dGdtxA0CE     = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtxA0C>>;
using dGdtxA1CE     = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtxA1C>>;
using dGdtxA2CE     = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtxA2C>>;

using dGdtxxA0CE    = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtxxA0C>>;
using dGdtxxA1CE    = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtxxA1C>>;
using dGdtxxA2CE    = ChebyshevEvaluatorBase<ChebyshevTraits<dGdtxxA2C>>;

using dGdttxA0CE    = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttxA0C>>;
using dGdttxA1CE    = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttxA1C>>;
using dGdttxA2CE    = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttxA2C>>;

using dGdttxxA0CE   = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttxxA0C>>;
using dGdttxxA1CE   = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttxxA1C>>;
using dGdttxxA2CE   = ChebyshevEvaluatorBase<ChebyshevTraits<dGdttxxA2C>>;
