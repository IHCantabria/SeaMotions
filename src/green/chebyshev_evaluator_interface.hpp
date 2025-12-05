
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
#include "chebyshev_traits.hpp"

/*************************************/
/******* Alias Scalar Classes ********/
/*************************************/
using R11CE         = ChebyshevEvaluatorBase<ChebyshevTraits<R11C>>;
using R11_dXCE      = ChebyshevEvaluatorBase<ChebyshevTraits<R11_dXC>>;


/*************************************/
/******* Alias Vector Classes ********/
/*************************************/
template<std::size_t N, int mode_loop>
using L1CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L1C>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using L1_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L1_dAC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using L1_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L1_dBC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using L2CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L2C>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using L3CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L3C>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using L3_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L3_dAC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using L3_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L3_dBC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M1CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M1C>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M1_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M1_dAC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M1_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M1_dBC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M2CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M2C>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M3CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M3C>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M3_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M3_dAC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using M3_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M3_dBC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using R11CEV        = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R11C>, N, mode_loop>;


template<std::size_t N, int mode_loop>
using R00_dXCEV     = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R00_dXC>, N, mode_loop>;

template<std::size_t N, int mode_loop>
using R11_dXCEV     = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R11_dXC>, N, mode_loop>;
