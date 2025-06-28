
#ifndef __chebyshev_evaluator_interface_hpp
#define __chebyshev_evaluator_interface_hpp

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
template<std::size_t N>
using L1CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L1C>, N>;

template<std::size_t N>
using L1_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L1_dAC>, N>;

template<std::size_t N>
using L1_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L1_dBC>, N>;

template<std::size_t N>
using L2CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L2C>, N>;

template<std::size_t N>
using L3CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L3C>, N>;

template<std::size_t N>
using L3_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L3_dAC>, N>;

template<std::size_t N>
using L3_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<L3_dBC>, N>;

template<std::size_t N>
using M1CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M1C>, N>;

template<std::size_t N>
using M1_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M1_dAC>, N>;

template<std::size_t N>
using M1_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M1_dBC>, N>;

template<std::size_t N>
using M2CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M2C>, N>;

template<std::size_t N>
using M3CEV         = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M3C>, N>;

template<std::size_t N>
using M3_dACEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M3_dAC>, N>;

template<std::size_t N>
using M3_dBCEV      = ChebyshevEvaluatorBaseVector<ChebyshevTraits<M3_dBC>, N>;

template<std::size_t N>
using R11CEV        = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R11C>, N>;

template<std::size_t N>
using R11_dXCEV     = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R11_dXC>, N>;




#endif