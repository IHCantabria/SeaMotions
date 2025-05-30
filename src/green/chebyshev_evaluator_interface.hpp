
#ifndef __chebyshev_evaluator_interface_hpp
#define __chebyshev_evaluator_interface_hpp

// Include local modules
#include "chebyshev_evaluator_base.hpp"
#include "chebyshev_traits.hpp"

/*************************************/
/******* Alias Scalar Classes ********/
/*************************************/
using R44CE         = ChebyshevEvaluatorBase<ChebyshevTraits<R44C>>;
using R44_dXCE      = ChebyshevEvaluatorBase<ChebyshevTraits<R44_dXC>>;


/*************************************/
/******* Alias Vector Classes ********/
/*************************************/
template<std::size_t N>
using R44CEV        = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R44C>, N>;

template<std::size_t N>
using R44_dXCEV     = ChebyshevEvaluatorBaseVector<ChebyshevTraits<R44_dXC>, N>;




#endif