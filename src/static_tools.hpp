
#pragma once

// Include local modules
#include "./math/math_tools.hpp"

#define IS_STATIC_LOOP                          mode_loop == 0
#define LOOP_DEF( N, code_block )               for ( std::size_t i=0; i<N; i++ ){ code_block }
#define STATIC_CLEAR( n, N, vec )               if constexpr( IS_STATIC_LOOP ) { (clear_vector<cusfloat, N>( vec )); } else { (clear_vector<cusfloat>( n, vec )); }
#define STATIC_COPY( n, N, vec )                if constexpr( IS_STATIC_LOOP ) { (copy_vector<cusfloat, N>( vec )); } else { (copy_vector<cusfloat>( n, vec )); }
#define STATIC_COND( cond, code_block )         if constexpr( cond ) { code_block }
#define STATIC_LOOP( n, N, code_block )         if constexpr( IS_STATIC_LOOP ) { LOOP_DEF( N, code_block ) } else { LOOP_DEF( n, code_block ) }
#define ONLY_FCN                                mode_f == 1
#define ONLY_FCNDR                              mode_dfdr == 1
#define ONLY_FCNDZ                              mode_dfdz == 1