
#pragma once

// Include local modules
#include "./math/math_tools.hpp"

#define IS_STATIC_LOOP                                      mode_loop == 1
#define IS_LID_PANEL                                        mode_fslid == 1
#define LOOP_DEF( N, code_block )                           for ( std::size_t i=0; i<N; i++ ){ code_block }
#define MPI_TIME_EXEC( block, var_stg )                     { cusfloat t0=MPI_Wtime( ); block cusfloat t1=MPI_Wtime( ); var_stg=t1-t0; }
#define STATIC_CLEAR( n, N, vec )                           if constexpr( IS_STATIC_LOOP ) { (clear_vector<cusfloat, N>( vec )); } else { (clear_vector<cusfloat>( n, vec )); }
#define STATIC_COPY( n, N, vec_1, vec_2 )                   if constexpr( IS_STATIC_LOOP ) { (copy_vector<cusfloat, N>( vec_1, vec_2 )); } else { (copy_vector<cusfloat>( n, vec_1, vec_2 )); }
#define STATIC_COND( cond, code_block )                     if constexpr( cond ) { code_block }
#define STATIC_COND_2( cond, code_block_1, code_block_2 )   if constexpr( cond ) { code_block_1 } else { code_block_2 }
#define STATIC_LOOP( n, N, code_block )                     if constexpr( IS_STATIC_LOOP ) { LOOP_DEF( N, code_block ) } else { LOOP_DEF( n, code_block ) }
#define ONLY_FCN                                            mode_f == 1
#define ONLY_FCNDR                                          mode_dfdr == 1
#define ONLY_FCNDZ                                          mode_dfdz == 1
#define ONLY_PF                                             mode_pf == 1