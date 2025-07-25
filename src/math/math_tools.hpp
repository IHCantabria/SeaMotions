
#ifndef __math_tools_hpp
#define __math_tools_hpp


// Include general usage libraries
#include <functional>
#include <string>

// Include local modules
#include "../config.hpp"

// Define local type to compactify expressions
typedef signed long long sll_type;

//////////////////////////////////////////////
////// MATHEMATICAL CONSTANTS BLOCK //////////
//////////////////////////////////////////////
const cusfloat PI = 3.141592653589793;
const cusfloat EULERGAMMA = 0.577215664901533;


//////////////////////////////////////////////
////////// MACRO DEFINITION BLOCK ////////////
//////////////////////////////////////////////
#define POW2S(x) ((x)*(x))


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
                                            cusfloat    angfreq_to_freq( 
                                                                                    cusfloat angfreq 
                                                                        );
                                                                    
                                            int         assert_complex_equality(
                                                                                    cuscomplex  u, 
                                                                                    cuscomplex  v, 
                                                                                    cusfloat    epsilon
                                                                                );
        
template<typename T>                inline  int         assert_scalar_equality(         
                                                                                    T           &u, 
                                                                                    T           &v, 
                                                                                    T           epsilon
                                                                                );

template<typename T>                inline  int         assert_scalar_equality(         
                                                                                    T           &u, 
                                                                                    T           &v, 
                                                                                    T           abs_eps,
                                                                                    T           rel_eps
                                                                                );

template<typename T>                inline  int         assert_vector_equality(         
                                                                                    int         N, 
                                                                                    T*          u, 
                                                                                    T*          v, 
                                                                                    cusfloat    epsilon
                                                                                );

template<typename T>                inline  int         assert_vector_equality(         
                                                                                    int N, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    cusfloat abs_eps,
                                                                                    cusfloat rel_eps
                                                                                );

template<typename T>                inline  int         assert_vector_equality(
                                                                                    int N, 
                                                                                    int* u, 
                                                                                    int* v, 
                                                                                    int epsilon
                                                                                );
        
                                            void        bisection(
                                                                                    std::function<cusfloat(cusfloat)> f_def, 
                                                                                    cusfloat a, 
                                                                                    cusfloat b, 
                                                                                    cusfloat abs_prec, 
                                                                                    cusfloat rel_prec, 
                                                                                    int max_iter, 
                                                                                    bool verbose,
                                                                                    cusfloat &sol, 
                                                                                    int &info
                                                                    );

                                            cusfloat    check_zero_eps( 
                                                                                    cusfloat value,
                                                                                    cusfloat eps
                                                                    );
        
template<typename T>                inline  void        clear_vector(                    
                                                                                    int n, 
                                                                                    T* vec 
                                                                    );

template<typename T, std::size_t N> inline  void        clear_vector(                    
                                                                                    T* vec 
                                                                    );
        
template<typename T>                inline  void        copy_vector(                    
                                                                                    int n, 
                                                                                    T* reference_vector, 
                                                                                    T* target_vector
                                                                    );

template<typename T, std::size_t N> inline  void        copy_vector(
                                                                                    T* reference_vector, 
                                                                                    T* target_vector
                                                                   );

                                            void        conj_vector(
                                                                                    int             n,
                                                                                    cuscomplex*     u,
                                                                                    cuscomplex*     v
                                                                    );

                                            cusfloat    cos3_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                    );

                                            cusfloat    cos2sin_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                        );

                                            cusfloat    cossin2_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                        );

template<typename T>                inline  void        cross(
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                            );

                                            cusfloat    deg_to_rad(    
                                                                                    cusfloat deg 
                                                                    );

template<typename T>                inline  T           eucledian_dist(
                                                                                    int     np,
                                                                                    T*      v0,
                                                                                    T*      v1
                                                                        );
        
                                            sll_type    factorial(
                                                                                    int n
                                                                );
                                
                                            cusfloat    freq_to_angfreq( 
                                                                                    cusfloat freq 
                                                                        );

template<typename T>                        inline  T*  generate_empty_vector(          
                                                                                    int size
                                                                            );

template<typename T>                        inline  T   l2_norm( 
                                                                                    int N,
                                                                                    T* x,
                                                                                    T* y
                                                                );

                                            void        newton_raphson(
                                                                                    std::function<cusfloat(cusfloat)> f_def, 
                                                                                    std::function<cusfloat(cusfloat)> f_der_def,
                                                                                    cusfloat x0, 
                                                                                    cusfloat fabs_tol, 
                                                                                    cusfloat xrel_tol, 
                                                                                    int max_iter, 
                                                                                    bool verbose, 
                                                                                    cusfloat &sol, 
                                                                                    int &info
                                                                        );

                                            cusfloat    period_to_angfreq( 
                                                                                    cusfloat period 
                                                                        );
        
template<typename T>                inline  T           pow2s(
                                                                                    T x
                                                            );
        
template<typename T>                inline  T           pow3s(
                                                                                    T x
                                                                );
        
template<typename T>                inline  void        print_matrix(
                                                                                    int num_rows, 
                                                                                    int num_cols, 
                                                                                    T* mat, 
                                                                                    int precision, 
                                                                                    int align, 
                                                                                    int scient_not
                                                                    );
        
template<typename T>                inline  void        print_vector(
                                                                                    int n, 
                                                                                    T* v, 
                                                                                    int mode, 
                                                                                    int precision
                                                                    );
        
template <typename T>               inline  int         sign(
                                                                                    T val
                                                                );

                                            cusfloat    sin_alpha( 
                                                                                    int         alpha,
                                                                                    cusfloat    theta
                                                                );

                                            cusfloat    sin3_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                        );
        
template<typename T>                inline  void        slice_vector(
                                                                                    T* parent_vector, 
                                                                                    int first_pos, 
                                                                                    int second_pos, 
                                                                                    T* slice_vector
                                                                    );
        
template<typename T>                inline  void        sv_add(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_add(                                
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T s, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        sv_cbrt(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_div(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        svs_div(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T s,
                                                                                    T* w
                                                                );

template<typename T>                inline  T           sv_dot(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T* v
                                                                );

template<typename T>                inline  void        sv_dot_vm(
                                                                                    int n,
                                                                                    T*  vec,
                                                                                    T*  mat,
                                                                                    T*  vec_out
                                                                );
        
template<typename T>                inline  void        sv_inv(
                                                                                    int n, 
                                                                                    T s, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_mod(                     int n, 
                                                                                    T* u, 
                                                                                    T &mod
                                                                );
        
template<typename T>                inline  void        sv_mult(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_mult(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T s, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        sv_pow(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T* v,
                                                                                    T* w
                                                                );

template<typename T>                inline  void        svs_pow(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T s,
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_pow2(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        sv_pow3(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_sqrt(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_sub(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_sub(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T s, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_sub(
                                                                                    int n, 
                                                                                    T s, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );


#include "math_tools.txx"

#endif