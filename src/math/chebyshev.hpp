
#pragma once

#include "../config.hpp"

cusfloat    chebyshev_poly(             
                                        int n, 
                                        cusfloat x
                            );

cusfloat    chebyshev_poly_raw(         
                                        int n, 
                                        cusfloat x
                                );

cusfloat    chebyshev_poly_der_raw(     int n,
                                        cusfloat x
                                    );

void        chebyshev_poly_roots(       
                                        int num_points,
                                        cusfloat* roots
                                );

void        get_gauss_chebyshev(        
                                        int num_points, 
                                        cusfloat* weights, 
                                        cusfloat* roots
                                );

void        chebyshev_poly_upto_order(  
                                        std::size_t max_order, 
                                        cusfloat x, 
                                        cusfloat* results 
                                    );


#include "chebyshev.txx"