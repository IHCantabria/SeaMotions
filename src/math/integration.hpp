
#ifndef __integration_math_hpp
#define __integration_math_hpp

// Include general usage libraries
#include <functional>

// Include local modules
#include "../config.hpp"


                                cuscomplex  complex_integration(
                                                                            std::function <cuscomplex(cuscomplex)> f_def,
                                                                            cuscomplex a,
                                                                            cuscomplex b,
                                                                            cusfloat tol
                                                                );
                                            
                                cuscomplex  complex_integration_path(
                                                                            std::function <cuscomplex(cuscomplex)> f_def,
                                                                            int num_segments,
                                                                            cuscomplex* waypoints,
                                                                            cusfloat tol,
                                                                            bool close_path,
                                                                            bool verbose
                                                                    );

template<typename Functor>      cusfloat    romberg_quadrature(
                                                                            Functor f, 
                                                                            cusfloat a, 
                                                                            cusfloat b, 
                                                                            double precision
                                                                );

// Include templates definition
#include "integration.txx"

#endif