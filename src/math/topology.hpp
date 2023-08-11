
#ifndef __topology_hpp
#define __topology_hpp

// Include local modules
#include "../config.hpp"


cusfloat    jacobi_det_2d(
                                    int         np,
                                    cusfloat*   xn,
                                    cusfloat*   yn,
                                    cusfloat    xi,
                                    cusfloat    eta
                        );

void        jacobi_mat_2d(
                                    int         np,
                                    cusfloat*   xn,
                                    cusfloat*   yn,
                                    cusfloat    xi,
                                    cusfloat    eta,
                                    cusfloat*   jac_mat
                        );

void        shape_fcn_2d( 
                                    int         np, 
                                    cusfloat    xi, 
                                    cusfloat    eta,
                                    cusfloat*   N
                           );

void        shape_fcn_der_2d(
                                    int         np,
                                    cusfloat    xi,
                                    cusfloat    eta,
                                    cusfloat*   dNdxi,
                                    cusfloat*   dNdeta
                        );

void        shape_fcn_deta_2d( 
                                    int         np, 
                                    cusfloat    xi, 
                                    cusfloat    eta,
                                    cusfloat*   N
                            );

void        shape_fcn_dxi_2d( 
                                    int         np, 
                                    cusfloat    xi, 
                                    cusfloat    eta,
                                    cusfloat*   N
                        );

void        simplex_to_cartesian( 
                                   cusfloat    x, 
                                   cusfloat    y, 
                                   cusfloat&   xm,
                                   cusfloat&   ym
                               );

#endif