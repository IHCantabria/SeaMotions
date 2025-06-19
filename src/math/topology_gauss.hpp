
#pragma once


template<int NV, int NGP>   inline  void    jacobi_det_2d_gp(
                                                                    cusfloat*   xn,
                                                                    cusfloat*   yn,
                                                                    cusfloat*   jac_det
                                                            );


template<int NV, int NGP>   inline  void    jacobi_mat_2d_gp(
                                                                    cusfloat*   xn,
                                                                    cusfloat*   yn,
                                                                    cusfloat*   jac_mat
                                                            );

                                                    
template<int NV, int NGP>   inline  void    shape_fcn_2d_gp( 
                                                                    cusfloat*   N
                                                            );


template<int NV, int NGP>   inline  void    shape_fcn_der_2d_gp(
                                                                    cusfloat*   dNdxi,
                                                                    cusfloat*   dNdeta
                                                                );


template<int NV, int NGP>   inline  void    shape_fcn_deta_2d_gp( 
                                                                    cusfloat*   N
                                                                );


template<int NV, int NGP>   inline  void    shape_fcn_dxi_2d_gp( 
                                                                    cusfloat*   N
                                                                );


#include "topology_gauss.txx"