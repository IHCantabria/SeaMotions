
#include "L2.hpp"

int                     L2C::dims = 1;
int                     L2C::num_intervals = 2;
cusfloat                L2C::interval_bounds[3] = {1.000000E-16, 1.000000E-03, 1.000000E+00};

int                     L2C::num_points[2] = {2, 8};
int                     L2C::num_points_cum[3] = {0, 2, 10};
int                     L2C::max_size_fold = 0;

bool                    L2C::x_log_scale[2] = {1, 1};
cusfloat                L2C::x_map_scale[2] = {0.000000E+00, 0.000000E+00};
cusfloat                L2C::x_map_scale_log[2] = {0.000000E+00, 0.000000E+00};
cusfloat                L2C::x_max[2] = {1.000000E-03, 1.000000E+00};
cusfloat                L2C::x_min[2] = {1.000000E-16, 1.000000E-03};
cusfloat                L2C::x_min_l10[2] = {0.000000E+00, 0.000000E+00};

int                     L2C::num_c = 10;
cusfloat                L2C::c[10] = {
                                   9.2439613664761549E+00,  // C[0]
                                   -7.4837300998278042E+00,  // C[1]
                                   -7.2607736162641739E-02,  // C[2]
                                   -1.8899219189923933E+00,  // C[3]
                                   -7.6792902850869413E-02,  // C[4]
                                   -2.0074581194437644E-02,  // C[5]
                                   1.0327991262758091E-03,  // C[6]
                                   4.5060307126852869E-03,  // C[7]
                                   3.1849598412797718E-03,  // C[8]
                                   1.6658728268160985E-03,  // C[9]
                                };
cusfloat                L2C::cf[0] = {
                            };
cusfloat                L2C::cf2[0] = {
                            };
int                     L2C::ncx[10] = {
                                   0,  // ncx[0]
                                   1,  // ncx[1]
                                   0,  // ncx[2]
                                   1,  // ncx[3]
                                   2,  // ncx[4]
                                   3,  // ncx[5]
                                   4,  // ncx[6]
                                   5,  // ncx[7]
                                   6,  // ncx[8]
                                   7,  // ncx[9]
                            };
int                     L2C::ncxf[0] = {
                            };
int                     L2C::ncxf2[0] = {
                            };
