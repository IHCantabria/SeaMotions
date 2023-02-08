
#ifndef __pulsating_hpp
#define __pulsating_hpp

#include "../../src/config.hpp"

void domain_inf_fit(cusfloat x, cusfloat y, cusfloat &xl, cusfloat &yl, cusfloat &jac);
cusfloat eval_chebyshev_fit(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y);
cusfloat eval_chebyshev_fit_dx(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y);
void get_inf_domain_bounds(cusfloat x, cusfloat y, cusfloat &x0, cusfloat &x1, cusfloat &y0, cusfloat &y1);
cusfloat expint_inf_depth_num(cusfloat X, cusfloat Y);
cusfloat expint_inf_depth_num_dx(cusfloat X, cusfloat Y);
cusfloat expint_inf_depth_num_dxt(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_dx(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_dy(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num_dx(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num_dy(cusfloat X, cusfloat Y);

#endif