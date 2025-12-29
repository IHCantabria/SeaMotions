
#
# Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
#
# This file is part of SeaMotions Software.
#
# SeaMotions is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SeaMotions is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#

# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
import numpy as np
from numpy import array, ceil, cos, exp, linspace, log, log10, mod, ndarray, pi, sign, sqrt, tan, zeros
from numpy import abs as np_abs
from scipy.integrate import quad
from scipy.special import expi, jv, kn, roots_laguerre, struve, yv

# Import local modules
from hydlib.waves.airy import w2k


def bisection(ft, a, b, abs_err=1e-6, max_iter=100, verbose=False):
    # Initialize method
    fa = ft(a)
    fb = ft(b)
    c = 0.0

    # Iterate to find the solution
    count_iter = 0
    while True:
        # Calculate center solution
        c = (a+b)/2.0
        fc = ft(c)

        if verbose:
            print(f"Iter: {count_iter} - a: {a} - b: {b} - c: {c} - fa: {fa} - fb: {fb} - fc: {fc}")

        # Check for convergence
        if (abs(fc) < abs_err):
            break

        if count_iter > max_iter:
            print(
                "Warning: Bisection method could not find the solution "
                + " with the precision required."
            )
            break

        count_iter += 1

        # Update to new interval
        if fc*fa > 0:
            fa = fc
            a = c
        else:
            fb = fc
            b = c

    return c


def _check_quad_bounds( X: float, Y: float, int_value: list ) -> None:
    if np_abs( int_value[1] ) > 1e-7:
            raise ValueError( f"Error Exceed! - X: {X:0.3E} - Y: {Y:0.3E} - Int.Value: {int_value[0]:0.2E} - Error: {int_value[1]:0.2E}" )


def complex_quadrature(func: Callable, a: float, b: float, **kwargs):
    # Define functions to parametrize the segment
    g = lambda t: a + t*(b-a)
    gp = lambda t: b-a

    # Define auxiliary funcitons to compute the numerical
    # integral
    def real_func(x):
        return (func(g(x))*gp(x)).real
    def imag_func(x):
        return (func(g(x))*gp(x)).imag
    
    # Add precision keyword arguments
    kwargs[ "limit" ]   = 10000
    kwargs[ "epsabs" ]  = 1e-14
    kwargs[ "epsrel" ]  = 1e-14
    
    # Calculate numerical functions
    real_integral = quad(real_func, 0.0, 1.0, **kwargs)
    imag_integral = quad(imag_func, 0.0, 1.0, **kwargs)

    _check_quad_bounds( a, b, real_integral )
    _check_quad_bounds( a, b, imag_integral )

    # Compose output tuple
    out = (
            complex(real_integral[0], imag_integral[0]),
            real_integral[1],
            imag_integral[1]
            )

    return out


def complex_quadrature_line(func: Callable, way_points: ndarray, close_path=False, **kwargs):
    out = [0+0j, 0.0, 0.0]
    num_points = way_points.shape[0]
    num_points_ref = num_points
    if not close_path:
        num_points_ref -= 1

    for i in range(num_points_ref):
        # Calculate integral along the segment i
        j = int(mod(i+1, num_points))
        out_i = complex_quadrature(func, way_points[i], way_points[j], **kwargs)

        # Sum segements solution
        out[0] = out[0] + out_i[0]
        out[1] = out[1] + out_i[1]
        out[2] = out[2] + out_i[2]

    return out


def fuh(u: ndarray, H: float) -> ndarray:
    return (u+H)/((u-H)-(u+H)*exp(-2*u))


def fxy(X: float, Y: float, only_int=False) -> float:
    # Calculate integral value
    ref_value = 8000
    if Y < ref_value:
        int_value = quad(lambda t: wave_term_expint_def(X, Y, t), 0, Y, limit=10000, epsabs=1e-12, epsrel=1e-14 )
        _check_quad_bounds( X, Y, int_value )
        
    else:
        num_points = int(ceil(Y/ref_value))+2
        int_points = linspace(0, Y, num_points)
        int_value = [0.0, 0.0]
        for i in range(num_points-1):
            int_value_i = quad(lambda t: wave_term_expint_def(X, Y, t), int_points[i], int_points[i+1], limit=10000, epsabs=1e-12, epsrel=1e-14)
            _check_quad_bounds( int_points[i], int_points[i+1], int_value_i )
            int_value[0] = int_value[0] + int_value_i[0]
            int_value[1] = int_value[1] + int_value_i[1]

    # Calculate function value
    if only_int:
        fun_val = int_value[0]

    else:
        fun_val = (
                    -pi*exp(-Y)*(struve(0, X)+yv(0, X))
                    -2*int_value[0]
                    )

    return fun_val


def fxy_polar( R: float, theta: float, only_int=False ) -> float:
    X = R * np.cos( theta )
    Y = R * np.sin( theta )

    return np.log10( fxy( X, Y, only_int=only_int ) )


def fxy_dx(X: float, Y: float, only_int=False) -> float:
    # Calculate the numerical value
    ref_value = 150
    if Y < 1e-8:
        int_value = [0, 0]
    elif Y < ref_value:
        if X < 0.5:
            num_points = 70
            int_points = 10**linspace(-30, log10(Y), num_points)
            # int_points = linspace(0, Y, num_points)
            int_value = [0.0, 0.0]
            for i in range(num_points-1):
                int_value_i = quad(lambda t: wave_term_expint_def_dxt(X, Y, t), int_points[i], int_points[i+1], limit=1000000, epsabs=1e-12, epsrel=1e-14)
                _check_quad_bounds( int_points[i], int_points[i+1], int_value_i )
                int_value[0] = int_value[0] + int_value_i[0]
                int_value[1] = int_value[1] + int_value_i[1]
        else:
            num_points = 10
            int_points = linspace(0, Y, num_points)
            int_value = [0.0, 0.0]
            for i in range(num_points-1):
                int_value_i = quad(lambda t: wave_term_expint_def_dxt(X, Y, t), int_points[i], int_points[i+1], limit=10000, epsabs=1e-12, epsrel=1e-14)
                _check_quad_bounds( int_points[i], int_points[i+1], int_value_i )
                int_value[0] = int_value[0] + int_value_i[0]
                int_value[1] = int_value[1] + int_value_i[1]
    else:
        num_points = int(ceil(Y/ref_value))+20
        int_points = linspace(0, Y, num_points)
        int_value = [0.0, 0.0]
        for i in range(num_points-1):
            int_value_i = quad(lambda t: wave_term_expint_def_dxt(X, Y, t), int_points[i], int_points[i+1], limit=10000, epsabs=1e-12, epsrel=1e-14)
            _check_quad_bounds( int_points[i], int_points[i+1], int_value_i )
            int_value[0] = int_value[0] + int_value_i[0]
            int_value[1] = int_value[1] + int_value_i[1]

    # Calculate function value
    if only_int:
        fun_val = int_value[0]
    else:
        fun_val = (
                    - 2.0/X*int_value[0]
                    + (2/X)*Y/sqrt(X**2.0+Y**2.0)
                    + pi*exp(-Y)*(yv(1, X)+struve(1, X)-2/pi)
                    )

    return fun_val


def fxy_dx_polar( R: float, theta: float, only_int=False ) -> float:
    X = R * np.cos( theta )
    Y = R * np.sin( theta )

    return fxy_dx( X, Y, only_int=only_int )


def fxy_dy(X: float, Y: float) -> float:
    return -(2/sqrt(X**2.0+Y**2.0)+fxy(X, Y))


def G_integral_input_params(R: float, z: float, zeta: float, T: float, h: float) -> list:
    # Calculate required input parameters
    g = 9.81
    w = 2*pi/T
    K = w**2.0/g

    # Calculate derivative input parameters
    A = R/h
    H = K*h
    X = K*R
    v = zeros((6, ))
    v[0] = np_abs(z-zeta)
    v[1] = z+zeta+2*h
    v[2] = np_abs(z+zeta)
    v[3] = z-zeta+2*h
    v[4] = zeta-z+2*h
    v[5] = z+zeta+4*h
    B = v/h
    r = sqrt(R**2.0+v**2.0)

    return g, w, K, A, H, X, v, B, r


def G_integral_dA_input_params(R: float, z: float, zeta: float, T: float, h: float) -> list:
    # Define the extra input params for the dA derivative
    dA_dR = 1/h
    
    # Extend input parameter list
    return_args = list(G_integral_input_params(R, z, zeta, T, h))
    return_args.append(dA_dR)

    return return_args


def G_integral_dB_input_params(R: float, z: float, zeta: float, T: float, h: float) -> list:
    # Define the extra input params for the dB derivative
    dv_dz = zeros((6, ))
    dv_dz[0] = sign(z-zeta)
    dv_dz[1] = 1
    dv_dz[2] = sign(z+zeta)
    dv_dz[3] = 1
    dv_dz[4] = -1
    dv_dz[5] = 1
    dB_dz = dv_dz/h

    # Extend input parameter list
    return_args = list(G_integral_input_params(R, z, zeta, T, h))
    return_args.append(dB_dz)

    return return_args


def G_integral(R: float, z: float, zeta: float, T: float, h: float) -> float:
    # Get derivative input arguments
    g, w, K, A, H, X, v, B, r = G_integral_input_params(R, z, zeta, T, h)

    # Include radius sumation to green function
    G = (1/r).sum()

    k0 = w2k(w, h, method="bisection")
    u0 = k0*h
    
    # Calculate real part
    if B[1] <= 1:
        G += (1/h)*(G1(A, B[0], H)+G1(A, B[1], H))
    else:
        G += (1/h)*(G1(A, B[0], H)+G2(A, B[1], H)) + K*fxy(X, K*v[2])
    
    # Calculate john series
    expsum = (
            + exp(-u0*(2+B[0]))                
            + exp(-u0*(2-B[0]))                
            + exp(-u0*(2+B[1]))                
            + exp(-u0*(2-B[1]))          
            )
    Ck = -pi*(u0+H)/(1+exp(-2*u0)*(2*(u0+H)-1.0))
    gi = Ck*jv(0, u0*A)*expsum/h

    G = G + gi*1j

    return G


def G_integral_dR(R: float, z: float, zeta: float, T: float, h: float) -> float:
    # Get derivative input arguments
    g, w, K, A, H, X, v, B, r, dA_dR = G_integral_dA_input_params(R, z, zeta, T, h)

    # Include radius sumation to green function
    G = (-R/r**3.0).sum()

    k0 = w2k(w, h, method="bisection")
    u0 = k0*h
    
    # Calculate real part
    if B[1] <= 1:
        G += (1/h)*(G1_dA(A, B[0], H)+G1_dA(A, B[1], H))*dA_dR
    else:
        G += (1/h)*(G1_dA(A, B[0], H)+G2_dA(A, B[1], H))*dA_dR + K**2.0*fxy_dx(X, K*v[2])
    
    # Calculate john series
    expsum = (
            + exp(-u0*(2+B[0]))                
            + exp(-u0*(2-B[0]))                
            + exp(-u0*(2+B[1]))                
            + exp(-u0*(2-B[1]))          
            )
    Ck = -pi*(u0+H)/(1+exp(-2*u0)*(2*(u0+H)-1.0))
    gi = -u0*Ck*jv(1, u0*A)*expsum/h*dA_dR

    G = G + gi*1j

    return G


def G_integral_dz(R: float, z: float, zeta: float, T: float, h: float) -> float:
    # Get derivative input arguments
    g, w, K, A, H, X, v, B, r, dB_dz = G_integral_dB_input_params(R, z, zeta, T, h)

    # Include radius sumation to green function
    G = (-v*dB_dz*h/r**3.0).sum()

    k0 = w2k(w, h, method="bisection")
    u0 = k0*h
    
    # Calculate real part
    if B[1] <= 1:
        G += (1/h)*(G1_dB(A, B[0], H)*dB_dz[0]+G1_dB(A, B[1], H)*dB_dz[1])
    else:
        G += (1/h)*(G1_dB(A, B[0], H)*dB_dz[0]+G2_dB(A, B[1], H)*dB_dz[1]) + K**2.0*fxy_dy(X, K*v[2])*sign(z+zeta)

    # Calculate john series
    expsum = (
            + exp(-u0*(2+B[0]))*dB_dz[0]              
            - exp(-u0*(2-B[0]))*dB_dz[0]       
            + exp(-u0*(2+B[1]))*dB_dz[1]        
            - exp(-u0*(2-B[1]))*dB_dz[1]  
            )
    Ck = -pi*(u0+H)/(1+exp(-2*u0)*(2*(u0+H)-1.0))
    gi = -u0*Ck*jv(0, u0*A)*expsum/h

    G = G + gi*1j

    return G


def G_series(R: float, z: float, zeta: float, T: float, h: float) -> float:
    # Calculate required input parameters
    g = 9.81
    w = 2*pi/T
    K = w**2.0/g
    H = K*h

    # Calculate derivative input parameters
    v = zeros((6, ))
    v[0] = np_abs(z-zeta)
    v[1] = z+zeta+2*h
    v[2] = np_abs(z+zeta)
    v[3] = z-zeta+2*h
    v[4] = zeta-z+2*h
    v[5] = z+zeta+4*h
    B = v/h

    # Calculate real wave number of the dispersion relation
    k0 = w2k(w, h, method="bisection")

    # Calculate imaginary wave numbers for the dispersion relation
    nk = 30
    ki_root = w2ki(w, h, nk)

    a = -2*pi*(k0**2.0)*(exp(-k0*v[2:]).sum())*(yv(0, k0*R)+jv(0,k0*R)*1j)
    b = ((k0**2.0-K**2.0)*h+K)*(1+exp(-2*k0*h))**2.0
    G = a/b

    Gs = 0.0
    for i in range(nk):
        Cn = (ki_root[i]**2.0+K**2.0)/((ki_root[i]**2.0+K**2.0)*h-K)
        a = Cn*cos(ki_root[i]*(z+h))*cos(ki_root[i]*(zeta+h))*kn(0, ki_root[i]*R)
        gi = 4*a
        Gs += gi

        if np_abs(gi) < 1e-7:
            break
    G += Gs

    return G


def G1(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = L1(A, B, H)[0] + L2(H)[0]
    else:
        int_value = (
            L3(A, B, H)[0]
            -2/sqrt(A**2.0+(2+B)**2.0)
            -2/sqrt(A**2.0+(2-B)**2.0)
            )
    return int_value.real


def G1_dA(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = L1_dA(A, B, H)[0] + L2_dA()
    else:
        int_value = (
            L3_dA(A, B, H)[0]
            +2*A/(A**2.0+(2+B)**2.0)**(3.0/2.0)
            +2*A/(A**2.0+(2-B)**2.0)**(3.0/2.0)
            )
    return int_value.real


def G1_dB(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = L1_dB(A, B, H)[0] + L2_dB()
    else:
        int_value = (
            L3_dB(A, B, H)[0]
            +2*(2+B)/(A**2.0+(2+B)**2.0)**(3.0/2.0)
            -2*(2-B)/(A**2.0+(2-B)**2.0)**(3.0/2.0)
            )
    return int_value.real


def G2(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = M1(A, B, H)[0]+M2(H)[0]
    else:
        int_value = (
            M3(A, B, H)[0]
            -2/sqrt(A**2.0+(2+B)**2.0)
        )

    return int_value.real


def G2_dA(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = M1_dA(A, B, H)[0]+M2_dA()
    else:
        int_value = (
            M3_dA(A, B, H)[0]
            +2*A/(A**2.0+(2+B)**2.0)**(3.0/2.0)
        )

    return int_value.real


def G2_dB(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = M1_dB(A, B, H)[0]+M2_dB()
    else:
        int_value = (
            M3_dB(A, B, H)[0]
            +2*(2+B)/(A**2.0+(2+B)**2.0)**(3.0/2.0)
        )

    return int_value.real


def guh(u: ndarray, H: float) -> ndarray:
    return (u+H)**2.0/((u-H)**2.0-(u**2.0-H**2.0)*exp(-2*u))


def L0(A: float, B: float) -> float:
    # Define way point for the integration
    way_points = array([
                            1e-16+0j,
                            1e1+0j,
                            1e2+0j,
                            1e3+0j,
                            1e4+0j,
                            1e5+0j
                        ])

    # Integrate function
    l0_def_dummy = lambda u: L0_def(A, B, u)
    int_value = complex_quadrature_line(l0_def_dummy, way_points)

    return int_value


def L0_dA(A: float, B: float) -> float:
    # Define way point for the integration
    way_points = array([
                            1e-16+0j,
                            1e1+0j,
                            1e2+0j,
                            1e3+0j,
                            1e4+0j,
                            1e5+0j
                        ])

    # Integrate function
    l0_da_def_dummy = lambda u: L0_dA_def(A, B, u)
    int_value = complex_quadrature_line(l0_da_def_dummy, way_points)

    return int_value


def L0_dB(A: float, B: float) -> float:
    # Define way point for the integration
    way_points = array([
                            1e-16+0j,
                            1e1+0j,
                            1e2+0j,
                            1e3+0j,
                            1e4+0j,
                            1e5+0j
                        ])

    # Integrate function
    l0_db_def_dummy = lambda u: L0_dB_def(A, B, u)
    int_value = complex_quadrature_line(l0_db_def_dummy, way_points)

    return int_value


def L0_def(A: float, B: float, u: ndarray) -> ndarray:
    return (1/(1-exp(-2*u)))*((exp(-u*(4+B))+exp(-u*(4-B)))*jv(0, u*A)-2*exp(-2*u)) - log( 2.0 )


def L0_dA_def(A: float, B: float, u: ndarray) -> ndarray:
    return (-u/(1-exp(-2*u)))*(exp(-u*(4+B))+exp(-u*(4-B)))*jv(1, u*A)


def L0_dB_def(A: float, B: float, u: ndarray) -> ndarray:
    return (u/(1-exp(-2*u)))*(-exp(-u*(4+B))+exp(-u*(4-B)))*jv(0, u*A)


def L1(A: float, B: float, H: float) -> float:
    # Define way point for the integration
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j
    ])

    # Integrate function
    if H <= 1e-6:
        int_value = list(L1_H0_endo(A, B))
        int_value[0] += complex(0.0, 0.0)
    else:
        l1_def_dummy = lambda u: L1_def(A, B, H, u)
        int_value = complex_quadrature_line(l1_def_dummy, way_points)

    return int_value


def L1_dA(A: float, B: float, H: float) -> float:
    # Define way point for the integration
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j
    ])

    # Integrate function
    if H <= 1e-17:
        raise ValueError("L1_dA Not Implemented for H <= 1e-17")
    else:
        l1_def_dummy = lambda u: L1_dA_def(A, B, H, u)
        int_value = complex_quadrature_line(l1_def_dummy, way_points)

    return int_value


def L1_dB(A: float, B: float, H: float) -> float:
    # Define way point for the integration
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j
    ])

    # Integrate function
    if H <= 1e-17:
        raise ValueError("L1_dB Not Implemented for H <= 1e-17")
    else:
        l1_def_dummy = lambda u: L1_dB_def(A, B, H, u)
        int_value = complex_quadrature_line(l1_def_dummy, way_points)

    return int_value


def L1_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return (fuh(u, H)-1)*((exp(-u*(2+B))+exp(-u*(2-B)))*jv(0, u*A)-2*exp(-2*u))


def L1_dA_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return -u*(fuh(u, H)-1)*((exp(-u*(2+B))+exp(-u*(2-B)))*jv(1, u*A))


def L1_dB_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return -u*(fuh(u, H)-1)*((exp(-u*(2+B))-exp(-u*(2-B)))*jv(0, u*A))


def L1_def_H0(A: float, B: float, H: float, u: ndarray) -> ndarray:
    u = np_abs(u)
    return ((exp(-u*(4+B))+exp(-u*(4-B)))*jv(0, u*A)-2*exp(-4*u))/(1-exp(-2*u))


def L1_H0(A: float, B: float, H: float) -> ndarray:
    # Define way point for the integration
    way_points = array([
        -1e6+0j,
        -1e5+0j,
        -1e4+0j,
        -1e3+0j,
        -1e2+0j,
        -20+0j,
        0.0+7j,
        20+0j,
        1e2+0j,
        1e3+0j,
        1e4+0j,
        1e5+0j,
        1e6+0j,
        0.0+20j,
        -1e6+0j
    ])

    # Integrate function
    l1_def_dummy = lambda u: L1_def_H0(A, B, H, u)
    int_value = complex_quadrature_line(l1_def_dummy, way_points)
    int_value[0] = int_value[0]/2.0 + 1.0+0j

    return int_value


def L1_H0_endo(A: float, B: float) -> float:
    # Get Gauss-Laguerre roots and weigths for the numerical
    # integration
    gl_roots, gl_weights = roots_laguerre(20)

    # Define function components
    fx = lambda u: (exp(-u*(3+B))+exp(-u*(3-B)))*jv(0, u*A)-2*exp(-3*u)
    gx = lambda u: 1-exp(-2*u)
    gxp = lambda u: 2*exp(-2*u)
    fxv = lambda u: (exp(-(u-1)*(3+B))+exp(-(u-1)*(3-B)))*jv(0, (u-1)*A)-2*exp(-3*(u-1))
    gxpv = lambda u: 2*exp(-2*(u-1))

    # Define principal functions
    u0 = 0.0
    u0v = 1.0
    Fx = lambda u: fx(u)/gx(u)
    Fl = lambda u: fx(u0)/((u-u0)*gxp(u0))
    Fp = lambda u: Fx(u)-Fl(u)

    # Integrate principal function using Gauss-Laguerre
    I1 = (gl_weights*Fp(gl_roots)).sum()
    I2 = -exp(-u0v)*expi(u0v)*fxv(u0v)/gxpv(u0v)

    return complex((I1+I2), 0.0), 0.0, 0.0


def L1_pole(A: float, B: float, H: float, u0: float) -> float:
    a = (u0+H)*((exp(-u0*(2+B))+exp(-u0*(2-B)))*jv(0, u0*A)-2.0*exp(-2.0*u0))
    b = 1+(2.0*(u0+H)-1.0)*exp(-2*u0)
    return a/b*pi


def L2(H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
        1e7*H1+0j
    ])

    # Integrate function
    if H<=1e-6:
        int_value = (complex(-log(H)/2-log(2)-1, 0.0), 0.0, 0.0)
    else:
        f_def_dummy = lambda u: L2_def(H, u)
        int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def L2_dA() -> float:
    return 0.0


def L2_dB() -> float:
    return 0.0


def L2_def(H: float, u: ndarray) -> ndarray:
    return (fuh(u, H)-1)*2.0*exp(-2*u)


def L2_pole(H: float, u0: float) -> float:
    return (u0+H)*2.0*exp(-2.0*u0)/(1+(2.0*(u0+H)-1.0)*exp(-2*u0))*pi


def L3(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: L3_def(A, B, H, u)
    int_value   = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def L3_dA(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: L3_dA_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def L3_dB(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: L3_dB_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def L3_def(A: float, B: float, H: float, u: ndarray) -> float:
    return (fuh(u, H)+1)*(exp(-u*(2+B))+exp(-u*(2-B)))*jv(0, u*A)


def L3_dA_def(A: float, B: float, H: float, u: ndarray) -> float:
    return -u*(fuh(u, H)+1)*(exp(-u*(2+B))+exp(-u*(2-B)))*jv(1, u*A)


def L3_dB_def(A: float, B: float, H: float, u: ndarray) -> float:
    return -u*(fuh(u, H)+1)*(exp(-u*(2+B))-exp(-u*(2-B)))*jv(0, u*A)


def L3_pole(A: float, B: float, H: float, u0: float) -> float:
    a = (u0+H)*(exp(-u0*(2+B))+exp(-u0*(2-B)))*jv(0, u0*A)
    b = 1+(2.0*(u0+H)-1.0)*exp(-2*u0)
    return a/b*pi


def Linf(A: float, B: float) -> float:
    # Define way point for the integration
    way_points = array([
                            0.0+0j,
                            1e1+0j,
                            1e2+0j,
                            1e3+0j,
                            1e4+0j,
                            1e5+0j
                        ])

    # Integrate function
    linf_def_dummy = lambda u: Linf_def(A, B, u)
    int_value = complex_quadrature_line(linf_def_dummy, way_points)

    return int_value


def Linf_dA(A: float, B: float) -> float:
    # Define way point for the integration
    way_points = array([
                            0.0+0j,
                            1e1+0j,
                            1e2+0j,
                            1e3+0j,
                            1e4+0j,
                            1e5+0j
                        ])

    # Integrate function
    linf_da_def_dummy = lambda u: Linf_dA_def(A, B, u)
    int_value = complex_quadrature_line(linf_da_def_dummy, way_points)

    return int_value


def Linf_dB(A: float, B: float) -> float:
    # Define way point for the integration
    way_points = array([
                            0.0+0j,
                            1e1+0j,
                            1e2+0j,
                            1e3+0j,
                            1e4+0j,
                            1e5+0j
                        ])

    # Integrate function
    linf_db_def_dummy = lambda u: Linf_dB_def(A, B, u)
    int_value = complex_quadrature_line(linf_db_def_dummy, way_points)

    return int_value


def Linf_def(A: float, B: float, u: ndarray) -> ndarray:
    return (1/(1+exp(-2*u)))*(exp(-u*(4+B))+exp(-u*(4-B)))*jv(0, u*A)


def Linf_dA_def(A: float, B: float, u: ndarray) -> ndarray:
    return (-u/(1+exp(-2*u)))*(exp(-u*(4+B))+exp(-u*(4-B)))*jv(1, u*A)


def Linf_dB_def(A: float, B: float, u: ndarray) -> ndarray:
    return (u/(1+exp(-2*u)))*(-exp(-u*(4+B))+exp(-u*(4-B)))*jv(0, u*A)


def L_minus(A: float, B: float) -> float:
    # Get Gauss-Laguerre points
    gl_roots, gl_weights = roots_laguerre(20)

    # Integrate L_minus function using Gauss-Laguerre
    int_value = (gl_weights*L_minus_def(A, B, gl_roots)).sum()

    return int_value


def L_minus_def(A: float, B: float, u: ndarray) -> ndarray:
    return (exp(-u*(4-B))+exp(-u*(4+B)))*jv(0, u*A)/(1+exp(-2*u))


def M1(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: M1_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M1_dA(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: M1_dA_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M1_dB(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: M1_dB_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M1_def(A: float, B: float, H: float, u: ndarray) -> float:
    c0 = (fuh(u, H)-1)*(exp(-u*(2+B))*jv(0, u*A)-exp(-3*u))
    c1 = guh(u, H)*(exp(-u*(4-B))*jv(0, u*A)-exp(-3*u))
    return c0+c1


def M1_dA_def(A: float, B: float, H: float, u: ndarray) -> float:
    c0 = (fuh(u, H)-1)*(exp(-u*(2+B))*jv(1, u*A))
    c1 = guh(u, H)*(exp(-u*(4-B))*jv(1, u*A))
    return -u*(c0+c1)


def M1_dB_def(A: float, B: float, H: float, u: ndarray) -> float:
    c0 = (fuh(u, H)-1)*(exp(-u*(2+B))*jv(0, u*A))
    c1 = guh(u, H)*(exp(-u*(4-B))*jv(0, u*A))
    return -u*(c0-c1)


def M1_pole(A: float, B: float, H: float, u0: float) -> float:
    # Add u0 fuh
    a = (u0+H)*(exp(-u0*(2+B))*jv(0, u0*A)-exp(-3*u0))
    b = 1+(2*(u0+H)-1)*exp(-2*u0)
    u0_f = a/b

    # Add u0 guh
    a = (exp(-u0*(4-B))*jv(0, u0*A)-exp(-3*u0))*(u0+H)**2.0
    b = 2*(u0-H)+2*exp(-2*u0)*(u0**2.0-u0-H**2.0)
    u0_g = a/b

    # Add H guh
    hg = 2*(exp(-H)-exp(-H*(2-B))*jv(0, H*A))


    return (u0_f+u0_g+hg)*pi


def M2(H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])

    # Integrate function
    f_def_dummy = lambda u: M2_def(H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M2_dA() -> float:
    return 0.0


def M2_dB() -> float:
    return 0.0


def M2_def(H: float, u: ndarray) -> ndarray:
    return (fuh(u, H) - 1.0 + guh(u, H))*exp(-3*u)


def M2_pole(H: float, u0: float) -> float:
    # Add u0 fuh
    a = (u0+H)*exp(-3*u0)
    b = 1+(2*(u0+H)-1)*exp(-2*u0)
    u0_f = a/b

    # Add u0 guh
    a = exp(-3*u0)*(u0+H)**2.0
    b = 2*(u0-H)+2*exp(-2*u0)*(u0**2.0-u0-H**2.0)
    u0_g = a/b

    # Add h guh
    hg = -H*exp(-H)

    return (u0_f+u0_g+hg)*pi


def M3(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])
    
    # Integrate function
    f_def_dummy = lambda u: M3_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M3_dA(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])
    
    # Integrate function
    f_def_dummy = lambda u: M3_dA_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M3_dB(A: float, B: float, H: float) -> float:
    H1 = H+1
    way_points = array([
        0.0+0j,
        H+1j,
        H1+0j,
        1e2*H1+0j,
        1e3*H1+0j,
        1e4*H1+0j,
        1e5*H1+0j,
        1e6*H1+0j,
    ])
    
    # Integrate function
    f_def_dummy = lambda u: M3_dB_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


def M3_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return ((fuh(u, H)+1)*exp(-u*(2+B))+guh(u, H)*exp(-u*(4-B)))*jv(0, u*A)


def M3_dA_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return -u*((fuh(u, H)+1)*exp(-u*(2+B))+guh(u, H)*exp(-u*(4-B)))*jv(1, u*A)


def M3_dB_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return -u*((fuh(u, H)+1)*exp(-u*(2+B))-guh(u, H)*exp(-u*(4-B)))*jv(0, u*A)


def M3_pole(A: float, B: float, H: float, u0: float) -> float:
    # Add u0 fuh
    a = (u0+H)*(exp(-u0*(2+B))*jv(0, u0*A))
    b = 1+(2*(u0+H)-1)*exp(-2*u0)
    u0_f = a/b

    # Add u0 guh
    a = exp(-u0*(4-B))*jv(0, u0*A)*(u0+H)**2.0
    b = 2*(u0-H)+2*exp(-2*u0)*(u0**2.0-u0-H**2.0)
    u0_g = a/b

    # Add h guh
    hg = -2*H*exp(-H*(2-B))*jv(0, H*A)

    return (u0_f+u0_g+hg)*pi


def wave_term_expint_def(X: float, Y: float, t: ndarray) -> ndarray:
    return exp(t-Y)/sqrt(X**2.0+t**2.0)


def wave_term_expint_def_dx(X: float, Y: float, t: ndarray) -> ndarray:
    return exp(t-Y)*X/(X**2.0+t**2.0)**(3.0/2.0)

def wave_term_expint_def_dxt(X: float, Y: float, t: ndarray) -> ndarray:
    return exp(t-Y)*t/(X**2.0+t**2.0)**(1.0/2.0)


def w2ki(w, h, n, g=9.81, abs_err=1e-6, max_iter=1000):
    # Calculate nu
    nu = w**2.0/g

    # Define relation dispersion equation
    kn = zeros((n, ))
    f_dispersion = lambda u: nu+u/h*tan(u)
    for i in range(1, n+1):
        # Define limits
        a = (i-0.5)*pi+1e-12
        b = i*pi

        # Find root using bisection method
        kn[i-1] = bisection(f_dispersion, a, b, abs_err=abs_err, max_iter=max_iter, verbose=False)/h

    return kn
