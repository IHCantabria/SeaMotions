
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
import matplotlib.pyplot as plt
from numpy import array, cos, cosh, exp, linspace, log, log10, meshgrid, mod, ndarray, pi, sin, sinh, sqrt, tan, tanh, zeros
from numpy import abs as np_abs
from scipy.integrate import quad
from scipy.special import expi, jv, kn, roots_chebyt, roots_laguerre, struve, yv

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
    
    # Calculate numerical functions
    real_integral = quad(real_func, 0.0, 1.0, **kwargs)
    imag_integral = quad(imag_func, 0.0, 1.0, **kwargs)

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


def fxy(X: float, Y: float) -> float:
    # Calculate integral value
    int_value = quad(lambda t: wave_term_expint_def(X, Y, t), 0, Y)

    # Caculate function value
    fun_val = (
                -pi*exp(-Y)*(struve(0, X)+yv(0, X))
                -2*int_value[0]
                )

    return fun_val


def G_integral(R: float, z: float, zeta: float, T: float, h: float) -> float:
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
    Ck = -2*pi*k0**2.0/((k0**2.0-K**2.0)*h+K)/(1+exp(-2*k0*h))**2.0
    gi = Ck*jv(0, k0*R)*(exp(-k0*v[2:])).sum()

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


def G2(A: float, B: float, H: float) -> float:
    if H <= 1:
        int_value = M1(A, B, H)[0]+M2(H)[0]
    else:
        int_value = (
            M3(A, B, H)[0]
            -2/sqrt(A**2.0+(2+B)**2.0)
        )

    return int_value.real


def guh(u: ndarray, H: float) -> ndarray:
    return (u+H)**2.0/((u-H)**2.0-(u**2.0-H**2.0)*exp(-2*u))


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
    if H <= 1e-6:
        raise ValueError("L1_dA Not Implemented for H <= 1e-6")
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
    if H <= 1e-6:
        raise ValueError("L1_dB Not Implemented for H <= 1e-6")
    else:
        l1_def_dummy = lambda u: L1_dB_def(A, B, H, u)
        int_value = complex_quadrature_line(l1_def_dummy, way_points)

    return int_value


def L1_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return (fuh(u, H)-1)*((exp(-u*(2+B))+exp(-u*(2-B)))*jv(0, u*A)-2*exp(-2*u))


def L1_dA_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return -u*(fuh(u, H)-1)*((exp(-u*(2+B))+exp(-u*(2-B)))*jv(1, u*A))


def L1_dB_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return -u*(fuh(u, H)-1)*((exp(-u*(2+B))+exp(-u*(2-B)))*jv(0, u*A))


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

    # u = linspace(0, 100, 1000)
    # y = Fp(u)
    # print(f"A: {A:0.3f} - B: {B:0.3f}")
    # plt.plot(u, Fx(u))
    # plt.plot(u, Fl(u))
    # plt.plot(u, Fp(u))
    # plt.show()

    # Integrate principal function using Gauss-Laguerre
    I1 = (gl_weights*Fp(gl_roots)).sum()
    I2 = -exp(-u0v)*expi(u0v)*fxv(u0v)/gxpv(u0v)
    # print("I1: ", I1)
    # print("I2: ", I2)
    # print("-exp(-u0): ", -exp(u0))
    # print("exp1(0): ", exp1(u0))
    # print("fx(u0): ", fx(u0))
    # print("gxp(u0): ", gxp(u0))

    # print("Endo: ", (I1+I2))

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
    # if H <= 1001:
    f_def_dummy = lambda u: L3_def(A, B, H, u)
    int_value = complex_quadrature_line(f_def_dummy, way_points)
    # else:
        # int_value = (complex(L_minus(A, B), 0.0), 0.0, 0.0)

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
    return -u*(c0+c1)


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
    return -u*((fuh(u, H)+1)*exp(-u*(2+B))+guh(u, H)*exp(-u*(4-B)))*jv(0, u*A)


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


if __name__ == "__main__":
    print(L1(1.0, 1.0, 1.0))