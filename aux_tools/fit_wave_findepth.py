
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
import matplotlib.pyplot as plt
from numpy import array, cos, cosh, exp, linspace, log, log10, meshgrid, mod, ndarray, pi, sin, sqrt, tan, tanh, zeros
from numpy import abs as np_abs
from scipy.integrate import quad
from scipy.special import expi, jv, kn, roots_chebyt, roots_laguerre, struve, yv

# Import local modules
from hydlib.waves.airy import w2k


def compare_integral_series()->None:
    # Define input parameters
    h = 10.0
    A = 1.0
    H = array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0])
    # H = array([1e-6])
    z = -0.2*h
    zeta = -0.2*h

    gi = zeros((H.shape[0], ))
    gs = zeros((H.shape[0], ))

    # R = A*h
    # K = H[-1]/h
    # w = sqrt(K*9.81)
    # T = 2*pi/w
    # print("G_series: ", G_series(R, z, zeta, T, h).real)

    for i in range(H.shape[0]):
        # Define derivative input parameters
        R = A*h
        K = H[i]/h
        w = sqrt(K*9.81)
        T = 2*pi/w

        gi[i] = G_integral(R, z, zeta, T, h).real
        gs[i] = G_series(R, z, zeta, T, h).real

    fig = plt.figure()
    ax0 = fig.add_subplot(311)
    ax1 = fig.add_subplot(312)
    ax2 = fig.add_subplot(313)
    ax0.plot(log10(H), gi, label=f"gi")
    ax1.plot(log10(H), gs, label=f"gs")
    ax2.plot(log10(H), gi-gs, label=f"gi-gs")
    ax0.set_xlabel("H [m]")
    ax1.set_xlabel("H [m]")
    ax2.set_xlabel("H [m]")

    plt.show()


def compare_L3_defs() -> None:
    # Get Gauss-Laguerre points
    gl_roots, gl_weights = roots_laguerre(20)
    
    # Define function properties
    A = 1.0
    B = 0.5
    H = 10000
    u = linspace(0, 1000, 100000)
    y1 = L3_def(A, B, H, u)
    y2 = L_minus_def(A, B, u)

    # Perform integrations using Gauss-Laguerre
    int_L3 = (gl_weights*L3_def(A, B, H, gl_roots)).sum()
    # int_L3 = L3(A, B, H)[0].real
    int_Lm = (gl_weights*L_minus_def(A, B, gl_roots)).sum()

    # Plot Graphs
    plt.plot(u, y1, label=f"L3 - I: {int_L3}")
    plt.plot(u, y2, label=f"Lm - I: {int_Lm}")
    plt.legend()
    plt.show()


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


def fit_residual_3d(f_residual: Callable,
                fit_props: FitProperties,
                region_name: str,
                show_figs = False
                )->None:
    # Define parametric space limits
    num_cross_sections = 5
    x_max = fit_props.x_max
    x_min = fit_props.x_min
    y_max = fit_props.y_max
    y_min = fit_props.y_min
    cheby_order = fit_props.cheby_order
    cheby_tol = fit_props.cheby_tol

    # Define parametric space for function fit
    num_x_fit = 150
    num_y_fit = 150
    x_fit_poly = roots_chebyt(num_x_fit)[0]
    y_fit_poly = roots_chebyt(num_y_fit)[0]
    x_fit = fit_props.x_map_lin(x_fit_poly)
    y_fit = fit_props.y_map_lin(y_fit_poly)
    Xfp,Yfp = meshgrid(x_fit_poly, y_fit_poly, indexing="ij")
    Xf,Yf = meshgrid(x_fit, y_fit, indexing="ij")

    # Calculate residual function over the fit points
    Ff = zeros((num_x_fit, num_y_fit))
    for i in range(num_x_fit):
        for j in range(num_y_fit):
            Ff[i, j] = f_residual(x_fit[i], y_fit[j])

    # Calculate polynomial fit coefficients
    C = fit_chebyshev_2d(Xfp.ravel(), Yfp.ravel(), Ff.ravel(), cheby_order, fit_props.x_map, fit_props.y_map)
    
    # Filter coefficients to the precision required
    n_cheby = linspace(0, cheby_order-1, cheby_order, dtype=int)
    NCX,NCY = meshgrid(n_cheby, n_cheby, indexing="ij")
    pos = np_abs(C) > cheby_tol
    C_filter = C[pos]
    NCX_filter = (NCX.ravel())[pos]
    NCY_filter = (NCY.ravel())[pos]

    # Define parametric space
    num_x = 150
    num_y = 100
    xe = linspace(-1.0, 1.0, num_x)
    ye = linspace(-1.0, 1.0, num_y)
    x = fit_props.x_map_lin(xe)
    y = fit_props.y_map_lin(ye)
    Xe,Ye = meshgrid(xe, ye, indexing="ij")
    X,Y = meshgrid(fit_props.x_map_fcn(x), fit_props.y_map_fcn(y), indexing="ij")
    X0 = x_min
    Y0 = y_min

    # Calculate residual function over the evaluation points
    Fr = zeros((num_x, num_y))
    for i in range(num_x):
        for j in range(num_y):
            Fr[i, j] = f_residual(x[i], y[j])
    
    Ffe = eval_chebyshev_filter_3d(Xe.ravel(), Ye.ravel(), fit_props.x_map, fit_props.y_map, C_filter, NCX_filter, NCY_filter)
    Ffe = Ffe.reshape(num_x, num_y)
    
    y0 = linspace(fit_props.y_min, fit_props.y_max, num_cross_sections)
    FY0 = zeros((num_cross_sections, num_x))
    for i in range(num_cross_sections):
        for j in range(num_x):
            FY0[i, j] = f_residual(x[j], y0[i])

    x0 = linspace(fit_props.y_min, fit_props.y_max, num_cross_sections)
    FX0 = zeros((num_cross_sections, num_y))
    for i in range(num_cross_sections):
        for j in range(num_y):
            FX0[i, j] = f_residual(x0[i], y[j])

    # Print error statistics
    err_tol = 1e-6
    err_abs = np_abs(Ffe-Fr)
    pos = err_abs > err_tol
    err_over_thr = pos.sum()/err_abs.shape[0]/err_abs.shape[1]*100
    print("Error Statistics:")
    print(f" -> Maximum: {err_abs.max()}")
    print(f" -> Minimum: {err_abs.min()}")
    print(f" -> Mean: {err_abs.mean()}")
    print(f" -> Values over threshold [{err_tol:0.1E}]: {err_over_thr} %")

    return C_filter, NCX_filter, NCY_filter


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
    v[1] = np_abs(z+zeta+2*h)
    v[2] = np_abs(z+zeta)
    v[3] = np_abs(z-zeta+2*h)
    v[4] = np_abs(zeta-z+2*h)
    v[5] = np_abs(z+zeta+4*h)
    B = v/h
    r = sqrt(R**2.0+v**2.0)

    # Include radius sumation to green function
    G = (1/r).sum()

    k0 = w2k(w, h, method="bisection")
    
    print(f"R: {R} - k0: {k0} - K: {K} - H: {H} - B2: {B[1]}")
    if B[1] <= 1:
        G += (1/h)*(G1(A, B[0], H)+G1(A, B[1], H))
    else:
        G += (1/h)*(G1(A, B[0], H)+G2(A, B[1], H)) + K*fxy(X, K*v[2])

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

    # print("T: ", T)
    # print("w: ", w)
    # print("h: ", h)
    # print("ki_root: ", ki_root)
    # plt.plot(ki_root*h, tan(ki_root*h))
    # plt.plot(ki_root*h, -K/ki_root)
    # plt.show()

    # Caculate Green function real part
    # print(f"R: {R} - k0: {k0} - K: {K} - H: {H} - B2: {B[1]}")
    # print(f"y0: {yv(0, k0*R)} - j0: {jv(0,k0*R)}")
    # print("cosh: ", cosh(k0*(z+h))*cosh(k0*(zeta+h)))
    # a = -2*pi*(k0**2.0-K**2.0)*(cosh(k0*(z+h))*cosh(k0*(zeta+h)))*(yv(0, k0*R)+jv(0,k0*R)*1j)
    a = -2*pi*(k0**2.0)*(exp(-k0*v[2:]).sum())*(yv(0, k0*R)+jv(0,k0*R)*1j)
    b = ((k0**2.0-K**2.0)*h+K)*(1+exp(-2*k0*h))**2.0
    # print("a: ", a)
    # print("b: ", b)
    G = a/b

    Gs = 0.0
    for i in range(nk):
        Cn = (ki_root[i]**2.0+K**2.0)/((ki_root[i]**2.0+K**2.0)*h-K)
        a = Cn*cos(ki_root[i]*(z+h))*cos(ki_root[i]*(zeta+h))*kn(0, ki_root[i]*R)
        # print(f"Iter: {i} - term_G: {4*a/b} - kn(n): {ki_root[i]} - R: {R} - Cn: {Cn} - K0: {kn(0, ki_root[i]*R)}")
        gi = 4*a
        Gs += gi

        if np_abs(gi) < 1e-7:
            break
    G += Gs

    return G


def G1(A: float, B: float, H: float) -> float:
    if H <= 1:
        print("L1")
        print(L1(A, B, H)[0])
        print("L2")
        print(L2(H)[0])
        int_value = L1(A, B, H)[0] + L2(H)[0]
    else:
        int_value = (
            L3(A, B, H)[0]
            -2/sqrt(A**2.0+(2+B)**2.0)
            -2/sqrt(A**2.0+(2-B)**2.0)
            )
    return int_value


def G2(A: float, B: float, H: float) -> float:
    if H <= 1:
        print("M1")
        print(M1(A, B, H)[0])
        print("M2")
        print(M2(H)[0])
        int_value = M1(A, B, H)[0]+M2(H)[0]
    else:
        int_value = (
            M3(A, B, H)[0]
            -2/sqrt(A**2.0+(2+B)**2.0)
        )

    return int_value


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
    if H <= 1e-5:
        int_value = list(L1_H0_endo(A, B))
        int_value[0] += complex(0.0, 0.0)
    else:
        l1_def_dummy = lambda u: L1_def(A, B, H, u)
        int_value = complex_quadrature_line(l1_def_dummy, way_points)
    print(f"--> L1 - A: {A} - B: {B} - C: {H} - Int_value: {int_value[0]}")

    return int_value


def L1_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return (fuh(u, H)-1)*((exp(-u*(2+B))+exp(-u*(2-B)))*jv(0, u*A)-2*exp(-2*u))


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
    if H<=1e-3:
        int_value = (complex(-log(H)/2-log(2)-1, 0.0), 0.0, 0.0)
    else:
        f_def_dummy = lambda u: L2_def(H, u)
        int_value = complex_quadrature_line(f_def_dummy, way_points)

    return int_value


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


def L3_def(A: float, B: float, H: float, u: ndarray) -> float:
    return (fuh(u, H)+1)*(exp(-u*(2+B))+exp(-u*(2-B)))*jv(0, u*A)


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


def M1_def(A: float, B: float, H: float, u: ndarray) -> float:
    c0 = (fuh(u, H)-1)*(exp(-u*(2+B))*jv(0, u*A)-exp(-3*u))
    c1 = guh(u, H)*(exp(-u*(4-B))*jv(0, u*A)-exp(-3*u))
    return c0+c1


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


def M3_def(A: float, B: float, H: float, u: ndarray) -> ndarray:
    return ((fuh(u, H)+1)*exp(-u*(2+B))+guh(u, H)*exp(-u*(4-B)))*jv(0, u*A)


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


def plot_H_values()->None:
    T = linspace(0.1, 100, int(1e3))
    h = linspace(0.1, 3000, int(1e3))
    TM, HM = meshgrid(T, h)
    w = 2*pi/TM
    K = w**2.0/9.81

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cplt = ax.contourf(K*HM)
    plt.colorbar(cplt, ax=ax)
    ax.set_xlabel("T [s]")
    ax.set_ylabel("H [m]")

    plt.show()


def plot_l1()->None:
    # Define parameter matrix
    A = linspace(0, 1, 10)
    B = array([0.0, 0.5, 1.0, 1.5, 1.7])
    # H = array([0.0, 0.5, 1.0, 10.0, 20.0, 50.0])
    H = array([0.0, 0.005, 0.01, 0.05, 0.1, 0.2])
    
    num_a = A.shape[0]
    num_b = B.shape[0]
    num_h = H.shape[0]

    # u = linspace(-25, 25, int(2e3))
    # y = L1_def_H0(0, 0, 0, u)
    # plt.plot(u, y)
    # plt.show()
    # print(L1_H0_endo(A[0], B[0], H[0])[0])
    # print(L1_H0_endo(A[0], B[1])[0])

    # Allocate space for the solution
    l1 = zeros((A.shape[0], B.shape[0], H.shape[0]))
    for i in range(num_h):
        for j in range(num_b):
            for k in range(num_a):
                # print(f"A[{k:d}]: {A[k]} - B[{j:d}]: {B[j]} - H[{i:d}]: {H[i]}")
                l1[k, j, i] = L1(A[k], B[j], H[i])[0].real
                # print(f"L1: {l1[k, j, i]}")
    # for i in range(3, 6):
    #     for j in range(num_b):
    #         for k in range(num_a):
    #             # print(f"A[{k:d}]: {A[k]} - B[{j:d}]: {B[j]} - H[{i:d}]: {H[i]}")
    #             l1[k, j, i] = L1_H0_endo(A[k], B[j])[0].real
                # print(f"L1: {l1[k, j, i]}")

    # Plot values
    fig = plt.figure()
    ax = []
    for i in range(num_h):
        axi = fig.add_subplot(2,3,i+1)
        axi.set_title(f"H: {H[i]}")
        ax.append(axi)

    for i in range(num_h):
        for j in range(num_b):
            ax[i].plot(A, l1[:, j, i], label=f"B={B[j]:0.1f}")
        ax[i].legend()

    plt.show()


def plot_integral_parametric_L1() -> None:
    # Define function to plot
    def plot_fun(ax, ft: Callable, X: ndarray, Y: ndarray, Z: ndarray, Ylabel: str, logx=False) -> None:
        x_data = zeros((X.shape[0], ))
        for zi in range(Z.shape[0]):
            for yi in range(Y.shape[0]):
                for xi in range(X.shape[0]):
                    x_data[xi] = ft(X[xi], Y[yi], Z[zi])[0].real
                ax[zi].plot(log10(X), x_data, label=Ylabel+f": {Y[yi]:0.1f}")
            ax[zi].legend()
            if logx:
                ax[zi].set_xscale("log")
        
    # Define parametric space
    
    # Draw figures For X = A
    # num_a = 10
    # num_b = 6
    # num_h = 9
    # A = linspace(0.0, 1.0, num_a)
    # B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])

    # fig = plt.figure()
    # ax = []
    # for i in range(num_h):
    #     ax.append(fig.add_subplot(3,3,i+1))
    #     ax[i].set_title(f"H: {H[i]}")
    # plot_fun(ax, lambda X, Y, Z: L1(X, Y, Z), A, B, H, "B")


    # Draw figures for X = B
    # num_a = 6
    # num_b = 10
    # num_h = 9
    # A = linspace(0.0, 1.0, num_a)
    # B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])

    # fig = plt.figure()
    # ax = []
    # for i in range(num_h):
    #     ax.append(fig.add_subplot(3,3,i+1))
    #     ax[i].set_title(f"H: {H[i]}")
    # plot_fun(ax, lambda X, Y, Z: L1(Y, X, Z), B, A, H, "A")

    # Draw figures for X = H
    num_a = 9
    num_b = 6
    num_h = 9
    A = linspace(0.0, 1.0, num_a)
    B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])
    H = 10**linspace(-6, -2, 30)

    fig = plt.figure()
    ax = []
    for i in range(num_a):
        ax.append(fig.add_subplot(3,3,i+1))
        ax[i].set_title(f"A: {A[i]}")
    plot_fun(ax, lambda X, Y, Z: L1(Z, Y, X), H, B, A, "B")

    plt.show()


def plot_integral_parametric_L1_derA() -> None:
    # Define function to plot
    def plot_fun(ax, ft: Callable, X: ndarray, Y: ndarray, Z: ndarray, der_chn: int, Ylabel: str, logx=False) -> None:
        eps = 1e-6
        x_data = zeros((X.shape[0], ))
        for zi in range(Z.shape[0]):
            for yi in range(Y.shape[0]):
                for xi in range(X.shape[0]):
                    if der_chn == 0:
                        x_data[xi] = (
                                        ft(X[xi]+eps, Y[yi], Z[zi])[0].real
                                        - ft(X[xi]-eps, Y[yi], Z[zi])[0].real
                                    )/2/eps
                    elif der_chn == 1:
                        x_data[xi] = (
                                        ft(X[xi], Y[yi]+eps, Z[zi])[0].real
                                        - ft(X[xi], Y[yi]-eps, Z[zi])[0].real
                                    )/2/eps
                    elif der_chn == 2:
                        x_data[xi] = (
                                        ft(X[xi], Y[yi], Z[zi]+eps)[0].real
                                        - ft(X[xi], Y[yi], Z[zi]-eps)[0].real
                                    )/2/eps
                ax[zi].plot(X, x_data, label=Ylabel+f": {Y[yi]:0.1f}")
            ax[zi].legend()
            if logx:
                ax[zi].set_xscale("log")
        
    # Define parametric space
    
    # Draw figures For X = A
    # num_a = 10
    # num_b = 6
    # num_h = 9
    # A = linspace(0.0, 1.0, num_a)
    # B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])

    # fig = plt.figure()
    # ax = []
    # for i in range(num_h):
    #     ax.append(fig.add_subplot(3,3,i+1))
    #     ax[i].set_title(f"H: {H[i]}")
    # plot_fun(ax, lambda X, Y, Z: L1(X, Y, Z), A, B, H, 0, "B")


    # Draw figures for X = B
    # num_a = 6
    # num_b = 10
    # num_h = 9
    # A = linspace(0.0, 1.0, num_a)
    # B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])

    # fig = plt.figure()
    # ax = []
    # for i in range(num_h):
    #     ax.append(fig.add_subplot(3,3,i+1))
    #     ax[i].set_title(f"H: {H[i]}")
    # plot_fun(ax, lambda X, Y, Z: L1(Y, X, Z), B, A, H, 1, "A")

    # Draw figures for X = H
    num_a = 9
    num_b = 6
    num_h = 9
    A = linspace(0.0, 1.0, num_a)
    B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])
    H = 10**linspace(-1, 1, 30)

    fig = plt.figure()
    ax = []
    for i in range(num_a):
        ax.append(fig.add_subplot(3,3,i+1))
        ax[i].set_title(f"A: {A[i]}")
    plot_fun(ax, lambda X, Y, Z: L1(Z, Y, X), H, B, A, 2, "B", logx=True)

    plt.show()


def plot_integral_parametric_L1_derB() -> None:
    # Define function to plot
    def plot_fun(ax, ft: Callable, X: ndarray, Y: ndarray, Z: ndarray, der_chn: int, Ylabel: str, logx=False) -> None:
        eps = 1e-6
        x_data = zeros((X.shape[0], ))
        for zi in range(Z.shape[0]):
            for yi in range(Y.shape[0]):
                for xi in range(X.shape[0]):
                    if der_chn == 0:
                        x_data[xi] = (
                                        ft(X[xi]+eps, Y[yi], Z[zi])[0].real
                                        - ft(X[xi]-eps, Y[yi], Z[zi])[0].real
                                    )/2/eps
                    elif der_chn == 1:
                        x_data[xi] = (
                                        ft(X[xi], Y[yi]+eps, Z[zi])[0].real
                                        - ft(X[xi], Y[yi]-eps, Z[zi])[0].real
                                    )/2/eps
                    elif der_chn == 2:
                        x_data[xi] = (
                                        ft(X[xi], Y[yi], Z[zi]+eps)[0].real
                                        - ft(X[xi], Y[yi], Z[zi]-eps)[0].real
                                    )/2/eps
                ax[zi].plot(X, x_data, label=Ylabel+f": {Y[yi]:0.1f}")
            ax[zi].legend()
            if logx:
                ax[zi].set_xscale("log")
        
    # Define parametric space
    
    # Draw figures For X = A
    # num_a = 10
    # num_b = 6
    # num_h = 9
    # A = linspace(0.0, 1.0, num_a)
    # B = linspace(0.0, 1.0, num_b)
    # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])

    # fig = plt.figure()
    # ax = []
    # for i in range(num_h):
    #     ax.append(fig.add_subplot(3,3,i+1))
    #     ax[i].set_title(f"H: {H[i]}")
    # plot_fun(ax, lambda X, Y, Z: L1(X, Y, Z), A, B, H, 1, "B")


    # Draw figures for X = B
    num_a = 6
    num_b = 10
    num_h = 9
    A = linspace(0.0, 1.0, num_a)
    B = linspace(0.0, 1.0, num_b)
    H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])

    fig = plt.figure()
    ax = []
    for i in range(num_h):
        ax.append(fig.add_subplot(3,3,i+1))
        ax[i].set_title(f"H: {H[i]}")
    plot_fun(ax, lambda X, Y, Z: L1(Y, X, Z), B, A, H, 1, "A")

    # Draw figures for X = H
    # num_a = 9
    # num_b = 6
    # num_h = 9
    # A = linspace(0.0, 1.0, num_a)
    # B = linspace(0.0, 1.0, num_b)
    # # H = array([1e-6, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3])
    # H = 10**linspace(-1, 1, 30)

    # fig = plt.figure()
    # ax = []
    # for i in range(num_a):
    #     ax.append(fig.add_subplot(3,3,i+1))
    #     ax[i].set_title(f"A: {A[i]}")
    # plot_fun(ax, lambda X, Y, Z: L1(Z, Y, X), H, B, A, 2, "B", logx=True)

    plt.show()


def test_complex_quadrature()->None:
    waypoints = array([
                        -1e6+0.0j,
                        -1e3+0.0j,
                        -1e2+0.0j,
                        -1.0+0.0j,
                         0.0+1.0j,
                         1.0+0.0j,
                         1e2+0.0j,
                         1e3+0.0j,
                         1e6+0.0j
                        ])
    print(complex_quadrature_line(lambda x: 1/x, waypoints))


def test_gauss_laguerre()->None:
    # Get Gauss-Laguerre points
    gl_roots, gl_weights = roots_laguerre(20)

    # Define function to integrate
    fint = lambda x: exp(-x)/x


def test_pole_values()->None:
    # Define problem parameters
    g = 9.81
    h = 10
    A = 1
    B = 1
    H = 1

    # Calculate problem derivative parameters
    nu = H/h
    w = sqrt(nu*g)
    k0 = w2k(w, h, method="bisection")
    u0 = k0*h

    # Calculate L1 using numerical rules
    l1_num, err_real, err_imag = L1(A, B, H)

    # Calculate L2 using numerical rules
    l2_num, err_real, err_imag = L2(H)

    # Calculate L3 using numerical rules
    l3_num, err_real, err_imag = L3(A, B, H)

    # Calculate M1 using numerical rules
    m1_num, err_real, err_imag = M1(A, B, H)

    # Calculate M2 using numerical rules
    m2_num, err_real, err_imag = M2(H)

    # Calculate M3 using numerical rules
    m3_num, err_real, err_imag = M3(A, B, H)

    # Calculate L1 pole
    l1_pole = L1_pole(A, B, H, u0)

    # Calculate L2 pole
    l2_pole = L2_pole(H, u0)

    # Calculate L3 pole
    l3_pole = L3_pole(A, B, H, u0)

    # Calculate M1 pole
    m1_pole = M1_pole(A, B, H, u0)

    # Calculate M2 pole
    m2_pole = M2_pole(H, u0)

    # Calculate M3 pole
    m3_pole = M3_pole(A, B, H, u0)

    # Plot results
    print(f"L1_num: {l1_num} - L1_pole: {l1_pole} - Poles Diff: {l1_num.imag+l1_pole}")
    print(f"L2_num: {l2_num} - L2_pole: {l2_pole} - Poles Diff: {l2_num.imag+l2_pole}")
    print(f"L3_num: {l3_num} - L3_pole: {l3_pole} - Poles Diff: {l3_num.imag+l3_pole}")
    print(f"M1_num: {m1_num} - M1_pole: {m1_pole} - Poles Diff: {m1_num.imag+m1_pole}")
    print(f"M2_num: {m2_num} - M2_pole: {m2_pole} - Poles Diff: {m2_num.imag+m2_pole}")
    print(f"M3_num: {m3_num} - M3_pole: {m3_pole} - Poles Diff: {m3_num.imag+m3_pole}")


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
    # test_pole_values()
    # plot_l1()
    # compare_integral_series()
    # plot_H_values()
    # compare_L3_defs()
    # print(M1(1.0, 1.8, 1e3))
    # print(M2(1e-6))
    # plot_integral_parametric_L1()
    # plot_integral_parametric_L1_derA()
    # plot_integral_parametric_L1_derB()
    # print("L1: ", L1(0.0, 2.0, 0.01))
    # u = linspace(0, 30, 1000)
    # plt.plot(u, L1_def(0.0, 2.0, 0.01, u))
    # plt.show()
    # print(L1(0.725, 0.5, 1.0))
    # from numpy import sin, tanh
    # x = linspace(-5, 5, 100)
    # y = (tanh(x)+1)/2+sin(x*pi/5)/10
    # plt.plot(x, y)
    # plt.show()
    test_complex_quadrature()