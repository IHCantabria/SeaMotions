
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import linspace, log10, meshgrid, ndarray, zeros
from numpy import abs as np_abs
from numpy.linalg import solve as np_solve
from scipy.special import eval_chebyt, roots_chebyt


class FitProperties:

    def __init__(self) -> None:
        self.num_x = 20
        self.num_x_fit = 20
        self.num_y = 20
        self.num_y_fit = 20
        self.x_log_scale = False
        self.x_max = 0
        self.x_map_fcn = None
        self.x_min = 0.0
        self.y_log_scale = False
        self.y_max = 0.0
        self.y_map_fcn = None
        self.y_min = 0.0
        self.cheby_order = 0
        self.cheby_tol = 0.0

    def x_map(self, x: ndarray)->ndarray:
        return self.x_map_fcn(self.x_map_lin(x))
        # return self.x_map_fcn(x)

    def x_map_lin(self, x: ndarray)->ndarray:
        if self.x_log_scale:
            return 10**((x+1)*(log10(self.x_max)-log10(self.x_min))/2 + log10(self.x_min))
        else:
            return (x+1)*(self.x_max-self.x_min)/2 + self.x_min
    
    def y_map(self, y: ndarray)->ndarray:
        # return self.y_map_fcn(self.y_map_lin(y))
        return self.y_map_fcn(y)
    
    def y_map_lin(self, y: ndarray)->ndarray:
        if self.y_log_scale:
            return 10**((y+1)*(log10(self.y_max)-log10(self.y_min))/2 + log10(self.y_min))
        else:
            return (y+1)*(self.y_max-self.y_min)/2 + self.y_min


def fit_integral_2d(f_residual: Callable,
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
    num_x_fit = fit_props.num_x_fit
    num_y_fit = fit_props.num_y_fit
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
    C = fit_chebyshev_2d(Xfp.ravel(), Yfp.ravel(), Ff.ravel(), cheby_order)
    
    # Filter coefficients to the precision required
    n_cheby = linspace(0, cheby_order-1, cheby_order, dtype=int)
    NCX,NCY = meshgrid(n_cheby, n_cheby, indexing="ij")
    pos = np_abs(C) > cheby_tol
    C_filter = C[pos]
    NCX_filter = (NCX.ravel())[pos]
    NCY_filter = (NCY.ravel())[pos]

    # Define parametric space
    num_x = fit_props.num_x
    num_y = fit_props.num_y
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
    
    Ffe = eval_chebyshev_2d_filter(Xe.ravel(), Ye.ravel(), C_filter, NCX_filter, NCY_filter)
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

    # Plot residual function
    fig = plt.figure()
    fig.suptitle(region_name + f" - Total Coeffs: {C_filter.shape[0]:d}")
    ax0 = fig.add_subplot(231)
    ax1 = fig.add_subplot(232, projection='3d')
    ax4 = fig.add_subplot(233)
    ax2 = fig.add_subplot(234)
    ax3 = fig.add_subplot(235)
    ax5 = fig.add_subplot(236)

    cfx = ax0.contourf(X, Y, Fr)
    ax0.set_xlabel("X")
    ax0.set_ylabel("Y")
    if fit_props.x_log_scale:
        ax0.set_xscale("log")
    if fit_props.y_log_scale:
        ax0.set_yscale("log")
    plt.colorbar(cfx, ax=ax0)

    psf = ax1.plot_surface(X, Y, Fr, cmap="jet")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    plt.colorbar(psf, ax=ax1)

    ax2.set_title(f"Y = {y_min:0.1E}")
    for i in range(num_cross_sections):
        ax2.plot(fit_props.x_map_fcn(x), FY0[i, :], label=f"Y = {y0[i]:0.2E}")
    ax2.legend()
    if fit_props.x_log_scale:
        ax2.set_xscale("log")
    ax2.set_xlabel("X")

    ax3.set_title(f"X = {x_min:0.1E}")
    for i in range(num_cross_sections):
        ax3.plot(fit_props.y_map_fcn(y), FX0[i, :], label=f"X = {x0[i]:0.2E}")
    ax3.legend()
    if fit_props.y_log_scale:
        ax3.set_xscale("log")
    ax3.set_xlabel("Y")

    ax4.set_title("log10(|Ffe-Fr|)")
    cntf = ax4.contourf(X, Y, log10(err_abs), cmap="jet")
    ax4.set_xlabel("X")
    ax4.set_ylabel("Y")
    plt.colorbar(cntf, ax=ax4)

    ax5.set_title("log10(|Ffe-Fr|/|Fr|)")
    cntf = ax5.contourf(X, Y, log10(np_abs(Ffe-Fr)/np_abs(Fr)))
    ax5.set_xlabel("X")
    ax5.set_ylabel("Y")
    plt.colorbar(cntf, ax=ax5)

    cheby_order_ticks = linspace(0, cheby_order-1, cheby_order)
    fig = plt.figure()
    fig.suptitle("Fit coefficient matrix")
    ax = fig.add_subplot(111)
    imc = ax.imshow(log10(C.reshape(cheby_order, cheby_order).T), origin="lower")
    ax.set_xlabel("X")
    ax.set_xticks(cheby_order_ticks)
    ax.set_ylabel("Y")
    ax.set_yticks(cheby_order_ticks)
    plt.colorbar(imc, ax=ax)

    if show_figs:
        plt.show()

    return C_filter, NCX_filter, NCY_filter


def eval_chebyshev_1d(x: ndarray, n: int, c: ndarray)->ndarray:
    sol = 0.0
    for i in range(n):
        sol += c[i]*eval_chebyt(i, x)

    return sol


def eval_chebyshev_2d(x: ndarray, y: ndarray, n: int, c: ndarray)->ndarray:
    sol = 0.0
    for i in range(n):
        for j in range(n):
            sol += c[i*n+j]*eval_chebyt(i, x)*eval_chebyt(j, y)

    return sol


def eval_chebyshev_2d_filter(x: ndarray, y: ndarray, c: ndarray, ncx: ndarray, ncy: ndarray)->ndarray:
    sol = 0.0
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)*eval_chebyt(ncy[i], y)

    return sol


def eval_chebyshev_3d_filter(x: ndarray, y: ndarray, z: ndarray, c: ndarray, ncx: ndarray, ncy: ndarray,
                            ncz: ndarray)->ndarray:
    sol = 0.0
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)*eval_chebyt(ncy[i], y)*eval_chebyt(ncz[i], z)

    return sol


def fit_chebyshev_1d(x: ndarray, f: ndarray, n: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, n))
    for i in range(n):
            A[:, i] = eval_chebyt(i, x)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C


def fit_chebyshev_2d(x: ndarray, y: ndarray, f: ndarray, n: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, n*n))
    for i in range(n):
        for j in range(n):
            A[:, i*n+j] = eval_chebyt(i, x)*eval_chebyt(j, y)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C


def fit_chebyshev_3d(x: ndarray, y: ndarray, z: ndarray, f: ndarray, nx: int, ny: int, nz: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, nx*ny*nz))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                A[:, i*(ny*nz)+j*nz+k] = eval_chebyt(i, x)*eval_chebyt(j, y)*eval_chebyt(k, z)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C