
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import exp, flip, linspace, log, log10, meshgrid, ndarray, pi, sqrt, sign, zeros
from numpy import abs as np_abs
from scipy.special import jv, roots_chebyt, yv

# Import local modules
from base_integrals import fxy
from fit_cheby import eval_chebyshev, eval_chebyshev_filter, fit_chebyshev


class FitProperties:

    def __init__(self) -> None:
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


def fit_residual_region_test()->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 3.0
    fit_props.x_min = 0.0
    fit_props.y_max = 4.0
    fit_props.y_min = 0.0
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_0, fit_props, "region_test")


def fit_residual_region_11_subregion_11(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max = 1e-4
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = True
    fit_props.y_max = 1e-4
    fit_props.y_min = 1e-12
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_11, fit_props, "region_11", show_figs=show_figs)


def fit_residual_region_11_subregion_12(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max = 1e-4
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-4
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_11, fit_props, "region_12", show_figs=show_figs)


def fit_residual_region_11_subregion_21(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-4
    fit_props.y_log_scale = True
    fit_props.y_max = 1e-4
    fit_props.y_min = 1e-12
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_11, fit_props, "region_21", show_figs=show_figs)


def fit_residual_region_11_subregion_22(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 3
    fit_props.x_min = 1e-4
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-4
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_11, fit_props, "region_22", show_figs=show_figs)


def fit_residual_region_21_subregion_11(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = True
    fit_props.y_max = 1e-1
    fit_props.y_min = 1e-12
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_12, fit_props, "region_12_subregion_11", show_figs=show_figs)


def fit_residual_region_21_subregion_12(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = True
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-1
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_12, fit_props, "region_12_subregion_21", show_figs=show_figs)


def fit_residual_region_21_subregion_21(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 150.0
    fit_props.x_min = 50.0
    fit_props.y_log_scale = True
    fit_props.y_max = 1e-1
    fit_props.y_min = 1e-12
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_12, fit_props, "region_12_subregion_11", show_figs=show_figs)


def fit_residual_region_12(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = True
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-4
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_residual(residual_region_12, fit_props, "region_12", show_figs=show_figs)


def residual_region_0(X: float, Y: float)->float:
    return X**2.0+Y**2.0


def residual_region_11(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        -2*jv(0, X)*log(R+Y)
        -(pi*yv(0, X)-2*jv(0, X)*log(X))
        )*exp(-Y)
    
    return (fxy(X, Y)-c0)*exp(Y)/R


def residual_region_12(X: float, Y: float)->float:
    return fxy(X, Y)+2*pi*exp(-Y)*yv(0, X)


def fit_residual(f_residual: Callable,
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
    C = fit_chebyshev(Xfp.ravel(), Yfp.ravel(), Ff.ravel(), cheby_order, fit_props.x_map, fit_props.y_map)
    
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
    
    Ffe = eval_chebyshev_filter(Xe.ravel(), Ye.ravel(), fit_props.x_map, fit_props.y_map, C_filter, NCX_filter, NCY_filter)
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


if __name__ == "__main__":
    # fit_residual_region_test()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_12()
    # fit_residual_region_11_subregion_21()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    fit_residual_region_21_subregion_21(show_figs=True)