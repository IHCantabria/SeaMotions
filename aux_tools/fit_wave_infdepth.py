
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
import matplotlib.pyplot as plt
from numpy import array, exp, log, log10, pi, polyval, ndarray, sqrt, zeros
from numpy import abs as np_abs
from scipy.special import eval_legendre, jv, legendre, roots_chebyt, yv

# Import local modules
from base_integrals import fxy, fxy_dx
from fit_cheby import FitProperties, fit_integral_2d


def factorial(n: int)->int:
    if n < 2:
        return 1
    else:
        f = 1
        for i in range(2, n+1):
            f *= i
        
        return f


def fit_residual_region_test()->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-12
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11, fit_props, "region_test", show_figs=True)


def fit_residual_region_11(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-12, 1e-3, 3.0])
    y = array([1e-12, 1e-3, 4.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_11, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_11", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_11A_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0])
    y = array([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-1, 4.0])
    # x = array([1e-8, 1e-5])
    # y = array([1e-8, 1e-4])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = True
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_11A_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_11_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_11B_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1.0, 3.0])
    y = array([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-1, 4.0])
    # x = array([1e-8, 1e-5])
    # y = array([1e-8, 1e-4])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_11B_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_11_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_12(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-12, 1e-3, 3.0])
    y = array([4.0, 12.0, 50.0, 200.0, 1000.0, 2000.0, 10000, 20000.0, 100000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_12, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_12", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_12_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3.0])
    y = array([4.0, 12.0, 50.0, 200.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_12_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_12_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_21(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 12.0, 50.0, 200.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    y = array([1e-12, 1e-3, 4.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_21, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_21", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_21_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 12.0, 50.0, 200.0, 500.0])
    y = array([1e-12, 1e-3, 4.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_21_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_21_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_22(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    y = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_22, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_22", log_scale=False, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_22_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 6.0, 12.0, 30.0, 50.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    y = array([3.0, 6.0, 12.0, 30.0, 50.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = 1e-8
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on Region: ", region_name)
            results.append(fit_integral_2d(residual_region_22_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_22_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def residual_region_0(X: float, Y: float)->float:
    return X**2.0+Y**2.0


def residual_region_11(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        -2*jv(0, X)*log(R+Y)
        -(pi*yv(0, X)-2*jv(0, X)*log(X))
        )*exp(-Y)
    
    return (fxy(X, Y)-c0)*exp(Y)/R


def residual_region_11A_dx(X: float, Y: float)->float:
    return log10(fxy_dx(X, Y, only_int=True))


def residual_region_11B_dx(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        +2*jv(1, X)*log(R+Y)
        -2*jv(0, X)*X/R/(R+Y)
        +pi*yv(1, X)
        -2*jv(1, X)*log(X)
        +2*jv(0, X)/X
        +X*residual_region_11(X, Y)/R
        )*exp(-Y)
    
    return (fxy_dx(X, Y)-c0)*exp(Y)/R


def residual_region_12(X: float, Y: float)->float:
    return fxy(X, Y)


def residual_region_12_dx(X: float, Y: float)->float:
    return log10(fxy_dx(X, Y, only_int=True))


def residual_region_21(X: float, Y: float)->float:
    return fxy(X, Y)+2*pi*exp(-Y)*yv(0, X)


def residual_region_21_dx(X: float, Y: float)->float:
    return fxy_dx(X, Y)-2*pi*exp(-Y)*yv(1, X)


def residual_region_22(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    costh = Y/R
    f = fxy(X, Y) + 2*pi*exp(-Y)*yv(0, X)
    for i in range(3):
        f -= factorial(i)*eval_legendre(i, costh)/R**(i+1)
    return f


def residual_region_22_dx(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    costh = Y/R
    costhd = -X*Y/R**3.0
    f = fxy_dx(X, Y) - 2*pi*exp(-Y)*yv(1, X)
    for i in range(3):
        poly = legendre(i)
        polyd = poly.deriv()
        f -= factorial(i)*(
                            polyval(polyd, costh)*costhd/R**(i+1)
                            -
                            eval_legendre(i, costh)*2*(i+1)*X/R**(i+2)
                            )
    return f


def plot_results_summary(x: ndarray, y: ndarray, results: list, region_name: str,
                        log_scale=False, show_figs=False)->None:
    def _plot_heat_map(axi, data, log_scale_i)->None:
        imis = axi.imshow(data, origin="lower", cmap="jet", extent=plot_extent)
        axi.set_xticks(x)
        axi.set_yticks(y)
        plt.colorbar(imis, ax=axi)

        if log_scale_i:
            axi.set_xscale("log")
            axi.set_yscale("log")

    # Arrange results in matrix order to plot as 
    # heat map
    num_x = x.shape[0]-1
    num_y = y.shape[0]-1
    err_max = zeros((num_x, num_y))
    err_mean = zeros((num_x, num_y))
    err_min = zeros((num_x, num_y))
    err_over_thr = zeros((num_x, num_y))
    num_coeffs = zeros((num_x, num_y))
    for i in range(num_x):
        for j in range(num_y):
            err_max[i, j] = results[i*num_y+j][3].max_err
            err_mean[i, j] = results[i*num_y+j][3].mean_err
            err_min[i, j] = results[i*num_y+j][3].min_err
            err_over_thr[i, j] = results[i*num_y+j][3].err_over_thr
            num_coeffs[i, j] = results[i*num_y+j][3].num_coeffs

    # Plot results
    plot_extent = [
                    x.min(),
                    x.max(),
                    y.min(),
                    y.max()
                    ]

    fig = plt.figure()
    fig.suptitle(region_name)
    ax0 = fig.add_subplot(231)
    ax0.set_title("Max. Error")
    ax1 = fig.add_subplot(234)
    ax1.set_title("Min. Error")
    ax2 = fig.add_subplot(232)
    ax2.set_title("Mean. Error")
    ax3 = fig.add_subplot(235)
    ax3.set_title("Values Over Threshold")
    ax4 = fig.add_subplot(233)
    ax4.set_title("Num. Coefficients")

    _plot_heat_map(ax0, err_max.T, log_scale)
    _plot_heat_map(ax1, err_min.T, log_scale)
    _plot_heat_map(ax2, err_mean.T, log_scale)
    _plot_heat_map(ax3, err_over_thr.T, log_scale)
    _plot_heat_map(ax4, num_coeffs.T, log_scale)

    if show_figs:
        plt.show()


if __name__ == "__main__":
    # fit_residual_region_test()

    # fit_residual_region_11_subregion_11(show_figs=True)
    # fit_residual_region_11_subregion_12(show_figs=True)
    # fit_residual_region_11_subregion_21(show_figs=True)
    # fit_residual_region_11_subregion_22(show_figs=True)

    # fit_residual_region_11_subregion_11_dx(show_figs=True)
    # fit_residual_region_11_subregion_12_dx(show_figs=True)
    # fit_residual_region_11_subregion_21_dx(show_figs=True)
    # fit_residual_region_11_subregion_22_dx(show_figs=True)

    # fit_residual_region_12_subregion_11(show_figs=True)
    # fit_residual_region_12_subregion_21(show_figs=True)
    # fit_residual_region_12_subregion_12(show_figs=True)
    # fit_residual_region_12_subregion_22(show_figs=True)
    # fit_residual_region_12_subregion_13(show_figs=True)
    # fit_residual_region_12_subregion_23(show_figs=True)
    # fit_residual_region_12_subregion_14(show_figs=True)
    # fit_residual_region_12_subregion_15(show_figs=True)
    # fit_residual_region_12_subregion_24(show_figs=True)

    # fit_residual_region_12_subregion_11_dx(show_figs=True)
    # fit_residual_region_12_subregion_12_dx(show_figs=True)
    # fit_residual_region_12_subregion_13_dx(show_figs=True)
    # fit_residual_region_12_subregion_14_dx(show_figs=True)
    # fit_residual_region_12_subregion_21_dx(show_figs=True)
    # fit_residual_region_12_subregion_22_dx(show_figs=True)
    # fit_residual_region_12_subregion_23_dx(show_figs=True)
    # fit_residual_region_12_subregion_24_dx(show_figs=True)

    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    # fit_residual_region_21_subregion_21(show_figs=True)
    # fit_residual_region_21_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_31(show_figs=True)
    # fit_residual_region_21_subregion_32(show_figs=True)
    # fit_residual_region_21_subregion_41(show_figs=True)
    # fit_residual_region_21_subregion_42(show_figs=True)

    # fit_residual_region_21_subregion_11_dx(show_figs=True)
    # fit_residual_region_21_subregion_12_dx(show_figs=True)
    # fit_residual_region_21_subregion_21_dx(show_figs=True)
    # fit_residual_region_21_subregion_22_dx(show_figs=True)
    # fit_residual_region_21_subregion_31_dx(show_figs=True)
    # fit_residual_region_21_subregion_32_dx(show_figs=True)
    # fit_residual_region_21_subregion_41_dx(show_figs=True)
    # fit_residual_region_21_subregion_42_dx(show_figs=True)

    # fit_residual_region_11(show_figs=True)
    fit_residual_region_11A_dx(show_figs=False)
    # fit_residual_region_11B_dx(show_figs=True)
    # fit_residual_region_22(show_figs=True)
    # fit_residual_region_12_dx(show_figs=True)
    # fit_residual_region_22_dx(show_figs=True)
    # fit_residual_region_22_subregion_11(show_figs=True)

    # fit_residual_region_12()
    # fit_residual_region_11_subregion_21()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    # fit_residual_region_21_subregion_21(show_figs=True)

    # X = 1e-3
    # Y = 1e-3
    # ri = residual_region_11_dx(X, Y)
    # print("residual_region_11_dx: ", ri)