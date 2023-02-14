
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
from numpy import array, exp, log,pi, sqrt
from numpy import abs as np_abs
from scipy.special import eval_legendre, jv, roots_chebyt, yv

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


def fit_residual_region_11_subregion_11(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11, fit_props, "region_11_11", show_figs=show_figs)


def fit_residual_region_11_subregion_11_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-4
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11_dx, fit_props, "region_11_11_dx", show_figs=show_figs)


def fit_residual_region_11_subregion_12(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-4
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-4
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11, fit_props, "region_11_12", show_figs=show_figs)


def fit_residual_region_11_subregion_12_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-4
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11_dx, fit_props, "region_11_12_dx", show_figs=show_figs)


def fit_residual_region_11_subregion_21(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-4
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-4
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11, fit_props, "region_11_21", show_figs=show_figs)


def fit_residual_region_11_subregion_21_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11_dx, fit_props, "region_11_21_dx", show_figs=show_figs)


def fit_residual_region_11_subregion_22(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_max = 3
    fit_props.x_min = 1e-4
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-4
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11, fit_props, "region_11_22", show_figs=show_figs)


def fit_residual_region_11_subregion_22_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11_dx, fit_props, "region_11_22_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_11(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 12.0
    fit_props.y_min = 4.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_11", show_figs=show_figs)


def fit_residual_region_12_subregion_11_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-1
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 12.0
    fit_props.y_min = 4.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12_dx, fit_props, "region_12_11_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_12(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 50.0
    fit_props.y_min = 12.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_12", show_figs=show_figs)


def fit_residual_region_12_subregion_12_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-1
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 50.0
    fit_props.y_min = 12.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12_dx, fit_props, "region_12_12_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_13(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 200.0
    fit_props.y_min = 50.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_13", show_figs=show_figs)


def fit_residual_region_12_subregion_13_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-1
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 200.0
    fit_props.y_min = 50.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12_dx, fit_props, "region_12_13_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_14(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-3
    fit_props.x_min = 1e-12
    fit_props.y_log_scale = False
    fit_props.y_max = 1000.0
    fit_props.y_min = 200.0
    fit_props.cheby_tol = 1e-14
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_14", show_figs=show_figs)


def fit_residual_region_12_subregion_14_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-1
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 1000.0
    fit_props.y_min = 200.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    # fit_integral_2d(residual_region_12_dx, fit_props, "region_12_14_dx", show_figs=show_figs)
    print("Return zero coefficients by hand.")


def fit_residual_region_12_subregion_21(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 12.0
    fit_props.y_min = 4.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_21", show_figs=show_figs)


def fit_residual_region_12_subregion_21_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-1
    fit_props.y_log_scale = False
    fit_props.y_max = 12.0
    fit_props.y_min = 4.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12_dx, fit_props, "region_12_21_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_22(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 50.0
    fit_props.y_min = 12.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_22", show_figs=show_figs)


def fit_residual_region_12_subregion_22_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-1
    fit_props.y_log_scale = False
    fit_props.y_max = 50.0
    fit_props.y_min = 12.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12_dx, fit_props, "region_12_22_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_23(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 200.0
    fit_props.y_min = 50.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_23", show_figs=show_figs)


def fit_residual_region_12_subregion_23_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-1
    fit_props.y_log_scale = False
    fit_props.y_max = 200.0
    fit_props.y_min = 50.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12_dx, fit_props, "region_12_23_dx", show_figs=show_figs)


def fit_residual_region_12_subregion_24(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 1000.0
    fit_props.y_min = 200.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_12, fit_props, "region_12_24", show_figs=show_figs)


def fit_residual_region_12_subregion_24_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-1
    fit_props.y_log_scale = False
    fit_props.y_max = 1000.0
    fit_props.y_min = 200.0
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    # fit_integral_2d(residual_region_12_dx, fit_props, "region_12_24_dx", show_figs=show_figs)
    print("Return zero coefficients by hand.")


def fit_residual_region_22(show_figs=False):
    # Define boundary values
    x = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0])
    y = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0])
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"
            fit_residual_region_22_subregion_11(x[i], x[i+1], y[j], y[j+1], region_name, show_figs=True)


def fit_residual_region_22_dx(show_figs=False):
    # Define boundary values
    x = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0])
    y = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0])
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"X: {x[i]:0.1f} - {x[i+1]:0.1f} | Y: {y[j]:0.1f} - {y[j+1]:0.1f}"
            fit_residual_region_22_subregion_11_dx(x[i], x[i+1], y[j], y[j+1], region_name, show_figs=True)


def fit_residual_region_22_subregion_11(x_min: float, x_max: float, y_min: float, y_max: float,
                                        region_name: str, show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = x_max
    fit_props.x_min = x_min
    fit_props.y_log_scale = False
    fit_props.y_max = y_max
    fit_props.y_min = y_min
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_22, fit_props, region_name, show_figs=show_figs)


def fit_residual_region_22_subregion_11_dx(x_min: float, x_max: float, y_min: float, y_max: float,
                                            region_name: str, show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = x_max
    fit_props.x_min = x_min
    fit_props.y_log_scale = False
    fit_props.y_max = y_max
    fit_props.y_min = y_min
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_22_dx, fit_props, region_name, show_figs=show_figs)


def fit_residual_region_21_subregion_11(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_11", show_figs=show_figs)


def fit_residual_region_21_subregion_11_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_11_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_12(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_12", show_figs=show_figs)


def fit_residual_region_21_subregion_12_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 12.0
    fit_props.x_min = 3.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-4
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_12_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_21(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 50.0
    fit_props.x_min = 12.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_21", show_figs=show_figs)


def fit_residual_region_21_subregion_21_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 50.0
    fit_props.x_min = 12.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_21_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_22(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 50.0
    fit_props.x_min = 12.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_22", show_figs=show_figs)


def fit_residual_region_21_subregion_22_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 50.0
    fit_props.x_min = 12.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_22_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_31(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 200.0
    fit_props.x_min = 50.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_31", show_figs=show_figs)


def fit_residual_region_21_subregion_31_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 200.0
    fit_props.x_min = 50.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_31_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_32_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 200.0
    fit_props.x_min = 50.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_32_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_32(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 200.0
    fit_props.x_min = 50.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_32", show_figs=show_figs)


def fit_residual_region_21_subregion_41(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1000.0
    fit_props.x_min = 200.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_41", show_figs=show_figs)


def fit_residual_region_21_subregion_41_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1000.0
    fit_props.x_min = 200.0
    fit_props.y_log_scale = False
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e-12
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_41_dx", show_figs=show_figs)


def fit_residual_region_21_subregion_42(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1000.0
    fit_props.x_min = 200.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21, fit_props, "region_21_42", show_figs=show_figs)


def fit_residual_region_21_subregion_42_dx(show_figs=False)->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.x_log_scale = False
    fit_props.x_max = 1000.0
    fit_props.x_min = 200.0
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_tol = 1e-8
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_21_dx, fit_props, "region_21_42_dx", show_figs=show_figs)


def residual_region_0(X: float, Y: float)->float:
    return X**2.0+Y**2.0


def residual_region_11(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        -2*jv(0, X)*log(R+Y)
        -(pi*yv(0, X)-2*jv(0, X)*log(X))
        )*exp(-Y)
    
    return (fxy(X, Y)-c0)*exp(Y)/R


def residual_region_11_dx(X: float, Y: float)->float:
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
    return fxy_dx(X, Y)


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
    f = fxy(X, Y) + 2*pi*exp(-Y)*yv(0, X)
    for i in range(3):
        f -= factorial(i)*eval_legendre(i, costh)/R**(i+1)
    return f


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

    fit_residual_region_22(show_figs=True)
    # fit_residual_region_22_subregion_11(show_figs=True)

    # fit_residual_region_12()
    # fit_residual_region_11_subregion_21()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    # fit_residual_region_21_subregion_21(show_figs=True)