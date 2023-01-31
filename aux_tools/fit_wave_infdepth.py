
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import exp, log,pi, sqrt
from numpy import abs as np_abs
from scipy.special import jv, roots_chebyt, yv

# Import local modules
from base_integrals import fxy
from fit_cheby import FitProperties, fit_integral_2d


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
    fit_integral_2d(residual_region_0, fit_props, "region_test")


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
    fit_integral_2d(residual_region_11, fit_props, "region_11", show_figs=show_figs)


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
    fit_integral_2d(residual_region_11, fit_props, "region_12", show_figs=show_figs)


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
    fit_integral_2d(residual_region_11, fit_props, "region_21", show_figs=show_figs)


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
    fit_integral_2d(residual_region_11, fit_props, "region_22", show_figs=show_figs)


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
    fit_integral_2d(residual_region_12, fit_props, "region_12_subregion_11", show_figs=show_figs)


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
    fit_integral_2d(residual_region_12, fit_props, "region_12_subregion_21", show_figs=show_figs)


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
    fit_integral_2d(residual_region_12, fit_props, "region_12_subregion_11", show_figs=show_figs)


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
    fit_integral_2d(residual_region_12, fit_props, "region_12", show_figs=show_figs)


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


if __name__ == "__main__":
    # fit_residual_region_test()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_12()
    # fit_residual_region_11_subregion_21()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    fit_residual_region_21_subregion_21(show_figs=True)