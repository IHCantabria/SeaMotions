
# Import general usage libraries
import h5py
import os

# Import general usage scientific libraries
from numpy import argsort, array, concatenate, log, ndarray
from numpy import abs as np_abs

# Import local modules
from base_integrals import L1, L1_dA, L1_dB, L2, L3, L3_dA, L3_dB, M1, M1_dA, M1_dB, M2, M3, M3_dA, M3_dB
from fit_cheby import fit_integral_1d, fit_integral_2d, fit_integral_3d, FitProperties


def data_to_dict(
                fit_props: FitProperties,
                cheby_coeffs: ndarray,
                ncx: ndarray,
                ncy: ndarray,
                ncz: ndarray
                ):
    # Generate dict with the results data
    data = {
            "fit_props": fit_props,
            "cheby_coeffs": cheby_coeffs,
            "ncx": ncx,
            "ncy": ncy,
            "ncz": ncz
            }
    
    return data


def fit_L1()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-16, 0.1]
    intervals_data.append(fit_L1_P0())

    # Calculate coefficients for H = [0.1, 1.0]
    intervals_data.append(fit_L1_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L1_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_L1_dA()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-16, 0.1]
    intervals_data.append(fit_L1_dA_P0())

    # Calculate coefficients for H = [0.1, 1.0]
    intervals_data.append(fit_L1_dA_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L1_dA_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_L1_dA_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 0.1
    fit_props.z_min = 1e-16
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L1_dA(x, y, z)[0].real, fit_props, "L1_dA_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L1_dA_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = False
    fit_props.z_max = 1.0
    fit_props.z_min = 0.1
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L1_dA(x, y, z)[0].real, fit_props, "L1_dA_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L1_dB()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-16, 0.1]
    intervals_data.append(fit_L1_dB_P0())

    # Calculate coefficients for H = [0.1, 1.0]
    intervals_data.append(fit_L1_dB_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L1_dB_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_L1_dB_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 0.1
    fit_props.z_min = 1e-16
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L1_dB(x, y, z)[0].real, fit_props, "L1_dB_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L1_dB_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = False
    fit_props.z_max = 1.0
    fit_props.z_min = 0.1
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L1_dB(x, y, z)[0].real, fit_props, "L1_dB_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L1_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 0.1
    fit_props.z_min = 1e-16
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L1(x, y, z)[0].real, fit_props, "L1_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L1_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = False
    fit_props.z_max = 1.0
    fit_props.z_min = 0.1
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L1(x, y, z)[0].real, fit_props, "L1_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L1_Afix_Bfix()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.num_x = 100
    fit_props.num_x_fit = 100
    fit_props.x_log_scale = True
    fit_props.x_max = 3.0
    fit_props.x_min = 1e-12
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.cheby_order = 50
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    A = 0.1
    B = 0.1

    # Launch fit
    fit_integral_1d(lambda z: L1(A, B, z)[0].real, fit_props, f"L1_A: {A:0.2E} - B: {B:0.2E}", show_figs=True)


def fit_L1_Hfix()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    H = 30.0

    # Launch fit
    fit_integral_2d(lambda x, y: L1(x, y, H)[0].real, fit_props, f"L1_H: {H:0.2E}", show_figs=True)


def fit_L2()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-12, 1.0]
    intervals_data.append(fit_L2_P0())

    # Calculate coefficients for H = [1.0, 50.0]
    intervals_data.append(fit_L2_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L2_database.h5")
    write_intervals(database_path, intervals_data, 1)


def fit_L2_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max = 0.001
    fit_props.x_min = 1e-16
    fit_props.cheby_order_x = 50
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx = fit_integral_1d(lambda z: L2(z)[0].real, fit_props, f"L2", show_figs=True)

    return data_to_dict(fit_props, cheby_coeffs, ncx, array([]), array([]))


def fit_L2_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max = 1.0
    fit_props.x_min = 0.001
    fit_props.cheby_order_x = 50
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx = fit_integral_1d(lambda z: L2(z)[0].real, fit_props, f"L2", show_figs=True)

    return data_to_dict(fit_props, cheby_coeffs, ncx, array([]), array([]))


def fit_L3()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1.0, 12.0]
    intervals_data.append(fit_L3_P0())

    # Calculate coefficients for H = [12.0, 50.0]
    intervals_data.append(fit_L3_P1())

    # Calculate coefficients for H = [50.0, 1000.0]
    intervals_data.append(fit_L3_P2())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L3_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_L3_dA()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1.0, 12.0]
    intervals_data.append(fit_L3_dA_P0())

    # Calculate coefficients for H = [12.0, 50.0]
    intervals_data.append(fit_L3_dA_P1())

    # Calculate coefficients for H = [50.0, 1000.0]
    intervals_data.append(fit_L3_dA_P2())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L3_dA_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_L3_dA_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 12.0
    fit_props.z_min = 1.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3_dA(x, y, z)[0].real, fit_props, "L3_dA_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_dA_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = False
    fit_props.z_max = 50.0
    fit_props.z_min = 12.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3_dA(x, y, z)[0].real, fit_props, "L3_dA_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_dA_P2()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 1000.0
    fit_props.z_min = 50.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3_dA(x, y, z)[0].real, fit_props, "L3_dA_P2")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_dB()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1.0, 12.0]
    intervals_data.append(fit_L3_dB_P0())

    # Calculate coefficients for H = [12.0, 50.0]
    intervals_data.append(fit_L3_dB_P1())

    # Calculate coefficients for H = [50.0, 1000.0]
    intervals_data.append(fit_L3_dB_P2())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "L3_dB_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_L3_dB_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 12.0
    fit_props.z_min = 1.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3_dB(x, y, z)[0].real, fit_props, "L3_dB_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_dB_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = False
    fit_props.z_max = 50.0
    fit_props.z_min = 12.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3_dB(x, y, z)[0].real, fit_props, "L3_dB_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_dB_P2()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 1000.0
    fit_props.z_min = 50.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3_dB(x, y, z)[0].real, fit_props, "L3_dB_P2")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_Afix_Bfix()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.num_x = 100
    fit_props.num_x_fit = 100
    fit_props.x_log_scale = True
    fit_props.x_max = 30.0
    fit_props.x_min = 1.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.cheby_order_x = 50
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    A = 1.0
    B = 1.0

    # Launch fit
    fit_integral_1d(lambda z: L3(A, B, z)[0].real, fit_props, f"L3_A: {A:0.2E} - B: {B:0.2E}", show_figs=True)


def fit_L3_Hfix()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    H = 50.0

    # Launch fit
    fit_integral_2d(lambda x, y: L3(x, y, H)[0].real, fit_props, f"L3_H: {H:0.2E}", show_figs=True)


def fit_L3_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 12.0
    fit_props.z_min = 1.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3(x, y, z)[0].real, fit_props, "L3_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = False
    fit_props.z_max = 50.0
    fit_props.z_min = 12.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3(x, y, z)[0].real, fit_props, "L3_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_L3_P2()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 1.0
    fit_props.y_min = 0.0
    fit_props.z_log_scale = True
    fit_props.z_max = 1000.0
    fit_props.z_min = 50.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: L3(x, y, z)[0].real, fit_props, "L3_P2")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-16, 0.1]
    intervals_data.append(fit_M1_P0())

    # Calculate coefficients for H = [0.1, 1.0]
    intervals_data.append(fit_M1_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M1_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_M1_dA()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-16, 0.1]
    intervals_data.append(fit_M1_dA_P0())

    # Calculate coefficients for H = [0.1, 1.0]
    intervals_data.append(fit_M1_dA_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M1_dA_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_M1_dA_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = True
    fit_props.z_max = 0.1
    fit_props.z_min = 1e-16
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M1_dA(x, y, z)[0].real, fit_props, "M1_dA_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1_dA_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 1.0
    fit_props.z_min = 0.1
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M1_dA(x, y, z)[0].real, fit_props, "M1_dA_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1_dB()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-16, 0.1]
    intervals_data.append(fit_M1_dB_P0())

    # Calculate coefficients for H = [0.1, 1.0]
    intervals_data.append(fit_M1_dB_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M1_dB_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_M1_dB_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = True
    fit_props.z_max = 0.1
    fit_props.z_min = 1e-16
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M1_dB(x, y, z)[0].real, fit_props, "M1_dB_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1_dB_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 1.0
    fit_props.z_min = 0.1
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M1_dB(x, y, z)[0].real, fit_props, "M1_dB_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = True
    fit_props.z_max = 0.1
    fit_props.z_min = 1e-16
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M1(x, y, z)[0].real, fit_props, "M1_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 1.0
    fit_props.z_min = 0.1
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M1(x, y, z)[0].real, fit_props, "M1_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M1_Hfix()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    H = 1.5

    # Launch fit
    fit_integral_2d(lambda x, y: M1(x, y, H)[0].real, fit_props, f"M1_H: {H:0.2E}", show_figs=True)


def fit_M2()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1e-12, 1.0]
    intervals_data.append(fit_M2_P0())

    # Calculate coefficients for H = [1.0, 50.0]
    intervals_data.append(fit_M2_P1())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M2_database.h5")
    write_intervals(database_path, intervals_data, 1)


def fit_M2_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max = 0.001
    fit_props.x_min = 1e-16
    fit_props.cheby_order_x = 50
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx = fit_integral_1d(lambda z: M2(z)[0].real, fit_props, f"M2", show_figs=True)

    return data_to_dict(fit_props, cheby_coeffs, ncx, array([]), array([]))


def fit_M2_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max = 1.0
    fit_props.x_min = 0.001
    fit_props.cheby_order_x = 50
    fit_props.cheby_tol = 1e-7
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx = fit_integral_1d(lambda z: M2(z)[0].real, fit_props, f"M2", show_figs=True)

    return data_to_dict(fit_props, cheby_coeffs, ncx, array([]), array([]))


def fit_M3()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1.0, 12.0]
    intervals_data.append(fit_M3_P0())

    # Calculate coefficients for H = [12.0, 50.0]
    intervals_data.append(fit_M3_P1())

    # Calculate coefficients for H = [50.0, 1000.0]
    intervals_data.append(fit_M3_P2())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M3_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_M3_dA()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1.0, 12.0]
    intervals_data.append(fit_M3_dA_P0())

    # Calculate coefficients for H = [12.0, 50.0]
    intervals_data.append(fit_M3_dA_P1())

    # Calculate coefficients for H = [50.0, 1000.0]
    intervals_data.append(fit_M3_dA_P2())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M3_dA_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_M3_dA_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 12.0
    fit_props.z_min = 1.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3_dA(x, y, z)[0].real, fit_props, "M3_dA_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_dA_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 50.0
    fit_props.z_min = 12.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3_dA(x, y, z)[0].real, fit_props, "M3_dA_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_dA_P2()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = True
    fit_props.z_max = 1000.0
    fit_props.z_min = 50.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3_dA(x, y, z)[0].real, fit_props, "M3_dA_P2")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_dB()->None:
    # Initialize object to storage the data
    intervals_data = []

    # Calculate coefficients for H = [1.0, 12.0]
    intervals_data.append(fit_M3_dB_P0())

    # Calculate coefficients for H = [12.0, 50.0]
    intervals_data.append(fit_M3_dB_P1())

    # Calculate coefficients for H = [50.0, 1000.0]
    intervals_data.append(fit_M3_dB_P2())

    # Save coefficients data into a database
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, "M3_dB_database.h5")
    write_intervals(database_path, intervals_data, 3)


def fit_M3_dB_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 12.0
    fit_props.z_min = 1.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3_dB(x, y, z)[0].real, fit_props, "M3_dB_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_dB_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 50.0
    fit_props.z_min = 12.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3_dB(x, y, z)[0].real, fit_props, "M3_dB_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_dB_P2()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = True
    fit_props.z_max = 1000.0
    fit_props.z_min = 50.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3_dB(x, y, z)[0].real, fit_props, "M3_dB_P2")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_P0()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 12.0
    fit_props.z_min = 1.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3(x, y, z)[0].real, fit_props, "M3_P0")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_P1()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = False
    fit_props.z_max = 50.0
    fit_props.z_min = 12.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3(x, y, z)[0].real, fit_props, "M3_P1")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def fit_M3_P2()->None:
    # Define parametric space and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 1.0
    fit_props.x_min = 0.0
    fit_props.y_max = 2.0
    fit_props.y_min = 1.0
    fit_props.z_log_scale = True
    fit_props.z_max = 1000.0
    fit_props.z_min = 50.0
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_order_z = 50
    fit_props.cheby_tol = 1e-8
    fit_props.fit_points_to_order()

    # Launch fit
    cheby_coeffs, ncx, ncy, ncz = fit_integral_3d(lambda x, y, z: M3(x, y, z)[0].real, fit_props, "M3_P2")

    return data_to_dict(fit_props, cheby_coeffs, ncx, ncy, ncz)


def write_intervals(database_path: str, intervals_data: list, dims: int) -> None:
    if dims == 1:
        dim_label = "x"
        fold_chn = None
    elif dims == 2:
        dim_label = "y"
        fold_chn = "ncx"
    elif dims == 3:
        dim_label = "z"
        fold_chn = "ncy"

    # Get intervals min position
    z_min = []
    for int_data in intervals_data:
        z_min.append(getattr(int_data["fit_props"], f"{dim_label}_min"))
    z_min = array(z_min)

    # Sort intervals by x_min values
    pos = argsort(z_min)
    intervals_sort = []
    for i in range(pos.shape[0]):
        intervals_sort.append(intervals_data[pos[i]])

    # Get intervals bounds
    intervals_bounds = concatenate((z_min[pos], array([getattr(intervals_sort[-1]["fit_props"], f"{dim_label}_max")])))

    # Get number of points per interval
    num_points = []
    for i in range(len(intervals_sort)):
        num_points.append(intervals_sort[i]["cheby_coeffs"].shape[0])
    num_points_cum = [0]
    for i,iv in enumerate(num_points):
        num_points_cum.append(num_points_cum[i]+iv)

    # Count folded points
    count_fold = 0
    count_fold_i = 0
    max_size_fold = 0
    if fold_chn is not None:
        for i in range(len(intervals_sort)):
            diff = np_abs(intervals_sort[i][fold_chn][1:]-intervals_sort[i][fold_chn][:-1])
            count_fold_i = (diff > 0.1).sum() + 1
            count_fold += count_fold_i
            if count_fold_i > max_size_fold:
                max_size_fold = count_fold_i

    # Open file unit
    fid = h5py.File(database_path, "w")

    # Storage data
    fid.create_dataset("dims", data=dims)
    fid.create_dataset("count_fold", data=count_fold)
    fid.create_dataset("intervals_bounds", data=intervals_bounds)
    fid.create_dataset("max_size_fold", data=max_size_fold)
    fid.create_dataset("num_points", data=num_points)
    fid.create_dataset("num_points_cum", data=num_points_cum)
    for i in range(len(intervals_sort)):
        # Create new group
        gp_int = fid.create_group(f"I{i:d}")

        # Fill intervals data
        gp_int.create_dataset("cheby_coeffs", data=intervals_sort[i]["cheby_coeffs"])
        gp_int.create_dataset("ncx", data=intervals_sort[i]["ncx"])
        gp_int.create_dataset("x_max", data=getattr(intervals_sort[i]["fit_props"], "x_max"))
        gp_int.create_dataset("x_min", data=getattr(intervals_sort[i]["fit_props"], "x_min"))
        gp_int.create_dataset("x_log_scale", data=getattr(intervals_sort[i]["fit_props"], "x_log_scale"))
        if dims >= 2:
            gp_int.create_dataset("ncy", data=intervals_sort[i]["ncy"])
            gp_int.create_dataset("y_max", data=getattr(intervals_sort[i]["fit_props"], "y_max"))
            gp_int.create_dataset("y_min", data=getattr(intervals_sort[i]["fit_props"], "y_min"))
            gp_int.create_dataset("y_log_scale", data=getattr(intervals_sort[i]["fit_props"], "y_log_scale"))
            if dims == 3:
                gp_int.create_dataset("ncz", data=intervals_sort[i]["ncz"])
                gp_int.create_dataset("z_max", data=getattr(intervals_sort[i]["fit_props"], "z_max"))
                gp_int.create_dataset("z_min", data=getattr(intervals_sort[i]["fit_props"], "z_min"))
                gp_int.create_dataset("z_log_scale", data=getattr(intervals_sort[i]["fit_props"], "z_log_scale"))


    # Close file unit
    fid.close()


if __name__ == "__main__":
    # fit_L1_Hfix()
    # fit_M1_Hfix()
    # fit_L1_Afix_Bfix()
    # fit_L1()
    # fit_L3()
    # fit_L3_Hfix()
    # fit_L3_Afix_Bfix()
    # fit_L1_P0()
    # fit_L1_P3()
    # fit_L2()
    # fit_L2_P1()
    # file_path = r"E:\sergio\developments\SeaMotions\aux_tools\test_coeffs.py"
    # fid = h5py.File(file_path, "w")
    # fid.create_dataset("C_filter", data=C_filter)
    # fid.create_dataset("NCX_filter", data=NCX_filter)
    # fid.create_dataset("NCY_filter", data=NCY_filter)
    # fid.create_dataset("NCZ_filter", data=NCZ_filter)
    # fid.close()
    # fit_M2()
    # fit_L3()
    # fit_L3_dA()
    # fit_L3_dB()
    # fit_L3_P2()
    # fit_M3()
    # fit_L1_dA()
    # fit_L1_dB()
    # fit_L3_dA()
    # fit_L3_dB()
    # fit_M1_dA()
    # fit_M1_dB()
    # fit_M3_dA()
    fit_M3_dB()