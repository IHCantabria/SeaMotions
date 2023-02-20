
# Import general usage libraries
import os

# Import general usage scientific libraries
from numpy import array, concatenate, ndarray, zeros

# Import local modules
from fit_wave_infdepth import (fit_residual_region_11, fit_residual_region_11_dx, fit_residual_region_11A_dx, fit_residual_region_11B_dx,
                                fit_residual_region_12, fit_residual_region_12_dx,
                                fit_residual_region_21, fit_residual_region_21_dx,
                                fit_residual_region_22, fit_residual_region_22_dx)


REF_COL = 24


def generate_region_11(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_11(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R11")


def generate_region_11_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_11_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R11_dX")


def generate_region_11A_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_11A_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R11A_dX")


def generate_region_11B_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_11B_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R11B_dX")


def generate_region_12(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_12(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R12")


def generate_region_12_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_12_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R12_dX")


def generate_region_21(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_21(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R21")


def generate_region_21_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_21_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R21_dX")


def generate_region_22(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_22(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R22")


def generate_region_22_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    x, y, fit_results = fit_residual_region_22_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module(x, y, fit_results, folder_path, "R22_dX")


def str_start_col(str0: str, str1: str, num_col: int)->str:
    l0 = len(str0)
    if l0 >= num_col:
        raise ValueError("First string is bigger than the reference column.")

    for i in range(num_col-l0):
        str0 += " "
    
    return str0+str1


def write_coeffs_module(x: ndarray, y: ndarray, results: list, folder_path: str, module_name: str)->None:
    # Get number of cumulative coefficients
    num_points_cum = [0]
    num_results = len(results)
    for i in range(num_results):
        num_points_cum.append(num_points_cum[i]+results[i][0].shape[0])
    num_points_cum = array(num_points_cum)
    num_points = num_points_cum[1:]

    # Get cumulative data
    cheby_coeffs = zeros((num_points_cum[-1], ))
    ncx = zeros((num_points_cum[-1], ))
    ncy = zeros((num_points_cum[-1], ))
    for i in range(num_results):
        idx0 = num_points_cum[i]
        idx1 = num_points_cum[i+1]
        cheby_coeffs[idx0:idx1] = results[i][0]
        ncx[idx0:idx1] = results[i][1]
        ncy[idx0:idx1] = results[i][2]

    # Get intervals size
    num_intervals_x = x.shape[0]-1
    num_intervals_y = y.shape[0]-1

    # Get X log scale intervals data
    x_log_scale = zeros((num_intervals_x, ))
    for i in range(num_intervals_x):
        x_log_scale[i] = results[i*num_intervals_y][4].x_log_scale

    # Get Y log scale intervals data
    y_log_scale = zeros((num_intervals_y, ))
    for i in range(num_intervals_y):
        y_log_scale[i] = results[i][4].y_log_scale


    # Open file unit to storage the coefficients in the 
    # database
    header_file_path = os.path.join(folder_path, module_name+".hpp")
    fid = open(header_file_path, "w")

    # Add header guard
    fid.writelines(f"#ifndef __{module_name}_coeffs_hpp\n")
    fid.writelines(f"#define __{module_name}_coeffs_hpp\n\n")

    # Add includes
    fid.writelines('#include "../../config.hpp"\n')
    fid.writelines("\n")

    # Open namespace field
    fid.writelines(f"namespace {module_name}" + "C\n{\n")

    fid.writelines(str_start_col("    extern int", "dims;\n", REF_COL))
    fid.writelines(str_start_col("    extern int", "num_intervals_x;\n", REF_COL))
    fid.writelines(str_start_col("    extern int", "num_intervals_y;\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "interval_bounds_x[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "interval_bounds_y[];\n", REF_COL))

    # Save number of points
    fid.writelines(str_start_col("    extern int", "num_points[];\n", REF_COL))
    fid.writelines(str_start_col("    extern int", "num_points_cum[];\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col("    extern bool", "x_log_scale[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "x_map_scale[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "x_map_scale_log[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "x_max[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "x_min[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "x_min_l10[];\n", REF_COL))

    fid.writelines(str_start_col("    extern bool", "y_log_scale[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "y_map_scale[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "y_map_scale_log[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "y_max[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "y_min[];\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "y_min_l10[];\n", REF_COL))
    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("    extern int", "num_c;\n", REF_COL))
    fid.writelines(str_start_col("    extern cusfloat", "c[];\n", REF_COL))

    # Write polynomials coefficients
    fid.writelines(str_start_col("    extern int", "ncx[];\n", REF_COL))
    fid.writelines(str_start_col("    extern int", "ncy[];\n", REF_COL))

    # Close namespace field
    fid.writelines("}\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close file unit for header file
    fid.close()

    # Open file unit for implementation file
    fid = open(os.path.join(folder_path, module_name+".cpp"), "w")

    # Include module header file
    int_name_sp = f"{module_name}C"
    fid.writelines(f'\n#include "{module_name}.hpp"\n\n')

    # Save number of intervals
    interval_str_x = ", ".join(f"{i:0.6E}" for i in x)
    interval_str_y = ", ".join(f"{i:0.6E}" for i in y)
    fid.writelines(str_start_col(f"int", f"{int_name_sp}::num_intervals_x = {num_intervals_x:d};\n", REF_COL))
    fid.writelines(str_start_col(f"int", f"{int_name_sp}::num_intervals_y = {num_intervals_y:d};\n", REF_COL))
    fid.writelines(str_start_col(f"cusfloat", f"{int_name_sp}::interval_bounds_x[{x.shape[0]:d}] = " + "{" + interval_str_x + "};\n\n", REF_COL))
    fid.writelines(str_start_col(f"cusfloat", f"{int_name_sp}::interval_bounds_y[{y.shape[0]:d}] = " + "{" + interval_str_y + "};\n\n", REF_COL))

    # Save number of points
    num_points_str = ", ".join(f"{i:d}" for i in num_points)
    num_points_cum_str = ", ".join(f"{i:d}" for i in num_points_cum)
    fid.writelines(str_start_col("int", f"{int_name_sp}::num_points[{num_points.shape[0]}] = " + "{" + num_points_str + "};\n", REF_COL))
    fid.writelines(str_start_col("int", f"{int_name_sp}::num_points_cum[{num_points_cum.shape[0]}] = " + "{" + num_points_cum_str + "};\n", REF_COL))

    # Save interval bounds
    zeros_map = zeros((num_intervals_x, ))
    write_vector_line(fid, x_log_scale, "x_log_scale", "bool", int_name_sp)
    write_vector_line(fid, zeros_map, "x_map_scale", "cusfloat", int_name_sp)
    write_vector_line(fid, zeros_map, "x_map_scale_log", "cusfloat", int_name_sp)
    write_vector_line(fid, x[1:], "x_max", "cusfloat", int_name_sp)
    write_vector_line(fid, x[:-1], "x_min", "cusfloat", int_name_sp)
    write_vector_line(fid, zeros_map, "x_min_l10", "cusfloat", int_name_sp)

    zeros_map = zeros((num_intervals_y, ))
    write_vector_line(fid, y_log_scale, "y_log_scale", "bool", int_name_sp)
    write_vector_line(fid, zeros_map, "y_map_scale", "cusfloat", int_name_sp)
    write_vector_line(fid, zeros_map, "y_map_scale_log", "cusfloat", int_name_sp)
    write_vector_line(fid, y[1:], "y_max", "cusfloat", int_name_sp)
    write_vector_line(fid, y[:-1], "y_min", "cusfloat", int_name_sp)
    write_vector_line(fid, zeros_map, "y_min_l10", "cusfloat", int_name_sp)
    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("int", f"{int_name_sp}::num_c = {cheby_coeffs.shape[0]};\n", REF_COL))
    write_vector(fid, cheby_coeffs, "c", "cusfloat", int_name_sp)

    # Write polynomials coefficients
    write_vector(fid, ncx, "ncx", "int", int_name_sp)
    write_vector(fid, ncy, "ncy", "int", int_name_sp)

    # Close coefficients file unit
    fid.close()


def write_vector(fid, field: ndarray, field_tag: str, var_type: str, int_name: str)->None:
    fid.writelines(str_start_col(f"{var_type}", f"{int_name}::{field_tag}[{field.shape[0]}] = " + "{\n", REF_COL))
    if var_type == "int":
        for i, iv in enumerate(field):
            fid.writelines(f"                                          {int(iv):d},  // {field_tag}[{i}]\n")
    elif var_type == "cusfloat":
        for i, iv in enumerate(field):
            fid.writelines(f"                                          {iv:0.16f},  // {field_tag}[{i}]\n")
    fid.writelines(f"                            " + "};\n")


def write_vector_line(fid, data: ndarray, field_name: str, var_type: str, int_name: str):
    if var_type == "cusfloat":
        field_data_str = ", ".join(f"{i:0.6E}" for i in data)
    elif var_type == "bool":
        field_data_str = ", ".join(f"{int(i):d}" for i in data)
    else:
        raise ValueError("Variable type not recognized!.")
    fid.writelines(str_start_col(f"{var_type}", f"{int_name}::{field_name}[{data.shape[0]}] = " + "{" + field_data_str + "};\n", REF_COL))


if __name__ == "__main__":
    this_path = os.path.dirname(os.path.abspath(__file__))
    folder_path = os.path.join(os.path.dirname(this_path), "src", "green", "inf_depth_coeffs")
    # generate_region_11(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_11_dx(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_11A_dx(folder_path, show_figs=False, show_summary_fig=True)
    generate_region_11B_dx(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_12(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_12_dx(folder_path, show_figs=True, show_summary_fig=True)
    # generate_region_21(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_21_dx(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_22(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_22_dx(folder_path, show_figs=False, show_summary_fig=True)