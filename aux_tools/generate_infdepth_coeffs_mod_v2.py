
# Import general usage libraries
import os

# Import general usage scientific libraries
import numpy as np
from numpy import arange, array, concatenate, ndarray, zeros

# Import local modules
from fit_wave_infdepth_v2 import (fit_residual_region_11, fit_residual_region_11_dx, fit_residual_region_11A_dx, fit_residual_region_11B_dx,
                                fit_residual_region_12, fit_residual_region_12_dx,
                                fit_residual_region_21, fit_residual_region_21_dx,
                                fit_residual_region_22, fit_residual_region_22_dx)
from fit_wave_infdepth_v2 import fit_residual_region_11_adaptive, fit_residual_region_11_dx_adaptive, fit_residual_region_AA_dx_adaptive
from fit_tools import RefLevel, write_coeffs_module_adaptive_2d_only_header


REF_COL = 60


def generate_region_11(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    # x, y, fit_results = fit_residual_region_11(show_summary_fig=show_summary_fig, show_figs=show_figs)
    fit_region = fit_residual_region_11_adaptive(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    # write_coeffs_module_adaptive( fit_region, folder_path, "R44" )
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "R11" )
    # write_coeffs_module(x, y, fit_results, folder_path, "R11")


def generate_region_11_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    # x, y, fit_results = fit_residual_region_11_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)
    fit_region = fit_residual_region_11_dx_adaptive(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    # write_coeffs_module(x, y, fit_results, folder_path, "R11_dX")
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "R11_dX" )


def generate_region_AA_dx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    # x, y, fit_results = fit_residual_region_11_dx(show_summary_fig=show_summary_fig, show_figs=show_figs)
    fit_region = fit_residual_region_AA_dx_adaptive(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    # write_coeffs_module(x, y, fit_results, folder_path, "R11_dX")
    write_coeffs_module_adaptive_only_header( fit_region, folder_path, "RAA_dX" )


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
    header_file_path = os.path.join(folder_path, module_name+"_newcode.hpp")
    fid = open(header_file_path, "w")

    # Add header guard
    fid.writelines(f"#ifndef __{module_name}_newcode_coeffs_hpp\n")
    fid.writelines(f"#define __{module_name}_newcode_coeffs_hpp\n\n")

    # Add includes
    fid.writelines('#include "../../config.hpp"\n')
    fid.writelines("\n")

    # Open namespace field
    fid.writelines("namespace newcode{\n")
    fid.writelines(f"namespace {module_name}" + "C{\n")

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
    fid.writelines("}\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close file unit for header file
    fid.close()

    # Open file unit for implementation file
    fid = open(os.path.join(folder_path, module_name+"_newcode.cpp"), "w")

    # Include module header file
    int_name_sp = f"{module_name}C"
    fid.writelines(f'\n#include "{module_name}_newcode.hpp"\n\n')

    fid.writelines("namespace newcode{\n")

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

    fid.writelines("}\n")


    # Close coefficients file unit
    fid.close()


def write_coeffs_module_adaptive( ref_level: RefLevel, folder_path: str, module_name: str)->None:
    # Get number of cumulative coefficients
    num_points_cum = ref_level.get_num_cheby_coeffs( )

    # Get cumulative data
    cum_coeffs      = ref_level.get_cheby_coeffs( )
    cheby_coeffs    = cum_coeffs[:, 0]
    ncx             = cum_coeffs[:, 1]
    ncy             = cum_coeffs[:, 2]

    # Get intervals size
    max_ref_level   = ref_level.get_max_level( )
    x_min           = ref_level.fit_props.x_min
    x_max           = ref_level.fit_props.x_max
    y_min           = ref_level.fit_props.y_min
    y_max           = ref_level.fit_props.y_max

    # Get log scales
    x_log_scale     = ref_level.fit_props.x_log_scale
    y_log_scale     = ref_level.fit_props.y_log_scale
    fcn_log_scale   = ref_level.fit_props.fcn_log_scale

    # Calculate hash table
    dx  = ( x_max - x_min ) / 2**(max_ref_level)
    x   = arange( x_min, x_max+dx, dx )
    xm  = ( x[1:] + x[:-1] ) / 2.0

    dy  = ( y_max - y_min ) / 2**(max_ref_level)
    y   = arange( y_min, y_max+dy, dy )
    ym  = ( y[1:] + y[:-1] ) / 2.0

    blocks_np           = xm.shape[0] * ym.shape[0]
    blocks_start        = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np    = np.zeros( ( blocks_np, ), dtype=int )
    count               = 0
    for i in range( xm.shape[0] ):
        xmi = xm[i]
        for j in range( ym.shape[0] ):
            ymi = ym[j]
            bs, bc                  = ref_level.get_start_index( np.array( [ xmi, ymi ] ) )
            blocks_start[count]     = bs
            blocks_coeffs_np[count] = bc
            count                   += 1

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
    fid.writelines("namespace newcode{\n")
    fid.writelines(f"namespace {module_name}" + "C{\n")

    fid.writelines(str_start_col("    extern const int", "max_ref_level;\n", REF_COL))

    # Save number of points
    fid.writelines(str_start_col("    extern const int", "blocks_start[];\n", REF_COL))
    fid.writelines(str_start_col("    extern const int", "blocks_coeffs_np[];\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col("    extern const bool", "fcn_log_scale;\n", REF_COL))

    fid.writelines(str_start_col("    extern const bool", "x_log_scale;\n", REF_COL))
    fid.writelines(str_start_col("    extern const cusfloat", "x_max;\n", REF_COL))
    fid.writelines(str_start_col("    extern const cusfloat", "x_min;\n", REF_COL))

    fid.writelines(str_start_col("    extern const bool", "y_log_scale;\n", REF_COL))
    fid.writelines(str_start_col("    extern const cusfloat", "y_max;\n", REF_COL))
    fid.writelines(str_start_col("    extern const cusfloat", "y_min;\n", REF_COL))
    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("    extern const int", "num_c;\n", REF_COL))
    fid.writelines(str_start_col("    extern const cusfloat", "c[];\n", REF_COL))

    # Write polynomials coefficients
    fid.writelines(str_start_col("    extern const int", "ncx[];\n", REF_COL))
    fid.writelines(str_start_col("    extern const int", "ncy[];\n", REF_COL))

    # Close namespace field
    fid.writelines("}\n")
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

    fid.writelines("namespace newcode{\n")

    # Save number of intervals
    fid.writelines(str_start_col(f"const int", f"{int_name_sp}::max_ref_level = {max_ref_level:d};\n", REF_COL))

    # Save number of points
    blocks_start_str = ", ".join(f"{i:d}" for i in blocks_start)
    fid.writelines(str_start_col(f"const int", f"{int_name_sp}::blocks_start[{blocks_np}] = " + "{" + blocks_start_str + "};\n", REF_COL))
    blocks_coeffs_np_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np)
    fid.writelines(str_start_col(f"const int", f"{int_name_sp}::blocks_coeffs_np[{blocks_np}] = " + "{" + blocks_coeffs_np_str + "};\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col(f"const bool", f"{int_name_sp}::fcn_log_scale = {fcn_log_scale:d};\n", REF_COL))

    fid.writelines(str_start_col(f"const bool", f"{int_name_sp}::x_log_scale = {x_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"const cusfloat", f"{int_name_sp}::x_max = {x_max:0.6f};\n", REF_COL))
    fid.writelines(str_start_col(f"const cusfloat", f"{int_name_sp}::x_min = {x_min:0.6f};\n", REF_COL))

    fid.writelines(str_start_col(f"const bool", f"{int_name_sp}::y_log_scale = {y_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"const cusfloat", f"{int_name_sp}::y_max = {y_max:0.6f};\n", REF_COL))
    fid.writelines(str_start_col(f"const cusfloat", f"{int_name_sp}::y_min = {y_min:0.6f};\n", REF_COL))

    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("const int", f"{int_name_sp}::num_c = {cheby_coeffs.shape[0]};\n", REF_COL))
    write_vector(fid, cheby_coeffs, "c", "cusfloat", int_name_sp)

    # Write polynomials coefficients
    write_vector(fid, ncx, "ncx", "int", int_name_sp)
    write_vector(fid, ncy, "ncy", "int", int_name_sp)

    fid.writelines("}\n")


    # Close coefficients file unit
    fid.close()


def write_coeffs_module_adaptive_only_header( ref_level: RefLevel, folder_path: str, module_name: str)->None:
    # Get number of cumulative coefficients
    num_points_cum = ref_level.get_num_cheby_coeffs( )

    # Get cumulative data
    cum_coeffs      = ref_level.get_cheby_coeffs( )
    cheby_coeffs    = cum_coeffs[:, 0]
    ncx             = cum_coeffs[:, 1]
    ncy             = cum_coeffs[:, 2]
    max_cheby_order = cum_coeffs[:, 1:].max( )

    # Get intervals size
    max_ref_level   = ref_level.get_max_level( )
    x_min           = ref_level.fit_props.x_min
    x_max           = ref_level.fit_props.x_max
    y_min           = ref_level.fit_props.y_min
    y_max           = ref_level.fit_props.y_max

    # Get log scales
    x_log_scale     = ref_level.fit_props.x_log_scale
    y_log_scale     = ref_level.fit_props.y_log_scale
    fcn_log_scale   = ref_level.fit_props.fcn_log_scale

    # Calculate hash table
    intervals_np    = 2**(max_ref_level)
    dx              = ( x_max - x_min ) / intervals_np
    x               = arange( x_min, x_max+dx, dx )
    xm              = ( x[1:] + x[:-1] ) / 2.0

    dy              = ( y_max - y_min ) / intervals_np
    y               = arange( y_min, y_max+dy, dy )
    ym              = ( y[1:] + y[:-1] ) / 2.0

    blocks_np           = xm.shape[0] * ym.shape[0]
    blocks_start        = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np    = np.zeros( ( blocks_np, ), dtype=int )
    x_min_vec           = np.zeros( ( blocks_np, ), dtype=float )
    x_max_vec           = np.zeros( ( blocks_np, ), dtype=float )
    y_min_vec           = np.zeros( ( blocks_np, ), dtype=float )
    y_max_vec           = np.zeros( ( blocks_np, ), dtype=float )
    count               = 0
    for i in range( xm.shape[0] ):
        xmi = xm[i]
        for j in range( ym.shape[0] ):
            ymi = ym[j]
            bs, bc, _, _, x_min_i, x_max_i, y_min_i, y_max_i, _, _  = ref_level.get_start_index( np.array( [ xmi, ymi ] ) )
            blocks_start[count]                                     = bs
            blocks_coeffs_np[count]                                 = bc
            x_min_vec[count]                                        = x_min_i
            x_max_vec[count]                                        = x_max_i
            y_min_vec[count]                                        = y_min_i
            y_max_vec[count]                                        = y_max_i
            count                                                   += 1

    dx_vec = x_max_vec - x_min_vec
    dy_vec = y_max_vec - y_min_vec
        
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
    fid.writelines(f"struct {module_name}" + "C\n{\n")

    # Save number of intervals
    fid.writelines(str_start_col(f"static constexpr int", f"max_ref_level = {max_ref_level:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"intervals_np = {intervals_np:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"max_cheby_order = {int(max_cheby_order):d};\n", REF_COL))

    # Save number of points
    blocks_start_str = ", ".join(f"{i:d}" for i in blocks_start)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_start[{blocks_np}] = " + "{" + blocks_start_str + "};\n", REF_COL))
    blocks_coeffs_np_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_coeffs_np[{blocks_np}] = " + "{" + blocks_coeffs_np_str + "};\n", REF_COL))
    x_min_str = ", ".join(f"{v:0.3E}" for v in x_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_region[{blocks_np}] = " + "{" + x_min_str + "};\n", REF_COL))
    x_max_str = ", ".join(f"{v:0.3E}" for v in x_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_region[{blocks_np}] = " + "{" + x_max_str + "};\n", REF_COL))
    dx_str = ", ".join(f"{v:0.3E}" for v in dx_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_region[{blocks_np}] = " + "{" + dx_str + "};\n", REF_COL))
    y_min_str = ", ".join(f"{v:0.3E}" for v in y_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_min_region[{blocks_np}] = " + "{" + y_min_str + "};\n", REF_COL))
    y_max_str = ", ".join(f"{v:0.3E}" for v in y_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_max_region[{blocks_np}] = " + "{" + y_max_str + "};\n", REF_COL))
    dy_str = ", ".join(f"{v:0.3E}" for v in dy_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dy_region[{blocks_np}] = " + "{" + dy_str + "};\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col(f"static constexpr bool", f"fcn_log_scale = {fcn_log_scale:d};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"x_log_scale = {x_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_global = {x_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_global = {x_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_min_region = {dx:0.16f};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"y_log_scale = {y_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_max_global = {y_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_min_global = {y_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dy_min_region = {dy:0.16f};\n", REF_COL))

    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("static constexpr std::size_t", f"num_cf = {cheby_coeffs.shape[0]};\n", REF_COL))
    write_vector_header(fid, cheby_coeffs, "cf", "cusfloat")

    # Write polynomials coefficients
    write_vector_header(fid, ncx, "ncxf", "std::size_t")
    write_vector_header(fid, ncy, "ncyf", "std::size_t")


    # Close namespace field
    fid.writelines("};\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close file unit for header file
    fid.close()


def write_vector(fid, field: ndarray, field_tag: str, var_type: str, int_name: str)->None:
    fid.writelines(str_start_col(f"alignas(FLOATING_PRECISION) const {var_type}", f"{int_name}::{field_tag}[{field.shape[0]}] = " + "{\n", REF_COL))
    if var_type == "int":
        for i, iv in enumerate(field):
            fid.writelines(f"                                          {int(iv):d},  // {field_tag}[{i}]\n")
    elif var_type == "cusfloat":
        for i, iv in enumerate(field):
            fid.writelines(f"                                          {iv:0.16f},  // {field_tag}[{i}]\n")
    fid.writelines(f"                            " + "};\n")


def write_vector_header( fid, field: ndarray, field_tag: str, var_type: str, keyword="constexpr" )->None:
    fid.writelines(str_start_col(f"alignas(FLOATING_PRECISION) static {keyword:s} {var_type}", f"{field_tag}[{field.shape[0]}] = " + "{\n", REF_COL))
    if var_type == "std::size_t":
        for i, iv in enumerate(field):
            fid.writelines(f"                                                                {int(iv):d},  // {field_tag}[{i}]\n")
    elif var_type == "cusfloat":
        for i, iv in enumerate(field):
            fid.writelines(f"                                                                {iv:0.16f},  // {field_tag}[{i}]\n")
    fid.writelines(f"                                                  " + "};\n")


def write_vector_line(fid, data: ndarray, field_name: str, var_type: str, int_name: str):
    if var_type == "cusfloat":
        field_data_str = ", ".join(f"{i:0.6E}" for i in data)
    elif var_type == "bool":
        field_data_str = ", ".join(f"{int(i):d}" for i in data)
    else:
        raise ValueError("Variable type not recognized!.")
    fid.writelines(str_start_col(f"{var_type}", f"{int_name}::{field_name}[{data.shape[0]}] = " + "{" + field_data_str + "};\n", REF_COL))


if __name__ == "__main__":
    import sys

    folder_path = sys.argv[1]
    option      = sys.argv[2]
    if option == "region_11":
        generate_region_11(folder_path, show_figs=False, show_summary_fig=True)
    
    elif option == "region_11_dx":
        generate_region_11_dx(folder_path, show_figs=False, show_summary_fig=True)
    
    else:
        raise ValueError( f"Option value: {option} - not valid" )
    
    # this_path = os.path.dirname(os.path.abspath(__file__))
    # folder_path = os.path.join(this_path, "0_databases", "0_infinite_water_depth", "Prec_1Em6")
    # generate_region_11(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_11_dx(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_AA_dx(folder_path, show_figs=False, show_summary_fig=True)
    # generate_region_11A_dx(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_11B_dx(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_12(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_12_dx(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_21(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_21_dx(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_22(folder_path, show_figs=False, show_summary_fig=False)
    # generate_region_22_dx(folder_path, show_figs=False, show_summary_fig=False)