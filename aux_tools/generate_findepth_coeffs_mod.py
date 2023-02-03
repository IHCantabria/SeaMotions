
# Import general usage libraries
import h5py
import os

# Import general usage scientific libraries
from numpy import array, concatenate, ndarray, zeros


REF_COL = 25


def generate_coeffs_modules(database_path: str, file_path: str, int_name: str)->None:
    # Open database file unit
    fid_db = h5py.File(database_path, "r")
    dims = fid_db["dims"][()]

    # Open file unit to storage the coefficients in the 
    # database
    fid = open(file_path, "w")

    # Add header guard
    fid.writelines(f"#ifndef __{int_name}_coeffs_hpp\n")
    fid.writelines(f"#define __{int_name}_coeffs_hpp\n\n")

    # Add includes
    fid.writelines('#include "../../config.hpp"\n')
    fid.writelines("\n")

    # Open namespace field
    fid.writelines(f"namespace {int_name}" + "C\n{\n")

    # Save number of intervals
    intervals_bounds = fid_db["intervals_bounds"][:]
    num_intervals = intervals_bounds.shape[0]-1
    interval_str = ", ".join(f"{i:0.6E}" for i in intervals_bounds)
    fid.writelines(str_start_col("    int", f"num_intervals = {num_intervals:d};\n", REF_COL))
    fid.writelines(str_start_col("    cusfloat", f"interval_bounds[{intervals_bounds.shape[0]:d}] = " + "{" + interval_str + "};\n\n", REF_COL))

    # Save number of points
    max_size_fold = fid_db["max_size_fold"][()]
    num_points = fid_db["num_points"][:]
    num_points_cum = fid_db["num_points_cum"][:]
    num_points_str = ", ".join(f"{i:d}" for i in num_points)
    num_points_cum_str = ", ".join(f"{i:d}" for i in num_points_cum)
    fid.writelines(str_start_col("    int", f"num_points[{num_points.shape[0]}] = " + "{" + num_points_str + "};\n", REF_COL))
    fid.writelines(str_start_col("    int", f"num_points_cum[{num_points_cum.shape[0]}] = " + "{" + num_points_cum_str + "};\n", REF_COL))
    fid.writelines(str_start_col("    int", f"max_size_fold = {max_size_fold:d};\n\n", REF_COL))

    # Save interval bounds
    zeros_map = zeros((num_intervals, ))
    write_interval_bounds(fid, fid_db, "x_log_scale", num_intervals, "bool")
    write_vector_line(fid, zeros_map, "x_map_scale", "cusfloat")
    write_vector_line(fid, zeros_map, "x_map_scale_log", "cusfloat")
    write_interval_bounds(fid, fid_db, "x_max", num_intervals, "cusfloat")
    write_interval_bounds(fid, fid_db, "x_min", num_intervals, "cusfloat")
    write_vector_line(fid, zeros_map, "x_min_l10", "cusfloat")
    if dims >=2:
        write_interval_bounds(fid, fid_db, "y_log_scale", num_intervals, "bool")
        write_vector_line(fid, zeros_map, "y_map_scale", "cusfloat")
        write_vector_line(fid, zeros_map, "y_map_scale_log", "cusfloat")
        write_interval_bounds(fid, fid_db, "y_max", num_intervals, "cusfloat")
        write_interval_bounds(fid, fid_db, "y_min", num_intervals, "cusfloat")
        write_vector_line(fid, zeros_map, "y_min_l10", "cusfloat")
        if dims == 3:
            write_interval_bounds(fid, fid_db, "z_log_scale", num_intervals, "bool")
            write_vector_line(fid, zeros_map, "z_map_scale", "cusfloat")
            write_vector_line(fid, zeros_map, "z_map_scale_log", "cusfloat")
            write_interval_bounds(fid, fid_db, "z_max", num_intervals, "cusfloat")
            write_interval_bounds(fid, fid_db, "z_min", num_intervals, "cusfloat")
            write_vector_line(fid, zeros_map, "z_min_l10", "cusfloat")
    fid.writelines("\n")

    # Write chebyshev polynomials
    cheby_coeffs = interval_to_vector(fid_db, "cheby_coeffs", intervals_bounds.shape[0]-1)
    fid.writelines(str_start_col("    int", f"num_c = {cheby_coeffs.shape[0]};\n", REF_COL))
    fid.writelines(str_start_col("    cusfloat", f"c[{cheby_coeffs.shape[0]}] = " + "{\n", REF_COL))
    for i, iv in enumerate(cheby_coeffs):
        fid.writelines(f"                                   {iv:0.16E},  // C[{i}]\n")
    fid.writelines(f"                                " + "};\n")
    cf = zeros((max_size_fold, ))
    write_vector(fid, cf, "cf", "cusfloat")

    # Write polynomials coefficients
    ncx = interval_to_vector(fid_db, "ncx", num_intervals)
    ncxf = zeros((max_size_fold, ))
    write_vector(fid, ncx, "ncx", "int")
    write_vector(fid, ncxf, "ncxf", "int")
    if dims >= 2:
        ncy = interval_to_vector(fid_db, "ncy", num_intervals)
        ncyf = zeros((max_size_fold, ))
        write_vector(fid, ncy, "ncy", "int")
        write_vector(fid, ncyf, "ncyf", "int")
        if dims == 3:
            ncz = interval_to_vector(fid_db, "ncz", num_intervals)
            write_vector(fid, ncz, "ncz", "int")
        else:
            raise ValueError(f"Number of dimensions: {dims} not available.")

    # Close namespace field
    fid.writelines("}\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close coefficients file unit
    fid.close()

    # Close database file unit
    fid_db.close()


def generate_test_database()->None:
    # Define file path for the database
    this_path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(this_path, "test_coeffs.h5")

    # Loads coefficients
    fid = h5py.File(file_path, "r")
    C_filter = fid["C_filter"][:]
    NCX_filter = fid["NCX_filter"][:]
    NCY_filter = fid["NCY_filter"][:]
    NCZ_filter = fid["NCZ_filter"][:]
    fid.close()

    # Create database
    intervals = array([1e-12, 1e-1, 3.0, 10.0, 50.0])
    database_file_path = os.path.join(this_path, "test_database.h5")
    fid = h5py.File(database_file_path, "w")
    fid.create_dataset("intervals", data=intervals)
    for i in range(intervals.shape[0]):
        gp_int = fid.create_group(f"I{i:d}")
        gp_int.create_dataset("C_filter", data=C_filter)
        gp_int.create_dataset("NCX_filter", data=NCX_filter)
        gp_int.create_dataset("NCY_filter", data=NCY_filter)
        gp_int.create_dataset("NCZ_filter", data=NCZ_filter)
    
    fid.close()


def interval_to_vector(fid_db, field_name: str, num_intervals: int, dims=1)->ndarray:
    data = array([])
    for i in range(num_intervals):
        if dims == 0:
            data = concatenate((data, array([fid_db[f"I{i:d}"][field_name][()]])))
        elif dims == 1:
            data = concatenate((data, fid_db[f"I{i:d}"][field_name]))
        else:
            raise ValueError("Dimensions error.")

    return data


def str_start_col(str0: str, str1: str, num_col: int)->str:
    l0 = len(str0)
    if l0 >= num_col:
        raise ValueError("First string is bigger than the reference column.")

    for i in range(num_col-l0):
        str0 += " "
    
    return str0+str1


def write_interval_bounds(fid, fid_db, field_name: str, num_intervals: int, var_type: str)->None:
    field_data = interval_to_vector(fid_db, field_name, num_intervals, dims=0)
    if var_type == "cusfloat":
        field_data_str = ", ".join(f"{i:0.6E}" for i in field_data)
    elif var_type == "bool":
        field_data_str = ", ".join(f"{int(i):d}" for i in field_data)
    else:
        raise ValueError("Variable type not recognized!.")
    fid.writelines(str_start_col(f"    {var_type}", f"{field_name}[{field_data.shape[0]}] = " + "{" + field_data_str + "};\n", REF_COL))


def write_vector(fid, field: ndarray, field_tag: str, var_type: str)->None:
    fid.writelines(str_start_col(f"    {var_type}", f"{field_tag}[{field.shape[0]}] = " + "{\n", REF_COL))
    for i, iv in enumerate(field):
        fid.writelines(f"                                   {int(iv):d},  // {field_tag}[{i}]\n")
    fid.writelines(f"                            " + "};\n")


def write_vector_line(fid, data: ndarray, field_name: str, var_type: str):
    if var_type == "cusfloat":
        field_data_str = ", ".join(f"{i:0.6E}" for i in data)
    elif var_type == "bool":
        field_data_str = ", ".join(f"{int(i):d}" for i in data)
    else:
        raise ValueError("Variable type not recognized!.")
    fid.writelines(str_start_col(f"    {var_type}", f"{field_name}[{data.shape[0]}] = " + "{" + field_data_str + "};\n", REF_COL))


if __name__ == "__main__":
    # generate_test_database()
    int_name = "L2"
    this_path = os.path.dirname(os.path.abspath(__file__))
    database_path = os.path.join(this_path, f"{int_name}_database.h5")
    file_path = os.path.join(os.path.dirname(this_path), "src", "green", "fin_depth_coeffs", f"{int_name}.hpp")
    generate_coeffs_modules(database_path, file_path, int_name)