#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <array>
#include <vector>

#include <hdf5.h>

#include "../../src/green/triple_hankel_far_field.hpp"

template <typename T>
T read_attribute(hid_t file, const char* name, hid_t h5_type)
{
    hid_t attr = H5Aopen( file, name, H5P_DEFAULT );
    if ( attr < 0 ) 
    {
        throw std::runtime_error(std::string("Could not open attribute: ") + name);
    }

    T value{ };
    if ( H5Aread( attr, h5_type, &value ) < 0)
    {
        H5Aclose( attr );
        throw std::runtime_error(std::string("Could not read attribute: ") + name);
    }

    H5Aclose( attr );
    return value;
}

template <typename T>
T read_attribute_or_default(hid_t file, const char* name, hid_t h5_type, T default_value)
{
    htri_t exists = H5Aexists(file, name);
    if (exists < 0) {
        throw std::runtime_error(std::string("Could not check attribute existence: ") + name);
    }
    if (exists == 0) {
        return default_value;
    }
    return read_attribute<T>(file, name, h5_type);
}

template <typename T>
std::vector<T> read_dataset_1d(hid_t file, const char* dataset_name, hid_t h5_type)
{
    hid_t dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);
    if (dataset < 0) {
        throw std::runtime_error(std::string("Could not open dataset: ") + dataset_name);
    }

    hid_t space = H5Dget_space(dataset);
    if (space < 0) {
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Could not get dataspace: ") + dataset_name);
    }

    int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 1) {
        H5Sclose(space);
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Unexpected rank for dataset: ") + dataset_name);
    }

    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);
    std::vector<T> values(static_cast<size_t>(dim));

    if (H5Dread(dataset, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
        H5Sclose(space);
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Could not read dataset: ") + dataset_name);
    }

    H5Sclose(space);
    H5Dclose(dataset);
    return values;
}

template <typename T>
std::vector<T> read_dataset_6d(
    hid_t file,
    const char* dataset_name,
    hid_t h5_type,
    std::array<hsize_t, 6>& dims)
{
    hid_t dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);
    if (dataset < 0) {
        throw std::runtime_error(std::string("Could not open dataset: ") + dataset_name);
    }

    hid_t space = H5Dget_space(dataset);
    if (space < 0) {
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Could not get dataspace: ") + dataset_name);
    }

    int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 6) {
        H5Sclose(space);
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Unexpected rank for dataset: ") + dataset_name);
    }

    std::array<hsize_t, 6> tmp{0, 0, 0, 0, 0, 0};
    H5Sget_simple_extent_dims(space, tmp.data(), nullptr);
    dims = tmp;

    size_t total_size = 1;
    for (hsize_t d : dims) {
        total_size *= static_cast<size_t>(d);
    }

    std::vector<T> values(total_size);
    if (H5Dread(dataset, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
        H5Sclose(space);
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Could not read dataset: ") + dataset_name);
    }

    H5Sclose(space);
    H5Dclose(dataset);
    return values;
}

size_t flatten_index(const std::array<hsize_t, 6>& dims, const std::array<hsize_t, 6>& idx)
{
    return (((((static_cast<size_t>(idx[0]) * static_cast<size_t>(dims[1]) + static_cast<size_t>(idx[1])) * static_cast<size_t>(dims[2]) + static_cast<size_t>(idx[2])) * static_cast<size_t>(dims[3]) + static_cast<size_t>(idx[3])) * static_cast<size_t>(dims[4]) + static_cast<size_t>(idx[4])) * static_cast<size_t>(dims[5]) + static_cast<size_t>(idx[5]));
}

std::vector<hsize_t> build_validation_indices(hsize_t dim, bool full_mode)
{
    std::vector<hsize_t> out;
    if (dim == 0) {
        return out;
    }

    if (full_mode) {
        out.resize(static_cast<size_t>(dim));
        for (hsize_t i = 0; i < dim; ++i) {
            out[static_cast<size_t>(i)] = i;
        }
        return out;
    }

    // Quick mode: sample representative locations per dimension.
    const hsize_t mid = dim / 2;
    out.push_back(0);
    if (mid != 0 && mid < dim) {
        out.push_back(mid);
    }
    if (dim > 1) {
        out.push_back(dim - 1);
    }

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}


int main(int argc, char* argv[]) 
{
    if ( argc < 2 )
    {
        std::cerr << "Usage: test_triple_hankel_mod_database_cmd <path_to_hdf5_db> [quick|full]\n";
        return 1;
    }

    bool full_mode = false;
    if (argc >= 3) {
        const std::string mode = argv[2];
        if (mode == "full") {
            full_mode = true;
        } else if (mode == "quick") {
            full_mode = false;
        } else {
            std::cerr << "Unknown validation mode: " << mode << "\n";
            std::cerr << "Allowed modes: quick (default), full\n";
            return 1;
        }
    }

    const std::string   file_path   = argv[1];
    hid_t               file        = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( file < 0 )
    {
        std::cerr << "Could not open database file: " << file_path << "\n";
        return 1;
    }

    try 
    {
        TripleHankelIO options;
        options.hkind0          = read_attribute<int>(file, "hkind0", H5T_NATIVE_INT);
        options.hkind1          = read_attribute<int>(file, "hkind1", H5T_NATIVE_INT);
        options.r_min           = static_cast<cusfloat>(read_attribute<double>(file, "r_min", H5T_NATIVE_DOUBLE));
        options.rtol            = static_cast<cusfloat>(read_attribute<double>(file, "rtol", H5T_NATIVE_DOUBLE));
        options.atol            = static_cast<cusfloat>(read_attribute<double>(file, "atol", H5T_NATIVE_DOUBLE));
        options.max_segments    = read_attribute<int>(file, "max_segments", H5T_NATIVE_INT);
        options.min_segments    = read_attribute<int>(file, "min_segments", H5T_NATIVE_INT);
        options.segment_cycles  = static_cast<cusfloat>(read_attribute<double>(file, "segment_cycles", H5T_NATIVE_DOUBLE));
        options.rotated_tail_cycles = static_cast<cusfloat>(
            read_attribute_or_default<double>(
                file,
                "rotated_tail_cycles",
                H5T_NATIVE_DOUBLE,
                static_cast<double>(options.rotated_tail_cycles)));
        options.verbose = false;

        std::vector<int> orders_n1 = read_dataset_1d<int>(file, "orders_n1", H5T_NATIVE_INT);
        std::vector<int> orders_n2 = read_dataset_1d<int>(file, "orders_n2", H5T_NATIVE_INT);
        std::vector<int> orders_n3 = read_dataset_1d<int>(file, "orders_n3", H5T_NATIVE_INT);
        std::vector<double> alpha = read_dataset_1d<double>(file, "alpha", H5T_NATIVE_DOUBLE);
        std::vector<double> beta = read_dataset_1d<double>(file, "beta", H5T_NATIVE_DOUBLE);
        std::vector<double> gamma = read_dataset_1d<double>(file, "gamma", H5T_NATIVE_DOUBLE);

        std::array<hsize_t, 6> dims_real{0, 0, 0, 0, 0, 0};
        std::array<hsize_t, 6> dims_imag{0, 0, 0, 0, 0, 0};
        std::vector<double> values_mod_real = read_dataset_6d<double>(file, "values_mod_real", H5T_NATIVE_DOUBLE, dims_real);
        std::vector<double> values_mod_imag = read_dataset_6d<double>(file, "values_mod_imag", H5T_NATIVE_DOUBLE, dims_imag);

        if (dims_real != dims_imag) {
            throw std::runtime_error("values_mod_real and values_mod_imag shapes differ.");
        }

        if (
            static_cast<size_t>(dims_real[0]) != orders_n1.size() ||
            static_cast<size_t>(dims_real[1]) != orders_n2.size() ||
            static_cast<size_t>(dims_real[2]) != orders_n3.size() ||
            static_cast<size_t>(dims_real[3]) != alpha.size() ||
            static_cast<size_t>(dims_real[4]) != beta.size() ||
            static_cast<size_t>(dims_real[5]) != gamma.size()) {
            throw std::runtime_error("Dataset dimensions do not match coordinate arrays.");
        }

        const std::vector<hsize_t> i0_vals = build_validation_indices(dims_real[0], full_mode);
        const std::vector<hsize_t> i1_vals = build_validation_indices(dims_real[1], full_mode);
        const std::vector<hsize_t> i2_vals = build_validation_indices(dims_real[2], full_mode);
        const std::vector<hsize_t> i3_vals = build_validation_indices(dims_real[3], full_mode);
        const std::vector<hsize_t> i4_vals = build_validation_indices(dims_real[4], full_mode);
        const std::vector<hsize_t> i5_vals = build_validation_indices(dims_real[5], full_mode);

        const double abs_tol = 1e-6;
        const double rel_tol = 5e-2;
        const size_t total_points =
            i0_vals.size() * i1_vals.size() * i2_vals.size() * i3_vals.size() * i4_vals.size() * i5_vals.size();
        size_t fail_count = 0;
        double max_abs_diff = 0.0;
        double max_rel_diff = 0.0;
        std::array<hsize_t, 6> worst_idx{0, 0, 0, 0, 0, 0};
        cuscomplex worst_ref(0.0, 0.0);
        cuscomplex worst_got(0.0, 0.0);

        size_t done = 0;
        std::cout << "Validation mode: " << (full_mode ? "full" : "quick") << "\n";
        for (hsize_t i0 : i0_vals) {
            for (hsize_t i1 : i1_vals) {
                for (hsize_t i2 : i2_vals) {
                    for (hsize_t i3 : i3_vals) {
                        for (hsize_t i4 : i4_vals) {
                            for (hsize_t i5 : i5_vals) {
                                const std::array<hsize_t, 6> idx{i0, i1, i2, i3, i4, i5};
                                const size_t flat = flatten_index(dims_real, idx);

                                const int n1 = orders_n1[static_cast<size_t>(i0)];
                                const int n2 = orders_n2[static_cast<size_t>(i1)];
                                const int n3 = orders_n3[static_cast<size_t>(i2)];
                                const cusfloat a = static_cast<cusfloat>(alpha[static_cast<size_t>(i3)]);
                                const cusfloat b = static_cast<cusfloat>(beta[static_cast<size_t>(i4)]);
                                const cusfloat c = static_cast<cusfloat>(gamma[static_cast<size_t>(i5)]);

                                const cuscomplex ref(values_mod_real[flat], values_mod_imag[flat]);
                                const cuscomplex got = integrate_triple_hankel_mod(n1, n2, n3, a, b, c, options);

                                const double abs_diff = std::abs(ref - got);
                                const double rel_diff = abs_diff / std::max(static_cast<double>(std::abs(ref)), 1e-14);

                                if (abs_diff > max_abs_diff) {
                                    max_abs_diff = abs_diff;
                                    max_rel_diff = rel_diff;
                                    worst_idx = idx;
                                    worst_ref = ref;
                                    worst_got = got;
                                }

                                if (abs_diff > abs_tol + rel_tol * std::abs(ref)) {
                                    ++fail_count;
                                }

                                ++done;
                                if (done % 100 == 0 || done == total_points) {
                                    std::cout << "Validation progress: " << done << "/" << total_points << "\r";
                                    std::cout.flush();
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout << "\nValidated " << total_points << " database points in "
              << (full_mode ? "full" : "quick") << " mode.\n";
        std::cout << "Max abs diff: " << max_abs_diff << "\n";
        std::cout << "Max rel diff: " << max_rel_diff << "\n";

        if (fail_count > 0) {
            std::cerr << "Validation failures: " << fail_count << "/" << total_points << "\n";
            std::cerr << "Worst mismatch at index ["
                      << worst_idx[0] << ", " << worst_idx[1] << ", " << worst_idx[2] << ", "
                      << worst_idx[3] << ", " << worst_idx[4] << ", " << worst_idx[5] << "]\n";
            std::cerr << "Reference: " << worst_ref << "\n";
            std::cerr << "Computed : " << worst_got << "\n";
            throw std::runtime_error(std::string("Triple Hankel database validation failed in ") + (full_mode ? "full" : "quick") + " mode");
        }

        H5Fclose(file);
        return 0;

    } 
    catch ( const std::exception& e ) 
    {
        std::cerr << e.what() << "\n";
        H5Fclose(file);
        return 1;

    }
}
