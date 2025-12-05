
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// Include local modules
#include "sparse_tools.hpp"


CSRMatrix* convert_mkl_to_csrmatrix(sparse_matrix_t* mkl_mat)
{
    // Get data from mkl matrix format in csr_matrix format
    // using three vectors: rows_cum, col_indexes, values
    sparse_index_base_t indexing;
    int cols_np = 0;
    int rows_np = 0;
    int* rows_index_cum = nullptr;
    int* rows_index_cum_dummy = nullptr;
    int* col_index = nullptr;
    double* values = nullptr;
    mkl_sparse_d_export_csr(
                            *mkl_mat,                   // Source
                            &indexing,                  // Indexing
                            &rows_np,                   // Rows
                            &cols_np,                   // Columns
                            &rows_index_cum,            // Rows Start
                            &rows_index_cum_dummy,      // Rows End
                            &col_index,                 // col_indx
                            &values                     // values
                            );

    // Create CSRMatrix to storage the data
    CSRMatrix* csr_mat = new CSRMatrix(
                                        rows_np,
                                        rows_index_cum[rows_np],
                                        rows_index_cum,
                                        col_index,
                                        values
                                        );

    return csr_mat;
}