
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


#pragma once

// Include general usage libraries

// Include local modules
#include "../math_tools.hpp"


struct CSRMatrix
{
    // Define class attributes
    int*        col_index       = nullptr;
    int*        row_index_cum   = nullptr;
    int         rows_np         = 0;
    cusfloat*   values          = nullptr;
    int         values_np       = 0;

    // Define class constructor
    CSRMatrix(int rows_np_in, int values_np_in)
    {
        this->col_index     = generate_empty_vector<int>(values_np_in);
        this->row_index_cum = generate_empty_vector<int>(rows_np_in+1);
        this->rows_np       = rows_np_in;
        this->values        = generate_empty_vector<cusfloat>(values_np_in);
        this->values_np     = values_np_in;
    }

    CSRMatrix(int rows_np_in, cusfloat* dense_mat_in)
    {
        // Storage the number of rows
        this->rows_np = rows_np_in;

        // Get the number of non zero-values per row
        cusfloat ref_value    = 0.0;
        this->row_index_cum = generate_empty_vector<int>(rows_np+1);
        for (int i=0; i<this->rows_np; i++)
        {
            for (int j=0; j<this->rows_np; j++)
            {
                if (assert_scalar_equality(dense_mat_in[i*rows_np+j], ref_value, EPS_PRECISION) == 0)
                {
                    this->row_index_cum[i+1] += 1;
                }
            }
        }
        for (int i=0; i<this->rows_np; i++)
        {
            this->row_index_cum[i+1] = this->row_index_cum[i+1] + this->row_index_cum[i];
        }

        // Get total number of non-zero values
        this->values_np = this->row_index_cum[this->rows_np];

        // Get the column indexes and the non-zero values and  
        // storage them into their proper vectors
        if (this->values_np > 0)
        {
            this->col_index = generate_empty_vector<int>(this->values_np);
            this->values    = generate_empty_vector<cusfloat>(this->values_np);
            int count       = 0;
            for (int i=0; i<this->rows_np; i++)
            {
                for (int j=0; j<this->rows_np; j++)
                {
                    if (assert_scalar_equality(dense_mat_in[i*rows_np+j], ref_value, EPS_PRECISION) == 0)
                    {
                        this->col_index[count]  = j;
                        this->values[count]     = dense_mat_in[i*rows_np+j];
                        count++;
                    }
                }
            }
        }
        else
        {
            // Generate matrixes to storage a single 0 value
            // It is done in this way to prevent undefined behaviour
            // due to the creation of zero length matrixes.
            this->col_index = generate_empty_vector<int>(1);
            this->values    = generate_empty_vector<cusfloat>(1);

            // Put value on first column of the last row
            this->row_index_cum[this->rows_np] = 1;
            
            // Update number of vales variable
            this->values_np = 1;
        }
    }

    CSRMatrix(int rows_np_in, int values_np_in, int* rows_cum, int* col_index, cusfloat* values)
    {
        // Storage the number of rows
        this->rows_np = rows_np_in;

        // Storage the number of non-zero values
        this->values_np = values_np_in;

        // Allocate space for the matrix data
        this->col_index     = generate_empty_vector<int>(this->values_np);
        this->row_index_cum = generate_empty_vector<int>(this->rows_np+1);
        this->values        = generate_empty_vector<cusfloat>(this->values_np);

        // Copy cumulative row values
        for (int i=0; i<this->rows_np+1; i++)
        {
            this->row_index_cum[i] = rows_cum[i];
        }

        // Copy column and values
        for (int i=0; i<this->values_np; i++)
        {
            this->col_index[i]  = col_index[i];
            this->values[i]     = values[i];
        }
    }

    ~CSRMatrix(void)
    {
        mkl_free(this->col_index);
        mkl_free(this->row_index_cum);
        mkl_free(this->values);
    }

    // Define class methods
    void print(void)
    {
        std::cout << "CSRMatrix - Contents:"  << std::endl;
        std::cout << " - Number of Rows: " << this->rows_np << std::endl;
        std::cout << " - Number of Values: " << this->values_np << std::endl;
        std::cout << " - Sparse Non-Zero Values:" << std::endl;
        print_vector(this->values_np, this->values, 0, 6);
        std::cout << " - Column Indexes:" << std::endl;
        print_vector(this->values_np, this->col_index, 0, 6);
        std::cout << " - Cumulative Row Indexes:" << std::endl;
        print_vector(this->rows_np+1, this->row_index_cum, 0, 0);
    }
    
};
