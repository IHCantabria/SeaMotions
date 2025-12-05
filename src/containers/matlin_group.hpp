
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
#include <fstream>
#include <string>

// Include local modules
#include "../config.hpp"
#include "../math/math_tools.hpp"
#include "../tools.hpp"


template <typename T>
struct MatLinGroup
{
private:
    // Define class private attributes
    cusfloat    _cog[3]             = { 0.0, 0.0, 0.0 };
    int         _dims_np            = 3;

    // Define class private methods
    void    _load_field_points(
                                    std::string fipath
                                );

public:
    // Define class public attributes
    cusfloat*   cog_to_field_points = nullptr;
    int         end_col             = 0;
    int         end_row             = 0;
    cusfloat*   field_points        = nullptr;
    int*        field_points_cnp    = nullptr;
    int         field_points_nb     = 0;
    int         field_points_np     = 0;
    T*          field_values        = nullptr;
    int         field_values_np     = 0;
    int         fields_np           = 0;
    bool        is_sysmat_field     = false;
    int         start_col           = 0;
    int         start_row           = 0;
    T*          sysmat              = nullptr;
    int         sysmat_ncols        = 0;
    int         sysmat_nrows        = 0;
    T*          sysmat_steady       = nullptr;

    // Define class constructors and destructor
    MatLinGroup(
                    int         sysmat_nrows_in,
                    int         sysmat_ncols_in,
                    int         field_points_nb_in,
                    int         fields_np_in,
                    int         start_row_in,
                    int         end_row_in,
                    int         start_col_in,
                    int         end_col_in,
                    bool        is_sysmat_field_in
                );

    MatLinGroup(
                    std::string fp_path_in,
                    cusfloat*   cog,
                    int         sysmat_ncols_in,
                    int         fields_np_in,
                    int         start_col_in,
                    int         end_col_in,
                    bool        is_sysmat_field_in
                );

    ~MatLinGroup(
                    void
                );

    // Define class methods
    void clear_field_values( void );

    void clear_sysmat( void );

};


template<typename T>
void MatLinGroup<T>::clear_field_values(
                                            void
                                        )
{
    clear_vector( this->sysmat_nrows * this->fields_np, this->field_values );
}


template<typename T>
void MatLinGroup<T>::clear_sysmat(
                                            void
                                        )
{
    clear_vector( this->field_values_np, this->field_values );
    clear_vector( this->sysmat_nrows * this->sysmat_ncols, this->sysmat );
}




template<typename T>
void MatLinGroup<T>::_load_field_points(
                                            std::string fipath
                                        )
{
    // Open file unit
    std::ifstream infile( fipath );
    CHECK_FILE_UNIT_STATUS( infile, fipath );

    // Loop over lines to count the number of points
    int         count_lines = 0;
    std::string line( "" );
    while ( getline( infile, line ) )
    {
        count_lines++;
    }

    this->field_points_np = count_lines;

    // Rewind file
    infile.clear( );
    infile.seekg( 0 );

    // Allocate space to storage field points
    this->cog_to_field_points   = generate_empty_vector<cusfloat>( this->_dims_np * this->field_points_np );
    this->field_points          = generate_empty_vector<cusfloat>( this->_dims_np * this->field_points_np );

    // Loop over lines to get the field points position
    for ( int i=0; i<this->field_points_np; i++ )
    {
        // Get field points position
        infile >> this->field_points[3*i];
        infile >> this->field_points[3*i+1];
        infile >> this->field_points[3*i+2];

        // Calculate distance from the cog the ith field point
        sv_sub( 
                    3,
                    &(this->field_points[3*i]),
                    this->_cog,
                    &(this->cog_to_field_points[3*i])
                );
    }

    // Close file unit
    infile.close( );
}


template <typename T>
MatLinGroup<T>::MatLinGroup(
                                int     sysmat_nrows_in,
                                int     sysmat_ncols_in,
                                int     field_points_nb_in,
                                int     fields_np_in,
                                int     start_row_in,
                                int     end_row_in,
                                int     start_col_in,
                                int     end_col_in,
                                bool    is_sysmat_field_in
                            )
{
    // Storage input arguments
    this->end_col           = end_col_in;
    this->end_row           = end_row_in;
    this->field_points_nb   = field_points_nb_in;
    this->field_points_np   = sysmat_nrows_in;
    this->fields_np         = fields_np_in;
    this->is_sysmat_field   = is_sysmat_field_in;
    this->start_col         = start_col_in;
    this->start_row         = start_row_in;
    this->sysmat_nrows      = sysmat_nrows_in;
    this->sysmat_ncols      = sysmat_ncols_in;

    // Define size dependent constant attributes
    this->field_values_np   = this->fields_np * this->sysmat_nrows;

    // Allocate space for the system matrixes
    this->sysmat            = generate_empty_vector<T>( this->sysmat_nrows * this->sysmat_ncols );
    this->sysmat_steady     = generate_empty_vector<T>( this->sysmat_nrows * this->sysmat_ncols );

    // Allocate space for the field points and the field values
    this->cog_to_field_points   = generate_empty_vector<cusfloat>( this->_dims_np * this->sysmat_nrows );
    this->field_points          = generate_empty_vector<cusfloat>( this->_dims_np * this->sysmat_nrows );
    this->field_points_cnp      = generate_empty_vector<int>( this->field_points_nb+1 );
    this->field_values          = generate_empty_vector<T>( this->field_values_np );
}


template<typename T>
MatLinGroup<T>::MatLinGroup(
                                std::string fp_path_in,
                                cusfloat*   cog_in,
                                int         sysmat_ncols_in,
                                int         fields_np_in,
                                int         start_col_in,
                                int         end_col_in,
                                bool        is_sysmat_field_in
                            )
{
    // Storage input arguments
    this->end_col           = end_col_in;
    this->is_sysmat_field   = is_sysmat_field_in;
    this->fields_np         = fields_np_in;
    this->start_col         = start_col_in;
    this->sysmat_ncols      = sysmat_ncols_in;

    copy_vector( 3, cog_in, this->_cog );

    // Read data from file
    this->_load_field_points( fp_path_in );

    // Define the remaining of the matrix dimensions
    this->end_row           = this->field_points_np;
    this->field_points_nb   = 1;
    this->start_row         = 0;
    this->sysmat_nrows      = this->field_points_np;

    // Allocate space for the system matrixes
    this->sysmat            = generate_empty_vector<T>( this->sysmat_nrows * this->sysmat_ncols );
    this->sysmat_steady     = generate_empty_vector<T>( this->sysmat_nrows * this->sysmat_ncols );

    // Allocate space for the field points and the field values
    this->field_points_cnp      = generate_empty_vector<int>( this->field_points_nb+1 );
    this->field_values          = generate_empty_vector<T>( this->fields_np * this->sysmat_nrows );
}


template <typename T>
MatLinGroup<T>::~MatLinGroup(
                            void
                        )
{
    mkl_free( this->cog_to_field_points );
    mkl_free( this->field_points );
    mkl_free( this->field_points_cnp );
    mkl_free( this->field_values );
    mkl_free( this->sysmat );
    mkl_free( this->sysmat_steady );
}


// Define short types to refer in a handly way to the MatLinGroup family members
typedef     MatLinGroup<cusfloat>       MLGFloat;
typedef     MatLinGroup<cuscomplex>     MLGCmpx;
