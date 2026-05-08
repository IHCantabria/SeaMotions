
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

// Include general usage scientific libraries
#include "mkl.h"
#include <cstddef>
#include <sstream>
#include <stdexcept>

// Include local modules
#include "../../config.hpp"
#include "data_fit_error.hpp"
#include "../math_interface.hpp"
#include "../math_tools.hpp"


template<typename T>
struct Interp1D
{
private:
    // Define class attributes
    T*          _spline_coeffs      = nullptr;  // Coefficients of the spline generated
    int         _spline_coeffs_np   = 0;        // Number of coefficients of the spline generated
    int         _data_np            = 0;        // Number point of interpolation points
    T           _max_x              = 0.0;      // Maximum x value of the input series
    T           _max_y              = 0.0;      // Maximum y value of the input series
    T           _min_x              = 0.0;      // Minimum x value of the input series
    T           _min_y              = 0.0;      // Minimum y value of the input series
    DFTaskPtr*  _task               = nullptr;  // Interpolation task object

    // Define class methods
    void _check_x_monotonic_ascendent( int data_np, T* x )
    {
        for ( int i=1; i<data_np; i++ )
        {
            if ( x[i] < x[i-1] )
            {
                _THROW_CLASS_ERROR( "X data series is not monotonically ascendent." );
            }
        }
    }

    void _get_data_bounding_box( int data_np, T* x_data, T* y_data )
    {
        // Get X maximum and minimum values. This field should be ordered monotonically
        // ascendent to the mininm and maximum values should be the first and the last one 
        // respectivelly.
        this->_max_x = x_data[data_np-1];
        this->_min_x = x_data[0];

        // Loop over data to get the Y minimum and maximum values
        this->_max_y = -1e302;
        this->_min_y = 1e302;
        for ( int i=0; i<data_np; i++ )
        {
            // Check for Y minimum values
            if ( y_data[i] < this->_min_y )
            {
                this->_min_y = y_data[i];
            }

            // Check for Y maximum values
            if ( y_data[i] > this->_max_y )
            {
                this->_max_y = y_data[i];
            }
        }
    }

public:
    Interp1D( 
                int data_np, 
                T*  x_data, 
                T*  y_data 
            )
    {
        /**
         * @brief Linear Piecewise interpolation 1D
         * 
         * \param   data_np     Number of data points
         * \param   x_data      X data values
         * \param   y_data      Y data values
         * 
         */


        // Storage the required input arguments into the class attributes
        this->_data_np = data_np;


        ///////////////////////////////////////////////////////////////////////
        /******************* Perform initial data checks *********************/
        ///////////////////////////////////////////////////////////////////////
        this->_check_x_monotonic_ascendent( data_np, x_data );

        // Get data bouding box
        this->_get_data_bounding_box( data_np, x_data, y_data );

        ///////////////////////////////////////////////////////////////////////
        /*************** Allocate requested memory locations *****************/
        ///////////////////////////////////////////////////////////////////////
        // Allocate space for the linear spline coefficients
        this->_spline_coeffs_np = ( this->_data_np - 1 ) * DF_PP_LINEAR;
        this->_spline_coeffs    = generate_empty_vector<T>( this->_spline_coeffs_np );

        // Create local variable to handle mkl routines status
        int errcode;

        ///////////////////////////////////////////////////////////////////////
        /***************** Create task for 1D interpolation ******************/
        ///////////////////////////////////////////////////////////////////////
        this->_task = new DFTaskPtr;
        errcode     = dfNewTask1D<T>( 
                                        this->_task, 
                                        data_np, 
                                        x_data, 
                                        DF_NON_UNIFORM_PARTITION, 
                                        1, 
                                        y_data, 
                                        DF_NO_HINT
                                    );
        DF_CHECK_STATUS( errcode );

        /////////////////////////////////////////////////////////////////////////
        /***** Edit task parameters for natural cubic spline construction ******/
        /////////////////////////////////////////////////////////////////////////
        errcode = dfEditPPSpline1D<T>( 
                                        *(this->_task), 
                                        DF_PP_LINEAR, 
                                        DF_PP_DEFAULT, 
                                        DF_NO_BC, 
                                        0,
                                        DF_NO_IC, 
                                        0, 
                                        this->_spline_coeffs, 
                                        DF_NO_HINT 
                                    );
        DF_CHECK_STATUS( errcode );

        /////////////////////////////////////////////////////////////////////////
        /**********  Construct natural cubic spline using STD method ***********/
        /////////////////////////////////////////////////////////////////////////
        errcode =  dfConstruct1D<T>( 
                                        *(this->_task), 
                                        DF_PP_SPLINE, 
                                        DF_METHOD_STD 
                                    );
        DF_CHECK_STATUS( errcode );

    }

    ~Interp1D( void )
    {
        // Delete data fitting task
        int errcode = dfDeleteTask( this->_task );
        DF_CHECK_STATUS( errcode );
        delete this->_task;

        // Delete spline coefficients memory
        mkl_free( this->_spline_coeffs );

    }

    // Define class methods
    T get_x_max( void )
    {
        return this->_max_x;
    }

    T operator()( T x_new ) const
    {
        // Configure input for the interpolation
        int     ndorder     = 1;
        int     dorder[1]   = { 0 };
        T  x_inp[1]    = { x_new };
        T  y_res       = 0.0;

        // Perform interpolation
        int errcode =   dfInterpolate1D<T>( 
                                            *(this->_task),
                                            DF_INTERP,
                                            DF_METHOD_PP,
                                            1, 
                                            x_inp, 
                                            DF_NON_UNIFORM_PARTITION, 
                                            ndorder,
                                            dorder, 
                                            0, 
                                            &y_res, 
                                            DF_MATRIX_STORAGE_ROWS, 
                                            0
                                        );
        DF_CHECK_STATUS( errcode );


        return y_res;
    }

    void operator()( int points_np, T* x_new, T* y_new ) const
    { 
        int ndorder = 1;
        int dorder[1] = {0};
        int errcode = dfInterpolate1D<T>(
                                            *(this->_task),
                                            DF_INTERP,
                                            DF_METHOD_PP,
                                            points_np,              // Number of interpolation points
                                            x_new,                  // Pointer      
                                            DF_NON_UNIFORM_PARTITION,
                                            ndorder,
                                            dorder,
                                            0,
                                            y_new,
                                            DF_MATRIX_STORAGE_ROWS,
                                            0
                                        );
        DF_CHECK_STATUS(errcode);

    }

    std::size_t memory_bytes( void ) const
    {
        return static_cast<std::size_t>( this->_spline_coeffs_np ) * sizeof( T );
    }

};
