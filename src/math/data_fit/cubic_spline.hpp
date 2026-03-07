
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

// Include local modules
#include "data_fit_error.hpp"
#include "../math_interface.hpp"
#include "../math_tools.hpp"


template<typename T>
struct CubicSplineInterp
{
    // Define class attributes
    T*              _spline_coeffs      = nullptr;  // Coefficients of the spline generated
    int             _spline_coeffs_np   = 0;        // Number of coefficients of the spline generated
    DFTaskPtr*      _task               = nullptr;  // Data Fitting task descriptor


    // Define class contructos and destructors
    CubicSplineInterp(
                        int     data_np, // Number of sampling points in the data
                        T* x_data,
                        T* y_data
                    )
    {
        // Create local variable to handle mkl routines status
        int errcode;

        // Calculate dimensions of spline coefficients
        this->_spline_coeffs_np = (data_np-1)*DF_PP_CUBIC;
        this->_spline_coeffs    = generate_empty_vector<T>(this->_spline_coeffs_np);

        ///////////////////////////////////////////////////////////////////////
        /***************** Create task for 1D interpolation ******************/
        ///////////////////////////////////////////////////////////////////////
        this->_task = new DFTaskPtr;
        errcode = dfNewTask1D<T>(
                                    this->_task,
                                    data_np, 
                                    x_data,
                                    DF_NON_UNIFORM_PARTITION,
                                    1,
                                    y_data,
                                    DF_NO_HINT 
                                );
        DF_CHECK_STATUS(errcode);

        /////////////////////////////////////////////////////////////////////////
        /***** Edit task parameters for natural cubic spline construction ******/
        /////////////////////////////////////////////////////////////////////////
        errcode = dfEditPPSpline1D<T>( 
                                        *(this->_task),
                                        DF_PP_CUBIC,
                                        DF_PP_NATURAL,
                                        DF_BC_NOT_A_KNOT,
                                        nullptr,
                                        DF_NO_IC,
                                        0,
                                        this->_spline_coeffs,
                                        DF_NO_HINT
                                    );
        DF_CHECK_STATUS(errcode);

        /////////////////////////////////////////////////////////////////////////
        /**********  Construct natural cubic spline using STD method ***********/
        /////////////////////////////////////////////////////////////////////////
        errcode = dfConstruct1D(
                                    *(this->_task),
                                    DF_PP_SPLINE,
                                    DF_METHOD_STD
                                );
        DF_CHECK_STATUS(errcode);

    }

    CubicSplineInterp( CubicSplineInterp& )
    {
        std::cerr << "It is not allowed make copies of the struct CubicSplineInterp" << std::endl;
        exit(800);
    }

    ~CubicSplineInterp( void )
    {
        // Delete task object
        delete this->_task;

        // Delete heap memory associated to vectors/matrixes allocated along 
        // the object
        mkl_free(this->_spline_coeffs);
    }

    // Define class methods
    T operator()( T x_new )  const
    { 
        T x_new_inp[1] = {x_new};
        int ndorder = 1;
        int dorder[1] = {0};
        T y_new = 0;
        int errcode = dfInterpolate1D(
                                        *(this->_task),
                                        DF_INTERP,
                                        DF_METHOD_PP,
                                        1,              // Number of interpolation points
                                        x_new_inp,     // Pointer      
                                        DF_NON_UNIFORM_PARTITION,
                                        ndorder,
                                        dorder,
                                        0,
                                        &y_new,
                                        DF_MATRIX_STORAGE_ROWS,
                                        0
                                        );
        DF_CHECK_STATUS(errcode);

        return y_new;
    }

    void operator()( int points_np, T* x_new, T* y_new ) const
    { 
        int ndorder = 1;
        int dorder[1] = {0};
        int errcode = dfInterpolate1D(
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
};
