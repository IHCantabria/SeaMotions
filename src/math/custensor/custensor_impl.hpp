
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

// Include local modules
#include "base.hpp"
#include "binary_expression.hpp"
#include "../math_tools.hpp"


namespace cut
{
    /********************************************/
    /************** MODULE MACROS ***************/
    /********************************************/
    #define CUSTENSOR_OPERATOR( OPERATOR, OP_TYPE )                                                                     \
    template<typename R, typename = typename std::enable_if<!std::is_arithmetic<R>::value>::type>                       \
    BinaryExpression<CusTensor, R, T, OP_TYPE> operator OPERATOR ( const R& other )                                     \
    {                                                                                                                   \
        return BinaryExpression<CusTensor, R, T, OP_TYPE>( *this, other );                                              \
    }                                                                                                                   \
                                                                                                                        \
    template<typename R, typename = typename std::enable_if<std::is_arithmetic<R>::value>::type>                        \
    BinaryExpression<CusTensor, ScalarExpression<T>, T, OP_TYPE> operator OPERATOR ( const R& other )                   \
    {                                                                                                                   \
        return BinaryExpression<CusTensor, ScalarExpression<T>, T, OP_TYPE>( *this, ScalarExpression( other ) );        \
    }                                                                                                                   \

    /********************************************/
    /************* Custensor Class **************/
    /********************************************/
    template<typename T>
    class CusTensor
    {
    private:
        // Define private variables
        T*          _data       = nullptr;
        bool        _is_heap    = false;
        vector_st   _shape;
        size_t      _size       = 0;
        vector_st   _strides;

        // Define private functions
        void _allocate_memory( void )
        {
            this->_data     = generate_empty_vector<T>( this->_size );
            this->_is_heap  = true;
        }

        void _compute_strides( void )
        {
            // Resize strides vector to allocate tensor dimensions
            this->_strides.resize( this->_shape.size( ) );

            // Pre-compute tensor index strides to save computational
            // time when accessing tensor indexes
            // --> Row Major Scheme is assumed
            size_t stride = 1;
            for ( int i = static_cast<int>( _shape.size( ) ) - 1; i >= 0; --i )
            {
                this->_strides[i]   = stride;
                stride              *= _shape[i];
            }
        }

    public:
        /******** Define class constructor *******/
        CusTensor( ) = default;

        CusTensor( int N )
        {
            // Storage input arguments
            this->_shape    = { static_cast<size_t>( N ) };
            this->_size     = N;

            // Compute strides
            this->_compute_strides( );

            // Allocate memory space
            this->_allocate_memory( );
        }

        CusTensor( size_t N )
        {
            // Storage input arguments
            this->_shape    = { N };
            this->_size     = N;

            // Compute strides
            this->_compute_strides( );

            // Allocate memory space
            this->_allocate_memory( );
        }

        CusTensor( const vector_st& shape )
        {
            // Storage input arguments
            this->_shape = shape;
            
            // Compute matrix size
            this->_size = 1;
            for ( size_t s : shape ) this->_size *= s;

            // Compute tensor strides
            this->_compute_strides( );

            // Allocate memory to stroage the data
            this->_allocate_memory( );
        }

        // Move constructor
        CusTensor( CusTensor&& other ) noexcept
        {
            // Move data into current class
            this->_data( other._data );
            this->_shape( std::move( other._shape ) );
            this->_strides( std::move( other._strides ) );
            this->_size( other._size );
            this->_is_heap( other._is_heap );

            // Unlink data from other class and reset it
            other._data     = nullptr;
            other._is_heap  = false;
            other._size     = 0;
        }

        // Copy assignment operator
        CusTensor& operator=( const CusTensor& other )
        {
            if ( this != &other ) 
            {
                if (_is_heap)
                {
                    mkl_free( this->_data );
                }

                this->_shape    = other._shape;
                this->_strides  = other._strides;
                this->_size     = other._size;
                this->_allocate_memory( );
                std::copy( other._data, other._data + _size, _data );
            }

            return *this;
        }

        // Move assignement operator
        CusTensor& operator=( CusTensor&& other ) noexcept
        {
            if ( this != &other ) 
            {
                if (_is_heap)
                {
                    mkl_free(_data);
                }

                // Move data into current class
                this->_data     = other._data;
                this->_shape    = std::move(other._shape);
                this->_strides  = std::move(other._strides);
                this->_size     = other._size;
                this->_is_heap  = other._is_heap;

                // Unlink and reset from previous class
                other._data = nullptr;
                other._is_heap = false;
                other._size = 0;

            }

            return *this;
        }

        template<typename Expr>
        CusTensor( const Expr& expr )
        {
            this->operator=( expr );
        }

        template<>
        CusTensor( const CusTensor& other )
        {
            this->_operator=( other );
        }

        ~CusTensor( )
        {
            if ( this->_is_heap )
            {
                mkl_free( this->_data );
            }
        }

        // Define class methods
        vector_st shape( void ) const
        {
            return this->_shape;
        }

        size_t  size( void ) const
        {
            return this->_size;
        }

        // Define class operators
        T& operator [ ] ( size_t i )
        {
            return this->_data[i];
        }

        const T operator [ ] ( size_t i ) const 
        {
            return this->_data[i];
        }

        template<typename... Args>
        T& operator( )( Args... args )
        {
            std::array<size_t, sizeof...( Args )> idx{ static_cast<size_t>( args )... };
            assert( idx.size( ) == _shape.size( ) );
            size_t offset = 0; 
            for ( size_t i = 0; i < _shape.size( ); i++ )
            {
                offset += idx[i] * _strides[i];
            }

            return _data[offset];
        }

        template<typename... Args>
        const T operator( )(Args... args) const 
        {
            std::array<size_t, sizeof...( Args )> idx{ static_cast<size_t>( args )... };
            assert( idx.size( ) == this->_shape.size( ) );
            size_t offset = 0;
            for ( size_t i = 0; i < this->_shape.size( ); i++ )
            {
                offset += idx[i] * this->_strides[i];
            }

            return this->_data[offset];
        }

        template<typename Expr>
        CusTensor& operator=( const Expr& expr )
        {
            // Check for size and shape inconsistency
            if ( this->_is_heap )
            {
                if ( this->_shape != expr.shape( ) )
                {
                    CUSTENSOR_SHAPE_MISMATCH( this->_shape, expr.shape( ) ) 
                }
            }
            else
            {
                this->_size     = expr.size( );
                this->_shape    = expr.shape( );
                this->_compute_strides( );
                this->_allocate_memory( );
            }

            // Evaluate expression
            expr.evaluate( this->_data );

            return *this;
        }

        template<>
        CusTensor& operator=( const CusTensor& other )
        {
            // Check for size and shape inconsistency
            if ( this->_is_heap )
            {
                if ( this->_shape != other.shape( ) )
                {
                    CUSTENSOR_SHAPE_MISMATCH( this->_shape, other.shape( ) )
                }
            }
            else
            {
                this->_size     = other.size( );
                this->_shape    = other.shape( );
                this->_compute_strides( );
                this->_allocate_memory( this->_size );
            }

            // Copy data
            for ( size_t i=0; i<this->_size; i++ )
            {
                this->_data[i] = other.data[i];
            }

            return *this;
        }

        // Define custensor operators
        CUSTENSOR_OPERATOR( +, AddOp  )
        CUSTENSOR_OPERATOR( -, SubOp  )
        CUSTENSOR_OPERATOR( *, MultOp )
        CUSTENSOR_OPERATOR( /, DivOp  )

        void print( void ) const
        {
            if ( this->_is_heap )
            {
                std::cout << this->_data[0];
                for ( size_t i=1; i<this->_size; i++ )
                {
                    std::cout << " " << this->_data[i];
                }
                std::cout << std::endl;
            }
            else
            {
                std::cout << "WARNING: there is no data allocated to be printed into the screen" << std::endl;
            }
        }

    };
}