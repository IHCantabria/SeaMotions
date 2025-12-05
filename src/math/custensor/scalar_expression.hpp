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
#include "../../config.hpp"


namespace cut // Namespace Custom Tensor
{
    template<typename T>
    class ScalarExpression
    {
    private:
        T   _value;

    public:
        // Define Constructor
        ScalarExpression( T value_in ): _value( value_in ) { }

        // Define class methods
        vector_st shape( void ) const
        {
            return { 1 };
        }

        size_t  size( void ) const
        { 
            return 1;
        }

        // Define class operators
        const T operator[ ] ( size_t  ) const 
        { 
            return this->_value; 
        }

    };
}