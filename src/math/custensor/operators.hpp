
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


namespace cut // Namespace Custom Tensor
{
    struct AddOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        { 
            return a + b;
        }
    };


    struct SubOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        {
            return a - b;
        }
    };


    struct MultOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        {
            return a * b;
        }
    };


    struct DivOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        {
            return a / b;
        }
    };


    struct ExpOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a )
        {
            return std::exp( a );
        }
    };
}