
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


template<typename T>
void read_binary( std::ifstream &in, T &value )
{
    in.read( reinterpret_cast<char*>( &value ), sizeof( T ) );
}

template<typename T>
void read_binary_array( std::ifstream &in, T *data, size_t count ) 
{
    in.read( reinterpret_cast<char*>( data ), sizeof( T ) * count );
}