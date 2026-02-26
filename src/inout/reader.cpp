
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

// Include general usage libraries
#include <fstream>
#include <sstream>
#include <string>

// Include local modules
#include "reader.hpp"


std::string read_channel_name( std::ifstream& infile )
{
    // Generate local auxiliar variables
    std::string aux_str, line, channel_name;

    // Read line from file unit
    std::getline(infile, line);

    // Get signal value and name from the line
    // read
    std::istringstream iss(line);
    iss >> aux_str >> channel_name >> aux_str;

    return channel_name;
}


void skip_header( 
                    std::ifstream&  infile, 
                    int&            line_count, 
                    int             np 
                )
{
    std::string line;
    for ( int i=0; i<np; i++)
    {
        std::getline( infile,  line );
    }
    line_count += np;
}