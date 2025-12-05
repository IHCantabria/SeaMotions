
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
#include <iostream>

// Include local modules
#include "shape_functions.hpp"


int     dofs_rectangular_region(
                                    int poly_order
                                )
{
    int dofs_total = 0;
    if ( poly_order == 0 )
    {
        dofs_total = 1;
    }
    else if ( poly_order > 0 )
    {
        dofs_total = ( poly_order + 1 ) * ( poly_order + 1 );
    }
    else
    {
        std::cerr << "ERROR - dofs_rectangular_region" << std::endl;
        std::cerr << "Polynomial order: " << poly_order << " not valid." << std::endl;
        throw std::runtime_error( "" );
    }

    return dofs_total;
}


int     dofs_triangular_region(
                                    int poly_order
                                )
{
    int dofs_total = 0;
    if ( poly_order == 0 )
    {
        dofs_total = 1;
    }
    else if ( poly_order > 0 )
    {
        for ( int i=1; i<poly_order+2; i++ )
        {
            dofs_total += i;
        }
    }
    else
    {
        std::cerr << "ERROR - dofs_rectangular_region" << std::endl;
        std::cerr << "Polynomial order: " << poly_order << " not valid." << std::endl;
        throw std::runtime_error( "" );
    }

    return dofs_total;
}


cusfloat shape_functions(
                            int p_order,
                            int q_order,
                            cusfloat ,
                            cusfloat 
                        )
{
    cusfloat val = 0.0;
    if ( p_order == 0 && q_order == 0 )
    {
        val = 1.0;
    }
    else
    {
        std::cerr << "ERROR - shape_functions" << std::endl;
        std::cerr << "Not implemented yet!" << std::endl;
        throw std::runtime_error( "" );
    }

    return val;
}