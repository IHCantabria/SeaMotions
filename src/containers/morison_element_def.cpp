
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


// Include local modules
#include "morison_element_def.hpp"


MorisonElementDef::MorisonElementDef( 
                                        std::ifstream& file,
                                        cusfloat*      cog 
                                    )
{
    // Storage of the center of gravity of the element
    this->cog[0] = cog[0];
    this->cog[1] = cog[1];
    this->cog[2] = cog[2];
    
    // Parse input file to fill the rest of the element properties
    this->_parse_input_file( file );

}


void MorisonElementDef::_parse_input_file( 
                                                std::ifstream& file 
                                            )
{
    // Read start position of the element
    file >> this->start_pos[0] >> this->start_pos[1] >> this->start_pos[2];

    // Read azimuth angle of the element
    file >> this->azimuth_angle;

    // Read declination angle of the element
    file >> this->declination_angle;

    // Read rotation angle along the direction axis
    file >> this->rotation_angle;

    // Read length of the element
    file >> this->length;

    // Read number of divisons of the element
    file >> this->divisions_np;

    // Read Keulegan-Carpenter characteristic length
    file >> this->kc_length;

    // Read axial area
    file >> this->area_axial;

    // Read transversal area
    file >> this->area_transversal;

    // Read vertical area
    file >> this->area_vertical;

    // Read axis added mass coefficient file name
    file >> this->axis_added_mass_coeff_file;

    // Read transversal added mass coefficient file name
    file >> this->transversal_added_mass_coeff_file;

    // Read vertical added mass coefficient file name
    file >> this->vertical_added_mass_coeff_file;

    // Read axis drag coefficient file name
    file >> this->axis_drag_coeff_file;

    // Read transversal drag coefficient file name
    file >> this->transversal_drag_coeff_file;

    // Read vertical drag coefficient file name
    file >> this->vertical_drag_coeff_file;

}