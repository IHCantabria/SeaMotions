
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
#include "field_mesh_data_config.hpp"


FieldMeshDataConfig::FieldMeshDataConfig( 
                                            FieldPointsDef* field_points_def
                                        )
{
    this->mesh_file_path  = field_points_def->mesh_finame;
    this->body_name       = field_points_def->mesh_body_name;
    this->out_components  = field_points_def->out_components;
    this->out_potential   = field_points_def->out_potential;
    this->out_pressure    = field_points_def->out_pressure;
    this->out_velocity    = field_points_def->out_velocity;
}