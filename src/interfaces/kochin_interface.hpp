
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
#include "../config.hpp"
#include "../containers/source_node.hpp"


struct KochinInterface
{
private:
    // Define class private attributes
    cuscomplex  _cf_l           = cuscomplex( 0.0, 0.0 );
    cusfloat    _ep_l           = 0.0;
    bool        _is_cos         = false;
    int         _l_order        = 0;
    SourceNode* _source         = nullptr;
    cusfloat    _water_depth    = 0.0;
    cusfloat    _wave_num       = 0.0;

public:
    // Defincde class constructors and destructor
    KochinInterface(
                        SourceNode* source,
                        cusfloat    wave_number,
                        cusfloat    water_depth,
                        int         l_order,
                        bool        is_cos
                    );

    // Define class public attributes
    cuscomplex  operator()( 
                                cusfloat    xi,
                                cusfloat    eta,
                                cusfloat    x,
                                cusfloat    y,
                                cusfloat    z
                        );
            
    void        set_is_cos(
                                bool        is_cos
                            );

    void        set_l_order( 
                                int         l_order 
                           );

    void        set_source(
                               SourceNode* source
                           );

};
