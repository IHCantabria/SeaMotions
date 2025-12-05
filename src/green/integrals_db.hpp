
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
#include "pulsating_fin_depth_cheby.hpp"
#include "pulsating_inf_depth_cheby.hpp"


class IntegralsDb
{
private:
    // Declare class methods
    void _error_copy_constructor( void );
    
public:
    // Declare class attributes
    bool is_build = false;
    P3*         l1;
    P3*         l1_da;
    P3*         l1_db;
    P3*         l2;
    P3*         l3;
    P3*         l3_da;
    P3*         l3_db;
    P3*         m1;
    P3*         m1_da;
    P3*         m1_db;
    P3*         m2;
    P3*         m3;
    P3*         m3_da;
    P3*         m3_db;
    R11*        r11;
    R11A_dX*    r11a_dx;
    R11B_dX*    r11b_dx;
    R12*        r12;
    R12_dX*     r12_dx;
    R21*        r21;
    R21_dX*     r21_dx;
    R22*        r22;
    R22_dX*     r22_dx;

    // Declare class constructors and destructor
    IntegralsDb( void );

    IntegralsDb( IntegralsDb& );

    IntegralsDb( IntegralsDb* );

    ~IntegralsDb(void);

    // Declare class methods
    void fold_h(cusfloat H);
    
};


void build_integrals_db(IntegralsDb &idb);
