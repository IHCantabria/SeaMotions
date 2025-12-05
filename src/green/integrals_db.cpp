
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
#include <cassert>
#include <exception>

// Include local modules
#include "integrals_db.hpp"
#include "pulsating_fin_depth_cheby.hpp"
#include "pulsating_inf_depth_cheby.hpp"


void IntegralsDb::_error_copy_constructor(
                                            void
                                        )
{
    std::cerr << "Invoking copy constructor for class IntegralsDB" << std::endl;
    throw std::runtime_error( "" );
}


void IntegralsDb::fold_h(cusfloat H)
{
    // Check if the database is already build
    assert(this->is_build && "Not build Integrals database.");

    // Fold all integrals with h
    this->l1->fold_h(H);
    this->l1_da->fold_h(H);
    this->l1_db->fold_h(H);
    this->l2->calculate_h_1D(H);
    this->l3->fold_h(H);
    this->l3_da->fold_h(H);
    this->l3_db->fold_h(H);
    this->m1->fold_h(H);
    this->m1_da->fold_h(H);
    this->m1_db->fold_h(H);
    this->m2->calculate_h_1D(H);
    this->m3->fold_h(H);
    this->m3_da->fold_h(H);
    this->m3_db->fold_h(H);

}


IntegralsDb::IntegralsDb( void )
{
    ////////////////////////////////////////////////
    ////// Load infinite water depth integrals//////
    ////////////////////////////////////////////////

    // Load R11 integral
    this->r11 = new R11();
    set_data_r11(this->r11);

    // Load R11A_dX integral
    this->r11a_dx = new R11A_dX();
    set_data_r11a_dx(this->r11a_dx);

    // Load R11B_dX integral
    this->r11b_dx = new R11B_dX();
    set_data_r11b_dx(this->r11b_dx);

    // Load R12 integral
    this->r12 = new R12();
    set_data_r12(this->r12);

     // Load R12_dX integral
    this->r12_dx = new R12_dX();
    set_data_r12_dx(this->r12_dx);

    // Load R21 integral
    this->r21 = new R21();
    set_data_r21(this->r21);

    // Load R21_dX integral
    this->r21_dx = new R21_dX();
    set_data_r21_dx(this->r21_dx);

    // Load R22 integral
    this->r22 = new R22();
    set_data_r22(this->r22);

    // Load R22_dX integral
    this->r22_dx = new R22_dX();
    set_data_r22_dx(this->r22_dx);

    ////////////////////////////////////////////////
    //////// Load finite water depth integrals//////
    ////////////////////////////////////////////////

    // Load L1 integral
    this->l1 = new P3();
    set_data_l1(this->l1);

    // Load L1_dA integral
    this->l1_da = new P3();
    set_data_l1_da(this->l1_da);

    // Load L1_dB integral
    this->l1_db = new P3();
    set_data_l1_db(this->l1_db);

    // Load L2 integral
    this->l2 = new P3();
    set_data_l2(this->l2);

    // Load L3 integral
    this->l3 = new P3();
    set_data_l3(this->l3);

    // Load L3_dA integral
    this->l3_da = new P3();
    set_data_l3_da(this->l3_da);

    // Load L3_dB integral
    this->l3_db = new P3();
    set_data_l3_db(this->l3_db);

    // Load M1 integral
    this->m1 = new P3();
    set_data_m1(this->m1);

    // Load M1_dA integral
    this->m1_da = new P3();
    set_data_m1_da(this->m1_da);

    // Load M1_dB integral
    this->m1_db = new P3();
    set_data_m1_db(this->m1_db);

    // Load M2 integral
    this->m2 = new P3();
    set_data_m2(this->m2);

    // Load M3 integral
    this->m3 = new P3();
    set_data_m3(this->m3);

    // Load M3_dA integral
    this->m3_da = new P3();
    set_data_m3_da(this->m3_da);

    // Load M3_dB integral
    this->m3_db = new P3();
    set_data_m3_db(this->m3_db);

    // Set build flag to true
    this->is_build = true;
}


IntegralsDb::IntegralsDb( IntegralsDb& )
{
    this->_error_copy_constructor( );
}


IntegralsDb::IntegralsDb( IntegralsDb* )
{
    this->_error_copy_constructor( );
}


IntegralsDb::~IntegralsDb(void)
{
    if (this->is_build)
    {
        delete this->l1;
        delete this->l1_da;
        delete this->l1_db;
        delete this->l2;
        delete this->l3;
        delete this->l3_da;
        delete this->l3_db;
        delete this->m1;
        delete this->m1_da;
        delete this->m1_db;
        delete this->m2;
        delete this->m3;
        delete this->m3_da;
        delete this->m3_db;
        delete this->r11;
        delete this->r11a_dx;
        delete this->r11b_dx;
        delete this->r12;
        delete this->r12_dx;
        delete this->r21;
        delete this->r21_dx;
        delete this->r22;
        delete this->r22_dx;
    }
}


void build_integrals_db(IntegralsDb &idb)
{
    ////////////////////////////////////////////////
    ////// Load infinite water depth integrals//////
    ////////////////////////////////////////////////

    // Load R11 integral
    idb.r11 = new R11();
    set_data_r11(idb.r11);

    // Load R11A_dX integral
    idb.r11a_dx = new R11A_dX();
    set_data_r11a_dx(idb.r11a_dx);

    // Load R11B_dX integral
    idb.r11b_dx = new R11B_dX();
    set_data_r11b_dx(idb.r11b_dx);

    // Load R12 integral
    idb.r12 = new R12();
    set_data_r12(idb.r12);

     // Load R12_dX integral
    idb.r12_dx = new R12_dX();
    set_data_r12_dx(idb.r12_dx);

    // Load R21 integral
    idb.r21 = new R21();
    set_data_r21(idb.r21);

    // Load R21_dX integral
    idb.r21_dx = new R21_dX();
    set_data_r21_dx(idb.r21_dx);

    // Load R22 integral
    idb.r22 = new R22();
    set_data_r22(idb.r22);

    // Load R22_dX integral
    idb.r22_dx = new R22_dX();
    set_data_r22_dx(idb.r22_dx);

    ////////////////////////////////////////////////
    //////// Load finite water depth integrals//////
    ////////////////////////////////////////////////

    // Load L1 integral
    idb.l1 = new P3();
    set_data_l1(idb.l1);

    // Load L1_dA integral
    idb.l1_da = new P3();
    set_data_l1_da(idb.l1_da);

    // Load L1_dB integral
    idb.l1_db = new P3();
    set_data_l1_db(idb.l1_db);

    // Load L2 integral
    idb.l2 = new P3();
    set_data_l2(idb.l2);

    // Load L3 integral
    idb.l3 = new P3();
    set_data_l3(idb.l3);

    // Load L3_dA integral
    idb.l3_da = new P3();
    set_data_l3_da(idb.l3_da);

    // Load L3_dB integral
    idb.l3_db = new P3();
    set_data_l3_db(idb.l3_db);

    // Load M1 integral
    idb.m1 = new P3();
    set_data_m1(idb.m1);

    // Load M1_dA integral
    idb.m1_da = new P3();
    set_data_m1_da(idb.m1_da);

    // Load M1_dB integral
    idb.m1_db = new P3();
    set_data_m1_db(idb.m1_db);

    // Load M2 integral
    idb.m2 = new P3();
    set_data_m2(idb.m2);

    // Load M3 integral
    idb.m3 = new P3();
    set_data_m3(idb.m3);

    // Load M3_dA integral
    idb.m3_da = new P3();
    set_data_m3_da(idb.m3_da);

    // Load M3_dB integral
    idb.m3_db = new P3();
    set_data_m3_db(idb.m3_db);

    // Set build flag to true
    idb.is_build = true;

}
