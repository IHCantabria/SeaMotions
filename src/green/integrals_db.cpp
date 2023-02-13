
// Include general usage libraries
#include <cassert>

// Include local modules
#include "integrals_db.hpp"
#include "pulsating_fin_depth_cheby.hpp"


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


void build_integrals_db(IntegralsDb &idb)
{
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
    idb.l2 = new P3();
    set_data_l2(idb.l2);

    // Load L3 integral
    idb.l3 = new P3();
    set_data_l3(idb.l3);

    // Load L3_dA integral
    idb.l3_da = new P3();
    set_data_l3(idb.l3_da);

    // Load L3_dB integral
    idb.l3_db = new P3();
    set_data_l3(idb.l3_db);

    // Load M1 integral
    idb.m1 = new P3();
    set_data_m1(idb.m1);

    // Load M1_dA integral
    idb.m1_da = new P3();
    set_data_m1(idb.m1_da);

    // Load M1_dB integral
    idb.m1_db = new P3();
    set_data_m1(idb.m1_db);

    // Load M2 integral
    idb.m2 = new P3();
    set_data_m2(idb.m2);

    // Load M3 integral
    idb.m3 = new P3();
    set_data_m3(idb.m3);

    // Load M3_dA integral
    idb.m3_da = new P3();
    set_data_m3(idb.m3_da);

    // Load M3_dB integral
    idb.m3_db = new P3();
    set_data_m3(idb.m3_db);

    // Set build flag to true
    idb.is_build = true;

}
