
#ifndef __integrals_db
#define __integrals_db

// Include local modules
#include "pulsating_fin_depth_cheby.hpp"


class IntegralsDb
{
public:
    bool is_build = false;
    P3* l1;
    P3* l1_da;
    P3* l1_db;
    P3* l2;
    P3* l3;
    P3* l3_da;
    P3* l3_db;
    P3* m1;
    P3* m1_da;
    P3* m1_db;
    P3* m2;
    P3* m3;
    P3* m3_da;
    P3* m3_db;

    ~IntegralsDb(void)
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
        }
    }
};


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


#endif