
#ifndef __integrals_db
#define __integrals_db

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


void build_integrals_db(IntegralsDb &idb);

#endif