
#ifndef __integrals_db
#define __integrals_db

// Include local modules
#include "pulsating_fin_depth_cheby.hpp"
#include "pulsating_inf_depth_cheby.hpp"


class IntegralsDb
{
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
    IntegralsDb( ) = default;

    IntegralsDb( void );

    ~IntegralsDb(void);

    // Declare class methods
    void fold_h(cusfloat H);
    
};


void build_integrals_db(IntegralsDb &idb);

#endif