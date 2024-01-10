
#ifndef __output_hpp
#define __output_hpp

// Include general usage libraries
#include <string>

// Include local modules
#include "../hydrostatics.hpp"
#include "input.hpp"


// Define output channels name
#define  _DN_ADDED_MASS         "added_mass"
#define  _DN_DAMPING_RAD        "damping_rad"
#define  _DN_DIFFRAC            "diffraction_force"
#define  _DN_DIFFRAC_MAG        "diffraction_force_mag"
#define  _DN_DIFFRAC_PHA        "diffraction_force_pha"
#define  _DN_FREQS              "frequencies"
#define  _DN_FK                 "froude_krylov_force"
#define  _DN_FK_MAG             "froude_krylov_force_mag"
#define  _DN_FK_PHA             "froude_krylov_force_pha"
#define  _DN_HEADS              "headings"
#define  _DN_HYDSTIFF           "hydstiffness"
#define  _DN_MDRIFT             "mean_drift"
#define  _DN_MDRIFT_WL          "mean_drift_wl"
#define  _DN_MDRIFT_BERN        "mean_drift_bern"
#define  _DN_MDRIFT_ACC         "mean_drift_acc"
#define  _DN_MDRIFT_MOM         "mean_drift_mom"
#define  _DN_MDRIFT_MAG         "mean_drift_mag"
#define  _DN_MDRIFT_PHA         "mean_drift_pha"
#define  _DN_MDRIFT_WL_MAG      "mean_drift_wl_mag"
#define  _DN_MDRIFT_WL_PHA      "mean_drift_wl_pha"
#define  _DN_MDRIFT_BERN_MAG    "mean_drift_bern_mag"
#define  _DN_MDRIFT_BERN_PHA    "mean_drift_bern_pha"
#define  _DN_MDRIFT_ACC_MAG     "mean_drift_acc_mag"
#define  _DN_MDRIFT_ACC_PHA     "mean_drift_acc_pha"
#define  _DN_MDRIFT_MOM_MAG     "mean_drift_mom_mag"
#define  _DN_MDRIFT_MOM_PHA     "mean_drift_mom_pha"
#define  _DN_QTF_DIFF_MAG       "qtf_diff_mag"
#define  _DN_QTF_DIFF_PHA       "qtf_diff_pha"
#define  _DN_QTF_SUM_MAG        "qtf_sum_mag"
#define  _DN_QTF_SUM_PHA        "qtf_sum_pha"
#define  _DN_RAO                "rao"
#define  _DN_RAO_MAG            "rao_mag"
#define  _DN_RAO_PHA            "rao_pha"
#define  _DN_STRUCT_MASS        "mass"
#define  _DN_WEX                "wave_exciting"
#define  _DN_WEX_MAG            "wave_exciting_mag"
#define  _DN_WEX_PHA            "wave_exciting_pha"

// Define datasets shape length
const hsize_t _DS_HF_NP = 1;
const hsize_t _DS_HM_NP = 5;
const hsize_t _DS_MH_NP = 3;
const hsize_t _DS_QF_NP = 6;
const hsize_t _DS_WX_NP = 4;


struct Output
{
private:
    // Define class attributes
    hsize_t         _ds_fh[_DS_HF_NP]   = { 0 };
    hsize_t         _ds_hm[_DS_HM_NP]   = { 0, 0, 0, 0, 0 };
    hsize_t         _ds_mh[_DS_MH_NP]   = { 0, 0, 0 };
    hsize_t         _ds_qf[_DS_QF_NP]   = { 0, 0, 0, 0, 0, 0 };
    hsize_t         _ds_wx[_DS_WX_NP]   = { 0, 0, 0, 0 };
    Input*          _input              = nullptr;
    std::string     _results_fipath     = "";
    
public:
    // Define class attributes

    // Define class constructors and destructor
    Output(
                Input* input
            );

    // Define class methods
    void    save_frequencies(
                                            cusfloat*   freqs
                            );
    
    void    save_headings(
                                            cusfloat*   heads
                            );

    void    save_hydstiffness(
                                            Hydrostatics** hydrostatics
                                );

    void    save_hydromechanics_format(
                                            int         freq_index,
                                            std::string channel_name,
                                            cusfloat*   added_mass
                                        );

    void    save_structural_mass( 
                                            void 
                                );

    void    qtf_format(
                                            std::string channel_name,
                                            cuscomplex* forces
                        );

    void    save_wave_exciting_format(
                                            int         freq_index,
                                            std::string channel_name,
                                            cuscomplex* forces
                                        );
};


#endif