
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

// Include general usage libraries
#include <string>

// Include local modules
#include "../hydrostatics.hpp"
#include "input.hpp"


// Define output channels name
#define  _DN_ADDED_MASS         "added_mass"
#define  _DN_ADDED_MASS_HF      "added_mass_hf"
#define  _DN_ADDED_MASS_LF      "added_mass_lf"
#define  _DN_DAMPING_RAD        "damping_rad"
#define  _DN_DAMPING_RAD_HF     "damping_rad_hf"
#define  _DN_DAMPING_RAD_LF     "damping_rad_lf"
#define  _DN_DIFFRAC            "diffraction_force"
#define  _DN_DIFFRAC_HF         "diffraction_force_hf"
#define  _DN_DIFFRAC_LF         "diffraction_force_lf"
#define  _DN_DIFFRAC_MAG        "diffraction_force_mag"
#define  _DN_DIFFRAC_HF_MAG     "diffraction_force_hf_mag"
#define  _DN_DIFFRAC_LF_MAG     "diffraction_force_lf_mag"
#define  _DN_DIFFRAC_PHA        "diffraction_force_pha"
#define  _DN_DIFFRAC_HF_PHA     "diffraction_force_hf_pha"
#define  _DN_DIFFRAC_LF_PHA     "diffraction_force_lf_pha"
#define  _DN_ELEMS              "elements"
#define  _DN_FREQS              "frequencies"
#define  _DN_FK                 "froude_krylov_force"
#define  _DN_FK_HF              "froude_krylov_force_hf"
#define  _DN_FK_LF              "froude_krylov_force_lf"
#define  _DN_FK_MAG             "froude_krylov_force_mag"
#define  _DN_FK_HF_MAG          "froude_krylov_force_hf_mag"
#define  _DN_FK_LF_MAG          "froude_krylov_force_lf_mag"
#define  _DN_FK_PHA             "froude_krylov_force_pha"
#define  _DN_FK_HF_PHA          "froude_krylov_force_hf_pha"
#define  _DN_FK_LF_PHA          "froude_krylov_force_lf_pha"
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
#define  _DN_NODES_X            "nodes_x"
#define  _DN_NODES_Y            "nodes_y"
#define  _DN_NODES_Z            "nodes_z"
#define  _DN_POT_INT            "potential"
#define  _DN_POT_INT_MAG        "potential_mag"
#define  _DN_POT_INT_PHA        "potential_pha"
#define  _DN_PRESS_INT          "pressure"
#define  _DN_PRESS_INT_MAG      "pressure_mag"
#define  _DN_PRESS_INT_PHA      "pressure_pha"
#define  _DN_QTF_DIFF_MAG       "qtf_diff_mag"
#define  _DN_QTF_DIFF_PHA       "qtf_diff_pha"
#define  _DN_QTF_DIFF_ACC_MAG   "qtf_diff_acc_mag"
#define  _DN_QTF_DIFF_ACC_PHA   "qtf_diff_acc_pha"
#define  _DN_QTF_DIFF_BERN_MAG  "qtf_diff_bern_mag"
#define  _DN_QTF_DIFF_BERN_PHA  "qtf_diff_bern_pha"
#define  _DN_QTF_DIFF_MOM_MAG   "qtf_diff_mom_mag"
#define  _DN_QTF_DIFF_MOM_PHA   "qtf_diff_mom_pha"
#define  _DN_QTF_DIFF_SOP_MAG   "qtf_diff_sop_mag"
#define  _DN_QTF_DIFF_SOP_PHA   "qtf_diff_sop_pha"
#define  _DN_QTF_DIFF_WL_MAG    "qtf_diff_wl_mag"
#define  _DN_QTF_DIFF_WL_PHA    "qtf_diff_wl_pha"
#define  _DN_QTF_SUM_MAG        "qtf_sum_mag"
#define  _DN_QTF_SUM_PHA        "qtf_sum_pha"
#define  _DN_QTF_SUM_ACC_MAG    "qtf_sum_acc_mag"
#define  _DN_QTF_SUM_ACC_PHA    "qtf_sum_acc_pha"
#define  _DN_QTF_SUM_BERN_MAG   "qtf_sum_bern_mag"
#define  _DN_QTF_SUM_BERN_PHA   "qtf_sum_bern_pha"
#define  _DN_QTF_SUM_MOM_MAG    "qtf_sum_mom_mag"
#define  _DN_QTF_SUM_MOM_PHA    "qtf_sum_mom_pha"
#define  _DN_QTF_SUM_SOP_MAG    "qtf_sum_sop_mag"
#define  _DN_QTF_SUM_SOP_PHA    "qtf_sum_sop_pha"
#define  _DN_QTF_SUM_WL_MAG     "qtf_sum_wl_mag"
#define  _DN_QTF_SUM_WL_PHA     "qtf_sum_wl_pha"
#define  _DN_RAO                "rao"
#define  _DN_RAO_MAG            "rao_mag"
#define  _DN_RAO_PHA            "rao_pha"
#define  _DN_SRC_INT            "source_intensity"
#define  _DN_SRC_INT_MAG        "source_intensity_mag"
#define  _DN_SRC_INT_PHA        "source_intensity_pha"
#define  _DN_STRUCT_MASS        "mass"
#define  _DN_WEX                "wave_exciting"
#define  _DN_WEX_HF             "wave_exciting_hf"
#define  _DN_WEX_LF             "wave_exciting_lf"
#define  _DN_WEX_MAG            "wave_exciting_mag"
#define  _DN_WEX_HF_MAG         "wave_exciting_hf_mag"
#define  _DN_WEX_LF_MAG         "wave_exciting_lf_mag"
#define  _DN_WEX_PHA            "wave_exciting_pha"
#define  _DN_WEX_HF_PHA         "wave_exciting_hf_pha"
#define  _DN_WEX_LF_PHA         "wave_exciting_lf_pha"
#define  _GN_MESH               "/mesh"

// Define datasets shape length
const hsize_t _DS_EL_NP     = 2;
const hsize_t _DS_FD_NP     = 3;
const hsize_t _DS_HF_NP     = 1;
const hsize_t _DS_HM_NP     = 5;
const hsize_t _DS_HMA_NP    = 4;
const hsize_t _DS_MH_NP     = 3;
const hsize_t _DS_ND_NP     = 1;
const hsize_t _DS_QF_NP     = 6;
const hsize_t _DS_WX_NP     = 4;
const hsize_t _DS_WXA_NP    = 3;


struct Output
{
private:
    // Define class attributes
    hsize_t         _ds_f[_DS_HF_NP]    = { 0 };                    // Dataset shape for frequencies
    hsize_t         _ds_fd[_DS_FD_NP]   = { 0, 0, 0 };              // Dataset shape for field data (sources, potential, pressure)
    hsize_t         _ds_h[_DS_HF_NP]    = { 0 };                    // Dataset shape for headings
    hsize_t         _ds_hm[_DS_HM_NP]   = { 0, 0, 0, 0, 0 };        // Dataset shape for hydromechanics values
    hsize_t         _ds_hma[_DS_HMA_NP] = { 0, 0, 0, 0 };           // Dataset shape for hydromechanic asymptotic values
    hsize_t         _ds_mh[_DS_MH_NP]   = { 0, 0, 0 };              // Dataset shape for mass and hydrostatic stiffness matrixes
    hsize_t         _ds_qf[_DS_QF_NP]   = { 0, 0, 0, 0, 0, 0 };     // Dataset shape for QTF values
    hsize_t         _ds_wx[_DS_WX_NP]   = { 0, 0, 0, 0 };           // Dataset shape for wave exciting forces ( diffraction, Froude-Krylov, total )
    hsize_t         _ds_wxa[_DS_WXA_NP] = { 0, 0, 0 };              // Dataset shape for wave exciting forces asymptotic values
    Input*          _input              = nullptr;
    std::string     _results_fipath     = "";
    int             _total_panels_np    = 0;
    
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

    void    save_hydromechanics_asympt_format(
                                                std::string channel_name,
                                                cusfloat*   added_mass
                                            );

    void    save_mesh(
                                                void
                        );

    void    save_fields_data(
                                                int         freq_index,
                                                std::string channel_name,
                                                cuscomplex* source_intensity
                                        );

    void    save_structural_mass( 
                                                void 
                                );

    void    save_qtf_format(
                                                std::string channel_name,
                                                cuscomplex* forces
                        );

    void    save_wave_exciting_format(
                                                int         freq_index,
                                                std::string channel_name,
                                                cuscomplex* forces
                                        );

    void    save_wave_exciting_asympt_format(
                                                std::string channel_name,
                                                cuscomplex* forces
                                        );

};
