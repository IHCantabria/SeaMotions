
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
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Import local modules
#include "../config.hpp"
#include "../containers/body_def.hpp"
#include "../containers/field_points_def.hpp"
#include "../mesh/mesh.hpp"

struct CliOptions;


struct Input
{
public:
    // Define class attributes
    std::vector<cusfloat>           angfreqs            ;
    int                             angfreqs_np         = 0;
    std::vector<cusfloat>           wave_numbers        ;
    BodyDef**                       bodies              = nullptr;
    std::vector<std::string>        bodies_finame       ;
    int                             bodies_np           = 0;
    std::string                     case_fopath         = "";
    const int                       dofs_np             = 6;
    std::vector<FieldPointsDef>     field_points        ;
    std::vector<std::string>        field_points_finame ;
    std::size_t                     field_points_np     = 0;
    std::string                     folder_path         = "";
    int                             gauss_order         = 2;
    cusfloat                        gfdn_abs_err        = 1e-3;
    cusfloat                        gfdn_rel_err        = 1e-3;
    cusfloat                        grav_acc            = 0.0;
    cusfloat*                       freqs               = nullptr;
    std::string                     freqs_unit          = "";
    int                             kochin_np           = 30;
    std::vector<cusfloat>           heads               ;
    int                             heads_np            = 0;
    std::string                     heads_units         = "";
    bool                            is_calc_mdrift      = false;
    bool                            is_block_adaption   = false;
    bool                            is_bodies           = false;
    bool                            is_fast_solver      = false;
    bool                            is_fs_qtf           = false;
    bool                            is_log_sin_ana      = false;
    bool                            is_wl_points        = false;
    bool                            out_diffrac         = false;
    bool                            out_fk              = false;
    bool                            out_hydmech         = false;
    bool                            out_hydstiff        = false;
    bool                            out_potential       = false;
    bool                            out_pressure        = false;
    bool                            out_mdrift          = false;
    bool                            out_mesh            = false;
    bool                            out_morison         = false;
    bool                            out_qtf             = false;
    bool                            out_qtf_comp        = false;
    QTFSOModelE                     out_qtf_so_model    = QTFSOModelE::PINKSTER;
    std::string                     qtf_pairs_mode                              = "FULL";
    std::vector<std::pair<cusfloat,cusfloat>> qtf_pair_values_in                ;
    std::vector<std::pair<int,int>> qtf_freq_pairs                              ;
    bool                            qtf_use_all_pairs                           = true;
    bool                            out_raos            = false;
    bool                            out_sources         = false;
    bool                            out_struct_mass     = false;
    bool                            out_wex             = false;
    bool                            out_memory_report   = false;
    int                             poly_order          = 0;
    cusfloat                        pot_abs_err         = 1e-3;
    cusfloat                        pot_rel_err         = 1e-3;
    cusfloat                        press_abs_err       = 1e-3;
    cusfloat                        press_rel_err       = 1e-3;
    bool                            use_field_points    = false;
    cusfloat                        water_density       = 0.0;
    cusfloat                        water_depth         = 0.0;
    cusfloat                        wave_amplitude      = 1.0;
    cusfloat                        wl_det_prec         = 0.0;
    int                             qtf_annuli_np       = 6;
    int                             qtf_annuli_theta_np = 64;

    // Define class constructors and destructors
    Input( ) = default;
    explicit Input( const std::string& folder_path );
    ~Input( void );

    // Define class methods
    void    configure( 
                                    void
                    );

    void    load(
                                    const std::string& folder_path
                );

    void    apply_cli_options(
                                    const CliOptions& options
                );

    int     gauss_np_factor_1d( 
                                    void
                                );

    int     gauss_np_factor_2d( 
                                    void
                                );

    void    print(     
                                    void 
                );

    void    get_wave_quality_params(
                                    cusfloat& min_wavelength,
                                    cusfloat& max_wave_number
                ) const;

    /**
     * @brief Return whether the (i,j) frequency pair must be evaluated
     *        during the second-order (QTF) calculation.
     *
     * When the user requested `QTFPairsMode = FULL` all pairs are
     * computed. Under `CUSTOM` mode, only the pairs explicitly listed
     * by the user are evaluated. Uncomputed cells of the QTF output
     * matrix remain zero.
     *
     * @param i  Row index into the angular-frequency array.
     * @param j  Column index into the angular-frequency array.
     * @return   True if the pair must be computed.
     */
    bool    should_compute_qtf_pair(
                                    int i,
                                    int j
                ) const;

private:
    void    read_case(
                                    const std::string& folder_path
                );

    void    read_bodies(
                                    const std::string& folder_path
                );

    void    read_body(
                                    const std::string& folder_path,
                                    const std::string& target_file,
                                    BodyDef* body
                );

    void    read_field_points( 
                                    const std::string& folder_path 
                                );

};
