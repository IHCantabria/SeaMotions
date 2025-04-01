
#ifndef __input_hpp
#define __input_hpp

// Include general usage libraries
#include <iostream>
#include <string>

// Import local modules
#include "../config.hpp"
#include "../containers/body_def.hpp"
#include "../mesh/mesh.hpp"


struct Input
{
public:
    // Define class attributes
    std::vector<cusfloat>       angfreqs            ;
    int                         angfreqs_np         = 0;
    BodyDef**                   bodies              = nullptr;
    std::vector<std::string>    bodies_finame       ;
    int                         bodies_np           = 0;
    std::string                 case_fopath         = "";
    const int                   dofs_np             = 6;
    std::string                 folder_path         = "";
    int                         gauss_order         = 0;
    cusfloat                    gfdn_abs_err        = 0.0;
    cusfloat                    gfdn_rel_err        = 0.0;
    cusfloat                    grav_acc            = 0.0;
    cusfloat*                   freqs               = nullptr;
    std::string                 freqs_unit          = "";
    int                         kochin_np           = 0;
    std::vector<cusfloat>       heads               ;
    int                         heads_np            = 0;
    std::string                 heads_units         = "";
    bool                        is_calc_mdrift      = false;
    bool                        is_block_adaption   = false;
    bool                        is_bodies           = false;
    bool                        is_fast_solver      = false;
    bool                        is_fs_qtf           = false;
    bool                        is_log_sin_ana      = false;
    bool                        is_wl_points        = false;
    bool                        out_diffrac         = false;
    bool                        out_fk              = false;
    bool                        out_hydmech         = false;
    bool                        out_hydstiff        = false;
    bool                        out_potential       = false;
    bool                        out_pressure        = false;
    bool                        out_mdrift          = false;
    bool                        out_mesh            = false;
    bool                        out_qtf             = false;
    bool                        out_qtf_comp        = false;
    int                         out_qtf_so_model    = 0;
    bool                        out_raos            = false;
    bool                        out_sources         = false;
    bool                        out_struct_mass     = false;
    bool                        out_wex             = false;
    int                         poly_order          = 0;
    cusfloat                    pot_abs_err         = 0.0;
    cusfloat                    pot_rel_err         = 0.0;
    cusfloat                    press_abs_err       = 0.0;
    cusfloat                    press_rel_err       = 0.0;
    cusfloat                    water_density       = 0.0;
    cusfloat                    water_depth         = 0.0;
    cusfloat                    wave_amplitude      = 1.0;
    cusfloat                    wl_det_prec         = 0.0;

    // Define class constructors and destructors
    ~Input( void );

    // Define class methods
    void    configure( 
                                    void
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

};

#endif