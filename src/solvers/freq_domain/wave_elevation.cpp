
// Include local modules
#include "../../config.hpp"
#include "../../inout/input.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../mesh/mesh_group.hpp"
#include "wave_elevation.hpp"


void    calculate_relative_wave_elevation_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     sources,
                                                cuscomplex*     pot_steady,
                                                cuscomplex*     raos,
                                                cusfloat*       field_points,
                                                int             field_points_np
                                            )
{
    /****************************************************/
    /************ Calculate wave elevation **************/
    /****************************************************/

    // Define parallel chunks
    int start_col   = 0;
    int end_col     = 0;
    mpi_config->get_1d_bounds( 
                                    mesh_gp->panels_wl_tnp, 
                                    start_col, 
                                    end_col
                                );
    int ipm_cols_np = start_col - end_col;

    

}