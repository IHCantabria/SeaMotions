
#ifndef __fd_wave_elevation_hpp
#define __fd_wave_elevation_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"


void    calculate_relative_wave_elevation_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MLGCmpx*        pot_gp,
                                                cuscomplex*     potpanel_total,
                                                cusfloat        ang_freq,
                                                cuscomplex*     raos,
                                                cuscomplex*     rel_wave_elevation
                                            );


void    calculate_wave_elevation_lin(
                                                cuscomplex*     pot_total,
                                                int             pot_total_np,
                                                cusfloat        ang_freq,
                                                cusfloat        grav_acc,
                                                cuscomplex*     wave_elevation
                                    );

#endif