
#ifndef __fd_wave_elevation_hpp
#define __fd_wave_elevation_hpp

void    calculate_relative_wave_elevation_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                cusfloat*       cog_to_fp,
                                                int*            fp_cnp,
                                                int             fp_nb,
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