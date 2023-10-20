
#ifndef __freq_solver_tools_hpp
#define __freq_solver_tools_hpp

// Include local modules
#include "../../containers/mpi_config.hpp"
#include "../../interfaces/hmf_interface.hpp"
#include "../../interfaces/grfdn_interface.hpp"
#include "../../interfaces/gwfdn_interface.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/input.hpp"
#include "../../inout/output.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"


// Define alias
typedef ScalapackSolver<cuscomplex> SclCmpx;


void    calculate_diffraction_forces_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                HMFInterface*   hmf_interf,
                                                cusfloat        w,
                                                cuscomplex*     wave_diffrac
                                        );


void    calculate_freq_domain_coeffs(
                                                MpiConfig*      mpi_config,
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                    );


void    calculate_froude_krylov(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     froude_krylov
                                );

void    calculate_global_hydstiffness(
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                cusfloat*       hydstiffness
                                    );

void    calculate_global_structural_mass(
                                                Input*          input,
                                                cusfloat*       structural_mass_p0
                                        );

void    calculate_hydromechanic_coeffs_lin( 
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panels_pot,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad
                                            );

void    calculate_hydromechanic_coeffs_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                HMFInterface*   hmf_interf,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad
                                        );

void    calculate_influence_potential_steady(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     inf_pot_mat
                                            );

void    calculate_influence_potential_total(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     inf_pot_steady,
                                                cuscomplex*     inf_pot_total
                                        );

void    calculate_panel_potentials_lin(
                                                Input*          input,
                                                cuscomplex*     inf_pot_mat,
                                                int             rows_np,
                                                int             cols_np,
                                                int             start_col,
                                                cuscomplex*     sources,
                                                cuscomplex*     panel_pot
                                        );

void    calculate_panel_potentials_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     all_sources,
                                                cusfloat        ang_freq
                                        );

void    calculate_sources_intensity(
                                                Input*          input,
                                                SclCmpx*        scl,
                                                MeshGroup*      mesh,
                                                GWFDnInterface* gwf_interf,
                                                cusfloat        w,
                                                cuscomplex*     sysmat_steady,
                                                cuscomplex*     sysmat,
                                                cuscomplex*     sources_int
                                   );

void    calculate_sources_sysmat_steady(
                                                Input*          input,
                                                SclCmpx*        scl,
                                                MeshGroup*      mesh_gp,
                                                GRFDnInterface* grf_interf,
                                                cuscomplex*     sysmat
                                        );

void    calculate_raos(
                                                Input*          input,
                                                cusfloat*       structural_mass,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad,
                                                cusfloat*       hydstiffness,
                                                cuscomplex*     wave_diffrac,
                                                cuscomplex*     froude_krylov,
                                                cusfloat        ang_freq,
                                                cuscomplex*     rao
                        );

#endif