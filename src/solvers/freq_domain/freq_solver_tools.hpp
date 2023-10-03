
#ifndef __freq_solver_tools_hpp
#define __freq_solver_tools_hpp

// Include local modules
#include "../../containers/mpi_config.hpp"
#include "../../interfaces/hmf_interface.hpp"
#include "../../interfaces/gwfdn_interface.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/input.hpp"
#include "../../inout/output.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"


// Define alias
typedef ScalapackSolver<cuscomplex> SclCmpx;


void    calculate_diffraction_forces(
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

void    calculate_hydromechanic_coeffs(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            HMFInterface*   hmf_interf,
                                            cusfloat        ang_freq,
                                            cusfloat*       added_mass,
                                            cusfloat*       damping_rad
                                        );

void    calculate_panel_potentials(
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
                                            GWFDnInterface* green_interf,
                                            cusfloat        w,
                                            cuscomplex*     sysmat,
                                            cuscomplex*     sources_int
                                   );

#endif