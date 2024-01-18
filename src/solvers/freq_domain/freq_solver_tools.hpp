
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


void    calculate_freq_domain_coeffs(
                                                MpiConfig*      mpi_config,
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
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


void    freq_domain_linear_solver(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                MeshGroup*      mesh_fs_qtf_gp,
                                                SclCmpx*        scl,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                );


void    freq_domain_nonlinear_solver(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                SclCmpx*        scl,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                    );

#endif