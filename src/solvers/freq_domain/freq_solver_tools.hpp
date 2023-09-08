
#ifndef __freq_solver_tools_hpp
#define __freq_solver_tools_hpp

// Include local modules
#include "../../containers/mpi_config.hpp"
#include "../../green_interfaces/gwf_interface.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/input.hpp"
#include "../../inout/output.hpp"
#include "../../math/scalapack_solver.hpp"


// Define alias
typedef ScalapackSolver<cuscomplex> SclCmpx;


void    calculate_freq_domain_coeffs(
                                        MpiConfig*      mpi_config,
                                        Input*          input,
                                        Hydrostatics*   hydrostatics,
                                        Output*         output
                                    );

void    calculate_sources_intensity(
                                        SclCmpx*        scl,
                                        Mesh*           mesh,
                                        GWFInterface*   green_interf,
                                        cusfloat        w,
                                        cusfloat*       cog,
                                        cuscomplex*     sysmat,
                                        cuscomplex*     sources_int
                                   );

#endif