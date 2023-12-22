
#ifndef __panel_fields_hpp
#define __panel_fields_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"

// Function declaration
template<typename T, typename V>    inline  void    calculate_fields_lin(
                                                                                    Input*          input,
                                                                                    MpiConfig*      mpi_config,
                                                                                    MeshGroup*      mesh_gp,
                                                                                    T*              field_func,
                                                                                    V*              fk_fcn,
                                                                                    cusfloat        ang_freq,
                                                                                    cuscomplex*     intensities,
                                                                                    cuscomplex*     raos,
                                                                                    MLGCmpx*        field_gp,
                                                                                    cuscomplex*     field_fk_p0,
                                                                                    cuscomplex*     field_raddif_p0,
                                                                                    cuscomplex*     field_total
                                                                        );

                                    inline  void    calculate_fields_raddif_lin(
                                                                                    Input*          input,
                                                                                    cuscomplex*     intensities,
                                                                                    MLGCmpx*        pot_gp
                                                                                );


template<typename T>                inline  void    calculate_influence_field_mat(
                                                                                    Input*          input,
                                                                                    MeshGroup*      mesh_gp,
                                                                                    T*              field_funct,
                                                                                    MLGCmpx*        field_gp
                                                                                );


                                    inline  void    calculate_pertubation_field(
                                                                                    Input*          input,
                                                                                    cuscomplex*     raos,
                                                                                    cusfloat        ang_freq,
                                                                                    MLGCmpx*        pot_gp,
                                                                                    cuscomplex*     pot_raddif,
                                                                                    cuscomplex*     pot_total
                                                                                );


                                    inline  void    calculate_raddiation_field(
                                                                                    Input*          input,
                                                                                    cuscomplex*     raos,
                                                                                    cusfloat        ang_freq,
                                                                                    MLGCmpx*        pot_gp,
                                                                                    cuscomplex*     pot_raddif,
                                                                                    cuscomplex*     pot_total
                                                                                );


                                    inline  void    calculate_total_field(
                                                                                    Input*          input,
                                                                                    cusfloat        ang_freq,
                                                                                    MLGCmpx*        pot_gp,
                                                                                    cuscomplex*     raos,
                                                                                    cuscomplex*     pot_fk,
                                                                                    cuscomplex*     pot_raddif,
                                                                                    cuscomplex*     pot_total
                                                                            );

// Include module function definition
#include "panel_fields.txx"


#endif