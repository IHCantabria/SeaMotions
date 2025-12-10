
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
#include <string>

// Include local modules
#include "gz_point.hpp"
#include "../../initial_stability.hpp"
#include "../../mesh/rigid_body_mesh.hpp"
#include "stab_input.hpp"


/*****************************************************************/
/*********** Auxiliar declarations for StabOutput ****************/
/*****************************************************************/

struct ScalarField 
{
    const char*                                         name;
    std::function< cusfloat ( std::size_t ) >           getter;
};


struct VectorField 
{
    const char*                                         name;
    std::function<const cusfloat* ( std::size_t ) >     getter;
};


template<typename T>
inline void _extract_hydrostats_scalar( 
                                            std::size_t heel_np,
                                            std::size_t draft_np,
                                            T           extract_fcn,
                                            cusfloat*   amb    
                                        )
{
    std::size_t index   = 0;
    std::size_t count   = 0;
    for ( std::size_t i=0; i<heel_np; i++ )
    {
        for ( std::size_t j=0; j<draft_np; j++ )
        {
            index           = i * draft_np + j;
            amb[count]      = extract_fcn( index );
            count++;
        }
    }
}


template<typename T>
inline void _extract_hydrostats_vector( 
                                            std::size_t heel_np,
                                            std::size_t draft_np,
                                            T           extract_fcn,
                                            cusfloat*   amb    
                                        )
{
    std::size_t index   = 0;
    std::size_t count   = 0;
    for ( std::size_t i=0; i<heel_np; i++ )
    {
        for ( std::size_t j=0; j<draft_np; j++ )
        {
            for ( std::size_t k=0; k<3; k++ )
            {
                index           = i * draft_np + j;
                amb[count]      = extract_fcn( index )[k];
                count++;
            }
        }
    }
}


// Define output channels name
#define  _DN_AREA_WL            "area_wl"
#define  _DN_AREA_IXX_WL        "area_ixx_wl"
#define  _DN_AREA_IYY_WL        "area_iyy_wl"
#define  _DN_AREA_MX_WL         "area_mx_wl"
#define  _DN_AREA_MY_WL         "area_my_wl"
#define  _DN_BMX                "bmx"
#define  _DN_BMY                "bmy"
#define  _DN_COB                "cob"
#define  _DN_DISPLACEMENT       "displacement"
#define  _DN_DRAFT              "drafts"
#define  _DN_GMX                "gmx"
#define  _DN_GMY                "gmy"
#define  _DN_GZ                 "gz"
#define  _DN_HEEL               "heelings"
#define  _DN_KMX                "kmx"
#define  _DN_KMY                "kmy"
#define  _DN_VOLUME             "volume"
#define  _DN_VOLUME_MX          "volume_mx"
#define  _DN_VOLUME_MY          "volume_my"
#define  _DN_VOLUME_MZ          "volume_mz"
#define  _GN_EQ                 "eq"
#define  _GN_GZ                 "gz"
#define  _GN_HYDROSTATS         "hydrostatics"


// Define datasets shape length
const hsize_t _DS_1D_NP     = 1;    // Dataset shape length for 1D vectors as input drafts or heelings.
const hsize_t _DS_GZ_NP     = 2;    // GZ curves
const hsize_t _DS_HS_S_NP   = 3;    // Hydrostatics scalar fields
const hsize_t _DS_HS_V_NP   = 4;    // Hydrostatics vector fields


/*****************************************************************/
/*************** StabOutput template declarations ****************/
/*****************************************************************/
struct StabOutput
{
private:
    /* Defince class local type aliases  */
    using   GZPVec      = std::vector< GZPoint >;
    using   InitStabVec = std::vector< InitialStability<NUM_GP, RigidBodyMesh> >;

    /* Define class private attributes */
    hsize_t         _ds_gz_h[_DS_1D_NP]     = { 0 };            // Dataset shape vector for GZ heeling points
    hsize_t         _ds_gz[_DS_GZ_NP]       = { 0, 0 };         // Dataset shape vector for GZ curves: Axis | Heelings
    hsize_t         _ds_hs_d[_DS_1D_NP]     = { 0 };            // Dataset shape vector for hydrostatics draft points
    hsize_t         _ds_hs_h[_DS_1D_NP]     = { 0 };            // Dataset shape vector for hydrostatics heeling points
    hsize_t         _ds_hs_s[_DS_HS_S_NP]   = { 0, 0, 0 };      // Dataset shape vector for Hydrostatic scalar values: Axis | Drafts | Heelings
    hsize_t         _ds_hs_v[_DS_HS_V_NP]   = { 0, 0, 0, 0 };   // Dataset shape vector for Hydrostatic vector values: Axis | Drafts | Heelings | VectorDims
    StabInput*      _input                  = nullptr;          // Pointer to stabilit input system to have access to all the case configurations
    std::string     _results_fipath         = "";               // Output file path

    /* Define class private methods */

    /**
     * @brief   Create GZ field in the output file
     */
    void    _create_gz_field(
                                        H5::Group&          gz_gp,
                                        const char*         field_name
                            );

    /**
     * @brief   Create hydrostatics scalar field in the output file
     */
    void    _create_hs_scalar_field(
                                        H5::Group&  gp,
                                        const char* ds_name
                                    );

    /**
     * @brief   Create hydrostatics vector field in the output file
     */
    void    _create_hs_vector_field(
                                        H5::Group&  gp,
                                        const char* ds_name
                                    );

    /**
     * @brief   Save input 1D vector data
     */
    void    _save_1D_input_data(
                                        H5::Group&              hs_gp,
                                        const char*             field_name,
                                        hsize_t*                dset_dims,
                                        std::vector<cusfloat>&  vec
                                );

    /**
     * @brief   Save hydrostatics scalar field into the output file
     */
    template<typename Getter>
    void    _save_hs_scalar_field(
                                        H5::Group&  gp,
                                        int         axis_id,
                                        const char* ds_name,
                                        cusfloat*   buffer,
                                        Getter&&    getter
                                );

    /**
     * @brief   Save hydrostatics scalar field into the output file
     */
    template<typename Getter>
    void    _save_hs_vector_field(
                                        H5::Group&  gp,
                                        int         axis_id,
                                        const char* ds_name,
                                        cusfloat*   buffer,
                                        Getter&&    getter
                                );
    
public:
    // Define class constructors and destructor
    StabOutput(
                    StabInput* input
                );

    // Define class methods

    /**
    * @brief   Save hydrostatics results into the output file
    * 
    * @param    axis_id     Axis identifier ( 0: Roll, 1: Pitch )
    * @param    hydrostats  Vector of InitialStability objects containing hydrostatics results      
    * 
    **/
    void    save_hydrostatics(
                                        int             axis_id,
                                        InitStabVec&    hydrostats
                                );

    /**
     * @brief   Save GZ curves into the output file
     * 
     * @param    lc_name     Loading condition name
     * @param    gz_points   Vector of GZPoint objects containing GZ curve points
     * 
     */
    void    save_gz(
                                        std::string&    lc_name,
                                        GZPVec&         gz_points
                                );
    
};

// Include template definitions
#include "stab_output.txx"