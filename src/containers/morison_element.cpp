
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

// Include general usage libraries
#include <numeric>


// Include local modules
#include "../math/euler_transforms.hpp"
#include "morison_element.hpp"
#include "../raos.hpp"


void    MorisonElement::_calculate_field_points( 
                                                    MorisonElementDef*  morison_def_
                                                )
{
    // Calculate rotation matrix for the Morison element based on the input angles
    cusfloat rot_mat[9];

    euler_rpy( 
                morison_def_->rotation_angle, 
                morison_def_->declination_angle, 
                morison_def_->azimuth_angle, 
                rot_mat 
            );

    // Rotate unit vector (1, 0, 0) to obtain the direction of the element axis
    this->_x_axis_l[0] = rot_mat[0];
    this->_x_axis_l[1] = rot_mat[3];
    this->_x_axis_l[2] = rot_mat[6];

    this->_y_axis_l[0] = rot_mat[1];
    this->_y_axis_l[1] = rot_mat[4];
    this->_y_axis_l[2] = rot_mat[7];

    this->_z_axis_l[0] = rot_mat[2];
    this->_z_axis_l[1] = rot_mat[5];
    this->_z_axis_l[2] = rot_mat[8];

    // Calculate divison segement length
    this->_div_len     = morison_def_->length / static_cast<cusfloat>( morison_def_->divisions_np );

    // Calculate nodes positions along the element axis
    cut::CusTensor<cusfloat> nodes_pos( { morison_def_->divisions_np + 1, 3 } );
    for ( std::size_t i = 0; i < morison_def_->divisions_np + 1; i++ )
    {
        nodes_pos( i, 0 ) = morison_def_->start_pos[0] + this->_x_axis_l[0] * this->_div_len * static_cast<cusfloat>( i );
        nodes_pos( i, 1 ) = morison_def_->start_pos[1] + this->_x_axis_l[1] * this->_div_len * static_cast<cusfloat>( i );
        nodes_pos( i, 2 ) = morison_def_->start_pos[2] + this->_x_axis_l[2] * this->_div_len * static_cast<cusfloat>( i );
    }

    // Calculate field points positions by replicating the nodes positions for each field point
    this->_field_points_l = cut::CusTensor<cusfloat>( { morison_def_->divisions_np, 3 } );
    this->_field_points_g = cut::CusTensor<cusfloat>( { morison_def_->divisions_np, 3 } );
    for ( std::size_t i = 0; i < morison_def_->divisions_np; i++ )
    {
        this->_field_points_l( i, 0 ) = ( nodes_pos( i, 0 ) + nodes_pos( i+1, 0 ) ) / 2.0;
        this->_field_points_l( i, 1 ) = ( nodes_pos( i, 1 ) + nodes_pos( i+1, 1 ) ) / 2.0;
        this->_field_points_l( i, 2 ) = ( nodes_pos( i, 2 ) + nodes_pos( i+1, 2 ) ) / 2.0;

        this->_field_points_g( i, 0 ) = this->_field_points_l( i, 0 ) + morison_def_->cog[0];
        this->_field_points_g( i, 1 ) = this->_field_points_l( i, 1 ) + morison_def_->cog[1];
        this->_field_points_g( i, 2 ) = this->_field_points_l( i, 2 ) + morison_def_->cog[2];
    }

}


void    MorisonElement::calculate_hydrodynamic_forces( 
                                                        cusfloat                    ang_freq,
                                                        cuscomplex*                 raos,
                                                        cut::CusTensor<cuscomplex>& inertial_force,
                                                        cut::CusTensor<cuscomplex>& drag_force
                                                    ) const
{
    // Define local variables for the calculation of hydrodynamic forces
    std::array<cuscomplex, 3>   drag_force_l        = { 0.0, 0.0, 0.0 };            // Drag force at the field point in local coordinates
    std::array<cuscomplex, 3>   fluid_acc_g         = { 0.0, 0.0, 0.0 };            // Fluid acceleration at the field point in global coordinates
    std::array<cuscomplex, 3>   fluid_acc_l         = { 0.0, 0.0, 0.0 };            // Fluid acceleration at the field point in local coordinates
    std::array<cuscomplex, 3>   inertial_force_l    = { 0.0, 0.0, 0.0 };            // Inertial force at the field point in local coordinates
    cusfloat                    wave_period         = 2.0 * PI / ang_freq;          // Wave period
    std::array<cuscomplex, 3>   rao_disp_fp         = { 0.0, 0.0, 0.0 };            // RAO displacement at the field point
    std::array<cuscomplex, 3>   rao_vel_fp          = { 0.0, 0.0, 0.0 };            // RAO velocity at the field point
    std::array<cuscomplex, 3>   rao_acc_fp          = { 0.0, 0.0, 0.0 };            // RAO acceleration at the field point
    std::array<cuscomplex, 3>   rel_acc_g           = { 0.0, 0.0, 0.0 };            // Relative acceleration at the field point in global coordinates
    std::array<cuscomplex, 3>   rel_acc_l           = { 0.0, 0.0, 0.0 };            // Relative acceleration at the field point in local coordinates
    std::array<cuscomplex, 3>   rel_vel_g           = { 0.0, 0.0, 0.0 };            // Relative velocity at the field point in global coordinates
    std::array<cuscomplex, 3>   rel_vel_l           = { 0.0, 0.0, 0.0 };            // Relative velocity at the field point in local coordinates
    cusfloat                    rhow                = this->_input->water_density;  // Water density (local copy for faster access and clear notation)

    for ( std::size_t ih = 0; ih < this->_input->heads_np; ih++ )
    {
        // Get field point velocities
        const cut::CusTensor<cuscomplex>* vel_x = this->_rdd_morison->template get_field_data<FieldTypeE::VELOCITY_X, FieldComponentE::TOTAL, ComplexDataTypeE::COMPLEX>( ih );
        const cut::CusTensor<cuscomplex>* vel_y = this->_rdd_morison->template get_field_data<FieldTypeE::VELOCITY_Y, FieldComponentE::TOTAL, ComplexDataTypeE::COMPLEX>( ih );
        const cut::CusTensor<cuscomplex>* vel_z = this->_rdd_morison->template get_field_data<FieldTypeE::VELOCITY_Z, FieldComponentE::TOTAL, ComplexDataTypeE::COMPLEX>( ih );

        for ( std::size_t i = 0; i < this->_rdd_morison->get_size_local( ); i++ )
        {
            // Calculate RAO displacement at the field point
            calculate_rao_disp( 
                                    raos,
                                &(this->_field_points_l.data( )[3*i]), 
                                rao_disp_fp.data( ) 
                            );

            // Calculate velocity and acceleration raos at the field point
            for ( std::size_t j = 0; j < 3; j++ )
            {
                rao_vel_fp[j] = cuscomplex( 0.0, - ang_freq ) * rao_disp_fp[j];
                rao_acc_fp[j] = cuscomplex( 0.0, - ang_freq ) * rao_vel_fp[j];
            }

            // Calculate relative velocity and acceleration at the field point
            fluid_acc_g[0] = cuscomplex( 0.0, - ang_freq ) * vel_x->data( )[i];
            fluid_acc_g[1] = cuscomplex( 0.0, - ang_freq ) * vel_y->data( )[i];
            fluid_acc_g[2] = cuscomplex( 0.0, - ang_freq ) * vel_z->data( )[i];

            rel_vel_g[0]   = vel_x->data( )[i] - rao_vel_fp[0];
            rel_vel_g[1]   = vel_y->data( )[i] - rao_vel_fp[1];
            rel_vel_g[2]   = vel_z->data( )[i] - rao_vel_fp[2];

            rel_acc_g[0]   = fluid_acc_g[0] - rao_acc_fp[0];
            rel_acc_g[1]   = fluid_acc_g[1] - rao_acc_fp[1];
            rel_acc_g[2]   = fluid_acc_g[2] - rao_acc_fp[2];

            // Project relative velocity and acceleration onto the local axes of the Morison element
            fluid_acc_l[0] = std::inner_product( this->_x_axis_l.begin( ), this->_x_axis_l.end( ), fluid_acc_g.begin( ), cuscomplex( 0.0 ) );
            fluid_acc_l[1] = std::inner_product( this->_y_axis_l.begin( ), this->_y_axis_l.end( ), fluid_acc_g.begin( ), cuscomplex( 0.0 ) );
            fluid_acc_l[2] = std::inner_product( this->_z_axis_l.begin( ), this->_z_axis_l.end( ), fluid_acc_g.begin( ), cuscomplex( 0.0 ) );

            rel_vel_l[0]   = std::inner_product( this->_x_axis_l.begin( ), this->_x_axis_l.end( ), rel_vel_g.begin( ), cuscomplex( 0.0 ) );
            rel_vel_l[1]   = std::inner_product( this->_y_axis_l.begin( ), this->_y_axis_l.end( ), rel_vel_g.begin( ), cuscomplex( 0.0 ) );
            rel_vel_l[2]   = std::inner_product( this->_z_axis_l.begin( ), this->_z_axis_l.end( ), rel_vel_g.begin( ), cuscomplex( 0.0 ) );

            rel_acc_l[0]   = std::inner_product( this->_x_axis_l.begin( ), this->_x_axis_l.end( ), rel_acc_g.begin( ), cuscomplex( 0.0 ) );
            rel_acc_l[1]   = std::inner_product( this->_y_axis_l.begin( ), this->_y_axis_l.end( ), rel_acc_g.begin( ), cuscomplex( 0.0 ) );
            rel_acc_l[2]   = std::inner_product( this->_z_axis_l.begin( ), this->_z_axis_l.end( ), rel_acc_g.begin( ), cuscomplex( 0.0 ) );

            // Calculate Keulegan-Carpenter number at the field point
            cusfloat velm  = std::sqrt( std::norm( rel_vel_l[0] ) + std::norm( rel_vel_l[1] ) + std::norm( rel_vel_l[2] ) );
            cusfloat kc    = velm * wave_period / this->_kc_length;

            // Calculate drag force at the field point in local coordinates
            drag_force_l[0] = 0.5 * rhow * this->_cd_axial->eval( kc ) * this->_area_axial_cd * this->_div_len * rel_vel_l[0] * std::abs( rel_vel_l[0] );
            drag_force_l[1] = 0.5 * rhow * this->_cd_trans->eval( kc ) * this->_area_trans_cd * this->_div_len * rel_vel_l[1] * std::abs( rel_vel_l[1] );
            drag_force_l[2] = 0.5 * rhow * this->_cd_vert->eval( kc )  * this->_area_vert_cd  * this->_div_len * rel_vel_l[2] * std::abs( rel_vel_l[2] );

            // Calculate inertial force at the field point in local coordinates
            cusfloat ca_axial = this->_ca_axial->eval( kc );
            cusfloat ca_trans = this->_ca_trans->eval( kc );
            cusfloat ca_vert  = this->_ca_vert->eval( kc );

            inertial_force_l[0] = rhow * this->_area_axial_ca * this->_div_len * ( ( 1.0 + ca_axial ) * fluid_acc_l[0] - ca_axial * rel_acc_l[0] );
            inertial_force_l[1] = rhow * this->_area_trans_cd * this->_div_len * ( ( 1.0 + ca_trans ) * fluid_acc_l[1] - ca_trans * rel_acc_l[1] );
            inertial_force_l[2] = rhow * this->_area_vert_cd  * this->_div_len * ( ( 1.0 + ca_vert )  * fluid_acc_l[2] - ca_vert  * rel_acc_l[2] );

            // Calculate drag and inertial forces at the field point in global coordinates by projecting the local forces onto the global axes of the Morison element
            for ( std::size_t j = 0; j < 3; j++ )
            {
                drag_force( j, i )     += drag_force_l[0] * this->_x_axis_l[j] + drag_force_l[1] * this->_y_axis_l[j] + drag_force_l[2] * this->_z_axis_l[j];
                inertial_force( j, i ) += inertial_force_l[0] * this->_x_axis_l[j] + inertial_force_l[1] * this->_y_axis_l[j] + inertial_force_l[2] * this->_z_axis_l[j];
            }

            cross( this->_field_points_l.data( ) + 3*i,     drag_force.data( ) + 6*ih,     drag_force.data( ) + 6*ih + 3 );
            cross( this->_field_points_l.data( ) + 3*i, inertial_force.data( ) + 6*ih, inertial_force.data( ) + 6*ih + 3 );

        }
    }
}


void    MorisonElement::_initialize_element( 
                                                MpiConfig*          mpi_config_,
                                                Input*              input_,
                                                MorisonElementDef*  morison_def_
                                            )
{
    // Store Morison element definition
    this->_area_axial_ca    = morison_def_->area_axial_ca;
    this->_area_axial_cd    = morison_def_->area_axial_cd;
    this->_area_trans_cd    = morison_def_->area_transversal_cd;
    this->_area_vert_cd     = morison_def_->area_vertical_cd;
    this->_kc_length        = morison_def_->kc_length;

    // Initialize interpolation objects for added mass coefficients
    this->_ca_axial     = new MorisonCoeffCurve(
                                                    input_->case_fopath,
                                                    morison_def_->axis_added_mass_coeff_file
                                                );

    this->_ca_trans     = new MorisonCoeffCurve(
                                                    input_->case_fopath,
                                                    morison_def_->transversal_added_mass_coeff_file
                                                );

    this->_ca_vert      = new MorisonCoeffCurve(
                                                    input_->case_fopath,
                                                    morison_def_->vertical_added_mass_coeff_file
                                                );

    // Initialize interpolation objects for drag coefficients
    this->_cd_axial     = new MorisonCoeffCurve(
                                                    input_->case_fopath,
                                                    morison_def_->axis_drag_coeff_file
                                                );

    this->_cd_trans     = new MorisonCoeffCurve(
                                                    input_->case_fopath,
                                                    morison_def_->transversal_drag_coeff_file
                                                );

    this->_cd_vert      = new MorisonCoeffCurve(
                                                    input_->case_fopath,
                                                    morison_def_->vertical_drag_coeff_file
                                                );

    // Initialize raddiff data for Morison element
    this->_rdd_morison  = new RadDiffData<cuscomplex, RDDMorisonConfig>(
                                                                            mpi_config_,
                                                                            this->_field_points_g.data( ),
                                                                            morison_def_->divisions_np,
                                                                            input_->angfreqs_np,
                                                                            input_->heads_np,
                                                                            input_->dofs_np
                                                                        );

}


MorisonElement::MorisonElement( 
                                    MpiConfig*          mpi_config_,
                                    Input*              input_,
                                    MorisonElementDef*  morison_def_
                                )
{
    // Calculate field points values for the Morison element
    this->_calculate_field_points( morison_def_ );

    // Initialize Morison element attributes and interpolation objects
    this->_initialize_element( mpi_config_, input_, morison_def_ );

}


MorisonElement::~MorisonElement( )
{
    if ( this->_ca_axial != nullptr )
    {
        delete this->_ca_axial;
        this->_ca_axial = nullptr;
    }

    if ( this->_ca_trans != nullptr )
    {
        delete this->_ca_trans;
        this->_ca_trans = nullptr;
    }

    if ( this->_ca_vert != nullptr )
    {
        delete this->_ca_vert;
        this->_ca_vert = nullptr;
    }

    if ( this->_cd_axial != nullptr )
    {
        delete this->_cd_axial;
        this->_cd_axial = nullptr;
    }

    if ( this->_cd_trans != nullptr )
    {
        delete this->_cd_trans;
        this->_cd_trans = nullptr;
    }

    if ( this->_cd_vert != nullptr )
    {
        delete this->_cd_vert;
        this->_cd_vert = nullptr;
    }

    if ( this->_rdd_morison != nullptr )
    {
        delete this->_rdd_morison;
        this->_rdd_morison = nullptr;
    }

}