
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
    this->_x_axis_l.resize( 3 );
    this->_x_axis_l[0] = rot_mat[0];
    this->_x_axis_l[1] = rot_mat[3];
    this->_x_axis_l[2] = rot_mat[6];

    this->_y_axis_l.resize( 3 );
    this->_y_axis_l[0] = rot_mat[1];
    this->_y_axis_l[1] = rot_mat[4];
    this->_y_axis_l[2] = rot_mat[7];

    this->_z_axis_l.resize( 3 );
    this->_z_axis_l[0] = rot_mat[2];
    this->_z_axis_l[1] = rot_mat[5];
    this->_z_axis_l[2] = rot_mat[8];

    // Calculate divison segement length
    cusfloat div_length = morison_def_->length / static_cast<cusfloat>( morison_def_->divisions_np );

    // Calculate nodes positions along the element axis
    cut::CusTensor<cusfloat> nodes_pos( { morison_def_->divisions_np + 1, 3 } );
    for ( std::size_t i = 0; i < morison_def_->divisions_np + 1; i++ )
    {
        nodes_pos( i, 0 ) = morison_def_->start_pos[0] + this->_x_axis_l[0] * div_length * static_cast<cusfloat>( i );
        nodes_pos( i, 1 ) = morison_def_->start_pos[1] + this->_x_axis_l[1] * div_length * static_cast<cusfloat>( i );
        nodes_pos( i, 2 ) = morison_def_->start_pos[2] + this->_x_axis_l[2] * div_length * static_cast<cusfloat>( i );
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
    std::array<cuscomplex, 3> rao_disp_fp = { 0.0, 0.0, 0.0 };
    std::array<cuscomplex, 3> rao_vel_fp  = { 0.0, 0.0, 0.0 };
    std::array<cuscomplex, 3> rao_acc_fp  = { 0.0, 0.0, 0.0 };

    for ( std::size_t ih = 0; ih < this->_input->heads_np; ih++ )
    {
        for ( std::size_t i = 0; i < this->_field_points_np; i++ )
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

            // Calculate inertial and drag forces at the field point
            const cut::CusTensor<cusfloat>* vel_x = this->_rdd_morison->template get_field_data<FieldComponentE::TOTAL, FieldTypeE::VELOCITY_X>( ih );
            const cut::CusTensor<cusfloat>* vel_y = this->_rdd_morison->template get_field_data<FieldComponentE::TOTAL, FieldTypeE::VELOCITY_Y>( ih );
            const cut::CusTensor<cusfloat>* vel_z = this->_rdd_morison->template get_field_data<FieldComponentE::TOTAL, FieldTypeE::VELOCITY_Z>( ih );

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
    this->_area_axial   = morison_def_->area_axial;
    this->_area_trans   = morison_def_->area_transversal;
    this->_area_vert    = morison_def_->area_vertical;
    this->_kc_length    = morison_def_->kc_length;

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