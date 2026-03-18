
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
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

// Include local modules
#include "input.hpp"
#include "reader.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"
#include "../tools.hpp"
#include "../version.hpp"
#include "../waves/wave_dispersion_base_fo.hpp"

// Manage namespaces
namespace fs = std::filesystem;


Input::Input( const std::string& folder_path )
{
    load( folder_path );
}


void Input::load( const std::string& folder_path )
{
    // Save case folder path
    this->case_fopath = folder_path;
    this->folder_path = folder_path;

    // Read case.input.dat file
    this->read_case( folder_path );

    // Read bodies
    this->read_bodies( folder_path );

    // Read field points
    this->read_field_points( folder_path );

    // Configure inputs
    this->configure( );
}


void Input::read_bodies( const std::string& folder_path )
{
    // Allocate space for the bodies description
    this->bodies = new BodyDef*[this->bodies_np];

    // Loop over body definitions to load them
    for ( int i=0; i<this->bodies_np; i++ )
    {
        // Create new body instance
        this->bodies[i] = new BodyDef;

        // Load body data
        read_body(
                    folder_path,
                    this->bodies_finame[i],
                    this->bodies[i]
                );

    }

    // Set bodies as read
    this->is_bodies = true;
}


void Input::read_body(
                            const std::string& folder_path,
                            const std::string& target_file,
                            BodyDef* body
                    )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Generate case.input.dat file path
    fs::path folder_path_( folder_path );
    fs::path file_name( target_file );
    fs::path file_path = folder_path_ / file_name;

    // Open file unit
    std::ifstream infile;
    infile.open( file_path, std::ios_base::in );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Read file header line
    skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

    //////////////////////////////////////////////
    /*********** General Properties *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read Is Fix
    target_signal   = "IsFix";
    read_signal     = read_channel_value( infile, body->is_fix );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************* Mass Properties **************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read mass
    target_signal   = "BodyMass";
    read_signal     = read_channel_value( infile, body->mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGX
    target_signal   = "COGX";
    read_signal     = read_channel_value( infile, body->cog[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGY
    target_signal   = "COGY";
    read_signal     = read_channel_value( infile, body->cog[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGZ
    target_signal   = "COGZ";
    read_signal     = read_channel_value( infile, body->cog[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read switch to define how the intertia is defined, by the 
    // radius of gyrantion or by the inertial matrix components
    target_signal   = "IBR";
    read_signal     = read_channel_value( infile, body->interia_by_rad );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXX
    target_signal   = "IXX";
    read_signal     = read_channel_value( infile, body->inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXY
    target_signal   = "IXY";
    read_signal     = read_channel_value( infile, body->inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXZ
    target_signal   = "IXZ";
    read_signal     = read_channel_value( infile, body->inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYY
    target_signal   = "IYY";
    read_signal     = read_channel_value( infile, body->inertia[3] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYZ
    target_signal   = "IYZ";
    read_signal     = read_channel_value( infile, body->inertia[4] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IZZ
    target_signal   = "IZZ";
    read_signal     = read_channel_value( infile, body->inertia[5] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RXX
    target_signal   = "RXX";
    read_signal     = read_channel_value( infile, body->rad_inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RYY
    target_signal   = "RYY";
    read_signal     = read_channel_value( infile, body->rad_inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RZZ
    target_signal   = "RZZ";
    read_signal     = read_channel_value( infile, body->rad_inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    if ( body->interia_by_rad )
    {
        body->inertia[0] = body->mass * pow2s( body->rad_inertia[0] );
        body->inertia[3] = body->mass * pow2s( body->rad_inertia[1] );
        body->inertia[5] = body->mass * pow2s( body->rad_inertia[2] );
    }

    //////////////////////////////////////////////
    /************* Mesh Description *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read mesh body name
    target_signal   = "BodyName";
    read_signal     = read_channel_value( infile, body->mesh_body_name );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read mesh file name
    target_signal   = "MeshFile";
    read_signal     = read_channel_value( infile, body->mesh_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read usage of internal lid
    target_signal   = "LidType";
    read_signal     = read_channel_value( infile, body->lid_type );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read usage of external lid
    target_signal   = "UseExtLid";
    read_signal     = read_channel_value( infile, body->use_ext_lid );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read external lid damping factor
    target_signal   = "ExtLidDF";
    read_signal     = read_channel_value( infile, body->ext_lid_damp_f );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /********** Morison Coefficients ************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read water line points detection precision
    target_signal   = "UseMC";
    read_signal     = read_channel_value( infile, body->use_morison_elements );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Loop over morison elements files in the morison_elements folder and 
    // load them if the user has defined the use of morison elements
    if ( !body->use_morison_elements )
    {
        body->morison_elements_names.clear( );
        body->morison_elements_np = 0;
    }
    else
    {
        body->morison_elements_names.clear( );
        fs::path morison_dir = folder_path_ / "morison_elements";

        auto has_me_extension = []( const std::string& name )
        {
            const std::string ext = ".me.dat";
            if ( name.size( ) < ext.size( ) )
            {
                return false;
            }
            return name.compare( name.size( ) - ext.size( ), ext.size( ), ext ) == 0;
        };

        if ( fs::exists( morison_dir ) && fs::is_directory( morison_dir ) )
        {
            for ( const auto& entry : fs::directory_iterator( morison_dir ) )
            {
                if ( !entry.is_regular_file( ) )
                {
                    continue;
                }
                const std::string filename = entry.path( ).filename( ).string( );
                if ( has_me_extension( filename ) )
                {
                    body->morison_elements_names.push_back( filename );
                }
            }
        }

        std::sort( body->morison_elements_names.begin( ), body->morison_elements_names.end( ) );
        body->morison_elements_np = body->morison_elements_names.size( );

        // Load morison elements
        body->morison_elements.reserve( body->morison_elements_np );
        for ( std::size_t i=0; i<body->morison_elements_np; i++ )
        {
            body->morison_elements.emplace_back(
                                                    MorisonElementDef(
                                                                            morison_dir.string( ),
                                                                            body->morison_elements_names[i],
                                                                            body->cog
                                                                        )
                                                );
        }
    }

    // Load mesh for the current body
    std::vector<Mesh*> _mesh_total;
    fs::path mesh_foname_{ MESH_FOLDER_NAME };
    fs::path mesh_finame_{ body->mesh_finame };
    fs::path mesh_fipath_ = folder_path_ / mesh_foname_ / mesh_finame_;
    body->mesh          = new   Mesh( 
                                        mesh_fipath_.string( ),
                                        body->mesh_body_name,
                                        body->cog,
                                        body->is_fix,
                                        PanelTypeE::DIFFRAC                                
                                    );
    body->is_mesh       = true;
    body->mesh_items_np = 1;
    _mesh_total.push_back( body->mesh );
    
    // Write mesh
    fs::path mesh_foname_1_{ RESULTS_MESH_FOLDER_NAME };
    fs::path results_foname_{ RESULTS_FOLDER_NAME };
    fs::path results_mesh_fipath_ = folder_path_ / results_foname_;
    fs::path plot_mesh_fipath_ = results_mesh_fipath_ / mesh_foname_1_;
    
    if ( !fs::exists( results_mesh_fipath_) )
    {
        fs::create_directory( results_mesh_fipath_ );
    }

    if ( !fs::exists( plot_mesh_fipath_) )
    {
        fs::create_directory( plot_mesh_fipath_ );
    }

    body->mesh->write( plot_mesh_fipath_.string( ) );

    // Load lid if it is required by the user
    if ( body->lid_type == 1 )
    {  
        // Define lid name
        std::stringstream ss;
        ss << body->mesh_body_name << "_lid";
        
        // Load lid mesh
        body->mesh_int_lid  = new Mesh(
                                            mesh_fipath_.string( ),
                                            ss.str( ),
                                            body->cog,
                                            body->is_fix,
                                            PanelTypeE::INT_LID
                                        );
        body->is_mesh_int_lid   = true;
        body->mesh_items_np++;
        _mesh_total.push_back( body->mesh_int_lid );

    }

    // Load external lid if it is required by the user
    if ( body->use_ext_lid )
    {
        // Define lid name
        std::stringstream sse;
        sse << body->mesh_body_name << "_ext_lid";
        
        // Load lid mesh
        body->mesh_ext_lid  = new Mesh(
                                            mesh_fipath_.string( ),
                                            sse.str( ),
                                            body->cog,
                                            body->is_fix,
                                            PanelTypeE::EXT_LID
                                        );
        body->mesh_items_np++;

        // Apply external lid damping factor
        for ( std::size_t i=0; i<static_cast<std::size_t>( body->mesh_ext_lid->elems_np ); i++ )
        {
            body->mesh_ext_lid->panels[i]->set_ext_lid_damp_f( body->ext_lid_damp_f );
        }

        // Append to total mesh vector
        _mesh_total.push_back( body->mesh_ext_lid );
    }

    // Load QTF free surface
    if ( this->out_qtf_so_model > 0 )
    {
        // Define QTF free surface name
        std::stringstream ss;
        ss << body->mesh_body_name << "_fs_qtf";

        // Load QTF free surface mesh
        body->mesh_fs_qtf       = new Mesh(
                                            mesh_fipath_.string( ),
                                            ss.str( ),
                                            body->cog,
                                            body->is_fix,
                                            PanelTypeE::QTF_LID
                                        );

        // Calculate radius of the fs mesh
        body->mesh_fs_qtf->calculate_fs_radius( );

        body->is_mesh_fs_qtf    = true;
        body->mesh_items_np++;
    }

    // Compose total mesh
    body->mesh_total    = new Mesh(
                                        body->mesh->name, 
                                        _mesh_total, 
                                        body->cog, 
                                        body->is_fix
                                    );
    body->is_mesh_total = true;

    body->mesh_total->write( plot_mesh_fipath_.string( ), "total" );

    // Close file unit
    infile.close();

}


void Input::read_case( const std::string& folder_path )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_file     = "case.input.dat";
    std::string target_signal   = "";
    std::string _version        = "";

    // Generate case.input.dat file path
    fs::path folder_path_( folder_path );
    fs::path file_name( "case.input.dat" );
    fs::path file_path = folder_path_ / file_name;

    // Open file unit
    std::ifstream infile;
    infile.open( file_path, std::ios_base::in );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Read file header line
    skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

    //////////////////////////////////////////////
    /************** Solver Controls *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read QTF second order potential model
    target_signal   = "QTFSOModel";
    read_signal     = read_channel_value( infile, this->out_qtf_so_model );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read water line points detection precision
    target_signal   = "WLDetPrec";
    read_signal     = read_channel_value( infile, this->wl_det_prec );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************** Body Definition *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read bodies name definition
    target_signal = "BodyFN";
    read_channel_list(
                            infile,
                            target_file,
                            target_signal,
                            line_count,
                            this->bodies_finame
                        );
    this->bodies_np = this->bodies_finame.size( );

    //////////////////////////////////////////////
    /********* Field Points Definition **********/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read water line points detection precision
    target_signal   = "UseFP";
    read_signal     = read_channel_value( infile, this->use_field_points );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read bodies name definition
    target_signal = "FieldPN";
    read_channel_list(
                            infile,
                            target_file,
                            target_signal,
                            line_count,
                            this->field_points_finame
                        );

    if ( !this->use_field_points )
    {
        this->field_points_finame.clear( );
    }

    this->field_points_np = this->field_points_finame.size( );

    //////////////////////////////////////////////
    /************** Output Channels *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read output flag for diffraction force
    target_signal   = "OutDiffrac";
    read_signal     = read_channel_value( infile, this->out_diffrac );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for Froude-Krylov force
    target_signal   = "OutFK";
    read_signal     = read_channel_value( infile, this->out_fk );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for hydromechanic coefficients
    target_signal   = "OutHydMech";
    read_signal     = read_channel_value( infile, this->out_hydmech );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for hydrostatic stiffness matrix
    target_signal   = "OutHydStiff";
    read_signal     = read_channel_value( infile, this->out_hydstiff );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for potential over the panels
    target_signal   = "OutPot";
    read_signal     = read_channel_value( infile, this->out_potential );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for pressure over the panels
    target_signal   = "OutPress";
    read_signal     = read_channel_value( infile, this->out_pressure );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for mean drift
    target_signal   = "OutMDrift";
    read_signal     = read_channel_value( infile, this->out_mdrift );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for mesh
    target_signal   = "OutMesh";
    read_signal     = read_channel_value( infile, this->out_mesh );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for Morison elements forces
    target_signal   = "OutMorison";
    read_signal     = read_channel_value( infile, this->out_morison );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for quadratic transfer functions
    target_signal   = "OutQTF";
    read_signal     = read_channel_value( infile, this->out_qtf );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for quadratic transfer functions components
    target_signal   = "OutQTFComp";
    read_signal     = read_channel_value( infile, this->out_qtf_comp );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for response amplitude operator
    target_signal   = "OutRAOs";
    read_signal     = read_channel_value( infile, this->out_raos );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for sources intensity over the panels
    target_signal   = "OutSources";
    read_signal     = read_channel_value( infile, this->out_sources );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for structural mass
    target_signal   = "OutStMass";
    read_signal     = read_channel_value( infile, this->out_struct_mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for wave excitation forces
    target_signal   = "OutWex";
    read_signal     = read_channel_value( infile, this->out_wex );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    
    //////////////////////////////////////////////
    /************** Site Conditions *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read water density
    target_signal   = "RhoW";
    read_signal     = read_channel_value( infile, this->water_density );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    
    // Read gravitational acceleration
    target_signal   = "GravAcc";
    read_signal     = read_channel_value( infile, this->grav_acc );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read end simulation time
    target_signal   = "WaterDepth";
    read_signal     = read_channel_value( infile, this->water_depth );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /****************** Heading *****************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read headings units
    target_signal   = "HeadUnits";
    read_signal     = read_channel_value( infile, this->heads_units );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read list contraction
    read_list_contraction(
                                infile,
                                line_count,
                                target_file,
                                this->heads
                            );
    this->heads_np = this->heads.size( );

    //////////////////////////////////////////////
    /*************** Frequencies ****************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read headings units
    target_signal   = "FreqUnit";
    read_signal     = read_channel_value( infile, this->freqs_unit );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read list contraction
    read_list_contraction(
                                infile,
                                line_count,
                                target_file,
                                this->angfreqs
                            );
    this->angfreqs_np = this->angfreqs.size( );
 
    // Close file unit
    infile.close();
}

void Input::read_field_points( 
                                    const std::string& folder_path 
                                )
{
    // Loop over field points definitions to load them
    this->field_points.reserve( this->field_points_np );
    for ( std::size_t i=0; i<this->field_points_np; i++ )
    {
        // Create new field points instance
        this->field_points.emplace_back( 
                                            FieldPointsDef(
                                                                folder_path,
                                                                this->field_points_finame[i]
                                                            ) 
                                        );

    }
    
}

void Input::configure( void )
{
    /**********************************************************/
    /****************** Check input flags *********************/
    /**********************************************************/

    // Check if fast mode is configured properly
    if ( 
            this->poly_order > 0 
            &&
            this->is_fast_solver
        )
    {
        std::cerr << "Using Fast Mode with high order polynomials ";
        std::cerr << "( > 0) is not allowed!" << std::endl;
        throw std::runtime_error( "" );
    }

    // Check if it is necessary to calculate Mean Drift values
    // This option will allow the frequency calculation of the 
    // QTF depedent signals
    this->is_calc_mdrift = ( this->out_mdrift || this->out_qtf );

    // Check if it is necessary to load the free surface 
    // QTF mesh
    this->is_fs_qtf = this->out_qtf_so_model > 0 ? true: false;

    /**********************************************************/
    /************** Check headings input units ****************/
    /**********************************************************/
    if ( this->heads_units.compare( "deg" ) == 0 )
    {
        for ( int i=0; i<this->heads_np; i++ )
        {
            this->heads[i] = deg_to_rad( this->heads[i] );
        }
    }
    else if ( this->heads_units.compare( "rad" ) != 0 )
    {
        std::cout << std::endl;
        std::cout << "ERROR - INPUT:" << std::endl;
        std::cout << "HeadUnits: " << this->heads_units << " is not a valid parameter." << std::endl;
        std::cout << "Valid heading units are: deg | rad." << std::endl;
        std::cout << std::endl;
        throw std::runtime_error( "" );
    }

    /**********************************************************/
    /************** Check input frequencies units *************/
    /**********************************************************/
    if ( this->freqs_unit.compare( "period" ) == 0 )
    {
        for ( int i=0; i<this->angfreqs_np; i++ )
        {
            this->angfreqs[i] = period_to_angfreq( this->angfreqs[i] );
        }
    }
    else if ( this->freqs_unit.compare( "freq" ) == 0 )
    {
        for ( int i=0; i<this->angfreqs_np; i++ )
        {
            this->angfreqs[i] = freq_to_angfreq( this->angfreqs[i] );
        }
    }
    else if ( this->freqs_unit.compare( "angfreq" ) != 0 )
    {
        std::cout << std::endl;
        std::cout << "ERROR - INPUT:" << std::endl;
        std::cout << "FreqUnit: " << this->freqs_unit << " is not a valid parameter." << std::endl;
        std::cout << "Valid frequency units are: period | freq | angfreq." << std::endl;
        std::cout << std::endl;
        throw std::runtime_error( "" );
    }

    /**********************************************************/
    /**** Sort frequencies from the lowest to the highest *****/
    /**********************************************************/
    int* sort_keys  = generate_empty_vector<int>( this->angfreqs_np );
    int  info       = 0;
    lasrt2<cusfloat>( 
                        "I",
                        &this->angfreqs_np,
                        this->angfreqs.data(),
                        sort_keys,
                        &info
                    );

    if ( info != 0 )
    {
        std::cerr << "ERROR - lasrt2" << std::endl;
        std::cerr << "Sort algorithm could not sort the input angular frequencies.";
        std::cerr << " - Error Code: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    mkl_free( sort_keys );

    /**********************************************************/
    /********** Create a vector for the frequencies ***********/
    /**********************************************************/
    this->freqs = generate_empty_vector<cusfloat>( this->angfreqs_np );
    for ( int i=0; i<this->angfreqs_np; i++ )
    {
        this->freqs[i] = angfreq_to_freq( this->angfreqs[i] );
    }

    /**********************************************************/
    /********** Detect points over the free surface ***********/
    /**********************************************************/
    this->is_wl_points = this->is_calc_mdrift;
    if ( 
            this->is_wl_points
        )
    {
        // Detect points over the free surface
        for ( int i=0; i<this->bodies_np; i++ )
        {
            this->bodies[i]->mesh_total->detect_wl_points( this->wl_det_prec );
        }
    }

    /**********************************************************/
    /**************** Calculate source nodes ******************/
    /**********************************************************/
    for ( int i=0; i<this->bodies_np; i++ )
    {
        this->bodies[i]->mesh_total->define_source_nodes(
                                                        this->poly_order,
                                                        this->bodies[i]->cog
                                                    );
    }

}


int  Input::gauss_np_factor_1d( void )
{
    int gf = this->gauss_order;
    // if ( !this->is_block_adaption )
    // {
    //     gf = 1;
    // }

    return gf;
}


int  Input::gauss_np_factor_2d( void )
{
    int gf = pow2s( this->gauss_order );
    // if ( !this->is_block_adaption )
    // {
    //     gf = 1;
    // }

    return gf;
}

Input::~Input( void )
{
    if ( this->is_bodies )
    {
        // Delete frequencies containers
        mkl_free( this->freqs );

        // Delete BodyDef object instances
        for ( int i=0; i<this->bodies_np; i++ )
        {
            delete this->bodies[i];
        }
        
        // Delete vector of BodyDef pointers
        delete [] this->bodies;

    }
}


void Input::print( void )
{
    std::cout << std::endl;
    std::cout << "SOLVER CONTROLS: " << std::endl;
    std::cout << " - Block Adaption: " << this->is_block_adaption << std::endl;
    std::cout << " - Fast Solver: " << this->is_fast_solver << std::endl;
    std::cout << " - Gauss Order: " << this->gauss_order << std::endl;
    std::cout << " - GFDnAbsErr: " << this->gfdn_abs_err << std::endl;
    std::cout << " - GFDnRelErr: " << this->gfdn_rel_err << std::endl;
    std::cout << " - KochinNC: " << this->kochin_np << std::endl;
    std::cout << " - LogSingAna: " << this->is_log_sin_ana << std::endl;
    std::cout << " - PolyOrder: " << this->poly_order << std::endl;
    std::cout << " - PotAbsErr: " << this->pot_abs_err << std::endl;
    std::cout << " - PotRelErr: " << this->pot_rel_err << std::endl;
    std::cout << " - PressAbsErr: " << this->press_abs_err << std::endl;
    std::cout << " - PressRelErr: " << this->press_rel_err << std::endl;
    std::cout << " - QTFSOModel: " << this->out_qtf_so_model << std::endl;
    std::cout << " - WLDetPrec: " << this->wl_det_prec << std::endl;

    std::cout << std::endl;
    std::cout << "BODY DEFINITION: " << std::endl;
    for ( int i=0; i<this->bodies_np; i++ )
    {
        std::cout << " - Body " << i << ": " << this->bodies_finame[i] << std::endl;
        this->bodies[i]->print( );
    }

    std::cout << std::endl;
    std::cout << "OUTPUT CHANNELS: " << std::endl;
    std::cout << " - OutDiffrac: " << this->out_diffrac << std::endl;
    std::cout << " - OutFK: " << this->out_fk << std::endl;
    std::cout << " - OutHydMech: " << this->out_hydmech << std::endl;
    std::cout << " - OutHydStiff: " << this->out_hydstiff << std::endl;
    std::cout << " - OutPress: " << this->out_pressure << std::endl;
    std::cout << " - OutMDrift: " << this->out_mdrift << std::endl;
    std::cout << " - OutMesh: " << this->out_mesh << std::endl;
    std::cout << " - OutQTF: " << this->out_qtf << std::endl;
    std::cout << " - OutQTFComp: " << this->out_qtf_comp << std::endl;
    std::cout << " - OutRAOs: " << this->out_raos << std::endl;
    std::cout << " - OutSources: " << this->out_sources << std::endl;
    std::cout << " - OutStMass: " << this->out_struct_mass << std::endl;
    std::cout << " - OutWex: " << this->out_wex << std::endl;

    std::cout << std::endl;
    std::cout << "SITE CONDITIONS: " << std::endl;
    std::cout << " - Water Density: " << this->water_density << std::endl;
    std::cout << " - Grav. Acc: " << this->grav_acc << std::endl;
    std::cout << " - Water Depth: " << this->water_depth << std::endl;

    std::cout << std::endl;
    std::cout << "HEADINGS: " << std::endl;
    print_vector( this->heads_np, this->heads.data( ), 1, 6 );

    std::cout << std::endl;
    std::cout << "ANGULAR FREQUENCIES: " << std::endl;
    print_vector( this->angfreqs_np, this->angfreqs.data( ), 1, 6 );
    std::cout << std::endl;
}


void Input::get_wave_quality_params(
                                        cusfloat& min_wavelength,
                                        cusfloat& max_wave_number
                                    ) const
{
    min_wavelength = 0.0;
    max_wave_number = 0.0;

    if ( this->angfreqs_np <= 0 || this->water_depth <= 0.0 || this->grav_acc <= 0.0 )
    {
        return;
    }

    cusfloat min_lambda = 1e30;
    cusfloat max_k = 0.0;
    for ( int i = 0; i < this->angfreqs_np; i++ )
    {
        const cusfloat k = w2k( this->angfreqs[i], this->water_depth, this->grav_acc );
        if ( k > 0.0 )
        {
            const cusfloat lambda = 2.0 * PI / k;
            if ( lambda < min_lambda )
            {
                min_lambda = lambda;
            }
            if ( k > max_k )
            {
                max_k = k;
            }
        }
    }

    if ( min_lambda < 1e29 )
    {
        min_wavelength = min_lambda;
    }
    max_wave_number = max_k;
}