
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


// Include standard headers
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>

// Include local modules
#include "hdf5_time_series_exporter.hpp"


// HDF5TimeSeriesExporterBase::HDF5TimeSeriesExporterBase(
//                                                             std::string folder_path,
//                                                             Mesh*       mesh
//                                                         )
//     :   _mesh( mesh ),
//         _folder_path( std::move( folder_path ) )
// {
//     this->_mesh_file  = this->_folder_path + "/mesh.h5";
//     this->_field_file = this->_folder_path + "/fields.h5";
//     this->_xdmf_file  = this->_folder_path + "/data.xdmf2";
// }


// std::size_t HDF5TimeSeriesExporterBase::nodes_count() const noexcept
// {
//     return static_cast<std::size_t>( this->_mesh ? this->_mesh->nodes_np : 0 );
// }


// HDF5TimeSeriesExporterBase::FieldDimensions
// HDF5TimeSeriesExporterBase::register_field_metadata(
//                                                         const   std::string&    field_name,
//                                                                 hsize_t         comps_np
//                                                     )
// {
//     FieldDimensions dims = { 1, this->nodes_count(), comps_np };
//     this->_field_names.push_back( field_name );
//     this->_field_dims[field_name] = dims;
//     return dims;
// }


// void HDF5TimeSeriesExporterBase::record_time_value( cusfloat append_value )
// {
//     this->_time.push_back( append_value );
//     this->_steps++;
// }


// void HDF5TimeSeriesExporterBase::write_mesh_payload()
// {
//     hid_t file = H5Fcreate(
//                                 this->_mesh_file.c_str(),
//                                 H5F_ACC_TRUNC,
//                                 H5P_DEFAULT,
//                                 H5P_DEFAULT
//                             );

//     H5Gcreate( file, "/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

//     std::size_t nodes_np = static_cast<std::size_t>( this->_mesh->nodes_np );
//     cusfloat* node_view_data = generate_empty_vector<cusfloat>( nodes_np * 3 );
//     for ( std::size_t i = 0; i < nodes_np; ++i )
//     {
//         node_view_data[3*i + 0] = this->_mesh->x[i];
//         node_view_data[3*i + 1] = this->_mesh->y[i];
//         node_view_data[3*i + 2] = this->_mesh->z[i];
//     }

//     write_hdf5_dataset_compressed(
//                                     file,
//                                     "/Mesh/Points",
//                                     cusfloat_h5,
//                                     { static_cast<hsize_t>( nodes_np ), 3 },
//                                     node_view_data
//                                 );

//     mkl_free( node_view_data );

//     std::vector<int> mixed;
//     int ernl = this->_mesh->enrl;
//     for ( size_t i = 0; i < static_cast<size_t>( this->_mesh->elems_np ); ++i )
//     {
//         if ( this->_mesh->elems_type[i] == 5 )
//         {
//             mixed.push_back( 4 );
//             mixed.push_back( this->_mesh->elems[ernl*i + 1] );
//             mixed.push_back( this->_mesh->elems[ernl*i + 2] );
//             mixed.push_back( this->_mesh->elems[ernl*i + 3] );
//             this->_tri_np++;
//         }
//         else if ( this->_mesh->elems_type[i] == 9 )
//         {
//             mixed.push_back( 5 );
//             mixed.push_back( this->_mesh->elems[ernl*i + 1] );
//             mixed.push_back( this->_mesh->elems[ernl*i + 2] );
//             mixed.push_back( this->_mesh->elems[ernl*i + 3] );
//             mixed.push_back( this->_mesh->elems[ernl*i + 4] );
//             this->_quads_np++;
//         }
//     }

//     write_hdf5_dataset_compressed(
//                                     file,
//                                     "/Mesh/MixedCells",
//                                     int_h5,
//                                     { static_cast<hsize_t>( mixed.size() ) },
//                                     mixed.data()
//                                 );

//     H5Fclose( file );
// }


// void HDF5TimeSeriesExporterBase::create_fields_container()
// {
//     hid_t file = H5Fcreate(
//                                 this->_field_file.c_str(),
//                                 H5F_ACC_TRUNC,
//                                 H5P_DEFAULT,
//                                 H5P_DEFAULT
//                             );

//     H5Gcreate( file, "/Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//     H5Fclose( file );
// }


// void HDF5TimeSeriesExporterBase::write_xdmf_payload() const
// {
//     std::string attr_type;
//     std::ofstream x( this->_xdmf_file );

//     x << "<?xml version=\"1.0\" ?>\n";
//     x << "<Xdmf Version=\"2.0\">\n";
//     x << "  <Domain>\n";
//     x << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

//     for ( size_t t = 0; t < this->_steps; ++t )
//     {
//         x << "      <Grid Name=\"Step" << t << "\" GridType=\"Uniform\">\n";
//         x << "        <Time Value=\"" << this->_time[t] << "\"/>\n";

//         if ( t == 0 )
//         {
//             x << "        <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << this->_mesh->elems_np << "\" BaseOffset=\"0\">\n";
//             x << "          <DataItem Format=\"HDF\" Dimensions=\"" << ( this->_tri_np * 4 + this->_quads_np * 5 ) << "\">\n";
//             x << "            mesh.h5:/Mesh/MixedCells\n";
//             x << "          </DataItem>\n";
//             x << "        </Topology>\n";

//             x << "        <Geometry GeometryType=\"XYZ\">\n";
//             x << "          <DataItem Format=\"HDF\" Dimensions=\"" << this->_mesh->nodes_np << " 3\">\n";
//             x << "            mesh.h5:/Mesh/Points\n";
//             x << "          </DataItem>\n";
//             x << "        </Geometry>\n";
//         }
//         else
//         {
//             x << "        <Topology Reference=\"/Xdmf/Domain/Grid/Grid[1]/Topology\"/>\n";
//             x << "        <Geometry Reference=\"/Xdmf/Domain/Grid/Grid[1]/Geometry\"/>\n";
//         }

//         for ( const auto& entry : this->_field_dims )
//         {
//             if ( entry.second[2] == 1 )
//             {
//                 attr_type = "Scalar";
//             }
//             else if ( entry.second[2] == 3 )
//             {
//                 attr_type = "Vector";
//             }
//             else
//             {
//                 throw std::runtime_error( "Unsupported number of components in field '" + entry.first + "'. Only scalar (1) and vector (3) fields are supported." );
//             }

//             x << "        <Attribute Name=\"" << entry.first << "\" AttributeType=\"" << attr_type << "\" Center=\"Node\">\n";
//             x << "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << this->_mesh->nodes_np << " " << entry.second[2] << "\">\n";
//             x << "            <DataItem Dimensions=\"3 3\" Format=\"XML\">\n";
//             x << "              " << t << " 0 0   1 1 1   1 " << this->_mesh->nodes_np << " " << entry.second[2] << "\n";
//             x << "            </DataItem>\n";
//             x << "            <DataItem Format=\"HDF\" Dimensions=\"" << this->_steps << " " << this->_mesh->nodes_np << " " << entry.second[2] << "\">\n";
//             x << "              fields.h5:/Fields/" << entry.first << "\n";
//             x << "            </DataItem>\n";
//             x << "          </DataItem>\n";
//             x << "        </Attribute>\n";
//         }

//         x << "      </Grid>\n";
//     }

//     x << "    </Grid>\n";
//     x << "  </Domain>\n";
//     x << "</Xdmf>\n";
// }


void HDF5TimeSeriesExporter::add_field(
                                            std::string             field_name,
                                            hsize_t                 samples_multiplier,
                                            hsize_t                 components_np
                                        )
{
    const auto write_dims = this->register_field_metadata( field_name, samples_multiplier, components_np );
    std::vector<hsize_t> initial_dims = { 0, write_dims[1], write_dims[2] };
    std::vector<hsize_t> max_dims     = { H5S_UNLIMITED, write_dims[1], write_dims[2] };

    hid_t file = H5Fopen( this->_field_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
    create_hdf5_extendable_dataset(
                                    file,
                                    std::string( "/Fields/" ) + field_name,
                                    cusfloat_h5,
                                    initial_dims,
                                    max_dims
                                );
    H5Fclose( file );
}


void HDF5TimeSeriesExporter::append_step(
                                            std::string                 field_name,
                                            cut::CusTensor<cusfloat>*   field_data
                                        )
{
    // Create file unit
    hid_t file = H5Fopen(
                            this->_field_file.c_str( ), 
                            H5F_ACC_RDWR, 
                            H5P_DEFAULT
                        );

    std::cout << "Appending field '" << field_name << "' for step " << this->_steps << " with dimensions: ";
    for ( auto& d : this->_field_dims[field_name] ) std::cout << d << " ";
    std::cout << std::endl;

    // Append field for current timestep
    append_to_hdf5_dataset(
                                file, 
                                std::string( "/Fields/" ) + field_name,
                                cusfloat_h5,
                                this->_field_dims[field_name],
                                field_data->data( )
                            );

    // Close file unit
    H5Fclose(file);
    
}


void HDF5TimeSeriesExporter::append_time( cusfloat append_value )
{
    this->record_time_value( append_value );
}


HDF5TimeSeriesExporter::HDF5TimeSeriesExporter(
                                                    const   std::string&    folder_path,
                                                            Mesh*           mesh
                                                )
    : HDF5TimeSeriesExporterBase( folder_path, mesh )
{
    this->write_mesh();
    this->initialize_time_series();
}


void HDF5TimeSeriesExporter::initialize_time_series(  )
{
    this->create_fields_container();
}


void HDF5TimeSeriesExporter::write_mesh( void )
{
    this->write_mesh_payload();
}


void HDF5TimeSeriesExporter::write_xdmf( void )
{
    this->write_xdmf_payload();
}


inline void example_hdf5_time_series_exporter( void )
{
    // Define output folder and mesh parameters
    std::string folder      = "path/to/output/folder";

    // Create mesh
    std::string mesh_fipath = "path/to/mesh/file.poly";
    cusfloat cog[3]         = {0.0, 0.0, 0.0};

    Mesh mesh( 
                    mesh_fipath,
                    "surface",
                    cog,
                    false,
                    PanelTypeE::DIFFRAC,
                    false
                );


    // Create exporter class
    HDF5TimeSeriesExporter exporter(
                                        folder,
                                        &mesh
                                    );

    

    // Add required fields to exporter (dataset names and dimensions must match the data passed to append_step())
    exporter.add_field( "Pressure", 1, 1 );
    exporter.add_field( "Velocity", 1, 3 );

    // Create auxiliary variables for synthetic solution
    const int nSteps    = 5;
    cusfloat T          = 2.0;
    cusfloat w          = 2.0 * PI / T;
    cusfloat k          = pow2s( w ) / 9.81;
    cusfloat time       = 0.0;
    cusfloat dt         = 0.1;

    cut::CusTensor<cusfloat> pressure( static_cast<std::size_t>( mesh.nodes_np ) );
    cut::CusTensor<cusfloat> velocity( { static_cast<std::size_t>( mesh.nodes_np ), 3 } );

    // Loop over time steps, fill fields with synthetic solution, and append to exporter
    for (int t = 0; t < nSteps; ++t)
    {
        time = t * dt;
        for ( std::size_t i=0; i<static_cast<std::size_t>(mesh.nodes_np); i++ )
        {
            pressure[i] = std::cos( k * mesh.x[i] - w * time ) * std::cos( k * mesh.y[i] - w * time );
            velocity(i, 0) = -k * std::sin( k * mesh.x[i] - w * time ) * std::cos( k * mesh.y[i] - w * time );
            velocity(i, 1) = -k * std::cos( k * mesh.x[i] - w * time ) * std::sin( k * mesh.y[i] - w * time );
            velocity(i, 2) = 0.0;

            
        }
        exporter.append_time( t * dt );
        exporter.append_step( "Pressure",  &pressure );
        exporter.append_step( "Velocity",  &velocity );
    }

    // Write XDMF companion file referencing the HDF5 datasets
    exporter.write_xdmf( );

    std::cout << "Compressed HDF5 time‑series written. Open data.xdmf in ParaView.\n";

}