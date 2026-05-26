
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
#include <cstdint>
#include <fstream>
#include <utility>
#include <vector>

// Include local modules
#include "../containers/logger.hpp"
#include "vtu.hpp"


// -----------------------------
// Binary VTU writer
// -----------------------------
bool    write_vtu_binary_appended(
                                    const std::string       &filename,
                                    const std::size_t       nodes_np,
                                    cusfloat*               nodes_x,
                                    cusfloat*               nodes_y,
                                    cusfloat*               nodes_z,
                                    const std::size_t       elems_np,
                                    const int               enrl,
                                    int*                    elements,
                                    int*                    elements_type
                                )
{
    // Open file unit
    std::ofstream out( filename, std::ios::binary );
    if ( !out.is_open( ) ) 
    {
        Logger logger;
        logger.error( "Output mesh file for binary VTU could not be created!\n" );
        return false;
    }

    //--------------------------------------------------------------------
    // Precompute binary blocks and their sizes
    //--------------------------------------------------------------------

    // Points
    const uint32_t points_bytes = nodes_np * 3 * sizeof( cusfloat );
    std::vector<cusfloat> points_bin;
    points_bin.reserve( nodes_np * 3 );
    for ( std::size_t i=0; i<nodes_np; i++ ) 
    {
        points_bin.push_back( nodes_x[i] );
        points_bin.push_back( nodes_y[i] );
        points_bin.push_back( nodes_z[i] );
    }

    // Connectivity
    std::vector<int32_t> conn_bin;
    for ( std::size_t i=0; i<elems_np; i++ )
    {
        for ( int j=0; j<elements[ i*enrl+0 ]; j++ )
        {
            conn_bin.push_back( elements[ i*enrl + 1 + j ] );   // zero based
        }
    }

    const uint32_t conn_bytes = conn_bin.size( ) * sizeof( int32_t );

    // Offsets
    std::vector<int32_t> offs_bin;
    offs_bin.reserve( elems_np );
    int32_t offset = 0;
    for ( std::size_t i=0; i<elems_np; i++ ) 
    {
        offset += elements[ i*enrl + 0 ];
        offs_bin.push_back( offset );
    }
    const uint32_t offs_bytes = offs_bin.size() * sizeof( int32_t );

    // Types
    std::vector<uint8_t> type_bin;
    type_bin.reserve( elems_np );
    for ( std::size_t i=0; i<elems_np; i++ )
    {
        type_bin.push_back( elements_type[i] );
    }

    const uint32_t type_bytes = type_bin.size( ) * sizeof( uint8_t );

    //--------------------------------------------------------------------
    // Compute offsets inside appended section
    //--------------------------------------------------------------------
    uint32_t offset_points  = 0;    
    uint32_t offset_conn    = offset_points + sizeof( uint32_t ) + points_bytes;
    uint32_t offset_offs    = offset_conn   + sizeof( uint32_t ) + conn_bytes;
    uint32_t offset_types   = offset_offs   + sizeof( uint32_t ) + offs_bytes;

    //--------------------------------------------------------------------
    // Write XML header
    //--------------------------------------------------------------------
    out <<
R"(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">
  <UnstructuredGrid>
    <Piece NumberOfPoints=")" << nodes_np <<
R"(" NumberOfCells=")" << elems_np <<
R"(">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset=")" << offset_points << R"("/>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="appended" offset=")" << offset_conn << R"("/>
        <DataArray type="Int32" Name="offsets" format="appended" offset=")" << offset_offs << R"("/>
        <DataArray type="UInt8" Name="types" format="appended" offset=")" << offset_types << R"("/>
      </Cells>
    </Piece>
  </UnstructuredGrid>
  <AppendedData encoding="raw">_)";

    //--------------------------------------------------------------------
    // Write appended binary blocks
    //--------------------------------------------------------------------

    // Points
    out.write( (char*)&points_bytes, sizeof( uint32_t ) );
    out.write( (char*)points_bin.data( ), points_bytes  );

    // Connectivity
    out.write( (char*)&conn_bytes, sizeof( uint32_t )   );
    out.write( (char*)conn_bin.data( ), conn_bytes      );

    // Offsets
    out.write( (char*)&offs_bytes, sizeof( uint32_t )   );
    out.write( (char*)offs_bin.data( ), offs_bytes      );

    // Types
    out.write( (char*)&type_bytes, sizeof( uint32_t )   );
    out.write( (char*)type_bin.data( ), type_bytes      );

    //--------------------------------------------------------------------
    // Close VTU
    //--------------------------------------------------------------------
    out << "</AppendedData>\n</VTKFile>\n";
    out.close( );
    return true;
}


// -----------------------------
// ASCII VTU writer
// -----------------------------
bool    write_vtu_ascii(
                            const std::string       &filename,
                            const std::size_t       nodes_np,
                            cusfloat*               nodes_x,
                            cusfloat*               nodes_y,
                            cusfloat*               nodes_z,
                            const std::size_t       elems_np,
                            const int               enrl,
                            int*                    elements,
                            int*                    elements_type
                        )
{
    // Open file unit
    std::ofstream out( filename );
    if ( !out.is_open( ) ) 
    {
        Logger logger;
        logger.error( "Output mesh file for binary VTU could not be created!\n" );
        return false;
    }

    // Print header
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << nodes_np
        << "\" NumberOfCells=\"" << elems_np << "\">\n";

    // Points
    out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for ( std::size_t i=0; i<nodes_np; i++ )
    {
        out << nodes_x[i] << " " << nodes_y[i] << " " << nodes_z[i] << "\n";
    }
    out << "</DataArray>\n</Points>\n";

    // Connectivity
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for ( std::size_t i=0; i<elems_np; i++ )
    {
        for ( int j=0; j<elements[ i*enrl + 0 ]; j++ )
        {
            out << elements[ i*enrl + 1 + j ] << " ";
        }
        out << "\n";
    }
    out << "</DataArray>\n";

    // Offsets (cumulative node counts)
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for ( std::size_t i=0; i<elems_np; i++ ) 
    {
        offset += elements[ i*enrl + 0 ];
        out << offset << "\n";
    }
    out << "</DataArray>\n";

    // VTK cell types
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( std::size_t i=0; i<elems_np; i++ )
    {
        out << elements_type[i] << "\n";
    }
    out << "</DataArray>\n";

    out << "</Cells>\n";
    out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";

    return true;
}


// -----------------------------------------------------------------------
// Compile-time helper: VTK type-name string for cusfloat
// -----------------------------------------------------------------------
static constexpr const char* vtk_float_type_name( )
{
    return ( sizeof( cusfloat ) == 4 ) ? "Float32" : "Float64";
}


// -----------------------------------------------------------------------
// write_vtu_panel_pressure  — binary-appended UnstructuredGrid + CellData
// -----------------------------------------------------------------------
bool    write_vtu_panel_pressure(
                                    const std::string&      filename,
                                    std::size_t             n_nodes,
                                    const cusfloat*         nodes_x,
                                    const cusfloat*         nodes_y,
                                    const cusfloat*         nodes_z,
                                    std::size_t             n_cells,
                                    const int32_t*          connectivity,
                                    const int32_t*          offsets,
                                    const uint8_t*          types,
                                    const cusfloat*         pressure,
                                    const cusfloat*         phi_dt_comp,
                                    const cusfloat*         kinetic_comp,
                                    const cusfloat*         hydrostatic_comp
                                )
{
    std::ofstream out( filename, std::ios::binary );
    if ( !out.is_open( ) )
    {
        Logger logger;
        logger.error( "Could not create VTU pressure file: " + filename + "\n" );
        return false;
    }

    const char* ftype = vtk_float_type_name( );

    // ----------------------------------------------------------------
    // Binary block sizes
    // ----------------------------------------------------------------
    const uint32_t pts_bytes  = static_cast<uint32_t>( n_nodes  * 3 * sizeof( cusfloat ) );
    const uint32_t conn_bytes = static_cast<uint32_t>( n_nodes      * sizeof( int32_t  ) );
    const uint32_t offs_bytes = static_cast<uint32_t>( n_cells      * sizeof( int32_t  ) );
    const uint32_t type_bytes = static_cast<uint32_t>( n_cells      * sizeof( uint8_t  ) );
    const uint32_t pres_bytes = static_cast<uint32_t>( n_cells      * sizeof( cusfloat ) );

    // ----------------------------------------------------------------
    // Appended-data section offsets (each block starts with its uint32 length)
    // ----------------------------------------------------------------
    const uint32_t off_pts        = 0;
    const uint32_t off_conn       = off_pts        + sizeof( uint32_t ) + pts_bytes;
    const uint32_t off_offs       = off_conn       + sizeof( uint32_t ) + conn_bytes;
    const uint32_t off_type       = off_offs       + sizeof( uint32_t ) + offs_bytes;
    const uint32_t off_pres       = off_type       + sizeof( uint32_t ) + type_bytes;
    const uint32_t off_phi_dt     = off_pres       + sizeof( uint32_t ) + pres_bytes;
    const uint32_t off_kinetic    = off_phi_dt     + sizeof( uint32_t ) + pres_bytes;
    const uint32_t off_hydrostatic= off_kinetic    + sizeof( uint32_t ) + pres_bytes;

    // ----------------------------------------------------------------
    // XML header
    // ----------------------------------------------------------------
    out <<
R"(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">
  <UnstructuredGrid>
    <Piece NumberOfPoints=")" << n_nodes << R"(" NumberOfCells=")" << n_cells << R"(">
      <Points>
        <DataArray type=")" << ftype << R"(" NumberOfComponents="3" format="appended" offset=")" << off_pts << R"("/>
      </Points>
      <Cells>
        <DataArray type="Int32"  Name="connectivity" format="appended" offset=")" << off_conn << R"("/>
        <DataArray type="Int32"  Name="offsets"      format="appended" offset=")" << off_offs << R"("/>
        <DataArray type="UInt8"  Name="types"        format="appended" offset=")" << off_type << R"("/>
      </Cells>
      <CellData Scalars="Pressure">
        <DataArray type=")" << ftype << R"(" Name="Pressure"             NumberOfComponents="1" format="appended" offset=")" << off_pres        << R"("/>
        <DataArray type=")" << ftype << R"(" Name="Pressure_PhiDt"       NumberOfComponents="1" format="appended" offset=")" << off_phi_dt      << R"("/>
        <DataArray type=")" << ftype << R"(" Name="Pressure_Kinetic"     NumberOfComponents="1" format="appended" offset=")" << off_kinetic     << R"("/>
        <DataArray type=")" << ftype << R"(" Name="Pressure_Hydrostatic" NumberOfComponents="1" format="appended" offset=")" << off_hydrostatic << R"("/>
      </CellData>
    </Piece>
  </UnstructuredGrid>
  <AppendedData encoding="raw">_)";

    // ----------------------------------------------------------------
    // Points block: interleaved x,y,z per node
    // ----------------------------------------------------------------
    {
        std::vector<cusfloat> pts_buf;
        pts_buf.reserve( n_nodes * 3 );
        for ( std::size_t i = 0; i < n_nodes; i++ )
        {
            pts_buf.push_back( nodes_x[i] );
            pts_buf.push_back( nodes_y[i] );
            pts_buf.push_back( nodes_z[i] );
        }
        out.write( reinterpret_cast<const char*>( &pts_bytes ), sizeof( uint32_t ) );
        out.write( reinterpret_cast<const char*>( pts_buf.data( ) ), pts_bytes );
    }

    // Connectivity block
    out.write( reinterpret_cast<const char*>( &conn_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( connectivity ), conn_bytes );

    // Offsets block
    out.write( reinterpret_cast<const char*>( &offs_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( offsets ), offs_bytes );

    // Types block
    out.write( reinterpret_cast<const char*>( &type_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( types ), type_bytes );

    // Pressure (cell data) block
    out.write( reinterpret_cast<const char*>( &pres_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( pressure ), pres_bytes );

    // Pressure_PhiDt block
    out.write( reinterpret_cast<const char*>( &pres_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( phi_dt_comp ), pres_bytes );

    // Pressure_Kinetic block
    out.write( reinterpret_cast<const char*>( &pres_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( kinetic_comp ), pres_bytes );

    // Pressure_Hydrostatic block
    out.write( reinterpret_cast<const char*>( &pres_bytes ), sizeof( uint32_t ) );
    out.write( reinterpret_cast<const char*>( hydrostatic_comp ), pres_bytes );

    // ----------------------------------------------------------------
    // Close
    // ----------------------------------------------------------------
    out << "\n  </AppendedData>\n</VTKFile>\n";
    out.close( );
    return true;
}


// -----------------------------------------------------------------------
// write_pvd  — ParaView Data collection file
// -----------------------------------------------------------------------
bool    write_pvd(
                    const std::string&                                  filename,
                    const std::vector<std::pair<double, std::string>>&  timesteps
                )
{
    std::ofstream out( filename );
    if ( !out.is_open( ) )
    {
        Logger logger;
        logger.error( "Could not create PVD file: " + filename + "\n" );
        return false;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";

    for ( const auto& [t, vtu_path] : timesteps )
    {
        out << "    <DataSet timestep=\"" << t
            << "\" group=\"\" part=\"0\" file=\"" << vtu_path << "\"/>\n";
    }

    out << "  </Collection>\n";
    out << "</VTKFile>\n";
    out.close( );
    return true;
}