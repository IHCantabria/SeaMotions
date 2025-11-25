
// Include general usage libraries
#include <fstream>

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
    std::cout << "FilePath: " << filename << std::endl;
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