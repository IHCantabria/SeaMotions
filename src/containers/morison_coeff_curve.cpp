
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
#include <filesystem>
#include <fstream>

// Include local modules
#include "morison_coeff_curve.hpp"


cusfloat MorisonCoeffCurve::eval( cusfloat x ) const
{
    return (*this->_interp)( x );
}


MorisonCoeffCurve::MorisonCoeffCurve( 
                                        std::string foname,
                                        std::string finame
                                    )
{
    // Read curve data from file
    this->_parse_input_file( foname, finame );

    // Generate interpolation object
    this->_interp = new Interp1D<cusfloat>( 
                                                this->_x_data.size( ),
                                                this->_x_data.data( ), 
                                                this->_y_data.data( ) 
                                            );
    
}


MorisonCoeffCurve::~MorisonCoeffCurve(  )
{
    if ( this->_interp != nullptr )
    {
        delete this->_interp;
        this->_interp = nullptr;
    }
}


void MorisonCoeffCurve::_parse_input_file( 
                                            std::string foname,
                                            std::string finame
                                        )
{
    // Alias namespaces
    namespace fs = std::filesystem;

    // Load curve data
    fs::path curve_foname_( foname );
    fs::path curve_finame_( finame );
    fs::path curve_fipath_ = curve_foname_ / curve_finame_;

    // Open file unit
    std::ifstream file( curve_fipath_.string( ) );

    if ( !file.is_open( ) )
    {
        std::cerr << "Error opening file: " << curve_fipath_ << std::endl;
        throw std::runtime_error( "Error opening file" );
    }

    // Loop over file lines
    std::size_t line_count_ = 0;
    while (  file.good( ) &&
            !file.eof( )
        )
    {
        // Increment line count
        ++line_count_;
    }

    // Remove header line from line count
    line_count_ -= 1;

    // Rewind file unit to read data
    file.clear( );
    file.seekg( 0, std::ios::beg );

    // Allocate data vectors
    cut::CusTensor<cusfloat> x_data_( line_count_ );
    cut::CusTensor<cusfloat> y_data_( line_count_ );

    // Loop over file lines to read data
    std::string dummy_string_ = "";
    file >> dummy_string_; // Skip header line
    file >> dummy_string_; // Skip header line

    for ( std::size_t i = 0; i < line_count_; ++i )
    {
        file >> this->_x_data[i] >> this->_y_data[i];
    }

    // Close file unit
    file.close( );

    // Store data in class attributes
    this->_x_data = std::move( x_data_ );
    this->_y_data = std::move( y_data_ );

}


std::size_t MorisonCoeffCurve::memory_bytes( void ) const
{
    std::size_t total = sizeof( MorisonCoeffCurve );
    total += this->_x_data.size( ) * sizeof( cusfloat );
    total += this->_y_data.size( ) * sizeof( cusfloat );
    if ( this->_interp != nullptr )
    {
        total += this->_interp->memory_bytes( );
    }
    return total;
}


void MorisonCoeffCurve::append_memory_report(
                                                std::vector<MemoryReportEntry>& entries,
                                                const std::string& prefix
                                            ) const
{
    add_memory_entry( entries, memory_report_path( prefix, "object" ), sizeof( MorisonCoeffCurve ) );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "x_data" ),
                        this->_x_data.size( ) * sizeof( cusfloat )
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "y_data" ),
                        this->_y_data.size( ) * sizeof( cusfloat )
                    );
    if ( this->_interp != nullptr )
    {
        add_memory_entry(
                            entries,
                            memory_report_path( prefix, "interp_spline_coeffs" ),
                            this->_interp->memory_bytes( )
                        );
    }
}