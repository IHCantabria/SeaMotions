
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
#include <iomanip>
#include <mpi.h>
#include <ostream>


class MpiTimer 
{
public:
    // sync: whether to call MPI_Barrier on stop
    // sync_on_print: whether operator<< should call stop() internally
    MpiTimer( MPI_Comm comm = MPI_COMM_WORLD, bool sync = true, bool sync_on_print = true )
        : _comm( comm ), _sync( sync ), _sync_on_print( sync_on_print )
    {
        if ( _sync ) MPI_Barrier( _comm );

        this->_start    = MPI_Wtime();
        this->_stopped  = false;
        this->_elapsed  = 0.0;
    }

    // manually stop the timer (performs barrier if _sync==true)
    double stop( void )
    {
        if ( !_stopped )
        {
            if ( _sync ) MPI_Barrier( _comm );
            this->_elapsed = MPI_Wtime( ) - _start;
            this->_stopped = true;
        }
        return this->_elapsed;
    }

    // elapsed without stopping
    double elapsed( void) const 
    {
        if ( this->_stopped ) return this->_elapsed;
        return MPI_Wtime( ) - _start;
    }

private:
    MPI_Comm    _comm;
    bool        _sync;
    bool        _sync_on_print;
    double      _start;
    bool        _stopped;
    double      _elapsed;

    friend std::ostream& operator<<(std::ostream& os, MpiTimer& t);
};

// operator<< optionally performs stop() before printing
inline std::ostream& operator<<( std::ostream& os, MpiTimer& t ) 
{
    if ( t._sync_on_print ) 
    {
        t.stop( );  // barrier + elapsed calculation
    }

    os << std::setw( 15 ) << std::fixed << std::setprecision( 3 ) << t.elapsed( ) << " s";

    return os;
}
