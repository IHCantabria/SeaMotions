
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
