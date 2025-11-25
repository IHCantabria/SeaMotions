
// Include local modules
#include "../config.hpp"
#include "logger.hpp"


const char* Logger::_color( Level level ) const 
{
    switch ( level ) 
    {
        case Level::INFO:    return TermColor::WHITE;
        case Level::WARNING: return TermColor::YELLOW;
        case Level::ERROR:   return TermColor::RED;
        case Level::TASK:    return TermColor::GREEN;
    }

    return TermColor::WHITE;
}


void Logger::error( const std::string &msg )
{ 
    this->log( Level::ERROR, msg );
}


void Logger::info( const std::string &msg )
{ 
    this->log( Level::INFO, msg );
}


void Logger::log( Level level, const std::string& msg ) 
{
    if ( this->_is_root )
    {
        std::cout   << this->_color( level ) << this->_prefix( level )
                    << msg << TermColor::RESET << std::endl << std::flush;
    }
}


Logger::Logger( void )
{
    // Get current process by using the MPI system
    // Get current process rank
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Check if it is the root process
    this->_is_root = proc_rank == MPI_ROOT_PROC_ID;
}


Logger::Logger( bool is_root )
{
    // Storage input arguments
    this->_is_root = is_root;
}


Logger::Logger( bool is_root, int time_colum )
{
    // Storage input arguments
    this->_is_root      = is_root;
    this->_time_column  = time_colum;
}


std::string Logger::_prefix( Level level ) const 
{
    switch ( level ) 
    {
        case Level::INFO:    return "[INFO]: ";
        case Level::WARNING: return "[WARNING]: ";
        case Level::ERROR:   return "[ERROR]: ";
        case Level::TASK:    return "[TASK]: ";
    }

    return "";
}


void Logger::set_time_column( int col ) 
{ 
    this->_time_column = col; 
}


void Logger::task( const std::string &msg ) 
{
    if ( this->_is_root )
    {
        this->_last_msg = msg;
        std::cout << TermColor::GREEN << this->_prefix( Level::TASK );
        std::cout << msg << TermColor::RESET << std::flush;
    }
}


void Logger::task_time( const std::string& msg, MpiTimer& timer )
{
    // Get new message as last message to append the time
    this->_last_msg = msg;

    // Task time
    this->task_time( timer );
}


void Logger::task_time( MpiTimer& timer )
{
    // Get time elapsed. Take care with the timer as 
    // it request syncrhonization of all processes,
    // it couldn't be only in the root process call.
    std::ostringstream time_out;
    time_out << " - Elapsed Time: " << timer;
    std::string time_str = time_out.str( );

    // Reject
    if ( this->_is_root )
    {
        // Alignment padding
        int padding = this->_time_column - static_cast<int>( this->_last_msg.size( ) );
        if ( padding < 1 ) padding = 1;
    
        std::cout << std::string( padding, ' ' );
        std::cout << TermColor::GREEN << " --> Done!" << TermColor::RESET;
        std::cout << TermColor::CYAN << time_str << TermColor::RESET;
        std::cout << "\n" << std::flush;
    }

}


void Logger::warn( const std::string &msg )
{ 
    this->log( Level::WARNING, msg );
}
