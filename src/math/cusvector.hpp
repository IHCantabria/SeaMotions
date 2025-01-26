
#ifndef __cusvector_hpp
#define __cusvector_hpp


// Define macros
#define CHECK_VECTOR_BOUNDS( x, n )                                                                 \
{                                                                                                   \
    if ( x < 0 )                                                                                    \
    {                                                                                               \
        std::cerr << "Requested index is less than 0" << std::endl;                                 \
        std::cerr << "  -> Requested Idx: " << x << std::endl;                                      \
        throw std::exception( );                                                                    \
    }                                                                                               \
                                                                                                    \
    if ( x > n )                                                                                    \
    {                                                                                               \
        std::cerr << "Requested index is higher than vector length." << std::endl;                  \
        std::cerr << "  -> Requested Idx: " << x << " - Upper limit: " << n << std::endl;           \
        throw std::exception( );                                                                    \
    }                                                                                               \
}                                                                                                   \


// Declare class structure
template<class T>
class CusVector
{
private:
    T*      _data   = nullptr;
    int     _np     = 0;

public:
    // Declare constructors and destructor
    CusVector( int n );

    ~CusVector( );

    // Declare operators
    T& operator [] ( int idx );
};

// Include template definitions
#include "cusvector.txx"

#endif