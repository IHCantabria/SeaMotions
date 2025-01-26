
// Include local modules
#include "math_tools.hpp"
#include "cusvector.hpp"

// Define class constructors and destructor
template<class T>
CusVector<T>::CusVector( int n )
{
    this->_np   = n;
    this->_data = generate_empty_vector<T>( this->_np );
}


template<class T>
CusVector<T>::~CusVector( )
{
    // Deallocate associated heap memory
    mkl_free( this->_data );
}


// Define class operators
template<class T>
T& CusVector<T>::operator []( int idx )
{
    #ifdef DEBUG_BUILD
    CHECK_VECTOR_BOUNDS( idx, this->_np )
    #endif
    return this->_data[ idx ];
}