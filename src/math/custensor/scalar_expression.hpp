
#pragma once

// Include local modules
#include "../../config.hpp"


namespace cut // Namespace Custom Tensor
{
    template<typename T>
    class ScalarExpression
    {
    private:
        T   _value;

    public:
        // Define Constructor
        ScalarExpression( T value_in ): _value( value_in ) { }

        // Define class methods
        vector_st shape( void ) const
        {
            return { 1 };
        }

        size_t  size( void ) const
        { 
            return 1;
        }

        // Define class operators
        const T operator[ ] ( size_t  ) const 
        { 
            return this->_value; 
        }

    };
}