
#pragma once

// Include local modules
#include "base.hpp"
#include "scalar_expression.hpp"
#include "operators.hpp"


namespace cut // Namespace Custom Tensor
{
    /********************************************/
    /************** MODULE MACROS ***************/
    /********************************************/
    #define BINARY_EXPR_OPERATOR( OPERATOR, OP_TYPE )                                                                       \
    template<typename RR>                                                                                               \
    BinaryExpression<BinaryExpression<L, R, T, Op>, RR, T, OP_TYPE> operator OPERATOR ( const RR& other ) const         \
    {                                                                                                                   \
        return BinaryExpression<BinaryExpression<L, R, T, Op>, RR, T, OP_TYPE>( *this, other );                         \
    }                                                                                                                   \

    /********************************************/
    /********* Binary Expression Class **********/
    /********************************************/
    template<typename L, typename R, typename T, typename Op>
    class BinaryExpression
    {
    private:
        // If L/R are ScalarExpression<T>, store them by value to own the temporary.
        // Otherwise, store a const reference.
        using LIsScalar = typename std::is_same<L, ScalarExpression<T>>;
        using RIsScalar = typename std::is_same<R, ScalarExpression<T>>;
        using LStore    = typename std::conditional<LIsScalar::value, ScalarExpression<T>, const L&>::type;
        using RStore    = typename std::conditional<RIsScalar::value, ScalarExpression<T>, const R&>::type;

        LStore  _lhs;
        RStore  _rhs;

        static LStore make_lstore( const L& l ) { return l; }
        static RStore make_rstore( const R& r ) { return r; }

    public:
        // Define class constructor
        BinaryExpression( const L& l_in, const R& r_in )
            : _lhs( make_lstore( l_in ) )
            , _rhs( make_rstore( r_in ) )
        { }

        // Define class methods
        void    evaluate( T* mstore ) const
        {
            for ( size_t i=0; i<this->size( ); i++ )
            {
                mstore[i] = Op::apply( this->_lhs[i], this->_rhs[i] );
            }
        }

        vector_st shape( void ) const
        {
            vector_st rshape = this->_lhs.shape( );
            if ( 
                    !( 
                        ( this->_lhs.size( ) == 1 ) 
                        || 
                        ( this->_rhs.size( ) == 1 ) 
                    ) 
                )
            {
                if ( this->_lhs.shape( ) != this->_rhs.shape( ) )
                {
                    CUSTENSOR_SHAPE_MISMATCH( this->_lhs.shape( ), this->_rhs.shape( ) )
                }
            }
            else
            {
                if ( this->_rhs.size( ) > this->_lhs.size( ) )
                {
                    rshape = this->_rhs.shape( );
                }
            }

            return rshape;
        }

        size_t  size( void ) const
        {
            return std::max( this->_lhs.size( ), this->_rhs.size( ) );
        }

        // Define class operators
        const T operator[ ]( size_t i ) const
        {
            return Op::apply( this->_lhs[i], this->_rhs[i] );
        }

        // Define aritmetic operators
        BINARY_EXPR_OPERATOR( +, AddOp  )
        BINARY_EXPR_OPERATOR( -, SubOp  )
        BINARY_EXPR_OPERATOR( *, MultOp )
        BINARY_EXPR_OPERATOR( /, DivOp  )

    };
}