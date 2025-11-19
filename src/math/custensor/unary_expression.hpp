
#pragma once

// Include local modules
#include "binary_expression.hpp"


namespace cut // Namespace Custom Tensor
{
    // Define a macro to generate both overloads for a given operator and operation type.
    #define UNARY_EXPR_OPERATOR( OPERATOR, OP_TYPE )                                                                        \
    template<typename RR, typename = typename std::enable_if<!std::is_arithmetic<RR>::value>::type>                     \
    BinaryExpression<UnaryExpression<R, T, Op>, RR, T, OP_TYPE> operator OPERATOR (const RR& other)                     \
    {                                                                                                                   \
        return BinaryExpression<UnaryExpression<R, T, Op>, RR, T, OP_TYPE>(*this, other);                               \
    }                                                                                                                   \
                                                                                                                        \
    template<typename RR, typename = typename std::enable_if<std::is_arithmetic<RR>::value>::type>                      \
    BinaryExpression<UnaryExpression<R, T, Op>, ScalarExpression<T>, T, OP_TYPE> operator OPERATOR (const RR& other)    \
    {                                                                                                                   \
        return BinaryExpression<UnaryExpression<R, T, Op>, ScalarExpression<T>, T, OP_TYPE>(                            \
            *this, ScalarExpression<T>(other));                                                                         \
    }                                                                                                                 \

    /********************************************/
    /********** Unary Expression Class **********/
    /********************************************/
    template<typename R, typename T, typename Op>
    class UnaryExpression
    {
    private:
        using RIsScalar = typename std::is_same<R, ScalarExpression<T>>;
        using RStore    = typename std::conditional<RIsScalar::value, ScalarExpression<T>, const R&>::type;

        RStore  _expr;

        static RStore make_rstore( const R& r ) { return r; }

    public:
        // Define constructors and destructor
        UnaryExpression( const R& r ): _expr( make_rstore( r ) ) { }

        // Define class methods
        void    evaluate( T* mstore ) const
        {
            for ( size_t i=0; i<this->_expr.size( ); i++ )
            {
                mstore[i] = Op::apply( this->_expr[i] );
            }
        }

        vector_st shape( void ) const
        {
            return this->_expr.shape( );
        }

        size_t  size( ) const
        {
            return this->_expr.size( );
        }

        // Define operators
        const T operator [ ] ( size_t i ) const
        {
            return Op::apply( this->_expr[i] );
        }

        // Define unary operatiors
        UNARY_EXPR_OPERATOR( +, AddOp  )
        UNARY_EXPR_OPERATOR( -, SubOp  )
        UNARY_EXPR_OPERATOR( *, MultOp )
        UNARY_EXPR_OPERATOR( /, DivOp  )

    };
}