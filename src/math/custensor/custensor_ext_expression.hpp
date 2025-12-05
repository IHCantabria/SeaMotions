
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

// Include local modules
#include "scalar_expression.hpp"
#include "custensor_impl.hpp"
#include "unary_expression.hpp"

/********************************************/
/************** MODULE MACROS ***************/
/********************************************/

/* UNARY OPERATORS MACRO DEFINITION */
#define UNARY_EXT_OPERATOR( FCN_NAME, OP_TYPE )                                                                                 \
template<typename T>                                                                                                            \
cut::UnaryExpression<cut::ScalarExpression<T>, T, OP_TYPE> FCN_NAME( T scalar )                                                 \
{                                                                                                                               \
    return cut::UnaryExpression<cut::ScalarExpression<T>, T, OP_TYPE>( cut::ScalarExpression( scalar ) );                       \
}                                                                                                                               \
                                                                                                                                \
template<typename T>                                                                                                            \
cut::UnaryExpression<cut::CusTensor<T>, T, OP_TYPE> FCN_NAME( const cut::CusTensor<T>& lhs )                                         \
{                                                                                                                               \
    return cut::UnaryExpression<cut::CusTensor<T>, T, OP_TYPE>( lhs );                                                          \
}                                                                                                                               \
                                                                                                                                \
template<typename L, typename R, typename T, typename Op>                                                                       \
cut::UnaryExpression<cut::BinaryExpression<L, R, T, Op>, T, OP_TYPE> FCN_NAME( const cut::BinaryExpression<L, R, T, Op>& lhs )  \
{                                                                                                                               \
    return cut::UnaryExpression<cut::BinaryExpression<L, R, T, Op>, T, OP_TYPE>( lhs );                                         \
}                                                                                                                               \


/* BINARY OPERATORS MACRO DEFINITION */
#define BINARY_EXT_OPERATOR( OPERATOR, OP_TYPE )                                                                                \
template <typename R, typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>                       \
cut::BinaryExpression<cut::ScalarExpression<T>, R, T, OP_TYPE> operator OPERATOR (T scalar, const R& other)                     \
{                                                                                                                               \
    return cut::BinaryExpression<cut::ScalarExpression<T>, R, T, OP_TYPE>( cut::ScalarExpression( scalar ), other );            \
}                                                                                                                               \


/********************************************/
/*********** Define operators  **************/
/********************************************/
/* Define Unary operators */
UNARY_EXT_OPERATOR( nda_exp, cut::ExpOp )


/* Define binary operators */
BINARY_EXT_OPERATOR( +, cut::AddOp  )
BINARY_EXT_OPERATOR( -, cut::SubOp  )
BINARY_EXT_OPERATOR( *, cut::MultOp )
BINARY_EXT_OPERATOR( /, cut::DivOp  )