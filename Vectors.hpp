/**
 *
 * This file is for educational purposes only,
 * the matrix operations are not optimized in any way,
 * please do not use it in a real production environment.
 *
 */

#pragma once

#include "Matrix.hpp"

template <usize nDims, typename Type = double>
using Vector = Matrix<nDims, 1, Type>;

template <usize nDims, typename Type = double>
Type dot(const Vector<nDims, Type> & lhs, const Vector<nDims, Type> & rhs)
{
    return (lhs.transpose() * rhs)(1);
}

template <usize nDims, typename Type = double>
Matrix<nDims, nDims, Type> cross(const Vector<nDims, Type> & lhs, const Vector<nDims, Type> & rhs)
{
    return lhs * rhs.transpose();
}
