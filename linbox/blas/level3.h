/* linbox/blas/level3.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * BLAS Level 3 public interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_H
#define __BLAS_LEVEL3_H

#include "linbox/blas/context.h"
#include "linbox/matrix/matrix-traits.h"
#include "linbox/blas/level3-ll.tcc"

namespace LinBox
{

/** This namespace contains the level 3 BLAS interface */
namespace BLAS3
{

/// @name Operations on matrices
//@{

/** Copy A into B
 *
 * A and B must both be defined over the same field.
 *
 * @param ctx @ref Context object for calculation
 * @param A Origin matrix
 * @param B Destination matrix
 * @returns Reference to B
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &copy (Context<Field, Modules> &ctx, const Matrix1 &A, Matrix2 &B)
	{ return _copy (ctx.F, ctx.M, A, B); }

/** Scale the matrix A by a scalar, A <- a A
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element scalar
 * @param A Matrix, to be replaced by result of computation
 * @returns Reference to A
 */

template <class Field, class Modules, class Matrix>
Matrix &scal (Context<Field, Modules> &ctx, const typename Field::Element &a, Matrix &A)
	{ return _scal (ctx.F, ctx.M, a, A); }

/** Matrix-scalar axpy, B <- a A + B
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element scalar
 * @param A Matrix, to be replaced by result of computation
 * @returns Reference to A
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B)
	{ return _axpy (ctx.F, ctx.M, a, A, B); }

/** General matrix-matrix multiply, C <- a AB + b C
 *
 * Not all combinations of iterator-types are available.
 *
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param B Matrix B
 * @param b Field::Element scalar b
 * @param C Matrix C, to be replaced by result of calculation
 * @returns Reference to y
 */

template <class Field, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C)
	{ return _gemm (ctx.F, ctx.M, a, A, B, b, C); }

/** Triangular matrix-matrix multiply, B <- a AB, where A is triangular
 *
 * A must be square
 *
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param B Matrix B
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &trmm (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return _trmm (ctx.F, ctx.M, a, A, B, type, diagIsOne); }

/** Triangular matrix-matrix solve, B <- a A^-1 B, where A is triangular
 *
 * A must be square. If diagIsOne is false and there are zeros on the
 * diagonal, an exception is thrown.
 *
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param B Matrix B
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &trsm (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return _trsm (ctx.F, ctx.M, a, A, B, type, diagIsOne); }

/** Permute rows of A, A <- PA, where P is a permutation
 *
 * @param P_begin beginning iterator of permutation
 * @param P_end ending iterator of permutation
 * @param A Matrix A, to be replaced by result of computation
 * @returns Reference to A
 */

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_rows (Context<Field, Modules> &ctx, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return _permute_rows (ctx.F, ctx.M, P_begin, P_end, A); }

/** Permute columns of A, A <- AP, where P is a permutation
 *
 * @param P_begin beginning iterator of permutation
 * @param P_end ending iterator of permutation
 * @param A Matrix A, to be replaced by result of computation
 * @returns Reference to A
 */

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_cols (Context<Field, Modules> &ctx, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return _permute_cols (ctx.F, ctx.M, P_begin, P_end, A); }

//@} Operations on matrices

/// @name Queries on matrices
//@{

/** Test whether A and B are equal
 *
 * @param ctx @ref Context object for calculation
 * @param A First matrix
 * @param B Second matrix
 * @returns true if equal, false otherwise
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
bool equal (Context<Field, Modules> &ctx, const Matrix1 &A, const Matrix2 &B)
	{ return _equal (ctx.F, ctx.M, A, B); }

/** Test whether A is the zero matrix
 *
 * @param ctx @ref Context object for calculation
 * @param A Matrix
 * @returns true if A is zero, false otherwise
 */

template <class Field, class Modules, class Matrix>
bool is_zero (Context<Field, Modules> &ctx, const Matrix &A)
	{ return _is_zero (ctx.F, ctx.M, A); }

//@} Queries on matrices

/// @name I/O of matrices
//@{

/** Read the given matrix from the given stream
 *
 * @param ctx Context-object
 * @param is Input-stream
 * @param A Matrix into which to place result
 * @param format Format in which matrix is, default FORMAT_DETECT, which tries to guess
 * @returns Reference to is
 */

template <class Field, class Modules, class Matrix>
std::istream &read (Context<Field, Modules> &ctx, std::istream &is, Matrix &A, FileFormatTag format = FORMAT_DETECT)
	{ return _read (ctx.F, ctx.M, is, A, format); }

/** Write the given matrix to the given stream
 *
 * @param ctx Context-object
 * @param os Output-stream
 * @param A Matrix to be written
 * @param format Format in which to write matrix, default FORMAT_PRETTY
 * @returns Reference to os
 */

template <class Field, class Modules, class Matrix>
std::ostream &write (Context<Field, Modules> &ctx, std::ostream &os, const Matrix &A, FileFormatTag format = FORMAT_PRETTY)
	{ return _write (ctx.F, ctx.M, os, A, format); }

//@} I/O of matrices

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
