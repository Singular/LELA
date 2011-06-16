/* linbox/blas/level3-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_GENERIC_H
#define __BLAS_LEVEL3_GENERIC_H

#include <algorithm>

#include "linbox/blas/context.h"
#include "linbox/matrix/traits.h"
#include "linbox/matrix/io.h"

namespace LinBox
{

namespace BLAS3
{

template <class Ring, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Ring &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

// FIXME: Not yet implemented (need generic way to attach entry to a vector first)
template <class Ring, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Ring &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag);

// FIXME: Not yet implemented (need generic way to attach entry to a vector first)
template <class Ring, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Ring &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Ring &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Ring, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Ring &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return copy_impl (F, M, A, B, MatrixCategories::RowMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Modules, class Matrix>
Matrix &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixCategories::RowMatrixTag);

template <class Ring, class Modules, class Matrix>
Matrix &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixCategories::ColMatrixTag);

template <class Ring, class Modules, class Matrix>
Matrix &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixCategories::RowColMatrixTag)
	{ return scal_impl (F, M, a, A, MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return axpy_impl (F, M, a, A, B, MatrixCategories::RowMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Ring &F, GenericModule &M,
		    const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Ring &F, GenericModule &M,
		    const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Ring &F, GenericModule &M,
		    const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Ring &F, GenericModule &M,
		    const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Ring &F, GenericModule &M,
		    const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return gemm_impl (F, M, a, A, B, b, C, MatrixCategories::RowMatrixTag (), MatrixCategories::ColMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Ring &F, GenericModule &M,
		    const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
		    MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return gemm_impl (F, M, a, A, B, b, C, MatrixCategories::RowMatrixTag (), MatrixCategories::ColMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::RowMatrixTag);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::ColMatrixTag);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::RowColMatrixTag)
	{ return permute_rows_impl (F, M, P_begin, P_end, A, MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::RowMatrixTag);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::ColMatrixTag);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::RowColMatrixTag)
	{ return permute_cols_impl (F, M, P_begin, P_end, A, MatrixCategories::ColMatrixTag ()); }

template <class Ring, class Matrix1, class Matrix2>
bool equal_impl (const Ring &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Ring, class Matrix1, class Matrix2>
bool equal_impl (const Ring &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Ring, class Matrix1, class Matrix2>
bool equal_impl (const Ring &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return equal_impl (F, M, A, B, MatrixCategories::RowMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Matrix>
bool is_zero_impl (const Ring &F, GenericModule &M, const Matrix &A, MatrixCategories::RowMatrixTag);

template <class Ring, class Matrix>
bool is_zero_impl (const Ring &F, GenericModule &M, const Matrix &A, MatrixCategories::ColMatrixTag);

template <class Ring, class Matrix>
bool is_zero_impl (const Ring &F, GenericModule &M, const Matrix &A, MatrixCategories::RowColMatrixTag)
	{ return is_zero_impl (F, M, A, MatrixCategories::RowMatrixTag ()); }

template <class Ring, class Modules, class Matrix>
std::istream &read_impl (const Ring &F, Modules &M, std::istream &is, Matrix &A, FileFormatTag format)
	{ MatrixReader<Ring> reader (F); return reader.read (is, A, format); }

template <class Ring, class Modules, class Matrix>
std::ostream &write_impl (const Ring &F, Modules &M, std::ostream &os, const Matrix &A, FileFormatTag format)
	{ MatrixWriter<Ring> writer (F); return writer.write (os, A, format); }

} // namespace BLAS3

} // namespace LinBox

#include "linbox/blas/level3-generic.tcc"

#endif // __BLAS_LEVEL3_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
