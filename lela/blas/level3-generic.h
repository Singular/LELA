/* lela/blas/level3-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_GENERIC_H
#define __BLAS_LEVEL3_GENERIC_H

#include <algorithm>

#include "lela/blas/context.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/io.h"
#include "lela/blas/level3-ll.h"

namespace LELA
{

namespace BLAS3
{

template <class Ring>
class _copy<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Vector, class Iterator>
	static void append_entries_spec (const Ring &R, const Vector &v, size_t idx, Iterator begin, Iterator end,
					 VectorRepresentationTypes::Dense);

	template <class Vector, class Iterator>
	static void append_entries_spec (const Ring &R, const Vector &v, size_t idx, Iterator begin, Iterator end,
					 VectorRepresentationTypes::Sparse);

	template <class Vector, class Iterator>
	static void append_entries_spec (const Ring &R, const Vector &v, size_t idx, Iterator begin, Iterator end,
					 VectorRepresentationTypes::Dense01);

	template <class Vector, class Iterator>
	static void append_entries_spec (const Ring &R, const Vector &v, size_t idx, Iterator begin, Iterator end,
					 VectorRepresentationTypes::Sparse01);

	template <class Vector, class Iterator>
	static void append_entries_spec (const Ring &R, const Vector &v, size_t idx, Iterator begin, Iterator end,
					 VectorRepresentationTypes::Hybrid01);

	template <class Vector, class Iterator>
	static inline void append_entries (const Ring &R, const Vector &v, size_t idx, Iterator begin, Iterator end)
		{ append_entries_spec (R, v, idx, begin, end, typename VectorTraits<Ring, Vector>::RepresentationType ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Row, MatrixIteratorTypes::Row);

	// FIXME: Not yet implemented (need generic way to attach entry to a vector first)
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Row, MatrixIteratorTypes::Col);

	// FIXME: Not yet implemented (need generic way to attach entry to a vector first)
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Col, MatrixIteratorTypes::Row);

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Col, MatrixIteratorTypes::Col);

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol)
		{ return copy_impl (F, M, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Row ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Row, MatrixIteratorTypes::RowCol)
		{ return copy_impl (F, M, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Row ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::RowCol, MatrixIteratorTypes::Row)
		{ return copy_impl (F, M, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Row ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Col, MatrixIteratorTypes::RowCol)
		{ return copy_impl (F, M, A, B, MatrixIteratorTypes::Col (), MatrixIteratorTypes::Col ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::RowCol, MatrixIteratorTypes::Col)
		{ return copy_impl (F, M, A, B, MatrixIteratorTypes::Col (), MatrixIteratorTypes::Col ()); }

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B)
		{ return copy_impl (F, M, A, B, typename Matrix1::IteratorType (), typename Matrix2::IteratorType ()); }
};

template <class Ring>
class _scal<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix>
	static Matrix &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixIteratorTypes::Row);

	template <class Modules, class Matrix>
	static Matrix &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixIteratorTypes::Col);

	template <class Modules, class Matrix>
	static Matrix &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixIteratorTypes::RowCol)
		{ return scal_impl (F, M, a, A, MatrixIteratorTypes::Row ()); }

public:
	template <class Modules, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A)
		{ return scal_impl (F, M, a, A, typename Matrix::IteratorType ()); }
};

template <class Ring>
class _axpy<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Row, MatrixIteratorTypes::Row);

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::Col, MatrixIteratorTypes::Col);

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
				   MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol)
		{ return axpy_impl (F, M, a, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Row ()); }

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B)
		{ return axpy_impl (F, M, a, A, B, typename Matrix1::IteratorType (), typename Matrix2::IteratorType ()); }
};

template <class Ring>
class _gemm<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::Row, MatrixIteratorTypes::Row, MatrixIteratorTypes::Row);

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::Col, MatrixIteratorTypes::Col, MatrixIteratorTypes::Col);

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::Row, MatrixIteratorTypes::Col, MatrixIteratorTypes::Generic);

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::Col, MatrixIteratorTypes::Row, MatrixIteratorTypes::Generic);

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol)
		{ return gemm_impl (F, M, a, A, B, b, C, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Col (), MatrixIteratorTypes::Generic ()); }

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C,
				    typename Matrix1::IteratorType (),
				    typename Matrix2::IteratorType (),
				    typename Matrix3::IteratorType ()); }
};

template <class Ring>
class _trmm<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne);

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixIteratorTypes::Row)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   VectorRepresentationTypes::Dense)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   VectorRepresentationTypes::Dense01)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixIteratorTypes::Col)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne, typename VectorTraits<Ring, typename Matrix2::Col>::RepresentationType ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixIteratorTypes::RowCol)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne, typename VectorTraits<Ring, typename Matrix2::Col>::RepresentationType ()); }

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne, typename Matrix2::IteratorType ()); }
};

template <class Ring>
class _trsm<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne);

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixIteratorTypes::Row)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   VectorRepresentationTypes::Dense)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   VectorRepresentationTypes::Dense01)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixIteratorTypes::Col)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne, typename VectorTraits<Ring, typename Matrix2::Col>::RepresentationType ()); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixIteratorTypes::RowCol)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne, typename VectorTraits<Ring, typename Matrix2::Col>::RepresentationType ()); }

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne, typename Matrix2::IteratorType ()); }
};

template <class Ring>
class _permute_rows<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Row);

	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Col);

	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::RowCol)
		{ return permute_rows_impl (F, M, P_begin, P_end, A, MatrixIteratorTypes::Row ()); }

public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return permute_rows_impl (F, M, P_begin, P_end, A, typename Matrix::IteratorType ()); }
};

template <class Ring>
class _permute_cols<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Row);

	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Col);

	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::RowCol)
		{ return permute_cols_impl (F, M, P_begin, P_end, A, MatrixIteratorTypes::Col ()); }

public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return permute_cols_impl (F, M, P_begin, P_end, A, typename Matrix::IteratorType ()); }
};

template <class Ring>
class _equal<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static bool equal_impl (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
				MatrixIteratorTypes::Generic, MatrixIteratorTypes::Generic);

	template <class Modules, class Matrix1, class Matrix2>
	static bool equal_impl (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
				MatrixIteratorTypes::Row, MatrixIteratorTypes::Row);

	template <class Modules, class Matrix1, class Matrix2>
	static bool equal_impl (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
				MatrixIteratorTypes::Col, MatrixIteratorTypes::Col);

	template <class Modules, class Matrix1, class Matrix2>
	static bool equal_impl (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
				MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol)
		{ return equal_impl (F, M, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Row ()); }

public:
	template <class Modules, class Matrix1, class Matrix2>
	static bool op (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B)
		{ return equal_impl (F, M, A, B, typename Matrix1::IteratorType (), typename Matrix2::IteratorType ()); }
};

template <class Ring>
class _is_zero<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Matrix>
	static bool is_zero_impl (const Ring &F, Modules &M, const Matrix &A, MatrixIteratorTypes::Row);

	template <class Modules, class Matrix>
	static bool is_zero_impl (const Ring &F, Modules &M, const Matrix &A, MatrixIteratorTypes::Col);

	template <class Modules, class Matrix>
	static bool is_zero_impl (const Ring &F, Modules &M, const Matrix &A, MatrixIteratorTypes::RowCol)
		{ return is_zero_impl (F, M, A, MatrixIteratorTypes::Row ()); }

public:
	template <class Modules, class Matrix>
	static bool op (const Ring &F, Modules &M, const Matrix &A)
		{ return is_zero_impl (F, M, A, typename Matrix::IteratorType ()); }
};

template <class Ring>
class _read<Ring, typename GenericModule<Ring>::Tag>
{
public:
	template <class Modules, class Matrix>
	static std::istream &op (const Ring &F, Modules &M, std::istream &is, Matrix &A, FileFormatTag format = FORMAT_DETECT)
		{ MatrixReader<Ring> reader (F); return reader.read (is, A, format); }
};

template <class Ring>
class _write<Ring, typename GenericModule<Ring>::Tag>
{
public:
	template <class Modules, class Matrix>
	static std::ostream &op (const Ring &F, Modules &M, std::ostream &os, const Matrix &A, FileFormatTag format)
		{ MatrixWriter<Ring> writer (F); return writer.write (os, A, format); }
};

} // namespace BLAS3

} // namespace LELA

#include "lela/blas/level3-generic.tcc"

#include "lela/matrix/io.tcc"

#endif // __BLAS_LEVEL3_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
