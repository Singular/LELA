/* linbox/blas/level3-ll.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_LL_TCC
#define __BLAS_LEVEL3_LL_TCC

#include "linbox/blas/context.h"
#include "linbox/matrix/matrix-traits.h"
#include "linbox/blas/level3-ll.h"

namespace LinBox
{

/** This namespace contains the level 3 BLAS interface */
namespace BLAS3
{

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_copy (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B)
	{ return copy_impl (F, M, A, B,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix>
Matrix &_scal (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A)
	{ return scal_impl (F, M, a, A, typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_axpy (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B)
	{ return axpy_impl (F, M, a, A, B,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
	{ return gemm_impl (F, M, a, A, B, b, C,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix3>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trmm (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return trmm_impl (F, M, a, A, B, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trsm (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return trsm_impl (F, M, a, A, B, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &_permute_rows (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return permute_rows_impl (F, M, P_begin, P_end, A,
				    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &_permute_cols (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return permute_cols_impl (F, M, P_begin, P_end, A,
				     typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix1, class Matrix2>
bool _equal (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B)
	{ return equal_impl (F, M, A, B,
			     typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			     typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix>
bool _is_zero (const Ring &F, Modules &M, const Matrix &A)
	{ return is_zero_impl (F, M, A, typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

template <class Ring, class Modules, class Matrix>
std::istream &_read (const Ring &F, Modules &M, std::istream &is, Matrix &A, FileFormatTag format)
	{ return read_impl (F, M, is, A, format); }

template <class Ring, class Modules, class Matrix>
std::ostream &_write (const Ring &F, Modules &M, std::ostream &os, const Matrix &A, FileFormatTag format)
	{ return write_impl (F, M, os, A, format); }

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_LL_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
