/* linbox/blas/level3-m4ri.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_M4RI_H
#define __BLAS_LEVEL3_M4RI_H

#ifndef __LINBOX_HAVE_M4RI
#  error "This header file requires that LinBox be configured with libm4ri enabled. Please ensure that libm4ri is properly installed and re-run configure."
#endif

#include <m4ri/m4ri.h>

#include "linbox/blas/context.h"
#include "linbox/blas/level3-ll.h"
#include "linbox/matrix/traits.h"
#include "linbox/matrix/m4ri-matrix.h"
#include "linbox/matrix/submatrix.h"

#ifdef __BLAS_LEVEL3_H
#  warning "linbox/blas/level3.h has already been included by this point. The M4RI-specialisations may not work."
#endif // __BLAS_LEVEL3_H

namespace LinBox
{

namespace BLAS3
{

template <>
class _copy<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const GF2 &F, Modules &M, const Matrix1 &A, Matrix2 &B)
		{ return _copy<GF2, M4RIModule::Tag::Parent>::op (F, M, A, B); }

	template <class Modules>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M, const M4RIMatrixBase &A, M4RIMatrixBase &B)
		{ mzd_copy (B._rep, A._rep); return B; }

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M, const DenseMatrix<bool> &A, DenseMatrix<bool> &B)
		{ op (F, M, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M, const M4RISubmatrix &A, DenseMatrix<bool> &B)
		{ op (F, M, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }

	template <class Modules>
	static M4RISubmatrix &op (const GF2 &F, Modules &M, const M4RISubmatrix &A, M4RISubmatrix &B)
		{ op (F, M, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }
};

template <>
class _scal<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, bool a, Matrix &A)
		{ return _scal<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A); }

	template <class Modules>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M, bool a, M4RIMatrixBase &A);

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M, bool a, DenseMatrix<bool> &A)
		{ op (F, M, a, (M4RIMatrixBase &) A); return A; }
};

template <>
class _axpy<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B)
		{ return _axpy<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A, B); }

	template <class Modules>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M, bool a, const M4RIMatrixBase &A, M4RIMatrixBase &B);

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B)
		{ op (F, M, a, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }
};

template <>
class _gemm<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const GF2 &F, Modules &M, bool a, const Matrix1 &A, const Matrix2 &B, bool b, Matrix3 &C)
		{ return _gemm<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A, B, b, C); }

	template <class Modules>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M,
				   bool a, const M4RIMatrixBase &A, const M4RIMatrix &B, bool b, M4RIMatrixBase &C);

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M,
				      bool a, const DenseMatrix<bool> &A, const DenseMatrix<bool> &B, bool b, DenseMatrix<bool> &C)
		{ op (F, M, a, (const M4RIMatrixBase &) A, (const M4RIMatrix &) B, b, (M4RIMatrixBase &) C); return C; }

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M,
				      bool a, const M4RISubmatrix &A, const DenseMatrix<bool> &B, bool b, DenseMatrix<bool> &C)
		{ op (F, M, a, (const M4RIMatrixBase &) A, (const M4RIMatrix &) B, b, (M4RIMatrixBase &) C); return C; }

	template <class Modules>
	static M4RISubmatrix &op (const GF2 &F, Modules &M,
				  bool a, const M4RISubmatrix &A, const DenseMatrix<bool> &B, bool b, M4RISubmatrix &C)
		{ op (F, M, a, (const M4RIMatrixBase &) A, (const M4RIMatrix &) B, b, (M4RIMatrixBase &) C); return C; }
};

template <>
class _trsm<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trsm<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	template <class Modules>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M, bool a, const M4RIMatrixBase &A, M4RIMatrixBase &B, TriangularMatrixType type, bool diagIsOne);

	template <class Modules>
	static DenseMatrix<bool> &op (const GF2 &F, Modules &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B, TriangularMatrixType type, bool diagIsOne)
		{ op (F, (Modules &) M, a, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B, type, diagIsOne); return B; }
};

template <>
class _permute_rows<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return _permute_rows<GF2, M4RIModule::Tag::Parent>::op (F, M, P_begin, P_end, A); }

	template <class Modules, class Iterator>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, M4RIMatrixBase &A);
};

template <>
class _permute_cols<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return _permute_cols<GF2, M4RIModule::Tag::Parent>::op (F, M, P_begin, P_end, A); }

	template <class Modules, class Iterator>
	static M4RIMatrixBase &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, M4RIMatrixBase &A);
};

template <>
class _equal<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static bool op (const GF2 &F, Modules &M, const Matrix1 &A, const Matrix2 &B)
		{ return _equal<GF2, M4RIModule::Tag::Parent>::op (F, M, A, B); }

	template <class Modules>
	static bool op (const GF2 &F, Modules &M, const M4RIMatrix &A, const M4RIMatrix &B)
		{ return mzd_equal (A._rep, B._rep); }

	template <class Modules>
	static bool op (const GF2 &F, Modules &M, const DenseMatrix<bool> &A, const DenseMatrix<bool> &B)
		{ return op (F, M, (const M4RIMatrix &) A, (const M4RIMatrix &) B); }
};

template <>
class _is_zero<GF2, M4RIModule::Tag>
{
public:
	template <class Modules, class Matrix1>
	static bool op (const GF2 &F, Modules &M, const Matrix1 &A)
		{ return _is_zero<GF2, M4RIModule::Tag::Parent>::op (F, M, A); }

	bool op (const GF2 &F, M4RIModule &M, const M4RIMatrixBase &A)
		{ return mzd_is_zero (A._rep); }

	bool op (const GF2 &F, M4RIModule &M, const DenseMatrix<bool> &A)
		{ return op (F, M, (const M4RIMatrixBase &) A); }

	bool op (const GF2 &F, M4RIModule &M, const M4RISubmatrix &A)
		{ return op (F, M, (const M4RIMatrixBase &) A); }
};

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_M4RI_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
