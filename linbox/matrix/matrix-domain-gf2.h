/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/matrix-domain-gf2.h
 * Copyright 2010 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __MATRIX_MATRIX_DOMAIN_GF2_H
#define __MATRIX_MATRIX_DOMAIN_GF2_H

#include <iostream>
#include <time.h>

#include "linbox/field/gf2.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/dense-zero-one.h"
#include "linbox/matrix/submatrix.h"
#include "linbox/vector/sparse-subvector-hybrid.h"

namespace LinBox 
{

/** Matrix-vector and matrix-matrix arithmetic over GF2
 *
 * This class provides broad support for all combinations of matrix-
 * and vector-types, but uses only naive algorithms and is therefore
 * not tuned well for performance, especially for dense matrices.
 */
class MatrixDomainSupportGF2
{
protected:
	const GF2 &_F;
	const VectorDomain<GF2> &_VD;

	const bool _zero;
	const bool _one;
	const bool _neg_one;

	///
	MatrixDomainSupportGF2 (const GF2 &F) : _F (F), _VD (F), _zero (false), _one (true), _neg_one (true) {}

public:

	template <class Vector1, class Matrix, class Vector2>
	inline Vector2 &gemv (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y) const
		{ return gemvSpecialized (a, A, x, b, y, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	template <class Matrix, class Vector>
	inline Vector &trsv (const Matrix &A, Vector &x) const
		{ return trsvSpecialized (A, x, typename MatrixTraits<Matrix>::MatrixCategory (), typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Matrix1, class Matrix2>
	inline Matrix1 &copy (Matrix1 &B, const Matrix2 &A) const
		{ return copySpecialized (B, A,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	template <class Matrix1, class Matrix2>
	inline bool areEqual (const Matrix1 &A, const Matrix2 &B) const
		{ return areEqualSpecialized (B, A,
					      typename MatrixTraits<Matrix1>::MatrixCategory (),
					      typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	template <class Matrix>
	inline bool isZero (const Matrix &A) const
		{ return isZeroSpecialized (A, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	template <class Matrix>
	inline Matrix &scal (Matrix &A, const bool &a) const
		{ return scalSpecialized (A, a,
					  typename MatrixTraits<Matrix>::MatrixCategory ()); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpy (const bool &a, const Matrix1 &A, Matrix2 &B) const
		{ return axpySpecialized (a, A, B,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	template <class Matrix1, class Matrix2, class Matrix3>
	inline Matrix3 &gemm (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C) const
		{ return gemmSpecialized (a, A, B, b, C,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory (),
					  typename MatrixTraits<Matrix3>::MatrixCategory ()); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &trsm (const bool &a, const Matrix1 &A, Matrix2 &B) const
		{ return trsmSpecialized (a, A, B, typename MatrixTraits<Matrix1>::MatrixCategory (), typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	template <class Matrix1, class Blackbox, class Matrix2>
	inline Matrix2 &gebmm (Matrix2 &C, const Blackbox &A, const Matrix1 &B) const;

	template <class Matrix1, class Matrix2, class Blackbox>
	inline Matrix2 &gembm (Matrix2 &C, const Matrix2 &A, const Blackbox &B) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteRows (Matrix   &A,
				    Iterator  P_start,
				    Iterator  P_end) const
		{ return permuteRowsSpecialized (A, P_start, P_end,
						 typename MatrixTraits<Matrix>::MatrixCategory ()); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColumns (Matrix   &A,
				       Iterator  P_start,
				       Iterator  P_end) const
		{ return permuteColsSpecialized (A, P_start, P_end,
						 typename MatrixTraits<Matrix>::MatrixCategory ()); }

private:
	// Private scratch-pad. It must never contain state, and even
	// const-methods are free to modify it.
	DenseZeroOneMatrix<> _Cp;

	// Specialized function implementations
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColDense (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::HybridZeroOneVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::DenseZeroOneVectorTag,
				     VectorCategories::DenseZeroOneVectorTag) const
		{ return gemvColDense (a, A, x, b, y); } 
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::DenseZeroOneVectorTag,
				     VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::SparseZeroOneVectorTag,
				     VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				     VectorCategories::HybridZeroOneVectorTag,
				     VectorCategories::DenseZeroOneVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				  MatrixCategories::RowMatrixTag) const
		{ return gemvRowSpecialized (a, A, x, b, y, typename GF2VectorTraits<Vector1>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				  MatrixCategories::ColMatrixTag) const
		{ return gemvColSpecialized (a, A, x, b, y,
					    typename GF2VectorTraits<Vector1>::VectorCategory (),
					    typename GF2VectorTraits<Vector2>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
				  MatrixCategories::ZeroOneRowMatrixTag) const
		{ return gemvRowSpecialized (a, A, x, b, y, typename GF2VectorTraits<Vector1>::VectorCategory ()); }

	template <class Matrix, class Vector>
	Vector &trsvSpecialized (const Matrix &A, Vector &x, MatrixCategories::RowMatrixTag, VectorCategories::DenseZeroOneVectorTag) const;

	template <class Matrix1, class Matrix2> Matrix1 &copyRow (Matrix1 &B, const Matrix2 &A) const;
	template <class Matrix1, class Matrix2> Matrix1 &copyCol (Matrix1 &B, const Matrix2 &A) const;

	template <class Matrix1, class Matrix2>
	inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{ return copyRow (B, A); }
	template <class Matrix1, class Matrix2>
	inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{ return copyCol (B, A); }
	template <class Matrix1, class Matrix2>
	inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
					 MatrixCategories::ZeroOneRowMatrixTag,
					 MatrixCategories::ZeroOneRowMatrixTag) const
		{ return copyRow (B, A); }

	template <class Matrix1, class Matrix2> bool areEqualRow (const Matrix1 &A, const Matrix2 &B) const;
	template <class Matrix1, class Matrix2> bool areEqualCol (const Matrix1 &A, const Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{ return areEqualRow (A, B); }
	template <class Matrix1, class Matrix2>
	inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{ return areEqualCol (A, B); }
	template <class Matrix1, class Matrix2>
	inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
					 MatrixCategories::ZeroOneRowMatrixTag,
					 MatrixCategories::ZeroOneRowMatrixTag) const
		{ return areEqualRow (A, B); }

	template <class Matrix> bool isZeroRow (const Matrix &v) const;
	template <class Matrix> bool isZeroCol (const Matrix &v) const;

	template <class Matrix>
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::RowMatrixTag) const
		{ return isZeroRow (A); }
	template <class Matrix>
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::ColMatrixTag) const
		{ return isZeroCol (A); }
	template <class Matrix>
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::ZeroOneRowMatrixTag) const
		{ return isZeroRow (A); }

	template <class Matrix> Matrix &scalRow (Matrix &A, const bool &a) const;
	template <class Matrix> Matrix &scalCol (Matrix &A, const bool &a) const;

	template <class Matrix>
	Matrix &scalSpecialized (Matrix &A, const bool &a,
				 MatrixCategories::RowMatrixTag) const
		{ return scalRow (A, a); }
	template <class Matrix>
	Matrix &scalSpecialized (Matrix &A, const bool &a,
				 MatrixCategories::ColMatrixTag) const
		{ return scalCol (A, a); }
	template <class Matrix>
	Matrix &scalSpecialized (Matrix &A, const bool &a,
				 MatrixCategories::ZeroOneRowMatrixTag) const
		{ return scalRow (A, a); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpyRow (const bool &a, const Matrix1 &A, Matrix2 &B) const;
	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpyCol (const bool &a, const Matrix1 &A, Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpySpecialized (const bool &a, const Matrix1 &A, Matrix2 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{ return axpyRow (a, A, B); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpySpecialized (const bool &a, const Matrix1 &A, Matrix2 &B,
					 MatrixCategories::ZeroOneRowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{ return axpyRow (a, A, B); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpySpecialized (const bool &a, const Matrix1 &A, Matrix2 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{ return axpyCol (a, A, B); }

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowCol (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmColRowCol (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowRow (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C) const
		{ return gemmRowRowRowSpecialised (a, A, B, b, C, typename GF2VectorTraits<typename Matrix1::ConstRow>::VectorCategory ()); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmColColCol (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C) const;

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
					   VectorCategories::DenseZeroOneVectorTag) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
					   VectorCategories::SparseZeroOneVectorTag) const;
	template <class Matrix1, class Matrix3>
	Matrix3 &gemmRowRowRow (const bool &a, const Matrix1 &A, const Submatrix<DenseZeroOneMatrix<> > &B, const bool &b, Matrix3 &C) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
					   VectorCategories::HybridZeroOneVectorTag) const;
	template <class Matrix1, class Matrix3>
	Matrix3 &gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Submatrix<DenseZeroOneMatrix<> > &B, const bool &b, Matrix3 &C,
					   VectorCategories::HybridZeroOneVectorTag) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
					   VectorCategories::HybridZeroOneSequenceVectorTag) const;

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return gemmRowRowCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return gemmColRowCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
				  MatrixCategories::ZeroOneRowMatrixTag,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return gemmRowRowCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag) const
		{ return gemmRowRowRow (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return gemmColColCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
				  MatrixCategories::ZeroOneRowMatrixTag,
				  MatrixCategories::ZeroOneRowMatrixTag,
				  MatrixCategories::ZeroOneRowMatrixTag) const
		{ return gemmRowRowRow (a, A, B, b, C); }

	template <class Matrix1, class Matrix2>
	Matrix2 &trsmSpecialized (const bool &a, const Matrix1 &A, Matrix2 &B, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag) const;

	template <class Matrix1, class Matrix2>
	Matrix2 &trsmSpecialized (const bool &a, const Matrix1 &A, Matrix2 &B, MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsByRow (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsByCol (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::ZeroOneRowMatrixTag) const
		{ return permuteRowsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowMatrixTag) const
		{ return permuteRowsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowColMatrixTag) const
		{ return permuteRowsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::ColMatrixTag) const
		{ return permuteRowsByCol (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsByRow (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsByCol (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::ZeroOneRowMatrixTag) const
		{ return permuteColsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowMatrixTag) const
		{ return permuteColsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::ColMatrixTag) const
		{ return permuteColsByCol (A, P_start, P_end); }
};

#ifdef __LINBOX_HAVE_M4RI
} // namespace LinBox

#  include "linbox/matrix/matrix-domain-m4ri.h"

namespace LinBox {

template <>
class MatrixDomainSupport<GF2> : public MatrixDomainM4RI
{
    public:
	MatrixDomainSupport (const GF2 &F) : MatrixDomainM4RI (F) {}
};

template <>
class MatrixDomainSupport<const GF2> : public MatrixDomainM4RI
{
    public:
	MatrixDomainSupport (const GF2 &F) : MatrixDomainM4RI (F) {}
};
#else // !__LINBOX_HAVE_M4RI
template <>
class MatrixDomainSupport<GF2> : public MatrixDomainSupportGF2
{
    public:
	MatrixDomainSupport (const GF2 &F) : MatrixDomainSupportGF2 (F) {}
};

template <>
class MatrixDomainSupport<const GF2> : public MatrixDomainSupportGF2
{
    public:
	MatrixDomainSupport (const GF2 &F) : MatrixDomainSupportGF2 (F) {}
};
#endif // __LINBOX_HAVE_M4RI

} // namespace LinBox

#include "linbox/matrix/matrix-domain-gf2.inl"

#endif // __MATRIX_MATRIX_DOMAIN_GF2_H
