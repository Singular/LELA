/* linbox/matrix/matrix-domain.h
 * Copyright (C) 2002 Zhendong Wan, Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------------------------------
 * 2002-11-26  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Added detailed documentation, cleaned up the interface slightly, and added
 * support for matrix traits. Added read, write, neg, negin, axpy, and
 * matrix-vector and matrix-black box operations.
 * ------------------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_domain_H
#define __LINBOX_matrix_domain_H

#include <iostream>

#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-traits.h"

// For some specialisations
#include "linbox/matrix/transpose-matrix.h"

namespace LinBox
{

/** \brief Helper class to allow specializations of certain matrix-vector products
 *
 * This class implements a method gemvColDense that multiplies a
 * column-represented matrix by a dense vector
 */

// FIXME: Now that we have MatrixDomainSupport, do we need a separate
// MVProductDomain? What does that bring to us?

template <class Field>
class MVProductDomain
{
    public:
	typedef typename Field::Element Element;

	MVProductDomain () {}

    protected:
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColDense (const VectorDomain<Field> &VD, const typename Field::Element &alpha, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y) const;
};

/** \brief Class of basic matrix-arithmetic-operations
 *
 * This class contains definitions for a BLAS-like interface for a set
 * of matrix-matrix and matrix-vector operations. It is intended that
 * it be specialised to various fields. The class \ref MatrixDomain
 * below makes use of these operations to provide an interface at a
 * higher level.
 *
 * The intention behind this class is to make it easier to specialize
 * matrix-arithmetic to various fields without having to reimplement
 * functionality which is common to all fields.
 *
 * Developers using LinBox should declare an instance of \ref
 * MatrixDomain instead of this class.
 */
template <class Field>
class MatrixDomainSupportGeneric : public MVProductDomain<Field>
{
protected:
	typename Field::Element _zero;
	typename Field::Element _one;
	typename Field::Element _neg_one;

	const Field &_F;
	const VectorDomain<Field> _VD;

	///
	MatrixDomainSupportGeneric (const Field &F) : _F (F), _VD (F) {
		F.init (_zero, 0);
		F.init (_one, 1);
		F.init (_neg_one, -1);
	}

public:

	/*? @name Matrix-vector arithmetic operations
	 *
	 * These operations take a matrix satisfying the \ref
	 * DenseMatrix archetype and LinBox vectors as inputs. They
	 * involve matrix-vector product and matrix-vector AXPY. They
	 * are equivalent to BLAS level 2.
	 */

	/** General matrix-vector multiplication
	 * y <- alpha A x + b y
	 *
	 * @param alpha Input scalar alpha
	 * @param A Input matrix A
	 * @param x Input vector x
	 * @param b Input scalar b
	 * @param y Output vector y
	 * @returns Reference to y
	 */
	template <class Vector1, class Matrix, class Vector2>
	inline Vector2 &gemv (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y) const
		{ return gemvSpecialized (a, A, x, b, y, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Triangular solve with vector
	 * x <- A^{-1} x
	 *
	 * @param A Input matrix A, assumed to be upper triangular
	 * @param x Input vector x
	 * @returns Reference to x
	 */
	template <class Matrix, class Vector>
	inline Vector &trsv (const Matrix &A, Vector &x) const
		{ return trsvSpecialized (A, x, typename MatrixTraits<Matrix>::MatrixCategory (), typename VectorTraits<Vector>::VectorCategory ()); }

	/*? @name Matrix-matrix arithmetic operations
	 *
	 * These provide the equivalent of BLAS level 3 functionality
	 */

	/** Matrix copy
	 * B <- A
	 * Copy the contents of the matrix B to the matrix A
	 *
	 * Both matrices must support the same iterators, row or column.
	 * 
	 * @param B Matrix B
	 * @param A Matrix A
	 * @returns Reference to B
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix1 &copy (Matrix1 &B, const Matrix2 &A) const
		{ return copySpecialized (B, A,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix equality
	 * Test whether the matrices A and B are equal
	 * @param A Input vector
	 * @param B Input vector
	 * @returns true if and only if the matrices A and B are equal
	 */
	template <class Matrix1, class Matrix2>
	inline bool areEqual (const Matrix1 &A, const Matrix2 &B) const
		{ return areEqualSpecialized (B, A,
					      typename MatrixTraits<Matrix1>::MatrixCategory (),
					      typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix equality with zero
	 * @param A Input matrix
	 * @returns true if and only if the matrix A is zero
	 */
	template <class Matrix>
	inline bool isZero (const Matrix &A) const
		{ return isZeroSpecialized (A, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Matrix-scalar multiply
	 * A <- a * A
	 *
	 * Multiply A by the scalar element a.
	 *
	 * @param C Output matrix C
	 * @param B Input matrix B
	 * @param a Input scalar a
	 * @returns Reference to C
	 */
	template <class Matrix>
	inline Matrix &scal (Matrix &A, const typename Field::Element &a) const
		{ return scalSpecialized (A, a,
					  typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Matrix-axpy
	 * B <- a * A + B
	 *
	 * A and B must support the same (row- or column-) iterators
	 * 
	 * @param a Input scalar a
	 * @param A Input matrix A
	 * @param B Output matrix B
	 * @returns Reference to A
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpy (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B) const
		{ return axpySpecialized (a, A, B,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** General matrix-matrix multiply
	 * C <- a * A * B + b * C
	 *
	 * C must support both row and column iterators, and the vector
	 * representations must be dense. Examples of supported matrices are
	 * \ref DenseMatrixBase and \ref DenseSubmatrix.
	 *
	 * Either A or B, or both, may have limited iterators. However, either A
	 * must support row iterators or B must support column iterators. If
	 * both A and B lack support for an iterator (either row or column),
	 * then C must support the same type of iterator as A and B.
	 *
	 * @param a Input scalar a
	 * @param A Input matrix A
	 * @param B Input matrix B
	 * @param b Input scalar b
	 * @param C Output matrix C
	 * @returns Reference to C
	 */
	template <class Matrix1, class Matrix2, class Matrix3>
	inline Matrix3 &gemm (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C) const
		{ return gemmSpecialized (a, A, B, b, C,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory (),
					  typename MatrixTraits<Matrix3>::MatrixCategory ()); }

	/** Triangular solve with matrix
	 * B <- a A^{-1} B
	 *
	 * @param a Input scalar a
	 * @param A Input matrix A, assumed to be upper triangular
	 * @param B Input matrix B
	 * @returns Reference to B
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix2 &trsm (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B) const
		{ return trsmSpecialized (a, A, B, typename MatrixTraits<Matrix1>::MatrixCategory (), typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/*? @name Matrix permutations
	 * These operations permute the rows or columns of a matrix based on
	 * the given permutation. They are intended for use with Gauss-Jordan
	 * elimination
	 */

	/** Permutation
	 *
	 * A permutation is represented as a vector of pairs, each
	 * pair representing a transposition.
	 */
	typedef std::pair<unsigned int, unsigned int> Transposition;
	typedef std::vector<Transposition> Permutation;

	/** Permute the rows of the given matrix
	 *
	 * @param A Output matrix
	 * @param P_start Start of permutation
	 * @param P_end End of permutation
	 * @returns Reference to A
	 */
	template <class Matrix, class Iterator>
	inline Matrix &permuteRows (Matrix   &A,
				    Iterator  P_start,
				    Iterator  P_end) const
		{ return permuteRowsSpecialized (A, P_start, P_end,
						 typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Permute the columns of the given matrix
	 *
	 * @param A Output matrix
	 * @param P_start Start of permutation
	 * @param P_end End of permutation
	 * @returns Reference to A
	 */
	template <class Matrix, class Iterator>
	inline Matrix &permuteColumns (Matrix   &A,
				       Iterator  P_start,
				       Iterator  P_end) const
		{ return permuteColsSpecialized (A, P_start, P_end,
						 typename MatrixTraits<Matrix>::MatrixCategory ()); }

private:
	// Specialized function implementations
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseSequenceVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseAssociativeVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvRowSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseParallelVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::DenseVectorTag,
				     VectorCategories::DenseVectorTag) const
		{ return gemvColDense (_VD, a, A, x, b, y); } 
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseSequenceVectorTag,
				     VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseAssociativeVectorTag,
				     VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseParallelVectorTag,
				     VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvColSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				     VectorCategories::SparseParallelVectorTag,
				     VectorCategories::SparseParallelVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				  MatrixCategories::RowMatrixTag) const
		{ return gemvRowSpecialized (a, A, x, b, y, typename VectorTraits<Vector2>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				 MatrixCategories::ColMatrixTag) const
		{ return gemvColSpecialized (a, A, x, b, y,
					    typename VectorTraits<Vector1>::VectorCategory (),
					    typename VectorTraits<Vector2>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector2 &gemvSpecialized (const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
				 MatrixCategories::RowColMatrixTag) const
		{ return gemvRowSpecialized (a, A, x, b, y, typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Matrix, class Vector>
	Vector &trsvSpecialized (const Matrix &A, Vector &x, MatrixCategories::RowMatrixTag, VectorCategories::DenseVectorTag);

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
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag) const
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
				  MatrixCategories::RowColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
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
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::RowColMatrixTag) const
		{ return isZeroRow (A); }

	template <class Matrix> Matrix &scalRow (Matrix &A, const typename Field::Element &a) const;
	template <class Matrix> Matrix &scalCol (Matrix &A, const typename Field::Element &a) const;

	template <class Matrix>
	Matrix &scalSpecialized (Matrix &A, const typename Field::Element &a,
				 MatrixCategories::RowMatrixTag) const
		{ return scalRow (A, a); }
	template <class Matrix>
	Matrix &scalSpecialized (Matrix &A, const typename Field::Element &a,
				  MatrixCategories::ColMatrixTag) const
		{ return scalCol (A, a); }
	template <class Matrix>
	Matrix &scalSpecialized (Matrix &A, const typename Field::Element &a,
				  MatrixCategories::RowColMatrixTag) const
		{ return scalRow (A, a); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpyRow (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B) const;
	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpyCol (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpySpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{ return axpyRow (a, A, B); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpySpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{ return axpyCol (a, A, B); }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpySpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{ return axpyRow (a, A, B); }

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowCol (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmColRowCol (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmRowRowRow (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmColColCol (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C) const;

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::RowMatrixTag) const
		{ return gemmRowRowCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return gemmColRowCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
		{ return gemmRowRowCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag) const
		{ return gemmRowRowRow (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return gemmColColCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::RowColMatrixTag,
				  MatrixCategories::RowColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
		{ return gemmRowRowCol (a, A, B, b, C); }

	// These shouldn't be necessary, but for some reason gcc
	// doesn't want to instantiate the above methods when
	// TransposeMatrix is used. Ugh.
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const TransposeMatrix<Matrix1> &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::RowColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
		{ return gemmColColCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const Matrix1 &A, const TransposeMatrix<Matrix2> &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::RowColMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
		{ return gemmColColCol (a, A, B, b, C); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemmSpecialized (const typename Field::Element &a, const TransposeMatrix<Matrix1> &A, const TransposeMatrix<Matrix2> &B, const typename Field::Element &b, Matrix3 &C,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
		{ return gemmColColCol (a, A, B, b, C); }

	template <class Matrix1, class Matrix2>
	Matrix2 &trsmSpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

	template <class Matrix1, class Matrix2>
	Matrix2 &trsmSpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag);

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
					       MatrixCategories::RowColMatrixTag) const
		{ return permuteRowsByCol (A, P_start, P_end); }

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
					       MatrixCategories::RowColMatrixTag) const
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

/** Class which enables a plug-in-architecture for matrix-operations
 *
 * May be specialised for specific fields and according to the
 * presence of optional libraries.
 */
template <class Field>
class MatrixDomainSupport : public MatrixDomainSupportGeneric<Field>
{
    public:
	MatrixDomainSupport (const Field &F) : MatrixDomainSupportGeneric<Field> (F) {}
};

/** \brief Class of matrix arithmetic functions
 *
 * This class encapuslated matrix-matrix and matrix-vector operations, roughly
 * equivalent to BLAS levels 2 and 3. The arithmetic methods are parameterized
 * by matrix type so that they may be used the same way with sparse matrices,
 * dense matrices, and dense submatrices. Except where otherwise noted, they
 * require the matrix inputs to meet the \ref DenseMatrix archetype.
 *
 * These methods are specialized so that they can run efficiently with different
 * matrix representations. If a matrix has an efficient row iterator, but not an
 * efficient column iterator, a specialization that makes use of the former will
 * be selected. This allows a great deal of flexibility when dealing with sparse
 * matrix arithmetic.
 *
 * For all of the arithmetic operations that output matrices, it is assumed that
 * the output matrix has an efficient row iterator. In typical use, the output
 * matrix will be a \ref DenseMatrixBase or a \ref DenseSubmatrix, which has
 * efficient row and column iterators. In particular, one should not perform
 * these arithmetic operations outputting to a \ref SparseMatrixBase.
 *
 * There are other restrictions. See the method-specific documentation for more
 * details.
 */

template <class Field>
class MatrixDomain : public MatrixDomainSupport<Field>
{
    public:
	typedef std::pair<unsigned int, unsigned int> Transposition;
	typedef std::vector<Transposition> Permutation;

	///
	MatrixDomain (const Field &F) : MatrixDomainSupport<Field> (F) {}

	/** Retrieve the underlying field
	 * Return a reference to the field that this matrix domain
	 * object uses
	 * @returns reference to field
	 */
	const Field &field () const
		{ return MatrixDomainSupport<Field>::_F; }
	Field &field () 
		{ return MatrixDomainSupport<Field>::_F; }

	/*? @name Field-independent functionality
	 *
	 * These functions do not do any arithmetic and therefore are
	 * not designed to be specialised to different fields.
	 */

	/** Print matrix.
	 * @param  os  Output stream to which matrix is written.
	 * @param  A   Matrix.
	 * @returns reference to os.
	 */
	template <class Matrix>
	inline std::ostream &write (std::ostream &os, const Matrix &A) const
		{ return A.write (os, MatrixDomainSupport<Field>::_F); }

	/** Read matrix
	 * @param  is  Input stream from which matrix is read.
	 * @param  A   Matrix.
	 * @returns reference to is.
	 */
	template <class Matrix>
	inline std::istream &read (std::istream &is, Matrix &A) const
		{ return A.read (is, MatrixDomainSupport<Field>::_F); }
};

} // namespace LinBox

#include "linbox/matrix/matrix-domain.inl"

#endif // __LINBOX_matrix_domain_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
