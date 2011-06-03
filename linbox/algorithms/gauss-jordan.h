/* linbox/algorithms/gauss-jordan.h
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Gauss-Jordan elimination
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ALGORITHMS_GAUSS_JORDAN_H
#define __LINBOX_ALGORITHMS_GAUSS_JORDAN_H

#include <iostream>
#include <iomanip>
#include <cassert>

#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"
#include "linbox/blas/context.h"
#include "linbox/field/gf2.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"

namespace LinBox
{

template <class Field>
class Adaptor {
public:
	typedef typename Vector<Field>::Sparse SparseVector;
	typedef LinBox::SparseMatrix<typename Field::Element, SparseVector> SparseMatrix;
	typedef LinBox::DenseMatrix<typename Field::Element> DenseMatrix;
	static const size_t cutoff = 1;
};

template <>
class Adaptor<GF2> {
public:
	typedef BigEndian<uint64> Endianness;
	typedef Vector<GF2>::Hybrid SparseVector;
	typedef LinBox::SparseMatrix<bool, SparseVector> SparseMatrix;
	typedef LinBox::DenseMatrix<bool> DenseMatrix;
	static const size_t cutoff = WordTraits<SparseVector::word_type>::bits;
};

/**
 * \brief Implementation of Gauss-Jordan elimination which
 * does not permute columns in the input-matrix
 *
 * The pivoting-strategy which it uses is to find the row with
 * the least number of nonzero entries and a nonzero entry in
 * the pivot-column.
 *
 * This class depends on the sparse-vector-facilities of LinBox, see http://www.linalg.org/
 *
 * This class depends on the following template-parameters:
 *
 * @param Field The field in which arithmetic takes
 * place. Must satisfy the field-archetype in LinBox.
 */
template <class Field, class Modules = AllModules<Field> >
class GaussJordan
{
public:
	typedef typename Field::Element Element;
	typedef std::pair<uint32, uint32> Transposition;
	typedef std::vector<Transposition> Permutation;
	typedef typename Adaptor<Field>::SparseMatrix SparseMatrix;
	typedef typename Adaptor<Field>::DenseMatrix DenseMatrix;

private:
	Context<Field, Modules> &ctx;

	const size_t _cutoff;

	// Find a suitable pivot from the matrix A starting at
	// start_row. If no pivot can be found (i.e. all rows
	// from start_row onwards are already 0) return
	// -1. Otherwise fill in col with the pivot-column.
	template <class Matrix>
	int GetPivot (const Matrix &A, int start_row, size_t &col) const
		{ return GetPivotSpecialised (A, start_row, col, typename VectorTraits<Field, typename Matrix::Row>::VectorCategory ()); }

	// Find the first nonzero element in the given column
	// starting at the row of the same index. Return -1 if
	// none found.
	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorCategories::DenseVectorTag) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorCategories::DenseZeroOneVectorTag) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorCategories::SparseVectorTag) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorCategories::SparseZeroOneVectorTag) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorCategories::HybridZeroOneVectorTag) const;

	// Set the given matrix to the identity-matrix
	template <class Matrix>
	void SetIdentity (Matrix &U, size_t start_row = 0) const;

	struct CompareSecond
	{
		bool operator () (const Transposition &t1, const Transposition &t2) const { return t1.second < t2.second; }
	};

	// Round n up to the nearest multiple of m
	template <class T>
	T round_up (T n, T m) const
		{ return m * ((n + m - 1) / m); }

	// Internal recursive procedure for the kth indexed
	// Gauss-Jordan transform. Uses a divide and conquer
	// method to maximise use of fast
	// matrix-multiplication. See Chapter 2 of "Algorithms
	// for Matrix Canonical Forms", Ph.D thesis by Arne
	// Storjohann.
	void GaussJordanTransform (DenseMatrix  &A,
				   int           k,
				   Element       d_0,
				   DenseMatrix  &U,
				   Permutation  &P,
				   size_t       &r,
				   int          &h,
				   Element      &d,
				   DenseMatrix  &S,
				   DenseMatrix  &T) const;

	// Internal recursive procedure for the Gauss transform. Uses
	// a divide and conquer method to maximise use of fast
	// matrix-multiplication. See Chapter 2 of "Algorithms for
	// Matrix Canonical Forms", Ph.D thesis by Arne Storjohann.
	void GaussTransform (DenseMatrix             &A,
			     Element                  d_0,
			     typename DenseMatrix::SubmatrixType &U,
			     Permutation             &P,
			     size_t                  &r,
			     int                     &h,
			     Element                 &d) const;

	// Optimised version of VD.addin which can take
	// advantage of knowledge of where in v the entries of
	// w start
	template <class Vector>
	Vector &FastAxpy (Vector &v, const typename Field::Element &a, const Vector &w, size_t idx) const;

	bool testFastAxpyHybridVector () const;

	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
					      VectorCategories::DenseVectorTag) const;

	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
					      VectorCategories::SparseVectorTag) const;

	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
					      VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector>
	class PivotRowCompare
	{
	public:
		inline bool operator () (const typename Vector::index_type x, const Vector &v1) const
			{ return x < v1.front ().first; }
	};

	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
					      VectorCategories::HybridZeroOneVectorTag) const;

	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
					      VectorCategories::DenseZeroOneVectorTag) const;

public:
	/**
	 * \brief Constructor
	 *
	 * @param _F Field over which operations take place
	 */
	GaussJordan (Context<Field, Modules> &_ctx)
		: ctx (_ctx), _cutoff (Adaptor<Field>::cutoff)
		{}

	/**
	 * \brief Convert the matrix A into (possibly reduced)
	 * row-echelon form.
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in reduced
	 * row-echelon form, L is lower triangular, and P is a
	 * permutation.
	 *
	 * If A is invertible, then L will be the inverse of
	 * A.
	 *
	 * @param A Dense matrix A. Will be replaced by its reduced
	 * row-echelon form
	 *
	 * @param U Dense matrix into which to store the matrix
	 * U. Should be n x n, with n the row-dimension of A.
	 *
	 * @param P Permutation in which to store the
	 * row-permutations of A made by the choice of pivots.
	 *
	 * @param rank Integer into which to store the
	 * computed rank
	 *
	 * @param det Field-element into which to store the
	 * computed determinant
	 *
	 * @param reduced Whether the output should be in reduced
	 * row-echelon form or not
	 */
	void DenseRowEchelonForm (DenseMatrix &A,
				  DenseMatrix &U,
				  Permutation &P,
				  size_t      &rank,
				  Element     &det,
				  bool         reduced = false);

	/**
	 * \brief Compute the (reduced or non-reduced)
	 * row-echelon form of a matrix using a sparse
	 * algorithm
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in reduced
	 * row-echelon form, L is lower triangular, and P is a
	 * permutation.
	 *
	 * In comparison with @see RowEchelonForm, this
	 * version does not take advantage of fast
	 * matrix-multiplication and does not use a
	 * divide-and-conquer method. It also modifies the
	 * input-matrix A.
	 *
	 * The pivot-strategy is to find the row with the
	 * fewest elements. This seems to be the only sensible
	 * approach, given that we are not allowed to permute
	 * columns.
	 *
	 * @param A The sparse matrix whose reduced
	 * row-echelon form to compute. Will be replaced by
	 * its reduced row-echelon form during computation.
	 *
	 * @param L The dense matrix into which to store the
	 * matrix L as defined above. Should be of size n x n,
	 * with n the row-dimension of A.
	 *
	 * @param P The permutation into which to store the
	 * permutation P as defined above.
	 *
	 * @param rank An integer into which to store the
	 * computed rank of A.
	 *
	 * @param det A field-element into which to store the
	 * computed determinant of the submatrix of A formed
	 * by taking pivot-rows and -columns.
	 *
	 * @param reduced True if the routine should compute
	 * the reduced row-echelon form and false if it should
	 * only compute the (non-reduced) row-echelon form.
	 *
	 * @param compute_L True if the matrix L should be
	 * computed. If false, then L is ignored.
	 *
	 * @param start_row Start at this row. Intended for
	 * internal use.
	 */
	template <class Matrix1, class Matrix2>
	void StandardRowEchelonForm (Matrix1       &A,
				     Matrix2       &L,
				     Permutation   &P,
				     size_t        &rank,
				     Element       &det,
				     bool           reduced = false,
				     bool           compute_L = true,
				     size_t         start_row = 0) const;

	/** \brief Take a matrix of known rank in row-echelon
	 * form and convert to reduced row-echelon form
	 *
	 * @param A Input matrix A in row-echelon form;
	 * replaced by its reduced row-echelon form
	 *
	 * @param L Dense matrix L into which to store
	 * conversion-matrix
	 *
	 * @param compute_L bool, true if L should be
	 * computed; if false then L is left unchanged
	 *
	 * @param rank Rank of A; must be known a priori
	 * (though can be easily computed by scanning rows of
	 * A)
	 *
	 * @param start_row Pivot-rows begin at this
	 * row. Intended for internal use only.
	 *
	 * @return Reference to A
	 */
	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelon (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row = 0) const
		{ return ReduceRowEchelonSpecialised (A, L, compute_L, rank, start_row,
						      typename VectorTraits<Field, typename Matrix1::Row>::VectorCategory ()); }

	/** Run internal tests
	 */
	void RunTests () const;
};

} // namespace LinBox

#include "linbox/algorithms/gauss-jordan.tcc"

#endif // __LINBOX_ALGORITHMS_GAUSS_JORDAN_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
