/* linbox/algorithms/gauss-jordan.h
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Asymptotically fast Gauss-Jordan elimination
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ALGORITHMS_GAUSS_JORDAN_H
#define __LINBOX_ALGORITHMS_GAUSS_JORDAN_H

#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"
#include "linbox/blas/context.h"
#include "linbox/ring/gf2.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"
#include "linbox/algorithms/elimination.h"

namespace LinBox
{

template <class Ring> struct DefaultCutoff { static const size_t cutoff = 1; };
template <> struct DefaultCutoff<GF2> { static const size_t cutoff = WordTraits<Vector<GF2>::Hybrid::word_type>::bits; };

/** Implementation of asymptotically fast Gauss-Jordan elimination
 */
template <class Ring, class Modules = AllModules<Ring> >
class GaussJordan
{
public:
	typedef typename Ring::Element Element;
	typedef std::pair<uint32, uint32> Transposition;
	typedef std::vector<Transposition> Permutation;

private:
	Context<Ring, Modules> &ctx;
	Elimination<Ring, Modules> elim;

	const size_t _cutoff;

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
	template <class Matrix1, class Matrix2>
	void GaussJordanTransform (Matrix1      &A,
				   int           k,
				   Element       d_0,
				   Matrix2      &U,
				   Permutation  &P,
				   size_t       &r,
				   int          &h,
				   Element      &d,
				   DenseMatrix<typename Ring::Element> &S,
				   DenseMatrix<typename Ring::Element> &T) const;

	// Internal recursive procedure for the Gauss transform. Uses
	// a divide and conquer method to maximise use of fast
	// matrix-multiplication. See Chapter 2 of "Algorithms for
	// Matrix Canonical Forms", Ph.D thesis by Arne Storjohann.
	template <class Matrix1, class Matrix2>
	void GaussTransform (Matrix1     &A,
			     Element      d_0,
			     Matrix2     &U,
			     Permutation &P,
			     size_t      &r,
			     int         &h,
			     Element     &d) const;

public:
	/**
	 * \brief Constructor
	 *
	 * @param _F Ring over which operations take place
	 */
	GaussJordan (Context<Ring, Modules> &_ctx)
		: ctx (_ctx), elim (_ctx), _cutoff (DefaultCutoff<Ring>::cutoff)
		{}

	/**
	 * \brief Convert the matrix A into (possibly reduced)
	 * row-echelon form.
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=UPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in reduced
	 * row-echelon form, U is lower triangular, and P is a
	 * permutation.
	 *
	 * If A is invertible, then U will be the inverse of
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
	 * @param det Ring-element into which to store the
	 * computed determinant
	 *
	 * @param reduced Whether the output should be in reduced
	 * row-echelon form or not
	 */
	template <class Matrix1, class Matrix2>
	void RowEchelonForm (Matrix1     &A,
			     Matrix2     &U,
			     Permutation &P,
			     size_t      &rank,
			     Element     &det,
			     bool         reduced = false);
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
