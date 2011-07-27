/* lela/algorithms/pivot-strategy.tcc
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Pivot-strategies for elimination-methods
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LELA_ALGORITHMS_PIVOT_STRATEGY_TCC
#define __LELA_ALGORITHMS_PIVOT_STRATEGY_TCC

namespace LELA
{

template <class Ring, class Modules>
template <class Matrix>
bool DensePivotStrategy<Ring, Modules>::getPivot_spec (const Matrix &A, typename Ring::Element &x, size_t &row, size_t &col,
						       VectorRepresentationTypes::Dense) const
{
	lela_check (row < A.rowdim ());
	lela_check (col < A.coldim ());

	typename Matrix::ConstRowIterator i;

	size_t k;

	for (; col < A.coldim (); ++col) {
		for (i = A.rowBegin () + row, k = row; i != A.rowEnd (); ++i, ++k) {
			if (!ctx.F.isZero ((*i)[col])) {
				row = k;
				x = (*i)[col];
				return true;
			}
		}
	}

	return false;
}

template <class Ring, class Modules>
template <class Matrix>
bool DensePivotStrategy<Ring, Modules>::getPivot_spec (const Matrix &A, typename Ring::Element &x, size_t &row, size_t &col,
						      VectorRepresentationTypes::Dense01) const
{
	lela_check (row < A.rowdim ());
	lela_check (col < A.coldim ());

	typename Matrix::ConstRowIterator i;

	size_t k;

	for (; col < A.coldim (); ++col) {
		for (i = A.rowBegin () + row, k = row; i != A.rowEnd (); ++i, ++k) {
			if ((*i)[col]) {
				row = k;
				x = true;
				return true;
			}
		}
	}

	return false;
}

template <class Ring, class Modules>
template <class Matrix>
bool SparsePartialPivotStrategy<Ring, Modules>::getPivot_spec (const Matrix &A, typename Ring::Element &x, size_t &row, size_t &col,
							       VectorRepresentationTypes::Sparse) const
{
	lela_check (row < A.rowdim ());
	lela_check (col < A.coldim ());

	typename Matrix::ConstRowIterator i;

	size_t min_nonzero = 0xffffffffU, k, start_col = col;
	col = A.coldim ();

	for (i = A.rowBegin () + row, k = row; i != A.rowEnd (); ++i, ++k) {
		typename VectorTraits<Ring, typename Matrix::ConstRow>::ConstSubvectorType row_sub (*i, start_col, A.coldim ());

		if (!row_sub.empty ()) {
			if (row_sub.front ().first + start_col < col) {
				col = row_sub.front ().first + start_col;
				min_nonzero = row_sub.size ();
				x = row_sub.front ().second;
				row = k;
			}
			else if (row_sub.front ().first + start_col == col && row_sub.size () < min_nonzero) {
				min_nonzero = row_sub.size ();
				x = row_sub.front ().second;
				row = k;
			}
		}
	}

	return min_nonzero != 0xffffffffU;
}

template <class Ring, class Modules>
template <class Matrix>
bool SparsePartialPivotStrategy<Ring, Modules>::getPivot_spec (const Matrix &A, typename Ring::Element &x, size_t &row, size_t &col,
							       VectorRepresentationTypes::Sparse01) const
{
	lela_check (row < A.rowdim ());
	lela_check (col < A.coldim ());

	typename Matrix::ConstRowIterator i;

	size_t min_nonzero = 0xffffffffU, k, start_col = col;
	col = A.coldim ();

	for (i = A.rowBegin () + row, k = row; i != A.rowEnd (); ++i, ++k) {
		typename VectorTraits<Ring, typename Matrix::ConstRow>::ConstSubvectorType row_sub (*i, start_col, A.coldim ());

		if (!row_sub.empty ()) {
			if (row_sub.front () + start_col < col) {
				col = row_sub.front () + start_col;
				min_nonzero = row_sub.size ();
				x = true;
				row = k;
			}
			else if (row_sub.front () + start_col == col && row_sub.size () < min_nonzero) {
				min_nonzero = row_sub.size ();
				x = true;
				row = k;
			}
		}
	}

	return min_nonzero != 0xffffffffU;
}

template <class Ring, class Modules>
template <class Matrix>
bool SparsePartialPivotStrategy<Ring, Modules>::getPivot_spec (const Matrix &A, typename Ring::Element &x, size_t &row, size_t &col,
							       VectorRepresentationTypes::Hybrid01) const
{
	lela_check (row < A.rowdim ());
	lela_check (col < A.coldim ());

	typename Matrix::ConstRowIterator i;

	size_t min_blocks = 0xffffffffU, k, s, start_col = col;
	col = A.coldim ();

	for (i = A.rowBegin () + row, k = row; i != A.rowEnd (); ++i, ++k) {
		typename VectorTraits<Ring, typename Matrix::ConstRow>::ConstSubvectorType row_sub (*i, start_col, A.coldim ());

		if (!row_sub.empty () && row_sub.front ().first << WordTraits<typename Matrix::Row::word_type>::logof_size <= (int) col) {
			int idx = BLAS1::head (ctx, x, row_sub) + start_col;
			if ((size_t) idx < col) {
				col = idx;
				min_blocks = i->size ();
				row = k;
			}
			else if ((size_t) idx == col && (s = i->size ()) < min_blocks) {
				min_blocks = s;
				row = k;
			}
		}
	}

	return min_blocks != 0xffffffffU;
}

} // namespace LELA

#endif // __LELA_ALGORITHMS_PIVOT_STRATEGY_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
