/* lela/algorithms/pivot-strategy.h
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

#ifndef __LELA_ALGORITHMS_PIVOT_STRATEGY_H
#define __LELA_ALGORITHMS_PIVOT_STRATEGY_H

namespace LELA
{

/** Pivot-strategy
 *
 * A pivot-strategy decides which entry in a matrix to use as a pivot
 * in elimination-algorithms. This class encapsulates such a
 * strategy. It is provided as an example for documentation only.
 *
 * \ingroup algorithms
 */
class PivotStrategy 
{
public:
	/** Obtain a pivot from the matrix
	 *
	 * The parameters row and col should be preset to the starting
	 * row and column, respectively. They are then set to the
	 * location of the pivot in the matrix, guaranteed to be at
	 * least the starting-point. The entry at that position is
	 * stored in pivot.
	 *
	 * If no pivot was found (because e.g. the matrix below
	 * (row,col) is zero) then the function returns
	 * false. Otherwise it returns true.
	 *
	 * @param A Matrix from which to find pivot
	 * @param pivot Place to store pivot-entry (not deep-copied)
	 *
	 * @param row Row: at input equal to starting row, at
	 * conclusion equal to row where pivot was found
	 *
	 * @param col Column: at input equal to starting column, at
	 * conclusion equal to column where pivot was found
	 *
	 * @return true if a pivot was found, false otherwise
 	 */
	template <class Matrix>
	bool getPivot (const Matrix &A, typename Matrix::Element &pivot, size_t &row, size_t &col) const;
};

/** Dense pivot-strategy
 *
 * This pivot-strategy is appropriate for dense matrices. It just
 * scans the matrix column by column for the first nonzero entry.
 *
 * This strategy is only available for matrices with dense or dense
 * 0-1 rows.
 *
 * \ingroup algorithms
 */
template <class Ring, class Modules>
class DensePivotStrategy 
{
	Context<Ring, Modules> &ctx;

	template <class Matrix>
	bool getPivot_spec (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col,
			    VectorRepresentationTypes::Dense) const;

	template <class Matrix>
	bool getPivot_spec (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col,
			    VectorRepresentationTypes::Dense01) const;

public:
	/** Constructor
	 *
	 * _ctx Context-object to be used
	 */
	DensePivotStrategy (Context<Ring, Modules> &_ctx)
		: ctx (_ctx) 
	{}

	template <class Matrix>
	bool getPivot (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col) const
		{ return getPivot_spec (A, pivot, row, col, typename VectorTraits<Ring, typename Matrix::Row>::RepresentationType ()); }
};

/** Sparse pivot-strategy without column-permutations
 *
 * This pivot-strategy is useful for sparse matrices in
 * elimination-problems where column-permutation is not permitted. It
 * searches among rows with the leftmost nonzero entry for the row
 * with the fewest nonzero entries.
 *
 * This strategy is only available for matrices with sparse, sparse
 * 0-1, or hybrid 0-1 rows.
 *
 * \ingroup algorithms
 */
template <class Ring, class Modules>
class SparsePartialPivotStrategy 
{
	Context<Ring, Modules> &ctx;

	template <class Matrix>
	bool getPivot_spec (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col,
			    VectorRepresentationTypes::Sparse) const;

	template <class Matrix>
	bool getPivot_spec (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col,
			    VectorRepresentationTypes::Sparse01) const;

	template <class Matrix>
	bool getPivot_spec (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col,
			    VectorRepresentationTypes::Hybrid01) const;

public:
	/** Constructor
	 *
	 * _ctx Context-object to be used
	 */
	SparsePartialPivotStrategy (Context<Ring, Modules> &_ctx)
		: ctx (_ctx) 
	{}

	template <class Matrix>
	bool getPivot (const Matrix &A, typename Ring::Element &pivot, size_t &row, size_t &col) const
		{ return getPivot_spec (A, pivot, row, col, typename VectorTraits<Ring, typename Matrix::Row>::RepresentationType ()); }
};

/** Sensible default pivot-strategies for row-types */
template <class Ring, class Modules, class Row, class Trait = typename VectorTraits<Ring, Row>::RepresentationType>
struct DefaultPivotStrategy
	{ typedef PivotStrategy Strategy; };

template <class Ring, class Modules, class Row>
struct DefaultPivotStrategy<Ring, Modules, Row, VectorRepresentationTypes::Dense>
	{ typedef DensePivotStrategy<Ring, Modules> Strategy; };

template <class Ring, class Modules, class Row>
struct DefaultPivotStrategy<Ring, Modules, Row, VectorRepresentationTypes::Dense01>
	{ typedef DensePivotStrategy<Ring, Modules> Strategy; };

template <class Ring, class Modules, class Row>
struct DefaultPivotStrategy<Ring, Modules, Row, VectorRepresentationTypes::Sparse>
	{ typedef SparsePartialPivotStrategy<Ring, Modules> Strategy; };

template <class Ring, class Modules, class Row>
struct DefaultPivotStrategy<Ring, Modules, Row, VectorRepresentationTypes::Sparse01>
	{ typedef SparsePartialPivotStrategy<Ring, Modules> Strategy; };

template <class Ring, class Modules, class Row>
struct DefaultPivotStrategy<Ring, Modules, Row, VectorRepresentationTypes::Hybrid01>
	{ typedef SparsePartialPivotStrategy<Ring, Modules> Strategy; };

} // namespace LELA

#include "lela/algorithms/pivot-strategy.tcc"

#endif // __LELA_ALGORITHMS_PIVOT_STRATEGY_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
