/* lela/matrix/sparse.h
 * Copyright 2001-2002, 2011 Bradford Hovinen
 *           1999-2001 William J Turner,
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * 
 * --------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_SPARSE_H
#define __LELA_MATRIX_SPARSE_H

#ifndef _SP_BB_VECTOR_
#  include <vector>
#  define _SP_BB_VECTOR_ std::vector
#endif

#include <utility>
#include <iostream>
#include <algorithm>
#include <vector>

#include "lela/lela-config.h"
#include "lela/util/debug.h"
#include "lela/vector/traits.h"
#include "lela/matrix/traits.h"
#include "lela/vector/sparse-subvector.h"

namespace LELA
{

/** Sparse matrix
 *
 * See @ref MatrixArchetype for documentation on the interface
 *
 * @param Element Element type
 * @param Row     Vector type to use for rows of matrix
\ingroup matrix
 */
template <class _Element, class _Row = typename RawVector<_Element>::Sparse, class Trait = typename ElementVectorTraits<_Element, _Row>::RepresentationType>
class SparseMatrix
{
    public:

	/// @name @ref MatrixArchetype interface
	//@{

	typedef _Element Element;
	typedef _Row Row;
	typedef SparseMatrix<Element, Row, Trait> Self_t;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;
	typedef MatrixIteratorTypes::Row IteratorType; 
	typedef MatrixStorageTypes::Rows StorageType;

	typedef Submatrix<Self_t> SubmatrixType;
	typedef Submatrix<const Self_t> ConstSubmatrixType;
	typedef Submatrix<Self_t> AlignedSubmatrixType;
	typedef Submatrix<const Self_t> ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = 1;

	typedef Self_t ContainerType;

        SparseMatrix () : _m (0), _n (0) {}
        SparseMatrix (size_t m, size_t n) : _A (m), _m (m), _n (n) {}
	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	SparseMatrix (const SparseMatrix<Element, Row, Trait> &A);

    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType, Trait> &A);

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void setEntry (size_t i, size_t j, const Element &value);
	void eraseEntry (size_t i, size_t j);
	bool getEntry (Element &x, size_t i, size_t j) const;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;

	typedef MatrixRawIterator<RowIterator, VectorRepresentationTypes::Sparse> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, VectorRepresentationTypes::Sparse> ConstRawIterator;
    
	RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, Trait, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	Row      &operator [] (size_t i);
	ConstRow &operator [] (size_t i) const;

	template <class Vector>
	Vector &columnDensity (Vector &v) const;

	//@}

	/// @name Additional interfaces
	//@{

	/** Retreive number of elements in the matrix.
	 * @return integer number of elements of SparseMatrix matrix.
	 */
	size_t size () const
	{
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s += LELA::RawVector<_Element>::size (*it);

		return s;
        }

	/** Construct the transpose of this matrix and place it in the
	 * matrix given
	 *
	 * Unlike @ref TransposeMatrix, which just provides a view
	 * into an existing matrix, this function actually constructs
	 * a copy.
	 */
	SparseMatrix &transpose (SparseMatrix &AT) const;

	//@}

    protected:
	
	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrix;
};

/* Specialization for sparse sequence vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorRepresentationTypes::Sparse>
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse> Self_t;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;
	typedef MatrixIteratorTypes::Row IteratorType; 
	typedef MatrixStorageTypes::Rows StorageType;

	typedef Submatrix<Self_t> SubmatrixType;
	typedef Submatrix<const Self_t> ConstSubmatrixType;
	typedef Submatrix<Self_t> AlignedSubmatrixType;
	typedef Submatrix<const Self_t> ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = 1;

	typedef Self_t ContainerType;

        SparseMatrix () : _m (0), _n (0) {}

	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}

	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	SparseMatrix (const SparseMatrix &A)
		: _A (A._A), _m (A._m), _n (A._n) {}

    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType, VectorRepresentationTypes::Sparse> &A)
		: _A (A._m), _m (A._m), _n (A._n)
	{
		typename Rep::iterator meit = this->_A.begin ();
		typename SparseMatrix<Element, VectorType, VectorRepresentationTypes::Sparse>::Rep::const_iterator copit = A._A.begin();
		for (; meit != this->_A.end (); ++meit, ++copit)
			LELA::RawVector<Element>::convert (*meit, *copit);
        }

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const
	{ 
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s+= LELA::RawVector<_Element>::size (*it);

		return s;
        }

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void setEntry (size_t i, size_t j, const Element &value);
	void eraseEntry (size_t i, size_t j);
	bool getEntry (Element &x, size_t i, size_t j) const;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()       { return _A.begin (); }
	ConstRowIterator rowBegin () const { return _A.begin (); }
	RowIterator      rowEnd ()         { return _A.end (); }
	ConstRowIterator rowEnd () const   { return _A.end (); }

	typedef MatrixRawIterator<ConstRowIterator, VectorRepresentationTypes::Sparse> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorRepresentationTypes::Sparse, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	Row &getRow (size_t i) { return _A[i]; }
	Row &operator [] (size_t i) { return _A[i]; }
	ConstRow &operator [] (size_t i) const { return _A[i]; }

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrix &transpose (SparseMatrix &AT) const;

    protected:

	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrix;
};

} // namespace LELA

#include "lela/matrix/sparse.tcc"

#include "lela/matrix/sparse-zero-one.h"

#endif // __LELA_MATRIX_SPARSE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
