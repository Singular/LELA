/* linbox/matrix/sparse.h
 * Copyright 2001-2002 Bradford Hovinen
 *           1999-2001 William J Turner,
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * 
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.h to matrix/sparse.h
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
 *   - Named template argument _Element rather than Row; add a typedef to Element
 *   - Renamed RawIndexIterator as RawIndexedIterator, and adjusted to match
 *     interface in DenseMatrixBase
 * ------------------------------------
 * 2002-08-06  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed to sparse-base.h from sparse0-base.h
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrixBase parameterized only on the element type
 *   - New read/write implementations for SparseMatrixBase, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_sparse_H
#define __LINBOX_matrix_sparse_H

#ifndef _SP_BB_VECTOR_
#  include <vector>
#  define _SP_BB_VECTOR_ std::vector
#endif

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/matrix/matrix-traits.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/util/matrix-stream.h"

namespace LinBox
{
	
/** Sparse matrix container
 * This class acts as a generic row-wise container for sparse
 * matrices. It is designed to provide various methods to access the
 * entries of the matrix. In particular, it does not meet the black box
 * archetype; see \ref{SparseMatrix} for an appropriate sparse matrix
 * black box.
 *
 * @param Element Element type
 * @param Row     LinBox sparse vector type to use for rows of matrix
\ingroup matrix
 */
template <class _Element, class _Row = typename RawVector<_Element>::Sparse, class Trait = typename VectorTraits<_Row>::VectorCategory>
class SparseMatrix
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef typename _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrix<typename _Tp1::Element, _R1, Trait> other; };

	/** Constructor.
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
        SparseMatrix (size_t m, size_t n): _A(m), _m(m), _n(n) {};

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrix ( MatrixStream<Field>& ms );

	/** Constructor from a VectorStream
	 *
	 * Fills the row-vectors of the matrix with vectors from the
	 * stream. Stream must have finite size.
	 */
	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	/** Copy constructor.
	 */
	SparseMatrix (const SparseMatrix<Element, Row, Trait> &A);

	/** Convert constructor.
	 */
    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType, Trait> &A);

	/** Destructor. */
	~SparseMatrix () {}

	/** Retreive row dimension of the matrix.
	 * @return integer number of rows of SparseMatrix matrix.
	 */
	size_t rowdim () const { return _m; }

	/** Retreive column dimension of matrix.
	 * @return integer number of columns of SparseMatrix matrix.
	 */
	size_t coldim () const { return _n; }

	/** Retreive number of elements in the matrix.
	 * @return integer number of elements of SparseMatrix matrix.
	 */
	size_t size () const
	{
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s += LinBox::RawVector<_Element>::size (*it);

		return s;
        }

	/** Resize the matrix.
	 *
	 * Does not touch the existing row-vectors, so if the new
	 * column-dimension is less than the existing one, then the
	 * row-vectors may no longer be valid.
	 *
	 * @param m New row-dimension
	 * @param n New column-dimension
	 */
	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	/** Set an individual entry
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @value Value of the new entry
	 */
	void setEntry (size_t i, size_t j, const Element &value);

	/** Erase an individual entry from the matrix.
	 * If the entry doesn't exist, then takes no action.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 */
	void eraseEntry (size_t i, size_t j);

	/** Get a writeable reference to an entry in the matrix
	 * If there is no entry at the position (i, j), then a new entry
	 * with a value of zero is inserted and a reference  to it is
	 * returned.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j);

	/** Get a read-only individual entry from the matrix
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const;

	/** Get an entry and store it in the given value
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const;

	/** @name Columns of rows iterator
	 * The columns of row iterator gives each of the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in sparse sequence format
	 */

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;

	/** @name Raw iterator
	 * The raw iterator is a method for accessing all nonzero
	 * entries in the matrix in some unspecified order. This can be
	 * used, e.g. to reduce all matrix entries modulo a prime before
	 * passing the matrix into an algorithm.
	 */

	typedef MatrixRawIterator<RowIterator, VectorCategories::SparseSequenceVectorTag> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::SparseSequenceVectorTag> ConstRawIterator;
    
	RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd ()); }
	RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd ()); }
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd ()); }

	/** @name Index iterator
	 * The index iterator gives the row, column indices of all matrix
	 * elements in the same order as the raw iterator above. Its value type
	 * is an STL pair with the row and column indices, starting at 0, in the
	 * first and second positions, respectively.
	 */

	typedef MatrixRawIndexedIterator<ConstRowIterator, Trait, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd ()); }

	/** Retrieve a row as a writeable reference
	 * @param i Row index
	 */
	Row &getRow (size_t i);

	/** Retrieve a row as a writeable reference
	 * @param i Row index
	 */
	Row &operator [] (size_t i);

	/** Retrieve a row as a read-only reference
	 * @param i Row index
	 */
	ConstRow &operator [] (size_t i) const;

	/** Compute the column density, i.e. the number of entries per column
	 * @param v Vector in which to store column density 
	 */
	template <class Vector>
	Vector &columnDensity (Vector &v) const;

	/** Construct the transpose of this matrix and place it in the
	 * matrix given
	 */
	SparseMatrix &transpose (SparseMatrix &AT) const;

    protected:
	
	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrix;
};

/* Specialization for sparse sequence vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::SparseSequenceVectorTag >
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template <typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other>
        struct rebind
	{
		typedef SparseMatrix<typename _Tp1::Element, _R1, VectorCategories::SparseSequenceVectorTag> other;
	};

	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrix ( MatrixStream<Field>& ms );

	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	SparseMatrix (const SparseMatrix &A)
		: _A (A._A), _m (A._m), _n (A._n) {}

    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType, VectorCategories::SparseSequenceVectorTag> &A)
		: _A (A._m), _m (A._m), _n (A._n)
	{
		typename Rep::iterator meit = this->_A.begin ();
		typename SparseMatrix<Element, VectorType, VectorCategories::SparseSequenceVectorTag>::Rep::const_iterator copit = A._A.begin();
		for (; meit != this->_A.end (); ++meit, ++copit)
			LinBox::RawVector<Element>::convert (*meit, *copit);
        }

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const
	{ 
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s+= LinBox::RawVector<_Element>::size (*it);

		return s;
        }

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void           setEntry (size_t i, size_t j, const Element &value);
	void           eraseEntry (size_t i, size_t j);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
			{ return x = getEntry (i, j); }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()       { return _A.begin (); }
	ConstRowIterator rowBegin () const { return _A.begin (); }
	RowIterator      rowEnd ()         { return _A.end (); }
	ConstRowIterator rowEnd () const   { return _A.end (); }

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::SparseSequenceVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::SparseSequenceVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd ()); }

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

/* Specialization for sparse associative vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::SparseAssociativeVectorTag >
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template <typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
	{
		typedef SparseMatrix<typename _Tp1::Element, _R1, VectorCategories::SparseAssociativeVectorTag> other;
	};

	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix (const SparseMatrix &A)
		: _A (A._A), _m (A._m), _n (A._n) {}

    	template <class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType, VectorCategories::SparseAssociativeVectorTag> &A)
		: _A(A.m), _m (A._m), _n (A._n)
	{
		typename Rep::iterator meit = this->_A.begin ();
		typename SparseMatrix<Element, VectorType, VectorCategories::SparseAssociativeVectorTag>::Rep::const_iterator copit = A._A.begin ();

		for (; meit != this->_A.end (); ++meit, ++copit)
			LinBox::RawVector<Element>::convert (*meit, *copit);
        }

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrix ( MatrixStream<Field>& ms );

	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const
	{ 
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s += LinBox::RawVector<_Element>::size (*it);

		return s;
        }

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void           setEntry (size_t i, size_t j, const Element &value) { _A[i][j] = value; }
	Element       &refEntry (size_t i, size_t j)                       { return _A[i][j]; }
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const     { return x = _A[i][j];}

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()       { return _A.begin (); }
	ConstRowIterator rowBegin () const { return _A.begin (); }
	RowIterator      rowEnd ()         { return _A.end (); }
	ConstRowIterator rowEnd () const   { return _A.end (); }

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::SparseAssociativeVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::SparseAssociativeVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd ()); }

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

/* Specialization for sparse parallel vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::SparseParallelVectorTag >
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrix<typename _Tp1::Element, _R1, VectorCategories::SparseParallelVectorTag> other; };

	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix (const SparseMatrix &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
    
    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType, VectorCategories::SparseParallelVectorTag> &A)
		: _A(A._m), _m (A._m), _n (A._n)
	{
		typename Rep::iterator meit = this->_A.begin ();
		typename SparseMatrix<Element, VectorType, VectorCategories::SparseParallelVectorTag>::Rep::const_iterator copit = A._A.begin ();

		for (; meit != this->_A.end (); ++meit, ++copit)
			LinBox::RawVector<Element>::convert (*meit, *copit);
        }

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrix ( MatrixStream<Field>& ms );

	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const
	{ 
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s+= LinBox::RawVector<_Element>::size (*it);

		return s;
        }

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
		{ return x = getEntry (i, j); }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator      rowBegin ()       { return _A.begin (); }
	ConstRowIterator rowBegin () const { return _A.begin (); }
	RowIterator      rowEnd ()         { return _A.end (); }
	ConstRowIterator rowEnd () const   { return _A.end (); }

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::SparseParallelVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::SparseParallelVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd ()); }

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

template <class Element, class Row, class Trait>
struct MatrixTraits< SparseMatrix<Element, Row, Trait> >
{ 
	typedef SparseMatrix<Element, Row, Trait> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag MatrixCategory; 
};

template <class Element, class Row, class Trait>
struct MatrixTraits< const SparseMatrix<Element, Row, Trait> >
{ 
	typedef const SparseMatrix<Element, Row, Trait> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/sparse.inl"

#endif // __LINBOX_matrix_sparse_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
