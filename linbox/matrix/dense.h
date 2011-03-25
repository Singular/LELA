/* linbox/matrix/dense.h
 * Copyright (C) 2001 B. David Saunders,
 *               2001-2002 Bradford Hovinen,
 *               2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/dense-base.h to matrix/dense.h
 * --------------------------------------------------------
 * 2002-11-29  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Swap the order of arguments in read and write, so that it is consistent with
 * SparseMatrixBase
 * --------------------------------------------------------
 * 2002-10-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Rename ColOfRowsIterator as RowIterator; similarly with RowOfColsIterator
 * --------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Split out container/iterator functionality into DenseMatrixBase
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_matrix_dense_H
#define __LINBOX_matrix_dense_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/raw-iterator.h"
#include "linbox/linbox-config.h"

namespace LinBox
{

template <class Iterator, class ConstIterator>
class DenseMatrixRowIterator;

template <class Iterator, class ConstIterator>
class DenseMatrixColIterator;

/** Blackbox dense matrix template. This is a class of dense matrices
 * templatized by the entry type, the Element type of some {@link Fields field}.
 * The matrix is stored as a one dimensional STL vector of the elements, by rows. 
 * The interface provides for iteration over rows and over columns.
 *
 * The class LinBox::Dense builds on this base.
 *
 * Currently, only dense vectors are supported when doing matrix-vector applies.
 *
\ingroup matrix
 */
  
template <class _Element>
class DenseMatrix
{
    public:

	typedef _Element Element;
	typedef typename RawVector<Element>::Dense Rep;
        typedef DenseMatrix<_Element> Self_t;

	typedef DenseMatrixRowIterator<typename Rep::iterator, typename Rep::const_iterator> RowIterator;
	typedef DenseMatrixRowIterator<typename Rep::const_iterator, typename Rep::const_iterator> ConstRowIterator;

	typedef Subvector<typename Rep::iterator, typename Rep::const_iterator> Row;
	typedef Subvector<typename Rep::const_iterator> ConstRow;  

	typedef DenseMatrixColIterator<typename Rep::iterator, typename Rep::const_iterator> ColIterator;
	typedef DenseMatrixColIterator<typename Rep::const_iterator, typename Rep::const_iterator> ConstColIterator;

	typedef Subvector<Subiterator<typename Rep::iterator> > Col;
	typedef Subvector<Subiterator<typename Rep::const_iterator> > ConstCol;
	typedef Col Column;
	typedef ConstCol ConstColumn;

	template<typename _Tp1>
        struct rebind
        { 
            typedef DenseMatrix<typename _Tp1::Element> other; 
        };

	///
	DenseMatrix ()
		: _rows (0), _cols (0)
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (size_t m, size_t n)
		: _rep (m * n), _rows (m), _cols (n), _ptr(&_rep[0])
	{}

	///
	DenseMatrix (const DenseMatrix &M)
		: _rep (M._rep),_rows (M._rows), _cols (M._cols), _ptr(&_rep[0])
	{}

	DenseMatrix (VectorStream<Row> &vs);

	~DenseMatrix(){}
	///
	DenseMatrix& operator= (const DenseMatrix& M) {
		(*this)._rep  = M._rep;
		(*this)._rows = M._rows;
		(*this)._cols = M._cols;
		(*this)._ptr  = &_rep[0];
		return (*this);
	}

	/** Get a pointer on the storage of the elements
	 * @returns a pointer on Elements
	/todo What is this?
	 */
	Element* FullIterator() const {return const_cast<Element*>(&_rep[0]);}

	/** Get the number of rows in the matrix
	 * @returns Number of rows in matrix
	 */
	size_t rowdim () const
		{ return _rows; }

	/** Get the number of columns in the matrix
	 * @returns Number of columns in matrix
	 */
	size_t coldim () const
		{ return _cols; }

	/** Resize the matrix to the given dimensions
	 * The state of the matrix's entries after a call to this method is
	 * undefined
	 * @param m Number of rows
	 * @param n Number of columns
	 */
	void resize (size_t m, size_t n, const Element& val = Element())
	{
		_rows = m;
		_cols = n;
		_rep.resize (m * n, val);
	}

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _rep[i * _cols + j] = a_ij; }

	/** Does nothing. Provided for compatibility
	 * @param i
	 * @param j
	 */
	void eraseEntry (size_t i, size_t j) {}

	/** Copy the (i, j) entry into x, and return a reference to x.
	 *
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @returns true
	 */
	bool getEntry (Element &x, size_t i, size_t j) const
		{ x = _rep[i * _cols + j]; return true; }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	RowIterator rowBegin ();  
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;        
	ConstRowIterator rowEnd () const;

	/** @name Row of columns iterator
	 * The row of columns iterator traverses the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */

	ColIterator colBegin ();
	ColIterator colEnd ();
	ConstColIterator colBegin () const;    
	ConstColIterator colEnd () const;

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	typedef typename Rep::iterator RawIterator;
	typedef typename Rep::const_iterator ConstRawIterator;
    
	RawIterator rawBegin ();		  
	RawIterator rawEnd ();
	ConstRawIterator rawBegin () const;
	ConstRawIterator rawEnd () const;

	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::DenseVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd ()); }
    
	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 * @param i Row index
	 */
	Row operator[] (size_t i)
		{ return Row (_rep.begin () + i * _cols, _rep.begin () + i * _cols + _cols); }

	ConstRow operator[] (size_t i) const
		{ return Row (_rep.begin () + i * _cols, _rep.begin () + i * _cols + _cols); }

	/** Compute column density
	 */

	template <class Vector>
	Vector &columnDensity (Vector &v) const
		{ std::fill (v.begin (), v.end (), _rows); return v; }

    protected:

	std::vector<Element>  _rep;
	size_t                _rows, _cols;
	Element *             _ptr;
};

template <class Element>
struct MatrixTraits< DenseMatrix<Element> >
{ 
	typedef DenseMatrix<Element> MatrixType;
	typedef typename MatrixCategories::RowColMatrixTag MatrixCategory; 
};

template <class Element>
struct MatrixTraits< const DenseMatrix<Element> >
{ 
	typedef const DenseMatrix<Element> MatrixType;
	typedef typename MatrixCategories::RowColMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "dense.inl"

#endif // __LINBOX_matrix_dense_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
