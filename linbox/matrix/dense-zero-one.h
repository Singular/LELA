/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/dense-zero-one.h
 * Copyright 2010 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of DenseMatrixBase for zero-one dense matrices
 *
 * Evolved from dense.h
 */

#ifndef __MATRIX_DENSE_ZERO_ONE_H
#define __MATRIX_DENSE_ZERO_ONE_H

#include <vector>

#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector-word-aligned.h"
#include "linbox/vector/bit-subvector.h"
#include "linbox/matrix/submatrix.h"

namespace LinBox
{

template <class Iterator = BitVector<DefaultEndianness>::word_iterator,
	  class ConstIterator = typename BitVector<DefaultEndianness>::const_word_iterator,
	  class Endianness = DefaultEndianness>
class DenseZeroOneMatrix
{
    public:

	typedef bool Element;
	typedef BitVector<Endianness> Rep;
        typedef DenseZeroOneMatrix Self_t;

	template <class It1, class It2, class End>
	friend class DenseZeroOneMatrix;

	///
	DenseZeroOneMatrix ()
		: _rows (0), _cols (0)
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseZeroOneMatrix (size_t m, size_t n)
		: _rows (m), _cols (n),
		  _disp ((!(n & __LINBOX_POS_ALL_ONES)) ? (n >> __LINBOX_LOGOF_SIZE) : ((n >> __LINBOX_LOGOF_SIZE) + 1)),
		  _rep (m * ((!(n & __LINBOX_POS_ALL_ONES)) ? (n >> __LINBOX_LOGOF_SIZE) : ((n >> __LINBOX_LOGOF_SIZE) + 1)) * __LINBOX_BITSOF_LONG)
		{ _begin = _rep.wordBegin (); }

	/** Construct a word-aligned dense submatrix of a given dense matrix
	 */
	template <class It1, class It2>
	DenseZeroOneMatrix (DenseZeroOneMatrix<It1, It2> &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		_begin = M._begin + beg_row * M._disp + (beg_col >> __LINBOX_LOGOF_SIZE);
	}

	/** Version of above for const matrices
	 */
	template <class It1, class It2>
	DenseZeroOneMatrix (const DenseZeroOneMatrix<It1, It2> &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		_begin = M._begin + beg_row * M._disp + (beg_col >> __LINBOX_LOGOF_SIZE);
	}

	///
	DenseZeroOneMatrix (const DenseZeroOneMatrix &M)
		: _rep (M._rep), _rows (M._rows), _cols (M._cols), _begin (M._begin), _disp (M._disp)
	{}

	~DenseZeroOneMatrix(){}
	///
	template <class It1, class It2>
	DenseZeroOneMatrix& operator= (const DenseZeroOneMatrix<It1, It2>& M) {
		(*this)._rep  = M._rep;
		(*this)._rows = M._rows;
		(*this)._cols = M._cols;
		(*this)._disp  = M._disp;
		return (*this);
	}

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
	void resize (size_t m, size_t n)
	{
		_rows = m;
		_cols = n;
		_disp = (!(n & __LINBOX_POS_ALL_ONES)) ? (n >> __LINBOX_LOGOF_SIZE) : ((n >> __LINBOX_LOGOF_SIZE) + 1);
		_rep.resize (m * _disp * __LINBOX_BITSOF_LONG);
		_begin = _rep.wordBegin ();
	}

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 * @param F Field over which to read
	 */
	template <class Field>
	std::istream &read (std::istream &file, const Field &F);

	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 * @param F Field over which to write
	 */
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F) const;

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, bool a_ij)
		{ *(BitVectorIterator<Iterator, ConstIterator, Endianness> (_begin, 0) + (i * _disp * __LINBOX_BITSOF_LONG + j)) = a_ij; }

	/** Get a writeable reference to the entry in the (i, j) position.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @returns Reference to matrix entry
	 */
	BitVectorReference<Iterator, Endianness> refEntry (size_t i, size_t j)
		{ return *(BitVectorIterator<Iterator, ConstIterator, Endianness> (_begin, 0) + (i * _disp * __LINBOX_BITSOF_LONG + j)); }

	/** Get a read-only reference to the entry in the (i, j) position.
	 * @param i Row index
	 * @param j Column index
	 * @returns Const reference to matrix entry
	 */
	bool getEntry (size_t i, size_t j) const
		{ return *(BitVectorConstIterator<ConstIterator, Endianness> (_begin, 0) + (i * _disp * __LINBOX_BITSOF_LONG + j)); }

	/** Copy the (i, j) entry into x, and return a reference to x.
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @returns Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const
		{ x = *(BitVectorConstIterator<ConstIterator, Endianness> (_begin, 0) + (i * _disp * __LINBOX_BITSOF_LONG + j)); return x; }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	typedef BitSubvectorWordAligned<Iterator, ConstIterator, LittleEndian<typename Iterator::value_type> > Row;  
	typedef BitSubvectorWordAligned<ConstIterator, ConstIterator, LittleEndian<typename ConstIterator::value_type> > ConstRow;

	class RowIterator;
	class ConstRowIterator;

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;

	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 * @param i Row index
	 */
	Row operator[] (size_t i)
		{ return Row (_begin + _disp * i, _begin + _disp * i + (_cols >> __LINBOX_LOGOF_SIZE), _cols); }

	ConstRow operator[] (size_t i) const
		{ return ConstRow (_begin + _disp * i, _begin + _disp * i + (_cols >> __LINBOX_LOGOF_SIZE), _cols); }

	/** Compute column density
	 */

	template <class Vector>
	Vector &columnDensity (Vector &v) const
		{ std::fill (v.begin (), v.end (), _rows); return v; }

    protected:

	Rep    _rep;
	size_t _rows, _cols;
	size_t _disp;
	Iterator _begin;
};

template <class Iterator, class ConstIterator>
struct MatrixTraits< DenseZeroOneMatrix<Iterator, ConstIterator> >
{ 
	typedef DenseZeroOneMatrix<Iterator, ConstIterator> MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

template <class Iterator, class ConstIterator>
struct MatrixTraits< const DenseZeroOneMatrix<Iterator, ConstIterator> >
{ 
	typedef const DenseZeroOneMatrix<Iterator, ConstIterator> MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/dense-zero-one.inl"

#endif // __MATRIX_DENSE_ZERO_ONE_H
