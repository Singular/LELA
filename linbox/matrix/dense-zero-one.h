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
#include "linbox/matrix/raw-iterator.h"
#include "linbox/matrix/dense.h"

namespace LinBox
{

template <class Iterator, class ConstIterator, class Endianness>
class DenseZeroOneMatrixRowIterator
{
    public:
	typedef BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> Row;  
	typedef BitSubvectorWordAligned<ConstIterator, ConstIterator, Endianness> ConstRow;

	typedef Row value_type;

	typedef typename Row::word_iterator::difference_type difference_type;

	DenseZeroOneMatrixRowIterator (Iterator p, size_t words, size_t bit_len, size_t d)
		: _row (p, p + words, bit_len), _disp (d) {}

	DenseZeroOneMatrixRowIterator () {}

	DenseZeroOneMatrixRowIterator (const DenseZeroOneMatrixRowIterator &colp)
		: _row (colp._row), _disp (colp._disp) {}

	template <class It, class CIt>
	DenseZeroOneMatrixRowIterator (const DenseZeroOneMatrixRowIterator<It, CIt, Endianness> &colp)
		: _row (colp._row), _disp (colp._disp) {}
    
	template <class It, class CIt>
	DenseZeroOneMatrixRowIterator& operator = (const DenseZeroOneMatrixRowIterator<It, CIt, Endianness> &colp)
	{
		_row = colp._row;
		_disp = colp._disp;
		return *this;
	}
    
	DenseZeroOneMatrixRowIterator& operator ++ ()
	{
		_row = Row (_row.wordBegin () + _disp, _row.wordEnd () + _disp, _row.size ());
		return *this;
	}
    
	DenseZeroOneMatrixRowIterator  operator ++ (int)
	{
		DenseZeroOneMatrixRowIterator tmp (*this);
		++*this;
		return tmp;
	}
    
        DenseZeroOneMatrixRowIterator& operator -- ()
        {
                _row = Row (_row.wordBegin () - _disp, _row.wordEnd () - _disp, _row.size ());
                return *this;
        }

        DenseZeroOneMatrixRowIterator  operator -- (int)
        {
                DenseZeroOneMatrixRowIterator tmp (*this);
                --*this;
                return tmp;
        }

	DenseZeroOneMatrixRowIterator operator + (int i) const
		{ return DenseZeroOneMatrixRowIterator (const_cast<Row *> (&_row)->wordBegin () + _disp * i, _row.word_size (), _row.size (), _disp); }

	difference_type operator- (const DenseZeroOneMatrixRowIterator &c) const
	{
		return (c._row.wordBegin () - _row.wordBegin ()) / _disp;
	}

	DenseZeroOneMatrixRowIterator& operator += (int i)
	{
		_row = Row (_row.wordBegin () + _disp * i, _row.wordEnd () + _disp * i, _row.size ());
		return *this;
	}

	Row operator[] (int i) const
		{ return Row (const_cast<Row&> (_row).wordBegin () + _disp * i,
			      const_cast<Row&> (_row).wordEnd () + _disp * i, _row.size ()); }

	Row *operator -> ()
		{ return &_row; }

	Row &operator * ()
		{ return _row; }

	bool operator == (const DenseZeroOneMatrixRowIterator& c) const
		{ return (_row.wordBegin () == c._row.wordBegin ()) && (_row.wordEnd () == c._row.wordEnd ()) && (_disp == c._disp); }

	template <class It, class CIt>
	bool operator == (const DenseZeroOneMatrixRowIterator<It, CIt, Endianness> &c) const
		{ return (_row.wordBegin () != c._row.wordBegin ()) || (_row.wordEnd () != c._row.wordEnd ()) || (_disp != c._disp); }

	bool operator != (const DenseZeroOneMatrixRowIterator& c) const
		{ return (_row.wordBegin () != c._row.wordBegin ()) || (_row.wordEnd () != c._row.wordEnd ()) || (_disp != c._disp); }

	template <class It, class CIt>
	bool operator != (const DenseZeroOneMatrixRowIterator<It, CIt, Endianness> &c) const
		{ return (_row.wordBegin () != c._row.wordBegin ()) || (_row.wordEnd () != c._row.wordEnd ()) || (_disp != c._disp); }

    private:
	Row _row;
	size_t _disp;

	template <class It, class CIt, class E>
	friend class DenseZeroOneMatrixRowIterator;
};

template <class Iterator = BitVector<DefaultEndianness>::word_iterator,
	  class ConstIterator = typename BitVector<DefaultEndianness>::const_word_iterator,
	  class Endianness = DefaultEndianness>
class DenseZeroOneMatrix
{
    public:

	typedef bool Element;
	typedef BitVector<Endianness> Rep;
        typedef DenseZeroOneMatrix Self_t;
	typedef typename std::iterator_traits<Iterator>::value_type word_type;

	typedef DenseZeroOneMatrixRowIterator<Iterator, ConstIterator, Endianness> RowIterator;
	typedef DenseZeroOneMatrixRowIterator<ConstIterator, ConstIterator, Endianness> ConstRowIterator;

	typedef typename RowIterator::value_type Row;  
	typedef typename ConstRowIterator::value_type ConstRow;

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
		: _rows (m), _cols (n)
		{ init_disp_rep (m, n); _begin = _rep.wordBegin (); }

	/** Construct a word-aligned dense submatrix of a given dense matrix
	 */
	template <class It1, class It2>
	DenseZeroOneMatrix (DenseZeroOneMatrix<It1, It2> &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		_begin = M._begin + beg_row * M._disp + (beg_col >> WordTraits<word_type>::logof_size);
	}

	/** Construct a dense matrix, filling the rows in from a VectorStream. The stream's size must be finite
	 */
	DenseZeroOneMatrix (VectorStream<Row> &vs)
		: _rows (vs.size ()), _cols (vs.dim ())
	{
		init_disp_rep (vs.size (), vs.dim ());

		_begin = _rep.wordBegin ();

		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	/** Version of above for const matrices
	 */
	template <class It1, class It2>
	DenseZeroOneMatrix (const DenseZeroOneMatrix<It1, It2> &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		_begin = M._begin + beg_row * M._disp + (beg_col >> WordTraits<word_type>::logof_size);
	}

	///
	DenseZeroOneMatrix (const DenseZeroOneMatrix &M)
		: _rep (M._rep), _rows (M._rows), _cols (M._cols), _disp (M._disp), _begin (M._begin)
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
		_disp = (!(n & WordTraits<word_type>::pos_mask)) ? (n >> WordTraits<word_type>::logof_size) : ((n >> WordTraits<word_type>::logof_size) + 1);
		_rep.resize (m * _disp * WordTraits<word_type>::bits);
		_begin = _rep.wordBegin ();
	}

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, bool a_ij)
		{ *(BitVectorIterator<Iterator, ConstIterator, Endianness> (_begin, 0) + (i * _disp * WordTraits<word_type>::bits + j)) = a_ij; }

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
		{ x = *(BitVectorConstIterator<ConstIterator, Endianness> (_begin, 0) + (i * _disp * WordTraits<word_type>::bits + j)); return true; }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::DenseZeroOneVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd ()); }

	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::DenseZeroOneVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd ()); }

	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 * @param i Row index
	 */
	Row operator[] (size_t i)
		{ return Row (_begin + _disp * i, _begin + _disp * i + (_cols >> WordTraits<word_type>::logof_size), _cols); }

	ConstRow operator[] (size_t i) const
		{ return ConstRow (_begin + _disp * i, _begin + _disp * i + (_cols >> WordTraits<word_type>::logof_size), _cols); }

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

    private:

	void init_disp_rep (size_t m, size_t n)
	{
		_disp = n >> WordTraits<word_type>::logof_size;

		if (n & WordTraits<word_type>::pos_mask)
			++_disp;

		_rep.resize (m * _disp * WordTraits<word_type>::bits);
	}
};

#ifndef __LINBOX_HAVE_M4RI

template <>
class DenseMatrix<bool> : public DenseZeroOneMatrix<>
{
public:
	DenseMatrix ()
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (size_t m, size_t n)
		: DenseZeroOneMatrix<> (m, n)
	{}

	///
	DenseMatrix (const DenseMatrix &M)
		: DenseZeroOneMatrix<> (M)
	{}

	DenseMatrix (VectorStream<Row> &vs)
		: DenseZeroOneMatrix<> (vs)
	{}
};

#endif // !__LINBOX_HAVE_M4RI

template <class Iterator, class ConstIterator, class Endianness>
struct MatrixTraits< DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> >
{
	typedef DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

template <class Iterator, class ConstIterator, class Endianness>
struct MatrixTraits< const DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> >
{ 
	typedef const DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

template <>
struct MatrixTraits<DenseMatrix<bool> >
{ 
	typedef DenseMatrix<bool> MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

template <>
struct MatrixTraits<const DenseMatrix<bool> >
{ 
	typedef const DenseMatrix<bool> MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/dense-zero-one.inl"

#endif // __MATRIX_DENSE_ZERO_ONE_H
