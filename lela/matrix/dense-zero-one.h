/* lela/matrix/dense-zero-one.h
 * Copyright 2010 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of DenseMatrix for zero-one dense matrices
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_DENSE_ZERO_ONE_H
#define __LELA_MATRIX_DENSE_ZERO_ONE_H

#include <vector>

#include "lela/vector/bit-vector.h"
#include "lela/vector/bit-subvector-word-aligned.h"
#include "lela/vector/bit-subvector.h"
#include "lela/matrix/submatrix.h"
#include "lela/matrix/raw-iterator.h"
#include "lela/matrix/dense.h"

#ifdef __LELA_HAVE_M4RI
#  include "lela/matrix/m4ri-matrix.h"
#else // !__LELA_HAVE_M4RI

namespace LELA
{

template <class Iterator, class ConstIterator, class Endianness>
class Dense01MatrixRowIterator
{
    public:
	typedef BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> Row;  
	typedef BitSubvectorWordAligned<ConstIterator, ConstIterator, Endianness> ConstRow;

	typedef Row value_type;

	typedef typename std::iterator_traits<typename Row::word_iterator>::difference_type difference_type;

	Dense01MatrixRowIterator (Iterator p, size_t words, size_t bit_len, size_t d)
		: _row (p, p + words, bit_len), _disp (d) {}

	Dense01MatrixRowIterator () {}

	Dense01MatrixRowIterator (const Dense01MatrixRowIterator &colp)
		: _row (colp._row), _disp (colp._disp) {}

	template <class It, class CIt>
	Dense01MatrixRowIterator (const Dense01MatrixRowIterator<It, CIt, Endianness> &colp)
		: _row (colp._row), _disp (colp._disp) {}
    
	template <class It, class CIt>
	Dense01MatrixRowIterator& operator = (const Dense01MatrixRowIterator<It, CIt, Endianness> &colp)
	{
		_row = colp._row;
		_disp = colp._disp;
		return *this;
	}
    
	Dense01MatrixRowIterator &operator ++ ()
	{
		_row = Row (_row.word_begin () + _disp, _row.word_begin () + (_disp + _row.word_size ()), _row.size ());
		return *this;
	}
    
	Dense01MatrixRowIterator  operator ++ (int)
	{
		Dense01MatrixRowIterator tmp (*this);
		++*this;
		return tmp;
	}
    
        Dense01MatrixRowIterator &operator -- ()
        {
                _row = Row (_row.word_begin () - _disp, _row.word_begin () + (_row.word_size () - _disp), _row.size ());
                return *this;
        }

        Dense01MatrixRowIterator  operator -- (int)
        {
                Dense01MatrixRowIterator tmp (*this);
                --*this;
                return tmp;
        }

	Dense01MatrixRowIterator operator + (int i) const
		{ return Dense01MatrixRowIterator (const_cast<Row *> (&_row)->word_begin () + _disp * i, _row.word_size (), _row.size (), _disp); }

	difference_type operator - (const Dense01MatrixRowIterator &c) const
		{ return (c._row.word_begin () - _row.word_begin ()) / _disp; }

	Dense01MatrixRowIterator& operator += (int i)
	{
		_row = Row (_row.word_begin () + _disp * i, _row.word_end () + (_disp * i + _row.word_size ()), _row.size ());
		return *this;
	}

	Row operator[] (int i) const
		{ return Row (const_cast<Row&> (_row).word_begin () + _disp * i,
			      const_cast<Row&> (_row).word_begin () + (_disp * i + _row.word_size ()), _row.size ()); }

	Row *operator -> ()
		{ return &_row; }

	Row &operator * ()
		{ return _row; }

	bool operator == (const Dense01MatrixRowIterator& c) const
		{ return (_row.word_begin () == c._row.word_begin ()) && (_row.word_end () == c._row.word_end ()) && (_disp == c._disp); }

	template <class It, class CIt>
	bool operator == (const Dense01MatrixRowIterator<It, CIt, Endianness> &c) const
		{ return (_row.word_begin () != c._row.word_begin ()) || (_row.word_end () != c._row.word_end ()) || (_disp != c._disp); }

	bool operator != (const Dense01MatrixRowIterator& c) const
		{ return (_row.word_begin () != c._row.word_begin ()) || (_row.word_end () != c._row.word_end ()) || (_disp != c._disp); }

	template <class It, class CIt>
	bool operator != (const Dense01MatrixRowIterator<It, CIt, Endianness> &c) const
		{ return (_row.word_begin () != c._row.word_begin ()) || (_row.word_end () != c._row.word_end ()) || (_disp != c._disp); }

    private:
	Row _row;
	size_t _disp;

	template <class It, class CIt, class E>
	friend class Dense01MatrixRowIterator;
};

template <class Iterator = BitVector<DefaultEndianness<uint64> >::word_iterator,
	  class ConstIterator = typename BitVector<DefaultEndianness<uint64> >::const_word_iterator,
	  class Endianness = DefaultEndianness<uint64> >
class Dense01Matrix
{
    public:

	typedef bool Element;
	typedef BitVector<Endianness> Rep;
        typedef Dense01Matrix Self_t;
	typedef typename std::iterator_traits<Iterator>::value_type word_type;
	typedef MatrixIteratorTypes::Row IteratorType; 
	typedef MatrixStorageTypes::Rows StorageType; 

	typedef Submatrix<Self_t> SubmatrixType;
	typedef Submatrix<const Self_t> ConstSubmatrixType;
	typedef Self_t AlignedSubmatrixType;
	typedef const Self_t ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = WordTraits<typename Iterator::value_type>::bits;

	typedef Self_t ContainerType;

	typedef Dense01MatrixRowIterator<Iterator, ConstIterator, Endianness> RowIterator;
	typedef Dense01MatrixRowIterator<ConstIterator, ConstIterator, Endianness> ConstRowIterator;

	typedef typename RowIterator::value_type Row;  
	typedef typename ConstRowIterator::value_type ConstRow;

	template <class It1, class It2, class End>
	friend class Dense01Matrix;

	///
	Dense01Matrix ()
		: _rows (0), _cols (0)
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	Dense01Matrix (size_t m, size_t n)
		: _rows (m), _cols (n)
		{ init_disp_rep (m, n); _begin = _rep.word_begin (); }

	/** Construct a word-aligned dense submatrix of a given dense matrix
	 */
	template <class It, class CIt>
	Dense01Matrix (Dense01Matrix<It, CIt, Endianness> &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		lela_check (beg_col & WordTraits<word_type>::pos_mask == 0);

		_begin = M._begin + beg_row * M._disp + (beg_col >> WordTraits<word_type>::logof_size);
	}

	Dense01Matrix (Dense01Matrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		lela_check ((beg_col & WordTraits<word_type>::pos_mask) == 0);

		_begin = M._begin + beg_row * M._disp + (beg_col >> WordTraits<word_type>::logof_size);
	}

	/** Construct a dense matrix, filling the rows in from a VectorStream. The stream's size must be finite
	 */
	Dense01Matrix (VectorStream<Row> &vs)
		: _rows (vs.size ()), _cols (vs.dim ())
	{
		init_disp_rep (vs.size (), vs.dim ());

		_begin = _rep.word_begin ();

		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	/** Version of above for const matrices
	 */
	template <class It, class CIt>
	Dense01Matrix (const Dense01Matrix<It, CIt, Endianness> &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rows (m), _cols (n), _disp (M._disp)
	{
		_begin = M._begin + beg_row * M._disp + (beg_col >> WordTraits<word_type>::logof_size);
	}

	///
	Dense01Matrix (const Dense01Matrix &M)
		: _rep (M._rep), _rows (M._rows), _cols (M._cols), _disp (M._disp), _begin (M._begin)
	{}

	~Dense01Matrix(){}
	///
	template <class It1, class It2>
	Dense01Matrix& operator= (const Dense01Matrix<It1, It2>& M) {
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
		_begin = _rep.word_begin ();
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

	RowIterator rowBegin ()
		{ return RowIterator (_begin, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp); }
	RowIterator rowEnd ()
		{ return RowIterator (_begin + _rows * _disp, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp); }
	ConstRowIterator rowBegin () const
		{ return ConstRowIterator (_begin, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp); }
	ConstRowIterator rowEnd () const
		{ return ConstRowIterator (_begin + _rows * _disp, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp); }

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	typedef MatrixRawIterator<ConstRowIterator, VectorRepresentationTypes::Dense01> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorRepresentationTypes::Dense01, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

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

#ifndef __LELA_HAVE_M4RI

template <>
class DenseMatrix<bool> : public Dense01Matrix<>
{
public:
	DenseMatrix ()
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (size_t m, size_t n)
		: Dense01Matrix<> (m, n)
	{}

	///
	DenseMatrix (const DenseMatrix &M)
		: Dense01Matrix<> (M)
	{}

	DenseMatrix (VectorStream<Row> &vs)
		: Dense01Matrix<> (vs)
	{}
};

#endif // !__LELA_HAVE_M4RI

} // namespace LELA

#endif // __LELA_HAVE_M4RI

#endif // __LELA_MATRIX_DENSE_ZERO_ONE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
