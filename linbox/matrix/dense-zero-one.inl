/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/dense-zero-one.inl
 * Copyright 2010 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of DenseMatrixBase for zero-one dense matrices -- support-functions
 *
 * Evolved from dense.h
 */

#ifndef __MATRIX_DENSE_ZERO_ONE_INL
#define __MATRIX_DENSE_ZERO_ONE_INL

#include <iostream>

#include <linbox/matrix/dense-zero-one.h>

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

template <class Iterator, class ConstIterator, class Endianness>
typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::RowIterator
DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::rowBegin ()
{
	return RowIterator (_begin, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp);
}

template <class Iterator, class ConstIterator, class Endianness>
typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::RowIterator
DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::rowEnd ()
{
	return RowIterator (_begin + _rows * _disp, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp);
}

template <class Iterator, class ConstIterator, class Endianness>
typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::ConstRowIterator
DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::rowBegin () const
{
	return ConstRowIterator (_begin, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp);
}

template <class Iterator, class ConstIterator, class Endianness>
typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::ConstRowIterator
DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::rowEnd () const
{
	return ConstRowIterator (_begin + _rows * _disp, (!(_cols & WordTraits<word_type>::pos_mask)) ? (_cols >> WordTraits<word_type>::logof_size) : ((_cols >> WordTraits<word_type>::logof_size) + 1), _cols, _disp);
}

template <class Iterator, class ConstIterator, class Endianness>
class SubvectorFactory<DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> >
{
    public:
	typedef BitSubvector<BitVectorIterator<Iterator, ConstIterator, Endianness>, BitVectorConstIterator<ConstIterator, Endianness> > RowSubvector;
	typedef BitSubvector<BitVectorConstIterator<ConstIterator, Endianness>, BitVectorConstIterator<ConstIterator, Endianness> > ConstRowSubvector;

	RowSubvector MakeRowSubvector (Submatrix<DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> > &M,
				       typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::RowIterator &pos)
		{ return RowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> > &M,
						 typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::ConstRowIterator &pos)
		{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

template <class Iterator, class ConstIterator, class Endianness>
class SubvectorFactory<const DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> >
{
    public:
	typedef BitSubvector<BitVectorConstIterator<ConstIterator, Endianness>, BitVectorConstIterator<ConstIterator, Endianness> > RowSubvector;
	typedef BitSubvector<BitVectorConstIterator<ConstIterator, Endianness>, BitVectorConstIterator<ConstIterator, Endianness> > ConstRowSubvector;

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<const DenseZeroOneMatrix<Iterator, ConstIterator, Endianness> > &M,
						 typename DenseZeroOneMatrix<Iterator, ConstIterator, Endianness>::ConstRowIterator &pos)
		{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

}

#endif // __MATRIX_DENSE_ZERO_ONE_INL
