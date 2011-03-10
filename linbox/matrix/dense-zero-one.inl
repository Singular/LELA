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
class DenseZeroOneMatrixConstRowIterator
{
    public:
	typedef BitSubvectorWordAligned<ConstIterator, ConstIterator, Endianness> ConstRow;
	typedef ConstRow Row;

	typedef ConstRow value_type;

	typedef typename ConstRow::const_word_iterator::difference_type difference_type;

	DenseZeroOneMatrixConstRowIterator (ConstIterator p, size_t words, size_t bit_len, size_t d)
		: _row (p, p + words, bit_len), _disp (d) {}
    
	DenseZeroOneMatrixConstRowIterator () {}
    
	DenseZeroOneMatrixConstRowIterator (const DenseZeroOneMatrixConstRowIterator& colp)
		: _row (colp._row), _disp (colp._disp) {}
    
	DenseZeroOneMatrixConstRowIterator& operator = (const DenseZeroOneMatrixConstRowIterator& colp)
	{
		_row = colp._row;
		_disp = colp._disp;
		return *this;
	}

	DenseZeroOneMatrixConstRowIterator& operator --()
	{
		_row = ConstRow (_row.wordBegin () - _disp, _row.wordEnd () - _disp, _row.size ());
		return *this;
	}

	DenseZeroOneMatrixConstRowIterator  operator-- (int)
        {
                DenseZeroOneMatrixConstRowIterator tmp (*this);
                --*this;
                return tmp;
	}

	
	DenseZeroOneMatrixConstRowIterator& operator++ ()
	{
		_row = ConstRow (_row.wordBegin () + _disp, _row.wordEnd () + _disp, _row.size ());
		return *this;
	}

	DenseZeroOneMatrixConstRowIterator  operator++ (int)
	{
		DenseZeroOneMatrixConstRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	DenseZeroOneMatrixConstRowIterator operator+ (int i) const
		{ return ConstRowIterator (_row.wordBegin () + _disp * i, _row.word_size (), _row.size (), _disp); }

	difference_type operator- (const DenseZeroOneMatrixConstRowIterator &c) const
	{
		return (c._row.wordBegin () - _row.wordBegin ()) / _disp;
	}

	DenseZeroOneMatrixConstRowIterator& operator += (int i)
	{
		_row = ConstRow (_row.wordBegin () + _disp * i, _row.wordEnd () + _disp * i, _row.size ());
		return *this;
	}

	ConstRow operator[] (int i) const
		{ return ConstRow (_row.wordBegin () + _disp * i, _row.wordEnd () + _disp * i, _row.size ()); }

	const ConstRow* operator-> () const
		{ return &_row; }

	const ConstRow& operator* () const
		{ return _row; }

	bool operator == (const DenseZeroOneMatrixConstRowIterator& c) const
		{ return (_row.wordBegin () == c._row.wordBegin ()) && (_row.wordEnd () == c._row.wordEnd ()) && (_disp == c._disp); }
    
	bool operator!= (const DenseZeroOneMatrixConstRowIterator& c) const
		{ return (_row.wordBegin () != c._row.wordBegin ()) || (_row.wordEnd () != c._row.wordEnd ()) || (_disp != c._disp); }

    private:
	friend class DenseZeroOneMatrixRowIterator<Iterator, ConstIterator, Endianness>;

	ConstRow _row;
	size_t _disp;
};

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

	DenseZeroOneMatrixRowIterator (const DenseZeroOneMatrixRowIterator& colp)
		: _row (colp._row), _disp (colp._disp) {}
    
	DenseZeroOneMatrixRowIterator& operator = (const DenseZeroOneMatrixRowIterator& colp)
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

	Row* operator-> ()
		{ return &_row; }

	const Row* operator-> () const
		{ return &_row; }
    
	Row& operator* ()
		{ return _row; }

	const Row& operator* () const
		{ return _row; }
 
	bool operator == (const DenseZeroOneMatrixRowIterator& c) const
		{ return (_row.wordBegin () == c._row.wordBegin ()) && (_row.wordEnd () == c._row.wordEnd ()) && (_disp == c._disp); }

	bool operator != (const DenseZeroOneMatrixRowIterator& c) const
		{ return (_row.wordBegin () != c._row.wordBegin ()) || (_row.wordEnd () != c._row.wordEnd ()) || (_disp != c._disp); }

	bool operator != (const DenseZeroOneMatrixConstRowIterator<Iterator, ConstIterator, Endianness> &c) const
		{ return (_row.wordBegin () != c._row.wordBegin ()) || (_row.wordEnd () != c._row.wordEnd ()) || (_disp != c._disp); }

	operator DenseZeroOneMatrixConstRowIterator<Iterator, ConstIterator, Endianness> ()
		{ return DenseZeroOneMatrixConstRowIterator<Iterator, ConstIterator, Endianness> (_row.wordBegin (), _row.word_size (), _row.size (), _disp); }

    private:
	Row _row;
	size_t _disp;
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
