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
