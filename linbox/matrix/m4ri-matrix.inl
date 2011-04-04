/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/m4ri-matrix.inl
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Wrapper for libm4ri-matrices
 *
 * Evolved from dense.h
 */

#ifndef __LINBOX_MATRIX_M4RI_MATRIX_INL
#define __LINBOX_MATRIX_M4RI_MATRIX_INL

#include <iostream>

#include "linbox/matrix/dense-zero-one.h"

namespace LinBox
{

template <class Iterator, class ConstIterator, class MatrixPointer>
class M4RIMatrixRowIterator
{
    public:
	typedef BitSubvectorWordAligned<Iterator, ConstIterator, BigEndian<word> > Row;  
	typedef BitSubvectorWordAligned<ConstIterator, ConstIterator, BigEndian<word> > ConstRow;

	typedef Row value_type;

	typedef typename std::iterator_traits<typename Row::word_iterator>::difference_type difference_type;

	M4RIMatrixRowIterator (MatrixPointer M, size_t idx)
		: _M (M), _idx (idx) { make_row (); }
    
	M4RIMatrixRowIterator () {}
    
	M4RIMatrixRowIterator (const M4RIMatrixRowIterator &colp)
		: _M (colp._M), _idx (colp._idx), _row (colp._row) {}

	template <class It, class CIt, class MP>
	M4RIMatrixRowIterator (const M4RIMatrixRowIterator<It, CIt, MP> &colp)
		: _M (colp._M), _idx (colp._idx), _row (colp._row) {}
    
	template <class It, class CIt, class MP>
	M4RIMatrixRowIterator &operator = (const M4RIMatrixRowIterator<It, CIt, MP> &colp)
	{
		_M = colp._M;
		_idx = colp._idx;
		_row = colp._row;
		return *this;
	}
    
	M4RIMatrixRowIterator &operator ++ ()
	{
		++_idx;
		make_row ();
		return *this;
	}
    
	M4RIMatrixRowIterator operator ++ (int)
	{
		M4RIMatrixRowIterator tmp (*this);
		++*this;
		return tmp;
	}
    
        M4RIMatrixRowIterator &operator -- ()
        {
		--_idx;
		make_row ();
		return *this;
        }

        M4RIMatrixRowIterator operator -- (int)
        {
                M4RIMatrixRowIterator tmp (*this);
                --*this;
                return tmp;
        }

	M4RIMatrixRowIterator operator+ (int i) const
		{ return M4RIMatrixRowIterator (_M, _idx + i); }

	difference_type operator- (const M4RIMatrixRowIterator &c) const
		{ return c._idx - _idx; }

	M4RIMatrixRowIterator &operator += (int i)
	{
		_idx += i;
		make_row ();
		return *this;
	}

	Row operator [] (int i) const
		{ return Row (_M->_rep->rows[_idx + i], _M->_rep->rows[_idx + i] + _M->_rep->width, _M->_rep->ncols); }

	Row *operator -> ()
		{ return &_row; }

	Row &operator * ()
		{ return _row; }

	bool operator == (const M4RIMatrixRowIterator& c) const
		{ return _idx == c._idx; }

	template <class It, class CIt, class MP>
	bool operator == (const M4RIMatrixRowIterator<It, CIt, MP> &c) const
		{ return _idx == c._idx; }

	bool operator != (const M4RIMatrixRowIterator& c) const
		{ return _idx != c._idx; }

	template <class It, class CIt, class MP>
	bool operator != (const M4RIMatrixRowIterator<It, CIt, MP> &c) const
		{ return _idx != c._idx; }

    private:
	inline void make_row ()
		{ if (_M->_rep->rows != NULL) _row = Row (_M->_rep->rows[_idx], _M->_rep->rows[_idx] + _M->_rep->width, _M->_rep->ncols); }

	MatrixPointer _M;
	size_t _idx;
	Row _row;

	template <class It, class CIt, class MP>
	friend class M4RIMatrixRowIterator;
};

M4RIMatrix::M4RIMatrix (VectorStream<Row> &vs)
	: M4RIMatrixBase (mzd_init (vs.size (), vs.dim ()))
{
	for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
		vs >> *i;
}

inline M4RIMatrixBase::RowIterator M4RIMatrixBase::rowBegin ()
{
	return RowIterator (this, 0);
}

inline M4RIMatrixBase::RowIterator M4RIMatrixBase::rowEnd ()
{
	return RowIterator (this, rowdim ());
}

inline M4RIMatrixBase::ConstRowIterator M4RIMatrixBase::rowBegin () const
{
	return ConstRowIterator (this, 0);
}

inline M4RIMatrixBase::ConstRowIterator M4RIMatrixBase::rowEnd () const
{
	return ConstRowIterator (this, rowdim ());
}

inline M4RIMatrixBase::ConstRawIterator M4RIMatrixBase::rawBegin () const
	{ return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
inline M4RIMatrixBase::ConstRawIterator M4RIMatrixBase::rawEnd () const
	{ return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

inline M4RIMatrixBase::ConstRawIndexedIterator M4RIMatrixBase::rawIndexedBegin() const
	{ return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
inline M4RIMatrixBase::ConstRawIndexedIterator M4RIMatrixBase::rawIndexedEnd() const
	{ return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

SubvectorFactory<M4RIMatrixBase>::RowSubvector SubvectorFactory<M4RIMatrixBase>::MakeRowSubvector (Submatrix<M4RIMatrixBase> &M, M4RIMatrixBase::RowIterator &pos)
	{ return RowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

SubvectorFactory<M4RIMatrixBase>::ConstRowSubvector SubvectorFactory<M4RIMatrixBase>::MakeConstRowSubvector (const Submatrix<M4RIMatrixBase> &M, M4RIMatrixBase::ConstRowIterator &pos)
	{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

SubvectorFactory<const M4RIMatrixBase>::ConstRowSubvector SubvectorFactory<const M4RIMatrixBase>::MakeConstRowSubvector (const Submatrix<const M4RIMatrixBase> &M, M4RIMatrixBase::ConstRowIterator &pos)
	{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

}

#endif // __LINBOX_MATRIX_M4RI_MATRIX_INL
