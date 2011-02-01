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

class M4RIMatrix::ConstRowIterator
{
    public:
	typedef int difference_type;

	ConstRowIterator (const M4RIMatrix &M, size_t idx)
		: _M (&M), _idx (idx) { make_row (); }
    
	ConstRowIterator () {}
    
	ConstRowIterator (const ConstRowIterator& colp)
		: _M (colp._M), _idx (colp._idx), _row (colp._row) {}
    
	ConstRowIterator& operator = (const ConstRowIterator& colp)
	{
		_M = colp._M;
		_idx = colp._idx;
		_row = colp._row;
		return *this;
	}

	ConstRowIterator& operator --()
	{
		--_idx;
		make_row ();
		return *this;
	}

	ConstRowIterator  operator-- (int)
        {
                ConstRowIterator tmp (*this);
                --*this;
                return tmp;
	}

	
	ConstRowIterator& operator++ ()
	{
		++_idx;
		make_row ();
		return *this;
	}

	ConstRowIterator  operator++ (int)
	{
		ConstRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	ConstRowIterator operator+ (int i) const
		{ return ConstRowIterator (*_M, _idx + i); }

	difference_type operator- (const ConstRowIterator &c) const
		{ return c._idx - _idx; }

	ConstRowIterator& operator += (int i)
	{
		_idx += i;
		make_row ();
		return *this;
	}

	ConstRow operator[] (int i) const
		{ return ConstRow (_M->_rep->rows[_idx + i], _M->_rep->rows[_idx + i] + _M->_rep->width, _M->_rep->ncols); }

	const ConstRow* operator-> () const
		{ return &_row; }

	const ConstRow& operator* () const
		{ return _row; }

	bool operator == (const ConstRowIterator& c) const
		{ return _row.wordBegin () == c._row.wordBegin (); }
    
	bool operator!= (const ConstRowIterator& c) const
		{ return _row.wordBegin () != c._row.wordBegin (); }

    private:
	friend class M4RIMatrix::RowIterator;

	inline void make_row ()
	{
		_row = ConstRow (_M->_rep->rows[_idx], _M->_rep->rows[_idx] + _M->_rep->width, _M->_rep->ncols);
	}

	const M4RIMatrix *_M;
	size_t _idx;
	ConstRow _row;
};

class M4RIMatrix::RowIterator
{
    public:
	typedef int difference_type;

	RowIterator (M4RIMatrix &M, size_t idx)
		: _M (&M), _idx (idx) { make_row (); }
    
	RowIterator () {}
    
	RowIterator (const RowIterator& colp)
		: _M (colp._M), _idx (colp._idx), _row (colp._row) {}
    
	RowIterator& operator = (const RowIterator& colp)
	{
		_M = colp._M;
		_idx = colp._idx;
		_row = colp._row;
		return *this;
	}

	RowIterator& operator --()
	{
		--_idx;
		make_row ();
		return *this;
	}

	RowIterator  operator-- (int)
        {
                RowIterator tmp (*this);
                --*this;
                return tmp;
	}

	
	RowIterator& operator++ ()
	{
		++_idx;
		make_row ();
		return *this;
	}

	RowIterator  operator++ (int)
	{
		RowIterator tmp (*this);
		++*this;
		return tmp;
	}

	RowIterator operator+ (int i) const
		{ return RowIterator (*_M, _idx + i); }

	difference_type operator- (const RowIterator &c) const
		{ return c._idx - _idx; }

	RowIterator& operator += (int i)
	{
		_idx += i;
		make_row ();
		return *this;
	}

	Row operator[] (int i) const
		{ return Row (_M->_rep->rows[_idx + i], _M->_rep->rows[_idx + i] + _M->_rep->width, _M->_rep->ncols); }

	Row* operator-> ()
		{ return &_row; }

	Row& operator* ()
		{ return _row; }

	bool operator == (const RowIterator& c) const
		{ return _row.wordBegin () == c._row.wordBegin (); }
    
	bool operator == (const ConstRowIterator& c) const
		{ return _row.wordBegin () == c._row.wordBegin (); }
    
	bool operator!= (const RowIterator& c) const
		{ return _row.wordBegin () != c._row.wordBegin (); }

	bool operator!= (const ConstRowIterator& c) const
		{ return _row.wordBegin () != c._row.wordBegin (); }

	operator ConstRowIterator ()
		{ return ConstRowIterator (*_M, _idx); }

    private:

	inline void make_row ()
	{
		_row = Row (_M->_rep->rows[_idx], _M->_rep->rows[_idx] + _M->_rep->width, _M->_rep->ncols);
	}

	M4RIMatrix *_M;
	size_t _idx;
	Row _row;
};

inline M4RIMatrix::RowIterator M4RIMatrix::rowBegin ()
{
	return RowIterator (*this, 0);
}

inline M4RIMatrix::RowIterator M4RIMatrix::rowEnd ()
{
	return RowIterator (*this, rowdim ());
}

inline M4RIMatrix::ConstRowIterator M4RIMatrix::rowBegin () const
{
	return ConstRowIterator (*this, 0);
}

inline M4RIMatrix::ConstRowIterator M4RIMatrix::rowEnd () const
{
	return ConstRowIterator (*this, rowdim ());
}

template <class Field>
std::istream& M4RIMatrix::read (std::istream &file, const Field& field)
{
	size_t i, j;

	for (i = 0; i < rowdim (); ++i) {
		for (j = 0; j < coldim (); ++j) {
			file.ignore (1);
			field.read (file, refEntry (i, j));
		}
	}

	return file;
}

template <class Field>
std::ostream& M4RIMatrix::write (std::ostream &os, const Field &F) const
{
	ConstRowIterator p;

	for (p = rowBegin (); p != rowEnd (); ++p) {
		typename ConstRow::const_iterator pe;

		os << "  [ ";

		for (pe = p->begin (); pe != p->end (); ++pe) {
			if (F.isZero (*pe))
				os << '.';
			else
				os << '1';
			os << " ";
		}

		os << "]" << std::endl;
	}

	return os;
}

template <>
class SubvectorFactory<M4RIMatrix>
{
    public:
	typedef BitSubvector<BitVectorIterator<word *, const word *, BigEndian<word> >, BitVectorConstIterator<const word *, BigEndian<word> > > RowSubvector;
	typedef BitSubvector<BitVectorConstIterator<const word *, BigEndian<word> > > ConstRowSubvector;

	RowSubvector MakeRowSubvector (SubmatrixBase<M4RIMatrix> &M, M4RIMatrix::RowIterator &pos)
		{ return RowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

	ConstRowSubvector MakeConstRowSubvector (const SubmatrixBase<M4RIMatrix> &M, M4RIMatrix::ConstRowIterator &pos)
		{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

template <>
class SubvectorFactory<const M4RIMatrix>
{
    public:
	typedef BitSubvector<BitVectorConstIterator<const word *, BigEndian<word> > > RowSubvector;
	typedef BitSubvector<BitVectorConstIterator<const word *, BigEndian<word> > > ConstRowSubvector;

	ConstRowSubvector MakeConstRowSubvector (const SubmatrixBase<const M4RIMatrix> &M, M4RIMatrix::ConstRowIterator &pos)
		{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

}

#endif // __LINBOX_MATRIX_M4RI_MATRIX_INL
