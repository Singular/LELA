/* linbox/matrix/dense.inl
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
 * Move from blackbox/dense-base.inl to matrix/dense.inl
 * --------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Split out container/iterator functionality into DenseMatrixBase
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_matrix_dense_INL
#define __LINBOX_matrix_dense_INL

#include <iostream>
#include <cmath>

#include "linbox/matrix/dense.h"
#include "linbox/util/debug.h"

namespace LinBox
{

template <class Iterator, class ConstIterator>
class DenseMatrixRowIterator
{
    public:
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

	typedef Subvector<Iterator, ConstIterator> Row;
	typedef Subvector<ConstIterator, ConstIterator> ConstRow;

	DenseMatrixRowIterator (const Iterator &p, size_t len, size_t d)
		: _row (p, p + len), _dis (d) {}

	DenseMatrixRowIterator () {}

	DenseMatrixRowIterator (const DenseMatrixRowIterator &colp)
		: _row (colp._row), _dis (colp._dis) {}

	DenseMatrixRowIterator &operator = (const DenseMatrixRowIterator &colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}
    
	DenseMatrixRowIterator &operator ++ ()
	{
		_row = Row (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}
    
	DenseMatrixRowIterator operator ++ (int)
	{
		DenseMatrixRowIterator tmp (*this);
		++*this;
		return tmp;
	}
    
        DenseMatrixRowIterator &operator -- ()
        {
                _row = Row (_row.begin () - _dis, _row.end () - _dis);
                return *this;
        }

        DenseMatrixRowIterator operator -- (int)
        {
                DenseMatrixRowIterator tmp (*this);
                --*this;
                return tmp;
        }

	DenseMatrixRowIterator operator + (int i)
		{ return DenseMatrixRowIterator (_row.begin () + _dis * i, _row.size (), _dis); }

	DenseMatrixRowIterator &operator += (int i)
	{
		_row = Row (_row.begin () + _dis * i, _row.end () + _dis * i);
		return *this;
	}

	Row operator [] (int i) const
		{ return Row (const_cast<Row&> (_row).begin () + _dis * i,
			      const_cast<Row&> (_row).end () + _dis * i); }

	Row *operator -> ()
		{ return &_row; }

	Row &operator * ()
		{ return _row; }

	bool operator == (const DenseMatrixRowIterator &c) const
		{ return (_row.begin () == c._row.begin ()) && (_row.end () == c._row.end ()) && (_dis == c._dis); }

	template <class It, class CIt>
	bool operator == (const DenseMatrixRowIterator<It, CIt> &c) const
		{ return (_row.begin () == c._row.begin ()) && (_row.end () == c._row.end ()) && (_dis == c._dis); }

	bool operator != (const DenseMatrixRowIterator &c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

	template <class It, class CIt>
	bool operator != (const DenseMatrixRowIterator<It, CIt> &c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

    private:
	Row _row;
	size_t _dis;

	template <class It, class CIt>
	friend class DenseMatrixRowIterator;
};

template <class Iterator, class ConstIterator>
class DenseMatrixColIterator
{
    public:
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

	typedef Subvector<Subiterator<Iterator>, Subiterator<ConstIterator> > Col;
	typedef Subvector<Subiterator<ConstIterator> > ConstCol;

	DenseMatrixColIterator (Iterator p, size_t stride, size_t len)
		: _col (Subiterator<Iterator> (p, stride),
			Subiterator<Iterator> (p + len * stride, stride)), _stride (stride)
	{}
    
	DenseMatrixColIterator () {}
    
	DenseMatrixColIterator (const DenseMatrixColIterator &rowp)
		:_col (rowp._col){}
    
	DenseMatrixColIterator &operator = (const DenseMatrixColIterator &rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}
    
	const DenseMatrixColIterator &operator = (const DenseMatrixColIterator &rowp) const
	{
		const_cast<DenseMatrixColIterator*> (this)->_col = rowp._col;
		return *this;
	}
    
	DenseMatrixColIterator &operator++ ()
	{
		_col = Col (Subiterator<Iterator> (_col.begin ().operator-> () + 1, _stride),
			    Subiterator<Iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}
    
	DenseMatrixColIterator operator++ (int)
	{
		Col tmp (_col);
		this->operator++ ();
		return tmp;
	}
        
	DenseMatrixColIterator operator + (int i)
		{ return DenseMatrixColIterator (_col.begin ().operator-> () + i, _stride, _col.size ()); }

	DenseMatrixColIterator &operator += (int i)
	{
		_col = Col (Subiterator<Iterator> (_col.begin ().operator-> () + i, _stride),
			    Subiterator<Iterator> (_col.end ().operator-> () + i, _stride));
		return *this;
	}

	Col operator [] (int i) const
		{ return Col (Subiterator<Iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride), 
			      Subiterator<Iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride)); }

	Col *operator -> ()
		{ return &_col; }

	Col &operator * ()
		{ return _col; }

	bool operator == (const DenseMatrixColIterator &c) const
		{ return (_col.begin () == c._col.begin ()) && (_col.end () == c._col.end ()); }

	template <class It, class CIt>
	bool operator == (const DenseMatrixColIterator<It, CIt> &c) const
		{ return (_col.begin () == c._col.begin ()) && (_col.end () == c._col.end ()); }

	bool operator != (const DenseMatrixColIterator &c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }

	template <class It, class CIt>
	bool operator != (const DenseMatrixColIterator<It, CIt> &c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }

    private:

	Col _col;
	size_t _stride;

	template <class It, class CIt>
	friend class DenseMatrixColIterator;
};

template <class _Element>
template <class Field>
DenseMatrix<_Element>::DenseMatrix (MatrixStream<Field> &ms)
	:_rep (0), _rows (0), _cols (0), _ptr (NULL)
{
	if (!ms.getArray (_rep) || !ms.getRows (_rows) || !ms.getColumns (_cols) )
		throw ms.reportError (__FUNCTION__, __LINE__);
	_ptr = &_rep[0];
}

template <class _Element>
DenseMatrix<_Element>::DenseMatrix (VectorStream<DenseMatrix<_Element>::Row> &vs)
	: _rep (vs.size () * vs.dim ()), _rows (vs.size ()), _cols (vs.dim ()), _ptr (&_rep[0])
{
	for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
		vs >> *i;
}

/// entry access raw view.  Size m*n vector in C (row major) order.
template <class Element>
typename DenseMatrix<Element>::RawIterator DenseMatrix<Element>::rawBegin ()
	{ return _rep.begin (); }  

template <class Element>
typename DenseMatrix<Element>::RawIterator DenseMatrix<Element>::rawEnd ()
	{ return _rep.end (); }
  
template <class Element>
typename DenseMatrix<Element>::ConstRawIterator DenseMatrix<Element>::rawBegin () const
	{ return _rep.begin (); }  

template <class Element>
typename DenseMatrix<Element>::ConstRawIterator DenseMatrix<Element>::rawEnd () const
	{ return _rep.end (); }

template <class Element>
typename DenseMatrix<Element>::RowIterator DenseMatrix<Element>::rowBegin ()
	{ return RowIterator (_rep.begin (), _cols, _cols); }

template <class Element>
typename DenseMatrix<Element>::RowIterator DenseMatrix<Element>::rowEnd ()
	{ return RowIterator (_rep.end (), _cols, _cols); }
  
template <class Element>
typename DenseMatrix<Element>::ConstRowIterator DenseMatrix<Element>::rowBegin () const
	{ return ConstRowIterator (_rep.begin (), _cols, _cols); }  

template <class Element>
typename DenseMatrix<Element>::ConstRowIterator DenseMatrix<Element>::rowEnd () const
	{return ConstRowIterator (_rep.end (), _cols, _cols); }
  
template <class Element>
typename DenseMatrix<Element>::ColIterator DenseMatrix<Element>::colBegin ()
	{ return  typename DenseMatrix<Element>::ColIterator (_rep.begin (), _cols, _rows); }

template <class Element>
typename DenseMatrix<Element>::ColIterator DenseMatrix<Element>::colEnd ()
	{ return  typename DenseMatrix<Element>::ColIterator (_rep.begin ()+_cols, _cols, _rows); }
  
template <class Element>
typename DenseMatrix<Element>::ConstColIterator DenseMatrix<Element>::colBegin () const
	{ return  typename DenseMatrix<Element>::ConstColIterator (_rep.begin (), _cols, _rows); }
  
template <class Element>
typename DenseMatrix<Element>::ConstColIterator DenseMatrix<Element>::colEnd () const
	{ return  typename DenseMatrix<Element>::ConstColIterator (_rep.begin ()+_cols, _cols, _rows); }

} // namespace LinBox

#endif // __LINBOX_matrix_dense_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
