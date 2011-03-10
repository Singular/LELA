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

template <class Element>
class DenseMatrix<Element>::ConstRowIterator
{
    public:
	typedef typename std::iterator_traits<typename Rep::const_iterator>::difference_type difference_type;

	ConstRowIterator (const typename Rep::const_iterator& p, size_t len, size_t d)
		: _row (p, p + len), _dis (d) {}
    
	ConstRowIterator () {}
    
	ConstRowIterator (const ConstRowIterator& colp)
		: _row (colp._row), _dis (colp._dis) {}
    
	ConstRowIterator& operator = (const ConstRowIterator& colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}

	ConstRowIterator& operator --()
	{
		_row = ConstRow (_row.begin () - _dis, _row.end () - _dis);
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
		_row = ConstRow (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}

	ConstRowIterator  operator++ (int)
	{
		ConstRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	ConstRowIterator operator+ (int i)
		{ return ConstRowIterator (_row.begin () + _dis * i, _row.size (), _dis); }

	ConstRowIterator& operator += (int i)
	{
		_row = ConstRow (_row.begin () + _dis * i, _row.end () + _dis * i);
		return *this;
	}

	ConstRow operator[] (int i) const
		{ return ConstRow (_row.begin () + _dis * i, _row.end () + _dis * i); }

	ConstRow* operator-> ()
		{ return &_row; }

	ConstRow& operator* ()
		{ return _row; }
    
	bool operator!= (const ConstRowIterator& c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

    private:
	ConstRow _row;
	size_t _dis;
};

template <class Element>
class DenseMatrix<Element>::RowIterator
{
    public:
	typedef typename std::iterator_traits<typename Rep::iterator>::difference_type difference_type;

	RowIterator (const typename Rep::iterator& p, size_t len, size_t d)
		: _row (p, p + len), _dis (d){}

	RowIterator () {}

	RowIterator (const RowIterator& colp)
		: _row (colp._row), _dis (colp._dis) {}
    
	RowIterator& operator = (const RowIterator& colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}
    
	RowIterator& operator ++ ()
	{
		_row = Row (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}
    
	RowIterator  operator ++ (int)
	{
		RowIterator tmp (*this);
		++*this;
		return tmp;
	}
    
        RowIterator& operator -- ()
        {
                _row = Row (_row.begin () - _dis, _row.end () - _dis);
                return *this;
        }

        RowIterator  operator -- (int)
        {
                RowIterator tmp (*this);
                --*this;
                return tmp;
        }

	RowIterator operator + (int i)
		{ return RowIterator (_row.begin () + _dis * i, _row.size (), _dis); }

	RowIterator& operator += (int i)
	{
		_row = Row (_row.begin () + _dis * i, _row.end () + _dis * i);
		return *this;
	}

	Row operator[] (int i) const
		{ return Row (const_cast<Row&> (_row).begin () + _dis * i,
			      const_cast<Row&> (_row).end () + _dis * i); }

	Row* operator-> ()
		{ return &_row; }
    
	Row& operator* ()
		{ return _row; }
 
	bool operator!= (const RowIterator& c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

	operator ConstRowIterator ()
		{ return ConstRowIterator (_row.begin (), _row.size (), _dis); }

    private:
	Row _row;
	size_t _dis;
};

template <class Element>
class DenseMatrix<Element>::ConstColIterator
{
    public:
	typedef typename std::iterator_traits<typename Rep::const_iterator>::difference_type difference_type;

	ConstColIterator (typename Rep::const_iterator p, size_t stride, size_t len)
		: _col (Subiterator<typename Rep::const_iterator> (p, stride),
			Subiterator<typename Rep::const_iterator> (p + len * stride, stride)), _stride (stride)
	{}
    
	ConstColIterator (const ConstCol& col, size_t stride)
		:_col (col), _stride (stride){}

	ConstColIterator () {}
    
	ConstColIterator (const ConstColIterator& rowp)
		:_col (rowp._col){}

	ConstColIterator& operator= (const ConstColIterator& rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}

	ConstColIterator& operator++ ()
	{
		_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + 1, _stride),
				 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}

	ConstColIterator  operator++ (int)
	{
		ConstColIterator old(*this);
		this->operator++ ();
		return old;
	}

	ConstColIterator operator + (int i)
		{ return ConstColIterator (_col.begin ().operator-> () + i, _stride, _col.size ()); }

	ConstColIterator& operator += (int i)
	{
		_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride),
				 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride));
		return *this;
	}

	ConstCol operator[] (int i) const
		{ return ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride), 
				   Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride)); }

	ConstCol* operator-> ()
		{ return &_col; }
 
	ConstCol& operator* ()
		{ return _col; }
    
	bool operator!= (const ConstColIterator& c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }
    
    private:
	ConstCol _col;
	size_t _stride;
};

template <class Element>
class DenseMatrix<Element>::ColIterator
{
    public:
	typedef typename std::iterator_traits<typename Rep::iterator>::difference_type difference_type;

	ColIterator (typename Rep::iterator p, size_t stride, size_t len)
		: _col (Subiterator<typename Rep::iterator> (p, stride),
			Subiterator<typename Rep::iterator> (p + len * stride, stride)), _stride (stride)
	{}
    
	ColIterator () {}
    
	ColIterator (const ColIterator& rowp)
		:_col (rowp._col){}
    
	ColIterator& operator= (const ColIterator& rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}
    
	const ColIterator& operator= (const ColIterator& rowp) const
	{
		const_cast<ColIterator*> (this)->_col = rowp._col;
		return *this;
	}
    
	ColIterator& operator++ ()
	{
		_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + 1, _stride),
			    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}
    
	ColIterator  operator++ (int)
	{
		Col tmp (_col);
		this->operator++ ();
		return tmp;
	}
        
	ColIterator operator + (int i)
		{ return ColIterator (_col.begin ().operator-> () + i, _stride, _col.size ()); }

	ColIterator& operator += (int i)
	{
		_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + i, _stride),
			    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + i, _stride));
		return *this;
	}

	Col operator[] (int i) const
		{ return Col (Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride), 
			      Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride)); }

	Col* operator-> ()
		{ return &_col; }

	Col& operator* ()
		{ return _col; }

	bool operator!= (const ColIterator& c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }
   
	operator ConstColIterator ()
	{
            return ConstColIterator (reinterpret_cast<ConstCol&> (_col) , _stride);
	}

    private:

	Col _col;
	size_t _stride;
};

template <class _Element>
template <class Field>
DenseMatrix<_Element>::DenseMatrix( MatrixStream<Field>& ms )
	:_rep(0), _rows(0), _cols(0), _ptr(NULL)
{
	if( !ms.getArray(_rep) || !ms.getRows(_rows) || !ms.getColumns(_cols) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	_ptr = &_rep[0];
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
