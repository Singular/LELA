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
#include "linbox/matrix/submatrix.h"
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

	template <class It, class CIt>
	DenseMatrixRowIterator (const DenseMatrixRowIterator<It, CIt> &colp)
		: _row (colp._row), _dis (colp._dis) {}

	DenseMatrixRowIterator &operator = (const DenseMatrixRowIterator &colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}

	template <class It, class CIt>
	DenseMatrixRowIterator &operator = (const DenseMatrixRowIterator<It, CIt> &colp)
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

	difference_type operator - (const DenseMatrixRowIterator &c) const
		{ return (_row.begin () - c._row.begin ()) / _dis; }

	template <class It, class CIt>
	difference_type operator - (const DenseMatrixRowIterator<It, CIt> &c) const
		{ return (_row.begin () - c._row.begin ()) / _dis; }

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
		: _col (rowp._col), _stride (rowp._stride) {}

	template <class It, class CIt>
	DenseMatrixColIterator (const DenseMatrixColIterator<It, CIt> &rowp)
		: _col (rowp._col), _stride (rowp._stride) {}

	DenseMatrixColIterator &operator = (const DenseMatrixColIterator &rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}

	template <class It, class CIt>
	DenseMatrixColIterator &operator = (const DenseMatrixColIterator<It, CIt> &rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
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

	difference_type operator - (const DenseMatrixColIterator &c) const
		{ return _col.begin () - c._col.begin (); }

	template <class It, class CIt>
	difference_type operator - (const DenseMatrixColIterator<It, CIt> &c) const
		{ return _col.begin () - c._col.begin (); }

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
DenseMatrix<_Element>::DenseMatrix (VectorStream<DenseMatrix<_Element>::Row> &vs)
	: _rep (vs.size () * vs.dim ()), _rows (vs.size ()), _cols (vs.dim ()), _disp (vs.dim ()), _ptr (&_rep[0])
{
	_rep_begin = _rep.begin ();
	_rep_end = _rep.end ();

	for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
		vs >> *i;
}

/// entry access raw view.  Size m*n vector in C (row major) order.
template <class Element>
typename DenseMatrix<Element>::RawIterator DenseMatrix<Element>::rawBegin ()
	{ return _rep_begin; }  

template <class Element>
typename DenseMatrix<Element>::RawIterator DenseMatrix<Element>::rawEnd ()
	{ return _rep_end; }
  
template <class Element>
typename DenseMatrix<Element>::ConstRawIterator DenseMatrix<Element>::rawBegin () const
	{ return _rep_begin; }  

template <class Element>
typename DenseMatrix<Element>::ConstRawIterator DenseMatrix<Element>::rawEnd () const
	{ return _rep_end; }

template <class Element>
typename DenseMatrix<Element>::RowIterator DenseMatrix<Element>::rowBegin ()
	{ return RowIterator (_rep_begin, _cols, _disp); }

template <class Element>
typename DenseMatrix<Element>::RowIterator DenseMatrix<Element>::rowEnd ()
	{ return RowIterator (_rep_begin + _rows * _disp, _cols, _disp); }
  
template <class Element>
typename DenseMatrix<Element>::ConstRowIterator DenseMatrix<Element>::rowBegin () const
	{ return ConstRowIterator (_rep_begin, _cols, _disp); }  

template <class Element>
typename DenseMatrix<Element>::ConstRowIterator DenseMatrix<Element>::rowEnd () const
	{return ConstRowIterator (_rep_begin + _rows * _disp, _cols, _disp); }
  
template <class Element>
typename DenseMatrix<Element>::ColIterator DenseMatrix<Element>::colBegin ()
	{ return  typename DenseMatrix<Element>::ColIterator (_rep_begin, _disp, _rows); }

template <class Element>
typename DenseMatrix<Element>::ColIterator DenseMatrix<Element>::colEnd ()
	{ return  typename DenseMatrix<Element>::ColIterator (_rep_begin + _cols, _disp, _rows); }
  
template <class Element>
typename DenseMatrix<Element>::ConstColIterator DenseMatrix<Element>::colBegin () const
	{ return  typename DenseMatrix<Element>::ConstColIterator (_rep_begin, _disp, _rows); }
  
template <class Element>
typename DenseMatrix<Element>::ConstColIterator DenseMatrix<Element>::colEnd () const
	{ return  typename DenseMatrix<Element>::ConstColIterator (_rep_begin + _cols, _disp, _rows); }

template <class Element>
class SubvectorFactory<DenseMatrix<Element> >
{
    public:
	typedef typename DenseMatrix<Element>::Row RowSubvector;
	typedef typename DenseMatrix<Element>::ConstRow ConstRowSubvector;
	typedef typename DenseMatrix<Element>::Col ColSubvector;
	typedef typename DenseMatrix<Element>::ConstCol ConstColSubvector;

	RowSubvector MakeRowSubvector (Submatrix<DenseMatrix<Element> > &M,
				       typename DenseMatrix<Element>::RowIterator &pos)
		{ return RowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<DenseMatrix<Element> > &M,
						 typename DenseMatrix<Element>::ConstRowIterator &pos)
		{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }

	ColSubvector MakeColSubvector (Submatrix<DenseMatrix<Element> > &M,
				       typename DenseMatrix<Element>::ColIterator &pos)
		{ return ColSubvector (Subiterator<typename DenseMatrix<Element>::Rep::iterator> (pos->begin ().operator-> (), M.parent ().coldim ()) + M.startRow (),
				       Subiterator<typename DenseMatrix<Element>::Rep::iterator> (pos->begin ().operator-> (), M.parent ().coldim ()) + (M.startRow () + M.rowdim ())); }
	ConstColSubvector MakeConstColSubvector (const Submatrix<DenseMatrix<Element> > &M,
						 typename DenseMatrix<Element>::ConstColIterator &pos)
		{ return ConstColSubvector (Subiterator<typename DenseMatrix<Element>::Rep::const_iterator> (pos->begin ().operator-> (), M.parent ().coldim ()) + M.startRow (),
					    Subiterator<typename DenseMatrix<Element>::Rep::const_iterator> (pos->begin ().operator-> (), M.parent ().coldim ()) + (M.startRow () + M.rowdim ())); }
};

template <class Element>
class SubvectorFactory<const DenseMatrix<Element> >
{
    public:
	typedef typename DenseMatrix<Element>::ConstRow RowSubvector;
	typedef typename DenseMatrix<Element>::ConstRow ConstRowSubvector;
	typedef typename DenseMatrix<Element>::ConstCol ColSubvector;
	typedef typename DenseMatrix<Element>::ConstCol ConstColSubvector;

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<const DenseMatrix<Element> > &M,
						 typename DenseMatrix<Element>::ConstRowIterator &pos)
		{ return ConstRowSubvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
	ConstColSubvector MakeConstColSubvector (const Submatrix<const DenseMatrix<Element> > &M,
						 typename DenseMatrix<Element>::ConstColIterator &pos)
		{ return ConstColSubvector (Subiterator<typename DenseMatrix<Element>::Rep::const_iterator> (pos->begin ().operator-> (), M.parent ().coldim ()) + M.startRow (),
					    Subiterator<typename DenseMatrix<Element>::Rep::const_iterator> (pos->begin ().operator-> (), M.parent ().coldim ()) + (M.startRow () + M.rowdim ())); }
};

} // namespace LinBox

#endif // __LINBOX_matrix_dense_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
