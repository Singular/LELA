/* lela/matrix/dense.h
 * Copyright 2001 B. David Saunders,
 *           2001-2002, 2011 Bradford Hovinen,
 *           2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@gmail.com>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 *
 * --------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_matrix_dense_H
#define __LELA_matrix_dense_H

#include <iostream>
#include <vector>
#include <fstream>

#include "lela/lela-config.h"
#include "lela/util/debug.h"
#include "lela/vector/subiterator.h"
#include "lela/vector/subvector.h"
#include "lela/vector/stream.h"
#include "lela/matrix/raw-iterator.h"
#include "lela/matrix/submatrix.h"
#include "lela/matrix/traits.h"

namespace LELA
{

/** Dense matrix
 *
 * See @ref MatrixArchetype for documentation on the interface
 *
\ingroup matrix
 */
  
template <class _Element>
class DenseMatrix
{
	template <class Iterator, class ConstIterator>
	class RowIteratorPT
	{
	public:
		typedef Subvector<Iterator, ConstIterator> Row;
		typedef Subvector<ConstIterator, ConstIterator> ConstRow;

		typedef std::random_access_iterator_tag iterator_category;
		typedef Row value_type;
		typedef Row &reference;
		typedef ConstRow &const_reference;
		typedef Row *pointer;
		typedef ConstRow *const_pointer;
		typedef size_t size_type;
		typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

		RowIteratorPT (const Iterator &p, size_t len, size_t d)
			: _row (p, p + len), _dis (d) {}

		RowIteratorPT () : _dis (0) {}

		RowIteratorPT (const RowIteratorPT &colp)
			: _row (colp._row), _dis (colp._dis) {}

		template <class It, class CIt>
		RowIteratorPT (const RowIteratorPT<It, CIt> &colp)
			: _row (colp._row), _dis (colp._dis) {}

		RowIteratorPT &operator = (const RowIteratorPT &colp)
		{
			_row = colp._row;
			_dis = colp._dis;
			return *this;
		}

		template <class It, class CIt>
		RowIteratorPT &operator = (const RowIteratorPT<It, CIt> &colp)
		{
			_row = colp._row;
			_dis = colp._dis;
			return *this;
		}
    
		RowIteratorPT &operator ++ ()
		{
			_row = Row (_row.begin () + _dis, _row.end () + _dis);
			return *this;
		}
    
		RowIteratorPT operator ++ (int)
		{
			RowIteratorPT tmp (*this);
			++*this;
			return tmp;
		}
    
		RowIteratorPT &operator -- ()
		{
			_row = Row (_row.begin () - _dis, _row.end () - _dis);
			return *this;
		}

		RowIteratorPT operator -- (int)
		{
			RowIteratorPT tmp (*this);
			--*this;
			return tmp;
		}

		RowIteratorPT operator + (int i)
			{ return RowIteratorPT (_row.begin () + _dis * i, _row.size (), _dis); }

		RowIteratorPT &operator += (int i)
		{
			_row = Row (_row.begin () + _dis * i, _row.end () + _dis * i);
			return *this;
		}

		difference_type operator - (const RowIteratorPT &c) const
			{ return (_row.begin () - c._row.begin ()) / _dis; }

		template <class It, class CIt>
		difference_type operator - (const RowIteratorPT<It, CIt> &c) const
			{ return (_row.begin () - c._row.begin ()) / _dis; }

		Row operator [] (int i) const
			{ return Row (const_cast<Row&> (_row).begin () + _dis * i,
				      const_cast<Row&> (_row).end () + _dis * i); }

		Row *operator -> ()
			{ return &_row; }

		Row &operator * ()
			{ return _row; }

		bool operator == (const RowIteratorPT &c) const
			{ return (_row.begin () == c._row.begin ()) && (_row.end () == c._row.end ()) && (_dis == c._dis); }

		template <class It, class CIt>
		bool operator == (const RowIteratorPT<It, CIt> &c) const
			{ return (_row.begin () == c._row.begin ()) && (_row.end () == c._row.end ()) && (_dis == c._dis); }

		bool operator != (const RowIteratorPT &c) const
			{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

		template <class It, class CIt>
		bool operator != (const RowIteratorPT<It, CIt> &c) const
			{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

	private:
		Row _row;
		size_t _dis;

		template <class It, class CIt>
		friend class RowIteratorPT;
	};

	template <class Iterator, class ConstIterator>
	class ColIteratorPT
	{
	public:
		typedef Subvector<Subiterator<Iterator>, Subiterator<ConstIterator> > Col;
		typedef Subvector<Subiterator<ConstIterator> > ConstCol;

		typedef std::random_access_iterator_tag iterator_category;
		typedef Col value_type;
		typedef Col &reference;
		typedef ConstCol &const_reference;
		typedef Col *pointer;
		typedef ConstCol *const_pointer;
		typedef size_t size_type;
		typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

		ColIteratorPT (Iterator p, size_t stride, size_t len)
			: _col (Subiterator<Iterator> (p, stride),
				Subiterator<Iterator> (p + len * stride, stride)), _stride (stride)
			{}
    
		ColIteratorPT () {}
    
		ColIteratorPT (const ColIteratorPT &rowp)
			: _col (rowp._col), _stride (rowp._stride) {}

		template <class It, class CIt>
		ColIteratorPT (const ColIteratorPT<It, CIt> &rowp)
			: _col (rowp._col), _stride (rowp._stride) {}

		ColIteratorPT &operator = (const ColIteratorPT &rowp)
		{
			_col = rowp._col;
			_stride = rowp._stride;
			return *this;
		}

		template <class It, class CIt>
		ColIteratorPT &operator = (const ColIteratorPT<It, CIt> &rowp)
		{
			_col = rowp._col;
			_stride = rowp._stride;
			return *this;
		}
    
		ColIteratorPT &operator++ ()
		{
			_col = Col (Subiterator<Iterator> (_col.begin ().operator-> () + 1, _stride),
				    Subiterator<Iterator> (_col.end ().operator-> () + 1, _stride));
			return *this;
		}
    
		ColIteratorPT operator++ (int)
		{
			Col tmp (_col);
			this->operator++ ();
			return tmp;
		}
        
		ColIteratorPT operator + (int i)
			{ return ColIteratorPT (_col.begin ().operator-> () + i, _stride, _col.size ()); }

		ColIteratorPT &operator += (int i)
		{
			_col = Col (Subiterator<Iterator> (_col.begin ().operator-> () + i, _stride),
				    Subiterator<Iterator> (_col.end ().operator-> () + i, _stride));
			return *this;
		}

		difference_type operator - (const ColIteratorPT &c) const
			{ return _col.begin () - c._col.begin (); }

		template <class It, class CIt>
		difference_type operator - (const ColIteratorPT<It, CIt> &c) const
			{ return _col.begin () - c._col.begin (); }

		Col operator [] (int i) const
			{ return Col (Subiterator<Iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride), 
				      Subiterator<Iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride)); }

		Col *operator -> ()
			{ return &_col; }

		Col &operator * ()
			{ return _col; }

		bool operator == (const ColIteratorPT &c) const
			{ return (_col.begin () == c._col.begin ()) && (_col.end () == c._col.end ()); }

		template <class It, class CIt>
		bool operator == (const ColIteratorPT<It, CIt> &c) const
			{ return (_col.begin () == c._col.begin ()) && (_col.end () == c._col.end ()); }

		bool operator != (const ColIteratorPT &c) const
			{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }

		template <class It, class CIt>
		bool operator != (const ColIteratorPT<It, CIt> &c) const
			{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }

	private:

		Col _col;
		size_t _stride;

		template <class It, class CIt>
		friend class ColIteratorPT;
	};

    public:

	/// @name @ref MatrixArchetype interface
	//@{

	typedef _Element Element;
	typedef typename RawVector<Element>::Dense Rep;
	typedef MatrixIteratorTypes::RowCol IteratorType;
	typedef MatrixStorageTypes::Dense StorageType;
        typedef DenseMatrix<_Element> Self_t;

	typedef Self_t SubmatrixType;
	typedef const Self_t ConstSubmatrixType;
	typedef Self_t AlignedSubmatrixType;
	typedef const Self_t ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = 1;

	typedef Self_t ContainerType;

	typedef Subvector<typename Rep::iterator, typename Rep::const_iterator> Row;
	typedef Subvector<typename Rep::const_iterator> ConstRow;  

	typedef Subvector<Subiterator<typename Rep::iterator> > Col;
	typedef Subvector<Subiterator<typename Rep::const_iterator> > ConstCol;
	typedef Col Column;
	typedef ConstCol ConstColumn;

	typedef RowIteratorPT<typename Rep::iterator, typename Rep::const_iterator> RowIterator;
	typedef RowIteratorPT<typename Rep::const_iterator, typename Rep::const_iterator> ConstRowIterator;

	typedef ColIteratorPT<typename Rep::iterator, typename Rep::const_iterator> ColIterator;
	typedef ColIteratorPT<typename Rep::const_iterator, typename Rep::const_iterator> ConstColIterator;

	template<typename _Tp1>
        struct rebind
        { 
		typedef DenseMatrix<typename _Tp1::Element> other; 
        };

	DenseMatrix ()
		: _rows (0), _cols (0), _disp (0), _start_row (0), _start_col (0)
	{}

	DenseMatrix (size_t m, size_t n)
		: _rep (m * n), _rows (m), _cols (n), _disp (n), _start_row (0), _start_col (0), _ptr (&_rep[0])
	{
		_rep_begin = _rep.begin ();
		_rep_end = _rep.end ();
	}

	DenseMatrix (const DenseMatrix &M)
		: _rep (M._rep), _rep_begin (M._rep_begin), _rep_end (M._rep_end), _rows (M._rows), _cols (M._cols), _disp (M._disp),
		  _start_row (M._start_row), _start_col (M._start_col)
	{
		if (!_rep.empty ()) {
			_rep_begin = _rep.begin ();
			_rep_end = _rep.end ();
			_ptr = &*_rep_begin;
		}
	}

	DenseMatrix (VectorStream<Row> &vs)
		: _rep (vs.size () * vs.dim ()), _rows (vs.size ()), _cols (vs.dim ()), _disp (vs.dim ()), _start_row (0), _start_col (0), _ptr (&_rep[0])
	{
		_rep_begin = _rep.begin ();
		_rep_end = _rep.end ();

		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	DenseMatrix &operator = (const DenseMatrix &M)
	{
		_rep = M._rep;
		_rep_begin = M._rep_begin;
		_rep_end = M._rep_end;
		_rows = M._rows;
		_cols = M._cols;
		_disp = M._disp;
		_start_row = M._start_row;
		_start_col = M._start_col;
		_ptr  = &*_rep_begin;

		if (!_rep.empty ()) {
			_rep_begin = _rep.begin ();
			_rep_end = _rep.end ();
		}

		return (*this);
	}

	size_t rowdim () const { return _rows; }
	size_t coldim () const { return _cols; }

	void resize (size_t m, size_t n)
	{
		_rows = m;
		_cols = n;
		_disp = n;
		_rep.resize (m * n);
		_rep_begin = _rep.begin ();
		_rep_end = _rep.end ();
	}

	void setEntry (size_t i, size_t j, const Element &a_ij) { _rep_begin[i * _disp + j] = a_ij; }
	void eraseEntry (size_t i, size_t j) {}
	bool getEntry (Element &x, size_t i, size_t j) const { x = _rep_begin[i * _disp + j]; return true; }

	RowIterator      rowBegin ()       { return RowIterator (_rep_begin, _cols, _disp); }
	RowIterator      rowEnd ()         { return RowIterator (_rep_begin + _rows * _disp, _cols, _disp); }
	ConstRowIterator rowBegin () const { return ConstRowIterator (_rep_begin, _cols, _disp); }  
	ConstRowIterator rowEnd ()   const { return ConstRowIterator (_rep_begin + _rows * _disp, _cols, _disp); }

	ColIterator      colBegin ()       { return ColIterator (_rep_begin, _disp, _rows); }
	ColIterator      colEnd ()         { return ColIterator (_rep_begin + _cols, _disp, _rows); }
	ConstColIterator colBegin () const { return ConstColIterator (_rep_begin, _disp, _rows); }
	ConstColIterator colEnd ()   const { return ConstColIterator (_rep_begin + _cols, _disp, _rows); }

	typedef MatrixRawIterator<RowIterator, VectorRepresentationTypes::Dense> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, VectorRepresentationTypes::Dense> ConstRawIterator;
    
	RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorRepresentationTypes::Dense, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }
    
	Row      operator[] (size_t i)       { return Row (_rep_begin + i * _disp, _rep_begin + i * _disp + _cols); }
	ConstRow operator[] (size_t i) const { return Row (_rep_begin + i * _disp, _rep_begin + i * _disp + _cols); }

	template <class Vector>
	Vector &columnDensity (Vector &v) const { std::fill (v.begin (), v.end (), _rows); return v; }

	//@}

	/// @name Additional interfaces
	//@{

	/** Construct a dense submatrix of the given matrix.
	 * @param  M parent-matrix
	 * @param  beg_row  Beginning row (indexed at 0)
	 * @param  beg_col  Beginning column (indexed at 0)
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (const DenseMatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rep_begin (M._rep_begin + (beg_row * M._disp + beg_col)),
		  _rep_end (M._rep_begin + ((beg_row + m) * M._disp + beg_col + n)),
		  _rows (m), _cols (n), _disp (M._disp), _start_row (M._start_row + beg_row), _start_col (M._start_col + beg_col)
		{ _ptr = &*_rep_begin; }

	/** Get the displacement from one row to the next
	 * @returns Number of words from one row to the next in memory
	 */
	size_t disp () const
		{ return _disp; }

	/** Return the starting row in the parent-matrix */
	size_t startRow () const { return _start_row; }

	/** Return the starting column in the parent-matrix */
	size_t startCol () const { return _start_col; }

	//@}

    protected:

	std::vector<Element> _rep;

	typename std::vector<Element>::iterator _rep_begin, _rep_end;

	size_t   _rows, _cols, _disp;
	size_t   _start_row, _start_col;
	Element *_ptr;
};

} // namespace LELA

// Pull in the specialisation for GF(2)
#include "lela/matrix/dense-zero-one.h"

#endif // __LELA_matrix_dense_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
