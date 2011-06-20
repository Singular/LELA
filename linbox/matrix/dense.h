/* linbox/matrix/dense.h
 * Copyright 2001 B. David Saunders,
 *           2001-2002 Bradford Hovinen,
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
 * See COPYING for license information
 */

#ifndef __LINBOX_matrix_dense_H
#define __LINBOX_matrix_dense_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/raw-iterator.h"
#include "linbox/matrix/submatrix.h"
#include "linbox/matrix/traits.h"

namespace LinBox
{

template <class Element>
struct DenseMatrixTag {};

/** Blackbox dense matrix template. This is a class of dense matrices
 * templatized by the entry type, the Element type of some {@link Rings ring}.
 * The matrix is stored as a one dimensional STL vector of the elements, by rows. 
 * The interface provides for iteration over rows and over columns.
 *
 * The class LinBox::Dense builds on this base.
 *
 * Currently, only dense vectors are supported when doing matrix-vector applies.
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

	typedef _Element Element;
	typedef typename RawVector<Element>::Dense Rep;
	typedef MatrixCategories::RowColMatrixTag MatrixCategory;
        typedef DenseMatrix<_Element> Self_t;
	typedef DenseMatrixTag<Element> Tag;

	typedef Self_t SubmatrixType;
	typedef const Self_t ConstSubmatrixType;

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

	///
	DenseMatrix ()
		: _rows (0), _cols (0), _disp (0), _start_row (0), _start_col (0)
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (size_t m, size_t n)
		: _rep (m * n), _rows (m), _cols (n), _disp (n), _start_row (0), _start_col (0), _ptr (&_rep[0])
	{
		_rep_begin = _rep.begin ();
		_rep_end = _rep.end ();
	}

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

	///
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

	///
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
		_disp = n;
		_rep.resize (m * n);
		_rep_begin = _rep.begin ();
		_rep_end = _rep.end ();
	}

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _rep_begin[i * _disp + j] = a_ij; }

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
		{ x = _rep_begin[i * _disp + j]; return true; }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	RowIterator rowBegin ()
		{ return RowIterator (_rep_begin, _cols, _disp); }
	RowIterator rowEnd ()
		{ return RowIterator (_rep_begin + _rows * _disp, _cols, _disp); }
	ConstRowIterator rowBegin () const
		{ return ConstRowIterator (_rep_begin, _cols, _disp); }  
	ConstRowIterator rowEnd () const
		{ return ConstRowIterator (_rep_begin + _rows * _disp, _cols, _disp); }

	/** @name Row of columns iterator
	 * The row of columns iterator traverses the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */

	ColIterator colBegin ()
		{ return ColIterator (_rep_begin, _disp, _rows); }
	ColIterator colEnd ()
		{ return ColIterator (_rep_begin + _cols, _disp, _rows); }
	ConstColIterator colBegin () const
		{ return ConstColIterator (_rep_begin, _disp, _rows); }
	ConstColIterator colEnd () const
		{ return ConstColIterator (_rep_begin + _cols, _disp, _rows); }

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	typedef typename Rep::iterator RawIterator;
	typedef typename Rep::const_iterator ConstRawIterator;
    
	RawIterator rawBegin ()
		{ return _rep_begin; }  
	RawIterator rawEnd ()
		{ return _rep_end; }
	ConstRawIterator rawBegin () const
		{ return _rep_begin; }  
	ConstRawIterator rawEnd () const
		{ return _rep_end; }

	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::DenseVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }
    
	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 * @param i Row index
	 */
	Row operator[] (size_t i)
		{ return Row (_rep_begin + i * _disp, _rep_begin + i * _disp + _cols); }

	ConstRow operator[] (size_t i) const
		{ return Row (_rep_begin + i * _disp, _rep_begin + i * _disp + _cols); }

	/** Compute column density
	 */

	template <class Vector>
	Vector &columnDensity (Vector &v) const
		{ std::fill (v.begin (), v.end (), _rows); return v; }

	/** Get the displacement from one row to the next
	 * @returns Number of words from one row to the next in memory
	 */
	size_t disp () const
		{ return _disp; }

	/** Return the starting row in the parent-matrix */
	size_t startRow () const { return _start_row; }

	/** Return the starting column in the parent-matrix */
	size_t startCol () const { return _start_col; }

    protected:

	std::vector<Element> _rep;

	typename std::vector<Element>::iterator _rep_begin, _rep_end;

	size_t   _rows, _cols, _disp;
	size_t   _start_row, _start_col;
	Element *_ptr;
};

template <class Element, class Submatrix>
class SubvectorFactory<DenseMatrixTag<Element>, Submatrix, typename DenseMatrix<Element>::RowIterator, DefaultSubvectorFactoryTrait>
{
    public:
	typedef typename DenseMatrix<Element>::RowIterator IteratorType;
	typedef typename DenseMatrix<Element>::Row Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

template <class Element, class Submatrix>
class SubvectorFactory<DenseMatrixTag<Element>, Submatrix, typename DenseMatrix<Element>::ConstRowIterator, DefaultSubvectorFactoryTrait>
{
    public:
	typedef typename DenseMatrix<Element>::ConstRowIterator IteratorType;
	typedef typename DenseMatrix<Element>::ConstRow Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

template <class Element, class Submatrix>
class SubvectorFactory<DenseMatrixTag<Element>, Submatrix, typename DenseMatrix<Element>::ColIterator, DefaultSubvectorFactoryTrait>
{
    public:
	typedef typename DenseMatrix<Element>::ColIterator IteratorType;
	typedef typename DenseMatrix<Element>::Col Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (Subiterator<typename DenseMatrix<Element>::Rep::iterator> (pos->begin ().operator-> (), M.parent ().disp ()) + M.startRow (),
				    Subiterator<typename DenseMatrix<Element>::Rep::iterator> (pos->begin ().operator-> (), M.parent ().disp ()) + (M.startRow () + M.rowdim ())); }
};

template <class Element, class Submatrix>
class SubvectorFactory<DenseMatrixTag<Element>, Submatrix, typename DenseMatrix<Element>::ConstColIterator, DefaultSubvectorFactoryTrait>
{
    public:
	typedef typename DenseMatrix<Element>::ConstColIterator IteratorType;
	typedef typename DenseMatrix<Element>::ConstCol Subvector;

	Subvector MakeSubvector (const Submatrix &M, IteratorType &pos)
		{ return Subvector (Subiterator<typename DenseMatrix<Element>::Rep::const_iterator> (pos->begin ().operator-> (), M.parent ().disp ()) + M.startRow (),
				    Subiterator<typename DenseMatrix<Element>::Rep::const_iterator> (pos->begin ().operator-> (), M.parent ().disp ()) + (M.startRow () + M.rowdim ())); }
};

} // namespace LinBox

#endif // __LINBOX_matrix_dense_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
