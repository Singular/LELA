/* linbox/matrix/m4ri-matrix.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Wrapper for libm4ri-matrices
 *
 * Evolved from dense.h
 */

#ifndef __LINBOX_MATRIX_M4RI_MATRIX_H
#define __LINBOX_MATRIX_M4RI_MATRIX_H

#include <vector>

#ifndef __LINBOX_HAVE_M4RI
#  error "This header file requires that LinBox be configured with libm4ri enabled. Please ensure that libm4ri is properly installed and re-run configure."
#endif

#include <m4ri/m4ri.h>

#include "linbox/vector/bit-subvector-word-aligned.h"
#include "linbox/vector/bit-iterator.h"
#include "linbox/vector/bit-subvector.h"
#include "linbox/matrix/submatrix.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/raw-iterator.h"

namespace LinBox
{

class M4RIMatrix;

/** Base class for wrapper for dense zero-one matrices in M4RI
 *
 * This is a base class providing just the M4RI-wrapper. Use @ref
 * M4RIMatrix instead.
 */
class M4RIMatrixBase
{
    public:

	typedef bool Element;
	typedef mzd_t *Rep;
        typedef M4RIMatrixBase Self_t;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory;
	struct Tag {};

	virtual ~M4RIMatrixBase () { if (_rep != NULL) mzd_free (_rep); }

	size_t rowdim () const
		{ return _rep->nrows; }

	size_t coldim () const
		{ return _rep->ncols; }

	void resize (size_t m, size_t n)
		// { if (_rep != NULL) mzd_free (_rep); _rep = mzd_init (m, n); }
		{ _rep = mzd_init (m, n); }

	void setEntry (size_t i, size_t j, bool a_ij)
		{ if (a_ij) mzd_write_bit (_rep, i, j, TRUE); else mzd_write_bit (_rep, i, j, FALSE); }

	void eraseEntry (size_t i, size_t j) {}

	bool getEntry (Element &x, size_t i, size_t j) const
		{ x = (mzd_read_bit (_rep, i, j) == TRUE) ? true : false; return true; }

	template <class Vector>
	Vector &columnDensity (Vector &v) const
		{ std::fill (v.begin (), v.end (), _rep->nrows); return v; }

	Rep    _rep;

    protected:

	M4RIMatrixBase (Rep rep) : _rep (rep) {}

    private:

	M4RIMatrixBase (const M4RIMatrixBase &M);
};

// Forward declaration
class M4RISubmatrix;

/** Wrapper for dense zero-one-matrices in M4RI
 *
 * This class wraps a dense matrix from M4RI into LinBox-matrix. It
 * provides row-iterators but not column-iterators.
 */
class M4RIMatrix : public M4RIMatrixBase
{
	template <class Iterator, class ConstIterator, class MatrixPointer>
	class RowIteratorPT
	{
	public:
		typedef BitSubvectorWordAligned<Iterator, ConstIterator, BigEndian<word> > Row;  
		typedef BitSubvectorWordAligned<ConstIterator, ConstIterator, BigEndian<word> > ConstRow;

		typedef Row value_type;

		typedef typename std::iterator_traits<typename Row::word_iterator>::difference_type difference_type;

		RowIteratorPT (MatrixPointer M, size_t idx)
			: _M (M), _idx (idx) { make_row (); }
    
		RowIteratorPT () {}
    
		RowIteratorPT (const RowIteratorPT &colp)
			: _M (colp._M), _idx (colp._idx), _row (colp._row) {}

		template <class It, class CIt, class MP>
		RowIteratorPT (const RowIteratorPT<It, CIt, MP> &colp)
			: _M (colp._M), _idx (colp._idx), _row (colp._row) {}
    
		template <class It, class CIt, class MP>
		RowIteratorPT &operator = (const RowIteratorPT<It, CIt, MP> &colp)
		{
			_M = colp._M;
			_idx = colp._idx;
			_row = colp._row;
			return *this;
		}
    
		RowIteratorPT &operator ++ ()
		{
			++_idx;
			make_row ();
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
			--_idx;
			make_row ();
			return *this;
		}

		RowIteratorPT operator -- (int)
		{
			RowIteratorPT tmp (*this);
			--*this;
			return tmp;
		}

		RowIteratorPT operator+ (int i) const
			{ return RowIteratorPT (_M, _idx + i); }

		difference_type operator- (const RowIteratorPT &c) const
			{ return c._idx - _idx; }

		RowIteratorPT &operator += (int i)
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

		bool operator == (const RowIteratorPT& c) const
			{ return _idx == c._idx; }

		template <class It, class CIt, class MP>
		bool operator == (const RowIteratorPT<It, CIt, MP> &c) const
			{ return _idx == c._idx; }

		bool operator != (const RowIteratorPT& c) const
			{ return _idx != c._idx; }

		template <class It, class CIt, class MP>
		bool operator != (const RowIteratorPT<It, CIt, MP> &c) const
			{ return _idx != c._idx; }

	private:
		inline void make_row ()
			{ if (_M->_rep->rows != NULL) _row = Row (_M->_rep->rows[_idx], _M->_rep->rows[_idx] + _M->_rep->width, _M->_rep->ncols); }

		MatrixPointer _M;
		size_t _idx;
		Row _row;

		template <class It, class CIt, class MP>
		friend class RowIteratorPT;
	};

public:
	typedef M4RISubmatrix SubmatrixType;
	typedef const M4RISubmatrix ConstSubmatrixType;
	typedef M4RIMatrix WordAlignedSubmatrixType;
	typedef const M4RIMatrix ConstWordAlignedSubmatrixType;
	
	typedef BitSubvectorWordAligned<word *, const word *, BigEndian<word> > Row;  
	typedef BitSubvectorWordAligned<const word *, const word *, BigEndian<word> > ConstRow;

	typedef RowIteratorPT<word *, const word *, M4RIMatrixBase *> RowIterator;
	typedef RowIteratorPT<const word *, const word *, const M4RIMatrixBase *> ConstRowIterator;

	M4RIMatrix ()
		: M4RIMatrixBase (NULL)
	{}

	M4RIMatrix (size_t m, size_t n)
		: M4RIMatrixBase (mzd_init (m, n)), _start_row (0), _start_col (0)
	{}

	M4RIMatrix (M4RIMatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: M4RIMatrixBase (mzd_init_window (M._rep, beg_row, beg_col, beg_row + m, beg_col + n)),
		  _start_row (M.startRow () + beg_row),
		  _start_col (M.startCol () + beg_col)
		{ linbox_check (_rep->offset == 0); }

	M4RIMatrix (const M4RIMatrix &M)
		: M4RIMatrixBase (mzd_init_window (M._rep, 0, 0, M.rowdim (), M.coldim ())),
		  _start_row (M.startRow ()), _start_col (M.startCol ())
	{}

	M4RIMatrix (VectorStream<Row> &vs)
		: M4RIMatrixBase (mzd_init (vs.size (), vs.dim ())), _start_row (0), _start_col (0)
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	M4RIMatrix &operator = (const M4RIMatrix &M)
	{
		_rep = M._rep;
		_start_row = M._start_row;
		_start_col = M._start_col;
		return *this;
	}

	inline RowIterator rowBegin ()
		{ return RowIterator (this, 0); }
	inline RowIterator rowEnd ()
		{ return RowIterator (this, rowdim ()); }
	inline ConstRowIterator rowBegin () const
		{ return ConstRowIterator (this, 0); }
	inline ConstRowIterator rowEnd () const
		{ return ConstRowIterator (this, rowdim ()); }

	Row operator[] (size_t i)
		{ return Row (_rep->rows[i], _rep->rows[i] + _rep->width, _rep->ncols); }

	ConstRow operator[] (size_t i) const
		{ return ConstRow (_rep->rows[i], _rep->rows[i] + _rep->width, _rep->ncols); }

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::DenseZeroOneVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const
		{ return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const
		{ return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::DenseZeroOneVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const
		{ return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const
		{ return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	size_t startRow () const { return _start_row; }
	size_t startCol () const { return _start_col; }

private:
	size_t _start_row, _start_col;
};

class M4RISubmatrix : public M4RIMatrixBase
{
	template <class Iterator, class MatrixPointer>
	class RowIteratorPT
	{
	public:

		typedef BitSubvector<Iterator> value_type;

		typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
		typedef typename Iterator::word_type word_type;

		RowIteratorPT (MatrixPointer M, size_t idx)
			: _M (M), _idx (idx), _last_offset ((M->_rep->offset + M->_rep->ncols) & WordTraits<word_type>::pos_mask)
			{ make_row (); }
    
		RowIteratorPT () {}
    
		RowIteratorPT (const RowIteratorPT &rowp)
			: _M (rowp._M), _idx (rowp._idx), _last_offset (rowp._last_offset), _row (rowp._row) {}

		template <class It, class MP>
		RowIteratorPT (const RowIteratorPT<It, MP> &rowp)
			: _M (rowp._M), _idx (rowp._idx), _last_offset (rowp._last_offset), _row (rowp._row) {}
    
		template <class It, class MP>
		RowIteratorPT &operator = (const RowIteratorPT<It, MP> &rowp)
		{
			_M = rowp._M;
			_idx = rowp._idx;
			_last_offset = rowp._last_offset;
			_row = rowp._row;
			return *this;
		}
    
		RowIteratorPT &operator ++ ()
		{
			++_idx;
			make_row ();
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
			--_idx;
			make_row ();
			return *this;
		}

		RowIteratorPT operator -- (int)
		{
			RowIteratorPT tmp (*this);
			--*this;
			return tmp;
		}

		RowIteratorPT operator+ (int i) const
			{ return RowIteratorPT (_M, _idx + i); }

		difference_type operator- (const RowIteratorPT &c) const
			{ return c._idx - _idx; }

		RowIteratorPT &operator += (int i)
		{
			_idx += i;
			make_row ();
			return *this;
		}

		value_type operator [] (int i) const
			{ return value_type (Iterator (_M->_rep->rows[_idx + i], _M->_rep->offset),
					     Iterator (_M->_rep->rows[_idx + i] + _M->_rep->width - ((_last_offset == 0) ? 0 : 1), _last_offset)); }

		value_type *operator -> ()
			{ return &_row; }

		value_type &operator * ()
			{ return _row; }

		bool operator == (const RowIteratorPT& c) const
			{ return _idx == c._idx; }

		template <class It, class MP>
		bool operator == (const RowIteratorPT<It, MP> &c) const
			{ return _idx == c._idx; }

		bool operator != (const RowIteratorPT& c) const
			{ return _idx != c._idx; }

		template <class It, class MP>
		bool operator != (const RowIteratorPT<It, MP> &c) const
			{ return _idx != c._idx; }

	private:
		inline void make_row ()
		{
			if (_M->_rep->rows != NULL)
				_row = value_type (Iterator (_M->_rep->rows[_idx], _M->_rep->offset),
						   Iterator (_M->_rep->rows[_idx] + _M->_rep->width - ((_last_offset == 0) ? 0 : 1), _last_offset));
		}

		MatrixPointer _M;
		size_t _idx, _last_offset;
		value_type _row;

		template <class It, class MP>
		friend class RowIteratorPT;
	};

public:
	typedef M4RISubmatrix SubmatrixType;
	typedef const M4RISubmatrix ConstSubmatrixType;
	
	typedef BitSubvector<BitVectorIterator<word *, const word *, BigEndian<word> > > Row;
	typedef BitSubvector<BitVectorIterator<const word *, const word *, BigEndian<word> > > ConstRow;

	typedef RowIteratorPT<BitVectorIterator<word *, const word *, BigEndian<word> >, M4RIMatrixBase *> RowIterator;
	typedef RowIteratorPT<BitVectorIterator<const word *, const word *, BigEndian<word> >, const M4RIMatrixBase *> ConstRowIterator;

	M4RISubmatrix ()
		: M4RIMatrixBase (NULL)
	{}

	M4RISubmatrix (const M4RIMatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: M4RIMatrixBase (mzd_init_window (M._rep, beg_row, beg_col, beg_row + m, beg_col + n)),
		  _start_row (M.startRow () + beg_row),
		  _start_col (M.startCol () + beg_col)
	{}

	M4RISubmatrix (const M4RISubmatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: M4RIMatrixBase (mzd_init_window (M._rep, beg_row, beg_col, beg_row + m, beg_col + n)),
		  _start_row (M.startRow () + beg_row),
		  _start_col (M.startCol () + beg_col)
	{}

	M4RISubmatrix (const M4RISubmatrix &M)
		: M4RIMatrixBase (mzd_init_window (M._rep, 0, 0, M.rowdim (), M.coldim ())),
		  _start_row (M.startRow ()),
		  _start_col (M.startCol ())
	{}

	M4RISubmatrix &operator = (const M4RISubmatrix &M)
	{
		_rep = M._rep;
		_start_row = M._start_row;
		_start_col = M._start_col;
		return *this;
	}

	inline RowIterator rowBegin ()
		{ return RowIterator (this, 0); }
	inline RowIterator rowEnd ()
		{ return RowIterator (this, rowdim ()); }
	inline ConstRowIterator rowBegin () const
		{ return ConstRowIterator (this, 0); }
	inline ConstRowIterator rowEnd () const
		{ return ConstRowIterator (this, rowdim ()); }

	Row operator[] (size_t i)
	{
		size_t last_offset = (_rep->offset + _rep->ncols) & WordTraits<word>::pos_mask;

		return Row (BitVectorIterator<word *, const word *, BigEndian<word> > (_rep->rows[i], _rep->offset),
			    BitVectorIterator<word *, const word *, BigEndian<word> > (_rep->rows[i] + _rep->width - ((last_offset == 0) ? 0 : 1), last_offset));
	}

	ConstRow operator[] (size_t i) const
	{
		size_t last_offset = (_rep->offset + _rep->ncols) & WordTraits<word>::pos_mask;

		return ConstRow (BitVectorIterator<const word *, const word *, BigEndian<word> > (_rep->rows[i], _rep->offset),
				 BitVectorIterator<const word *, const word *, BigEndian<word> > (_rep->rows[i] + _rep->width - ((last_offset == 0) ? 0 : 1), last_offset));
	}

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::DenseZeroOneVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const
		{ return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const
		{ return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::DenseZeroOneVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const
		{ return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const
		{ return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	size_t startRow () const { return _start_row; }
	size_t startCol () const { return _start_col; }

private:
	size_t _start_row, _start_col;
};

template <>
class DenseMatrix<bool> : public M4RIMatrix
{
public:
	DenseMatrix ()
	{}

	DenseMatrix (size_t m, size_t n)
		: M4RIMatrix (m, n)
	{}

	DenseMatrix (DenseMatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: M4RIMatrix (M, beg_row, beg_col, m, n)
	{}

	DenseMatrix (const DenseMatrix &M)
		: M4RIMatrix (M)
	{}

	DenseMatrix (VectorStream<Row> &vs)
		: M4RIMatrix (vs)
	{}
};

} // namespace LinBox

#endif // __LINBOX_MATRIX_M4RI_MATRIX_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
