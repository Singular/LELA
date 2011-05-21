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

/** Wrapper for dense zero-one matrices in M4RI
 */
class M4RIMatrixBase
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

	typedef bool Element;
	typedef mzd_t *Rep;
        typedef M4RIMatrixBase Self_t;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory;
	struct Tag {};

	typedef BitSubvectorWordAligned<word *, const word *, BigEndian<word> > Row;  
	typedef BitSubvectorWordAligned<const word *, const word *, BigEndian<word> > ConstRow;

	typedef RowIteratorPT<word *, const word *, M4RIMatrixBase *> RowIterator;
	typedef RowIteratorPT<const word *, const word *, const M4RIMatrixBase *> ConstRowIterator;

	virtual ~M4RIMatrixBase () { if (_rep != NULL) mzd_free (_rep); }

	/** Get the number of rows in the matrix
	 * @returns Number of rows in matrix
	 */
	size_t rowdim () const
		{ return _rep->nrows; }

	/** Get the number of columns in the matrix
	 * @returns Number of columns in matrix
	 */
	size_t coldim () const
		{ return _rep->ncols; }

	/** Resize the matrix to the given dimensions
	 * The state of the matrix's entries after a call to this method is
	 * undefined
	 * @param m Number of rows
	 * @param n Number of columns
	 */
	void resize (size_t m, size_t n)
		// { if (_rep != NULL) mzd_free (_rep); _rep = mzd_init (m, n); }
		{ _rep = mzd_init (m, n); }

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, bool a_ij)
		{ if (a_ij) mzd_write_bit (_rep, i, j, TRUE); else mzd_write_bit (_rep, i, j, FALSE); }

	/** Does nothing. Provided for compatibility
	 * @param i
	 * @param j
	 */
	void eraseEntry (size_t i, size_t j) {}

	/** Copy the (i, j) entry into x, and return a reference to x.
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @returns true
	 */
	bool getEntry (Element &x, size_t i, size_t j) const
		{ x = (mzd_read_bit (_rep, i, j) == TRUE) ? true : false; return true; }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	inline RowIterator rowBegin ()
		{ return RowIterator (this, 0); }
	inline RowIterator rowEnd ()
		{ return RowIterator (this, rowdim ()); }
	inline ConstRowIterator rowBegin () const
		{ return ConstRowIterator (this, 0); }
	inline ConstRowIterator rowEnd () const
		{ return ConstRowIterator (this, rowdim ()); }

	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 * @param i Row index
	 */
	Row operator[] (size_t i)
		{ return Row (_rep->rows[i], _rep->rows[i] + _rep->width, _rep->ncols); }

	ConstRow operator[] (size_t i) const
		{ return ConstRow (_rep->rows[i], _rep->rows[i] + _rep->width, _rep->ncols); }

	/** Compute column density
	 */

	template <class Vector>
	Vector &columnDensity (Vector &v) const
		{ std::fill (v.begin (), v.end (), _rep->nrows); return v; }

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

	Rep    _rep;

    protected:

	friend class RowIteratorPT<word *, const word *, M4RIMatrixBase *>;
	friend class RowIteratorPT<const word *, const word *, const M4RIMatrixBase *>;

	M4RIMatrixBase (Rep rep) : _rep (rep) {}

    private:

	M4RIMatrixBase (const M4RIMatrixBase &M);
};

class M4RIMatrix : public M4RIMatrixBase
{
public:
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory;
	
	///
	M4RIMatrix ()
		: M4RIMatrixBase (NULL)
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	M4RIMatrix (size_t m, size_t n)
		: M4RIMatrixBase (mzd_init (m, n))
	{}

	/** Construct a dense submatrix of a given dense matrix
	 */
	M4RIMatrix (M4RIMatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: M4RIMatrixBase (mzd_init_window (M._rep, beg_row, beg_col, beg_row + m, beg_col + n))
	{}

	///
	M4RIMatrix (const M4RIMatrix &M)
		: M4RIMatrixBase (mzd_init_window (M._rep, 0, 0, M.rowdim (), M.coldim ()))
	{}

	M4RIMatrix (VectorStream<Row> &vs)
		: M4RIMatrixBase (mzd_init (vs.size (), vs.dim ()))
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	M4RIMatrix &operator = (const M4RIMatrix &M)
	{
		_rep = M._rep;
		return *this;
	}

private:

	M4RIMatrix (Rep rep) : M4RIMatrixBase (rep) {}

	friend class Submatrix<M4RIMatrix>;
	friend class Submatrix<const M4RIMatrix>;
};

template <class Submatrix>
class SubvectorFactory<M4RIMatrixBase::Tag, Submatrix, M4RIMatrixBase::RowIterator, DefaultSubvectorFactoryTrait>
{
    public:
	typedef M4RIMatrixBase::RowIterator IteratorType;
	typedef BitSubvector<BitVectorIterator<word *, const word *, BigEndian<word> >, BitVectorConstIterator<const word *, BigEndian<word> > > Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

template <class Submatrix>
class SubvectorFactory<M4RIMatrixBase::Tag, Submatrix, M4RIMatrixBase::ConstRowIterator, DefaultSubvectorFactoryTrait>
{
    public:
	typedef M4RIMatrixBase::ConstRowIterator IteratorType;
	typedef BitSubvector<BitVectorConstIterator<const word *, BigEndian<word> >, BitVectorConstIterator<const word *, BigEndian<word> > > Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (pos->begin () + M.startCol (), pos->begin () + (M.startCol () + M.coldim ())); }
};

/* Specialisation of Submatrix to M4RIMatrix using facilities in M4RI */

template<>
class Submatrix<M4RIMatrix> : public Submatrix<M4RIMatrixBase>
{
    public:
	typedef M4RIMatrix Matrix;
	typedef MatrixCategories::ZeroOneRowMatrixTag Trait;
 
	Submatrix () : _rep (NULL) {}
	Submatrix (Matrix &M,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<M4RIMatrixBase> (M, row, col, rowdim, coldim),
		  _rep (mzd_init_window (M._rep, row, col, row + rowdim, col + coldim))
		{}

	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<M4RIMatrixBase> ((Submatrix<M4RIMatrixBase> &) SM, row, col, rowdim, coldim),
		  _rep (mzd_init_window (SM._rep._rep, row, col, row + rowdim, col + coldim))
		{}

    private:
	friend class Submatrix<const M4RIMatrix>;

	M4RIMatrix _rep;
};

template<>
class Submatrix<const M4RIMatrix> : public Submatrix<const M4RIMatrixBase>
{
    public:
	typedef const M4RIMatrix Matrix;
	typedef MatrixCategories::ZeroOneRowMatrixTag Trait;
 
	Submatrix () : _rep (NULL) {}
	Submatrix (const Matrix &M,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<const M4RIMatrixBase> (M, row, col, rowdim, coldim),
		  _rep (mzd_init_window (M._rep, row, col, row + rowdim, col + coldim))
		{}

	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<const M4RIMatrixBase> (SM, row, col, rowdim, coldim),
		  _rep (mzd_init_window (SM._rep._rep, row, col, row + rowdim, col + coldim))
		{}

	Submatrix (const Submatrix<M4RIMatrix> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<const M4RIMatrixBase> (SM, row, col, rowdim, coldim),
		  _rep (mzd_init_window (SM._rep._rep, row, col, row + rowdim, col + coldim))
		{}

    private:
	const M4RIMatrix _rep;
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

template <>
class Submatrix<DenseMatrix<bool> > : public Submatrix<M4RIMatrix>
{
    public:
	typedef DenseMatrix<bool> Matrix;
	typedef MatrixCategories::ZeroOneRowMatrixTag Trait;
 
	Submatrix () {}
	Submatrix (Matrix &M,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<M4RIMatrix> (M, row, col, rowdim, coldim)
		{}

	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<M4RIMatrix> (SM, row, col, rowdim, coldim)
		{}
};

template <>
class Submatrix<const DenseMatrix<bool> > : public Submatrix<const M4RIMatrix>
{
    public:
	typedef DenseMatrix<bool> Matrix;
	typedef MatrixCategories::ZeroOneRowMatrixTag Trait;
 
	Submatrix () {}
	Submatrix (const Matrix &M,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<const M4RIMatrix> (M, row, col, rowdim, coldim)
		{}

	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<const M4RIMatrix> (SM, row, col, rowdim, coldim)
		{}

	Submatrix (const Submatrix<DenseMatrix<bool> > &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: Submatrix<const M4RIMatrix> (SM, row, col, rowdim, coldim)
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
