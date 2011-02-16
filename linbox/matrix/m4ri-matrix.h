/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#include "linbox/matrix/matrix-domain-gf2.h"
#include "linbox/vector/bit-subvector-word-aligned.h"
#include "linbox/matrix/submatrix.h"

namespace LinBox
{

class MatrixDomainM4RI;

/** Specialisation of WordTraits for M4RI-words */
template <>
struct WordTraits<word> {
	typedef word Word;
	static const unsigned int bits = 64;
	static const unsigned int logof_size = 6;
	static const unsigned int pos_mask = 0x3F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 32);
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};

/** Wrapper for dense zero-one matrices in M4RI
 */
class M4RIMatrix
{
    public:

	typedef bool Element;
	typedef mzd_t *Rep;
        typedef M4RIMatrix Self_t;

	///
	M4RIMatrix ()
	{}

	/** Constructor.
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	M4RIMatrix (size_t m, size_t n)
		: _rep (mzd_init (m, n))
	{}

	/** Construct a dense submatrix of a given dense matrix
	 */
	M4RIMatrix (M4RIMatrix &M, size_t beg_row, size_t beg_col, size_t m, size_t n)
		: _rep (mzd_init_window (M._rep, beg_row, beg_col, beg_row + m, beg_col + n))
	{}

	///
	M4RIMatrix (const M4RIMatrix &M)
		: _rep (M._rep)
	{}

	~M4RIMatrix () { mzd_free (_rep); }

	///
	M4RIMatrix &operator = (const M4RIMatrix &M)
	{
		_rep = M._rep;
		return *this;
	}

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
	{
		_rep = mzd_init (m, n);
	}

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 * @param F Field over which to read
	 */
	template <class Field>
	std::istream &read (std::istream &file, const Field &F);

	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 * @param F Field over which to write
	 */
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const;

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, bool a_ij)
		{ if (a_ij) mzd_write_bit (_rep, i, j, TRUE); else mzd_write_bit (_rep, i, j, FALSE); }

	/** Get a writeable reference to the entry in the (i, j) position.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @returns Reference to matrix entry
	 */
	BitVectorReference<word *, BigEndian<word> > refEntry (size_t i, size_t j)
		{ return (*this)[i][j]; }

	/** Get a read-only reference to the entry in the (i, j) position.
	 * @param i Row index
	 * @param j Column index
	 * @returns Const reference to matrix entry
	 */
	bool getEntry (size_t i, size_t j) const
		{ return (mzd_read_bit (_rep, i, j) == TRUE) ? true : false; }

	/** Copy the (i, j) entry into x, and return a reference to x.
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @returns Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const
		{ x = (mzd_read_bit (_rep, i, j) == TRUE) ? true : false; return x; }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	typedef BitSubvectorWordAligned<word *, const word *, BigEndian<word> > Row;  
	typedef BitSubvectorWordAligned<const word *, const word *, BigEndian<word> > ConstRow;

	class RowIterator;
	class ConstRowIterator;

	inline RowIterator rowBegin ();
	inline RowIterator rowEnd ();
	inline ConstRowIterator rowBegin () const;
	inline ConstRowIterator rowEnd () const;

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

    protected:

	friend class MatrixDomainM4RI;

	Rep    _rep;
};

template <>
struct MatrixTraits<M4RIMatrix>
{ 
	typedef M4RIMatrix MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

template <>
struct MatrixTraits<const M4RIMatrix>
{ 
	typedef const M4RIMatrix MatrixType;
	typedef MatrixCategories::ZeroOneRowMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/m4ri-matrix.inl"

#endif // __LINBOX_MATRIX_M4RI_MATRIX_H
