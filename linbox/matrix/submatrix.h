/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/submatrix.h
 * Copyright 2001 B. David Saunders,
 *           2001-2002, 2010 Bradford Hovinen,
 *           2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@gmail.com>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * evolved from dense-submatrix.h by -bds, Zhendong Wan, Bradford Hovinen
 *
 * See COPYING for license information
 */

#ifndef __MATRIX_SUBMATRIX_H
#define __MATRIX_SUBMATRIX_H

#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-traits.h"

namespace LinBox
{

template <class Matrix, class Trait>
class SubmatrixBaseRowIterator;

template <class Matrix, class Trait>
class SubmatrixBaseConstRowIterator;

template <class Matrix, class Trait>
class SubmatrixBaseColIterator;

template <class Matrix, class Trait>
class SubmatrixBaseConstColIterator;

template <class Matrix, class Trait>
class SubmatrixBase;

/** Factory for subvectors
 *
 * This class allows construction of subvectors of row- and
 * column-vectors. It exists so that different subvector-classes with
 * possibly different constructor-interfaces -- as well as different
 * vector-representations -- can be used. It should be specialised for
 * the parent-matrix being used.
 */

template <class Matrix, class Trait = typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory>
class SubvectorFactory {
    public:
	class RowSubvector;
	class ConstRowSubvector;
	class ColSubvector;
	class ConstColSubvector;

	virtual RowSubvector MakeRowSubvector (SubmatrixBase<Matrix, Trait> &M, const typename Matrix::RowIterator &pos) = 0;
	virtual ConstRowSubvector MakeConstRowSubvector (const SubmatrixBase<Matrix, Trait> &M, const typename Matrix::ConstRowIterator &pos) = 0;

	virtual ColSubvector MakeColSubvector (SubmatrixBase<Matrix, Trait> &M, const typename Matrix::ColIterator &pos) = 0;
	virtual ConstColSubvector MakeConstColSubvector (const SubmatrixBase<Matrix, Trait> &M, const typename Matrix::ConstColIterator &pos) = 0;
};

/** Generic submatrix
 *
 * This is an archetype for submatrices. It is parameterised by the
 * class of containing matrix and specialisations provide support
 * for various containers.

\ingroup matrix
 */
template<class _Matrix, class Trait = typename MatrixIteratorTypes<typename MatrixTraits<_Matrix>::MatrixCategory>::MatrixCategory>
class SubmatrixBase
{
    public:
 
	class RawIterator;
	class ConstRawIterator;

	typedef _Matrix Matrix;
        typedef SubmatrixBase<Matrix, Trait> Self_t;
    
	typedef typename Matrix::Element Element;

	/** \brief 
	 *
	 * The row iterator gives the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */
	typedef SubmatrixBaseRowIterator<Matrix, Trait> RowIterator;
	typedef SubmatrixBaseConstRowIterator<Matrix, Trait> ConstRowIterator;
	typedef typename SubvectorFactory<Matrix, Trait>::RowSubvector Row;
	typedef const typename SubvectorFactory<Matrix, Trait>::ConstRowSubvector ConstRow;

	/** \brief
	 *
	 * The columns iterator gives the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */
	typedef SubmatrixBaseColIterator<Matrix, Trait> ColIterator;
	typedef SubmatrixBaseConstColIterator<Matrix, Trait> ConstColIterator;
	typedef typename SubvectorFactory<Matrix, Trait>::ColSubvector Col;
	typedef const typename SubvectorFactory<Matrix, Trait>::ConstColSubvector ConstCol;

	typedef Col Column;

	/** \brief
	 */
	SubmatrixBase () {}

	/** Constructor from an existing matrix and dimensions
	 * \param M Containing matrix in which to construct submatrix
	 * \param row Starting row
	 * \param col Starting column
	 * \param rowdim Row dimension
	 * \param coldim Column dimension
	 */
	SubmatrixBase (Matrix &M,
		       size_t row,
		       size_t col,
		       size_t rowdim,
		       size_t coldim)
		: _M (&M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
		{}

	/** Constructor from an existing submatrix and dimensions
	 * @param SM Constant reference to DenseSubmatrix from which to
	 *           construct submatrix
	 * @param row Starting row
	 * @param col Starting column
	 * @param rowdim Row dimension
	 * @param coldim Column dimension
	 */
	SubmatrixBase (const SubmatrixBase<Matrix, Trait> &SM,
		       size_t row,
		       size_t col,
		       size_t rowdim,
		       size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	/** Copy constructor
	 * @param _M SubmatrixBase to copy
	 */
	SubmatrixBase (const SubmatrixBase<Matrix, Trait> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	/** Assignment operator
	 * Assign the given submatrix to this one
	 * @param _M SubmatrixBase to assign
	 * @return Reference to this submatrix
	 */
	SubmatrixBase &operator = (const SubmatrixBase<Matrix, Trait> &SM)
	{
		_M = SM._M;
		_beg_row = SM._beg_row;
		_beg_col = SM._beg_col;
		_end_row = SM._end_row;
		_end_col = SM._end_col;
	}

	/** Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	size_t rowdim () const
		{ return _end_row - _beg_row; }

	/** Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	size_t coldim () const
		{ return _end_col - _beg_col; }

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 * @param field 
	 */
	template<class Field>
	std::istream& read (std::istream &file, const Field& field);
    
	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 * @param field
	 */
	template<class Field>
	std::ostream& write (std::ostream &os, const Field& field, bool mapleFormat = false) const;
	
	/** Set the entry at (i, j)
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	/** Get a writeable reference to an entry in the matrix
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j)
		{ return _M->refEntry (i + _beg_row, j + _beg_col); } 

	/** Get a read-only individual entry from the matrix
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const
		{ return _M->getEntry (i + _beg_row, j + _beg_col); } 

	/** Get an entry and store it in the given value
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j)
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;
 
	ColIterator colBegin ();
	ColIterator colEnd ();
	ConstColIterator colBegin () const;
	ConstColIterator colEnd () const;

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	class RawIterator;   
	class ConstRawIterator;
   
	RawIterator rawBegin ();     
	RawIterator rawEnd ();
	ConstRawIterator rawBegin () const;       
	ConstRawIterator rawEnd () const;  

	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

        class RawIndexedIterator;
        class ConstRawIndexedIterator;

        RawIndexedIterator rawIndexedBegin();
        RawIndexedIterator rawIndexedEnd();   
	ConstRawIndexedIterator rawIndexedBegin() const;
        ConstRawIndexedIterator rawIndexedEnd() const;   

	/** Access to the parent matrix */
	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	/** Start-row and -column within parent */
	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

	friend class SubmatrixBaseRowIterator<Matrix, MatrixCategories::RowMatrixTag>;
	friend class SubmatrixBaseConstRowIterator<Matrix, MatrixCategories::RowMatrixTag>;

    protected:

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

/// Specialisation for matrices indexed only by rows

template <class _Matrix>
class SubmatrixBase<_Matrix, MatrixCategories::RowMatrixTag>
{
    public:
 
	typedef _Matrix Matrix;
        typedef SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag> Self_t;
    
	typedef typename Matrix::Element Element;

	typedef SubmatrixBaseRowIterator<Matrix, MatrixCategories::RowMatrixTag> RowIterator;
	typedef SubmatrixBaseConstRowIterator<Matrix, MatrixCategories::RowMatrixTag> ConstRowIterator;
	typedef typename SubvectorFactory<Matrix, MatrixCategories::RowMatrixTag>::RowSubvector Row;
	typedef const typename SubvectorFactory<Matrix, MatrixCategories::RowMatrixTag>::ConstRowSubvector ConstRow;

	SubmatrixBase () {}

	SubmatrixBase (Matrix &M,
		       size_t row,
		       size_t col,
		       size_t rowdim,
		       size_t coldim)
		: _M (&M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
		{}

	SubmatrixBase (const SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag> &SM,
		       size_t row,
		       size_t col,
		       size_t rowdim,
		       size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	SubmatrixBase (const SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	SubmatrixBase &operator = (const SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag> &SM)
	{
		_M = SM._M;
		_beg_row = SM._beg_row;
		_beg_col = SM._beg_col;
		_end_row = SM._end_row;
		_end_col = SM._end_col;
	}

	size_t rowdim () const
		{ return _end_row - _beg_row; }

	size_t coldim () const
		{ return _end_col - _beg_col; }

	template<class Field>
	std::istream& read (std::istream &file, const Field& field);
    
	template<class Field>
	std::ostream& write (std::ostream &os, const Field& field, bool mapleFormat = false) const;
	
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	Element &refEntry (size_t i, size_t j)
		{ return _M->refEntry (i + _beg_row, j + _beg_col); } 

	const Element getEntry (size_t i, size_t j) const
		{ return _M->getEntry (i + _beg_row, j + _beg_col); } 

	Element &getEntry (Element &x, size_t i, size_t j)
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	inline RowIterator rowBegin ();
	inline RowIterator rowEnd ();
	inline ConstRowIterator rowBegin () const;
	inline ConstRowIterator rowEnd () const;

	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:
	friend class SubmatrixBaseRowIterator<Matrix, MatrixCategories::RowMatrixTag>;
	friend class SubmatrixBaseConstRowIterator<Matrix, MatrixCategories::RowMatrixTag>;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

/// Specialisation for matrices indexed only by columns

template <class _Matrix>
class SubmatrixBase<_Matrix, MatrixCategories::ColMatrixTag>
{
    public:
 
	typedef _Matrix Matrix;
        typedef SubmatrixBase<Matrix, MatrixCategories::ColMatrixTag> Self_t;
    
	typedef typename Matrix::Element Element;

	typedef SubmatrixBaseColIterator<Matrix, MatrixCategories::ColMatrixTag> ColIterator;
	typedef SubmatrixBaseConstColIterator<Matrix, MatrixCategories::ColMatrixTag> ConstColIterator;
	typedef typename SubvectorFactory<Matrix, MatrixCategories::ColMatrixTag>::ColSubvector Col;
	typedef const typename SubvectorFactory<Matrix, MatrixCategories::ColMatrixTag>::ConstColSubvector ConstCol;
	typedef Col Column;

	SubmatrixBase () {}

	SubmatrixBase (Matrix &M,
		       size_t row,
		       size_t col,
		       size_t rowdim,
		       size_t coldim)
		: _M (&M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
		{}

	SubmatrixBase (const SubmatrixBase<Matrix, MatrixCategories::ColMatrixTag> &SM,
		       size_t row,
		       size_t col,
		       size_t rowdim,
		       size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	SubmatrixBase (const SubmatrixBase<Matrix, MatrixCategories::ColMatrixTag> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	SubmatrixBase &operator = (const SubmatrixBase<Matrix, MatrixCategories::ColMatrixTag> &SM)
	{
		_M = SM._M;
		_beg_row = SM._beg_row;
		_beg_col = SM._beg_col;
		_end_row = SM._end_row;
		_end_col = SM._end_col;
	}

	size_t rowdim () const
		{ return _end_row - _beg_row; }

	size_t coldim () const
		{ return _end_col - _beg_col; }

	template<class Field>
	std::istream& read (std::istream &file, const Field& field);
    
	template<class Field>
	std::ostream& write (std::ostream &os, const Field& field, bool mapleFormat = false) const;
	
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	Element &refEntry (size_t i, size_t j)
		{ return _M->refEntry (i + _beg_row, j + _beg_col); } 

	const Element getEntry (size_t i, size_t j) const
		{ return _M->getEntry (i + _beg_row, j + _beg_col); } 

	Element &getEntry (Element &x, size_t i, size_t j)
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	inline ColIterator colBegin ();
	inline ColIterator colEnd ();
	inline ConstColIterator colBegin () const;
	inline ConstColIterator colEnd () const;

	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:
	friend class SubmatrixBaseColIterator<Matrix, MatrixCategories::ColMatrixTag>;
	friend class SubmatrixBaseConstColIterator<Matrix, MatrixCategories::ColMatrixTag>;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

template <class Matrix, class Trait>
struct MatrixTraits< SubmatrixBase<Matrix, Trait> >
{ 
	typedef SubmatrixBase<Matrix, Trait> MatrixType;
	typedef typename MatrixTraits<Matrix>::MatrixCategory MatrixCategory; 
};

template <class Matrix, class Trait>
struct MatrixTraits< const SubmatrixBase<Matrix, Trait> >
{ 
	typedef const SubmatrixBase<Matrix, Trait> MatrixType;
	typedef typename MatrixTraits<Matrix>::MatrixCategory MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/submatrix.inl"

#endif // __MATRIX_SUBMATRIX_H

