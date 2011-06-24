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
#include "linbox/matrix/traits.h"
#include "linbox/matrix/raw-iterator.h"

namespace LinBox
{

template <class Matrix, class SubvectorFactory, class Trait>
class Submatrix;

struct DefaultSubvectorFactoryTrait {};

/** Factory for subvectors
 *
 * This class allows construction of subvectors of row- and
 * column-vectors. It exists so that different subvector-classes with
 * possibly different constructor-interfaces -- as well as different
 * vector-representations -- can be used. It should be specialised for
 * the parent-matrix being used.
 */

template <class MatrixTag, class Submatrix, class _IteratorType, class SFTrait>
class SubvectorFactory {
    public:
	typedef SubvectorFactory<MatrixTag, Submatrix, _IteratorType, SFTrait> Self_t;
	typedef _IteratorType IteratorType;

	class Subvector;

	virtual Subvector MakeSubvector (Submatrix &M, const IteratorType &pos) = 0;
};

template <class Submatrix, class SubvectorFactory, class Trait>
class SubmatrixRowColIteratorPT {
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename SubvectorFactory::Subvector &reference;
	typedef typename SubvectorFactory::Subvector *pointer;
	typedef typename SubvectorFactory::Subvector value_type;
	typedef ptrdiff_t difference_type;

	SubmatrixRowColIteratorPT () {}

	SubmatrixRowColIteratorPT (Submatrix *M, typename SubvectorFactory::IteratorType pos)
		: _M (M), _pos (pos), _row_valid (false) {}

	SubmatrixRowColIteratorPT (const SubmatrixRowColIteratorPT &i) : _M (i._M), _pos (i._pos), _row (i._row), _row_valid (i._row_valid) {}

	SubmatrixRowColIteratorPT &operator = (const SubmatrixRowColIteratorPT &i)
	{
		_M = i._M;
		_pos = i._pos;
		_row = i._row;
		_row_valid = i._row_valid;
		return *this;
	}

	SubmatrixRowColIteratorPT &operator ++ () 
	{
		++_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixRowColIteratorPT operator ++ (int) 
	{
		SubmatrixRowColIteratorPT tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixRowColIteratorPT operator + (difference_type i)
		{ return SubmatrixRowColIteratorPT (_M, _pos + i); }

	SubmatrixRowColIteratorPT &operator += (difference_type i) 
	{
		_pos += i;
		_row_valid = false;
		return *this;
	}

	SubmatrixRowColIteratorPT &operator -- () 
	{
		--_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixRowColIteratorPT operator -- (int) 
	{
		SubmatrixRowColIteratorPT tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixRowColIteratorPT operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixRowColIteratorPT &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixRowColIteratorPT &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i); }

	const reference operator * ()
		{ update_row (); return _row; }

	const pointer operator -> ()
		{ update_row (); return &_row; }

	bool operator == (const SubmatrixRowColIteratorPT &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixRowColIteratorPT &c) const 
		{ return (_pos != c._pos); }

	template <class SM, class SF, class T>
	operator SubmatrixRowColIteratorPT<SM, SF, T> ()
		{ return SubmatrixRowColIteratorPT<SM, SF, T> (_M, _pos); }

private:
	template <class SM, class SF, class T>
	friend class SubmatrixRowColIteratorPT;

	Submatrix *_M;
	typename SubvectorFactory::IteratorType _pos;
	value_type _row;
	bool _row_valid;

	inline void update_row ()
	{
		if (!_row_valid) {
			_row = (SubvectorFactory ()).MakeSubvector (*_M, _pos);
			_row_valid = true;
		}
	}
};

/** Generic submatrix
 *
 * This is an archetype for submatrices. It is parameterised by the
 * class of containing matrix and specialisations provide support
 * for various containers.

\ingroup matrix
 */
template<class _Matrix, class SFTrait = DefaultSubvectorFactoryTrait, class Trait = typename MatrixIteratorTypes<typename MatrixTraits<_Matrix>::MatrixCategory>::MatrixCategory>
class Submatrix
{
    public:
 
	typedef _Matrix Matrix;
        typedef Submatrix<Matrix, SFTrait, Trait> Self_t;
	typedef typename Matrix::Tag Tag;

	typedef typename Matrix::SubmatrixType SubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstSubmatrixType;
	typedef typename Matrix::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename Matrix::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = _Matrix::rowAlign;
	static const size_t colAlign = _Matrix::colAlign;

	typedef Matrix ContainerType;

	typedef typename MatrixTraits<Matrix>::MatrixCategory MatrixCategory; 
    
	typedef typename Matrix::Element Element;

	/** \brief 
	 *
	 * The row iterator gives the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */
	typedef SubmatrixRowColIteratorPT<Self_t, SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::RowIterator, SFTrait>, Trait> RowIterator;
	typedef SubmatrixRowColIteratorPT<const Self_t, SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstRowIterator, SFTrait>, Trait> ConstRowIterator;
	typedef typename SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::RowIterator, SFTrait>::Subvector Row;
	typedef const typename SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstRowIterator, SFTrait>::Subvector ConstRow;

	/** \brief
	 *
	 * The columns iterator gives the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */
	typedef SubmatrixRowColIteratorPT<Self_t, SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::ColIterator, SFTrait>, Trait> ColIterator;
	typedef SubmatrixRowColIteratorPT<const Self_t, SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstColIterator, SFTrait>, Trait> ConstColIterator;
	typedef typename SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::ColIterator, SFTrait>::Subvector Col;
	typedef const typename SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstColIterator, SFTrait>::Subvector ConstCol;

	typedef Col Column;

	/** \brief
	 */
	Submatrix () {}

	/** Constructor from an existing matrix and dimensions
	 * \param M Containing matrix in which to construct submatrix
	 * \param row Starting row
	 * \param col Starting column
	 * \param rowdim Row dimension
	 * \param coldim Column dimension
	 */
	Submatrix (Matrix &M,
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
	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	// Parametrised version
	template <class M2>
	Submatrix (const Submatrix<M2, SFTrait, Trait> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	/** Copy constructor
	 * @param _M Submatrix to copy
	 */
	Submatrix (const Submatrix &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	/** Parametrised copy constructor
	 * @param _M Submatrix to copy
	 */
	template <class M2>
	Submatrix (const Submatrix<M2, SFTrait, Trait> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	/** Parametrised conversion
	 */
	template <class M2>
	operator Submatrix<M2, SFTrait, Trait> () const
		{ return Submatrix<M2, SFTrait, Trait> (static_cast<M2> (*_M), _beg_row, _beg_col, _end_row - _beg_row, _end_col - _beg_col); }

	/** Assignment operator
	 * Assign the given submatrix to this one
	 * @param _M Submatrix to assign
	 * @return Reference to this submatrix
	 */
	Submatrix &operator = (const Submatrix &SM)
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

	/** Set the entry at (i, j)
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	/** Erase an individual entry from the matrix.
	 * If the entry doesn't exist, then takes no action. If the underlying
	 * matrix is dense, then takes no action.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 */
	void eraseEntry (size_t i, size_t j)
		{ _M->eraseEntry (_beg_row + i, _beg_col + j); }

	/** Get an entry and store it in the given value
	 * If entry does not exist in matrix, x is left unchanged and false is returned
	 *
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return true if entry exists in matrix, false otherwise
	 */
	bool getEntry (Element &x, size_t i, size_t j) const
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	RowIterator rowBegin ()
		{ return RowIterator (this, _M->rowBegin () + _beg_row); }
	RowIterator rowEnd ()
		{ return RowIterator (this, _M->rowBegin () + _end_row); }
	ConstRowIterator rowBegin () const
		{ return ConstRowIterator (this, _M->rowBegin () + _beg_row); }
	ConstRowIterator rowEnd () const
		{ return ConstRowIterator (this, _M->rowBegin () + _end_row); }
 
	ColIterator colBegin ()
		{ return ColIterator (this, _M->colBegin () + _beg_col); }
	ColIterator colEnd ()
		{ return ColIterator (this, _M->colBegin () + _end_col); }
	ConstColIterator colBegin () const
		{ return ConstColIterator (this, _M->colBegin () + _beg_col); }
	ConstColIterator colEnd () const
		{ return ConstColIterator (this, _M->colBegin () + _end_col); }

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::VectorCategory> RawIterator;
	typedef RawIterator ConstRawIterator;

	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::VectorCategory, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	/** Access to the parent matrix */
	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	/** Start-row and -column within parent */
	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:

	template <class M2, class SFT2, class T2>
	friend class Submatrix;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

/// Specialisation for matrices indexed only by rows

template <class _Matrix, class SFTrait>
class Submatrix<_Matrix, SFTrait, MatrixCategories::RowMatrixTag>
{
    public:
 
	typedef _Matrix Matrix;
        typedef Submatrix<Matrix, SFTrait, MatrixCategories::RowMatrixTag> Self_t;
	typedef typename Matrix::Tag Tag;

	typedef typename Matrix::SubmatrixType SubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstSubmatrixType;
	typedef typename Matrix::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename Matrix::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = _Matrix::rowAlign;
	static const size_t colAlign = _Matrix::colAlign;

	typedef Matrix ContainerType;

	typedef typename MatrixTraits<Matrix>::MatrixCategory MatrixCategory; 
    
	typedef typename Matrix::Element Element;

	typedef SubmatrixRowColIteratorPT<Self_t, SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::RowIterator, SFTrait>, MatrixCategories::RowMatrixTag> RowIterator;
	typedef SubmatrixRowColIteratorPT<const Self_t, SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstRowIterator, SFTrait>, MatrixCategories::RowMatrixTag> ConstRowIterator;
	typedef typename SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::RowIterator, SFTrait>::Subvector Row;
	typedef const typename SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstRowIterator, SFTrait>::Subvector ConstRow;

	Submatrix () {}

	Submatrix (Matrix &M,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (&M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
		{}

	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	template <class M2>
	Submatrix (const Submatrix<M2, SFTrait, MatrixCategories::RowMatrixTag> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	template <class SFT2>
	Submatrix (const Submatrix<Matrix, SFT2, MatrixCategories::RowMatrixTag> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	Submatrix (const Submatrix &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	template <class M2>
	Submatrix (const Submatrix<M2, SFTrait, MatrixCategories::RowMatrixTag> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	Submatrix &operator = (const Submatrix &SM)
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

	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	void eraseEntry (size_t i, size_t j)
		{ _M->eraseEntry (_beg_row + i, _beg_col + j); }

	bool getEntry (Element &x, size_t i, size_t j) const
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	inline RowIterator rowBegin ()
		{ return RowIterator (this, _M->rowBegin () + _beg_row); }
	inline RowIterator rowEnd ()
		{ return RowIterator (this, _M->rowBegin () + _end_row); }
	inline ConstRowIterator rowBegin () const
		{ return ConstRowIterator (this, _M->rowBegin () + _beg_row); }
	inline ConstRowIterator rowEnd () const
		{ return ConstRowIterator (this, _M->rowBegin () + _end_row); }

	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::VectorCategory> RawIterator;
	typedef RawIterator ConstRawIterator;
   
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::VectorCategory, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:

	template <class M2, class SFT2, class T2>
	friend class Submatrix;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

/// Specialisation for matrices indexed only by columns

template <class _Matrix, class SFTrait>
class Submatrix<_Matrix, SFTrait, MatrixCategories::ColMatrixTag>
{
    public:
 
	typedef _Matrix Matrix;
        typedef Submatrix<Matrix, SFTrait, MatrixCategories::ColMatrixTag> Self_t;
	typedef typename Matrix::Tag Tag;

	typedef typename Matrix::SubmatrixType SubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstSubmatrixType;
	typedef typename Matrix::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename Matrix::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = _Matrix::rowAlign;
	static const size_t colAlign = _Matrix::colAlign;

	typedef Matrix ContainerType;

	typedef typename MatrixTraits<Matrix>::MatrixCategory MatrixCategory; 
    
	typedef typename Matrix::Element Element;

	typedef SubmatrixRowColIteratorPT<Self_t, SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::ColIterator, SFTrait>, MatrixCategories::ColMatrixTag> ColIterator;
	typedef SubmatrixRowColIteratorPT<const Self_t, SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstColIterator, SFTrait>, MatrixCategories::ColMatrixTag> ConstColIterator;
	typedef typename SubvectorFactory<typename Matrix::Tag, Self_t, typename Matrix::ColIterator, SFTrait>::Subvector Col;
	typedef const typename SubvectorFactory<typename Matrix::Tag, const Self_t, typename Matrix::ConstColIterator, SFTrait>::Subvector ConstCol;

	typedef Col Column;

	Submatrix () {}

	Submatrix (Matrix &M,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (&M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
		{}

	Submatrix (const Submatrix &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	template <class M2>
	Submatrix (const Submatrix<M2, SFTrait, MatrixCategories::ColMatrixTag> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	Submatrix (const Submatrix &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	template <class M2>
	Submatrix (const Submatrix<M2, SFTrait, MatrixCategories::ColMatrixTag> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	Submatrix &operator = (const Submatrix &SM)
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

	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	void eraseEntry (size_t i, size_t j)
		{ _M->eraseEntry (_beg_row + i, _beg_col + j); }

	bool getEntry (Element &x, size_t i, size_t j) const
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	ColIterator colBegin ()
		{ return ColIterator (this, _M->colBegin () + _beg_col); }
	ColIterator colEnd ()
		{ return ColIterator (this, _M->colBegin () + _end_col); }
	ConstColIterator colBegin () const
		{ return ConstColIterator (this, _M->colBegin () + _beg_col); }
	ConstColIterator colEnd () const
		{ return ConstColIterator (this, _M->colBegin () + _end_col); }

	typedef MatrixRawIterator<ConstColIterator, typename ElementVectorTraits<Element, Column>::VectorCategory> RawIterator;
	typedef RawIterator ConstRawIterator;
   
	ConstRawIterator rawBegin () const { return ConstRawIterator (colBegin (), 0, colEnd (), rowdim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (colEnd (), 0, colEnd (), rowdim ()); }

	typedef MatrixRawIndexedIterator<ConstColIterator, typename ElementVectorTraits<Element, Column>::VectorCategory, true> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (colBegin (), 0, colEnd (), rowdim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (colEnd (), coldim (), colEnd (), rowdim ()); }

	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:

	template <class M2, class SFT2, class T2>
	friend class Submatrix;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

} // namespace LinBox

#endif // __MATRIX_SUBMATRIX_H

