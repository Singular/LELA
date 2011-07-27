/* lela/matrix/submatrix.h
 * Copyright 2001 B. David Saunders,
 *           2001-2002, 2010 Bradford Hovinen,
 *           2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@gmail.com>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * evolved from dense-submatrix.h by -bds, Zhendong Wan, Bradford Hovinen
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_SUBMATRIX_H
#define __LELA_MATRIX_SUBMATRIX_H

#include <cstddef>

#include "lela/util/debug.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/raw-iterator.h"

namespace LELA
{

template <class Matrix, class AlignedTrait, class Trait>
class Submatrix;

/// Tag which indicates that the submatrix is not aligned
///
/// \ingroup matrix
struct ArbitrarySubmatrixTag {};

/// Tag which indicates that the submatrix is aligned
///
/// \ingroup matrix
struct AlignedSubmatrixTag {};

/// Tag which indicates that the submatrix is const
///
/// \ingroup matrix
struct RCConstIteratorTag {};

/// Tag which indicates that the submatrix is mutable
///
/// \ingroup matrix
struct RCMutableIteratorTag {};

/// Structure to select the type of submatrix to be used depending on
/// whether the parent submatrix is aligned or not
///
/// \ingroup matrix
template <class Matrix, class AlignedTrait>
struct SubSubmatrixTraits
{
	typedef typename Matrix::SubmatrixType AlignedSubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstAlignedSubmatrixType;
};

template <class Matrix>
struct SubSubmatrixTraits<Matrix, AlignedSubmatrixTag>
{
	typedef typename Matrix::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename Matrix::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;
};

/// Structure to select which subvector of a row to take
///
/// \ingroup matrix
template <class Row, class Element, class AlignedTrait, class ConstTrait>
struct SubmatrixRowTraits
	{ typedef typename ElementVectorTraits<Element, Row>::SubvectorType SubvectorType; };

template <class Row, class Element>
struct SubmatrixRowTraits<Row, Element, AlignedSubmatrixTag, RCMutableIteratorTag>
	{ typedef typename ElementVectorTraits<Element, Row>::AlignedSubvectorType SubvectorType; };

template <class Row, class Element>
struct SubmatrixRowTraits<Row, Element, ArbitrarySubmatrixTag, RCConstIteratorTag>
	{ typedef typename ElementVectorTraits<Element, Row>::ConstSubvectorType SubvectorType; };

template <class Row, class Element>
struct SubmatrixRowTraits<Row, Element, AlignedSubmatrixTag, RCConstIteratorTag>
	{ typedef typename ElementVectorTraits<Element, Row>::ConstAlignedSubvectorType SubvectorType; };

template <class IteratorType, class Element, class AlignedTrait, class ConstTrait>
class SubmatrixRowColIteratorPT {
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename SubmatrixRowTraits<typename IteratorType::value_type, Element, AlignedTrait, ConstTrait>::SubvectorType &reference;
	typedef typename SubmatrixRowTraits<typename IteratorType::value_type, Element, AlignedTrait, ConstTrait>::SubvectorType *pointer;
	typedef typename SubmatrixRowTraits<typename IteratorType::value_type, Element, AlignedTrait, ConstTrait>::SubvectorType  value_type;
	typedef ptrdiff_t difference_type;

	SubmatrixRowColIteratorPT () {}

	SubmatrixRowColIteratorPT (IteratorType pos, size_t begin_idx, size_t end_idx)
		: _pos (pos), _begin_idx (begin_idx), _end_idx (end_idx), _row_valid (false) {}

	SubmatrixRowColIteratorPT (const SubmatrixRowColIteratorPT &i) : _pos (i._pos), _begin_idx (i._begin_idx), _end_idx (i._end_idx), _row (i._row), _row_valid (i._row_valid) {}

	SubmatrixRowColIteratorPT &operator = (const SubmatrixRowColIteratorPT &i)
	{
		_pos = i._pos;
		_begin_idx = i._begin_idx;
		_end_idx = i._end_idx;
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
		{ return SubmatrixRowColIteratorPT (_pos + i, _begin_idx, _end_idx); }

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

	difference_type operator - (SubmatrixRowColIteratorPT i) const 
		{ return _pos - i._pos; }

	template <class I, class AT, class CT>
	difference_type operator - (SubmatrixRowColIteratorPT<I, Element, AT, CT> i) const 
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

	template <class I, class AT, class CT>
	operator SubmatrixRowColIteratorPT<I, Element, AT, CT> ()
		{ return SubmatrixRowColIteratorPT<I, Element, AT, CT> (_pos, _begin_idx, _end_idx); }

private:
	template <class I, class E, class AT, class CT>
	friend class SubmatrixRowColIteratorPT;

	IteratorType _pos;
	size_t _begin_idx;
	size_t _end_idx;
	value_type _row;
	bool _row_valid;

	inline void update_row ()
	{
		if (!_row_valid) {
			_row = value_type (*_pos, _begin_idx, _end_idx);
			_row_valid = true;
		}
	}
};

/** Generic submatrix
 *
 * This class permits the construction of a submatrix of an arbitrary
 * matrix. It should however not be used directly by the user; please
 * use one of the types
 * typename Matrix::SubmatrixType,
 * typename Matrix::ConstSubmatrixType,
 * typename Matrix::AlignedSubmatrixType,
 * typename Matrix::ConstAlignedSubmatrixType,
 * as appropriate.
 *
 * If defining one's own matrix-type, this class may be used in the
 * definition of the above submatrix-types. See e.g. SparseMatrix for
 * an example.
 *
 * @param _Matrix Parent matrix-type
 *
 * @param AlignedTrait Indicates whether the submatrix is aligned or
 * not. Should be ArbitrarySubmatrixTag for arbitrary submatrices or
 * AlignedSubmatrixTag for aligned submatrices.
 *
 * \ingroup matrix
 */
template<class _Matrix, class AlignedTrait = ArbitrarySubmatrixTag, class Trait = typename _Matrix::IteratorType>
class Submatrix
{
    public:
 
	typedef _Matrix Matrix;
        typedef Submatrix<Matrix, AlignedTrait, Trait> Self_t;

	typedef typename Matrix::SubmatrixType SubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstSubmatrixType;
	typedef typename SubSubmatrixTraits<Matrix, AlignedTrait>::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename SubSubmatrixTraits<Matrix, AlignedTrait>::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = _Matrix::rowAlign;
	static const size_t colAlign = _Matrix::colAlign;

	typedef Matrix ContainerType;

	typedef typename Matrix::IteratorType IteratorType; 
	typedef typename Matrix::StorageType StorageType; 
    
	typedef typename Matrix::Element Element;

	/** \brief 
	 *
	 * The row iterator gives the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */
	typedef SubmatrixRowColIteratorPT<typename Matrix::RowIterator, Element, AlignedTrait, RCMutableIteratorTag> RowIterator;
	typedef SubmatrixRowColIteratorPT<typename Matrix::ConstRowIterator, Element, AlignedTrait, RCConstIteratorTag> ConstRowIterator;
	typedef typename SubmatrixRowTraits<typename Matrix::Row, Element, AlignedTrait, RCMutableIteratorTag>::SubvectorType Row;
	typedef typename SubmatrixRowTraits<typename Matrix::ConstRow, Element, AlignedTrait, RCConstIteratorTag>::SubvectorType ConstRow;

	/** \brief
	 *
	 * The columns iterator gives the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */
	typedef SubmatrixRowColIteratorPT<typename Matrix::ColIterator, Element, AlignedTrait, RCMutableIteratorTag> ColIterator;
	typedef SubmatrixRowColIteratorPT<typename Matrix::ConstColIterator, Element, AlignedTrait, RCConstIteratorTag> ConstColIterator;
	typedef typename SubmatrixRowTraits<typename Matrix::Col, Element, AlignedTrait, RCMutableIteratorTag>::SubvectorType Col;
	typedef typename SubmatrixRowTraits<typename Matrix::ConstCol, Element, AlignedTrait, RCConstIteratorTag>::SubvectorType ConstCol;

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
	template <class M2, class AT2>
	Submatrix (const Submatrix<M2, AT2, Trait> &SM,
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
	Submatrix (const Submatrix<M2, AlignedTrait, Trait> &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	/** Parametrised conversion
	 */
	template <class M2>
	operator Submatrix<M2, AlignedTrait, Trait> () const
		{ return Submatrix<M2, AlignedTrait, Trait> (static_cast<M2> (*_M), _beg_row, _beg_col, _end_row - _beg_row, _end_col - _beg_col); }

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

	RowIterator rowBegin ()            { return RowIterator (_M->rowBegin () + _beg_row, _beg_col, _end_col); }
	RowIterator rowEnd   ()            { return RowIterator (_M->rowBegin () + _end_row, _beg_col, _end_col); }
	ConstRowIterator rowBegin () const { return ConstRowIterator (_M->rowBegin () + _beg_row, _beg_col, _end_col); }
	ConstRowIterator rowEnd   () const { return ConstRowIterator (_M->rowBegin () + _end_row, _beg_col, _end_col); }
 
	ColIterator colBegin ()            { return ColIterator (_M->colBegin () + _beg_col, _beg_row, _end_row); }
	ColIterator colEnd   ()            { return ColIterator (_M->colBegin () + _end_col, _beg_row, _end_row); }
	ConstColIterator colBegin () const { return ConstColIterator (_M->colBegin () + _beg_col, _beg_row, _end_row); }
	ConstColIterator colEnd   () const { return ConstColIterator (_M->colBegin () + _end_col, _beg_row, _end_row); }

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
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

	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
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

	template <class M2, class AT2, class T2>
	friend class Submatrix;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

/// Specialisation for matrices indexed only by rows

template <class _Matrix, class AlignedTrait>
class Submatrix<_Matrix, AlignedTrait, MatrixIteratorTypes::Row>
{
    public:
 
	typedef _Matrix Matrix;
        typedef Submatrix<Matrix, AlignedTrait, MatrixIteratorTypes::Row> Self_t;

	typedef typename Matrix::SubmatrixType SubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstSubmatrixType;
	typedef typename SubSubmatrixTraits<Matrix, AlignedTrait>::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename SubSubmatrixTraits<Matrix, AlignedTrait>::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = _Matrix::rowAlign;
	static const size_t colAlign = _Matrix::colAlign;

	typedef Matrix ContainerType;

	typedef typename Matrix::IteratorType IteratorType; 
	typedef typename Matrix::StorageType StorageType; 
    
	typedef typename Matrix::Element Element;

	typedef SubmatrixRowColIteratorPT<typename Matrix::RowIterator, Element, AlignedTrait, RCMutableIteratorTag> RowIterator;
	typedef SubmatrixRowColIteratorPT<typename Matrix::ConstRowIterator, Element, AlignedTrait, RCConstIteratorTag> ConstRowIterator;
	typedef typename SubmatrixRowTraits<typename Matrix::Row, Element, AlignedTrait, RCMutableIteratorTag>::SubvectorType Row;
	typedef typename SubmatrixRowTraits<typename Matrix::ConstRow, Element, AlignedTrait, RCConstIteratorTag>::SubvectorType ConstRow;

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

	template <class M2, class AT2>
	Submatrix (const Submatrix<M2, AT2, MatrixIteratorTypes::Row> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	Submatrix (const Submatrix &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	template <class M2>
	Submatrix (const Submatrix<M2, AlignedTrait, MatrixIteratorTypes::Row> &SM)
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

	RowIterator rowBegin ()            { return RowIterator (_M->rowBegin () + _beg_row, _beg_col, _end_col); }
	RowIterator rowEnd   ()            { return RowIterator (_M->rowBegin () + _end_row, _beg_col, _end_col); }
	ConstRowIterator rowBegin () const { return ConstRowIterator (_M->rowBegin () + _beg_row, _beg_col, _end_col); }
	ConstRowIterator rowEnd   () const { return ConstRowIterator (_M->rowBegin () + _end_row, _beg_col, _end_col); }

	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, typename Matrix::Row>::RepresentationType> RawIterator;
	typedef RawIterator ConstRawIterator;
   
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, typename Matrix::Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:

	template <class M2, class AT2, class T2>
	friend class Submatrix;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

/// Specialisation for matrices indexed only by columns

template <class _Matrix, class AlignedTrait>
class Submatrix<_Matrix, AlignedTrait, MatrixIteratorTypes::Col>
{
    public:
 
	typedef _Matrix Matrix;
        typedef Submatrix<Matrix, AlignedTrait, MatrixIteratorTypes::Col> Self_t;

	typedef typename Matrix::SubmatrixType SubmatrixType;
	typedef typename Matrix::ConstSubmatrixType ConstSubmatrixType;
	typedef typename SubSubmatrixTraits<Matrix, AlignedTrait>::AlignedSubmatrixType AlignedSubmatrixType;
	typedef typename SubSubmatrixTraits<Matrix, AlignedTrait>::ConstAlignedSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = _Matrix::rowAlign;
	static const size_t colAlign = _Matrix::colAlign;

	typedef Matrix ContainerType;

	typedef typename Matrix::IteratorType IteratorType; 
	typedef typename Matrix::StorageType StorageType; 
    
	typedef typename Matrix::Element Element;

	typedef SubmatrixRowColIteratorPT<typename Matrix::ColIterator, Element, AlignedTrait, RCMutableIteratorTag> ColIterator;
	typedef SubmatrixRowColIteratorPT<typename Matrix::ConstColIterator, Element, AlignedTrait, RCConstIteratorTag> ConstColIterator;
	typedef typename SubmatrixRowTraits<typename Matrix::Col, Element, AlignedTrait, RCMutableIteratorTag>::SubvectorType Col;
	typedef typename SubmatrixRowTraits<typename Matrix::ConstCol, Element, AlignedTrait, RCConstIteratorTag>::SubvectorType ConstCol;

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

	template <class M2, class AT2>
	Submatrix (const Submatrix<M2, AT2, MatrixIteratorTypes::Col> &SM,
		   size_t row,
		   size_t col,
		   size_t rowdim,
		   size_t coldim)
		: _M (SM._M), _beg_row (SM._beg_row + row), _end_row (SM._beg_row + row + rowdim), _beg_col (SM._beg_col + col), _end_col (SM._beg_col + col + coldim)
		{}

	Submatrix (const Submatrix &SM)
		: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col) {}

	template <class M2>
	Submatrix (const Submatrix<M2, AlignedTrait, MatrixIteratorTypes::Col> &SM)
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

	ColIterator colBegin ()            { return ColIterator (_M->colBegin () + _beg_col, _beg_row, _end_row); }
	ColIterator colEnd   ()            { return ColIterator (_M->colBegin () + _end_col, _beg_row, _end_row); }
	ConstColIterator colBegin () const { return ConstColIterator (_M->colBegin () + _beg_col, _beg_row, _end_row); }
	ConstColIterator colEnd   () const { return ConstColIterator (_M->colBegin () + _end_col, _beg_row, _end_row); }

	typedef MatrixRawIterator<ConstColIterator, typename ElementVectorTraits<Element, Column>::RepresentationType> RawIterator;
	typedef RawIterator ConstRawIterator;
   
	ConstRawIterator rawBegin () const { return ConstRawIterator (colBegin (), 0, colEnd (), rowdim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (colEnd (), 0, colEnd (), rowdim ()); }

	typedef MatrixRawIndexedIterator<ConstColIterator, typename ElementVectorTraits<Element, Column>::RepresentationType, true> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (colBegin (), 0, colEnd (), rowdim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (colEnd (), coldim (), colEnd (), rowdim ()); }

	Matrix &parent () { return *_M; }
	const Matrix &parent () const { return *_M; }

	inline size_t startRow () const { return _beg_row; }
	inline size_t startCol () const { return _beg_col; }

    protected:

	template <class M2, class AT2, class T2>
	friend class Submatrix;

	Matrix *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

} // namespace LELA

#endif // __LELA_MATRIX_SUBMATRIX_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
