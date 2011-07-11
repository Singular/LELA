/* linbox/matrix/transpose-submatrix.h
 * Copyright 2002 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *
 * Evolved from dense-base.h by Bradford Hovinen
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_MATRIX_TRANSPOSE_SUBMATRIX_H
#define __LINBOX_MATRIX_TRANSPOSE_SUBMATRIX_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/vector/stream.h"
#include "linbox/matrix/traits.h"
#include "linbox/matrix/transpose.h"

#undef _A

namespace LinBox
{

/** Submatrix of a transposed matrix
 * 
 * This is the submatrix-type for Transpose, making it possible to
 * take submatrices of transposed matrices.
 */
  
template <class Matrix, class Trait>
class TransposeSubmatrix
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, Trait> Self_t;

	typedef MatrixIteratorTypes::Generic IteratorType;
	typedef MatrixStorageTypes::Generic StorageType;

	typedef typename Matrix::SubmatrixType::ColIterator RowIterator;
	typedef typename Matrix::SubmatrixType::RowIterator ColIterator;
	typedef typename Matrix::SubmatrixType::ConstColIterator ConstRowIterator;
	typedef typename Matrix::SubmatrixType::ConstRowIterator ConstColIterator;

	typedef typename Matrix::SubmatrixType::Row Column;
	typedef typename Matrix::SubmatrixType::Row Col;
	typedef typename Matrix::SubmatrixType::Col Row;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef typename Matrix::ContainerType ContainerType;

	/** Constructor.
	 * @param  A  Underlying matrix of which to construct the transpose
	 */
	TransposeSubmatrix (const TransposeMatrix<Matrix> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}
	
	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeSubmatrix<M2, Trait> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M)
		: _A (M._A)
	{}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { _A.setEntry (j, i, a_ij); }
	inline void eraseEntry (size_t i, size_t j) { _A.eraseEntry (j, i); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline RowIterator rowBegin () { return _A.colBegin (); }
	inline RowIterator rowEnd () { return _A.colEnd (); }
	inline ConstRowIterator rowBegin () const { return _A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return _A.colEnd (); }

	inline ColIterator colBegin () { return _A.rowBegin (); }
	inline ColIterator colEnd () { return _A.rowEnd (); }
	inline ConstColIterator colBegin () const { return _A.rowBegin (); }
	inline ConstColIterator colEnd () const { return _A.rowEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

    protected:

	typename Matrix::SubmatrixType _A;
};

// Specialization for matrices that have both row and column iterators

template <class Matrix>
class TransposeSubmatrix<Matrix, MatrixIteratorTypes::RowCol>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, MatrixIteratorTypes::RowCol> Self_t;
	typedef MatrixIteratorTypes::RowCol IteratorType;
	typedef typename Matrix::StorageType StorageType;

	typedef typename Matrix::SubmatrixType::ColIterator RowIterator;
	typedef typename Matrix::SubmatrixType::RowIterator ColIterator;
	typedef typename Matrix::SubmatrixType::ConstColIterator ConstRowIterator;
	typedef typename Matrix::SubmatrixType::ConstRowIterator ConstColIterator;

	typedef typename Matrix::SubmatrixType::Row Column;
	typedef typename Matrix::SubmatrixType::Row Col;
	typedef typename Matrix::SubmatrixType::Col Row;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::SubmatrixType::ConstRow ConstColumn;
	typedef typename Matrix::SubmatrixType::ConstRow ConstCol;
	typedef typename Matrix::SubmatrixType::ConstCol ConstRow;

	typedef TransposeSubmatrix<Self_t> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeSubmatrix (const TransposeMatrix<Matrix> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeSubmatrix<M2, MatrixIteratorTypes::RowCol> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M) : _A (M._A) {}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { _A.setEntry (j, i, a_ij); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline RowIterator rowBegin () { return _A.colBegin (); }
	inline RowIterator rowEnd () { return _A.colEnd (); }
	inline ConstRowIterator rowBegin () const { return _A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return _A.colEnd (); }

	inline ColIterator colBegin () { return _A.rowBegin (); }
	inline ColIterator colEnd () { return _A.rowEnd (); }
	inline ConstColIterator colBegin () const { return _A.rowBegin (); }
	inline ConstColIterator colEnd () const { return _A.rowEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

    protected:

	typename Matrix::SubmatrixType _A;
};

// Specialization for matrices that have only row iterators

template <class Matrix>
class TransposeSubmatrix<Matrix, MatrixIteratorTypes::Row>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, MatrixIteratorTypes::Row> Self_t;
	typedef MatrixIteratorTypes::Col IteratorType;
	typedef typename Matrix::StorageType StorageType;

	typedef typename Matrix::SubmatrixType::RowIterator ColIterator;
	typedef typename Matrix::SubmatrixType::ConstRowIterator ConstColIterator;

	typedef typename Matrix::SubmatrixType::Row Column;
	typedef typename Matrix::SubmatrixType::Row Col;

	typedef typename Matrix::SubmatrixType::ConstRow ConstColumn;
	typedef typename Matrix::SubmatrixType::ConstRow ConstCol;

	typedef MatrixRawIterator<ColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType, true> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeSubmatrix (const TransposeMatrix<Matrix> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeSubmatrix<M2, MatrixIteratorTypes::Row> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M) : _A (M._A) {}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { _A.setEntry (j, i, a_ij); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline ColIterator colBegin () { return _A.rowBegin (); }
	inline ColIterator colEnd () { return _A.rowEnd (); }
	inline ConstColIterator colBegin () const { return _A.rowBegin (); }
	inline ConstColIterator colEnd () const { return _A.rowEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (colBegin (), 0, colEnd (), rowdim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (colEnd (), 0, colEnd (), rowdim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (colBegin (), 0, colEnd (), rowdim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (colEnd (), 0, colEnd (), rowdim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (colBegin (), 0, colEnd (), rowdim ()); }
	inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (colEnd (), coldim (), colEnd (), rowdim ()); }

    protected:

	typename Matrix::SubmatrixType _A;
};

// Specialization for matrices that have only column iterators

template <class Matrix>
class TransposeSubmatrix<Matrix, MatrixIteratorTypes::Col>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, MatrixIteratorTypes::Col> Self_t;
	typedef MatrixIteratorTypes::Row IteratorType;
	typedef typename Matrix::StorageType StorageType;

	typedef typename Matrix::SubmatrixType::ColIterator RowIterator;
	typedef typename Matrix::SubmatrixType::ConstColIterator ConstRowIterator;

	typedef typename Matrix::SubmatrixType::Col Row;
	typedef typename Matrix::SubmatrixType::ConstCol ConstRow;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeSubmatrix (TransposeMatrix<Matrix> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (A._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (TransposeSubmatrix<M2, MatrixIteratorTypes::Col> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: _A (M._A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M) : _A (M._A) {}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { _A.setEntry (j, i, a_ij); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline RowIterator rowBegin () { return _A.colBegin (); }
	inline RowIterator rowEnd () { return _A.colEnd (); }
	inline ConstRowIterator rowBegin () const { return _A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return _A.colEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

    protected:

	typename Matrix::SubmatrixType _A;
};

} // namespace LinBox

#endif // __LINBOX_MATRIX_TRANSPOSE_SUBMATRIX_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
