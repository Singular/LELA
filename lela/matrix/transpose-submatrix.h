/* lela/matrix/transpose-submatrix.h
 * Copyright 2002 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *
 * Submatrix of a transposed matrix
 *
 * --------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_TRANSPOSE_SUBMATRIX_H
#define __LELA_MATRIX_TRANSPOSE_SUBMATRIX_H

#include <iostream>
#include <vector>
#include <fstream>

#include "lela/vector/stream.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/transpose.h"

namespace LELA
{

/** Submatrix of a transposed matrix
 * 
 * This is the submatrix-type for Transpose, making it possible to
 * take submatrices of transposed matrices.
 *
 * \ingroup matrix
 */
  
template <class Matrix, class Submatrix, class Trait>
class TransposeSubmatrix
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, Submatrix, Trait> Self_t;

	typedef MatrixIteratorTypes::Generic IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Submatrix::ColIterator RowIterator;
	typedef typename Submatrix::RowIterator ColIterator;
	typedef typename Submatrix::ConstColIterator ConstRowIterator;
	typedef typename Submatrix::ConstRowIterator ConstColIterator;

	typedef typename Submatrix::Row Column;
	typedef typename Submatrix::Row Col;
	typedef typename Submatrix::Col Row;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t, Submatrix> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef typename Matrix::ContainerType ContainerType;

	/** Constructor.
	 * @param  A  Underlying matrix of which to construct the transpose
	 */
	TransposeSubmatrix (TransposeMatrix<Matrix> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}
	
	template <class M2>
	TransposeSubmatrix (TransposeMatrix<M2> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeSubmatrix<M2, Trait> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M)
		: A (M.A)
	{}

	inline size_t rowdim () const { return A.coldim (); }
	inline size_t coldim () const { return A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { A.setEntry (j, i, a_ij); }
	inline void eraseEntry (size_t i, size_t j) { A.eraseEntry (j, i); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return A.getEntry (x, j, i); }

	inline RowIterator rowBegin () { return A.colBegin (); }
	inline RowIterator rowEnd () { return A.colEnd (); }
	inline ConstRowIterator rowBegin () const { return A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return A.colEnd (); }

	inline ColIterator colBegin () { return A.rowBegin (); }
	inline ColIterator colEnd () { return A.rowEnd (); }
	inline ConstColIterator colBegin () const { return A.rowBegin (); }
	inline ConstColIterator colEnd () const { return A.rowEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	inline const Matrix &parent () const { return A; }

    protected:

	Submatrix A;
};

// Specialization for matrices that have both row and column iterators

template <class Matrix, class Submatrix>
class TransposeSubmatrix<Matrix, Submatrix, MatrixIteratorTypes::RowCol>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, Submatrix, MatrixIteratorTypes::RowCol> Self_t;
	typedef MatrixIteratorTypes::RowCol IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Submatrix::ColIterator RowIterator;
	typedef typename Submatrix::RowIterator ColIterator;
	typedef typename Submatrix::ConstColIterator ConstRowIterator;
	typedef typename Submatrix::ConstRowIterator ConstColIterator;

	typedef typename Submatrix::Row Column;
	typedef typename Submatrix::Row Col;
	typedef typename Submatrix::Col Row;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef typename Submatrix::ConstRow ConstColumn;
	typedef typename Submatrix::ConstRow ConstCol;
	typedef typename Submatrix::ConstCol ConstRow;

	typedef TransposeSubmatrix<Self_t, Submatrix> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeSubmatrix (TransposeMatrix<Matrix> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (TransposeMatrix<M2> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeMatrix<Matrix> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeSubmatrix<M2, MatrixIteratorTypes::RowCol> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M) : A (M.A) {}

	inline size_t rowdim () const { return A.coldim (); }
	inline size_t coldim () const { return A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { A.setEntry (j, i, a_ij); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return A.getEntry (x, j, i); }

	inline RowIterator rowBegin () { return A.colBegin (); }
	inline RowIterator rowEnd () { return A.colEnd (); }
	inline ConstRowIterator rowBegin () const { return A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return A.colEnd (); }

	inline ColIterator colBegin () { return A.rowBegin (); }
	inline ColIterator colEnd () { return A.rowEnd (); }
	inline ConstColIterator colBegin () const { return A.rowBegin (); }
	inline ConstColIterator colEnd () const { return A.rowEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	inline const Matrix &parent () const { return A; }

    protected:

	Submatrix A;
};

// Specialization for matrices that have only row iterators

template <class Matrix, class Submatrix>
class TransposeSubmatrix<Matrix, Submatrix, MatrixIteratorTypes::Row>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, Submatrix, MatrixIteratorTypes::Row> Self_t;
	typedef MatrixIteratorTypes::Col IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Submatrix::RowIterator ColIterator;
	typedef typename Submatrix::ConstRowIterator ConstColIterator;

	typedef typename Submatrix::Row Column;
	typedef typename Submatrix::Row Col;

	typedef typename Submatrix::ConstRow ConstColumn;
	typedef typename Submatrix::ConstRow ConstCol;

	typedef MatrixRawIterator<ColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType, true> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t, Submatrix> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeSubmatrix (TransposeMatrix<Matrix> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (TransposeMatrix<M2> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeMatrix<Matrix> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeSubmatrix<M2, MatrixIteratorTypes::Row> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M) : A (M.A) {}

	inline size_t rowdim () const { return A.coldim (); }
	inline size_t coldim () const { return A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { A.setEntry (j, i, a_ij); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return A.getEntry (x, j, i); }

	inline ColIterator colBegin () { return A.rowBegin (); }
	inline ColIterator colEnd () { return A.rowEnd (); }
	inline ConstColIterator colBegin () const { return A.rowBegin (); }
	inline ConstColIterator colEnd () const { return A.rowEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (colBegin (), 0, colEnd (), rowdim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (colEnd (), 0, colEnd (), rowdim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (colBegin (), 0, colEnd (), rowdim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (colEnd (), 0, colEnd (), rowdim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (colBegin (), 0, colEnd (), rowdim ()); }
	inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (colEnd (), coldim (), colEnd (), rowdim ()); }

	inline const Matrix &parent () const { return A; }

    protected:

	Submatrix A;
};

// Specialization for matrices that have only column iterators

template <class Matrix, class Submatrix>
class TransposeSubmatrix<Matrix, Submatrix, MatrixIteratorTypes::Col>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeSubmatrix<Matrix, Submatrix, MatrixIteratorTypes::Col> Self_t;
	typedef MatrixIteratorTypes::Row IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Submatrix::ColIterator RowIterator;
	typedef typename Submatrix::ConstColIterator ConstRowIterator;

	typedef typename Submatrix::Col Row;
	typedef typename Submatrix::ConstCol ConstRow;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t, Submatrix> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeSubmatrix (TransposeMatrix<Matrix> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (const TransposeMatrix<M2> &__A, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (__A.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (TransposeSubmatrix &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	template <class M2>
	TransposeSubmatrix (TransposeSubmatrix<M2, MatrixIteratorTypes::Col> &M, size_t row, size_t col, size_t rowdim, size_t coldim)
		: A (M.A, col, row, coldim, rowdim)
	{}

	TransposeSubmatrix (const TransposeSubmatrix &M) : A (M.A) {}

	inline size_t rowdim () const { return A.coldim (); }
	inline size_t coldim () const { return A.rowdim (); }

	inline void setEntry (size_t i, size_t j, const Element &a_ij) { A.setEntry (j, i, a_ij); }
	inline bool getEntry (Element &x, size_t i, size_t j) const { return A.getEntry (x, j, i); }

	inline RowIterator rowBegin () { return A.colBegin (); }
	inline RowIterator rowEnd () { return A.colEnd (); }
	inline ConstRowIterator rowBegin () const { return A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return A.colEnd (); }

	inline RawIterator      rawBegin ()       { return RawIterator      (rowBegin (), 0, rowEnd (), coldim ()); }
	inline RawIterator      rawEnd ()         { return RawIterator      (rowEnd (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	inline ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        inline ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	inline const Matrix &parent () const { return A; }

    protected:

	Submatrix A;
};

} // namespace LELA

#endif // __LELA_MATRIX_TRANSPOSE_SUBMATRIX_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
