/* lela/matrix/transpose.h
 * Copyright 2002 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Transposed matrix
 *
 * --------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_TRANSPOSE_H
#define __LELA_MATRIX_TRANSPOSE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "lela/vector/stream.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/raw-iterator.h"

namespace LELA
{

// Forward declaration
template <class Matrix, class Submatrix, class Trait = typename Matrix::IteratorType>
class TransposeSubmatrix;

// Struct to determine which StorageType to use
template <class ST> struct TransposeStorageType { typedef ST StorageType; };
template <> struct TransposeStorageType<MatrixStorageTypes::Dense> { typedef MatrixStorageTypes::DenseTranspose StorageType; };
template <> struct TransposeStorageType<MatrixStorageTypes::DenseTranspose> { typedef MatrixStorageTypes::Generic StorageType; };
template <> struct TransposeStorageType<MatrixStorageTypes::M4RI> { typedef MatrixStorageTypes::M4RITranspose StorageType; };
template <> struct TransposeStorageType<MatrixStorageTypes::M4RITranspose> { typedef MatrixStorageTypes::Generic StorageType; };

/** Transpose of a matrix
 * 
 * This class takes a matrix meeting the @ref{DenseMatrixBase} archetype and
 * switches the row and column iterators, giving the transpose of the original
 * matrix. It is generic with respect to the matrix given.
 * 
 * If the matrix given has limited iterators, then its transpose will have
 * limited iterators as well. In particular, if the matrix given has only row
 * iterators, then the transpose object will have only column iterators, and
 * vice versa.
 * 
 * This class differs from @ref{Transpose} in that it constructs a full matrix
 * representation, with row and/or column iterators. It does not include any
 * logic for matrix-vector products, and does not meet the
 * @ref{BlackboxArchetype} interface. Nor does it make such assumptions about
 * the matrix given.
 *
 * This class gives a constant matrix as output. It provides no iterators for
 * modification of the data in the matrix.
 * 
 * The input/output functionality of this class passes requests directly through
 * to the underlying matrix. In particular, the output will be the transpose of
 * the matrix expected and the input will expect the transpose of the matrix
 * given. Thus, it is not recommended to use TransposeMatrix for reading and
 * writing matrices, except for testing purposes.
 *
 * \ingroup matrix
 */
  
template <class Matrix, class Trait = typename Matrix::IteratorType>
class TransposeMatrix
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeMatrix<Matrix, Trait> Self_t;

	typedef MatrixIteratorTypes::Generic IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Matrix::ColIterator RowIterator;
	typedef typename Matrix::RowIterator ColIterator;
	typedef typename Matrix::ConstColIterator ConstRowIterator;
	typedef typename Matrix::ConstRowIterator ConstColIterator;

	typedef typename Matrix::Row Column;
	typedef typename Matrix::Row Col;
	typedef typename Matrix::Col Row;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t, typename Matrix::SubmatrixType> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef typename Matrix::ContainerType ContainerType;

	/** Constructor.
	 * @param  A  Underlying matrix of which to construct the transpose
	 */
	TransposeMatrix (Matrix &__A)
		: A (__A)
	{}

	TransposeMatrix (const TransposeMatrix &M)
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

	friend class TransposeSubmatrix<Matrix, typename Matrix::SubmatrixType, Trait>;
	friend class TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType, Trait>;

	Matrix &A;
};

// Specialization for matrices that have both row and column iterators

template <class Matrix>
class TransposeMatrix<Matrix, MatrixIteratorTypes::RowCol>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeMatrix<Matrix, MatrixIteratorTypes::RowCol> Self_t;
	typedef MatrixIteratorTypes::RowCol IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Matrix::ColIterator RowIterator;
	typedef typename Matrix::RowIterator ColIterator;
	typedef typename Matrix::ConstColIterator ConstRowIterator;
	typedef typename Matrix::ConstRowIterator ConstColIterator;

	typedef typename Matrix::Row Column;
	typedef typename Matrix::Row Col;
	typedef typename Matrix::Col Row;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::ConstRow ConstColumn;
	typedef typename Matrix::ConstRow ConstCol;
	typedef typename Matrix::ConstCol ConstRow;

	typedef TransposeSubmatrix<Self_t, typename Matrix::SubmatrixType> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeMatrix (Matrix &__A) : A (__A) {}
	TransposeMatrix (const TransposeMatrix &M) : A (M.A) {}

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

	friend class TransposeSubmatrix<Matrix, typename Matrix::SubmatrixType, MatrixIteratorTypes::RowCol>;
	friend class TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType, MatrixIteratorTypes::RowCol>;

	Matrix &A;
};

// Specialization for matrices that have only row iterators

template <class Matrix>
class TransposeMatrix<Matrix, MatrixIteratorTypes::Row>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeMatrix<Matrix, MatrixIteratorTypes::Row> Self_t;
	typedef MatrixIteratorTypes::Col IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Matrix::RowIterator ColIterator;
	typedef typename Matrix::ConstRowIterator ConstColIterator;

	typedef typename Matrix::Row Column;
	typedef typename Matrix::Row Col;

	typedef MatrixRawIterator<ColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstColIterator, typename ElementVectorTraits<Element, Col>::RepresentationType, true> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::ConstRow ConstColumn;
	typedef typename Matrix::ConstRow ConstCol;

	typedef TransposeSubmatrix<Self_t, typename Matrix::SubmatrixType> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	//TransposeMatrix () {}
	TransposeMatrix (Matrix &__A) : A (__A) {}
	TransposeMatrix (const TransposeMatrix &M) : A (M.A) {}

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

	friend class TransposeSubmatrix<Matrix, typename Matrix::SubmatrixType, MatrixIteratorTypes::Row>;
	friend class TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType, MatrixIteratorTypes::Row>;

	Matrix &A;
};

// Specialization for matrices that have only column iterators

template <class Matrix>
class TransposeMatrix<Matrix, MatrixIteratorTypes::Col>
{
    public:

	typedef typename Matrix::Element Element;
	typedef TransposeMatrix<Matrix, MatrixIteratorTypes::Col> Self_t;
	typedef MatrixIteratorTypes::Row IteratorType;
	typedef typename TransposeStorageType<typename Matrix::StorageType>::StorageType StorageType;

	typedef typename Matrix::ColIterator RowIterator;
	typedef typename Matrix::ConstColIterator ConstRowIterator;

	typedef typename Matrix::Col Row;
	typedef typename Matrix::ConstCol ConstRow;

	typedef MatrixRawIterator<RowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> RawIterator;
	typedef MatrixRawIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType> ConstRawIterator;
	typedef MatrixRawIndexedIterator<ConstRowIterator, typename ElementVectorTraits<Element, Row>::RepresentationType, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	typedef TransposeSubmatrix<Self_t, typename Matrix::SubmatrixType> SubmatrixType;
	typedef TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType> ConstSubmatrixType;
	typedef SubmatrixType AlignedSubmatrixType;
	typedef ConstSubmatrixType ConstAlignedSubmatrixType;

	static const size_t rowAlign = Matrix::colAlign;
	static const size_t colAlign = Matrix::rowAlign;

	typedef Matrix ContainerType;

	TransposeMatrix (Matrix &__A) : A (__A) {}
	TransposeMatrix (const TransposeMatrix &M) : A (M.A) {}

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

	friend class TransposeSubmatrix<Matrix, typename Matrix::SubmatrixType, MatrixIteratorTypes::Col>;
	friend class TransposeSubmatrix<const Matrix, typename Matrix::ConstSubmatrixType, MatrixIteratorTypes::Col>;

	const Matrix &A;
};

} // namespace LELA

#include "lela/matrix/transpose-submatrix.h"

#endif // __LELA_MATRIX_TRANSPOSE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
