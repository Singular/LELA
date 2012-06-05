/* lela/matrix/interface.h
 * Copyright 2001 B. David Saunders,
 *           2001-2002 Bradford Hovinen,
 *           2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * Borrowed from dense-base.h by Bradford Hovinen
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 *
 * This holds the "directly represented" matrix interface. It is provided here
 * only for reference; it does not provide any useful functionality. See the
 * other headers in this directory for useful classes.
 *
 * --------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_matrix_interface_H
#define __LELA_matrix_interface_H

#include <iostream>
#include <vector>
#include <fstream>

namespace LELA
{

/** Matrix interface
 *
 * This defines the interface which matrices should satisfy. It is
 * provided for documentation only.
 *
 * \ingroup matrix
 */
  
template <class _Element>
class MatrixInterface
{
    public:

	/** @name Type-definitions
	 */

	//@{

	typedef _Element Element;

	/// What iterators the matrix supports; see @ref MatrixIteratorTypes
	typedef typename MatrixIteratorTypes::Generic IteratorType;

	/// How the matrix is stores; see @ref MatrixStorageTypes
	typedef typename MatrixStorageTypes::Generic StorageType;

	/** Arbitrary writeable submatrix
	 *
	 * If an implementation cannot allow arbitrary writeable
	 * submatrices, it may make this type equal to @ref
	 * ConstSubmatrixType. For formal reasons it is not possible
	 * to omit it.
	 */
	class SubmatrixType;

	/// Arbitrary const submatrix
	class ConstSubmatrixType;

	/** Writeable aligned submatrix
	 *
	 * This type may be faster than @ref SubmatrixType with the
	 * cost that the choice of beginning and ending row-
	 * respectively column-indices may be restricted. Thus it
	 * should be used in for example recursive algorithms where
	 * the choice of submatrix is free to the developer.
	 *
	 * The implementation is permitted to specify that the
	 * beginning row- respectively column-indices be multiples of
	 * @ref rowAlign respectively @ref colAlign and the the ending
	 * indices be eiter coincident to the end of the matrix or
	 * multiples of these values.
	 */
	class AlignedSubmatrixType;

	/// Const version of @ref AlignedSubmatrixType
	class ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = 1;

	/** Container type
	 *
	 * This should specify a type with the same underlying
	 * storage-representation (e.g. dense, sparse, hybrid) as this
	 * matrix. This can then be used in generic code to construct
	 * a temporary matrix with the same underlying storage-type as
	 * some parametrised input.
	 *
	 * For most matrix-types, it suffices to set ContainerType to
	 * the matrix-type itself. For submatrices and transpose
	 * matrices, ContainerType should be the parent-matrix.
	 */
	class ContainerType;

	/// Row-vector-type
	class Row;
	class ConstRow;  

	/// Column-vector-type
	class Col;
	class ConstCol;

	/// Synonymous with @ref Col resp. @ref ConstCol
	typedef Col Column;
	typedef ConstCol ConstColumn;

	/* N.B. A matrix type may omit either one, but not both, of the
	 * following two iterator types. If one type is omitted, then certain
	 * restrictions on matrix-matrix arithmetic apply; see
	 * @ref{MatrixDomain}
	 */

	/** @name Iterator over rows
	 *
	 * Traverses the rows of the matrix in ascending
	 * order. Dereferencing the iterator yields a row-vector.
	 */

	class RowIterator;    
	class ConstRowIterator;

	/** @name Iterator over columns
	 *
	 * Traverses the columns of the matrix in ascending
	 * order. Dereferencing the iterator yields a column-vector.
	 */

	class ColIterator;
	class ConstColIterator;

	//@}

	/** Default constructor.
	 *
	 * Builds empty 0x0 matrix.
	 */
	MatrixInterface ();

	/** Constructor with size
	 * @param  m  row dimension
	 * @param  n  column dimension
	 *
	 * An implementation may omit this constructor if it makes no
	 * sense to provide it (e.g. with @ref Submatrix or @ref
	 * TransposeMatrix)
	 */
	MatrixInterface (size_t m, size_t n);

	/** Copy constructor
	 */
	MatrixInterface (const MatrixInterface &M);

	/** Constructor from a @ref VectorStream
	 *
	 * The @ref VectorStream supplies either row- or
	 * column-vectors for the matrix, depending on the
	 * implementation. Its parameters also specify the size of the
	 * matrix.
	 *
	 * @param vs @ref VectorStream from which to get vectors
	 */
	MatrixInterface (VectorStream<Row> &vs)
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	/** Operator =
	 */
	MatrixInterface& operator= (const MatrixInterface& M);

	/** Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	size_t rowdim () const;

	/** Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	size_t coldim () const;

	/** \brief Resize the matrix to the given dimensions
	 *
	 * The state of the matrix's entries after a call to this method is
	 * undefined.
	 *
	 * An implementation may omit this interface if it makes no
	 * sense (e.g. for @ref Submatrix or @ref TransposeMatrix)
	 *
	 * @param m Number of rows
	 * @param n Number of columns
	 */
	void resize (size_t m, size_t n);

	/** @name Access to matrix elements
	 */

	//@{

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij);

	/** Remove the entry at (i,j) from the matrix.
	 *
	 * This only has an effect if the matrix has a sparse
	 * representation, in which case it is equivalent to setting
	 * the entry at (i,j) to 0. If the matrix has a dense
	 * representation, then this function need do nothing --
	 * including setting the entry to 0.
	 *
	 * Thus, when an entry is to be set to 0 in generic code, it
	 * is necessary to call both @ref setEntry (i, j, zero) and
	 * eraseEntry (i, j).
	 *
	 * @param i
	 * @param j
	 */
	void eraseEntry (size_t i, size_t j);

	/** Copy the (i, j) entry into x, and return whether the entry exists in the matrix
	 *
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return true if the entry exists in the matrix, false otherwise
	 */
	bool getEntry (Element &x, size_t i, size_t j) const;

	/// Obtain @ref RowIterator for first row of the matrix
	RowIterator rowBegin ();  

	/// Obtain @ref RowIterator for one past the last row of the matrix
	RowIterator rowEnd ();

	/// Obtain @ref ConstRowIterator for first row of the matrix
	ConstRowIterator rowBegin () const;        

	/// Obtain @ref ConstRowIterator for one past the last row of the matrix
	ConstRowIterator rowEnd () const;

	/// Obtain @ref ColIterator for first column of the matrix
	ColIterator colBegin ();

	/// Obtain @ref ColIterator for one past the last column of the matrix
	ColIterator colEnd ();

	/// Obtain @ref ConstColIterator for first column of the matrix
	ConstColIterator colBegin () const;    

	/// Obtain @ref ConstColIterator for one past the last column of the matrix
	ConstColIterator colEnd () const;

	/** @name Raw iterator
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

	/** @name Raw Indexed iterator
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

        class RawIndexedIterator;
        typedef const RawIndexedIterator ConstRawIndexedIterator;

        RawIndexedIterator rawIndexedBegin();
        RawIndexedIterator rawIndexedEnd();   
	ConstRawIndexedIterator rawIndexedBegin() const;
        ConstRawIndexedIterator rawIndexedEnd() const;   
    
	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 *
	 * This may be omitted by an implementation if no Row type is available
	 *
	 * @param i Row index
	 */
	Row operator[] (size_t i);
	ConstRow operator[] (size_t i) const;

	//@}

	/** @name Computing matrix information
	 */

	//@{

	/** Compute the column density, i.e. the number of entries per column
	 */
	template <class Vector>
	Vector &columnDensity (Vector &v) const;

	/** Compute the transpose
	 */
	MatrixInterface &transpose (MatrixInterface &M) const;

	//@}

    protected:

	std::vector<Element>  _rep;
	size_t                _rows, _cols;
};

} // namespace LELA

#endif // __LELA_matrix_interface_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
