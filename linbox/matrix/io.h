/* linbox/matrix/io.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_MATRIX_IO_H
#define __LINBOX_MATRIX_IO_H

#include <iostream>

namespace LinBox
{

/// File-formats for matrix-output
enum FileFormatTag {
	FORMAT_DETECT, FORMAT_UNKNOWN, FORMAT_TURNER, FORMAT_ONE_BASED, FORMAT_GUILLAUME, FORMAT_MAPLE, FORMAT_MATLAB, FORMAT_SAGE, FORMAT_PRETTY
};

/// Exception thrown when the format cannot be detected
class UnrecognisedFormat {};

/// Exception class for invalid matrix input
class InvalidMatrixInput {};

/// Exception class for functions which haven't been implemented
class NotImplemented {};

/// Class to read a matrix from an istream
template <class Field>
class MatrixReader {
	const Field &_F;

public:
	/// Construct a new MatrixReader using the field F for element-input
	MatrixReader (const Field &F) : _F (F) {}

	/// Read a matrix from the istream is and place the data in the matrix A
	template <class Matrix>
	std::istream &read (std::istream &is, Matrix &A, FileFormatTag format = FORMAT_DETECT) const;

private:
	FileFormatTag detectFormat (std::istream &is) const;

	template <class Matrix>
	std::istream &readTurner (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readOneBased (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readGuillaume (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readMaple (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readMatlab (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readSage (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readPretty (std::istream &is, Matrix &A) const;

	template <class Vector>
	void appendEntry (Vector &v, size_t index, const typename Field::Element &a) const
		{ appendEntrySpecialised (v, index, a, typename VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::DenseVectorTag) const;

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::SparseVectorTag) const
		{ if (!_F.isZero (a)) v.push_back (typename std::iterator_traits<typename Vector::iterator>::value_type (index, a)); }

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::DenseZeroOneVectorTag) const;

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::SparseZeroOneVectorTag) const
		{ if (!_F.isZero (a)) v.push_back (index); }

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::HybridZeroOneVectorTag) const;
};

template <class Field>
class MatrixWriter {
	const Field &_F;

public:
	/// Construct a new MatrixWriter using the field F for element-output
	MatrixWriter (const Field &F) : _F (F) {}

	template <class Matrix>
	std::ostream &write (std::ostream &os, const Matrix &A, FileFormatTag format = FORMAT_PRETTY) const;

private:
	template <class Matrix>
	std::ostream &writeTurner (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeOneBased (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeGuillaume (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeMaple (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeMatlab (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeSage (std::ostream &os, const Matrix &A) const;	

	template <class Matrix>
	std::ostream &writePretty (std::ostream &os, const Matrix &A) const;
};

} // namespace LinBox

#include "linbox/matrix/io.inl"

#endif // __LINBOX_MATRIX_IO_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
