/* lela/matrix/io.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Support for reading and writing matrices to and from streams
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_IO_H
#define __LELA_MATRIX_IO_H

#include <iostream>

#include "lela/lela-config.h"

#ifdef __LELA_HAVE_LIBPNG
#  include <png.h>
#endif

#include "lela/vector/traits.h"
#include "lela/matrix/traits.h"

namespace LELA
{

/// @name Matrix I/O support
///
/// \ingroup matrix
//@{

/// File-formats for matrix-output
///
/// \ingroup matrix
enum FileFormatTag {
	FORMAT_DETECT, FORMAT_UNKNOWN, FORMAT_TURNER, FORMAT_ONE_BASED, FORMAT_DUMAS, FORMAT_MAPLE, FORMAT_MATLAB, FORMAT_SAGE, FORMAT_PRETTY,
#ifdef __LELA_HAVE_LIBPNG
	FORMAT_PNG
#endif // __LELA_HAVE_LIBPNG
};

/// Exception thrown when the data-format of a matrix for reading cannot be detected
///
/// \ingroup matrix
class UnrecognisedFormat {};

/// Exception class for invalid input when reading a matrix
///
/// \ingroup matrix
class InvalidMatrixInput {};

/// Class to read a matrix from an istream
///
/// \ingroup matrix
template <class Ring>
class MatrixReader {
	const Ring &_F;

public:
	/// Construct a new MatrixReader using the ring F for element-input
	MatrixReader (const Ring &F) : _F (F) {}

	/// Read a matrix from the istream is and place the data in the matrix A
	template <class Matrix>
	std::istream &read (std::istream &is, Matrix &A, FileFormatTag format = FORMAT_DETECT) const;

	static FileFormatTag detectFormat (std::istream &is);

private:
	template <class Matrix>
	std::istream &readTurner (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readOneBased (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readDumas (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readMaple (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readMatlab (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readSage (std::istream &is, Matrix &A) const;

	template <class Matrix>
	std::istream &readPretty (std::istream &is, Matrix &A) const;

	static bool isDumas (char *buf, std::streamsize n);
	static bool isTurner (char *buf, std::streamsize n);
	static bool isMaple (char *buf, std::streamsize n);
	static bool isMatlab (char *buf, std::streamsize n);
	static bool isSage (char *buf, std::streamsize n);
	static bool isPretty (char *buf, std::streamsize n);

#ifdef __LELA_HAVE_LIBPNG
	static const unsigned _png_sig_size = 8;

	static bool isPNG (std::istream &is);

	static void PNGReadData (png_structp pngPtr, png_bytep data, png_size_t length);

	template <class Vector>
	void readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
				      VectorRepresentationTypes::Dense) const
		{ throw NotImplemented (); }

	template <class Vector>
	void readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
				      VectorRepresentationTypes::Sparse) const
		{ throw NotImplemented (); }

	template <class Vector>
	void readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
				      VectorRepresentationTypes::Dense01) const;
	template <class Vector>
	void readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
				      VectorRepresentationTypes::Sparse01) const;
	template <class Vector>
	void readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
				      VectorRepresentationTypes::Hybrid01) const;
	template <class Vector>
	void readPNGBlock (Vector &v, png_byte x, size_t start, size_t stop) const
		{ readPNGBlockSpecialised (v, x, start, stop, typename VectorTraits<Ring, Vector>::RepresentationType ()); }

	template <class Vector>
	void readPNGRow (Vector &v, png_bytep row, size_t width) const;

	template <class Matrix>
	std::istream &readPNGSpecialised (std::istream &is, Matrix &A, MatrixIteratorTypes::Row) const;

	template <class Matrix>
	std::istream &readPNGSpecialised (std::istream &is, Matrix &A, MatrixIteratorTypes::Col) const
		{ throw NotImplemented (); }

	template <class Matrix>
	std::istream &readPNGSpecialised (std::istream &is, Matrix &A, MatrixIteratorTypes::RowCol) const
		{ return readPNGSpecialised (is, A, MatrixIteratorTypes::Row ()); }

	template <class Matrix>
	std::istream &readPNG (std::istream &is, Matrix &A) const
		{ return readPNGSpecialised (is, A, typename Matrix::IteratorType ()); }

#endif // __LELA_HAVE_LIBPNG	

	template <class Vector>
	void appendEntry (Vector &v, size_t index, const typename Ring::Element &a) const
		{ appendEntrySpecialised (v, index, a, typename VectorTraits<Ring, Vector>::RepresentationType ()); }

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Dense) const;

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Sparse) const
		{ if (!_F.isZero (a)) v.push_back (typename std::iterator_traits<typename Vector::iterator>::value_type (index, a)); }

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Dense01) const;

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Sparse01) const
		{ if (!_F.isZero (a)) v.push_back (index); }

	template <class Vector>
	void appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Hybrid01) const;
};

/// Class for writing a matrix to a stream
///
/// \ingroup matrix
template <class Ring>
class MatrixWriter {
	const Ring &_F;

public:
	/// Construct a new MatrixWriter using the ring F for element-output
	MatrixWriter (const Ring &F) : _F (F) {}

	template <class Matrix>
	std::ostream &write (std::ostream &os, const Matrix &A, FileFormatTag format = FORMAT_PRETTY) const;

private:
	template <class Matrix>
	std::ostream &writeTurner (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeOneBased (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeDumas (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeMaple (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeMatlab (std::ostream &os, const Matrix &A) const;

	template <class Matrix>
	std::ostream &writeSage (std::ostream &os, const Matrix &A) const;	

	template <class Matrix>
	std::ostream &writePretty (std::ostream &os, const Matrix &A) const;

#ifdef __LELA_HAVE_LIBPNG
	static void PNGWriteData (png_structp png_ptr, png_bytep data, png_size_t length);
	static void PNGFlush (png_structp png_ptr);

	template <class Vector>
	void copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Dense) const
		{ throw NotImplemented (); }

	template <class Vector>
	void copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Sparse) const
		{ throw NotImplemented (); }

	template <class Vector>
	void copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Dense01) const;

	template <class Vector>
	void copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Sparse01) const;

	template <class Vector>
	void copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Hybrid01) const;

	template <class Vector>
	void copyToPNGData (png_bytep data, const Vector &v, size_t width) const
		{ copyToPNGDataSpecialised (data, v, width, typename VectorTraits<Ring, Vector>::RepresentationType ()); }

	template <class Matrix>
	std::ostream &writePNGSpecialised (std::ostream &is, const Matrix &A, MatrixIteratorTypes::Row) const;

	template <class Matrix>
	std::ostream &writePNGSpecialised (std::ostream &is, const Matrix &A, MatrixIteratorTypes::Col) const
		{ throw NotImplemented (); }

	template <class Matrix>
	std::ostream &writePNGSpecialised (std::ostream &is, const Matrix &A, MatrixIteratorTypes::RowCol) const
		{ return writePNGSpecialised (is, A, MatrixIteratorTypes::Row ()); }

	template <class Matrix>
	std::ostream &writePNG (std::ostream &is, const Matrix &A) const
		{ return writePNGSpecialised (is, A, typename Matrix::IteratorType ()); }
#endif // __LELA_HAVE_LIBPNG
};

//@}

} // namespace LELA

#endif // __LELA_MATRIX_IO_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
