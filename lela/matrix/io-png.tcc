/* lela/matrix/io-png.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Support for reading and writing matrices in PNG-format
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_IO_PNG_TCC
#define __LELA_MATRIX_IO_PNG_TCC

#ifndef __LELA_HAVE_LIBPNG
#  error "This header file requires that LELA be configured with libpng enabled. Please ensure that libpng is properly installed and re-run configure."
#else
#  include <pngconf.h>
#  include <png.h>
#endif

#include "lela/matrix/io.h"
#include "lela/vector/traits.h"
#include "lela/vector/bit-iterator.h"
#include "lela/vector/hybrid.h"

#include <iostream>

/// compression level for z-lib default:(3-6), best: 9
const int COMPRESSION_LEVEL_PNG = 6;

namespace LELA
{

template <class Ring>
bool MatrixReader<Ring>::isPNG (std::istream &is)
{
	png_byte pngsig[_png_sig_size];

	is.read ((char*) pngsig, _png_sig_size);

	if (!is.good ()) return false;

	return png_sig_cmp (pngsig, 0, _png_sig_size) == 0;
}

template <class Ring>
void MatrixReader<Ring>::PNGReadData (png_structp png_ptr, png_bytep data, png_size_t length)
{
	png_voidp a = png_get_io_ptr (png_ptr);
	((std::istream *) a)->read ((char *) data, length);
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::readPNGRow (Vector &v, png_bytep row, size_t width) const
{
	size_t i;

	for (i = 0; i < width / (8 * sizeof (png_byte)); ++i)
		readPNGBlock (v, row[i], i * (8 * sizeof (png_byte)), 8 * sizeof (png_byte));

	if (width % (8 * sizeof (png_byte)))
		readPNGBlock (v, row[i], i * (8 * sizeof (png_byte)), width % (8 * sizeof (png_byte)));
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
						   VectorRepresentationTypes::Dense01) const
{
	typedef typename std::iterator_traits<typename Vector::word_iterator>::value_type word_type;

	size_t count;
	png_byte t;

	typename Vector::word_iterator i = v.word_begin () + (start >> WordTraits<word_type>::logof_size);
	word_type mask = Vector::Endianness::e_j (start & WordTraits<word_type>::pos_mask);

	for (count = 0, t = BigEndian<png_byte>::e_0; count < stop; t >>= 1, mask = Vector::Endianness::shift_right (mask, 1), ++count)
		if (!(x & t))
			*i |= mask;
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
						   VectorRepresentationTypes::Sparse01) const
{
	size_t idx;
	png_byte t;

	for (idx = 0, t = ~(((png_byte) -1) & (((png_byte) -1) >> 1)); idx < stop; ++idx, t >>= 1)
		if (!(x & t))
			v.push_back (start + idx);
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::readPNGBlockSpecialised (Vector &v, png_byte x, size_t start, size_t stop,
						   VectorRepresentationTypes::Hybrid01) const
{
	typename Vector::index_type idx = start >> WordTraits<typename Vector::word_type>::logof_size;
	size_t count;
	png_byte t;

	typename Vector::word_type mask = Vector::Endianness::e_j (start & WordTraits<typename Vector::word_type>::pos_mask);

	for (count = 0, t = ~(((png_byte) -1) & (((png_byte) -1) >> 1)); count < stop; t >>= 1, mask = Vector::Endianness::shift_right (mask, 1), ++count) {
		if (!(x & t)) {
			if (v.empty () || v.back ().first != idx)
				v.push_back (typename Vector::value_type (idx, mask));
			else
				v.back ().second |= mask;
		}
	}
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readPNGSpecialised (std::istream &is, Matrix &A, MatrixIteratorTypes::Row) const
{
	png_structp png_ptr;
	png_infop info_ptr, end_ptr;
	png_bytep data;
	png_bytepp row_pointers;

	png_uint_32 width, height;
	png_uint_32 bit_depth, colour_type;

	typename Matrix::RowIterator i;

	size_t idx, stride;

	// Check that we are actually reading a PNG
	if (!isPNG (is))
		throw InvalidMatrixInput ();

	// Set up data-structures
	png_ptr = png_create_read_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (png_ptr == NULL)
		throw LELAError ("Cannot initialise PNG-library");

	info_ptr = png_create_info_struct (png_ptr);

	if (info_ptr == NULL) {
		png_destroy_read_struct (&png_ptr, NULL, NULL);
		throw LELAError ("Cannot initialise PNG-library");
	}

	end_ptr = png_create_info_struct (png_ptr);

	if (end_ptr == NULL) {
		png_destroy_read_struct (&png_ptr, &info_ptr, NULL);
		throw LELAError ("Cannot initialise PNG-library");
	}

	if (setjmp (png_jmpbuf (png_ptr))) {
		png_destroy_read_struct (&png_ptr, &info_ptr, &end_ptr);
		throw LELAError ("Error reading PNG-data");
	}

	png_set_read_fn (png_ptr, (png_voidp) &is, PNGReadData);

	// Get header-information and make sure PNG is suitable
	png_set_sig_bytes (png_ptr, _png_sig_size); // We've already read the signature
	png_read_info (png_ptr, info_ptr);

	if ((bit_depth = png_get_bit_depth (png_ptr, info_ptr)) != 1) {
		std::ostringstream str;
		str << "Incorrect bit_depth " << bit_depth << " (should be 1)" << std::ends;
		throw LELAError (str.str ().c_str ());
	}

	if ((colour_type = png_get_color_type (png_ptr, info_ptr)) != 0) {
		std::ostringstream str;
		str << "Incorrect colour_type " << colour_type << " (should be 0)" << std::ends;
		throw LELAError (str.str ().c_str ());
	}

	width = png_get_image_width (png_ptr, info_ptr);
	height = png_get_image_height (png_ptr, info_ptr);

	// Allocate data for image
	stride = (width + WordTraits<png_byte>::bits - 1) / WordTraits<png_byte>::bits;
	data = new png_byte[height * stride];
	row_pointers = new png_bytep[height];

	for (idx = 0; idx < height; ++idx)
		row_pointers[idx] = data + idx * stride;

	// Read image from file
	png_read_image (png_ptr, row_pointers);

	// Load image-data into matrix
	A.resize (height, width);

	for (i = A.rowBegin (), idx = 0; i != A.rowEnd (); ++idx, ++i)
		readPNGRow (*i, row_pointers[idx], width);

	// Clean up
	png_destroy_read_struct (&png_ptr, &info_ptr, &end_ptr);

	delete[] data;
	delete[] row_pointers;

	return is;
}

template <class Ring>
void MatrixWriter<Ring>::PNGWriteData (png_structp png_ptr, png_bytep data, png_size_t length)
{
	png_voidp a = png_get_io_ptr (png_ptr);
	((std::ostream *) a)->write ((char *) data, length);
}

template <class Ring>
void MatrixWriter<Ring>::PNGFlush (png_structp png_ptr)
{
	png_voidp a = png_get_io_ptr (png_ptr);
	((std::ostream *) a)->flush ();
}

template <class Ring>
template <class Vector>
void MatrixWriter<Ring>::copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Dense01) const
{
	typename Vector::const_iterator i;
	size_t idx;

	for (i = v.begin (), idx = 0; i != v.end (); ++i, ++idx) {
		if (*i)
			data[idx / WordTraits<png_byte>::bits] &= ~BigEndian<png_byte>::e_j (idx & WordTraits<png_byte>::pos_mask);
		else
			data[idx / WordTraits<png_byte>::bits] |= BigEndian<png_byte>::e_j (idx & WordTraits<png_byte>::pos_mask);
	}
}

template <class Ring>
template <class Vector>
void MatrixWriter<Ring>::copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Sparse01) const
{
	typename Vector::const_iterator i;

	std::fill (data, data + (width + WordTraits<png_byte>::bits - 1) / WordTraits<png_byte>::bits, (png_byte) -1);

	for (i = v.begin (); i != v.end (); ++i)
		data[*i / WordTraits<png_byte>::bits] &= ~BigEndian<png_byte>::e_j (*i & WordTraits<png_byte>::pos_mask);
}

template <class Ring>
template <class Vector>
void MatrixWriter<Ring>::copyToPNGDataSpecialised (png_bytep data, const Vector &v, size_t width, VectorRepresentationTypes::Hybrid01) const
{
	typename Vector::const_iterator i;
	typename Vector::word_type t;
	size_t idx;

	std::fill (data, data + (width + WordTraits<png_byte>::bits - 1) / WordTraits<png_byte>::bits, (png_byte) -1);

	for (i = v.begin (); i != v.end (); ++i) {
		idx = i->first << WordTraits<typename Vector::word_type>::logof_size;

		for (t = Vector::Endianness::e_0; t != 0 && idx < width; t = Vector::Endianness::shift_right (t, 1), ++idx) {
			if (i->second & t)
				data[idx / WordTraits<png_byte>::bits] &= ~BigEndian<png_byte>::e_j (idx & WordTraits<png_byte>::pos_mask);
			else
				data[idx / WordTraits<png_byte>::bits] |= BigEndian<png_byte>::e_j (idx & WordTraits<png_byte>::pos_mask);
		}
	}
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writePNGSpecialised (std::ostream &os, const Matrix &A, MatrixIteratorTypes::Row) const
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_bytep data;

	typename Matrix::ConstRowIterator i;

	size_t idx;

	// Set up data-structures
	png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (png_ptr == NULL)
		throw LELAError ("Cannot initialise PNG-library");

	info_ptr = png_create_info_struct (png_ptr);

	if (info_ptr == NULL) {
		png_destroy_write_struct (&png_ptr, NULL);
		throw LELAError ("Cannot initialise PNG-library");
	}

	if (setjmp (png_jmpbuf (png_ptr))) {
		png_destroy_write_struct (&png_ptr, &info_ptr);
		throw LELAError ("Error writing PNG-data");
	}

	png_set_write_fn (png_ptr, (png_voidp) &os, PNGWriteData, PNGFlush);

	// Write the header
	png_set_compression_level(png_ptr, COMPRESSION_LEVEL_PNG);
	png_set_IHDR (png_ptr, info_ptr, A.coldim (), A.rowdim (), 1, 0, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info (png_ptr, info_ptr);

	// Allocate data-buffer
	data = new png_byte[(A.coldim () + WordTraits<png_byte>::bits - 1) / WordTraits<png_byte>::bits];

	// Write image-data to file
	for (i = A.rowBegin (), idx = 0; i != A.rowEnd (); ++idx, ++i) {
		copyToPNGData (data, *i, A.coldim ());
		png_write_row (png_ptr, data);
	}

	// Clean up
	png_write_end (png_ptr, NULL);
	png_destroy_write_struct (&png_ptr, &info_ptr);

	delete[] data;

	return os;
}

}

#endif // __LELA_MATRIX_IO_PNG_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
