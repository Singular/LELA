/* -*- mode: c++;-*-
 * from-png.inl
 * Read matrix from PNG-file
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#ifndef __FROM_PNG_INL
#define __FROM_PNG_INL

#include <cstdio>
#include <iostream>
#include <png.h>

namespace F4 {
	std::ostream &operator << (std::ostream &os, const ReaderErrorException &ex) {
		os << ex.error;
	}

	template <class Field>
	void SparseMatrixReader<Field>::readBlockSpecialised (typename SparseMatrixReader<Field>::SparseMatrix::Row &v, png_byte x, int start, int stop,
							      VectorCategories::SparseVectorTag)
	{
		int idx;
		png_byte t;

		for (idx = 0, t = ~(((png_byte) -1) & (((png_byte) -1) >> 1)); idx < stop; ++idx, t >>= 1)
			if (!(x & t))
				v.push_back (std::pair<size_t, typename Field::Element> (start + idx, F.one ()));
	}

	template <class Field>
	void SparseMatrixReader<Field>::readBlockSpecialised (typename SparseMatrixReader<Field>::SparseMatrix::Row &v, png_byte x, int start, int stop,
							      VectorCategories::SparseZeroOneVectorTag)
	{
		int idx;
		png_byte t;

		for (idx = 0, t = ~(((png_byte) -1) & (((png_byte) -1) >> 1)); idx < stop; ++idx, t >>= 1) {
			if (!(x & t)) {
				v.push_back (start + idx);
				++total_nonzero;
			}
		}
	}

	template <class Field>
	void SparseMatrixReader<Field>::readBlockSpecialised (typename SparseMatrixReader<Field>::SparseMatrix::Row &v, png_byte x, int start, int stop,
							      VectorCategories::HybridZeroOneVectorTag)
	{
		readBlockHybridSpecialised (v, x, start, stop, typename SparseMatrix::Row::second_type::Endianness ());
	}

	template <class Field>
	template <class Endianness>
	void SparseMatrixReader<Field>::readBlockHybridSpecialised (typename SparseMatrixReader<Field>::SparseMatrix::Row &v, png_byte x, int start, int stop, Endianness)
	{
		typedef typename std::iterator_traits<typename SparseMatrixReader<Field>::SparseMatrix::Row::first_type::iterator>::value_type index_type;
		typedef typename std::iterator_traits<typename SparseMatrixReader<Field>::SparseMatrix::Row::second_type::word_iterator>::value_type word_type;

		index_type idx = start >> WordTraits<word_type>::logof_size;
		int count;
		png_byte t;

		word_type mask = Endianness::e_j (start & WordTraits<word_type>::pos_mask);

		for (count = 0, t = ~(((png_byte) -1) & (((png_byte) -1) >> 1)); count < stop; t >>= 1, mask = Endianness::shift_right (mask, 1), ++count) {
			if (!(x & t)) {
				if (v.first.empty () || v.first.back () != idx) {
					v.first.push_back (idx);
					v.second.push_word_back (mask);
				}
				else {
					v.second.back_word () |= mask;
				}

				++total_nonzero;
			}
		}
	}

	template <class Field>
	void SparseMatrixReader<Field>::readRow (typename SparseMatrixReader<Field>::SparseMatrix::Row &v, png_bytep row, int width)
	{
		int i, j;

		for (i = 0; i < width / (8 * sizeof (png_byte)); ++i)
			readBlock (v, row[i], i * (8 * sizeof (png_byte)), 8 * sizeof (png_byte));

		if (width % (8 * sizeof (png_byte)))
			readBlock (v, row[i], i * (8 * sizeof (png_byte)), width % (8 * sizeof (png_byte)));
	}

	template <class Field>
	typename SparseMatrixReader<Field>::SparseMatrix *SparseMatrixReader<Field>::ReadFromPNG (char *filename) {
		png_byte header[8];
		const static int number = 8;
		bool is_png;

		png_structp png_ptr;
		png_infop info_ptr, end_ptr;
		png_bytepp row_pointers;

		png_uint_32 width, height;
		int bit_depth, colour_type, filter_method, compression_type, interlace_type;

		FILE *F = fopen (filename, "rb");

		if (!F) {
			std::string err = "Cannot open file ";
			err += filename;

			throw ReaderErrorException (err);
		}

		png_ptr = png_create_read_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

		if (png_ptr == NULL)
			throw ReaderErrorException ("Cannot initialise PNG-library");

		info_ptr = png_create_info_struct (png_ptr);

		if (info_ptr == NULL) {
			png_destroy_read_struct (&png_ptr, NULL, NULL);
			throw ReaderErrorException ("Cannot initialise PNG-library");
		}

		end_ptr = png_create_info_struct (png_ptr);

		if (end_ptr == NULL) {
			png_destroy_read_struct (&png_ptr, &info_ptr, NULL);
			throw ReaderErrorException ("Cannot initialise PNG-library");
		}

		if (setjmp (png_jmpbuf (png_ptr))) {
			png_destroy_read_struct (&png_ptr, &info_ptr, &end_ptr);

			std::string err ("File ");
			err += filename;
			err += " cannot be read by PNG-reader (is it really a PNG-file)?";

			throw ReaderErrorException (err);
		}

		png_init_io (png_ptr, F);
		png_read_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
		fclose (F);
		png_get_IHDR (png_ptr, info_ptr, &width, &height, &bit_depth, &colour_type, &filter_method, &compression_type, &interlace_type);

		if (bit_depth != 1) {
			std::cerr << "Incorrect bit_depth " << bit_depth << " (should be 1)" << std::endl;
			return NULL;
		}

		if (colour_type != 0) {
			std::cerr << "Incorrect colour_type " << colour_type << " (should be 0)" << std::endl;
			return NULL;
		}

		std::cout << "Image is " << height << " by " << width << std::endl;

		row_pointers = png_get_rows (png_ptr, info_ptr);

		SparseMatrix *matrix = new SparseMatrix (height, width);

		typename SparseMatrix::RowIterator i = matrix->rowBegin ();
		int idx;

		std::cout << "Placing PNG-data into matrix [          ]\b\b\b\b\b\b\b\b\b\b\b";
		std::cout.flush ();
		float prog = 0.1;

		total_nonzero = 0;

		for (idx = 0; idx < height; ++idx, ++i) {
			if ((float) idx / (float) height >= prog) {
				std::cout << '.';
				std::cout.flush ();
				prog += .1;
			}
			readRow (*i, row_pointers[idx], width);
		}

		std::cout << '.' << std::endl;

		std::cout << "Total non-zero entries: " << total_nonzero << " / " << width * height << " ("
			  << 100 * ((double) total_nonzero) / ((double) width * height) << "%)" << std::endl;

		png_destroy_read_struct (&png_ptr, &info_ptr, &end_ptr);

		return matrix;
	}
}

#endif // __FROM_PNG_INL
