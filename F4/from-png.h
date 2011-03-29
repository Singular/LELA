/* -*- mode: c++;-*-
 * from-png.h
 * Read matrix from PNG-file
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#ifndef __F4_FROM_PNG_H
#define __F4_FROM_PNG_H

#include <png.h>

#include <vector>
#include <string>

#include <linbox/matrix/sparse.h>
#include <linbox/matrix/dense.h>
#include <linbox/field/modular.h>

#include "gauss-jordan.h"

namespace F4 {
	using namespace LinBox;

	// FIXME Document this!!!
	class ReaderErrorException {
		const std::string error;

	public:
		ReaderErrorException (std::string _error) : error (_error) {}

		friend std::ostream &operator << (std::ostream &os, const ReaderErrorException &ex);
	};

	// FIXME Document this!!!
	template <class Field>
	class SparseMatrixReader {
	public:
		typedef typename GaussJordan<Field>::SparseMatrix SparseMatrix;
		typedef typename GaussJordan<Field>::DenseMatrix DenseMatrix;

	private:
		const Field &F;
		typename Field::Element one;

		size_t total_nonzero;

		void readBlock (typename SparseMatrix::Row &v, png_byte x, int start, int stop)
			{ readBlockSpecialised (v, x, start, stop, typename VectorTraits<Field, typename SparseMatrix::Row>::VectorCategory ()); }
		void readRow (typename SparseMatrix::Row &v, png_bytep row, int width);

		void readBlockSpecialised (typename SparseMatrix::Row &v, png_byte x, int start, int stop,
					   VectorCategories::SparseVectorTag);
		void readBlockSpecialised (typename SparseMatrix::Row &v, png_byte x, int start, int stop,
					   VectorCategories::SparseZeroOneVectorTag);
		void readBlockSpecialised (typename SparseMatrix::Row &v, png_byte x, int start, int stop,
					   VectorCategories::HybridZeroOneVectorTag);

		template <class Endianness>
		void readBlockHybridSpecialised (typename SparseMatrixReader<Field>::SparseMatrix::Row &v, png_byte x, int start, int stop,
						 Endianness);

	public:
		SparseMatrixReader (const Field &_F) : F (_F) {
			F.init (one, 1);
		}

		SparseMatrix *ReadFromPNG (char *filename);
	};
}

#include "from-png.inl"

#endif // __F4_FROM_PNG_H
