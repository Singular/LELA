/* -*- mode: c++;-*-
 * test-from-png.C
 * Test code to read matrix from PNG-file
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include <iostream>

#include "linbox/field/gf2.h"

#include "from-png.C"

using namespace LinBox;
using namespace F4;

typedef GF2 Field;

int main (int argc, char **argv)
{
	if (argc < 2) {
		std::cout << "Usage: test-from-png <png-input>" << std::endl;
		return 0;
	}

	Field F;
	SparseMatrixReader reader (F);

	SparseMatrixReader::SparseMatrix *A = reader.ReadFromPNG (argv[1]);

	if (A == NULL) {
		std::cerr << "Could not read file" << std::endl;
		return 1;
	}

	std::ofstream f ("matrix.out");
	A->write (f, F, FORMAT_GUILLAUME);

	delete A;
}
