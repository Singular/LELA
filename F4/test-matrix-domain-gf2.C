/* -*- mode: c++;-*-
 * test-matrix-domain-gf2.C
 * Tests for specialisation of MatrixDomain to GF2
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2011 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#define __LINBOX_BITVECTOR_WORD_TYPE uint8

#include "linbox/field/gf2.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/dense-zero-one.h"
#include "linbox/matrix/matrix-domain-gf2.h"
#include "linbox/matrix/sparse.h"
#include "linbox/matrix/sparse-zero-one.h"
#include "linbox/matrix/submatrix.h"

namespace F4Tests {

using namespace LinBox;

typedef GF2 Field;
typedef __LINBOX_BITVECTOR_WORD_TYPE Word;
typedef BigEndian<Word> Endianness;
typedef std::vector<size_t> SparseVector;
typedef std::pair<std::vector<uint16>, BitVector<Endianness> > HybridVector;
typedef LinBox::SparseMatrix<bool, SparseVector, VectorCategories::SparseZeroOneVectorTag> SparseMatrix;
typedef LinBox::SparseMatrix<bool, HybridVector, VectorCategories::HybridZeroOneVectorTag> HybridMatrix;
typedef LinBox::DenseZeroOneMatrix<BitVector<Endianness>::word_iterator, BitVector<Endianness>::const_word_iterator, Endianness> DenseMatrix;

void testGEMMSubmatrixDense ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	Field F;
	MatrixDomain<Field> MD (F);

	size_t n = 96;
	size_t m = 96;
	size_t l = 96;
	size_t k = 20;

	RandomDenseStream<Field, DenseMatrix::Row> A_stream (F, m, n);
	RandomDenseStream<Field, DenseMatrix::Row> B_stream (F, l, m);

	DenseMatrix A (A_stream), B (B_stream), C (n, l), Cp (n, k);

	DenseMatrix::RowIterator iter;

	for (iter = A.rowBegin (); iter != A.rowEnd (); ++iter)
		A_stream >> *iter;
	
	for (iter = B.rowBegin (); iter != B.rowEnd (); ++iter)
		B_stream >> *iter;

	MD.gemm (true, A, B, false, C);

	std::cout << __FUNCTION__ << ": Input matrix A: " << std::endl;
	MD.write (std::cout, A);

	std::cout << __FUNCTION__ << ": Input matrix B: " << std::endl;
	MD.write (std::cout, B);

	std::cout << __FUNCTION__ << ": Output matrix C:=AB: " << std::endl;
	MD.write (std::cout, C);

	for (int s = 0; s <= l - k; ++s) {
		std::cout << __FUNCTION__ << ": Testing offset " << s << "...";
		std::cout.flush ();

		Submatrix<DenseMatrix > Bp (B, 0, s, B.rowdim (), k);
		MD.gemm (true, A, Bp, false, Cp);

		Submatrix<DenseMatrix > Cpp (C, 0, s, C.rowdim (), k);

		if (MD.areEqual (Cp, Cpp))
			std::cout << "okay" << std::endl;
		else {
			std::cout << "NOT okay" << std::endl;

			std::cout << __FUNCTION__ << ": Computed submatrix with column-offset " << s << ": " << std::endl;
			MD.write (std::cout, Cp);

			std::cout << __FUNCTION__ << ": terminated (ERROR), offset " << s << std::endl;
			return;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;
}

void testGEMMSubmatrixHybrid ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	Field F;
	MatrixDomain<Field> MD (F);

	size_t n = 96;
	size_t m = 96;
	size_t l = 96;
	size_t k = 20;
	double p = 0.1;

#if 1
	RandomSparseStream<Field, SparseMatrix::Row, Field::RandIter, VectorCategories::SparseZeroOneVectorTag> A_stream (F, p, m, n);
	RandomSparseStream<Field, HybridMatrix::Row, Field::RandIter, VectorCategories::HybridZeroOneVectorTag> B_stream (F, p, l, m);

	SparseMatrix A (A_stream);
	HybridMatrix B (B_stream);
	DenseMatrix C (n, l), Cp (n, k);

	std::ofstream A_output ("A.out");
	MD.write (A_output, A, FORMAT_GUILLAUME);
	A_output.close ();

	std::ofstream B_output ("B.out");
	MD.write (B_output, B, FORMAT_GUILLAUME);
	B_output.close ();
#else // 0
	SparseMatrix A (n, m);
	HybridMatrix B (m, l);
	DenseMatrix C (n, l), Cp (n, k);

	std::ifstream A_input ("A.out");
	MD.read (A_input, A, FORMAT_GUILLAUME);
	A_input.close ();

	std::ifstream B_input ("B.out");
	MD.read (B_input, B, FORMAT_GUILLAUME);
	B_input.close ();
#endif // 0

	MD.gemm (true, A, B, false, C);

	std::cout << __FUNCTION__ << ": Input matrix A: " << std::endl;
	MD.write (std::cout, A);

	std::cout << __FUNCTION__ << ": Input matrix B: " << std::endl;
	MD.write (std::cout, B);

	std::cout << __FUNCTION__ << ": Output matrix C:=AB: " << std::endl;
	MD.write (std::cout, C);

	for (int s = 0; s <= l - k; ++s) {
		std::cout << __FUNCTION__ << ": Testing offset " << s << "...";
		std::cout.flush ();

		Submatrix<HybridMatrix> Bp (B, 0, s, B.rowdim (), k);
		MD.gemm (true, A, Bp, false, Cp);

		Submatrix<DenseMatrix> Cpp (C, 0, s, C.rowdim (), k);

		if (MD.areEqual (Cp, Cpp))
			std::cout << "okay" << std::endl;
		else {
			std::cout << "NOT okay" << std::endl;

			std::cout << __FUNCTION__ << ": Computed submatrix with column-offset " << s << ": " << std::endl;
			MD.write (std::cout, Cp);

			std::cout << __FUNCTION__ << ": terminated (ERROR), offset " << s << std::endl;
			return;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;
}

} // namespace F4Tests

int main (int argc, char **argv)
{
	F4Tests::testGEMMSubmatrixDense ();
	F4Tests::testGEMMSubmatrixHybrid ();
}
