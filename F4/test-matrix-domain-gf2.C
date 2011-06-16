/* -*- mode: c++;-*-
 * test-matrix-domain-gf2.C
 * Tests for specialisation of MatrixDomain to GF2
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2011 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include "linbox/ring/gf2.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/dense-zero-one.h"
#include "linbox/matrix/matrix-domain-gf2.h"
#include "linbox/matrix/sparse.h"
#include "linbox/matrix/sparse-zero-one.h"
#include "linbox/matrix/submatrix.h"

namespace F4Tests {

using namespace LinBox;

typedef Vector<GF2>::Dense::Endianness Endianness;
typedef DenseZeroOneMatrix<> DenseMatrix;

void testGEMMSubmatrixDense ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	GF2 F;
	MatrixDomain<GF2> MD (F);

	size_t n = 96;
	size_t m = 96;
	size_t l = 96;
	size_t k = 20;

	RandomDenseStream<GF2, DenseMatrix::Row> A_stream (F, m, n);
	RandomDenseStream<GF2, DenseMatrix::Row> B_stream (F, l, m);

	DenseMatrix A (A_stream), B (B_stream), C (n, l), Cp (n, k);

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

		Submatrix<DenseMatrix> Bp (B, 0, s, B.rowdim (), k);
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

void testGEMMSubmatrixHybrid ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	GF2 F;
	MatrixDomain<GF2> MD (F);

	size_t n = 96;
	size_t m = 96;
	size_t l = 96;
	size_t k = 20;
	double p = 0.1;

#if 1
	RandomSparseStream<GF2, Vector<GF2>::Sparse, GF2::RandIter, VectorCategories::SparseZeroOneVectorTag> A_stream (F, p, m, n);
	RandomSparseStream<GF2, Vector<GF2>::Hybrid, GF2::RandIter, VectorCategories::HybridZeroOneVectorTag> B_stream (F, p, l, m);

	SparseMatrix<GF2::Element, Vector<GF2>::Sparse> A (A_stream);
	SparseMatrix<GF2::Element, Vector<GF2>::Hybrid> B (B_stream);
	DenseMatrix C (n, l), Cp (n, k);

	std::ofstream A_output ("A.out");
	MD.write (A_output, A, FORMAT_GUILLAUME);
	A_output.close ();

	std::ofstream B_output ("B.out");
	MD.write (B_output, B, FORMAT_GUILLAUME);
	B_output.close ();
#else // 0
	SparseMatrix<GF2::Element, Vector<GF2>::Sparse> A (n, m);
	SparseMatrix<GF2::Element, Vector<GF2>::Hybrid> B (m, l);
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

		Submatrix<SparseMatrix<GF2::Element, Vector<GF2>::Hybrid> > Bp (B, 0, s, B.rowdim (), k);
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
