/* lela/tests/test-matrix.C
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic tests for matrices
 * 
 * ------------------------------------
 *
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/ring/gf2.h"
#include "lela/ring/modular.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"
#include "lela/blas/level3.h"
#include "lela/matrix/dense.h"
#include "lela/matrix/sparse.h"
#include "lela/matrix/dense-zero-one.h"
#include "lela/matrix/sparse-zero-one.h"
#include "lela/vector/stream.h"

#include "test-common.h"
#include "test-matrix.h"

template <class Field, class Matrix>
bool runAllTests (const Field &F, const char *text, Matrix &M, size_t n, size_t m)
{
	std::ostringstream str;
	str << "Running matrix-tests (" << text << ")" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool pass = true;

	Context<Field> ctx (F);

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << std::endl;
	BLAS3::write (ctx, report, M);

	pass = testRowdimColdim (F, M, n, m) && pass;
	pass = testRowColIterator (F, M) && pass;
	pass = testGetSetEntry (F, M) && pass;
	pass = testSetEraseEntry (F, M) && pass;
	pass = testRawIterator (F, M) && pass;
	pass = testSubmatrix (F, M) && pass;
	pass = testTransposeMatrix (F, M) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Field, class Matrix>
bool runDenseTests (const Field &F, const char *text, const Matrix &M, size_t n, size_t m)
{
	std::ostringstream str;
	str << "Running matrix-tests (" << text << ", dense tests only)" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool pass = true;

	Context<Field> ctx (F);

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << std::endl;
	BLAS3::write (ctx, report, M);

	pass = testRowdimColdim (F, M, n, m) && pass;
	pass = testRowColIterator (F, M) && pass;
	pass = testGetSetEntry (F, M) && pass;
	pass = testRawIterator (F, M) && pass;
	pass = testSubmatrix (F, M) && pass;
	pass = testTransposeMatrix (F, M) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long m = 10;
	static long n = 10;
	static long k = 3;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'm', "-m M", "Set row-dimension of matrix to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set column-dimension of matrix to N.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint32> Field;
	typedef Field::Element Element;

	Field F (q);

	commentator.start ("Matrix-test-suite", "Matrix");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	RandomDenseStream<Field, DenseMatrix<Element>::Row> stream1 (F, n, m);
	RandomSparseStream<Field, SparseMatrix<Element>::Row> stream2 (F, (double) k / (double) m, n, m);

	DenseMatrix<Element> M1 (stream1);
	SparseMatrix<Element> M2 (stream2);

	pass = runDenseTests (F, "dense GF(q)", M1, m, n) && pass;
	pass = runAllTests (F, "sparse row-wise GF(q)", M2, m, n) && pass;

	GF2 gf2;

	RandomDenseStream<GF2, DenseMatrix<GF2::Element>::Row> stream3 (gf2, n, m);
	RandomSparseStream<GF2, Vector<GF2>::Sparse> stream4 (gf2, (double) k / (double) m, m, n);
	RandomSparseStream<GF2, Vector<GF2>::Hybrid> stream5 (gf2, (double) k / (double) m, m, n);

	DenseMatrix<GF2::Element> M3 (stream3);
	SparseMatrix<GF2::Element, Vector<GF2>::Sparse> M4 (stream4);
	SparseMatrix<GF2::Element, Vector<GF2>::Hybrid> M5 (stream5);

	pass = runDenseTests (gf2, "dense GF(2)", M3, m, n) && pass;
	pass = runAllTests (gf2, "sparse row-wise GF(2)", M4, m, n) && pass;
	pass = runDenseTests (gf2, "hybrid row-wise GF(2)", M5, m, n) && pass;

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
