/* tests/test-matrix-domain-gf2.C
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 * Based on test-matrix-domain.C by same author
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test suite for MatrixDomain
 */

#include "linbox/util/commentator.h"
#include "linbox/blas/context.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"
#include "linbox/vector/stream.h"
#include "linbox/ring/gf2.h"

#include "test-common.h"
#include "test-blas-level1.h"
#include "test-blas-level2.h"
#include "test-blas-level3.h"

using namespace LinBox;

typedef GF2 Field;

template <class Modules>
bool testBLAS1 (Context<GF2, Modules> &ctx, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing VectorDomain <" << text << ">" << std::ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<GF2, Vector<GF2>::Dense> stream1 (ctx.F, n, iterations), stream2 (ctx.F, n, iterations);
	RandomSparseStream<GF2, Vector<GF2>::Sparse> stream3 (ctx.F, 0.1, n, iterations), stream4 (ctx.F, 0.1, n, iterations);
	RandomSparseStream<GF2, Vector<GF2>::Hybrid> stream5 (ctx.F, 0.1, n, iterations), stream6 (ctx.F, 0.1, n, iterations);

	if (!testCopyEqual (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "dense/sparse", stream1, stream4)) pass = false;  stream1.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "dense/hybrid", stream1, stream6)) pass = false;  stream1.reset (); stream6.reset ();
	if (!testCopyEqual (ctx, "sparse/dense", stream3, stream2)) pass = false;  stream3.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "sparse/hybrid", stream3, stream6)) pass = false; stream3.reset (); stream6.reset ();
	if (!testCopyEqual (ctx, "hybrid/dense", stream5, stream2)) pass = false;  stream5.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "hybrid/sparse", stream5, stream4)) pass = false; stream5.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "hybrid/hybrid", stream5, stream6)) pass = false; stream5.reset (); stream6.reset ();

	if (!testInequality (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testInequality (ctx, "dense/sparse", stream1, stream4)) pass = false;  stream1.reset (); stream4.reset ();
	if (!testInequality (ctx, "dense/hybrid", stream1, stream6)) pass = false;  stream1.reset (); stream6.reset ();
	if (!testInequality (ctx, "sparse/dense", stream3, stream2)) pass = false;  stream3.reset (); stream2.reset ();
	if (!testInequality (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();
	if (!testInequality (ctx, "sparse/hybrid", stream3, stream6)) pass = false; stream3.reset (); stream6.reset ();
	if (!testInequality (ctx, "hybrid/dense", stream5, stream2)) pass = false;  stream5.reset (); stream2.reset ();
	if (!testInequality (ctx, "hybrid/sparse", stream5, stream4)) pass = false; stream5.reset (); stream4.reset ();
	if (!testInequality (ctx, "hybrid/hybrid", stream5, stream6)) pass = false; stream5.reset (); stream6.reset ();

	if (!testNonzero (ctx, "dense", stream1)) pass = false;   stream1.reset ();
	if (!testNonzero (ctx, "sparse", stream3)) pass = false;  stream3.reset ();
	if (!testNonzero (ctx, "hybrid", stream5)) pass = false;  stream5.reset ();

	if (!testDotProduct (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testDotProduct (ctx, "dense/sparse", stream1, stream3)) pass = false;  stream1.reset (); stream3.reset ();
	if (!testDotProduct (ctx, "dense/hybrid", stream1, stream5)) pass = false;  stream1.reset (); stream5.reset ();
	if (!testDotProduct (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();

	if (!testScal (ctx, "dense", stream1)) pass = false;   stream1.reset ();
	if (!testScal (ctx, "sparse", stream3)) pass = false;  stream3.reset ();
	if (!testScal (ctx, "hybrid", stream5)) pass = false;  stream5.reset ();

	if (!testAXPY (ctx, "dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testAXPY (ctx, "sparse", stream3, stream4)) pass = false;  stream3.reset (); stream4.reset ();
	if (!testAXPY (ctx, "hybrid", stream5, stream6)) pass = false;  stream5.reset (); stream6.reset ();

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long l = 50;
	static long m = 50;
	static long n = 80;
	static long p = 60;
	static long k = 10;
	static int iterations = 1;

	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set row-dimension of matrix C and column-dimension of B to N.", TYPE_INT, &n },
		{ 'p', "-p P", "Set column-dimension of matrix C to P.", TYPE_INT, &p },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Field F;
	Context<GF2, GenericModule<GF2> > ctx (F);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("BLAS GenericModule, GF2-specialisation, test suite", "BLASGenericModuleGF2");

	if (!testBLAS1 (ctx, "GF2", l, iterations)) pass = false;

	RandomDenseStream<Field, Vector<Field>::Dense> stream_v1 (F, l, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v2 (F, m, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v3 (F, n, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v4 (F, p, 1);

	Vector<Field>::Dense v1 (l), v2 (m), v3 (n), v4 (p);
	stream_v1 >> v1;
	stream_v2 >> v2;
	stream_v3 >> v3;
	stream_v4 >> v4;

	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream11 (F, m, l);
	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream12 (F, n, m);
	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream13 (F, p, n);

	DenseMatrix<Field::Element> M1 (stream11);
	DenseMatrix<Field::Element> M2 (stream12);
	DenseMatrix<Field::Element> M3 (stream13);

	if (!testBLAS2 (ctx, "dense", M1, M2, v1, v2, iterations,
			DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "dense", M1, M2, M3, iterations,
			DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;

	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream21 (F, (double) k / (double) m, m, l);
	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream22 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream23 (F, (double) k / (double) p, p, n);

	SparseMatrix<Field::Element, Vector<Field>::Sparse> M4 (stream21);
	SparseMatrix<Field::Element, Vector<Field>::Sparse> M5 (stream22);
	SparseMatrix<Field::Element, Vector<Field>::Sparse> M6 (stream23);

	if (!testBLAS2 (ctx, "sparse", M4, M5, v1, v2, iterations,
			SparseMatrix<Field::Element, Vector<Field>::Sparse>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "sparse", M4, M5, M6, iterations,
			SparseMatrix<Field::Element, Vector<Field>::Sparse>::IteratorType ()))
		pass = false;

	RandomSparseStream<Field, Vector<Field>::Hybrid> stream31 (F, (double) k / (double) m, m, l);
	RandomSparseStream<Field, Vector<Field>::Hybrid> stream32 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Field, Vector<Field>::Hybrid> stream33 (F, (double) k / (double) p, p, n);

	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M7 (stream31);
	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M8 (stream32);
	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M9 (stream33);

	if (!testBLAS2 (ctx, "hybrid", M7, M8, v1, v2, iterations,
			SparseMatrix<Field::Element, Vector<Field>::Hybrid>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "hybrid", M7, M8, M9, iterations,
			SparseMatrix<Field::Element, Vector<Field>::Hybrid>::IteratorType ()))
		pass = false;

	RandomDenseStream<Field, Vector<Field>::Dense> stream_v1p (F, l - 8, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v2p (F, m - 8, 1);

	Vector<Field>::Dense v1p (l - 8), v2p (m - 8);
	stream_v1p >> v1p;
	stream_v2p >> v2p;

	DenseMatrix<Field::Element>::SubmatrixType M10 (M1, 4, 4, l - 8, m - 8);
	DenseMatrix<Field::Element>::SubmatrixType M11 (M2, 4, 4, m - 8, n - 8);
	DenseMatrix<Field::Element>::SubmatrixType M12 (M3, 4, 4, n - 8, p - 8);

	if (!testBLAS2Submatrix (ctx, "dense (submatrix)", M10, M11, v1p, v2p, iterations,
				 DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3Submatrix (ctx, "dense (submatrix)", M10, M11, M12, iterations,
				 DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;

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
