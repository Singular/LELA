/* tests/test-blas-zp-module.C
 * Copyright 2001, 2002, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test suite for BLAS-routines using ZpModule
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/util/commentator.h"
#include "lela/blas/context.h"
#include "lela/ring/modular.h"
#include "lela/matrix/dense.h"
#include "lela/matrix/sparse.h"
#include "lela/vector/stream.h"
#include "lela/matrix/transpose.h"

#include "test-common.h"
#include "test-blas-level1.h"
#include "test-blas-level2.h"
#include "test-blas-level3.h"

using namespace LELA;

template <class Element>
bool runTests (const integer &q, const char *text, long l, long m, long n, long p, long k, int iterations)
{
	bool pass = true;

	Modular<Element> F (q);
	Context<Modular<Element>, ZpModule<Element> > ctx (F);

	Context<Modular<Element>, GenericModule<Modular<Element> > > ctx_gen (F);

	ostringstream str;
	str << "Testing BLAS ZpModule with ring-type " << text << std::ends;

	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Working over ";
	F.write (report) << std::endl;

	if (!testBLAS1 (ctx, text, l, iterations)) pass = false;
	if (!testBLAS1RepsConsistency (ctx, text, l, iterations)) pass = false;

	RandomDenseStream<Modular<Element>, typename Vector<Modular<Element> >::Dense> stream_v1 (F, l, 1);
	RandomDenseStream<Modular<Element>, typename Vector<Modular<Element> >::Dense> stream_v2 (F, m, 1);
	RandomDenseStream<Modular<Element>, typename Vector<Modular<Element> >::Dense> stream_v3 (F, n, 1);
	RandomDenseStream<Modular<Element>, typename Vector<Modular<Element> >::Dense> stream_v4 (F, p, 1);

	typename Vector<Modular<Element> >::Dense v1 (l), v2 (m), v3 (n), v4 (p);
	stream_v1 >> v1;
	stream_v2 >> v2;
	stream_v3 >> v3;
	stream_v4 >> v4;

	RandomDenseStream<Modular<Element>, typename DenseMatrix<Element>::Row> stream11 (F, m, l);
	RandomDenseStream<Modular<Element>, typename DenseMatrix<Element>::Row> stream12 (F, n, m);
	RandomDenseStream<Modular<Element>, typename DenseMatrix<Element>::Row> stream13 (F, p, n);
	RandomDenseStream<Modular<Element>, typename DenseMatrix<Element>::Row> stream14 (F, m, m);

	DenseMatrix<Element> M1 (stream11);
	DenseMatrix<Element> M2 (stream12);
	DenseMatrix<Element> M3 (stream13);
	DenseMatrix<Element> M4 (stream14);

	if (!testBLAS2 (ctx, "dense", M1, M2, v1, v2,
			typename DenseMatrix<Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "dense", M1, M2, M3, M4,
			typename DenseMatrix<Element>::IteratorType ()))
		pass = false;

	RandomSparseStream<Modular<Element>, typename SparseMatrix<Element>::Row> stream21 (F, (double) k / (double) m, m, l);
	RandomSparseStream<Modular<Element>, typename SparseMatrix<Element>::Row> stream22 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Modular<Element>, typename SparseMatrix<Element>::Row> stream23 (F, (double) k / (double) p, p, n);
	RandomSparseStream<Modular<Element>, typename SparseMatrix<Element>::Row> stream24 (F, (double) k / (double) p, m, m);

	SparseMatrix<Element> M5 (stream21);
	SparseMatrix<Element> M6 (stream22);
	SparseMatrix<Element> M7 (stream23);
	SparseMatrix<Element> M8 (stream24);

	if (!testBLAS2 (ctx, "sparse row-wise", M5, M6, v1, v2,
			typename SparseMatrix<Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "sparse row-wise", M5, M6, M7, M8,
			typename SparseMatrix<Element>::IteratorType ()))
		pass = false;

	TransposeMatrix<SparseMatrix<Element> > M9 (M7);
	TransposeMatrix<SparseMatrix<Element> > M10 (M6);
	TransposeMatrix<SparseMatrix<Element> > M11 (M5);

	RandomSparseStream<Modular<Element>, typename SparseMatrix<Element>::Row> stream31 (F, (double) k / (double) n, n, n);

	SparseMatrix<Element> M12 (stream31);
	TransposeMatrix<SparseMatrix<Element> > M12T (M12);

	if (!testBLAS2 (ctx, "sparse column-wise", M9, M10, v4, v3,
			typename TransposeMatrix<SparseMatrix<Element> >::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "sparse column-wise", M9, M10, M11, M12T,
			typename TransposeMatrix<SparseMatrix<Element> >::IteratorType ()))
		pass = false;

	pass = testBLAS2ModulesConsistency(ctx, ctx_gen, text, m, n, k) && pass;
	pass = testBLAS2RepsConsistency(ctx, text, m, n, k) && pass;
	pass = testBLAS3ModulesConsistency(ctx, ctx_gen, text, m, n, p, k) && pass;
	pass = testBLAS3RepsConsistency(ctx, text, m, n, p, k) && pass;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long l = 50;
	static long m = 50;
	static long n = 30;
	static long p = 30;
	static long k = 10;
	static integer q_integer = 65521;
	static integer q_uint32 = 2147483647;
	static integer q_uint16 = 65521;
	static integer q_uint8 = 251;
	static integer q_float_small = 2039;
	static integer q_float_big = 4093;
	static integer q_double_small = 33554393;
	static integer q_double_big = 67108859;
	static int iterations = 1;

	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set row-dimension of matrix C and column-dimension of B to N.", TYPE_INT, &n },
		{ 'p', "-p P", "Set column-dimension of matrix C to P.", TYPE_INT, &p },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint8 modulus", TYPE_INTEGER, &q_uint8 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (7);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (7);

	commentator.start ("BLAS ZpModule test-suite", "ZpModule");

	pass = runTests<integer> (q_integer, "Modular<integer>", l, m, n, p, k, iterations) && pass;
	pass = runTests<uint32> (q_uint32, "Modular<uint32>", l, m, n, p, k, iterations) && pass;
	pass = runTests<uint16> (q_uint16, "Modular<uint16>", l, m, n, p, k, iterations) && pass;
	pass = runTests<uint8> (q_uint8, "Modular<uint8>", l, m, n, p, k, iterations) && pass;
	pass = runTests<float> (q_float_small, "Modular<float>", l, m, n, p, k, iterations) && pass;
	pass = runTests<float> (q_float_big, "Modular<float>", l, m, n, p, k, iterations) && pass;
	pass = runTests<double> (q_double_small, "Modular<double>", l, m, n, p, k, iterations) && pass;
	pass = runTests<double> (q_double_big, "Modular<double>", l, m, n, p, k, iterations) && pass;

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
