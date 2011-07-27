/* tests/benchmark-blas.C
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Benchmarks for BLAS-routines
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/util/commentator.h"
#include "lela/blas/context.h"
#include "lela/ring/gf2.h"
#include "lela/ring/modular.h"
#include "lela/matrix/dense.h"
#include "lela/matrix/sparse.h"
#include "lela/vector/stream.h"

#include "test-common.h"
#include "test-blas-level2.h"
#include "test-blas-level3.h"

using namespace LELA;

static long l = 5000;
static long n = 5000;
static long m = 5000;
static long k = 100;
static integer q_uint8 = 101U;
static integer q_uint32 = 2147483647U;
static integer q_float = 101U;
static integer q_double = 2000003U;
static int iterations = 1;
static bool enable_dense = true;
static bool enable_sparse = false;
static bool enable_hybrid = true;
static bool enable_gf2 = true;
static bool enable_uint8 = true;
static bool enable_uint32 = true;
static bool enable_float = true;
static bool enable_double = true;

template <class Ring, class Modules>
void runBenchmarks (Context<Ring, Modules> &ctx, const char *text1, const char *text2)
{
	typedef typename Ring::Element Element;

	std::ostringstream str;
	str << "Running benchmarks for " << text1 << " over " << text2 << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	NonzeroRandIter<Ring> r (ctx.F, typename Ring::RandIter (ctx.F));
	typename Ring::Element a, b;
	r.random (a);
	r.random (b);

	if (enable_dense) {
		RandomDenseStream<Ring, typename DenseMatrix<Element>::Row> stream1 (ctx.F, m, l);
		RandomDenseStream<Ring, typename DenseMatrix<Element>::Row> stream2 (ctx.F, n, m);
		RandomDenseStream<Ring, typename DenseMatrix<Element>::Row> stream3 (ctx.F, n, l);

		DenseMatrix<Element> M1 (stream1);
		DenseMatrix<Element> M2 (stream2);
		DenseMatrix<Element> M3 (stream3);

		commentator.start ("gemm (dense)", "gemm");
		BLAS3::gemm (ctx, a, M1, M2, b, M3);
		commentator.stop ("done");
	}

	if (enable_sparse) {
		RandomSparseStream<Ring, typename SparseMatrix<Element>::Row> stream1 (ctx.F, (double) k / (double) m, m, l);
		RandomSparseStream<Ring, typename SparseMatrix<Element>::Row> stream2 (ctx.F, (double) k / (double) n, n, m);
		RandomSparseStream<Ring, typename SparseMatrix<Element>::Row> stream3 (ctx.F, (double) k / (double) n, n, l);

		SparseMatrix<Element> M1 (stream1);
		SparseMatrix<Element> M2 (stream2);
		SparseMatrix<Element> M3 (stream3);

		commentator.start ("gemm (SparseVector)", "gemm");
		BLAS3::gemm (ctx, a, M1, M2, b, M3);
		commentator.stop ("done");
	}

	commentator.stop (MSG_DONE);
}

template <class Modules>
void runBenchmarks (Context<GF2, Modules> &ctx, const char *text)
{
	typedef GF2 Ring;
	typedef Ring::Element Element;

	std::ostringstream str;
	str << "Running benchmarks for " << text << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	if (enable_dense) {
		RandomDenseStream<Ring, DenseMatrix<bool>::Row> stream1 (ctx.F, m, l);
		RandomDenseStream<Ring, DenseMatrix<bool>::Row> stream2 (ctx.F, n, m);
		RandomDenseStream<Ring, DenseMatrix<bool>::Row> stream3 (ctx.F, n, l);

		DenseMatrix<bool> M1 (stream1);
		DenseMatrix<bool> M2 (stream2);
		DenseMatrix<bool> M3 (stream3);

		commentator.start ("gemm (dense)", "gemm");
		BLAS3::gemm (ctx, true, M1, M2, true, M3);
		commentator.stop ("done");
	}

	if (enable_sparse) {
		typedef LELA::SparseMatrix<bool> SparseMatrix;

		RandomSparseStream<Ring, SparseMatrix::Row> stream1 (ctx.F, (double) k / (double) m, m, l);
		RandomSparseStream<Ring, SparseMatrix::Row> stream2 (ctx.F, (double) k / (double) n, n, m);
		RandomSparseStream<Ring, SparseMatrix::Row> stream3 (ctx.F, (double) k / (double) n, n, k);

		SparseMatrix M1 (stream1);
		SparseMatrix M2 (stream2);
		SparseMatrix M3 (stream3);

		commentator.start ("gemm (dense)", "gemm");
		BLAS3::gemm (ctx, true, M1, M2, true, M3);
		commentator.stop ("done");
	}

	if (enable_hybrid) {
		typedef LELA::SparseMatrix<bool, Vector<GF2>::Hybrid> HybridMatrix;

		RandomHybridStream<Ring, HybridMatrix::Row> stream1 (ctx.F, (double) k / (double) m, m, l);
		RandomHybridStream<Ring, HybridMatrix::Row> stream2 (ctx.F, (double) k / (double) n, n, m);
		RandomHybridStream<Ring, HybridMatrix::Row> stream3 (ctx.F, (double) k / (double) n, n, l);

		HybridMatrix M1 (stream1);
		HybridMatrix M2 (stream2);
		HybridMatrix M3 (stream3);

		commentator.start ("gemm (dense)", "gemm");
		BLAS3::gemm (ctx, true, M1, M2, true, M3);
		commentator.stop ("done");
	}

	commentator.stop (MSG_DONE);
}

template <class Element>
void runBenchmarksForElement (const integer &q, const char *text)
{
	Modular<Element> F (q);

	Context<Modular<Element> > ctx_all (F);
	runBenchmarks (ctx_all, "AllModules", text);

	Context<Modular<Element>, ZpModule<Element> > ctx_zp (F);
	runBenchmarks (ctx_zp, "ZpModule", text);

	Context<Modular<Element>, GenericModule<Modular<Element> > > ctx_gen (F);
	runBenchmarks (ctx_gen, "GenericModule", text);
}

int main (int argc, char **argv)
{
	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set column-dimension of B to N.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'q', "-q Q", "Operate over the ring Z/Q for uint8 modulus.", TYPE_INTEGER, &q_uint8 },
		{ 'Q', "-Q Q", "Operate over the ring Z/Q for uint32 modulus.", TYPE_INTEGER, &q_uint32 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ 'd', "-d", "Enable dense tests", TYPE_NONE, &enable_dense },
		{ 's', "-s", "Enable sparse tests (SparseVector)", TYPE_NONE, &enable_sparse },
		{ '2', "-2", "Enable tests for GF(2)", TYPE_NONE, &enable_gf2 },
		{ 'b', "-b", "Enable tests for integers mod uint8", TYPE_NONE, &enable_uint8 },
		{ 'w', "-w", "Enable tests for integers mod uint32", TYPE_NONE, &enable_uint32 },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	if (!(enable_dense || enable_sparse || enable_hybrid)) {
		printHelpMessage (argv[0], args);
		return 0;
	}

	if (!(enable_gf2 || enable_uint8 || enable_uint32 || enable_float || enable_double)) {
		printHelpMessage (argv[0], args);
		return 0;
	}

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, true, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (6);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (6);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator.start ("BLAS3 benchmark suite", "BLAS3");

	if (enable_gf2) {
		GF2 F;

		Context<GF2> ctx_all (F);
		runBenchmarks (ctx_all, "AllModules", "GF2");

		Context<GF2> ctx_gen (F);
		runBenchmarks (ctx_gen, "GenericModule", "GF2");
	}

	if (enable_float)
		runBenchmarksForElement<float> (q_float, "Modular<float>");

	if (enable_double)
		runBenchmarksForElement<double> (q_double, "Modular<double>");

	if (enable_uint8)
		runBenchmarksForElement<uint8> (q_uint8, "Modular<uint8>");

	if (enable_uint32)
		runBenchmarksForElement<uint32> (q_uint32, "Modular<uint32>");

	commentator.stop ("done");

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
