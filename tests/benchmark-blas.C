/* tests/benchmark-blas.C
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Benchmarks for BLAS-routines
 */

#include "linbox/util/commentator.h"
#include "linbox/blas/context.h"
#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"
#include "linbox/vector/stream.h"

#include "test-common.h"
#include "test-blas-level2.h"
#include "test-blas-level3.h"

using namespace LinBox;

static long l = 5000;
static long n = 5000;
static long m = 5000;
static long k = 50;
static integer q_uint8 = 101U;
static integer q_uint32 = 2147483647U;
static int iterations = 1;
static bool enable_dense = false;
static bool enable_sparse = false;
static bool enable_gf2 = false;
static bool enable_uint8 = false;
static bool enable_uint32 = false;

template <class Field>
void runBenchmarks (const Field &F)
{
	typedef typename Field::Element Element;

	Context<Field> ctx (F);

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));
	typename Field::Element a;
	r.random (a);

	DenseMatrix<Element> C (l, n), D (l, m), D1 (l, m);

	if (enable_dense) {
		RandomDenseStream<Field, typename DenseMatrix<Element>::Row> stream1 (F, m, l);
		RandomDenseStream<Field, typename DenseMatrix<Element>::Row> stream2 (F, n, m);

		DenseMatrix<Element> M1 (stream1);
		DenseMatrix<Element> M2 (stream2);

		commentator.start ("gemm (dense)", "gemm");
		BLAS3::gemm (ctx, a, M1, M2, F.zero (), C);
		commentator.stop ("done");
	}

	if (enable_sparse) {
		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		RandomSparseStream<Field, typename SparseMatrix<Element>::Row> stream1 (F, (double) k / (double) m, m, l);
		RandomSparseStream<Field, typename SparseMatrix<Element>::Row> stream2 (F, (double) k / (double) n, n, m);

		SparseMatrix<Element> M1 (stream1);
		SparseMatrix<Element> M2 (stream2);

		TransposeMatrix<SparseMatrix<Element> > M2T (M2);

		SparseMatrix<Element> M2Tp (M2.coldim (), M2.rowdim ());

		M2.transpose (M2Tp);

		commentator.start ("gemm (SparseVector)", "gemm");
		BLAS3::gemm (ctx, a, M1, M2Tp, F.zero (), D1);
		commentator.stop ("done");

		commentator.start ("gemm (SparseVector transpose)", "gemm");
		BLAS3::gemm (ctx, a, M1, M2T, F.zero (), D);
		commentator.stop ("done");

		report << "Origial M2^T:" << std::endl;
		BLAS3::write (ctx, report, M2T);

		report << "Copy of M2^T:" << std::endl;
		BLAS3::write (ctx, report, M2Tp);

		report << "Original product:" << std::endl;
		BLAS3::write (ctx, report, D1);

		report << "Check product:" << std::endl;
		BLAS3::write (ctx, report, D);

		if (!BLAS3::equal (ctx, D, D1))
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Products are not equal!" << std::endl;
	}
}

#if 0

template <>
void runBenchmarks (const GF2 &F)
{
	typedef GF2 Field;
	typedef Field::Element Element;

	if (enable_dense) {
		RandomDenseStream<Field, Dense01Matrix<>::Row> stream1 (F, m, l);
		RandomDenseStream<Field, Dense01Matrix<>::Row> stream2 (F, n, m);

		Dense01Matrix<> M1 (stream1);
		Dense01Matrix<> M2 (stream2);

		testGemmCoeff (F, "dense", M1, M2);
	}

	if (enable_sparse) {
		typedef std::vector<uint32> Sparse01Vector;
		typedef LinBox::SparseMatrix<bool, Sparse01Vector, VectorRepresentationTypes::Sparse01> Sparse01Matrix;

		RandomSparseStream<Field, Sparse01Matrix::Row> stream1 (F, (double) k / (double) m, m, l);
		RandomSparseStream<Field, Sparse01Matrix::Row> stream2 (F, (double) k / (double) n, n, m);

		Sparse01Matrix M1 (stream1);
		Sparse01Matrix M2 (stream2);

		testGemmCoeff (F, "sparse (SparseVector)", M1, M2);
	}
}

#endif // Disabled

int main (int argc, char **argv)
{
	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set column-dimension of B to N.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint8 modulus.", TYPE_INTEGER, &q_uint8 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q_uint32 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ 'd', "-d", "Enable dense tests", TYPE_NONE, &enable_dense },
		{ 's', "-s", "Enable sparse tests (SparseVector)", TYPE_NONE, &enable_sparse },
		{ '2', "-2", "Enable tests for GF(2)", TYPE_NONE, &enable_gf2 },
		{ 'b', "-b", "Enable tests for integers mod uint8", TYPE_NONE, &enable_uint8 },
		{ 'w', "-w", "Enable tests for integers mod uint32", TYPE_NONE, &enable_uint32 },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	if (!(enable_dense || enable_sparse)) {
		printHelpMessage (argv[0], args);
		return 0;
	}

	if (!(enable_gf2 || enable_uint8 || enable_uint32)) {
		printHelpMessage (argv[0], args);
		return 0;
	}

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, true, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_IMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator.start ("Matrix domain benchmark suite", "MatrixDomain");

#if 0
	if (enable_gf2) {
		GF2 F;
		runBenchmarks (F);
	}
#endif // Disabled

	if (enable_uint8) {
		Modular<uint8> F (q_uint8);
		runBenchmarks (F);
	}

	if (enable_uint32) {
		Modular<uint32> F (q_uint32);
		runBenchmarks (F);
	}

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
