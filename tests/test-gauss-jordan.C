/* tests/test-gauss-jordan.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test for hybrid sparse-dense vector-format
 */

#include <iostream>

#include "test-common.h"

#include <linbox/blas/context.h>
#include <linbox/ring/gf2.h>
#include <linbox/ring/modular.h>
#include <linbox/matrix/dense.h>
#include <linbox/vector/stream.h>
#include <linbox/algorithms/gauss-jordan.h>

using namespace LinBox;

template <class Ring>
bool testDenseGJ (const Ring &F, size_t m, size_t n)
{
	commentator.start ("Testing GaussJordan::DenseRowEchelonForm", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	DenseMatrix<typename Ring::Element> R (m, n);
	DenseMatrix<typename Ring::Element> U (m, m);
	DenseMatrix<typename Ring::Element> UPA (m, n);

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> A_stream (F, n, m);

	DenseMatrix<typename Ring::Element> A (A_stream);

	typename GaussJordan<Ring>::Permutation P;

	Context<Ring> ctx (F);

	GaussJordan<Ring> GJ (ctx);
	size_t rank;
	typename Ring::Element det;

	BLAS3::copy (ctx, A, R);

	GJ.DenseRowEchelonForm (R, U, P, rank, det);

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "U = " << std::endl;
	BLAS3::write (ctx, report, U);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	report << "R = " << std::endl;
	BLAS3::write (ctx, report, R);

	BLAS3::permute_rows (ctx, P.begin (), P.end (), A);

	report << "PA = " << std::endl;
	BLAS3::write (ctx, report, A);

	BLAS3::gemm (ctx, F.one (), U, A, F.zero (), UPA);

	report << "UPA = " << std::endl;
	BLAS3::write (ctx, report, UPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!BLAS3::equal (ctx, UPA, R)) {
		error << "ERROR: UPA != R" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring>
bool testStandardGJ (const Ring &F, size_t m, size_t n, size_t k)
{
	commentator.start ("Testing GaussJordan::StandardRowEchelonForm", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Context<Ring> ctx (F);

	GaussJordan<Ring> GJ (ctx);

	RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> A_stream (F, (double) k / (double) n, n, m);

	SparseMatrix<typename Ring::Element> A (A_stream), Aorig (m, n);
	DenseMatrix<typename Ring::Element> U (m, m);
	DenseMatrix<typename Ring::Element> UPA (m, n);

	BLAS3::copy (ctx, A, Aorig);

	typename GaussJordan<Ring>::Permutation P;

	size_t rank;
	typename Ring::Element det;

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	GJ.StandardRowEchelonForm (A, U, P, rank, det, false, true);

	report << "R = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	report << "U = " << std::endl;
	BLAS3::write (ctx, report, U, FORMAT_PRETTY);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	BLAS3::permute_rows (ctx, P.begin (), P.end (), Aorig);
	report << "PA = " << std::endl;
	BLAS3::write (ctx, report, Aorig, FORMAT_PRETTY);

	BLAS3::gemm (ctx, F.one (), U, Aorig, F.zero (), UPA);
	report << "UPA = " << std::endl;
	BLAS3::write (ctx, report, UPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!BLAS3::equal (ctx, UPA, A)) {
		error << "UPA != R, not okay" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass1 = true, pass2 = true;

	static long m = 96;
	static long n = 96;
	static long k = 10;
	static integer q = 101U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'm', "-m M", "Set row-dimension of matrix A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set column-dimension of matrix A to N.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ 'q', "-q Q", "Operate over the ring ZZ/Q [1] for uint32 modulus.", TYPE_INTEGER, &q },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Modular<uint32> GFq (q);
	GF2 gf2;

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("GaussJordan test suite", "GaussJordan");

	std::ostringstream str;
	str << "Running tests over GF(" << q << ")" << std::ends;

	commentator.start (str.str ().c_str (), "GaussJordan");

	pass1 = testDenseGJ (GFq, m, n) && pass1;
	pass1 = testStandardGJ (GFq, m, n, k) && pass1;

	commentator.stop (MSG_STATUS (pass1));

	commentator.start ("Running tests over GF(2)", "GaussJordan");

	pass2 = testDenseGJ (gf2, m, n) && pass2;
	pass2 = testStandardGJ (gf2, m, n, k) && pass2;

	commentator.stop (MSG_STATUS (pass2));

	commentator.stop (MSG_STATUS (pass1 && pass2));

	return (pass1 && pass2) ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
