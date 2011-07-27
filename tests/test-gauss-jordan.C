/* tests/test-gauss-jordan.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test for Gauss-Jordan transform
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include <iostream>

#include "test-common.h"

#include <lela/blas/context.h>
#include <lela/ring/gf2.h>
#include <lela/ring/modular.h>
#include <lela/matrix/dense.h>
#include <lela/vector/stream.h>
#include <lela/algorithms/elimination.h>
#include <lela/algorithms/gauss-jordan.h>

using namespace LELA;

template <class Ring>
bool testGaussTransform (const Ring &F, size_t m, size_t n)
{
	commentator.start ("Testing GaussJordan::echelonize", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	DenseMatrix<typename Ring::Element> R (m, n);

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> A_stream (F, n, m);

	DenseMatrix<typename Ring::Element> A (A_stream), A_copy (m, n);

	typename GaussJordan<Ring>::Permutation P;

	Context<Ring> ctx (F);

	Elimination<Ring> elim (ctx);
	GaussJordan<Ring> GJ (ctx);
	size_t rank;
	typename Ring::Element det;

	BLAS3::copy (ctx, A, R);
	BLAS3::copy (ctx, A, A_copy);

	GJ.echelonize (R, P, rank, det);

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	report << "R, L = " << std::endl;
	BLAS3::write (ctx, report, R);

	BLAS3::permute_rows (ctx, P.begin (), P.end (), A);

	report << "PA = " << std::endl;
	BLAS3::write (ctx, report, A);

	typename DenseMatrix<typename Ring::Element>::SubmatrixType Rp (R, 0, 0, R.rowdim (), R.rowdim ());

	BLAS3::trmm (ctx, F.one (), Rp, A, LowerTriangular, true);

	report << "LPA = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	// Trick to eliminate part below diagonal so that equality-check works
	elim.move_L (R, R);

	if (!BLAS3::equal (ctx, A, R)) {
		error << "ERROR: LPA != R" << std::endl;
		pass = false;
	}

	elim.echelonize (A_copy, P, rank, det, false);

	report << "Result of Elimination::echelonize: " << std::endl;
	BLAS3::write (ctx, report, A_copy);

	if (!BLAS3::equal (ctx, A_copy, R)) {
		error << "ERROR: Results from Elimination and GaussJordan not equal" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring>
bool testGaussJordanTransform (const Ring &F, size_t m, size_t n)
{
	commentator.start ("Testing GaussJordan::echelonize_reduced", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	DenseMatrix<typename Ring::Element> R (m, n);

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> A_stream (F, n, m);

	DenseMatrix<typename Ring::Element> A (A_stream), L (m, m), LPA (m, n);

	typename GaussJordan<Ring>::Permutation P;

	Context<Ring> ctx (F);

	GaussJordan<Ring> GJ (ctx);
	size_t rank;
	typename Ring::Element det;

	BLAS3::copy (ctx, A, R);

	GJ.echelonize_reduced (R, L, P, rank, det);

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	report << "R = " << std::endl;
	BLAS3::write (ctx, report, R);

	report << "L = " << std::endl;
	BLAS3::write (ctx, report, L);

	BLAS3::permute_rows (ctx, P.begin (), P.end (), A);

	report << "PA = " << std::endl;
	BLAS3::write (ctx, report, A);

	BLAS3::scal (ctx, F.zero (), LPA);
	BLAS3::gemm (ctx, F.one (), L, A, F.zero (), LPA);

	report << "LPA = " << std::endl;
	BLAS3::write (ctx, report, LPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	// Trick to eliminate part below diagonal so that equality-check works
	Elimination<Ring> elim (ctx);
	elim.move_L (R, R);

	if (!BLAS3::equal (ctx, LPA, R)) {
		error << "ERROR: LPA != R" << std::endl;
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

	Modular<float> GFq (q);
	GF2 gf2;

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("GaussJordan test suite", "GaussJordan");

	std::ostringstream str;
	str << "Running tests over GF(" << q << ")" << std::ends;

	commentator.start (str.str ().c_str (), "GaussJordan");

	pass1 = testGaussTransform (GFq, m, n) && pass1;
	pass1 = testGaussJordanTransform (GFq, m, n) && pass1;

	commentator.stop (MSG_STATUS (pass1));

	commentator.start ("Running tests over GF(2)", "GaussJordan");

	pass2 = testGaussTransform (gf2, m, n) && pass2;
	pass2 = testGaussJordanTransform (gf2, m, n) && pass2;

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
