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

#include <linbox/field/modular.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/vector/stream.h>
#include <linbox/algorithms/gauss-jordan.h>

using namespace LinBox;

typedef GF2 Field;

bool testDenseGJ (const Field &F, size_t m, size_t n)
{
	commentator.start ("Testing GaussJordan::DenseRowEchelonForm", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	GaussJordan<Field>::DenseMatrix R (m, n);
	GaussJordan<Field>::DenseMatrix U (m, m);
	GaussJordan<Field>::DenseMatrix UPA (m, n);

	RandomDenseStream<Field, GaussJordan<Field>::DenseMatrix::Row> A_stream (F, n, m);

	GaussJordan<Field>::DenseMatrix A (A_stream);

	report << "Number of words per row: " << A.rowBegin ()->word_size () << std::endl;

	GaussJordan<Field>::Permutation P;

	GaussJordan<Field> GJ (F);
	size_t rank;
	Field::Element det;

	MatrixDomain<Field> MD (F);

	MD.copy (R, A);

	GJ.DenseRowEchelonForm (R, U, P, rank, det);

	report << "A = " << std::endl;
	MD.write (report, A);

	report << "U = " << std::endl;
	MD.write (report, U);

	report << "P = ";
	MD.writePermutation (report, P.begin (), P.end ()) << std::endl;

	report << "R = " << std::endl;
	MD.write (report, R);

	MD.permuteRows (A, P.begin (), P.end ());

	report << "PA = " << std::endl;
	MD.write (report, A);

	MD.gemm (F.one (), U, A, F.zero (), UPA);

	report << "UPA = " << std::endl;
	MD.write (report, UPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!MD.areEqual (UPA, R)) {
		error << "ERROR: UPA != R" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

bool testStandardGJ (const Field &F, size_t m, size_t n, size_t k)
{
	commentator.start ("Testing GaussJordan::StandardRowEchelonForm", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	MatrixDomain<Field> MD (F);
	GaussJordan<Field> GJ (F);

	typedef GaussJordan<Field>::SparseMatrix Matrix;

	RandomSparseStream<Field, Matrix::Row, Field::RandIter, VectorCategories::HybridZeroOneVectorTag> A_stream (F, (double) k / (double) n, n, m);

	Matrix A (A_stream), Aorig (m, n);
	GaussJordan<Field>::DenseMatrix U (m, m);
	GaussJordan<Field>::DenseMatrix UPA (m, n);

	MD.copy (Aorig, A);

	GaussJordan<Field>::Permutation P;

	size_t rank;
	Field::Element det;

	report << "A = " << std::endl;
	MD.write (report, A, FORMAT_PRETTY);

	GJ.StandardRowEchelonForm (A, U, P, rank, det, true, true);

	report << "R = " << std::endl;
	MD.write (report, A, FORMAT_PRETTY);

	report << "U = " << std::endl;
	MD.write (report, U);

	report << "P = ";
	MD.writePermutation (report, P.begin (), P.end ()) << std::endl;

	MD.permuteRows (Aorig, P.begin (), P.end ());
	report << "PA = " << std::endl;
	MD.write (report, Aorig, FORMAT_PRETTY);

	MD.gemm (F.one (), U, Aorig, F.zero (), UPA);
	report << "UPA = " << std::endl;
	MD.write (report, UPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!MD.areEqual (UPA, A)) {
		error << "UPA != R, not okay" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long m = 96;
	static long n = 96;
	static long k = 10;
	static int iterations = 1;

	static Argument args[] = {
		{ 'm', "-m M", "Set row-dimension of matrix A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set column-dimension of matrix A to N.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Field F;

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("GaussJordan test suite", "GaussJordan");

	pass = testDenseGJ (F, m, n) && pass;
	pass = testStandardGJ (F, m, n, k) && pass;

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
