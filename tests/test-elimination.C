/* tests/test-elimination.C
 * Copyright 2011 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test for naive Gaussian elimination
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
#include <lela/matrix/sparse.h>
#include <lela/vector/stream.h>
#include <lela/algorithms/elimination.h>

using namespace LELA;

template <class Ring, class Matrix>
bool testEchelonize (const Ring &F, const char *text, Matrix &A)
{
	std::ostringstream str;
	str << "Testing Elimination::echelonize for " << text << " matrices" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Context<Ring> ctx (F);

	Elimination<Ring> elim (ctx);

	typename Matrix::ContainerType PA (A.rowdim (), A.coldim ()), LPA (A.rowdim (), A.coldim ());

	BLAS3::copy (ctx, A, PA);

	typename Elimination<Ring>::Permutation P;

	size_t rank;
	typename Ring::Element det;

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	elim.echelonize (A, P, rank, det, true);

	report << "L, R = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	BLAS3::permute_rows (ctx, P.begin (), P.end (), PA);
	report << "PA = " << std::endl;
	BLAS3::write (ctx, report, PA, FORMAT_PRETTY);

	typename Matrix::ContainerType L (A.rowdim (), A.rowdim ());

	typename Matrix::ContainerType::RowIterator i_L;
	StandardBasisStream<Ring, typename Matrix::ContainerType::Row> s (ctx.F, A.rowdim ());

	for (i_L = L.rowBegin (); i_L != L.rowEnd (); ++i_L)
		s >> *i_L;

	elim.move_L (L, A);

	BLAS3::scal (ctx, F.zero (), LPA);
	BLAS3::gemm (ctx, F.one (), L, PA, F.zero (), LPA);
	report << "LPA = " << std::endl;
	BLAS3::write (ctx, report, LPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!BLAS3::equal (ctx, LPA, A)) {
		error << "LPA != R, not okay" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Matrix>
bool testEchelonizeReduced (const Ring &F, const char *text, Matrix &A)
{
	std::ostringstream str;
	str << "Testing Elimination::echelonize_reduced for " << text << " matrices" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Context<Ring> ctx (F);
	Elimination<Ring> elim (ctx);

	typename Matrix::ContainerType Acopy (A.rowdim (), A.coldim ()), LPA (A.rowdim (), A.coldim ()), L (A.rowdim (), A.rowdim ());

	BLAS3::copy (ctx, A, Acopy);

	typename Elimination<Ring>::Permutation P;

	size_t rank;
	typename Ring::Element det;

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	elim.echelonize_reduced (A, L, P, rank, det, true);

	report << "R = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	report << "L = " << std::endl;
	BLAS3::write (ctx, report, L, FORMAT_PRETTY);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	BLAS3::permute_rows (ctx, P.begin (), P.end (), Acopy);
	report << "PA = " << std::endl;
	BLAS3::write (ctx, report, Acopy, FORMAT_PRETTY);

	BLAS3::scal (ctx, F.zero (), LPA);
	BLAS3::gemm (ctx, F.one (), L, Acopy, F.zero (), LPA);

	report << "LPA = " << std::endl;
	BLAS3::write (ctx, report, LPA);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!BLAS3::equal (ctx, LPA, A)) {
		error << "LPA != R, not okay" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Matrix>
bool testPLUQ (const Ring &F, const char *text, Matrix &A)
{
	std::ostringstream str;
	str << "Testing Elimination::pluq for " << text << " matrices" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Context<Ring> ctx (F);

	Elimination<Ring> elim (ctx);

	typename Matrix::ContainerType Acopy (A.rowdim (), A.coldim ()), L (A.rowdim (), A.rowdim ());

	BLAS3::copy (ctx, A, Acopy);

	typename Elimination<Ring>::Permutation P, Q;

	size_t rank;
	typename Ring::Element det;

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	elim.pluq (A, P, Q, rank, det);

	report << "L, U = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	report << "P = ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;

	report << "Q = ";
	BLAS1::write_permutation (report, Q.begin (), Q.end ()) << std::endl;

	BLAS3::scal (ctx, ctx.F.zero (), L);
	elim.move_L (L, A);

	report << "L = " << std::endl;
	BLAS3::write (ctx, report, L, FORMAT_PRETTY);

	report << "U = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	BLAS3::trmm (ctx, F.one (), L, A, LowerTriangular, true);

	report << "LU = " << std::endl;
	BLAS3::write (ctx, report, A, FORMAT_PRETTY);

	BLAS3::permute_rows (ctx, P.begin (), P.end (), A);

	report << "PLU = " << std::endl;
	BLAS3::write (ctx, report, A);

	BLAS3::permute_cols (ctx, Q.begin (), Q.end (), A);

	report << "PLUQ = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "Computed rank = " << rank << std::endl;
	report << "Computed det = ";
	F.write (report, det);
	report << std::endl;

	if (!BLAS3::equal (ctx, A, Acopy)) {
		error << "PLUQ != A, not okay" << std::endl;

		BLAS3::axpy (ctx, ctx.F.minusOne (), A, Acopy);

		error << "Difference is:" << std::endl;
		BLAS3::write (ctx, error, Acopy);

		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass1 = true, pass2 = true;

	static long m = 100;
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

	commentator.start ("Elimination test suite", "Elimination");

	std::ostringstream str;
	str << "Running tests over GF(" << q << ")" << std::ends;

	commentator.start (str.str ().c_str (), "Elimination");

	RandomDenseStream<Modular<uint32>, DenseMatrix<uint32>::Row> A1_stream (GFq, n, m);
	RandomSparseStream<Modular<uint32>, SparseMatrix<uint32>::Row> A2_stream (GFq, (double) k / (double) n, n, m);

	DenseMatrix<uint32> A1 (A1_stream);
	SparseMatrix<uint32> A2 (A2_stream);

	pass1 = testEchelonize (GFq, "dense", A1) && pass1;
	pass1 = testEchelonize (GFq, "sparse", A2) && pass1;

	A1_stream.reset ();
	A2_stream.reset ();

	DenseMatrix<uint32> A3 (A1_stream);
	SparseMatrix<uint32> A4 (A2_stream);

	pass1 = testEchelonizeReduced (GFq, "dense", A3) && pass1;
	pass1 = testEchelonizeReduced (GFq, "sparse", A4) && pass1;

	A1_stream.reset ();
	A2_stream.reset ();

	DenseMatrix<uint32> A5 (A1_stream);
	SparseMatrix<uint32> A6 (A2_stream);

	pass1 = testPLUQ (GFq, "dense", A5) && pass1;
	pass1 = testPLUQ (GFq, "sparse", A6) && pass1;

	commentator.stop (MSG_STATUS (pass1));

	commentator.start ("Running tests over GF(2)", "Elimination");

	RandomDenseStream<GF2, DenseMatrix<bool>::Row> B1_stream (gf2, n, m);
	RandomSparseStream<GF2, Vector<GF2>::Sparse> B2_stream (gf2, (double) k / (double) n, n, m);
	RandomHybridStream<GF2, Vector<GF2>::Hybrid> B3_stream (gf2, (double) k / (double) n, n, m);

	DenseMatrix<bool> B1 (B1_stream);
	SparseMatrix<bool, Vector<GF2>::Sparse> B2 (B2_stream);
	SparseMatrix<bool, Vector<GF2>::Hybrid> B3 (B3_stream);

	pass2 = testEchelonize (gf2, "dense", B1) && pass2;
	pass2 = testEchelonize (gf2, "sparse", B2) && pass2;
	pass2 = testEchelonize (gf2, "hybrid", B3) && pass2;

	B1_stream.reset ();
	B2_stream.reset ();
	B3_stream.reset ();

	DenseMatrix<bool> B4 (B1_stream);
	SparseMatrix<bool, Vector<GF2>::Sparse> B5 (B2_stream);
	SparseMatrix<bool, Vector<GF2>::Hybrid> B6 (B3_stream);

	pass2 = testEchelonizeReduced (gf2, "dense", B4) && pass2;
	pass2 = testEchelonizeReduced (gf2, "sparse", B5) && pass2;
	pass2 = testEchelonizeReduced (gf2, "hybrid", B6) && pass2;

	B1_stream.reset ();
	B2_stream.reset ();
	DenseMatrix<bool> B7 (B1_stream);
	SparseMatrix<bool, Vector<GF2>::Sparse> B8 (B2_stream);

	pass2 = testPLUQ (gf2, "dense", B7) && pass2;
	pass2 = testPLUQ (gf2, "sparse", B8) && pass2;

	// Note: PLUQ not allowed on hybrid vectors

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
