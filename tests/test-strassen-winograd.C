/* tests/test-strassen-winograd.C
 * Copyright 2011 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test for implementation of Strassen-Winograd fast
 * matrix-multiplication
 */

#include <iostream>

#include "test-common.h"

#include <linbox/blas/context.h>
#include <linbox/ring/gf2.h>
#include <linbox/ring/modular.h>
#include <linbox/matrix/dense.h>
#include <linbox/vector/stream.h>
#include <linbox/algorithms/strassen-winograd.h>
#include <linbox/blas/level3.h>

using namespace LinBox;

template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
bool testMul (Context<Ring, Modules> &ctx, const Matrix1 &A, const Matrix2 &B, const Matrix3 &C)
{
	commentator.start ("Testing StrassenWinograd::gemm (b = 0)", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	typename Matrix3::ContainerType Cp1 (C.rowdim (), C.coldim ()), Cp2 (C.rowdim (), C.coldim ());

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "B = " << std::endl;
	BLAS3::write (ctx, report, B);

	typename Ring::Element a;
	NonzeroRandIter<Ring> nri (ctx.F, typename Ring::RandIter (ctx.F));
	nri.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	StrassenWinograd<GenericModule::Tag> sw;

	commentator.start ("Strassen-Winograd multiplication");

	sw.gemm (ctx.F, ctx.M, a, A, B, ctx.F.zero (), Cp1);

	commentator.stop (MSG_DONE);

	report << "(From Strassen-Winograd) a * A * B = " << std::endl;
	BLAS3::write (ctx, report, Cp1);

	commentator.start ("Classical multiplication");

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), Cp2);

	commentator.stop (MSG_DONE);

	report << "(From classical method) a * A * B = " << std::endl;
	BLAS3::write (ctx, report, Cp2);

	if (!BLAS3::equal (ctx, Cp1, Cp2)) {
		error << "ERROR: Results differ" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
bool testAddMul (Context<Ring, Modules> &ctx, const Matrix1 &A, const Matrix2 &B, const Matrix3 &C)
{
	commentator.start ("Testing StrassenWinograd::gemm (b != 0)", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	typename Matrix3::ContainerType Cp1 (C.rowdim (), C.coldim ()), Cp2 (C.rowdim (), C.coldim ());

	BLAS3::copy (ctx, C, Cp1);
	BLAS3::copy (ctx, C, Cp2);

	report << "A = " << std::endl;
	BLAS3::write (ctx, report, A);

	report << "B = " << std::endl;
	BLAS3::write (ctx, report, B);

	report << "C(1) = " << std::endl;
	BLAS3::write (ctx, report, Cp1);

	report << "C(2) = " << std::endl;
	BLAS3::write (ctx, report, Cp2);

	typename Ring::Element a, b;
	NonzeroRandIter<Ring> nri (ctx.F, typename Ring::RandIter (ctx.F));

	nri.random (a);
	nri.random (b);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	report << "Coefficient b = ";
	ctx.F.write (report, b) << std::endl;

	StrassenWinograd<GenericModule::Tag> sw;

	commentator.start ("Strassen-Winograd multiplication");

	sw.gemm (ctx.F, ctx.M, a, A, B, b, Cp1);

	commentator.stop (MSG_DONE);

	report << "(From Strassen-Winograd) a * A * B + b * C = " << std::endl;
	BLAS3::write (ctx, report, Cp1);

	commentator.start ("Classical multiplication");

	BLAS3::gemm (ctx, a, A, B, b, Cp2);

	commentator.stop (MSG_DONE);

	report << "(From classical method) a * A * B + b * C = " << std::endl;
	BLAS3::write (ctx, report, Cp2);

	if (!BLAS3::equal (ctx, Cp1, Cp2)) {
		error << "ERROR: Results differ" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass1 = true, pass2 = true;

	static long m = 128;
	static long k = 128;
	static long n = 128;
	static integer q = 101U;

	static Argument args[] = {
		{ 'm', "-m M", "Set row-dimension of matrix A to M.", TYPE_INT, &m },
		{ 'k', "-k K", "Set column-dimension of matrix A to K.", TYPE_INT, &k },
		{ 'n', "-n N", "Set column-dimension of matrix B to N.", TYPE_INT, &n },
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

	commentator.start ("Strassen-Winograd test suite", "Strassen-Winograd");

	std::ostringstream str;
	str << "Running tests over GF(" << q << ")" << std::ends;

	Context<Modular<uint32> > ctx_GFq (GFq);

	RandomDenseStream<Modular<uint32>, DenseMatrix<uint32>::Row> stream_A_gfq (GFq, m, k);
	RandomDenseStream<Modular<uint32>, DenseMatrix<uint32>::Row> stream_B_gfq (GFq, k, n);
	RandomDenseStream<Modular<uint32>, DenseMatrix<uint32>::Row> stream_C_gfq (GFq, m, n);

	DenseMatrix<uint32> A_gfq (stream_A_gfq), B_gfq (stream_B_gfq), C_gfq (stream_C_gfq);

	commentator.start (str.str ().c_str (), "Strassen-Winograd");

	pass1 = testMul (ctx_GFq, A_gfq, B_gfq, C_gfq) && pass1;
	pass1 = testAddMul (ctx_GFq, A_gfq, B_gfq, C_gfq) && pass1;

	commentator.stop (MSG_STATUS (pass1));

	commentator.start ("Running tests over GF(2)", "Strassen-Winograd");

	Context<GF2> ctx_gf2 (gf2);

	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream_A_gf2 (gf2, m, k);
	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream_B_gf2 (gf2, k, n);
	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream_C_gf2 (gf2, m, n);

	DenseMatrix<bool> A_gf2 (stream_A_gf2), B_gf2 (stream_B_gf2), C_gf2 (stream_C_gf2);

	pass1 = testMul (ctx_gf2, A_gf2, B_gf2, C_gf2) && pass1;
	pass1 = testAddMul (ctx_gf2, A_gf2, B_gf2, C_gf2) && pass1;

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
