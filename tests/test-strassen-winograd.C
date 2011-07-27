/* tests/test-strassen-winograd.C
 * Copyright 2011 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test for implementation of Strassen-Winograd fast
 * matrix-multiplication
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
#include <lela/algorithms/strassen-winograd.h>
#include <lela/blas/level3.h>

using namespace LELA;

template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
bool testMul (Context<Ring, Modules> &ctx, const Matrix1 &A, const Matrix2 &B, const Matrix3 &C)
{
	commentator.start ("Testing StrassenWinograd::gemm (b = 0)", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	typename Matrix3::ContainerType Cp1 (C.rowdim (), C.coldim ()), Cp2 (C.rowdim (), C.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), Cp1);
	BLAS3::scal (ctx, ctx.F.zero (), Cp2);

	reportUI << "A = " << std::endl;
	BLAS3::write (ctx, reportUI, A);

	reportUI << "B = " << std::endl;
	BLAS3::write (ctx, reportUI, B);

	typename Ring::Element a;
	NonzeroRandIter<Ring> nri (ctx.F, typename Ring::RandIter (ctx.F));
	nri.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	StrassenWinograd<typename GenericModule<Ring>::Tag> sw (1);

	commentator.start ("Strassen-Winograd multiplication");

	sw.gemm (ctx.F, ctx.M, a, A, B, ctx.F.zero (), Cp1);

	commentator.stop (MSG_DONE);

	reportUI << "(From Strassen-Winograd) a * A * B = " << std::endl;
	BLAS3::write (ctx, reportUI, Cp1);

	commentator.start ("Classical multiplication");

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), Cp2);

	commentator.stop (MSG_DONE);

	reportUI << "(From classical method) a * A * B = " << std::endl;
	BLAS3::write (ctx, reportUI, Cp2);

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
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	typename Matrix3::ContainerType Cp1 (C.rowdim (), C.coldim ()), Cp2 (C.rowdim (), C.coldim ());

	BLAS3::copy (ctx, C, Cp1);
	BLAS3::copy (ctx, C, Cp2);

	reportUI << "A = " << std::endl;
	BLAS3::write (ctx, reportUI, A);

	reportUI << "B = " << std::endl;
	BLAS3::write (ctx, reportUI, B);

	reportUI << "C(1) = " << std::endl;
	BLAS3::write (ctx, reportUI, Cp1);

	reportUI << "C(2) = " << std::endl;
	BLAS3::write (ctx, reportUI, Cp2);

	typename Ring::Element a, b;
	NonzeroRandIter<Ring> nri (ctx.F, typename Ring::RandIter (ctx.F));

	nri.random (a);
	nri.random (b);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	report << "Coefficient b = ";
	ctx.F.write (report, b) << std::endl;

	StrassenWinograd<typename GenericModule<Ring>::Tag> sw (1);

	commentator.start ("Strassen-Winograd multiplication");

	sw.gemm (ctx.F, ctx.M, a, A, B, b, Cp1);

	commentator.stop (MSG_DONE);

	reportUI << "(From Strassen-Winograd) a * A * B + b * C = " << std::endl;
	BLAS3::write (ctx, reportUI, Cp1);

	commentator.start ("Classical multiplication");

	BLAS3::gemm (ctx, a, A, B, b, Cp2);

	commentator.stop (MSG_DONE);

	reportUI << "(From classical method) a * A * B + b * C = " << std::endl;
	BLAS3::write (ctx, reportUI, Cp2);

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

	static long m = 100;
	static long k = 103;
	static long n = 99;
	static integer q = 101U;

	static Argument args[] = {
		{ 'm', "-m M", "Set row-dimension of matrix A to M.", TYPE_INT, &m },
		{ 'k', "-k K", "Set column-dimension of matrix A to K.", TYPE_INT, &k },
		{ 'n', "-n N", "Set column-dimension of matrix B to N.", TYPE_INT, &n },
		{ 'q', "-q Q", "Operate over the ring ZZ/Q [1] for float modulus.", TYPE_INTEGER, &q },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Modular<float> GFq (q);
	GF2 gf2;

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("Strassen-Winograd test suite", "Strassen-Winograd");

	std::ostringstream str;
	str << "Running tests over GF(" << q << ")" << std::ends;

	Context<Modular<float> > ctx_GFq (GFq);

	RandomDenseStream<Modular<float>, DenseMatrix<float>::Row> stream_A_gfq (GFq, k, m);
	RandomDenseStream<Modular<float>, DenseMatrix<float>::Row> stream_B_gfq (GFq, n, k);
	RandomDenseStream<Modular<float>, DenseMatrix<float>::Row> stream_C_gfq (GFq, n, m);

	DenseMatrix<float> A_gfq (stream_A_gfq), B_gfq (stream_B_gfq), C_gfq (stream_C_gfq);

	commentator.start (str.str ().c_str (), "Strassen-Winograd");

	pass1 = testMul (ctx_GFq, A_gfq, B_gfq, C_gfq) && pass1;
	pass1 = testAddMul (ctx_GFq, A_gfq, B_gfq, C_gfq) && pass1;

	commentator.stop (MSG_STATUS (pass1));

	commentator.start ("Running tests over GF(2)", "Strassen-Winograd");

	Context<GF2> ctx_gf2 (gf2);

	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream_A_gf2 (gf2, k, m);
	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream_B_gf2 (gf2, n, k);
	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream_C_gf2 (gf2, n, m);

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
