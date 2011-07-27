/* tests/test-blas-cblas-module.C
 * Copyright 2011 Florian Fischer
 *
 * Written by Florian Fischer <florian-fischer@gmx.net>
 *
 * Test suite for BLAS routines using CBLAS-module
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/lela-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "lela/util/commentator.h"
#include "lela/ring/gf2.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"
#include "lela/blas/level2.h"
#include "lela/blas/level3.h"
#include "lela/vector/stream.h"
#include "lela/matrix/dense.h"

#include "test-common.h"
#include "test-blas-level3.h"

using namespace std;
using namespace LELA;

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2>
bool testscalcblasConsistency  (LELA::Context<Ring, Modules1> &ctx1,
				LELA::Context<Ring, Modules2> &ctx2,
				const char *text,
				const typename Ring::Element &a,
				const Matrix1 &A1, 
				const Matrix2 &A2)
{

	ostringstream str;
	str << "Testing " << text << " scal consistency" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	bool pass = true;

	typename Matrix1::ContainerType A3 (A1.rowdim (), A1.coldim ());
	typename Matrix2::ContainerType A4 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1, A1, A3);
	BLAS3::copy(ctx1, A1, A4);

	report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A3);

	BLAS3::scal (ctx1, a, A3);

	report << "Matrix a A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A3);

	report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx2, report, A4);

	BLAS3::scal (ctx2, a, A4);

	report << "Matrix a A_2: " << std::endl;;
	BLAS3::write (ctx2, report, A4);

	if (!BLAS3::equal (ctx1, A3, A4))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a A_1 !=  a A_2 " << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2,  class Matrix3, class Matrix4>
bool testaxpycblasConsistency  (LELA::Context<Ring, Modules1> &ctx1,
				LELA::Context<Ring, Modules2> &ctx2,
				const char *text,
				const typename Ring::Element &a,
				const Matrix1 &A1, const Matrix2 &A2,
				const Matrix3 &A3, const Matrix4 &A4)
{

        ostringstream str;
        str << "Testing " << text << " axpy consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

        bool pass = true;

	typename Matrix3::ContainerType A6 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1, A1, A6);

        typename Matrix2::ContainerType A7 (A2.rowdim (), A2.coldim ());
        typename Matrix4::ContainerType A8 (A2.rowdim (), A2.coldim ());

	BLAS3::copy(ctx1, A2, A7);
	BLAS3::copy(ctx1, A2, A8);

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

        report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A7);

	BLAS3::axpy (ctx1, a, A1, A7);

        report << "Matrix a A_1 + A_2: " << std::endl;;
	BLAS3::write (ctx1, report, A7);

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A6);

        report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A8);

	BLAS3::axpy (ctx2, a, A6, A8);

        report << "Matrix a A_3 + A_4: " << std::endl;
	BLAS3::write (ctx2, report, A8);

        if (!BLAS3::equal (ctx1, A7, A8))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a A_1 + A_2  !=  a A_3 + A_4 " << std::endl;
		pass = false;
	}

        commentator.stop (MSG_STATUS (pass));

        return pass;
}


template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2,  class Matrix3, class Matrix4, class Matrix5, class Matrix6>
bool testgemmcblasConsistency  (LELA::Context<Ring, Modules1> &ctx1,
				LELA::Context<Ring, Modules2> &ctx2,
				const char *text,
				const typename Ring::Element &a, const typename Ring::Element &b,
				const Matrix1 &A1, const Matrix2 &A2,
				const Matrix3 &A3, const Matrix4 &A4,
				const Matrix5 &A5, const Matrix6 &A6)
{

        ostringstream str;
        str << "Testing " << text << " gemm consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

        bool pass = true;

        typename Matrix4::ContainerType A7 (A1.rowdim (), A1.coldim ());
	typename Matrix5::ContainerType A8 (A2.rowdim (), A2.coldim ());

	BLAS3::copy(ctx1, A1, A7);
	BLAS3::copy(ctx1, A2, A8);

        typename Matrix3::ContainerType A9 (A3.rowdim (), A3.coldim ());
	typename Matrix6::ContainerType A10 (A3.rowdim (), A3.coldim ());

	BLAS3::copy(ctx1, A3, A9);
	BLAS3::copy(ctx1, A3, A10);

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

        report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A2);

	report << "Coefficient b: ";
	ctx1.F.write (report, b) << std::endl;

	report << "Matrix A_3: "<< std::endl;
	BLAS3::write (ctx1, report, A9);

	BLAS3::gemm (ctx1, ctx1.F.one (), A1, A2, b, A9);

        report << "Matrix a A_1 A_2 + b A_3: " << std::endl;
	BLAS3::write (ctx1, report, A9);

        report << "Matrix A_4: "<< std::endl;
	BLAS3::write (ctx1, report, A7);

        report << "Matrix A_5: "<< std::endl;
	BLAS3::write (ctx1, report, A8);

        report << "Matrix A_6: "<< std::endl;
	BLAS3::write (ctx1, report, A10);

	BLAS3::gemm (ctx2, ctx2.F.one (), A7, A8, b, A10);

        report << "Matrix a A_4 A_5 + b A_6: " << std::endl;
	BLAS3::write (ctx2, report, A10);

        if (!BLAS3::equal (ctx1, A9, A10))
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a A_1 A_2 + b A_3  !=  a A_4 A_5 + b A_6 " << std::endl;
			pass = false;
		}

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring1, class Ring2, class Modules1, class Modules2>
bool testBLASModulesConsistency (const Ring1 &R, LELA::Context<Ring2, Modules1> &ctx1, LELA::Context<Ring2, Modules2> &ctx2, const char *text, size_t m, size_t n, size_t p) 
{
	std::ostringstream str;
	str << "Testing BLASModules consistency over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Ring1, typename DenseMatrix<typename Ring1::Element>::Row> stream11 (R, n, m);
	DenseMatrix<typename Ring1::Element> M1 (stream11);
	TransposeMatrix<DenseMatrix<typename Ring1::Element> > M2 (M1);

        RandomDenseStream<Ring1, typename DenseMatrix<typename Ring1::Element>::Row> stream111 (R, m, n);
        DenseMatrix<typename Ring1::Element> M3 (stream111);
        TransposeMatrix<DenseMatrix<typename Ring1::Element> > M4 (M3);
       
	typename Ring1::Element a, b;
        NonzeroRandIter<Ring1, typename Ring1::RandIter> r (R, typename Ring1::RandIter (R));

        r.random (a);
        r.random (b);

	pass = testscalcblasConsistency (ctx1, ctx2, "dense", a, M1, M1) && pass;
	pass = testscalcblasConsistency (ctx1, ctx2, "dense(col-wise)", a, M2, M2) && pass;

	pass = testaxpycblasConsistency (ctx1, ctx2, "dense /dense", a, M1, M1, M1, M1) && pass;
	pass = testaxpycblasConsistency (ctx1, ctx2, "dense(col-wise)/dense", a, M2, M3, M2, M3) &&pass;
	pass = testaxpycblasConsistency (ctx1, ctx2, "dense(col-wise)/dense(col-wise)", a, M2, M2, M2, M2) &&pass;

	RandomDenseStream<Ring1, typename DenseMatrix<typename Ring1::Element>::Row> stream12 (R, n, m);
	DenseMatrix<typename Ring1::Element> A1_dense (stream12);
	TransposeMatrix<DenseMatrix<typename Ring1::Element> > A1_trans (A1_dense);

	RandomDenseStream<Ring1, typename DenseMatrix<typename Ring1::Element>::Row> stream13 (R, p, n);
	DenseMatrix<typename Ring1::Element> A2_dense (stream13);
	TransposeMatrix<DenseMatrix<typename Ring1::Element> > A2_trans (A2_dense);

	RandomDenseStream<Ring1, typename DenseMatrix<typename Ring1::Element>::Row> stream14 (R, p, m);
	DenseMatrix<typename Ring1::Element> A3_dense (stream14);
	TransposeMatrix<DenseMatrix<typename Ring1::Element> > A3_trans (A3_dense);

	RandomDenseStream<Ring1, typename DenseMatrix<typename Ring1::Element>::Row> stream15 (R, m, p);
	DenseMatrix<typename Ring1::Element> A4_dense (stream15);
	TransposeMatrix<DenseMatrix<typename Ring1::Element> > A4_trans (A4_dense);

	pass = testgemmcblasConsistency (ctx1, ctx2, "dense/dense/dense", a, b, A1_dense, A2_dense, A3_dense, A1_dense, A2_dense, A3_dense) && pass;
	pass = testgemmcblasConsistency (ctx1, ctx2, "dense(col-wise)/dense/dense", a, b, A1_trans, A3_dense, A2_dense, A1_trans, A3_dense, A2_dense)  && pass;
	pass = testgemmcblasConsistency (ctx1, ctx2, "dense(col-wise)/dense(col-wise)/dense", a, b, A1_trans, A4_trans, A2_dense, A1_trans, A4_trans, A2_dense) && pass;
	pass = testgemmcblasConsistency (ctx1, ctx2, "dense/dense(col-wise)/dense", a, b, A3_dense, A2_trans, A1_dense, A3_dense, A2_trans, A1_dense) && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long l = 100;
	static long m = 121;
	static long n = 114;
	static long p = 137;
	static integer q_float = 101;
	static integer q_double = 65521;
	static int iterations = 1;

	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set row-dimension of matrix C and column-dimension of B to N.", TYPE_INT, &n },
		{ 'p', "-p P", "Set column-dimension of matrix C to P.", TYPE_INT, &p },
		{ 'q', "-q Q", "Operate over the ring Z/Q for float modulus", TYPE_INTEGER, &q_float },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (7);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (7);

	commentator.start ("BLAS BLASModule test-suite", "cBLASModule");

#ifdef __LELA_BLAS_AVAILABLE
	TypeWrapperRing<float> R_TW_float;
	Modular<float> F_float (q_float);
	Context<TypeWrapperRing<float>, BLASModule<float> > ctx_float (R_TW_float);
	Context<TypeWrapperRing<float>, GenericModule<TypeWrapperRing<float> > > ctx_float_gen (R_TW_float);

	TypeWrapperRing<double> R_TW_double;
	Modular<double> F_double (q_double);
        Context<TypeWrapperRing<double>, BLASModule<double> > ctx_double (R_TW_double);
        Context<TypeWrapperRing<double>, GenericModule<TypeWrapperRing<double> > > ctx_double_gen (R_TW_double);

	pass = testBLASModulesConsistency(F_float, ctx_float, ctx_float_gen, "Modular<float>", m, n, p) && pass;
        pass = testBLASModulesConsistency(F_double, ctx_double, ctx_double_gen, "Modular<double>", m, n, p) && pass;
#else // !__LELA_BLAS_AVAILABLE
	commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
		<< "LELA was not configured to use BLAS, so these tests cannot be run. Skipping." << std::endl;
#endif // __LELA_BLAS_AVAILABLE

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
