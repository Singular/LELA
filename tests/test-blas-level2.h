 /* tests/test-blas-level2.h
 * Copyright 2001, 2002, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic test suite for BLAS level 2 routines
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_BLAS_LEVEL2_H
#define __LELA_TESTS_TEST_BLAS_LEVEL2_H

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

using namespace std;
using namespace LELA;

typedef std::vector<std::pair<uint32, uint32> > Permutation;

template <class Matrix>
TransposeMatrix<Matrix> transpose (Matrix &M)
{
	return TransposeMatrix<Matrix> (M);
}

/* Test 6: ger and gemm
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix, class Vector1, class Vector2>
static bool testGerGemm (Context<Field, Modules> &ctx, const char *text, const Matrix &A, const Vector1 &u, const Vector2 &v)
{
	ostringstream str;

	str << "Testing " << text << " ger (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	bool ret = true;

	DenseMatrix<typename Field::Element> A2 (A.rowdim (), A.coldim ());

	typename Matrix::ContainerType A1 (A.rowdim (), A.coldim ());

	BLAS3::copy (ctx, A, A1);
	BLAS3::copy (ctx, A, A2);

	DenseMatrix<typename Field::Element> U (A.rowdim (), 1);
	DenseMatrix<typename Field::Element> V (1, A.coldim ());

	BLAS1::copy (ctx, u, *U.colBegin ());
	BLAS1::copy (ctx, v, *V.rowBegin ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A1);

	reportUI << "Input vector u: ";
	BLAS1::write (ctx, reportUI, u) << endl;

	reportUI << "Input vector v: ";
	BLAS1::write (ctx, reportUI, v) << endl;

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS2::ger (ctx, a, u, v, A1);

	reportUI << "Output matrix a * u * v^T + A: " << std::endl;
	BLAS3::write (ctx, reportUI, A1);

	BLAS3::gemm (ctx, a, U, V, ctx.F.one (), A2);

	reportUI << "Output matrix a * U * V + A: " << std::endl;
	BLAS3::write (ctx, reportUI, A2);

	if (!BLAS3::equal (ctx, A1, A2)) {
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported a * u * v^T + A != a * U * V + A" << endl;

		BLAS3::axpy (ctx, ctx.F.minusOne (), A1, A2);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx, error, A2);

		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 9: gemv and gemm
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixIteratorTypes::Generic, MatrixIteratorTypes::Col) 
{
	lela_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemv (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> ABgemm (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> ABgemv (A.rowdim (), B.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), ABgemm);
	BLAS3::scal (ctx, ctx.F.zero (), ABgemv);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A);

	reportUI << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	typename DenseMatrix<typename Field::Element>::ColIterator i_AB = ABgemv.colBegin ();
	typename Matrix2::ConstColIterator i_B = B.colBegin ();

	for (; i_AB != ABgemv.colEnd (); ++i_B, ++i_AB)
		BLAS2::gemv (ctx, a, A, *i_B, ctx.F.zero (), *i_AB);

	reportUI << "Output matrix AB (from gemv):" << endl;
	BLAS3::write (ctx, reportUI, ABgemv);

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), ABgemm);

	reportUI << "Output matrix AB (from gemm):" << endl;
	BLAS3::write (ctx, reportUI, ABgemm);

	if (!BLAS3::equal (ctx, ABgemv, ABgemm)) {
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported matrices from gemm and gemv are not equal" << endl;

		BLAS3::axpy (ctx, ctx.F.minusOne (), ABgemv, ABgemm);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx, error, ABgemm);

		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixIteratorTypes::Row, MatrixIteratorTypes::Generic) 
{
	lela_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemv (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> ABgemm (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> ABgemv (A.rowdim (), B.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), ABgemm);
	BLAS3::scal (ctx, ctx.F.zero (), ABgemv);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A);

	reportUI << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	typename DenseMatrix<typename Field::Element>::RowIterator i_AB = ABgemv.rowBegin ();
	typename Matrix1::ConstRowIterator i_A = A.rowBegin ();

	for (; i_AB != ABgemv.rowEnd (); ++i_A, ++i_AB)
		BLAS2::gemv (ctx, a, transpose (B), *i_A, ctx.F.zero (), *i_AB);

	reportUI << "Output matrix AB (from gemv):" << endl;
	BLAS3::write (ctx, reportUI, ABgemv);

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), ABgemm);

	reportUI << "Output matrix AB (from gemm):" << endl;
	BLAS3::write (ctx, reportUI, ABgemm);

	if (!BLAS3::equal (ctx, ABgemv, ABgemm)) {
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported matrices from gemm and gemv are not equal" << endl;

		BLAS3::axpy (ctx, ctx.F.minusOne (), ABgemv, ABgemm);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx, error, ABgemm);

		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixIteratorTypes::Generic, MatrixIteratorTypes::RowCol) 
{
	return testGemvGemmSpecialised (ctx, text, A, B, MatrixIteratorTypes::Generic (), MatrixIteratorTypes::Col ());
}

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixIteratorTypes::RowCol, MatrixIteratorTypes::Generic) 
{
	return testGemvGemmSpecialised (ctx, text, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Generic ());
}

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol) 
{
	return testGemvGemmSpecialised (ctx, text, A, B, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Generic ());
}

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemm (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B) 
{
	return testGemvGemmSpecialised (ctx, text, A, B, typename Matrix1::IteratorType (), typename Matrix2::IteratorType ());
}

/* Test 10: gemv with coefficients
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix, class Vector>
static bool testGemvCoeff (Context<Field, Modules> &ctx, const char *text, const Matrix &A, const Vector &v)
{
	ostringstream str;

	str << "Testing " << text << " gemv (coefficients)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element a, b, negainvb;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	bool ret = true;

	typename LELA::Vector<Field>::Dense aAv (A.rowdim ());

	BLAS1::scal (ctx, ctx.F.zero (), aAv);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A);

	reportUI << "Input vector v: ";
	BLAS1::write (ctx, reportUI, v) << endl;

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	ctx.F.write (report, a) << ", b = ";
	ctx.F.write (report, b) << std::endl;

	BLAS2::gemv (ctx, a, A, v, ctx.F.zero (), aAv);

	reportUI << "Output vector a * A * v: ";
	BLAS1::write (ctx, reportUI, aAv) << std::endl;

	ctx.F.inv (negainvb, a);
	ctx.F.mulin (negainvb, b);
	ctx.F.negin (negainvb);

	report << "Coefficient: -a^-1 b = ";
	ctx.F.write (report, negainvb) << std::endl;

	BLAS2::gemv (ctx, b, A, v, negainvb, aAv);

	reportUI << "Output vector w := b * A * v - b * a^-1 * aAv: ";
	BLAS1::write (ctx, reportUI, aAv) << std::endl;

	if (!BLAS1::is_zero (ctx, aAv)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS1 reported vector w is not zero" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 12: trsm and trsv
 *
 * B must be iterable by columns
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testTrsmTrsv (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, TriangularMatrixType type) 
{
	lela_check (A.coldim () == B.rowdim ());

	ostringstream str;

	const char *type_str[] = { "upper", "lower" };

	str << "Testing " << text << " trsm (consistency with trsv, " << type_str[type] << " triangular matrix)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	Matrix1 U (A.rowdim (), A.coldim ());

	BLAS3::copy (ctx, A, U);
	makeNonsingDiag (ctx.F, U, false);

	DenseMatrix<typename Field::Element> UinvBtrsm (U.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> UinvBtrsv (U.rowdim (), B.coldim ());

        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix U:" << endl;
	BLAS3::write (ctx, reportUI, U);

	reportUI << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B);

	BLAS3::copy (ctx, B, UinvBtrsm);
	BLAS3::copy (ctx, B, UinvBtrsv);

	typename DenseMatrix<typename Field::Element>::ColIterator i_UinvB = UinvBtrsv.colBegin ();

	for (; i_UinvB != UinvBtrsv.colEnd (); ++i_UinvB)
		BLAS2::trsv (ctx, U, *i_UinvB, type, false);

	reportUI << "Output matrix U^-1 * B (from trsv):" << endl;
	BLAS3::write (ctx, reportUI, UinvBtrsv);

	BLAS3::trsm (ctx, ctx.F.one (), U, UinvBtrsm, type, false);

	reportUI << "Output matrix U^-1 * B (from trsm):" << endl;
	BLAS3::write (ctx, reportUI, UinvBtrsm);

	if (!BLAS3::equal (ctx, UinvBtrsv, UinvBtrsm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported matrices from trsm and trsv are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Modules, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testBLAS2 (Context<Field, Modules> &ctx, const char *text,
		Matrix1 &M1, Matrix2 &M2,
		Vector1 &v1, Vector2 &v2,
		MatrixIteratorTypes::RowCol)
{
	ostringstream str;
	str << "Testing BLAS2 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testGerGemm (ctx, text, M1, v1, v2)) pass = false;
	if (!testGemvGemm (ctx, text, M1, M2)) pass = false;
	if (!testGemvCoeff (ctx, text, M1, v2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsmTrsv (ctx, text, M1, M2, LowerTriangular)) pass = false;
		if (!testTrsmTrsv (ctx, text, M1, M2, UpperTriangular)) pass = false;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testBLAS2 (Context<Field, Modules> &ctx, const char *text,
		Matrix1 &M1, Matrix2 &M2,
		Vector1 &v1, Vector2 &v2,
		MatrixIteratorTypes::Row) 
{
	ostringstream str;
	str << "Testing BLAS2 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

//	if (!testGerGemm (ctx, text, M1, v1, v2)) pass = false; // Needs ColIterator
	if (!testGemvGemm (ctx, text, M1, M2)) pass = false;
	if (!testGemvCoeff (ctx, text, M1, v2)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testBLAS2 (Context<Field, Modules> &ctx, const char *text,
		Matrix1 &M1, Matrix2 &M2,
		Vector1 &v1, Vector2 &v2,
		MatrixIteratorTypes::Col) 
{
	ostringstream str;
	str << "Testing BLAS2 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testGerGemm (ctx, text, M1, v1, v2)) pass = false;
	if (!testGemvGemm (ctx, text, M1, M2)) pass = false;
	if (!testGemvCoeff (ctx, text, M1, v2)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testBLAS2Submatrix (Context<Field, Modules> &ctx, const char *text,
			 const Matrix1 &M1, const Matrix2 &M2,
			 Vector1 &v1, Vector2 &v2,
			 MatrixIteratorTypes::Row) 
{
	ostringstream str;
	str << "Testing BLAS2 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

//	if (!testGerGemm (ctx, text, M1, v1, v2)) pass = false; // Needs ColIterator
	if (!testGemvGemm (ctx, text, M1, M2)) pass = false;
	if (!testGemvCoeff (ctx, text, M1, v2)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2, class Vector1, class Vector2, class Vector3, class Vector4>
bool testgemvConsistency (LELA::Context<Ring, Modules1> &ctx1,
			  LELA::Context<Ring, Modules2> &ctx2,
			  const char *text,
			  const Matrix1 &M1, 
			  const Vector1 &v1, const Vector2 &v2,
			  const Matrix2 &M2,
			  const Vector3 &v3, const Vector4 &v4)
{
	ostringstream str;
	str << "Testing " << text << " gemv consistency" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool pass = true;
	
	typename VectorTraits<Ring,Vector3>::ContainerType w1;
	typename VectorTraits<Ring,Vector4>::ContainerType w2;
	typename VectorTraits<Ring,Vector2>::ContainerType w3;

	VectorUtils::ensureDim<Ring, Vector3> (w1, M1.coldim ());
	VectorUtils::ensureDim<Ring, Vector4> (w2, M1.rowdim ());

	VectorUtils::ensureDim<Ring, Vector2> (w3, M1.rowdim ());

	BLAS1::copy (ctx1, v1, w1);
	BLAS1::copy (ctx1, v2, w2);
	BLAS1::copy (ctx1, v2, w3);

	typename Matrix2::ContainerType A (M1.rowdim (), M1.coldim ());

	BLAS3::copy(ctx1,M1,A);

	typename Ring::Element a,b;
	NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

	r.random (a);
	r.random (b);

	report << "Coefficient a: ";
	ctx1.F.write (report, a) << std::endl;

	report << "Coefficient b: ";
	ctx1.F.write (report, b) << std::endl;

	reportUI << "Vector v_1: ";
	BLAS1::write (ctx1, reportUI, v1) << std::endl;

	reportUI << "Vector v_2: ";
	BLAS1::write (ctx1, reportUI, v2) << std::endl;

	reportUI << "Matrix A: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A);

	BLAS2::gemv(ctx1, a, M1, v1, b, w3);

	reportUI << "Vector a Av_1 + bv_2: ";
	BLAS1::write (ctx1, reportUI, w3) << std::endl;

	reportUI << "Vector w_1: ";
	BLAS1::write (ctx2, reportUI, w1) << std::endl;

	reportUI << "Vector w_2: ";
	BLAS1::write (ctx2, reportUI, w2) << std::endl;

	BLAS2::gemv (ctx2, a, A, w1, b, w2);

	reportUI << "Vector a Aw_1 + bw_2: ";
	BLAS1::write (ctx1, reportUI, w2) << std::endl;

	if (!BLAS1::equal (ctx1, w3, w2))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a Av_1 + bv_2 != a Aw_1 + bw_2" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}	

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testtrmvConsistency  (LELA::Context<Ring, Modules1> &ctx1,
    			   LELA::Context<Ring, Modules2> &ctx2,
		           const char *text,
      			   const Matrix1 &A1, const Vector1 &x1,
      			   const Matrix2 &A2, const Vector2 &x2,
		           TriangularMatrixType type,
		           bool diagIsOne)
{
	ostringstream str;
	str << "Testing " << text << " trmv consistency" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool pass = true;

	typename VectorTraits<Ring,Vector1>::ContainerType w1;
	typename VectorTraits<Ring,Vector2>::ContainerType w2;

	VectorUtils::ensureDim<Ring, Vector1> (w1, A1.coldim ());
	VectorUtils::ensureDim<Ring, Vector2> (w2, A1.rowdim ());

	BLAS1::copy (ctx1, x1, w1);
	BLAS1::copy (ctx1, x1, w2);

	typename Matrix2::ContainerType A3 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1,A1,A3);

	report << "Vector x_1: ";
	BLAS1::write (ctx1, report, x1) << std::endl;

	report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

	BLAS2::trmv(ctx1, A1, w1, type, diagIsOne);

	report << "Vector A_1 x_1: ";
	BLAS1::write (ctx1, report, w1) << std::endl;

	report << "Vector w_1: ";
	BLAS1::write (ctx2, report, w1) << std::endl;

	report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A3);

	BLAS2::trmv (ctx2, A3, w2, type, diagIsOne);

	report << "Vector A_2 w_1: ";
	BLAS1::write (ctx1, report, w2) << std::endl;

	if (!BLAS1::equal (ctx1, w1, w2))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
    			<< "ERROR: A_1 x_1 !=  A_2 w_1 " << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testtrsvConsistency (LELA::Context<Ring, Modules1> &ctx1,
      			  LELA::Context<Ring, Modules2> &ctx2,
      			  const char *text,
      			  const Matrix1 &A1, const Vector1 &x1,
   			  const Matrix2 &A2, const Vector2 &x2,
      			  TriangularMatrixType type,
    			  bool diagIsOne)
{
	ostringstream str;
	str << "Testing " << text << " trsv consistency" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool pass = true;

	typename VectorTraits<Ring,Vector1>::ContainerType w1;
	typename VectorTraits<Ring,Vector2>::ContainerType w2;

	VectorUtils::ensureDim<Ring, Vector1> (w1, A1.coldim ());
	VectorUtils::ensureDim<Ring, Vector2> (w2, A1.rowdim ());

	BLAS1::copy (ctx1, x1, w1);
	BLAS1::copy (ctx1, x1, w2);

	typename Matrix2::ContainerType A3 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1,A1,A3);

	report << "Vector x_1: ";
	BLAS1::write (ctx1, report, x1) << std::endl;

	report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

	BLAS2::trsv(ctx1, A1, w1, type, diagIsOne);

	report << "Vector (A_1)^-1 x_1: ";
	BLAS1::write (ctx1, report, w1) << std::endl;

	report << "Vector w_1: ";
	BLAS1::write (ctx2, report, w1) << std::endl;

	report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A3);

	BLAS2::trsv (ctx2, A3, w2, type, diagIsOne);

	report << "Vector (A_2)^-1 w_1: ";
	BLAS1::write (ctx1, report, w2) << std::endl;

	if (!BLAS1::equal (ctx1, w1, w2))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
    			<< "ERROR: (A_1)^-1 x_1 !=  (A_2)^-1 w_1 " << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2, class Vector1, class Vector2, class Vector3, class Vector4>
bool testgerConsistency (LELA::Context<Ring, Modules1> &ctx1,
      			 LELA::Context<Ring, Modules2> &ctx2,
      			 const char *text,
      			 const Matrix1 &A1, 
			 const Vector1 &x1, const Vector2 &x2,
			 const Matrix2 &A2,
			 const Vector3 &y1, const Vector4 &y2)
{
	ostringstream str;
	str << "Testing " << text << " ger consistency" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool pass = true;

	typename VectorTraits<Ring,Vector3>::ContainerType w1;
	typename VectorTraits<Ring,Vector4>::ContainerType w2;

	VectorUtils::ensureDim<Ring, Vector3> (w1, A1.rowdim ());
	VectorUtils::ensureDim<Ring, Vector4> (w2, A1.coldim ());

	BLAS1::copy (ctx1, x1, w1);
	BLAS1::copy (ctx1, x2, w2);

	typename Ring::Element a;
	NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

	r.random (a);
        
	typename Matrix2::ContainerType A3 (A1.rowdim (), A1.coldim ());
        typename Matrix2::ContainerType A4 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1,A1,A3);
	BLAS3::copy(ctx1,A1,A4);

	reportUI << "Vector x_1: ";
	BLAS1::write (ctx1, reportUI, x1) << std::endl;

	reportUI << "Vector x_2: ";
	BLAS1::write (ctx1, reportUI, x2) << std::endl;

	reportUI << "Matrix A_1: " << std::endl;
	BLAS3::write (ctx1, reportUI, A3);

	BLAS2::ger(ctx1, a, x1, x2, A3);

	report << "Coefficient a: ";
	ctx1.F.write (report, a) << std::endl;

	reportUI << "Matrix A_1 + a x_1 y^T: " << std::endl;
	BLAS3::write (ctx1, reportUI, A3);

	reportUI << "Vector w_1: ";
	BLAS1::write (ctx2, reportUI, w1) << std::endl;

	reportUI << "Vector w_2: ";
	BLAS1::write (ctx2, reportUI, w2) << std::endl;

	reportUI << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A4);

	BLAS2::ger (ctx2, a, w1, w2, A4);

	reportUI << "Matrix A_2 + a x_2 y^T: " << std::endl;
	BLAS3::write (ctx1, reportUI, A4);

	if (!BLAS3::equal (ctx1, A3, A4))
	{
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
    			<< "ERROR: A_1 + a x_1 y^T  !=  A_2 + a x_2  y^T " << std::endl;

		BLAS3::axpy (ctx1, ctx1.F.minusOne (), A3, A4);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx1, error, A4);

		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2>
bool testBLAS2ModulesConsistency (LELA::Context<Ring, Modules1> &ctx1, LELA::Context<Ring, Modules2> &ctx2, const char *text, size_t m, size_t n, size_t k) 
{
	std::ostringstream str;
	str << "Testing BLAS2 consistency of modules over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;

	typename LELA::Vector<Ring>::Dense v1(n), v2(m);
	typename LELA::Vector<Ring>::Sparse w1, w2;

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream11 (ctx1.F, n, m); 
	DenseMatrix<typename Ring::Element> M1 (stream11);

	RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream21 (ctx1.F, (double) k / (double) n, n, m);
	SparseMatrix<typename Ring::Element> M2 (stream21);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M3 (M2);

	LELA::RandomDenseStream<Ring, typename LELA::Vector<Ring>::Dense> stream1 (ctx1.F, n), stream2 (ctx1.F, m);
	LELA::RandomSparseStream<Ring, typename LELA::Vector<Ring>::Sparse> stream3 (ctx1.F, (double) k / (double) n, n), stream4 (ctx1.F, (double) k / (double) n, m);

	stream1 >>  v1;
	stream2 >>  v2;
	stream3 >>  w1;
	stream4 >>  w2;

	pass = testgemvConsistency (ctx1, ctx2, "sparse(row-wise)/dense/dense", M2, v1, v2, M1, v1, v2) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "sparse(col-wise)/dense/dense", M3, v2, v1, M3, v2, v1) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "sparse(row-wise)/sparse/sparse", M2, w1, w2, M2, w1, w2) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "sparse(col-wise)/sparse/sparse", M3, w2, w1, M3, w2, w1) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "sparse(row-wise)/sparse/dense", M2, w1, v2, M2, w1, v2) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "sparse(col-wise)/sparse/dense", M3, w2, v1, M3, w2, v1) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "dense/sparse/dense", M1, w1, v2, M1, w1, v2) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "dense/sparse/sparse", M1, w1, w2, M1, w1, w2) && pass; 
	pass = testgemvConsistency (ctx1, ctx2, "dense/dense/sparse", M1, v1, w2, M1, v1, w2) && pass;

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream41 (ctx1.F, n, n); 
	DenseMatrix<typename Ring::Element> M4 (stream41);

	RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream51 (ctx1.F, (double) k / (double) n, n, n);
	SparseMatrix<typename Ring::Element> M5 (stream51);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M6 (M5);

	pass = testtrmvConsistency (ctx1, ctx2, "sparse(row-wise)(LT, diag = 1)", M5, v1, M5, v1, LowerTriangular, true) && pass;
	pass = testtrmvConsistency (ctx1, ctx2, "sparse(row-wise)(UT, diag = 1)", M5, v1, M5, v1, UpperTriangular, true) && pass;
	pass = testtrmvConsistency (ctx1, ctx2, "sparse(row-wise)(LT, diag != 1)", M5, v1, M5, v1, LowerTriangular, false) && pass;
	pass = testtrmvConsistency (ctx1, ctx2, "sparse(row-wise)(UT, diag != 1)", M5, v1, M5, v1, UpperTriangular, false) && pass;

	pass = testtrmvConsistency (ctx1, ctx2, "sparse(col-wise)(LT, diag = 1)", M6, v1, M6, v1, LowerTriangular, true) && pass;
	pass = testtrmvConsistency (ctx1, ctx2, "sparse(col-wise)(UT, diag = 1)", M6, v1, M6, v1, UpperTriangular, true) && pass;
	pass = testtrmvConsistency (ctx1, ctx2, "sparse(col-wise)(LT, diag != 1)", M6, v1, M6, v1, LowerTriangular, false) && pass;
	pass = testtrmvConsistency (ctx1, ctx2, "sparse(col-wise)(UT, diag != 1)", M6, v1, M6, v1, UpperTriangular, false) && pass;

	SparseMatrix<typename Ring::Element> M5p (n, n);
	BLAS3::copy (ctx1, M5, M5p);
	makeNonsingDiag (ctx1.F, M5p, false);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M6p (M5p);

	pass = testtrsvConsistency (ctx1, ctx2, "sparse(row-wise)", M5, v1, M5, v1, LowerTriangular, true) && pass;
	pass = testtrsvConsistency (ctx1, ctx2, "sparse(row-wise)", M5, v1, M5, v1, UpperTriangular, true) && pass;
	pass = testtrsvConsistency (ctx1, ctx2, "sparse(row-wise)", M5p, v1, M5p, v1, LowerTriangular, false) && pass;
	pass = testtrsvConsistency (ctx1, ctx2, "sparse(row-wise))", M5p, v1, M5p, v1, UpperTriangular, false) && pass;

	pass = testtrsvConsistency (ctx1, ctx2, "sparse(col-wise)(LT, diag = 1)", M6, v1, M6, v1, LowerTriangular, true) && pass;
	pass = testtrsvConsistency (ctx1, ctx2, "sparse(col-wise)(UT, diag = 1)", M6, v1, M6, v1, UpperTriangular, true) && pass;
	pass = testtrsvConsistency (ctx1, ctx2, "sparse(col-wise)(LT, diag != 1)", M6p, v1, M6p, v1, LowerTriangular, false) && pass;
	pass = testtrsvConsistency (ctx1, ctx2, "sparse(col-wise)(UT, diag != 1)", M6p, v1, M6p, v1, UpperTriangular, false) && pass;

	pass = testgerConsistency (ctx1, ctx2, "sparse(row-wise) /dense  /dense ", M2, v2, v1, M2, v2, v1) && pass;
        pass = testgerConsistency (ctx1, ctx2, "sparse(col-wise) /dense  /dense ", M3, v1, v2, M3, v1, v2) && pass;
        pass = testgerConsistency (ctx1, ctx2, "sparse(row-wise) /sparse /sparse", M2, w2, w1, M2, w2, w1) && pass;
        pass = testgerConsistency (ctx1, ctx2, "sparse(col-wise) /sparse /sparse", M3, w1, w2, M3, w1, w2) && pass;
        pass = testgerConsistency (ctx1, ctx2, "sparse(row-wise) /sparse /dense ", M2, v2, w1, M2, v2, w1) && pass;
        pass = testgerConsistency (ctx1, ctx2, "sparse(col-wise) /sparse /dense ", M3, v1, w2, M3, v1, w2) && pass;
        pass = testgerConsistency (ctx1, ctx2, "dense            /sparse /dense ", M1, v2, w1, M1, v2, w1) && pass;
        pass = testgerConsistency (ctx1, ctx2, "dense            /sparse /sparse", M1, w2, w1, M1, w2, w1) && pass;
	pass = testgerConsistency (ctx1, ctx2, "dense            /dense  /sparse", M1, w2, v1, M1, w2, v1) && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules>
bool testBLAS2RepsConsistency (LELA::Context<Ring, Modules> &ctx, const char *text, size_t m, size_t n, size_t k) 
{
	std::ostringstream str;
	str << "Testing BLAS2 consistency over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;

	typename LELA::Vector<Ring>::Dense v1(n), v2(m);
	typename LELA::Vector<Ring>::Sparse w1, w2;

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream11 (ctx.F, n, m); 
	DenseMatrix<typename Ring::Element> M1 (stream11);

	RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream21 (ctx.F, (double) k / (double) n, n, m);
	SparseMatrix<typename Ring::Element> M2 (stream21);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M3 (M2);

	LELA::RandomDenseStream<Ring, typename LELA::Vector<Ring>::Dense> stream1 (ctx.F, n), stream2 (ctx.F, m);
	LELA::RandomSparseStream<Ring, typename LELA::Vector<Ring>::Sparse> stream3 (ctx.F, (double) k / (double) n, n), stream4 (ctx.F, (double) k / (double) n, m);

	stream1 >>  v1;
	stream2 >>  v2;
	stream3 >>  w1;
	stream4 >>  w2;
	
	pass = testgemvConsistency (ctx, ctx, "sparse(row-wise)	/dense	/dense 		with dense/dense/dense", M2, v1, v2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(col-wise)	/dense	/dense 		with dense/dense/dense", M3, v2, v1, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(row-wise)	/sparse	/sparse 	with dense/dense/dense", M2, w1, w2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(col-wise)	/sparse	/sparse 	with dense/dense/dense", M3, w2, w1, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(row-wise)	/sparse	/dense		with dense/dense/dense", M2, w1, v2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(col-wise)	/sparse	/dense		with dense/dense/dense", M3, w2, v1, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "dense		/sparse	/dense		with dense/dense/dense", M1, w1, v2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "dense		/sparse	/sparse		with dense/dense/dense", M1, w1, w2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "dense            /dense  /sparse         with dense/dense/dense", M1, v1, w2, M1, v1, v1) && pass;

	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream41 (ctx.F, n, n); 
	DenseMatrix<typename Ring::Element> M4 (stream41);

	RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream51 (ctx.F, (double) k / (double) n, n, n);
	SparseMatrix<typename Ring::Element> M5 (stream51);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M6 (M5);

	pass = testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag = 1)", M5, v1, M4, v1, LowerTriangular, true) && pass;
	pass = testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag = 1)", M5, v1, M4, v1, UpperTriangular, true) && pass;	
	pass = testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag != 1)", M5, v1, M4, v1, LowerTriangular, false) && pass;
	pass = testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag != 1)", M5, v1, M4, v1, UpperTriangular, false) && pass;

	pass = testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag = 1)", M6, v1, M4, v1, LowerTriangular, true) && pass;
	pass = testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag = 1)", M6, v1, M4, v1, UpperTriangular, true) && pass;	
	pass = testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag != 1)", M6, v1, M4, v1, LowerTriangular, false) && pass;
	pass = testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag != 1)", M6, v1, M4, v1, UpperTriangular, false) && pass;

	SparseMatrix<typename Ring::Element> M5p (n, n);
	BLAS3::copy (ctx, M5, M5p);
	makeNonsingDiag (ctx.F, M5p, false);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M6p (M5p);

	pass = testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag = 1)", M5, v1, M4, v1, LowerTriangular, true) && pass;
	pass = testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag = 1)", M5, v1, M4, v1, UpperTriangular, true) && pass;	
	pass = testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag != 1)", M5p, v1, M4, v1, LowerTriangular, false) && pass;
	pass = testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag != 1)", M5p, v1, M4, v1, UpperTriangular, false) && pass;

	pass = testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag = 1)", M6, v1, M4, v1, LowerTriangular, true) && pass;
	pass = testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag = 1)", M6, v1, M4, v1, UpperTriangular, true) && pass;	
	pass = testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag != 1)", M6p, v1, M4, v1, LowerTriangular, false) && pass;
	pass = testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag != 1)", M6p, v1, M4, v1, UpperTriangular, false) && pass;

	pass = testgerConsistency (ctx, ctx, "sparse(row-wise) /dense  /dense          with dense/dense/dense", M2, v2, v1, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "sparse(col-wise) /dense  /dense          with dense/dense/dense", M3, v1, v2, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "sparse(row-wise) /sparse /sparse         with dense/dense/dense", M2, w2, w1, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "sparse(col-wise) /sparse /sparse         with dense/dense/dense", M3, w1, w2, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "sparse(row-wise) /sparse /dense          with dense/dense/dense", M2, v2, w1, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "sparse(col-wise) /sparse /dense          with dense/dense/dense", M3, v1, w2, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "dense            /sparse /dense          with dense/dense/dense", M1, v2, w1, M1, v1, v1) && pass;
        pass = testgerConsistency (ctx, ctx, "dense            /sparse /sparse         with dense/dense/dense", M1, w2, w1, M1, v1, v1) && pass;
	pass = testgerConsistency (ctx, ctx, "dense            /dense  /sparse         with dense/dense/dense", M1, w2, v1, M1, v1, v1) && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

#endif // __LELA_TESTS_TEST_BLAS_LEVEL2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
