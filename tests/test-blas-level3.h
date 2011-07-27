/* tests/test-blas-level3.h
 * Copyright 2001, 2002, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic test suite for BLAS Level 3 routines
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_BLAS_LEVEL3_H
#define __LELA_TESTS_TEST_BLAS_LEVEL3_H

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
#include "lela/matrix/sparse.h"

using namespace std;
using namespace LELA;

typedef std::vector<std::pair<uint32, uint32> > Permutation;

/* Construct a random nonsingular upper triangular matrix from a random matrix */

template <class Field, class Matrix>
static Matrix &makeUpperTriangular (Field &F, Matrix &A, bool nonsingular)
{
	size_t i, j;

	typename Field::Element a;

	for (i = 0; i < A.rowdim (); ++i) {
		for (j = 0; j < std::min (i, A.coldim ()); ++j) {
			if (A.getEntry (a, i, j) && !F.isZero (a)) {
				A.setEntry (i, j, F.zero ());
				A.eraseEntry (i, j);
			}
		}
	}

	if (nonsingular) {
		NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

		for (i = 0; i < std::min (A.rowdim (), A.coldim ()); ++i) {
			if (!A.getEntry (a, i, i) || F.isZero (a)) {
				r.random (a);
				A.setEntry (i, i, a);
			}
		}
	}

	return A;
}

/* Construct a random (possibly nonsingular) lower triangular matrix from a random matrix */

template <class Field, class Matrix>
static Matrix &makeLowerTriangular (Field &F, Matrix &A, bool nonsingular)
{
	size_t i, j;

	typename Field::Element a;

	for (i = 0; i < A.rowdim (); ++i) {
		for (j = i + 1; j < A.coldim (); ++j) {
			if (A.getEntry (a, i, j) && !F.isZero (a)) {
				A.setEntry (i, j, F.zero ());
				A.eraseEntry (i, j);
			}
		}
	}

	if (nonsingular) {
		NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

		for (i = 0; i < std::min (A.rowdim (), A.coldim ()); ++i) {
			if (!A.getEntry (a, i, i) || F.isZero (a)) {
				r.random (a);
				A.setEntry (i, i, a);
			}
		}
	}

	return A;
}

/* Test 1: Copy and equal
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
static bool testCopyEqual (Context<Field, Modules> &ctx, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix copy, equal" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	report << "Input matrix M:" << endl;
	BLAS3::write (ctx, report, M);

	BLAS3::copy (ctx, M, M1);

	report << "Output matrix M1:" << endl;
	BLAS3::write (ctx, report, M1);

	if (!BLAS3::equal (ctx, M1, M)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M and M1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 2: scal, axpy, is_zero
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
static bool testScalAxpyIsZero (Context<Field, Modules> &ctx, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix scal, axpy, is_zero" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element a, nega;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	bool ret = true;

	DenseMatrix<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	
	reportUI << "Input matrix M:" << endl;
	BLAS3::write (ctx, reportUI, M);

	r.random (a);
	ctx.F.neg (nega, a);

	report << "Coefficient: ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::copy (ctx, M, M1);
	BLAS3::scal (ctx, a, M1);

	reportUI << "Output matrix a * M:" << endl;
	BLAS3::write (ctx, reportUI, M1);

	BLAS3::axpy (ctx, nega, M, M1);

	reportUI << "Output matrix -a * M + a * M:" << endl;
	BLAS3::write (ctx, reportUI, M1);

	if (!BLAS3::is_zero (ctx, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported matrix M1 is not zero" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 3: gemm with coefficients
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemmCoeff (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B)
{
	lela_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemm (coefficients)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> C (A.rowdim (), B.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), C);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B);

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));
	typename Field::Element a, b;

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	ctx.F.write (report, a) << ", b = ";
	ctx.F.write (report, b) << std::endl;

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), C);

	reportUI << "Output matrix C := a * A * B:" << endl;
	BLAS3::write (ctx, reportUI, C);

	typename Field::Element negainvb;

	ctx.F.inv (negainvb, a);
	ctx.F.mulin (negainvb, b);
	ctx.F.negin (negainvb);

	report << "Coefficient: -a^-1 b = ";
	ctx.F.write (report, negainvb) << std::endl;

	BLAS3::gemm (ctx, b, A, B, negainvb, C);

	reportUI << "Output matrix D := b * A * B - b * a^-1 * C:" << endl;
	BLAS3::write (ctx, reportUI, C);

	if (!BLAS3::is_zero (ctx, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported matrix D is not zero" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 4: gemm is associative
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
static bool testGemmAssoc (Context<Field, Modules> &ctx, const char *text, const Matrix &A, const Matrix &B, const Matrix &C)
{
	lela_check (A.coldim () == B.rowdim ());
	lela_check (B.coldim () == C.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemm (associativity)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> AB (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> BC (B.rowdim (), C.coldim ());
	DenseMatrix<typename Field::Element> ABpC (A.rowdim (), C.coldim ());
	DenseMatrix<typename Field::Element> ApBC (A.rowdim (), C.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), AB);
	BLAS3::scal (ctx, ctx.F.zero (), BC);
	BLAS3::scal (ctx, ctx.F.zero (), ABpC);
	BLAS3::scal (ctx, ctx.F.zero (), ApBC);

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B);

	report << "Input matrix C:" << endl;
	BLAS3::write (ctx, report, C);

	BLAS3::gemm (ctx, ctx.F.one (), A, B, ctx.F.zero (), AB);

	report << "Output matrix A * B:" << endl;
	BLAS3::write (ctx, report, AB);

	BLAS3::gemm (ctx, ctx.F.one (), AB, C, ctx.F.zero (), ABpC);

	report << "Output matrix (A * B) * C:" << endl;
	BLAS3::write (ctx, report, ABpC);

	BLAS3::gemm (ctx, ctx.F.one (), B, C, ctx.F.zero (), BC);

	report << "Output matrix B * C:" << endl;
	BLAS3::write (ctx, report, AB);

	BLAS3::gemm (ctx, ctx.F.one (), A, BC, ctx.F.zero (), ApBC);

	report << "Output matrix A * (B * C):" << endl;
	BLAS3::write (ctx, report, ApBC);

	if (!BLAS3::equal (ctx, ABpC, ApBC)) {
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported (A * B) * C != A * (B * C)" << endl;

		BLAS3::axpy (ctx, ctx.F.minusOne (), ABpC, ApBC);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx, error, ApBC);

		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* test 5: gemm with identity-matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
static bool testGemmIdent (Context<Field, Modules> &ctx, const char *text, const Matrix &A)
{
	ostringstream str;

	str << "Testing " << text << " gemm (identity)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> AI (A.rowdim (), A.coldim ());
	DenseMatrix<typename Field::Element> IA (A.rowdim (), A.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), IA);
	BLAS3::scal (ctx, ctx.F.zero (), AI);

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> I_l_stream (ctx.F, A.rowdim ()), I_r_stream (ctx.F, A.coldim ());

	DenseMatrix<typename Field::Element> I_l (I_l_stream);
	DenseMatrix<typename Field::Element> I_r (I_r_stream);

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A);

	report << "Identity matrix I_l:" << endl;
	BLAS3::write (ctx, report, I_l);

	report << "Identity matrix I_r:" << endl;
	BLAS3::write (ctx, report, I_r);

	BLAS3::gemm (ctx, ctx.F.one (), I_l, A, ctx.F.zero (), IA);

	report << "Output matrix I * A:" << endl;
	BLAS3::write (ctx, report, IA);

	if (!BLAS3::equal (ctx, A, IA)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported A != I * A" << endl;
		ret = false;
	}

	BLAS3::gemm (ctx, ctx.F.one (), A, I_r, ctx.F.zero (), AI);

	report << "Output matrix A * I:" << endl;
	BLAS3::write (ctx, report, AI);

	if (!BLAS3::equal (ctx, AI, A)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported A * I != A" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 7a: trmm and gemm (upper triangular)
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testTrmmGemmUpper (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B)
{
	ostringstream str;

	str << "Testing " << text << " trmm (upper triangular, consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;

	typename Matrix1::ContainerType A1 (A.rowdim (), A.coldim ());
	typename Matrix2::ContainerType B1 (B.rowdim (), B.coldim ());
	typename Matrix2::ContainerType C (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, A1);
	BLAS3::copy (ctx, B, B1);

	BLAS3::scal (ctx, ctx.F.zero (), C);

	makeUpperTriangular (ctx.F, A1, false);

	reportUI << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A1);

	reportUI << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B1);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trmm (ctx, a, A, B1, UpperTriangular, false);

	reportUI << "Output matrix a * A * B (trmm): " << std::endl;
	BLAS3::write (ctx, reportUI, B1);

	BLAS3::gemm (ctx, a, A1, B, ctx.F.zero (), C);

	reportUI << "Output matrix a * A * B (gemm): " << std::endl;
	BLAS3::write (ctx, reportUI, C);

	if (!BLAS3::equal (ctx, B1, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported results from trmm and gemm are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 7b: trmm and gemm (lower triangular)
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testTrmmGemmLower (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B)
{
	ostringstream str;

	str << "Testing " << text << " trmm (lower triangular, consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;

	typename Matrix1::ContainerType A1 (A.rowdim (), A.coldim ());
	typename Matrix2::ContainerType B1 (B.rowdim (), B.coldim ());
	typename Matrix2::ContainerType C (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, A1);
	BLAS3::copy (ctx, B, B1);

	BLAS3::scal (ctx, ctx.F.zero (), C);

	makeLowerTriangular (ctx.F, A1, false);

	reportUI << "Input matrix A:" << endl;
	BLAS3::write (ctx, reportUI, A1);

	reportUI << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B1);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trmm (ctx, a, A, B1, LowerTriangular, false);

	reportUI << "Output matrix a * A * B (trmm): " << std::endl;
	BLAS3::write (ctx, reportUI, B1);

	BLAS3::gemm (ctx, a, A1, B, ctx.F.zero (), C);

	reportUI << "Output matrix a * A * B (gemm): " << std::endl;
	BLAS3::write (ctx, reportUI, C);

	if (!BLAS3::equal (ctx, B1, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported results from trmm and gemm are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 8: Row-echelon and gemm
 *
 * Return true on success and false on failure
 */

/* Elimination-step */

template <class Field, class Modules, class Matrix1, class Matrix2>
void eliminate (Context<Field, Modules> &ctx, Matrix1 &M, Matrix2 &pivotRow,
		size_t row, size_t col, size_t rowdim, size_t coldim) 
{
	typename Matrix1::ConstSubmatrixType pivotCol (M, row, col, rowdim, 1);
	typename Matrix1::SubmatrixType block (M, row, col, rowdim, coldim);

	BLAS2::ger (ctx, ctx.F.minusOne (), *pivotCol.colBegin (), *pivotRow.rowBegin (), block);
}

/* Dumb elimination code
 *
 * This computes the reduced row-echelon form of a matrix, placing the
 * transform in U and the row-echelon form in R. The output should
 * satisfy R = UA. Places rank in parameter rank.
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
void rowEchelon (Context<Field, Modules> &ctx, Matrix1 &U, Matrix1 &R, const Matrix2 &A, size_t &rank) 
{
	lela_check (U.rowdim () == A.rowdim ());
	lela_check (U.coldim () == A.rowdim ());
	lela_check (R.rowdim () == A.rowdim ());
	lela_check (R.coldim () == A.coldim ());

	typename Matrix1::ContainerType M (U.rowdim (), U.coldim () + R.coldim ());
	typename Matrix1::ContainerType::SubmatrixType M1 (M, 0, 0, R.rowdim (), R.coldim ());
	typename Matrix1::ContainerType::SubmatrixType M2 (M, 0, R.coldim (), U.rowdim (), U.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), M);

	StandardBasisStream<Field, typename Matrix1::ContainerType::SubmatrixType::Row> stream (ctx.F, U.coldim ());
	typename Matrix1::ContainerType::SubmatrixType::RowIterator ip = M2.rowBegin ();

	for (; ip != M2.rowEnd (); ++ip)
		stream >> *ip;

	BLAS3::copy (ctx, A, M1);

	unsigned int idx;
	typename Field::Element Mjj_inv, a;

	rank = 0;

	for (idx = 0; idx < M.rowdim (); ++idx) {
		if (!M.getEntry (a, idx, idx) || ctx.F.isZero (a)) {
			typename Matrix1::ContainerType::ColIterator col;
			typename Matrix1::ContainerType::Col::iterator i;
			unsigned int c_idx = idx + 1;

			col = M.colBegin () + idx;
			i = col->begin () + idx + 1;

			while (ctx.F.isZero (*i) && i != col->end ())
				++i, ++c_idx;

			if (i == col->end ())
				continue;
			else {
				typename Matrix1::ContainerType::RowIterator row1, row2;

				row1 = M.rowBegin () + idx;
				row2 = M.rowBegin () + c_idx;

				std::swap_ranges (row1->begin () + idx, row1->end (), row2->begin () + idx);
			}
		}

		M.getEntry (Mjj_inv, idx, idx);
		ctx.F.invin (Mjj_inv);
		typename Matrix1::ContainerType::SubmatrixType realPivotRow (M, idx, idx, 1, M.coldim () - idx);
		BLAS3::scal (ctx, Mjj_inv, realPivotRow);

		if (idx > 0)
			eliminate (ctx, M, realPivotRow, 0, idx, idx, M.coldim () - idx);

		if (idx < M.rowdim () - 1)
			eliminate (ctx, M, realPivotRow, idx + 1, idx,
				   M.rowdim () - idx - 1, M.coldim () - idx);

		++rank;
	}

	BLAS3::copy (ctx, M1, R);
	BLAS3::copy (ctx, M2, U);
}

template <class Field, class Modules, class Matrix>
static bool testGemmRowEchelon (Context<Field, Modules> &ctx, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " gemm (row-echelon)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	size_t rank;

	DenseMatrix<typename Field::Element> R (M.rowdim (), M.coldim ());
	DenseMatrix<typename Field::Element> U (M.rowdim (), M.rowdim ());
	DenseMatrix<typename Field::Element> UM (U.rowdim (), M.coldim ());

	BLAS3::scal (ctx, ctx.F.zero (), UM);

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> stream (ctx.F, M.rowdim ());

	DenseMatrix<typename Field::Element> I (UM.rowdim (), UM.coldim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	BLAS3::write (ctx, report, M);

	rowEchelon (ctx, U, R, M, rank);

	report << "Computed transform U:" << std::endl;
	BLAS3::write (ctx, report, U);

	report << "Computed row-echelon form R:" << std::endl;
	BLAS3::write (ctx, report, R);

	report << "Computed rank = " << rank << std::endl;

	BLAS3::gemm (ctx, ctx.F.one (), U, M, ctx.F.zero (), UM);

	report << "Computed product UM:" << endl;
	BLAS3::write (ctx, report, UM);

	if (!BLAS3::equal (ctx, UM, R)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported matrix UM != R" << endl;
		ret = false;
	}

	if (rank == M.rowdim () && rank == M.coldim ()) {
		report << "Full rank reported and matrix is square. Checking whether R is identity-matrix..." << std::endl;

		if (!BLAS3::equal (ctx, UM, I)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: BLAS3 reported matrix R is not identity" << endl;
			ret = false;
		}

		report << "Checking whether U is inverse..." << std::endl;

		BLAS3::gemm (ctx, ctx.F.one (), M, U, ctx.F.zero (), UM);

		report << "Computed product MU:" << endl;
		BLAS3::write (ctx, report, UM);

		if (!BLAS3::equal (ctx, UM, R)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: BLAS3 reported matrix MU is not identity" << endl;
			ret = false;
		}
	} else {
		report << "Reportedly M is not invertible, so remaining tests not meaningful" << std::endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 11a: trsm (lower triangular)
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
static bool testTrsmLower (Context<Field, Modules> &ctx, const char *text, const Matrix &A, const Matrix &B)
{
	ostringstream str;

	str << "Testing " << text << " trsm (lower triangular)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	typename Matrix::ContainerType U (A.rowdim (), A.coldim ());
	typename Matrix::ContainerType B1 (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, U);
	makeNonsingDiag (ctx.F, U, false);

	BLAS3::copy (ctx, B, B1);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix U:" << std::endl;
	BLAS3::write (ctx, reportUI, U);

	reportUI << "Input matrix B: " << std::endl;
	BLAS3::write (ctx, reportUI, B);

	typename Field::Element a, ainv;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);
	ctx.F.inv (ainv, a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trsm (ctx, a, U, B1, LowerTriangular, false);

	reportUI << "Output matrix U^-1 * B: " << std::endl;
	BLAS3::write (ctx, reportUI, B1);

	BLAS3::trmm (ctx, ainv, U, B1, LowerTriangular, false);

	reportUI << "Output matrix U U^-1 * B: " << std::endl;
	BLAS3::write (ctx, reportUI, B1);

	if (!BLAS3::equal (ctx, B, B1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: VectorxDomain reported B != U * U^-1 * B" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 11a: trsm (upper triangular)
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
static bool testTrsmUpper (Context<Field, Modules> &ctx, const char *text, const Matrix &A, const Matrix &B)
{
	ostringstream str;

	str << "Testing " << text << " trsm (upper triangular)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	typename Matrix::ContainerType U (A.rowdim (), A.coldim ());
	typename Matrix::ContainerType B1 (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, U);
	makeNonsingDiag (ctx.F, U, false);

	BLAS3::copy (ctx, B, B1);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix U:" << std::endl;
	BLAS3::write (ctx, reportUI, U);

	reportUI << "Input matrix B: " << std::endl;
	BLAS3::write (ctx, reportUI, B);

	typename Field::Element a, ainv;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);
	ctx.F.inv (ainv, a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trsm (ctx, a, U, B1, UpperTriangular, false);

	reportUI << "Output matrix U^-1 * B: " << std::endl;
	BLAS3::write (ctx, reportUI, B1);

	BLAS3::trmm (ctx, ainv, U, B1, UpperTriangular, false);

	reportUI << "Output matrix U U^-1 * B: " << std::endl;
	BLAS3::write (ctx, reportUI, B1);

	if (!BLAS3::equal (ctx, B, B1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: VectorxDomain reported B != U * U^-1 * B" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 13: trsm with coefficient
 *
 * B must be iterable by columns
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testTrsmCoeff (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B) 
{
	lela_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " trsm (coefficient)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	typename Matrix1::ContainerType U (A.rowdim (), A.coldim ());

	BLAS3::copy (ctx, A, U);
	makeUpperTriangular (ctx.F, U, true);

	DenseMatrix<typename Field::Element> aUinvB (U.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> UinvB (U.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
        ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	reportUI << "Input matrix U:" << endl;
	BLAS3::write (ctx, reportUI, U);

	reportUI << "Input matrix B:" << endl;
	BLAS3::write (ctx, reportUI, B);

	BLAS3::copy (ctx, B, aUinvB);
	BLAS3::copy (ctx, B, UinvB);

	typename Field::RandIter r (ctx.F);
	typename Field::Element a;

	r.random (a);

	report << "Coefficient: a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trsm (ctx, a, U, aUinvB, UpperTriangular, false);

	reportUI << "Output matrix a U^-1 * B:" << endl;
	BLAS3::write (ctx, reportUI, aUinvB);

	BLAS3::trsm (ctx, ctx.F.one (), U, UinvB, UpperTriangular, false);

	reportUI << "Output matrix U^-1 * B:" << endl;
	BLAS3::write (ctx, reportUI, UinvB);

	BLAS3::scal (ctx, a, UinvB);

	reportUI << "Output matrix a (U^-1 * B):" << endl;
	BLAS3::write (ctx, reportUI, UinvB);

	if (!BLAS3::equal (ctx, UinvB, aUinvB)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: BLAS3 reported a U^-1 * B != a (U^-1 * B)" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 14: permuteRows, permuteColumns
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
bool testPermutation (Context<Field, Modules> &ctx, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " permutations" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	MersenneTwister MT;

	Permutation P, Pinv;

	// Create a random row permutation
	for (unsigned int i = 0; i < M.rowdim (); ++i) {
		unsigned int row1, row2;

		do {
			row1 = MT.randomInt () % M.rowdim ();
			row2 = MT.randomInt () % M.rowdim ();
		} while (row1 == row2);

		P.push_back (Permutation::value_type (row1, row2));
	}

	// Construct the inverse of this transposition
	Pinv.resize (P.size ());
	std::copy (P.begin (), P.end (), Pinv.begin ());
	std::reverse (Pinv.begin (), Pinv.end ());

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	BLAS3::write (ctx, report, M);

	report << "Permutation P:    ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << endl;
	report << "Permutation P^-1: ";
	BLAS1::write_permutation (report, Pinv.begin (), Pinv.end ()) << endl;

	// Apply the permutation and then its inverse to a copy of M
	Matrix M1 (M);

	BLAS3::permute_rows (ctx, P.begin (), P.end (), M1);
	BLAS3::permute_rows (ctx, Pinv.begin (), Pinv.end (), M1);

	report << "Output matrix P^-1 PM:" << endl;
	BLAS3::write (ctx, report, M1);

	// Compare M and M1
	if (!BLAS3::equal (ctx, M, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: M != P^-1 PM" << endl;
		ret = false;
	}

	// Now we do exactly the same thing with columns
	BLAS3::copy (ctx, M, M1);
	P.clear ();
	Pinv.clear ();

	// Create a random column permutation
	for (unsigned int i = 0; i < M.coldim (); ++i) {
		unsigned int col1, col2;

		do {
			col1 = MT.randomInt () % M.coldim ();
			col2 = MT.randomInt () % M.coldim ();
		} while (col1 == col2);

		P.push_back (Permutation::value_type (col1, col2));
	}

	// Construct the inverse of this transposition
	Pinv.resize (P.size ());
	std::copy (P.begin (), P.end (), Pinv.begin ());
	std::reverse (Pinv.begin (), Pinv.end ());

	report << "Permutation P:    ";
	BLAS1::write_permutation (report, P.begin (), P.end ()) << endl;
	report << "Permutation P^-1: ";
	BLAS1::write_permutation (report, Pinv.begin (), Pinv.end ()) << endl;

	// Apply the permutation and then its inverse to a copy of M
	BLAS3::permute_cols (ctx, P.begin (), P.end (), M1);
	BLAS3::permute_cols (ctx, Pinv.begin (), Pinv.end (), M1);

	report << "Output matrix MPP^-1:" << endl;
	BLAS3::write (ctx, report, M1);

	// Compare M and M1
	if (!BLAS3::equal (ctx, M, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: M != MPP^-1" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Ring, class Modules, class Matrix>
bool testpermute_rows (LELA::Context<Ring, Modules> &ctx,
		       const char *text,
		       const Matrix &A1)
{

        ostringstream str;
        str << "Testing " << text << " permute_rows (does it do anything?)" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

	if (BLAS3::is_zero (ctx, A1)) {
		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_WARNING)
			<< "Input-matrix is 0. Skipping test." << std::endl;
		commentator.stop ("skipped");
		return true;
	}

	typename Matrix::ContainerType A2 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx, A1, A2);

	MersenneTwister MT;

	Permutation P;

	// Create a random row permutation
	for (unsigned int i = 0; i < A1.rowdim (); ++i) {
		unsigned int row1, row2;

		do {
			row1 = MT.randomInt () % A1.rowdim ();
			row2 = MT.randomInt () % A1.rowdim ();
		} while (row1 == row2);

		P.push_back (Permutation::value_type (row1, row2));
	}

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx, report, A1);

	BLAS3::permute_rows (ctx, P.begin (), P.end (), A2);

	report << "Matrix A_1 P: " << std::endl;
	BLAS3::write(ctx, report, A2);

        if (BLAS3::equal (ctx, A1, A2))
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A_1 =  A_1 P " << std::endl;
			pass = false;
		}

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring, class Modules, class Matrix>
bool testpermute_cols (LELA::Context<Ring, Modules> &ctx,
                       const char *text,
                       const Matrix &A1)
{

        ostringstream str;
        str << "Testing " << text << " permute_cols (does it do anything?)" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

        typename Matrix::ContainerType A2 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx, A1, A2);

        MersenneTwister MT;

        Permutation P;

        // Create a random row permutation
        for (unsigned int i = 0; i < A1.coldim (); ++i) {
                unsigned int col1,col2;

                do {
                        col1 = MT.randomInt () % A1.coldim ();
                        col2 = MT.randomInt () % A1.coldim ();
                } while (col1 == col2);

                P.push_back (Permutation::value_type (col1, col2));
        }

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx, report, A1);

	BLAS3::permute_cols (ctx, P.begin (), P.end (), A2);

        report << "Matrix A_1 P: " << std::endl;
	BLAS3::write(ctx, report, A2);

        if (BLAS3::equal (ctx, A1, A2))
                {
                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                << "ERROR: A_1 =  A_1 P " << std::endl;
                        pass = false;
                }

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

/* Test 15: read and write
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Matrix>
bool testReadWriteFormat (Context<Field, Modules> &ctx, const char *text, const Matrix &M, FileFormatTag format)
{
	static const char *format_names[] = 
		{ "detect", "unknown", "Turner", "one-based", "Dumas", "Maple", "Matlab", "Sage", "pretty", "PNG" };

	ostringstream str;
	str << "Testing " << text << " read/write (format " << format_names[format] << ")" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool pass = true;

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Matrix M1 (M.rowdim (), M.coldim ());

	ostringstream output;

	try {
		BLAS3::write (ctx, output, M, format);
	}
	catch (LELAError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Caught LELAError: " << e << endl;
		pass = false;
	}

#ifdef __LELA_HAVE_LIBPNG
	if (format != FORMAT_PNG)
#endif // __LELA_HAVE_LIBPNG
	report << "Matrix-output as string:" << std::endl << output.str ();

	istringstream input (output.str ());

	try {
		BLAS3::read (ctx, input, M1, format);
		report << "Matrix as read from " << format_names[format] << " format" << std::endl;
		BLAS3::write (ctx, report, M1);

		if (!BLAS3::equal (ctx, M, M1)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Matrix as read does not equal original" << endl;
			pass = false;
		}
	}
	catch (InvalidMatrixInput) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Caught InvalidMatrixInput while trying to read the matrix" << endl;
		pass = false;
	}
	catch (LELAError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Caught LELAError: " << e << endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testReadWrite (Context<Field, Modules> &ctx, const char *text, const Matrix &M)
{
	ostringstream str;
	str << "Testing " << text << " read/write" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool pass = true;

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << std::endl;
	BLAS3::write (ctx, report, M);

	FileFormatTag formats[] = { FORMAT_TURNER, FORMAT_DUMAS, FORMAT_MATLAB, FORMAT_PRETTY };

	for (size_t i = 0; i < sizeof (formats) / sizeof (FileFormatTag); ++i)
		pass = testReadWriteFormat (ctx, text, M, formats[i]) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Modules, class Matrix>
bool testReadWrite (Context<GF2, Modules> &ctx, const char *text, const Matrix &M)
{
	ostringstream str;
	str << "Testing " << text << " read/write" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool pass = true;

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << std::endl;
	BLAS3::write (ctx, report, M);

#ifdef __LELA_HAVE_LIBPNG
	FileFormatTag formats[] = { FORMAT_TURNER, FORMAT_DUMAS, FORMAT_MATLAB, FORMAT_PRETTY, FORMAT_PNG };
#else // !__LELA_HAVE_LIBPNG
	FileFormatTag formats[] = { FORMAT_TURNER, FORMAT_DUMAS, FORMAT_MATLAB, FORMAT_PRETTY };
#endif // __LELA_HAVE_LIBPNG

	for (size_t i = 0; i < sizeof (formats) / sizeof (FileFormatTag); ++i)
		pass = testReadWriteFormat (ctx, text, M, formats[i]) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2>
bool testscalConsistency  (LELA::Context<Ring, Modules1> &ctx1,
			   LELA::Context<Ring, Modules2> &ctx2,
			   const char *text,
			   const Matrix1 &A1, 
                           const Matrix2 &A2)
{

	ostringstream str;
	str << "Testing " << text << " scal consistency" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool pass = true;

	typename Matrix1::ContainerType A3 (A1.rowdim (), A1.coldim ());
	typename Matrix2::ContainerType A4 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1, A1, A3);
	BLAS3::copy(ctx1, A1, A4);

	typename Ring::Element a;
	NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

	r.random (a);

	reportUI << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A3);

	BLAS3::scal (ctx1, a, A3);

	report << "Coefficient a: "<< std::endl;
	ctx1.F.write(report, a);

	reportUI << "Matrix a A_1: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A3);

	reportUI << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx2, reportUI, A4);

	BLAS3::scal (ctx2, a, A4);

	reportUI << "Matrix a A_2: " << std::endl;;
	BLAS3::write (ctx2, reportUI, A4);

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
bool testaxpyConsistency  (LELA::Context<Ring, Modules1> &ctx1,
                           LELA::Context<Ring, Modules2> &ctx2,
                           const char *text,
                           const Matrix1 &A1, const Matrix2 &A2,
                           const Matrix3 &A3, const Matrix4 &A4)
{

        ostringstream str;
        str << "Testing " << text << " axpy consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

	typename Matrix3::ContainerType A6 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1, A1, A6);

        typename Matrix2::ContainerType A7 (A2.rowdim (), A2.coldim ());
        typename Matrix4::ContainerType A8 (A2.rowdim (), A2.coldim ());

	BLAS3::copy(ctx1, A2, A7);
	BLAS3::copy(ctx1, A2, A8);

        typename Ring::Element a;
        NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

        r.random (a);

        reportUI << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A1);

        reportUI << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A7);

	BLAS3::axpy (ctx1, a, A1, A7);

	report << "Coefficient a: "<< std::endl;
	ctx1.F.write(report, a);

        reportUI << "Matrix a A_1 + A_2: " << std::endl;;
	BLAS3::write (ctx1, reportUI, A7);

        reportUI << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A6);

        reportUI << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A8);

	BLAS3::axpy (ctx2, a, A6, A8);

        reportUI << "Matrix a A_3 + A_4: " << std::endl;
	BLAS3::write (ctx2, reportUI, A8);

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
bool testgemmConsistency  (LELA::Context<Ring, Modules1> &ctx1,
                           LELA::Context<Ring, Modules2> &ctx2,
                           const char *text,
                           const Matrix1 &A1, const Matrix2 &A2,
                           const Matrix3 &A3, const Matrix4 &A4,
			   const Matrix5 &A5, const Matrix6 &A6)
{

        ostringstream str;
        str << "Testing " << text << " gemm consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

        typename Matrix4::ContainerType A7 (A1.rowdim (), A1.coldim ());
	typename Matrix5::ContainerType A8 (A2.rowdim (), A2.coldim ());

	BLAS3::copy(ctx1, A1, A7);
	BLAS3::copy(ctx1, A2, A8);

        typename Matrix3::ContainerType A9 (A3.rowdim (), A3.coldim ());
	typename Matrix6::ContainerType A10 (A3.rowdim (), A3.coldim ());

	BLAS3::copy(ctx1, A3, A9);
	BLAS3::copy(ctx1, A3, A10);

        typename Ring::Element a, b;
        NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

        r.random (a);
	r.random (b);

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

        report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A2);

	report << "Matrix A_3: "<< std::endl;
	BLAS3::write (ctx1, report, A9);

	BLAS3::gemm (ctx1, a, A1, A2, b, A9);

	report << "Coefficient a: ";
	ctx1.F.write(report, a) << std::endl;

        report << "Coefficient b: ";
        ctx1.F.write(report, b) << std::endl;

        reportUI << "Matrix a A_1 A_2 + c A_3: " << std::endl;
	BLAS3::write (ctx1, reportUI, A9);

        reportUI << "Matrix A_4: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A7);

        reportUI << "Matrix A_5: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A8);

        reportUI << "Matrix A_6: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A10);

	BLAS3::gemm (ctx2, a, A7, A8, b, A10);

        reportUI << "Matrix a A_4 A_5 + b A_6: " << std::endl;
	BLAS3::write (ctx2, reportUI, A10);

        if (!BLAS3::equal (ctx1, A9, A10))
	{
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a A_1 A_2 + b A_3  !=  a A_4 A_5 + b A_6 " << std::endl;

		BLAS3::axpy (ctx1, ctx1.F.minusOne (), A9, A10);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx1, error, A10);

		pass = false;
	}

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2,  class Matrix3, class Matrix4>
bool testtrmmConsistency  (LELA::Context<Ring, Modules1> &ctx1,
                           LELA::Context<Ring, Modules2> &ctx2,
                           const char *text,
                           const Matrix1 &A1, const Matrix2 &A2,
                           const Matrix3 &A3, const Matrix4 &A4,
                           TriangularMatrixType type, bool diagIsOne)
{

        ostringstream str;
        str << "Testing " << text << " trmm consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

        typename Matrix3::ContainerType A5 (A1.rowdim (), A1.coldim ());
        typename Matrix2::ContainerType A6 (A2.rowdim (), A1.coldim ());
	typename Matrix4::ContainerType A7 (A2.rowdim (), A2.coldim ()); 

	BLAS3::copy(ctx1, A1, A5);
	BLAS3::copy(ctx1, A2, A6);
	BLAS3::copy(ctx1, A2, A7);

        typename Ring::Element a;
        NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

        r.random (a);

        reportUI << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A1);

        reportUI << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A6);

	BLAS3::trmm (ctx1, a, A1, A6, type, diagIsOne);

	report << "Coefficient a: ";
	ctx1.F.write (report, a) << std::endl;

        reportUI << "Matrix a A_1 A_2: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A6);

        reportUI << "Matrix A_3: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A5);

        reportUI << "Matrix A_4: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A7);

	BLAS3::trmm (ctx2, a, A5, A7, type, diagIsOne);

        reportUI << "Matrix a A_3 A_4: "<< std::endl;
	BLAS3::write (ctx2, reportUI, A7);

        if (!BLAS3::equal (ctx1, A6, A7))
        {
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a A_1 A_2  !=  a A_3 A_4 " << std::endl;

		BLAS3::axpy (ctx1, ctx1.F.minusOne (), A6, A7);

		error << "Difference is" << std::endl;
		BLAS3::write (ctx1, error, A7);

                pass = false;
        }

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2,  class Matrix3, class Matrix4>
bool testtrsmConsistency  (LELA::Context<Ring, Modules1> &ctx1,
                           LELA::Context<Ring, Modules2> &ctx2,
                           const char *text,
                           const Matrix1 &A1, const Matrix2 &A2,
                           const Matrix3 &A3, const Matrix4 &A4,
                           TriangularMatrixType type, bool diagIsOne)
{

        ostringstream str;
        str << "Testing " << text << " trsm consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

        typename Matrix3::ContainerType A5 (A1.rowdim (), A1.coldim ());
        typename Matrix2::ContainerType A6 (A2.rowdim (), A2.coldim ());
        typename Matrix4::ContainerType A7 (A2.rowdim (), A2.coldim ());

	BLAS3::copy(ctx1, A1, A5);
	BLAS3::copy(ctx1, A2, A6);
	BLAS3::copy(ctx1, A2, A7);

        typename Ring::Element a;
        NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

        r.random (a);

        reportUI << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A1);

        reportUI << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A6);

	BLAS3::trsm (ctx1, a, A1, A6, type, diagIsOne);

	report <<"Coefficient a: ";
	ctx1.F.write(report, a) << std::endl;

        reportUI << "Matrix a ((A_1)^-1) A_2: "<<std::endl;
	BLAS3::write (ctx1, reportUI, A6);

        reportUI << "Matrix A_3: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A5);

        reportUI << "Matrix A_4: "<< std::endl;
	BLAS3::write (ctx1, reportUI, A7);

	BLAS3::trsm (ctx2, a, A5, A7, type, diagIsOne);

        reportUI << "Matrix a ((A_3)^-1) A_4: "<< std::endl; 
	BLAS3::write (ctx2, reportUI, A7);

        if (!BLAS3::equal (ctx1, A6, A7))
        {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a ((A_1)^-1) A_2  !=  a ((A_3)^-1) A_4 " << std::endl;
                pass = false;
        }

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2>
bool testpermute_rowsConsistency  (LELA::Context<Ring, Modules1> &ctx1,
                           LELA::Context<Ring, Modules2> &ctx2,
                           const char *text,
                           const Matrix1 &A1, const Matrix2 &A2)
{

        ostringstream str;
        str << "Testing " << text << " permute_rows consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

	typename Matrix1::ContainerType A3 (A1.rowdim (), A1.coldim ());
        typename Matrix2::ContainerType A4 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1, A1, A3);
	BLAS3::copy(ctx1, A1, A4);

	MersenneTwister MT;

	Permutation P;

	// Create a random row permutation
	for (unsigned int i = 0; i < A1.rowdim (); ++i) {
		unsigned int row1, row2;

		do {
			row1 = MT.randomInt () % A1.rowdim ();
			row2 = MT.randomInt () % A1.rowdim ();
		} while (row1 == row2);

		P.push_back (Permutation::value_type (row1, row2));
	}

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

	BLAS3::permute_rows (ctx1, P.begin (), P.end (), A3);

	report << "Matrix A_1 P: " << std::endl;
	BLAS3::write(ctx1, report, A3);

	report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A4);

	BLAS3::permute_rows (ctx1, P.begin (), P.end (), A4);

        report << "A_2 P: "<< std::endl;
	BLAS3::write (ctx2, report, A4);

        if (!BLAS3::equal (ctx1, A3, A4))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A_1 P  !=  A_2 P " << std::endl;
		pass = false;
	}

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring, class Modules1, class Modules2, class Matrix1, class Matrix2>
bool testpermute_colsConsistency  (LELA::Context<Ring, Modules1> &ctx1,
				   LELA::Context<Ring, Modules2> &ctx2,
				   const char *text,
				   const Matrix1 &A1, const Matrix2 &A2)
{

        ostringstream str;
        str << "Testing " << text << " permute_cols consistency" << std::ends;
        commentator.start (str.str ().c_str ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

        bool pass = true;

	typename Matrix1::ContainerType A3 (A1.rowdim (), A1.coldim ());
        typename Matrix2::ContainerType A4 (A1.rowdim (), A1.coldim ());

	BLAS3::copy(ctx1, A1, A3);
	BLAS3::copy(ctx1, A1, A4);

	MersenneTwister MT;

	Permutation P;

	// Create a random column permutation
	for (unsigned int i = 0; i < A1.coldim (); ++i) {
		unsigned int col1, col2;

		do {
			col1 = MT.randomInt () % A1.coldim ();
			col2 = MT.randomInt () % A1.coldim ();
		} while (col1 == col2);

		P.push_back (Permutation::value_type (col1, col2));
	}

        report << "Matrix A_1: "<< std::endl;
	BLAS3::write (ctx1, report, A1);

	BLAS3::permute_cols (ctx1, P.begin (), P.end (), A3);

	report << "Matrix A_1 P: " << std::endl;
	BLAS3::write(ctx1, report, A3);

	report << "Matrix A_2: "<< std::endl;
	BLAS3::write (ctx1, report, A4);

	BLAS3::permute_cols (ctx1, P.begin (), P.end (), A4);

        report << "A_2 P: "<< std::endl;
	BLAS3::write (ctx2, report, A4);

        if (!BLAS3::equal (ctx1, A3, A4))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A_1 P  !=  A_2 P " << std::endl;
		pass = false;
	}

        commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3 (Context<Field, Modules> &ctx, const char *text,
		Matrix &M1, Matrix &M2, Matrix &M3, Matrix &M4,
		MatrixIteratorTypes::RowCol)
{
	ostringstream str;
	str << "Testing BLAS3 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;

	if (!testTrmmGemmUpper (ctx, text, M4, M2)) pass = false;
	if (!testTrmmGemmLower (ctx, text, M4, M2)) pass = false;
	if (!testTrsmLower (ctx, text, M4, M2)) pass = false;
	if (!testTrsmUpper (ctx, text, M4, M2)) pass = false;
	if (!testTrsmCoeff (ctx, text, M4, M2)) pass = false;

	if (!testPermutation (ctx, text, M1)) pass = false;
	if (!testpermute_rows(ctx, text, M1)) pass = false;
	if (!testpermute_cols(ctx, text, M1)) pass = false;
	if (!testReadWrite (ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3 (Context<Field, Modules> &ctx, const char *text,
		Matrix &M1, Matrix &M2, Matrix &M3, Matrix &M4,
		MatrixIteratorTypes::Row) 
{
	ostringstream str;
	str << "Testing BLAS3 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;

	if (!testTrmmGemmUpper (ctx, text, M4, M2)) pass = false;
	if (!testTrmmGemmLower (ctx, text, M4, M2)) pass = false;
	if (!testTrsmLower (ctx, text, M4, M2)) pass = false;
	if (!testTrsmUpper (ctx, text, M4, M2)) pass = false;
	if (!testTrsmCoeff (ctx, text, M4, M2)) pass = false;

//	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;  // Needs ColIterator

	if (!testPermutation (ctx, text, M1)) pass = false;
        if (!testpermute_rows(ctx, text, M1)) pass = false;
        if (!testpermute_cols(ctx, text, M1)) pass = false;
	if (!testReadWrite (ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3 (Context<Field, Modules> &ctx, const char *text,
		Matrix &M1, Matrix &M2, Matrix &M3, Matrix &M4,
		MatrixIteratorTypes::Col) 
{
	ostringstream str;
	str << "Testing BLAS3 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;

	if (!testTrmmGemmUpper (ctx, text, M4, M2)) pass = false;
	if (!testTrmmGemmLower (ctx, text, M4, M2)) pass = false;
	if (!testTrsmLower (ctx, text, M4, M2)) pass = false;
	if (!testTrsmUpper (ctx, text, M4, M2)) pass = false;

        if (!testpermute_rows(ctx, text, M1)) pass = false;
        if (!testpermute_cols(ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3Submatrix (Context<Field, Modules> &ctx, const char *text, // 
			 const Matrix &M1, const Matrix &M2, const Matrix &M3, const Matrix &M4,
			 MatrixIteratorTypes::Row) 
{
	ostringstream str;
	str << "Testing BLAS3 with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
//	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;  // Needs ColIterator

	if (!testTrsmCoeff (ctx, text, M4, M2)) pass = false;

	if (!testPermutation (ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2>
bool testBLAS3ModulesConsistency (LELA::Context<Ring, Modules1> &ctx1, LELA::Context<Ring, Modules2> &ctx2, const char *text, size_t m, size_t n, size_t p, size_t k) 
{
	std::ostringstream str;
	str << "Testing BLAS3 consistency of modules over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;
        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream11 (ctx1.F, n, m);
        DenseMatrix<typename Ring::Element> M1 (stream11);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream21 (ctx1.F, (double) k / (double) n, n, m);
        SparseMatrix<typename Ring::Element> M2 (stream21);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > M3 (M2);

	pass = testscalConsistency (ctx1, ctx2, "dense", M1, M1) && pass;
        pass = testscalConsistency (ctx1, ctx2, "sparse(row-wise)", M2, M2) && pass;
        pass = testscalConsistency (ctx1, ctx2, "sparse(col-wise)", M3, M3) && pass;

	pass = testaxpyConsistency (ctx1, ctx2, "dense /dense", M1, M1, M1, M1) && pass;
        pass = testaxpyConsistency (ctx1, ctx2, "sparse(row-wise)/dense", M2, M1, M2, M1) &&pass;
        pass = testaxpyConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)", M2, M2, M2, M2) &&pass;

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream12 (ctx1.F, n, m);
        DenseMatrix<typename Ring::Element> A1_dense (stream12);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream13 (ctx1.F, (double) k / (double) n, n, m);
        SparseMatrix<typename Ring::Element> A1_sparse (stream13);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A1_trans (A1_sparse);

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream14 (ctx1.F, p, n);
        DenseMatrix<typename Ring::Element> A2_dense (stream14);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream15 (ctx1.F, (double) k / (double) n, p, n);
        SparseMatrix<typename Ring::Element> A2_sparse (stream15);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A2_trans (A2_sparse);

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream16 (ctx1.F, p, m);
        DenseMatrix<typename Ring::Element> A3_dense (stream16);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream17 (ctx1.F, (double) k / (double) n, p, m);
        SparseMatrix<typename Ring::Element> A3_sparse (stream17);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A3_trans (A3_sparse);

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream18 (ctx1.F, m, p);
        DenseMatrix<typename Ring::Element> A4_dense (stream18);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream19 (ctx1.F, (double) k / (double) n, m, p);
        SparseMatrix<typename Ring::Element> A4_sparse (stream19);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A4_trans (A4_sparse);

        pass = testgemmConsistency (ctx1, ctx2, "sparse(row-wise)/dense/dense", A1_sparse, A2_dense, A3_dense, A1_sparse, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(col-wise)/dense/dense", A1_trans, A3_dense, A2_dense, A1_trans, A3_dense, A2_dense)  && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)/dense", A1_sparse, A2_sparse, A3_dense, A1_sparse, A2_sparse, A3_dense) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)/dense", A1_trans, A3_sparse, A2_dense, A1_trans, A3_sparse, A2_dense) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(col-wise)/dense", A3_sparse, A2_trans, A1_dense, A3_sparse, A2_trans, A1_dense) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(col-wise)/dense", A2_trans, A1_trans, A4_dense, A2_trans, A1_trans, A4_dense) && pass;

        pass = testgemmConsistency (ctx1, ctx2, "sparse(row-wise)/dense/sparse(row-wise)", A1_sparse, A2_dense, A3_sparse, A1_sparse, A2_dense, A3_sparse) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(col-wise)/dense/sparse(row-wise)", A3_trans, A1_dense, A2_sparse, A3_trans, A1_dense, A2_sparse) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)/sparse(row-wise)", A1_sparse, A2_sparse, A3_sparse, A1_sparse, A2_dense, A3_sparse) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)/sparse(row-wise)", A1_trans, A3_sparse, A2_sparse, A1_trans, A2_sparse, A3_sparse) && pass;
        pass = testgemmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(col-wise)/sparse(row-wise)", A3_sparse, A2_trans, A1_sparse, A3_sparse, A2_trans, A3_sparse) && pass;

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream41 (ctx1.F, n, n);
        DenseMatrix<typename Ring::Element> M4 (stream41);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream51 (ctx1.F, (double) k / (double) n, n, n);
        SparseMatrix<typename Ring::Element> M5 (stream51);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > M6 (M5);

        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(LT, diag = 1)", M5, M4, M5, M4, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(UT, diag = 1)", M5, M4, M5, M4,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(LT, diag != 1)", M5, M4, M5, M4, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(UT, diag != 1)", M5, M4, M5, M4, UpperTriangular, false) && pass;

        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(LT, diag = 1)", M6, M4, M6, M4, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(UT, diag = 1)", M6, M4, M6, M4,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(LT, diag != 1)", M6, M4, M6, M4, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(UT, diag != 1)", M6, M4, M6, M4, UpperTriangular, false) && pass;

        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(LT, diag = 1)", M5, M5, M5, M5, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(UT, diag = 1)", M5, M5, M5, M5,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(LT, diag != 1)", M5, M5, M5, M5, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(UT, diag != 1)", M5, M5, M5, M5, UpperTriangular, false) && pass;

        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(LT, diag = 1)", M6, M5, M6, M5, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(UT, diag = 1)", M6, M5, M6, M5,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(LT, diag != 1)", M6, M5, M6, M5, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(UT, diag != 1)", M6, M5, M6, M5, UpperTriangular, false) && pass;

        SparseMatrix<typename Ring::Element> M5p (n, n);
	BLAS3::copy (ctx1, M5, M5p);
        makeNonsingDiag (ctx1.F, M5p, false);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > M6p (M5p);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream81 (ctx1.F, (double) k / (double) n, m, n);
        SparseMatrix<typename Ring::Element> M5n (stream81);

        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(LT, diag = 1)", M5, M4, M5, M4, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(UT, diag = 1)", M5, M4, M5, M4,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(LT, diag != 1)", M5p, M4, M5p, M4, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/dense(UT, diag != 1)", M5p, M4, M5p, M4, UpperTriangular, false) && pass;

        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(LT, diag = 1)", M5, M5n, M5, M5n, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(UT, diag = 1)", M5, M5n, M5, M5n,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(LT, diag != 1)", M5p, M5, M5p, M5, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(row-wise)/sparse(row-wise)(UT, diag != 1)", M5p, M5, M5p, M5, UpperTriangular, false) && pass;

        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(LT, diag = 1)", M6, M4, M6, M4, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(UT, diag = 1)", M6, M4, M6, M4,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(LT, diag != 1)", M6p, M4, M6p, M4, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/dense(UT, diag != 1)", M6p, M4, M6p, M4, UpperTriangular, false) && pass;

        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(LT, diag = 1)", M6, M5n, M6, M5n, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(UT, diag = 1)", M6, M5n, M6, M5n,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(LT, diag != 1)", M6p, M5n, M6p, M5n, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx1, ctx2, "sparse(col-wise)/sparse(row-wise)(UT, diag != 1)", M6p, M5n, M6p, M5n, UpperTriangular, false) && pass;

        pass = testpermute_rowsConsistency(ctx1, ctx2, "sparse(row-wise)", M2, M2) && pass;
        pass = testpermute_rowsConsistency(ctx1, ctx2, "sparse(col-wise)", M3, M3) && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

        return pass;
}

template <class Ring, class Modules>
bool testBLAS3RepsConsistency (LELA::Context<Ring, Modules> &ctx, const char *text, size_t m, size_t n, size_t p, size_t k)
{
	std::ostringstream str;
        str << "Testing BLAS3 consistency over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

        bool pass = true;

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream11 (ctx.F, n, m);
        DenseMatrix<typename Ring::Element> M1 (stream11);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream21 (ctx.F, (double) k / (double) n, n, m);
        SparseMatrix<typename Ring::Element> M2 (stream21);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > M3 (M2);

        pass = testscalConsistency (ctx, ctx, "sparse(row-wise)           with dense", M2, M1) && pass;
        pass = testscalConsistency (ctx, ctx, "sparse(col-wise)           with dense", M3, M1) && pass;

	pass = testaxpyConsistency (ctx, ctx, "sparse(row-wise)/dense                with dense/dense", M2, M1, M1, M1) &&pass;
	//	pass = testaxpyConsistency (ctx, ctx, "sparse(col-wise)/dense                with dense/dense", M3, M1, M1, M1) &&pass;
        pass = testaxpyConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)     with dense/dense", M2, M2, M1, M1) &&pass;
	//	pass = testaxpyConsistency (ctx, ctx, "sparse(col-wise)/sparse(col-wise)     with dense/dense", M3, M3, M1, M1) &&pass;

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream12 (ctx.F, n, m);
        DenseMatrix<typename Ring::Element> A1_dense (stream12);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream13 (ctx.F, (double) k / (double) n, n, m);
        SparseMatrix<typename Ring::Element> A1_sparse (stream13);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A1_trans (A1_sparse);

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream14 (ctx.F, p, n);
        DenseMatrix<typename Ring::Element> A2_dense (stream14);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream15 (ctx.F, (double) k / (double) n, p, n);
        SparseMatrix<typename Ring::Element> A2_sparse (stream15);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A2_trans (A2_sparse);

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream16 (ctx.F, p, m);
        DenseMatrix<typename Ring::Element> A3_dense (stream16);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream17 (ctx.F, (double) k / (double) n, p, m);
        SparseMatrix<typename Ring::Element> A3_sparse (stream17);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A3_trans (A3_sparse);

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream18 (ctx.F, m, p);
        DenseMatrix<typename Ring::Element> A4_dense (stream18);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream19 (ctx.F, (double) k / (double) n, m, p);
        SparseMatrix<typename Ring::Element> A4_sparse (stream19);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > A4_trans (A4_sparse);

        pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/dense/dense                with dense/dense/dense",
				    A1_sparse, A2_dense, A3_dense, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/dense/dense                with dense/dense/dense",
				    A1_trans, A3_dense, A2_dense, A1_dense, A2_dense, A3_dense)  && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)/dense     with dense/dense/dense",
				    A1_sparse, A2_sparse, A3_dense, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)/dense     with dense/dense/dense",
				    A1_trans, A3_sparse, A2_dense, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/sparse(col-wise)/dense     with dense/dense/dense",
				    A3_sparse, A2_trans, A1_dense, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/sparse(col-wise)/dense     with dense/dense/dense",
				    A2_trans, A1_trans, A4_dense, A1_dense, A2_dense, A3_dense) && pass;

        pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/dense/sparse(row-wise)                with dense/dense/dense",
				    A1_sparse, A2_dense, A3_sparse, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/dense/sparse(row-wise)                with dense/dense/dense",
				    A3_trans, A1_dense, A2_sparse, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)/sparse(row-wise)     with dense/dense/dense",
				    A1_sparse, A2_sparse, A3_sparse, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)/sparse(row-wise)     with dense/dense/dense",
				    A1_trans, A3_sparse, A2_sparse, A1_dense, A2_dense, A3_dense) && pass;
        pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/sparse(col-wise)/sparse(row-wise)     with dense/dense/dense",
				    A3_sparse, A2_trans, A1_sparse, A1_dense, A2_dense, A3_dense) && pass;

        // pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/dense/sparse(col-wise)                with dense/dense/dense",
	// 			    A1_sparse, A2_dense, A4_trans, A1_dense, A2_dense, A3_dense) && pass;
        // pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/dense/sparse(col-wise)                with dense/dense/dense",
	// 			    A3_trans, A1_dense, A2_trans, A1_dense, A2_dense, A3_dense) && pass;
        // pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)/sparse(col-wise)     with dense/dense/dense",
	// 			    A3_trans, A1_sparse, A2_trans, A1_dense, A2_dense, A3_dense) && pass;
        // pass = testgemmConsistency (ctx, ctx, "sparse(row-wise)/sparse(col-wise)/sparse(col-wise)     with dense/dense/dense",
	// 			    A2_sparse, A3_trans, A1_trans, A1_dense, A2_dense, A3_dense) && pass;
        // pass = testgemmConsistency (ctx, ctx, "sparse(col-wise)/sparse(col-wise)/sparse(col-wise)     with dense/dense/dense",
	// 			    A2_trans, A1_trans, A3_trans, A1_dense, A2_dense, A3_dense) && pass;

        RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream41 (ctx.F, n, n);
        DenseMatrix<typename Ring::Element> M4 (stream41);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream51 (ctx.F, (double) k / (double) n, n, n);
        SparseMatrix<typename Ring::Element> M5 (stream51);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > M6 (M5);

        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (LT, diag = 1)", M5, M4, M4, M4, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (UT, diag = 1)", M5, M4, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (LT, diag != 1)", M5, M4, M4, M4, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (UT, diag != 1)", M5, M4, M4, M4, UpperTriangular, false) && pass;

        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (LT, diag = 1)", M6, M4, M4, M4, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (UT, diag = 1)", M6, M4, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (LT, diag != 1)", M6, M4, M4, M4, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (UT, diag != 1)", M6, M4, M4, M4, UpperTriangular, false) && pass;

        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (LT, diag = 1)", M5, M5, M4, M4, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (UT, diag = 1)", M5, M5, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (LT, diag != 1)", M5, M5, M4, M4, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (UT, diag != 1)", M5, M5, M4, M4, UpperTriangular, false) && pass;

        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (LT, diag = 1)", M6, M5, M4, M4, LowerTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (UT, diag = 1)", M6, M5, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (LT, diag != 1)", M6, M5, M4, M4, LowerTriangular, false) && pass;
        pass = testtrmmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (UT, diag != 1)", M6, M5, M4, M4, UpperTriangular, false) && pass;

        SparseMatrix<typename Ring::Element> M5p (n, n);
	BLAS3::copy (ctx, M5, M5p);
        makeNonsingDiag (ctx.F, M5p, false);

        TransposeMatrix<SparseMatrix<typename Ring::Element> > M6p (M5p);

        RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream81 (ctx.F, (double) k / (double) n, m, n);
        SparseMatrix<typename Ring::Element> M5n (stream81);

        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (LT, diag = 1)", M5, M4, M4, M4, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (UT, diag = 1)", M5, M4, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (LT, diag != 1)", M5p, M4, M4, M4, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/dense   with dense/dense (UT, diag != 1)", M5p, M4, M4, M4, UpperTriangular, false) && pass;

        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (LT, diag = 1)", M5, M5n, M4, M4, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (UT, diag = 1)", M5, M5n, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (LT, diag != 1)", M5p, M5, M4, M4, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(row-wise)/sparse(row-wise)   with dense/dense (UT, diag != 1)", M5p, M5, M4, M4, UpperTriangular, false) && pass;

        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (LT, diag = 1)", M6, M4, M4, M4, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (UT, diag = 1)", M6, M4, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (LT, diag != 1)", M6p, M4, M4, M4, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/dense   with dense/dense (UT, diag != 1)", M6p, M4, M4, M4, UpperTriangular, false) && pass;

        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (LT, diag = 1)", M6, M5n, M4, M4, LowerTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (UT, diag = 1)", M6, M5n, M4, M4,  UpperTriangular, true) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (LT, diag != 1)", M6p, M5n, M4, M4, LowerTriangular, false) && pass;
        pass = testtrsmConsistency (ctx, ctx, "sparse(col-wise)/sparse(row-wise)   with dense/dense (UT, diag != 1)", M6p, M5n, M4, M4, UpperTriangular, false) && pass;

	pass = testpermute_rowsConsistency(ctx, ctx, "sparse(row-wise)   with dense", M2, M1) && pass;
	pass = testpermute_rowsConsistency(ctx, ctx, "sparse(col-wise)   with dense", M3, M1) && pass;

        pass = testpermute_colsConsistency(ctx, ctx, "sparse(row-wise)   with dense", M2, M1) && pass;
        pass = testpermute_colsConsistency(ctx, ctx, "sparse(col-wise)   with dense", M3, M1) && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}


#endif // __LELA_TESTS_TEST_BLAS_LEVEL3_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
