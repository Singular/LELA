/* tests/test-blas-level3.h
 * Copyright 2001, 2002, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Generic test suite for BLAS Level 3 routines
 */

#ifndef __LINBOX_TESTS_TEST_BLAS_LEVEL3_H
#define __LINBOX_TESTS_TEST_BLAS_LEVEL3_H

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/ring/gf2.h"
#include "linbox/blas/context.h"
#include "linbox/blas/level1.h"
#include "linbox/blas/level2.h"
#include "linbox/blas/level3.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/dense.h"

using namespace std;
using namespace LinBox;

typedef std::vector<std::pair<uint32, uint32> > Permutation;

/* Force the input-matrix to have nonsingular entries on the diagonal */

template <class Field, class Matrix>
static Matrix &makeNonsingDiag (Field &F, Matrix &A)
{
	size_t i;

	typename Field::Element a;

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	for (i = 0; i < std::min (A.rowdim (), A.coldim ()); ++i) {
		if (!A.getEntry (a, i, i) || F.isZero (a)) {
			r.random (a);
			A.setEntry (i, i, a);
		}
	}

	return A;
}

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

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
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
	report << "Input matrix M:" << endl;
	BLAS3::write (ctx, report, M);

	r.random (a);
	ctx.F.neg (nega, a);

	report << "Coefficient: ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::copy (ctx, M, M1);
	BLAS3::scal (ctx, a, M1);

	report << "Output matrix a * M:" << endl;
	BLAS3::write (ctx, report, M1);

	BLAS3::axpy (ctx, nega, M, M1);

	report << "Output matrix -a * M + a * M:" << endl;
	BLAS3::write (ctx, report, M1);

	if (!BLAS3::is_zero (ctx, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M1 is not zero" << endl;
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
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemm (coefficients)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> C (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B);

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));
	typename Field::Element a, b;

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	ctx.F.write (report, a) << ", b = ";
	ctx.F.write (report, b) << std::endl;

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), C);

	report << "Output matrix C := a * A * B:" << endl;
	BLAS3::write (ctx, report, C);

	typename Field::Element negainvb;

	ctx.F.inv (negainvb, a);
	ctx.F.mulin (negainvb, b);
	ctx.F.negin (negainvb);

	report << "Coefficient: -a^-1 b = ";
	ctx.F.write (report, negainvb) << std::endl;

	BLAS3::gemm (ctx, b, A, B, negainvb, C);

	report << "Output matrix D := b * A * B - b * a^-1 * C:" << endl;
	BLAS3::write (ctx, report, C);

	if (!BLAS3::is_zero (ctx, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix D is not zero" << endl;
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
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (B.coldim () == C.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemm (associativity)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> AB (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> BC (B.rowdim (), C.coldim ());
	DenseMatrix<typename Field::Element> ABpC (A.rowdim (), C.coldim ());
	DenseMatrix<typename Field::Element> ApBC (A.rowdim (), C.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
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
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported (A * B) * C != A * (B * C)" << endl;
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

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> I_l_stream (ctx.F, A.rowdim ()), I_r_stream (ctx.F, A.coldim ());

	DenseMatrix<typename Field::Element> I_l (I_l_stream);
	DenseMatrix<typename Field::Element> I_r (I_r_stream);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
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
			<< "ERROR: MatrixDomain reported A != I * A" << endl;
		ret = false;
	}

	BLAS3::gemm (ctx, ctx.F.one (), A, I_r, ctx.F.zero (), AI);

	report << "Output matrix A * I:" << endl;
	BLAS3::write (ctx, report, AI);

	if (!BLAS3::equal (ctx, AI, A)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported A * I != A" << endl;
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

	bool ret = true;

	typename Matrix1::ContainerType A1 (A.rowdim (), A.coldim ());
	typename Matrix2::ContainerType B1 (B.rowdim (), B.coldim ());
	typename Matrix2::ContainerType C (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, A1);
	BLAS3::copy (ctx, B, B1);

	makeUpperTriangular (ctx.F, A1, false);

	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A1);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B1);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trmm (ctx, a, A, B1, UpperTriangular, false);

	report << "Output matrix a * A * B (trmm): " << std::endl;
	BLAS3::write (ctx, report, B1);

	BLAS3::gemm (ctx, a, A1, B, ctx.F.zero (), C);

	report << "Output matrix a * A * B (gemm): " << std::endl;
	BLAS3::write (ctx, report, C);

	if (!BLAS3::equal (ctx, B1, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported results from trmm and gemm are not equal" << endl;
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

	bool ret = true;

	typename Matrix1::ContainerType A1 (A.rowdim (), A.coldim ());
	typename Matrix2::ContainerType B1 (B.rowdim (), B.coldim ());
	typename Matrix2::ContainerType C (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, A1);
	BLAS3::copy (ctx, B, B1);

	makeLowerTriangular (ctx.F, A1, false);

	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A1);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B1);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trmm (ctx, a, A, B1, LowerTriangular, false);

	report << "Output matrix a * A * B (trmm): " << std::endl;
	BLAS3::write (ctx, report, B1);

	BLAS3::gemm (ctx, a, A1, B, ctx.F.zero (), C);

	report << "Output matrix a * A * B (gemm): " << std::endl;
	BLAS3::write (ctx, report, C);

	if (!BLAS3::equal (ctx, B1, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported results from trmm and gemm are not equal" << endl;
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
	linbox_check (U.rowdim () == A.rowdim ());
	linbox_check (U.coldim () == A.rowdim ());
	linbox_check (R.rowdim () == A.rowdim ());
	linbox_check (R.coldim () == A.coldim ());

	typename Matrix1::ContainerType M (U.rowdim (), U.coldim () + R.coldim ());
	typename Matrix1::ContainerType::SubmatrixType M1 (M, 0, 0, R.rowdim (), R.coldim ());
	typename Matrix1::ContainerType::SubmatrixType M2 (M, 0, R.coldim (), U.rowdim (), U.coldim ());

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

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> stream (ctx.F, M.rowdim ());

	DenseMatrix<typename Field::Element> I (UM.rowdim (), UM.coldim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
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
			<< "ERROR: MatrixDomain reported matrix UM != R" << endl;
		ret = false;
	}

	if (rank == M.rowdim () && rank == M.coldim ()) {
		report << "Full rank reported and matrix is square. Checking whether R is identity-matrix..." << std::endl;

		if (!BLAS3::equal (ctx, UM, I)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: MatrixDomain reported matrix R is not identity" << endl;
			ret = false;
		}

		report << "Checking whether U is inverse..." << std::endl;

		BLAS3::gemm (ctx, ctx.F.one (), M, U, ctx.F.zero (), UM);

		report << "Computed product MU:" << endl;
		BLAS3::write (ctx, report, UM);

		if (!BLAS3::equal (ctx, UM, R)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: MatrixDomain reported matrix MU is not identity" << endl;
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

	Matrix U (A.rowdim (), A.coldim ());
	Matrix B1 (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, U);
	makeNonsingDiag (ctx.F, U);

	BLAS3::copy (ctx, B, B1);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << std::endl;
	BLAS3::write (ctx, report, U);

	report << "Input matrix B: " << std::endl;
	BLAS3::write (ctx, report, B);

	typename Field::Element a, ainv;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);
	ctx.F.inv (ainv, a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trsm (ctx, a, U, B1, LowerTriangular, false);

	report << "Output matrix U^-1 * B: " << std::endl;
	BLAS3::write (ctx, report, B1);

	BLAS3::trmm (ctx, ainv, U, B1, LowerTriangular, false);

	report << "Output matrix U U^-1 * B: " << std::endl;
	BLAS3::write (ctx, report, B1);

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

	Matrix U (A.rowdim (), A.coldim ());
	Matrix B1 (B.rowdim (), B.coldim ());

	BLAS3::copy (ctx, A, U);
	makeNonsingDiag (ctx.F, U);

	BLAS3::copy (ctx, B, B1);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << std::endl;
	BLAS3::write (ctx, report, U);

	report << "Input matrix B: " << std::endl;
	BLAS3::write (ctx, report, B);

	typename Field::Element a, ainv;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);
	ctx.F.inv (ainv, a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trsm (ctx, a, U, B1, UpperTriangular, false);

	report << "Output matrix U^-1 * B: " << std::endl;
	BLAS3::write (ctx, report, B1);

	BLAS3::trmm (ctx, ainv, U, B1, UpperTriangular, false);

	report << "Output matrix U U^-1 * B: " << std::endl;
	BLAS3::write (ctx, report, B1);

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
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " trsm (coefficient)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	Matrix1 U (A.rowdim (), A.coldim ());

	BLAS3::copy (ctx, A, U);
	makeUpperTriangular (ctx.F, U, true);

	DenseMatrix<typename Field::Element> aUinvB (U.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> UinvB (U.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << endl;
	BLAS3::write (ctx, report, U);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B);

	BLAS3::copy (ctx, B, aUinvB);
	BLAS3::copy (ctx, B, UinvB);

	typename Field::RandIter r (ctx.F);
	typename Field::Element a;

	r.random (a);

	report << "Coefficient: a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS3::trsm (ctx, a, U, aUinvB, UpperTriangular, false);

	report << "Output matrix a U^-1 * B:" << endl;
	BLAS3::write (ctx, report, aUinvB);

	BLAS3::trsm (ctx, ctx.F.one (), U, UinvB, UpperTriangular, false);

	report << "Output matrix U^-1 * B:" << endl;
	BLAS3::write (ctx, report, UinvB);

	BLAS3::scal (ctx, a, UinvB);

	report << "Output matrix a (U^-1 * B):" << endl;
	BLAS3::write (ctx, report, UinvB);

	if (!BLAS3::equal (ctx, UinvB, aUinvB)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported a U^-1 * B != a (U^-1 * B)" << endl;
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

	MersenneTwister MT (time (NULL));

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
	catch (LinboxError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Caught LinboxError: " << e << endl;
		pass = false;
	}

#ifdef __LINBOX_HAVE_LIBPNG
	if (format != FORMAT_PNG)
#endif // __LINBOX_HAVE_LIBPNG
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
	catch (LinboxError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Caught LinboxError: " << e << endl;
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

#ifdef __LINBOX_HAVE_LIBPNG
	FileFormatTag formats[] = { FORMAT_TURNER, FORMAT_DUMAS, FORMAT_MATLAB, FORMAT_PRETTY, FORMAT_PNG };
#else // !__LINBOX_HAVE_LIBPNG
	FileFormatTag formats[] = { FORMAT_TURNER, FORMAT_DUMAS, FORMAT_MATLAB, FORMAT_PRETTY };
#endif // __LINBOX_HAVE_LIBPNG

	for (size_t i = 0; i < sizeof (formats) / sizeof (FileFormatTag); ++i)
		pass = testReadWriteFormat (ctx, text, M, formats[i]) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3 (Context<Field, Modules> &ctx, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       unsigned int iterations,
		       MatrixIteratorTypes::RowCol)
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
	if (!testTrmmGemmUpper (ctx, text, M1, M2)) pass = false;
	if (!testTrmmGemmLower (ctx, text, M1, M2)) pass = false;
	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsmLower (ctx, text, M1, M2)) pass = false;
		if (!testTrsmUpper (ctx, text, M1, M2)) pass = false;
		if (!testTrsmCoeff (ctx, text, M1, M2)) pass = false;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}

	if (!testPermutation (ctx, text, M1)) pass = false;
	if (!testReadWrite (ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3 (Context<Field, Modules> &ctx, const char *text,
		Matrix &M1, Matrix &M2, Matrix &M3,
		unsigned int iterations,
		MatrixIteratorTypes::Row) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
	if (!testTrmmGemmUpper (ctx, text, M1, M2)) pass = false;
	if (!testTrmmGemmLower (ctx, text, M1, M2)) pass = false;
//	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;  // Needs ColIterator

	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsmLower (ctx, text, M1, M2)) pass = false;
		if (!testTrsmUpper (ctx, text, M1, M2)) pass = false;
		if (!testTrsmCoeff (ctx, text, M1, M2)) pass = false;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}

	if (!testPermutation (ctx, text, M1)) pass = false;
	if (!testReadWrite (ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3 (Context<Field, Modules> &ctx, const char *text,
		Matrix &M1, Matrix &M2, Matrix &M3,
		unsigned int iterations,
		MatrixIteratorTypes::Col) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
	if (!testTrmmGemmUpper (ctx, text, M1, M2)) pass = false;
	if (!testTrmmGemmLower (ctx, text, M1, M2)) pass = false;
	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;

#if 0
	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsmLower (ctx, text, M1, M2)) pass = false;
		if (!testTrsmUpper (ctx, text, M1, M2)) pass = false;

		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M2 is not iterable by columns, so skipping tests of trsm" << std::endl;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}
#endif

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Modules, class Matrix>
bool testBLAS3Submatrix (Context<Field, Modules> &ctx, const char *text, // 
			 const Matrix &M1, const Matrix &M2, const Matrix &M3,
			 unsigned int iterations,
			 MatrixIteratorTypes::Row) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (ctx, text, M1)) pass = false;
	if (!testScalAxpyIsZero (ctx, text, M1)) pass = false;
	if (!testGemmCoeff (ctx, text, M1, M2)) pass = false;
	if (!testGemmAssoc (ctx, text, M1, M2, M3)) pass = false;
	if (!testGemmIdent (ctx, text, M1)) pass = false;
//	if (!testGemmRowEchelon (ctx, text, M1)) pass = false;  // Needs ColIterator

#if 0 // Disabled until we have a good candidate here
	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsmCoeff (ctx, text, M1, M2)) pass = false;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}
#endif

	if (!testPermutation (ctx, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

#endif // __LINBOX_TESTS_TEST_MATRIX_DOMAIN_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
