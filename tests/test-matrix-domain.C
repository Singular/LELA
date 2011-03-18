/* tests/test-matrix-domain.C
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test suite for MatrixDomain. This also effectively tests DenseMatrix,
 * SparseMatrixBase, and TransposeMatrix
 */

/* ERRORS:
 *
 * X- Ambiguous specializations with RolColMatrixTag
 *  - Can't use leftMulin or rightMulin on some tests
 *  - VectorDomain needs subin, addin, etc. with multiple vector representations
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/sparse.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/sparse.h"
#include "linbox/matrix/dense-submatrix.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

class SingularMatrix {};

template <class Matrix>
TransposeMatrix<Matrix> transpose (Matrix &M)
{
	return TransposeMatrix<Matrix> (M);
}

/* Test 1: Copy and areEqual
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testCopyEqual (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix copy, areEqual" << ends;
	commentator.start (str.str ().c_str (), "testCopyEqual");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.copy (M1, M);

	report << "Output matrix M1:" << endl;
	MD.write (report, M1);

	if (!MD.areEqual (M1, M)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M and M1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testCopyEqual");

	return ret;
}

/* Test 2: scal, axpy, isZero
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testScalAxpyIsZero (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " matrix scal, axpy, isZero" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element a, nega;

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	r.random (a);
	F.neg (nega, a);

	report << "Coefficient: ";
	F.write (report, a) << std::endl;

	MD.copy (M1, M);
	MD.scal (M1, a);

	report << "Output matrix a * M:" << endl;
	MD.write (report, M1);

	MD.axpy (nega, M, M1);

	report << "Output matrix -a * M + a * M:" << endl;
	MD.write (report, M1);

	if (!MD.isZero (M1)) {
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

template <class Field, class Matrix>
static bool testGemmCoeff (Field &F, const char *text, const Matrix &A, const Matrix &B)
{
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemm (coefficients)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element zero;

	F.init (zero, 0);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> C (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));
	typename Field::Element a, b;

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	F.write (report, a) << ", b = ";
	F.write (report, b) << std::endl;

	MD.gemm (a, A, B, zero, C);

	report << "Output matrix C := a * A * B:" << endl;
	MD.write (report, C);

	typename Field::Element negainvb;

	F.inv (negainvb, a);
	F.mulin (negainvb, b);
	F.negin (negainvb);

	report << "Coefficient: -a^-1 b = ";
	F.write (report, negainvb) << std::endl;

	MD.gemm (b, A, B, negainvb, C);

	report << "Output matrix D := b * A * B - b * a^-1 * C:" << endl;
	MD.write (report, C);

	if (!MD.isZero (C)) {
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

template <class Field, class Matrix>
static bool testGemmAssoc (Field &F, const char *text, const Matrix &A, const Matrix &B, const Matrix &C)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (B.coldim () == C.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemm (associativity)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element one, zero;

	F.init (one, 1);
	F.init (zero, 0);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> AB (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> BC (B.rowdim (), C.coldim ());
	DenseMatrix<typename Field::Element> ABpC (A.rowdim (), C.coldim ());
	DenseMatrix<typename Field::Element> ApBC (A.rowdim (), C.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	report << "Input matrix C:" << endl;
	MD.write (report, C);

	MD.gemm (one, A, B, zero, AB);

	report << "Output matrix A * B:" << endl;
	MD.write (report, AB);

	MD.gemm (one, AB, C, zero, ABpC);

	report << "Output matrix (A * B) * C:" << endl;
	MD.write (report, ABpC);

	MD.gemm (one, B, C, zero, BC);

	report << "Output matrix B * C:" << endl;
	MD.write (report, AB);

	MD.gemm (one, A, BC, zero, ApBC);

	report << "Output matrix A * (B * C):" << endl;
	MD.write (report, ApBC);

	if (!MD.areEqual (ABpC, ApBC)) {
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

template <class Field, class Matrix>
static bool testGemmIdent (Field &F, const char *text, const Matrix &A)
{
	ostringstream str;

	str << "Testing " << text << " gemm (identity)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element one, zero;

	F.init (one, 0);
	F.init (zero, 0);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> AI (A.rowdim (), A.coldim ());
	DenseMatrix<typename Field::Element> IA (A.rowdim (), A.coldim ());

	StandardBasisStream<Field, typename Matrix::Row> I_l_stream (F, A.rowdim ()), I_r_stream (F, A.coldim ());

	DenseMatrix<typename Field::Element> I_l (I_l_stream);
	DenseMatrix<typename Field::Element> I_r (I_r_stream);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Identity matrix I_l:" << endl;
	MD.write (report, I_l);

	report << "Identity matrix I_r:" << endl;
	MD.write (report, I_r);

	MD.gemm (one, I_l, A, zero, IA);

	report << "Output matrix I * A:" << endl;
	MD.write (report, IA);

	if (!MD.areEqual (A, IA)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported A != I * A" << endl;
		ret = false;
	}

	MD.gemm (one, A, I_r, zero, AI);

	report << "Output matrix A * I:" << endl;
	MD.write (report, AI);

	if (!MD.areEqual (AI, A)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported A * I != A" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 6: Row-echelon and gemm
 *
 * Return true on success and false on failure
 */

/* Elimination-step */

template <class Field, class Matrix1, class Matrix2>
void eliminate (MatrixDomain<Field> &MD, Matrix1 &M, Matrix2 &pivotRow,
		size_t row, size_t col, size_t rowdim, size_t coldim) 
{
	typename Field::Element zero, one, neg_one;
	MD.field ().init (one, 1);
	MD.field ().init (neg_one, -1);

	DenseMatrix<typename Matrix1::Element> pivotCol (rowdim, 1);
	DenseSubmatrix<typename Matrix1::Element> realPivotCol (M, row, col, rowdim, 1);
	DenseSubmatrix<typename Matrix1::Element> block (M, row, col, rowdim, coldim);

	DenseMatrix<typename Matrix1::Element> update (rowdim, coldim);

	MD.copy (pivotCol, realPivotCol);
	MD.scal (pivotCol, neg_one);
	MD.gemm (one, pivotCol, pivotRow, zero, update);
	MD.axpy (one, update, block);
}

/* Dumb elimination code
 *
 * This computes the reduced row-echelon form of a matrix, placing the
 * transform in U and the row-echelon form in R. The output should
 * satisfy R = UA. Places rank in parameter rank.
 */

template <class Field, class Matrix1, class Matrix2>
void rowEchelon (MatrixDomain<Field> &MD, Matrix1 &U, Matrix1 &R, const Matrix2 &A, size_t &rank) 
{
	linbox_check (U.rowdim () == A.rowdim ());
	linbox_check (U.coldim () == A.rowdim ());
	linbox_check (R.rowdim () == A.rowdim ());
	linbox_check (R.coldim () == A.coldim ());

	DenseMatrix<typename Matrix1::Element> M (U.rowdim (), U.coldim () + R.coldim ());
	DenseSubmatrix<typename Matrix1::Element> M1 (M, 0, 0, R.rowdim (), R.coldim ());
	DenseSubmatrix<typename Matrix1::Element> M2 (M, 0, R.coldim (), U.rowdim (), U.coldim ());

	StandardBasisStream<Field, typename DenseSubmatrix<typename Matrix1::Element>::Row> stream (MD.field (), U.coldim ());
	typename DenseSubmatrix<typename Matrix1::Element>::RowIterator ip = M2.rowBegin ();

	for (; ip != M2.rowEnd (); ++ip)
		stream >> *ip;

	MD.copy (M1, A);

	unsigned int idx;
	typename Field::Element Mjj_inv;

	rank = 0;

	for (idx = 0; idx < M.rowdim (); ++idx) {
		if (MD.field ().isZero (M.getEntry (idx, idx))) {
			typename DenseMatrix<typename Matrix1::Element>::ColIterator col;
			typename DenseMatrix<typename Matrix1::Element>::Col::iterator i;
			unsigned int c_idx = idx + 1;

			col = M.colBegin () + idx;
			i = col->begin () + idx + 1;

			while (MD.field ().isZero (*i) && i != col->end ())
				++i, ++c_idx;

			if (i == col->end ())
				continue;
			else {
				typename DenseMatrix<typename Matrix1::Element>::RowIterator row1, row2;

				row1 = M.rowBegin () + idx;
				row2 = M.rowBegin () + c_idx;

				std::swap_ranges (row1->begin () + idx, row1->end (), row2->begin () + idx);
			}
		}

		MD.field ().inv (Mjj_inv, M.getEntry (idx, idx));
		DenseSubmatrix<typename Matrix1::Element> realPivotRow (M, idx, idx, 1, M.coldim () - idx);
		MD.scal (realPivotRow, Mjj_inv);

		if (idx > 0)
			eliminate (MD, M, realPivotRow, 0, idx, idx, M.coldim () - idx);

		if (idx < M.rowdim () - 1)
			eliminate (MD, M, realPivotRow, idx + 1, idx,
				   M.rowdim () - idx - 1, M.coldim () - idx);

		++rank;
	}

	MD.copy (R, M1);
	MD.copy (U, M2);
}

template <class Field, class Matrix>
static bool testGemmRowEchelon (Field &F, const char *text, const Matrix &M) 
{
	typename Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	ostringstream str;

	str << "Testing " << text << " gemm (row-echelon)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	size_t rank;

	DenseMatrix<typename Field::Element> R (M.rowdim (), M.coldim ());
	DenseMatrix<typename Field::Element> U (M.rowdim (), M.rowdim ());
	DenseMatrix<typename Field::Element> UM (U.rowdim (), M.coldim ());

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrix<typename Field::Element> I (UM.rowdim (), UM.coldim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	rowEchelon (MD, U, R, M, rank);

	report << "Computed transform U:" << std::endl;
	MD.write (report, U);

	report << "Computed row-echelon form R:" << std::endl;
	MD.write (report, R);

	report << "Computed rank = " << rank << std::endl;

	MD.gemm (one, U, M, zero, UM);

	report << "Computed product UM:" << endl;
	MD.write (report, UM);

	if (!MD.areEqual (UM, R)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix UM != R" << endl;
		ret = false;
	}

	if (rank == M.rowdim () && rank == M.coldim ()) {
		report << "Full rank reported and matrix is square. Checking whether R is identity-matrix..." << std::endl;

		if (!MD.areEqual (UM, I)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: MatrixDomain reported matrix R is not identity" << endl;
			ret = false;
		}

		report << "Checking whether U is inverse..." << std::endl;

		MD.gemm (one, M, U, zero, UM);

		report << "Computed product MU:" << endl;
		MD.write (report, UM);

		if (!MD.areEqual (UM, R)) {
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

/* Test 7: gemv and gemm
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix1, class Matrix2>
static bool testGemvGemm (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixCategories::BlackboxTag, MatrixCategories::ColMatrixTag) 
{
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemv (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> ABgemm (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> ABgemv (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	typename Field::Element a, zero;

	F.init (zero, 0);

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	r.random (a);

	report << "Coefficient a = ";
	F.write (report, a) << std::endl;

	typename DenseMatrix<typename Field::Element>::ColIterator i_AB = ABgemv.colBegin ();
	typename Matrix2::ColIterator i_B = B.colBegin ();

	for (; i_AB != ABgemv.colEnd (); ++i_B, ++i_AB)
		MD.gemv (a, A, *i_B, zero, *i_AB);

	report << "Output matrix AB (from gemv):" << endl;
	MD.write (report, ABgemv);

	MD.gemm (a, A, B, zero, ABgemm);

	report << "Output matrix AB (from gemm):" << endl;
	MD.write (report, ABgemm);

	if (!MD.areEqual (ABgemv, ABgemm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices from gemm and gemv are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixCategories::RowMatrixTag, MatrixCategories::BlackboxTag) 
{
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemv (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> ABgemm (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> ABgemv (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	typename Field::Element a, zero;

	F.init (zero, 0);

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	r.random (a);

	report << "Coefficient a = ";
	F.write (report, a) << std::endl;

	typename DenseMatrix<typename Field::Element>::RowIterator i_AB = ABgemv.rowBegin ();
	typename Matrix1::ConstRowIterator i_A = A.rowBegin ();

	for (; i_AB != ABgemv.rowEnd (); ++i_A, ++i_AB)
		MD.gemv (a, transpose (B), *i_A, zero, *i_AB);

	report << "Output matrix AB (from gemv):" << endl;
	MD.write (report, ABgemv);

	MD.gemm (a, A, B, zero, ABgemm);

	report << "Output matrix AB (from gemm):" << endl;
	MD.write (report, ABgemm);

	if (!MD.areEqual (ABgemv, ABgemm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices from gemm and gemv are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixCategories::BlackboxTag, MatrixCategories::RowColMatrixTag) 
{
	return testGemvGemmSpecialised (F, text, A, B, MatrixCategories::BlackboxTag (), MatrixCategories::ColMatrixTag ());
}

template <class Field, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixCategories::RowColMatrixTag, MatrixCategories::BlackboxTag) 
{
	return testGemvGemmSpecialised (F, text, A, B, MatrixCategories::RowMatrixTag (), MatrixCategories::BlackboxTag ());
}

template <class Field, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag) 
{
	return testGemvGemmSpecialised (F, text, A, B, MatrixCategories::RowMatrixTag (), MatrixCategories::BlackboxTag ());
}

template <class Field, class Matrix1, class Matrix2>
static bool testGemvGemm (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B) 
{
	return testGemvGemmSpecialised (F, text, A, B,
					typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
					typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ());
}

/* Test 8: gemv with coefficients
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix, class Vector>
static bool testGemvCoeff (Field &F, const char *text, const Matrix &A, const Vector &v)
{
	ostringstream str;

	str << "Testing " << text << " gemv (coefficients)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element a, b, zero, negainvb;

	F.init (zero, 0);

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	std::vector<typename Field::Element> aAv (A.rowdim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Input vector v: ";
	VD.write (report, v) << endl;

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	F.write (report, a) << ", b = ";
	F.write (report, b) << std::endl;

	MD.gemv (a, A, v, zero, aAv);

	report << "Output vector a * A * v: ";
	VD.write (report, aAv) << std::endl;

	F.inv (negainvb, a);
	F.mulin (negainvb, b);
	F.negin (negainvb);

	report << "Coefficient: -a^-1 b = ";
	F.write (report, negainvb) << std::endl;

	MD.gemv (b, A, v, negainvb, aAv);

	report << "Output vector w := b * A * v - b * a^-1 * aAv: ";
	VD.write (report, aAv) << std::endl;

	if (!VD.isZero (aAv)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported vector w is not zero" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 9: trsv
 *
 * Return true on success and false on failure
 */

/* Construct a random nonsingular upper triangular matrix from a random matrix */

template <class Field, class Matrix>
static Matrix &makeUpperTriangular (Field &F, Matrix &A)
{
	typename Field::Element zero;

	F.init (zero, 0);

	size_t i, j;

	for (i = 0; i < A.rowdim (); ++i) {
		for (j = 0; j < std::min (i, A.coldim ()); ++j) {
			if (!F.isZero (A.getEntry (i, j))) {
				A.setEntry (i, j, zero);
				A.eraseEntry (i, j);
			}
		}
	}

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	for (i = 0; i < std::min (A.rowdim (), A.coldim ()); ++i)
		if (F.isZero (A.getEntry (i, i)))
			r.random (A.refEntry (i, i));

	return A;
}

template <class Field, class Matrix, class Vector>
static bool testTrsv (Field &F, const char *text, const Matrix &A, const Vector &v)
{
	ostringstream str;

	str << "Testing " << text << " trsv" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	typename Field::Element zero, one;

	F.init (zero, 0);
	F.init (one, 1);

	Matrix U (A.rowdim (), A.coldim ());

	MD.copy (U, A);
	makeUpperTriangular (F, U);

	typename LinBox::Vector<Field>::Dense Uinv_v (U.coldim ()), UUinv_v (U.rowdim ());
	VD.copy (Uinv_v, v);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << endl;
	MD.write (report, U);

	report << "Input vector v: ";
	VD.write (report, Uinv_v) << endl;

	MD.trsv (U, Uinv_v);

	report << "Output vector U^-1 * v: ";
	VD.write (report, Uinv_v) << std::endl;

	MD.gemv (one, U, Uinv_v, zero, UUinv_v);

	report << "Output vector U U^-1 * v: ";
	VD.write (report, UUinv_v) << std::endl;

	if (!VD.areEqual (v, UUinv_v)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: VectorxDomain reported v != U * U^-1 * v" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 10: trsm and trsv
 *
 * B must be iterable by columns
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix1, class Matrix2>
static bool testTrsmTrsv (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B) 
{
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " trsm (consistency with trsv)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	typename Field::Element one;

	F.init (one, 1);

	Matrix1 U (A.rowdim (), A.coldim ());

	MD.copy (U, A);
	makeUpperTriangular (F, U);

	DenseMatrix<typename Field::Element> UinvBtrsm (U.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> UinvBtrsv (U.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << endl;
	MD.write (report, U);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	MD.copy (UinvBtrsm, B);
	MD.copy (UinvBtrsv, B);

	typename DenseMatrix<typename Field::Element>::ColIterator i_UinvB = UinvBtrsv.colBegin ();

	for (; i_UinvB != UinvBtrsv.colEnd (); ++i_UinvB)
		MD.trsv (U, *i_UinvB);

	report << "Output matrix U^-1 * B (from trsv):" << endl;
	MD.write (report, UinvBtrsv);

	MD.trsm (one, U, UinvBtrsm);

	report << "Output matrix U^-1 * B (from trsm):" << endl;
	MD.write (report, UinvBtrsm);

	if (!MD.areEqual (UinvBtrsv, UinvBtrsm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices from trsm and trsv are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 11: trsm with coefficient
 *
 * B must be iterable by columns
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix1, class Matrix2>
static bool testTrsmCoeff (Field &F, const char *text, const Matrix1 &A, const Matrix2 &B) 
{
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " trsm (coefficient)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	typename Field::Element one;

	F.init (one, 1);

	Matrix1 U (A.rowdim (), A.coldim ());

	MD.copy (U, A);
	makeUpperTriangular (F, U);

	DenseMatrix<typename Field::Element> aUinvB (U.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> UinvB (U.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << endl;
	MD.write (report, U);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	typename Field::RandIter r (F);
	typename Field::Element a;

	r.random (a);

	report << "Coefficient: a = ";
	F.write (report, a) << std::endl;

	MD.trsm (a, U, aUinvB);

	report << "Output matrix a U^-1 * B:" << endl;
	MD.write (report, aUinvB);

	MD.trsm (one, U, UinvB);

	report << "Output matrix U^-1 * B:" << endl;
	MD.write (report, UinvB);

	MD.scal (UinvB, a);

	report << "Output matrix a (U^-1 * B):" << endl;
	MD.write (report, UinvB);

	if (!MD.areEqual (UinvB, aUinvB)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported a U^-1 * B != a (U^-1 * B)" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

std::ostream &reportPermutation
	(std::ostream &out,
	 const std::vector<std::pair<unsigned int, unsigned int> > &P) 
{
	std::vector<std::pair<unsigned int, unsigned int> >::const_iterator i;

	for (i = P.begin (); i != P.end (); ++i)
		out << "(" << i->first << " " << i->second << ")";

	return out;
}

template <class Field, class Matrix>
bool testPermutation (const Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	str << "Testing " << text << " permutations" << ends;
	commentator.start (str.str ().c_str (), "testPermutation");

	bool ret = true;

	MatrixDomain<Field> MD (F);
	MersenneTwister MT (time (NULL));

	typename MatrixDomain<Field>::Permutation P, Pinv;

	// Create a random row permutation
	for (unsigned int i = 0; i < M.rowdim (); ++i) {
		unsigned int row1, row2;

		do {
			row1 = MT.randomInt () % M.rowdim ();
			row2 = MT.randomInt () % M.rowdim ();
		} while (row1 == row2);

		P.push_back (typename MatrixDomain<Field>::Transposition (row1, row2));
	}

	// Construct the inverse of this transposition
	Pinv.resize (P.size ());
	std::copy (P.begin (), P.end (), Pinv.begin ());
	std::reverse (Pinv.begin (), Pinv.end ());

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Permutation P:    ";
	reportPermutation (report, P) << endl;
	report << "Permutation P^-1: ";
	reportPermutation (report, Pinv) << endl;

	// Apply the permutation and then its inverse to a copy of M
	Matrix M1 (M);

	MD.permuteRows (M1, P.begin (), P.end ());
 	MD.permuteRows (M1, Pinv.begin (), Pinv.end ());

	report << "Output matrix P^-1 PM:" << endl;
	MD.write (report, M1);

	// Compare M and M1
	if (!MD.areEqual (M, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: M != P^-1 PM" << endl;
		ret = false;
	}

	// Now we do exactly the same thing with columns
	P.clear ();
	Pinv.clear ();

	// Create a random column permutation
	for (unsigned int i = 0; i < M.coldim (); ++i) {
		unsigned int col1, col2;

		do {
			col1 = MT.randomInt () % M.coldim ();
			col2 = MT.randomInt () % M.coldim ();
		} while (col1 == col2);

		P.push_back (typename MatrixDomain<Field>::Transposition (col1, col2));
	}

	// Construct the inverse of this transposition
	Pinv.resize (P.size ());
	std::copy (P.begin (), P.end (), Pinv.begin ());
	std::reverse (Pinv.begin (), Pinv.end ());

	report << "Permutation P:    ";
	reportPermutation (report, P) << endl;
	report << "Permutation P^-1: ";
	reportPermutation (report, Pinv) << endl;

	// Apply the permutation and then its inverse to a copy of M
	MD.permuteColumns (M1, P.begin (), P.end ());
	MD.permuteColumns (M1, Pinv.begin (), Pinv.end ());

	report << "Output matrix MPP^-1:" << endl;
	MD.write (report, M1);

	// Compare M and M1
	if (!MD.areEqual (M, M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: M != MPP^-1" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testPermutation");

	return ret;
}

template <class Field, class Matrix, class Vector>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       Vector &v,
		       unsigned int iterations,
		       MatrixCategories::RowColMatrixTag)
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testScalAxpyIsZero (F, text, M1)) pass = false;
	if (!testGemmCoeff (F, text, M1, M2)) pass = false;
	if (!testGemmAssoc (F, text, M1, M2, M3)) pass = false;
	if (!testGemmRowEchelon (F, text, M1)) pass = false;
	if (!testGemvGemm (F, text, M1, M2)) pass = false;
	if (!testGemvCoeff (F, text, M1, v)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsv (F, text, M1, v)) pass = false;
		if (!testTrsmTrsv (F, text, M1, M2)) pass = false;
		if (!testTrsmCoeff (F, text, M1, M2)) pass = false;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}

	if (!testPermutation (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Matrix, class Vector>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       Vector &v,
		       unsigned int iterations,
		       MatrixCategories::RowMatrixTag) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testScalAxpyIsZero (F, text, M1)) pass = false;
	if (!testGemmCoeff (F, text, M1, M2)) pass = false;
	if (!testGemmAssoc (F, text, M1, M2, M3)) pass = false;
	if (!testGemmRowEchelon (F, text, M1)) pass = false;
	if (!testGemvGemm (F, text, M1, M2)) pass = false;
	if (!testGemvCoeff (F, text, M1, v)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsv (F, text, M1, v)) pass = false;
		if (!testTrsmTrsv (F, text, M1, M2)) pass = false;
		if (!testTrsmCoeff (F, text, M1, M2)) pass = false;
	} else {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Input matrix M1 is not square, so skipping tests of trsv and trsm" << std::endl;
	}

	if (!testPermutation (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Matrix, class Vector>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       Vector &v,
		       unsigned int iterations,
		       MatrixCategories::ColMatrixTag) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testScalAxpyIsZero (F, text, M1)) pass = false;
	if (!testGemmCoeff (F, text, M1, M2)) pass = false;
	if (!testGemmAssoc (F, text, M1, M2, M3)) pass = false;
	if (!testGemmRowEchelon (F, text, M1)) pass = false;
//	if (!testGemvGemm (F, text, M1, M2)) pass = false;    Not compiling because the compiler won't instantiate with TransposeMatrix. Hmmm....
	if (!testGemvCoeff (F, text, M1, v)) pass = false;

#if 0 // disabled because of TransposeMatrix-issues
	if (M1.rowdim () == M1.coldim ()) {
		if (!testTrsv (F, text, M1, v)) pass = false;

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

int main (int argc, char **argv)
{
	bool pass = true;

	static long l = 50;
	static long n = 50;
	static long m = 40;
	static long p = 30;
	static long k = 10;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set row-dimension of matrix C and column-dimension of B to N.", TYPE_INT, &n },
		{ 'p', "-p P", "Set column-dimension of matrix C to P.", TYPE_INT, &p },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint32> Field;
	typedef Field::Element Element;

	Field F (q);

	commentator.start("Matrix domain test suite", "MatrixDomain");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	RandomDenseStream<Field, Vector<Field>::Dense> stream_v1 (F, m, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v2 (F, n, 1);

	Vector<Field>::Dense v1 (m), v2 (n);
	stream_v1 >> v1;
	stream_v2 >> v2;

	RandomDenseStream<Field, DenseMatrix<Element>::Row> stream11 (F, m, l);
	RandomDenseStream<Field, DenseMatrix<Element>::Row> stream12 (F, n, m);
	RandomDenseStream<Field, DenseMatrix<Element>::Row> stream13 (F, p, n);

	DenseMatrix<Element> M1 (stream11);
	DenseMatrix<Element> M2 (stream12);
	DenseMatrix<Element> M3 (stream13);

	if (!testMatrixDomain (F, "dense", M1, M2, M3, v1, iterations,
			       MatrixTraits<DenseMatrix<Element> >::MatrixCategory ()))
		pass = false;

	RandomSparseStream<Field, SparseMatrix<Element>::Row> stream21 (F, (double) k / (double) m, m, l);
	RandomSparseStream<Field, SparseMatrix<Element>::Row> stream22 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Field, SparseMatrix<Element>::Row> stream23 (F, (double) k / (double) p, p, n);

	SparseMatrix<Element> M4 (stream21);
	SparseMatrix<Element> M5 (stream22);
	SparseMatrix<Element> M6 (stream23);

	if (!testMatrixDomain (F, "sparse row-wise", M4, M5, M6, v1, iterations,
			       MatrixTraits<SparseMatrix<Element> >::MatrixCategory ()))
		pass = false;

	TransposeMatrix<SparseMatrix<Element> > M7 (M6);
	TransposeMatrix<SparseMatrix<Element> > M8 (M5);
	TransposeMatrix<SparseMatrix<Element> > M9 (M4);

	if (!testMatrixDomain (F, "sparse column-wise", M7, M8, M9, v2, iterations,
			       MatrixTraits<TransposeMatrix<SparseMatrix<Element> > >::MatrixCategory ()))
		pass = false;

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
