
/* tests/test-matrix-domain.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
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
 * This computes the inverse of a nonsingular matrix. It throws SingularMatrix
 * in the case that no inverse exists
 */

template <class Field, class Matrix1, class Matrix2>
Matrix1 &inv (MatrixDomain<Field> &MD, Matrix1 &res, const Matrix2 &A) 
{
	linbox_check (res.coldim () == A.coldim ());
	linbox_check (res.rowdim () == A.rowdim ());

	DenseMatrix<typename Matrix1::Element> M (res.rowdim (), res.coldim () * 2);
	DenseSubmatrix<typename Matrix1::Element> M1 (M, 0, 0, res.rowdim (), res.coldim ());
	DenseSubmatrix<typename Matrix1::Element> M2 (M, 0, res.coldim (), res.rowdim (), res.coldim ());

	StandardBasisStream<Field, typename DenseSubmatrix<typename Matrix1::Element>::Row> stream (MD.field (), res.coldim ());
	typename DenseSubmatrix<typename Matrix1::Element>::RowIterator ip = M2.rowBegin ();

	for (; ip != M2.rowEnd (); ++ip)
		stream >> *ip;

	MD.copy (M1, A);

	unsigned int idx;
	typename Field::Element Mjj_inv;

	for (idx = 0; idx < M.rowdim (); ++idx) {
		if (MD.field ().isZero (M.getEntry (idx, idx))) {
			typename DenseMatrix<typename Matrix1::Element>::ColIterator col;
			typename DenseMatrix<typename Matrix1::Element>::Col::iterator i;
			unsigned int c_idx = idx + 1;

			col = M.colBegin () + idx;
			i = col->begin () + idx + 1;

			while (MD.field ().isZero (*i) && i != col->end ()) ++i, ++c_idx;

			if (i == col->end ())
				throw SingularMatrix ();
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
	}

	MD.copy (res, M2);
	return res;
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

/* Test 2: subin and isZero
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
	MD.axpy (nega, M, M1);

	report << "Output matrix M1:" << endl;
	MD.write (report, M1);

	if (!MD.isZero (M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M1 is not zero" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

/* Test 3: gemm, scal
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testGemm (Field &F, const char *text, const Matrix &A, const Matrix &B)
{
	ostringstream str;

	str << "Testing " << text << " matrix gemm" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	typename Field::Element a, b, zero, negainvb;

	F.init (zero, 0);

	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> C (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	MD.write (report, A);

	report << "Input matrix B:" << endl;
	MD.write (report, B);

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	F.write (report, a) << ", b = ";
	F.write (report, b) << std::endl;

	MD.gemm (a, A, B, zero, C);

	report << "Output matrix C := a * A * B:" << endl;
	MD.write (report, C);

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

/* Test 5: inv and mul
 *
 * Return true on success and false on failure
 */

// Version for square matrices

template <class Field, class Matrix>
static bool testInvMulSquare (Field &F, const char *text, const Matrix &M) 
{
	linbox_check (M.rowdim () == M.coldim ());

	typename Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	ostringstream str;

	str << "Testing " << text << " matrix multiplication (square)" << ends;
	commentator.start (str.str ().c_str (), "testInvMulSquare");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> Minv (M.rowdim (), M.rowdim ());
	DenseMatrix<typename Field::Element> M2 (M.rowdim (), M.rowdim ());

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrix<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	try {
		inv (MD, Minv, M);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.gemm (one, M, Minv, zero, M2);

	report << "Computed product M Minv:" << endl;
	MD.write (report, M2);

	if (!MD.areEqual (M2, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M Minv is not the identity" << endl;
		ret = false;
	}

	MD.gemm (one, Minv, M, zero, M2);

	report << "Computed product Minv M:" << endl;
	MD.write (report, M2);

	if (!MD.areEqual (M2, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMulSquare");

	return ret;
}

// Version for over-determined matrices

template <class Field, class Matrix>
static bool testInvMulOver (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.coldim () <= M.rowdim ());

	typename Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	ostringstream str;

	str << "Testing " << text << " matrix multiplication (over-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvMulOver");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> Minv (M.coldim (), M.coldim ());
	DenseMatrix<typename Field::Element> M2 (M.coldim (), M.coldim ());
	DenseMatrix<typename Field::Element> M3 (M.coldim (), M.coldim ());

	DenseMatrix<typename Field::Element> MTM (M.coldim (), M.coldim ());

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row>
		stream (F, M.coldim ());

	DenseMatrix<typename Field::Element> I (M.coldim (), M.coldim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.gemm (one, transpose (M), M, zero, MTM);

	try {
		inv (MD, Minv, MTM);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.gemm (one, transpose(M), M, zero, M2);
	MD.gemm (one, M2, Minv, zero, M3);

	report << "Computed product M^T M Minv:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M^T M Minv is not the identity" << endl;
		ret = false;
	}

	M2.resize (M.coldim (), M.rowdim ());

	MD.gemm (one, Minv, transpose (M), zero, M2);
	MD.gemm (one, M2, M, zero, M3);

	report << "Computed product Minv M^T M:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M^T M is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMulOver");

	return ret;
}

// Version for under-determined matrices

template <class Field, class Matrix>
static bool testInvMulUnder (Field &F, const char *text, Matrix &M) 
{
	linbox_check (M.rowdim () <= M.coldim ());

	typename Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	ostringstream str;

	str << "Testing " << text << " matrix multiplication (under-determined)" << ends;
	commentator.start (str.str ().c_str (), "testInvMulUnder");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> Minv (M.rowdim (), M.rowdim ());
	DenseMatrix<typename Field::Element> M2 (M.rowdim (), M.rowdim ());
	DenseMatrix<typename Field::Element> M3 (M.rowdim (), M.rowdim ());

	DenseMatrix<typename Field::Element> MMT (M.rowdim (), M.rowdim ());

	StandardBasisStream<Field, typename DenseMatrix<typename Field::Element>::Row> stream (F, M.rowdim ());

	DenseMatrix<typename Field::Element> I (M.rowdim (), M.rowdim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = I.rowBegin ();

	while (i != I.rowEnd ())
		stream >> *i++;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	MD.gemm (one, M, transpose (M), zero, MMT);

	try {
		inv (MD, Minv, MMT);
	}
	catch (SingularMatrix) {
		commentator.stop ("ok", (const char *) 0, "testInvMul");
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Matrix was found singular" << endl;
		return true;
	}

	report << "Computed inverse Minv:" << endl;
	MD.write (report, Minv);

	MD.gemm (one, M, transpose (M), zero, M2);
	MD.gemm (one, M2, Minv, zero, M3);

	report << "Computed product M M^T Minv:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix M M^T Minv is not the identity" << endl;
		ret = false;
	}

	M2.resize (M.rowdim (), M.coldim ());

	MD.gemm (one, Minv, M, zero, M2);
	MD.gemm (one, M2, transpose (M), zero, M3);

	report << "Computed product Minv M M^T:" << endl;
	MD.write (report, M3);

	if (!MD.areEqual (M3, I)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrix Minv M M^T is not the identity" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvMulUnder");

	return ret;
}

/* Test 9: m-v mul by e_i and sub
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testMVMulSub (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	typename Field::Element zero, one, neg_one;
	F.init (zero, 0);
	F.init (one, 1);
	F.init (neg_one, -1);

	str << "Testing " << text << " matrix-vector mul" << ends;
	commentator.start (str.str ().c_str (), "testMVMulSub");

	bool ret = true;

	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M:" << endl;
	MD.write (report, M);

	StandardBasisStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, M.rowdim ());
	typename LinBox::Vector<Field>::Dense v (M.coldim ());
	typename DenseMatrix<typename Field::Element>::ColIterator i = M1.colBegin ();

	for (; i != M1.colEnd (); ++i) {
		stream >> v;
		MD.gemv (one, M, v, zero, *i);
	}

	report << "Output matrix M1:" << endl;
	MD.write (report, M1);

	MD.axpy (neg_one, M, M1);

	if (!MD.isZero (M1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices M and M1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testMVMulSub");

	return ret;
}

/* Test 10: m-v axpy by e_i
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
static bool testMVAxpy (Field &F, const char *text, const Matrix &M) 
{
	ostringstream str;

	typename Field::Element zero, one, neg_one;
	F.init (zero, 0);
	F.init (one, 1);
	F.init (neg_one, -1);

	str << "Testing " << text << " matrix-vector axpy" << ends;
	commentator.start (str.str ().c_str (), "testMVAxpy");

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	DenseMatrix<typename Field::Element> M1 (M.rowdim (), M.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix M1:" << endl;
	MD.write (report, M1);

	StandardBasisStream<Field, typename LinBox::Vector<Field>::Dense> stream (F, M.rowdim ());
	typename LinBox::Vector<Field>::Dense v (M.coldim ()), w (M.rowdim ());
	typename DenseMatrix<typename Field::Element>::RowIterator i = M1.rowBegin ();

	VD.subin (w, w);

	for (; i != M1.rowEnd (); ++i) {
		stream >> v;
		MD.gemv (one, M, v, one, w);
	}

	report << "Output vector w:" << endl;
	VD.write (report, w);

	typename LinBox::Vector<Field>::Dense z (M.coldim (), one), w1 (M.rowdim ());

	MD.gemv (one, M, z, zero, w1);

	report << "Output vector w1:" << endl;
	VD.write (report, w1);

	if (!VD.areEqual (w, w1)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: VectorDomain reported vectors w and w1 are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testMVAxpy");

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

template <class Field, class Matrix>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       unsigned int iterations,
		       MatrixCategories::RowColMatrixTag)
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testScalAxpyIsZero (F, text, M1)) pass = false;
	if (!testGemm (F, text, M1, M2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testInvMulSquare (F, text, M1)) pass = false;
	}
	else if (M1.coldim () < M1.rowdim ()) {
		if (!testInvMulOver (F, text, M1)) pass = false;
	}
	else if (M1.rowdim () < M1.coldim ()) {
		if (!testInvMulUnder (F, text, M1)) pass = false;
	}

	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testPermutation (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Matrix>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       unsigned int iterations,
		       MatrixCategories::RowMatrixTag) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testScalAxpyIsZero (F, text, M1)) pass = false;
	if (!testGemm (F, text, M1, M2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testInvMulSquare (F, text, M1)) pass = false;
	}
	else if (M1.rowdim () < M1.coldim ()) {
		if (!testInvMulUnder (F, text, M1)) pass = false;
	}

	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;
	if (!testPermutation (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field, class Matrix>
bool testMatrixDomain (const Field &F, const char *text,
		       Matrix &M1, Matrix &M2, Matrix &M3,
		       unsigned int iterations,
		       MatrixCategories::ColMatrixTag) 
{
	ostringstream str;
	str << "Testing MatrixDomain with " << text << " matrices" << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	if (!testCopyEqual (F, text, M1)) pass = false;
	if (!testScalAxpyIsZero (F, text, M1)) pass = false;
	if (!testGemm (F, text, M1, M2)) pass = false;

	if (M1.rowdim () == M1.coldim ()) {
		if (!testInvMulSquare (F, text, M1)) pass = false;
	}
	else if (M1.coldim () < M1.rowdim ()) {
		if (!testInvMulOver (F, text, M1)) pass = false;
	}

	if (!testMVMulSub (F, text, M1)) pass = false;
	if (!testMVAxpy (F, text, M1)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long n = 50;
	static long m = 50;
	static long k = 10;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set row of test matrices to N.", TYPE_INT,     &n },
		{ 'm', "-m M", "Set column of test vectors to M.", TYPE_INT,     &m },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT,     &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
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

	RandomDenseStream<Field, DenseMatrix<Element>::Row> stream1 (F, m, n);
	RandomDenseStream<Field, DenseMatrix<Element>::Row> stream1p (F, m, m);

	DenseMatrix<Element> M1 (stream1); stream1.reset ();
	DenseMatrix<Element> M2 (stream1);
	DenseMatrix<Element> M3 (stream1p);

	if (!testMatrixDomain (F, "dense", M1, M2, M3, iterations,
			       MatrixTraits<DenseMatrix<Element> >::MatrixCategory ()))
		pass = false;

	RandomSparseStream<Field, SparseMatrix<Element>::Row> stream2 (F, (double) k / (double) n, m, n);
	RandomSparseStream<Field, SparseMatrix<Element>::Row> stream2p (F, (double) k / (double) n, m, m);

	SparseMatrix<Element> M4 (stream2); stream2.reset ();
	SparseMatrix<Element> M5 (stream2);
	SparseMatrix<Element> M6 (stream2p);

	if (!testMatrixDomain (F, "sparse row-wise", M4, M5, M6, iterations,
			       MatrixTraits<SparseMatrix<Element> >::MatrixCategory ()))
		pass = false;

	TransposeMatrix<SparseMatrix<Element> > M7 (M4);
	TransposeMatrix<SparseMatrix<Element> > M8 (M5);
	TransposeMatrix<SparseMatrix<Element> > M9 (M6);

	if (!testMatrixDomain (F, "sparse column-wise", M7, M8, M9, iterations,
			       MatrixTraits<TransposeMatrix<SparseMatrix<Element> > >::MatrixCategory ()))
		pass = false;

	typedef SparseVector<Element> Row;

	RandomSparseStream<Field, Row> stream3 (F, (double) k / (double) n, m, n);
	RandomSparseStream<Field, Row> stream3p (F, (double) k / (double) n, m, m);

	SparseMatrix<Element, Row> M10 (stream3); stream3.reset ();
	SparseMatrix<Element, Row> M11 (stream3);
	SparseMatrix<Element, Row> M12 (stream3p);

	if (!testMatrixDomain (F, "sparse row-wise (new rows)", M10, M11, M12, iterations,
			       MatrixTraits<SparseMatrix<Element, Row> >::MatrixCategory ()))
		pass = false;

	commentator.stop (MSG_STATUS (pass));
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
