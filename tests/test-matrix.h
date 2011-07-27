/* lela/tests/test-matrix.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic tests for matrices
 * 
 * ------------------------------------
 *
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_MATRIX_H
#define __LELA_TESTS_TEST_MATRIX_H

#include <iostream>
#include <sstream>

#include "lela/util/commentator.h"
#include "lela/blas/context.h"
#include "lela/blas/level3.h"
#include "lela/matrix/transpose.h"

using namespace LELA;

/* Test 1: getEntry and Row/ColIterator (reading)
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
bool testRowIterator (const Field &F, const Matrix &M)
{
	commentator.start ("Testing RowIterator", __FUNCTION__);

	bool pass = true;

	typename Matrix::ConstRowIterator i_M;
	size_t i, j;

	typename Field::Element a;

	for (i_M = M.rowBegin (), i = 0; i_M != M.rowEnd (); ++i_M, ++i) {
		if (i >= M.rowdim ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Matrix RowIterator did not end where it was supposed to (" << i << ")" << std::endl;
			pass = false;
			break;
		}

		for (j = 0; j < M.coldim (); ++j) {
			if (M.getEntry (a, i, j)) {
				typename Field::Element b;

				if (!VectorUtils::getEntry<typename Field::Element, typename Matrix::ConstRow> (*i_M, b, j)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: Entry at position (" << i << "," << j << ") is not present in row-vector (should be ";
					F.write (error, a) << ")" << std::endl;
					pass = false;
				}
				else if (!F.areEqual (a, b)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: Entry at position (" << i << "," << j << ") is ";
					F.write (error, b) << " (should be ";
					F.write (error, a) << ")" << std::endl;
					pass = false;
				}
			} else {
				typename Field::Element b;

				if (VectorUtils::getEntry<typename Field::Element, typename Matrix::ConstRow> (*i_M, b, j) && !F.isZero (b)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: getEntry reports entry at position (" << i << "," << j << ") does not exist, but ref reports ";
					F.write (error, b) << std::endl;
					pass = false;
				}
			}
		}
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Field, class Matrix>
bool testColIterator (const Field &F, const Matrix &M)
{
	commentator.start ("Testing ColIterator", __FUNCTION__);

	bool pass = true;

	typename Matrix::ConstColIterator j_M;
	size_t i, j;

	typename Field::Element a;

	for (j_M = M.colBegin (), j = 0; j_M != M.colEnd (); ++j_M, ++j) {
		if (j >= M.coldim ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Matrix ColIterator did not end where it was supposed to (" << j << ")" << std::endl;
			pass = false;
			break;
		}

		for (i = 0; i < M.rowdim (); ++i) {
			if (M.getEntry (a, i, j)) {
				typename Field::Element b;

				if (!VectorUtils::getEntry<typename Field::Element, typename Matrix::ConstCol> (*j_M, b, i)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: Entry at position (" << i << "," << j << ") is not present in column-vector (should be ";
					F.write (error, a) << ")" << std::endl;
					pass = false;
				}
				else if (!F.areEqual (a, b)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: Entry at position (" << i << "," << j << ") is ";
					F.write (error, b) << " (should be ";
					F.write (error, a) << ")" << std::endl;
					pass = false;
				}
			} else {
				typename Field::Element b;

				if (VectorUtils::getEntry<typename Field::Element, typename Matrix::ConstCol> (*j_M, b, i) && !F.isZero (b)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: getEntry reports entry at position (" << i << "," << j << ") does not exist, but ref reports ";
					F.write (error, b) << std::endl;
					pass = false;
				}
			}
		}
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Field, class Matrix>
bool testRowColIteratorSpecialised (const Field &F, const Matrix &M, MatrixIteratorTypes::Row)
{
	return testRowIterator (F, M);
}

template <class Field, class Matrix>
bool testRowColIteratorSpecialised (const Field &F, const Matrix &M, MatrixIteratorTypes::Col)
{
	return testColIterator (F, M);
}

template <class Field, class Matrix>
bool testRowColIteratorSpecialised (const Field &F, const Matrix &M, MatrixIteratorTypes::RowCol)
{
	bool pass = true;

	pass = testRowIterator (F, M) && pass;
	pass = testColIterator (F, M) && pass;

	return pass;
}

template <class Field, class Matrix>
bool testRowColIterator (const Field &F, const Matrix &M)
{
	return testRowColIteratorSpecialised (F, M, typename Matrix::IteratorType ());
}

/* Test 2: getEntry, setEntry
 *
 * Return true on success and false on failure
 */
template <class Field, class Matrix>
bool testGetSetEntry (const Field &F, const Matrix &M)
{
	commentator.start ("Testing getEntry, setEntry", __FUNCTION__);

	bool pass = true;

	Matrix M1 (M.rowdim (), M.coldim ());

	Context<Field> ctx (F);

	BLAS3::copy (ctx, M, M1);

	NonzeroRandIter <Field> r (F, typename Field::RandIter (F));
	typename Field::Element a, b;
	
	size_t i, j;

	for (i = 0; i < M.rowdim (); ++i) {
		for (j = 0; j < M.coldim (); ++j) {
			r.random (a);
			M1.setEntry (i, j, a);

			if (!M1.getEntry (b, i, j)) {
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Entry M_(" << i << "," << j << ") does not exist after setEntry to nonzero value" << std::endl;
				pass = false;
				continue;
			}

			if (!F.areEqual (a, b)) {
				std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
				error << "ERROR: (setEntry, getEntry) M_(" << i << "," << j << ") = ";
				F.write (error, b) << " (should be ";
				F.write (error, a) << ")" << std::endl;
				pass = false;
			}
		}
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

/* Test 3: setEntry, eraseEntry
 *
 * Should only be called on sparse matrices, since it will always fail
 * on dense matrices.
 *
 * Return true on success and false on failure
 */
template <class Field, class Matrix>
bool testSetEraseEntry (const Field &F, const Matrix &M)
{
	commentator.start ("Testing setEntry, eraseEntry", __FUNCTION__);

	bool pass = true;

	Matrix M1 (M.rowdim (), M.coldim ());

	Context<Field> ctx (F);

	BLAS3::copy (ctx, M, M1);

	NonzeroRandIter <Field> r (F, typename Field::RandIter (F));
	typename Field::Element a, b;

	size_t i, j;

	for (i = 0; i < M.rowdim (); ++i) {
		for (j = 0; j < M.coldim (); ++j) {
			r.random (a);
			M1.setEntry (i, j, a);
			M1.eraseEntry (i, j);

			if (M1.getEntry (b, i, j)) {
				std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
				error << "ERROR: M_(" << i << "," << j << ") = ";
				F.write (error, b) << " (should not exist)" << std::endl;
				error << "Entry had been set to ";
				F.write (error, a) << std::endl;
				pass = false;
			}
		}
	}

	if (!BLAS3::is_zero (ctx, M1)) {
		std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		error << "ERROR: MatrixDomain reports that M is not zero" << std::endl;
		error << "State of M:" << std::endl;
		BLAS3::write (ctx, error, M1);
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

/* Test 4: rowdim, coldim
 *
 * Return true on success and false on failure
 */
template <class Field, class Matrix>
bool testRowdimColdim (const Field &F, const Matrix &M, size_t rowdim, size_t coldim)
{
	commentator.start ("Testing rowdim, coldim", __FUNCTION__);

	bool pass = true;

	if (M.rowdim () != rowdim) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Incorrect row-dimension " << M.rowdim () << " (should be " << rowdim << ")" << std::endl;
		pass = false;
	}

	if (M.coldim () != coldim) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Incorrect column-dimension " << M.coldim () << " (should be " << coldim << ")" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

/* Test 5: RawIterator, RawIndexedIterator
 *
 * Return true on success and false on failure
 */
template <class Field, class Matrix>
bool testRawIterator (const Field &F, const Matrix &M)
{
	commentator.start ("Testing RawIterator, RawIndexedIterator", __FUNCTION__);

	bool pass = true;

	typename Matrix::ConstRawIterator i_elt;
	typename Matrix::ConstRawIndexedIterator i_idx;

	typename Field::Element a;

	bool occupied[M.rowdim ()][M.coldim ()];

	size_t i, j;

	for (i = 0; i < M.rowdim (); ++i)
		for (j = 0; j < M.coldim (); ++j)
			occupied[i][j] = false;

	for (i_idx = M.rawIndexedBegin (), i_elt = M.rawBegin (); i_idx != M.rawIndexedEnd (); ++i_idx, ++i_elt) {
		if (i_idx->first > M.rowdim ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Row-index " << i_idx->first << " out of range (0.." << M.rowdim () << ")" << std::endl;
			pass = false;
			continue;
		}

		if (i_idx->second > M.coldim ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Column-index " << i_idx->second << " out of range (0.." << M.coldim () << ")" << std::endl;
			pass = false;
			continue;
		}

		occupied[i_idx->first][i_idx->second] = true;

		if (!M.getEntry (a, i_idx->first, i_idx->second)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: M_(" << i_idx->first << "," << i_idx->second << ") listed in RawIndexedIterator but not present in matrix according to getEntry" << std::endl;
			pass = false;
		}
		else if (!F.areEqual (a, *i_elt)) {
			std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "ERROR: M_(" << i_idx->first << "," << i_idx->second << ") = ";
			F.write (error, a) << " (from getEntry) and ";
			F.write (error, *i_elt) << " (from RawIterator)" << std::endl;
			pass = false;
		}
	}

	for (i = 0; i < M.rowdim (); ++i) {
		for (j = 0; j < M.coldim (); ++j) {
			if (!occupied[i][j] && M.getEntry (a, i, j)) {
				std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
				error << "ERROR: M_(" << i << "," << j << ") = ";
				F.write (error, a) << " but was not hit by RawIndexedIterator" << std::endl;
				pass = false;
			}
		}
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

/* Test 6: Submatrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Matrix>
bool testSubmatrixDim (const Field &F, Matrix &M, size_t row_begin, size_t col_begin, size_t rowdim, size_t coldim)
{
	std::ostringstream str;
	str << "Testing Submatrix, start: (" << row_begin << ", " << col_begin << "), dimensions: (" << rowdim << ", " << coldim << ")" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool pass = true;

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	typename Matrix::ConstSubmatrixType Mp (M, row_begin, col_begin, rowdim, coldim);

	Context<Field> ctx (F);

	report << "Submatrix:" << std::endl;
	BLAS3::write (ctx, report, Mp);

	if (Mp.startRow () != row_begin) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Incorrect start row " << Mp.startRow () << " (should be " << row_begin << ")" << std::endl;
		pass = false;
	}

	if (Mp.startCol () != col_begin) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Incorrect start row " << Mp.startCol () << " (should be " << col_begin << ")" << std::endl;
		pass = false;
	}

	if (Mp.rowdim () != rowdim) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Incorrect row-dimension " << Mp.rowdim () << " (should be " << rowdim << ")" << std::endl;
		pass = false;
	}

	if (Mp.coldim () != coldim) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Incorrect column-dimension " << Mp.coldim () << " (should be " << coldim << ")" << std::endl;
		pass = false;
	}

	std::ostringstream str1;

	str1 << __FUNCTION__ << "_check1" << std::ends;
	commentator.start ("Checking entry-equality with getEntry", str1.str ().c_str ());

	bool pass1 = true;

	size_t i, j;

	typename Field::Element a, b;

	for (i = 0; i < rowdim; ++i) {
		for (j = 0; j < coldim; ++j) {
			if (Mp.getEntry (a, i, j)) {
				if (!M.getEntry (b, i + row_begin, j + col_begin)) {
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR: M'_(" << i << "," << j << ") exists, but M_(" << i + row_begin << "," << j + col_begin << ") does not" << std::endl;
					pass1 = false;
				}
				else if (!F.areEqual (a, b)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: M'_(" << i << "," << j << ") != M_(" << i + row_begin << "," << j + col_begin << ")" << std::endl
					      << "M'_(" << i << "," << j << "): ";
					F.write (error, a) << std::endl << "M_(" << i + row_begin << "," << j + col_begin << "): ";
					F.write (error, b) << std::endl;
					pass1 = false;
				}
			} else {
				if (M.getEntry (b, i + row_begin, j + col_begin)) {
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR: M'_(" << i << "," << j << ") does not exist, but M_(" << i + row_begin << "," << j + col_begin << ") does" << std::endl;
					pass1 = false;
				}
			}
		}
	}

	commentator.stop (MSG_STATUS (pass1), (const char *) 0, str1.str ().c_str ());

	pass = pass && pass1;

	pass = testRowColIterator (F, Mp) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Field, class Matrix>
bool testSubmatrix (const Field &F, Matrix &M)
{
	commentator.start ("Testing Submatrix", __FUNCTION__);

	bool pass = true;

	size_t i, j, n, m;

	for (i = 0; i < M.rowdim (); ++i)
		for (j = 0; j < M.coldim (); ++j)
			for (n = 0; n < M.rowdim () - i; ++n)
				for (m = 0; m < M.coldim () - j; ++m)
					pass = testSubmatrixDim (F, M, i, j, n, m) && pass;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

/* Test 7: TransposeMatrix
 *
 * Return true on success and false on failure
 */
template <class Field, class Matrix>
bool testTransposeMatrix (const Field &F, Matrix &M)
{
	commentator.start ("Testing TransposeMatrix", __FUNCTION__);

	bool pass = true;

	TransposeMatrix<Matrix> MT (M);

	size_t i, j;

	typename Field::Element a, b;

	for (i = 0; i < M.rowdim (); ++i) {
		for (j = 0; j < M.coldim (); ++j) {
			if (M.getEntry (a, i, j)) {
				if (!MT.getEntry (b, j, i)) {
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR: M_(" << i << "," << j << ") exists, but M^T_(" << j << "," << i << ") does not" << std::endl;
					pass = false;
				}
				else if (!F.areEqual (a, b)) {
					std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
					error << "ERROR: M_(" << i << "," << j << ") != M^T_(" << j << "," << i << ")" << std::endl
					      << "M_(" << i << "," << j << "): ";
					F.write (error, a) << std::endl << "M^T_(" << j << "," << i << "): ";
					F.write (error, b) << std::endl;
					pass = false;
				}
			} else {
				if (MT.getEntry (b, j, i)) {
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR: M_(" << i << "," << j << ") does not exist, but M^T_(" << j << "," << i << ") does" << std::endl;
					pass = false;
				}
			}
		}
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

#endif // __LELA_TESTS_TEST_MATRIX_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
