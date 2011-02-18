/* -*- mode: c++;-*-
 * test-gauss-jordan.C
 * Test for Gauss-Jordan elimination
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include <iostream>

#include <linbox/field/modular.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/vector/stream.h>

#include "gauss-jordan.h"

using namespace LinBox;
using namespace F4;

typedef GF2 Field;

std::ostream &writePermutation (std::ostream &os, GaussJordan<Field>::Permutation P)
{
	GaussJordan<Field>::Permutation::iterator i;

	for (i = P.begin (); i != P.end (); ++i)
		os << "(" << i->first << " " << i->second << ")";

	return os;
}

void testDenseGJ (const Field &F)
{
	std::cout << "Testing DenseRowEchelonForm..." << std::endl;

	Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	int m = 96;
	int n = 96;

	GaussJordan<Field>::DenseMatrix A (m, n);
	GaussJordan<Field>::DenseMatrix R (m, n);
	GaussJordan<Field>::DenseMatrix U (m, m);
	GaussJordan<Field>::DenseMatrix UPA (m, n);

	RandomDenseStream<Field, GaussJordan<Field>::DenseMatrix::Row> A_stream (F, n);
	GaussJordan<Field>::DenseMatrix::RowIterator iter;

	for (iter = A.rowBegin (); iter != A.rowEnd (); ++iter)
		A_stream >> *iter;

	std::cout << "Number of words per row: " << A.rowBegin ()->word_size () << std::endl;

	GaussJordan<Field>::Permutation P;

	GaussJordan<Field> GJ (F);
	size_t rank;
	Field::Element det;

	MatrixDomain<Field> MD (F);

	GJ.DenseRowEchelonForm (R, U, P, A, rank, det);

	std::cout << "A = " << std::endl;
	MD.write (std::cout, A, FORMAT_SAGE);

	std::cout << "U = " << std::endl;
	MD.write (std::cout, U);

	std::cout << "P = ";
	writePermutation (std::cout, P) << std::endl;

	std::cout << "R = " << std::endl;
	MD.write (std::cout, R);

	MD.permuteRows (A, P.begin (), P.end ());

	std::cout << "PA = " << std::endl;
	MD.write (std::cout, A);

	MD.gemm (one, U, A, zero, UPA);

	std::cout << "UPA = " << std::endl;
	MD.write (std::cout, UPA);

	std::cout << "Computed rank = " << rank << std::endl;
	std::cout << "Computed det = ";
	F.write (std::cout, det);
	std::cout << std::endl;

	if (MD.areEqual (UPA, R))
		std::cout << "UPA == R, okay" << std::endl;
	else
		std::cerr << "UPA != R, not okay" << std::endl;
}

void testSparseGJ (const Field &F)
{
	GaussJordan<Field> GJ (F);

	GJ.RunTests ();

	Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	std::cout << "Testing StandardRowEchelonForm..." << std::endl;

	typedef GaussJordan<Field>::SparseMatrix Matrix;

	int k = 3;
	int m = 96;
	int n = 96;

	Matrix A (m, n), Aorig (m, n);
	GaussJordan<Field>::DenseMatrix U (m, m);
	GaussJordan<Field>::DenseMatrix UPA (m, n);

	RandomSparseStream<Field, Matrix::Row, Field::RandIter, VectorCategories::HybridZeroOneVectorTag> A_stream (F, (double) k / (double) n, n, m);
	// RandomDenseStream<Field, GaussJordan<Field>::DenseMatrix::Row> A_stream (F, n);
	Matrix::RowIterator iter;

	for (iter = A.rowBegin (); iter != A.rowEnd (); ++iter)
		A_stream >> *iter;

	GaussJordan<Field>::Permutation P;

	size_t rank;
	Field::Element det;

	MatrixDomain<Field> MD (F);

	std::cout << "A = " << std::endl;
	A.write (std::cout, F, FORMAT_PRETTY);

	MD.copy (Aorig, A);

	GJ.StandardRowEchelonForm (A, U, P, rank, det, true, true);

	std::cout << "R = " << std::endl;
	A.write (std::cout, F, FORMAT_PRETTY);

	std::cout << "U = " << std::endl;
	MD.write (std::cout, U);

	std::cout << "P = ";
	writePermutation (std::cout, P) << std::endl;

	MD.permuteRows (Aorig, P.begin (), P.end ());
	std::cout << "PA = " << std::endl;
	Aorig.write (std::cout, F, FORMAT_PRETTY);

	MD.gemm (one, U, Aorig, zero, UPA);
	std::cout << "UPA = " << std::endl;
	MD.write (std::cout, UPA);

	std::cout << "Computed rank = " << rank << std::endl;
	std::cout << "Computed det = ";
	F.write (std::cout, det);
	std::cout << std::endl;

	if (MD.areEqual (UPA, A))
		std::cout << "UPA == R, okay" << std::endl;
	else
		std::cerr << "UPA != R, not okay" << std::endl;

	std::cout << "Finished testing StandardRowEchelonForm" << std::endl;
}

#if 0

void smallTestSparseTrSolve (const Field &F)
{
	typedef GaussJordan<Field>::SparseMatrix Matrix;
	
	Matrix U (2, 2);

	Field::Element zero, one;

	F.init (zero, 0);
	F.init (one, 1);

	U.setEntry (0, 0, one);
	U.setEntry (0, 1, one);
	U.setEntry (1, 1, one);

	GaussJordan<Field>::DenseMatrix R (2, 1);

	R.setEntry (0, 0, one);
	R.setEntry (1, 0, one);

	GaussJordan<Field> GJ (F);

	GJ.SparseTrSolve (R, U);

	MatrixDomain<Field> MD (F);
	std::cout << "Solution for R = (1, 1):" << std::endl;
	MD.write (std::cout, R);

	R.setEntry (0, 0, zero);
	R.setEntry (1, 0, one);

	GJ.SparseTrSolve (R, U);

	std::cout << "Solution for R = (0, 1):" << std::endl;
	MD.write (std::cout, R);

	Matrix V (3, 3);

	V.setEntry (0, 0, one);
	V.setEntry (0, 1, one);
	V.setEntry (0, 2, one);
	V.setEntry (1, 1, one);
	V.setEntry (1, 2, one);
	V.setEntry (2, 2, one);

	std::cout << "Matrix V =" << std::endl;
	V.write (std::cout, F);

	GaussJordan<Field>::DenseMatrix S (3, 1);

	S.setEntry (0, 0, one);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, one);

	GJ.SparseTrSolve (S, V);

	std::cout << "Solution for R = (1, 1, 1):" << std::endl;
	MD.write (std::cout, S);

	S.setEntry (0, 0, zero);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, one);

	GJ.SparseTrSolve (S, V);

	std::cout << "Solution for R = (0, 1, 1):" << std::endl;
	MD.write (std::cout, S);

	Matrix V1 (3, 3);

	V1.setEntry (0, 0, one);
	V1.setEntry (0, 1, one);
	V1.setEntry (1, 1, one);
	V1.setEntry (1, 2, one);
	V1.setEntry (2, 2, one);

	std::cout << "Matrix V =" << std::endl;
	V1.write (std::cout, F);

	S.setEntry (0, 0, zero);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, one);

	GJ.SparseTrSolve (S, V1);

	std::cout << "Solution for R = (0, 1, 1):" << std::endl;
	MD.write (std::cout, S);

	Matrix W (3, 3);

	W.setEntry (0, 1, one);
	W.setEntry (0, 2, one);
	W.setEntry (1, 2, one);

	std::cout << "Matrix W =" << std::endl;
	W.write (std::cout, F);

	S.setEntry (0, 0, zero);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, zero);

	GJ.SparseTrSolve (S, W);

	std::cout << "Solution for R = (0, 1, 0):" << std::endl;
	MD.write (std::cout, S);

	Matrix W1 (3, 3);

	W1.setEntry (0, 0, one);
	W1.setEntry (0, 2, one);
	W1.setEntry (1, 2, one);

	std::cout << "Matrix W =" << std::endl;
	W1.write (std::cout, F);

	S.setEntry (0, 0, one);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, zero);

	GJ.SparseTrSolve (S, W1);

	std::cout << "Solution for R = (1, 1, 0):" << std::endl;
	MD.write (std::cout, S);

	W1.setEntry (0, 1, one);

	std::cout << "Matrix W =" << std::endl;
	W1.write (std::cout, F);

	S.setEntry (0, 0, one);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, zero);

	GJ.SparseTrSolve (S, W1);

	std::cout << "Solution for R = (1, 1, 0):" << std::endl;
	MD.write (std::cout, S);

	Matrix W2 (3, 3);

	W2.setEntry (0, 0, one);
	W2.setEntry (0, 1, one);
	W2.setEntry (1, 2, one);

	std::cout << "Matrix W =" << std::endl;
	W2.write (std::cout, F);

	S.setEntry (0, 0, one);
	S.setEntry (1, 0, one);
	S.setEntry (2, 0, zero);

	GJ.SparseTrSolve (S, W2);

	std::cout << "Solution for R = (1, 1, 0):" << std::endl;
	MD.write (std::cout, S);
}

void testSparseTrSolve (const Field &F)
{
	std::cout << "Testing SparseTrSolve..." << std::endl;

	typedef GaussJordan<Field>::SparseMatrix Matrix;

	Field::Element zero, one;
	F.init (zero, 0);
	F.init (one, 1);

	int k = 3;
	int m = 10;
	int n = 10;
	int p = 1;

	// First, obtain a random sparse matrix and put it in row-echelon form
	Matrix U (m, n);
	GaussJordan<Field>::DenseMatrix T (m, m);

	RandomSparseStream<Field, Matrix::Row> U_stream (F, (double) k / (double) n, m, n);
	Matrix::RowIterator U_iter;

	for (U_iter = U.rowBegin (); U_iter != U.rowEnd (); ++U_iter)
		U_stream >> *U_iter;

	GaussJordan<Field>::Permutation P;

	GaussJordan<Field> GJ (F);
	size_t rank;
	Field::Element det;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	GJ.StandardRowEchelonForm (U, T, P, rank, det);

	// Next, obtain a random right-hand side
	GaussJordan<Field>::DenseMatrix R (n, p), Rorig (n, p);
	RandomDenseStream<Field, GaussJordan<Field>::DenseMatrix::Row> R_stream (F, n);

	GaussJordan<Field>::DenseMatrix::RowIterator R_iter;

	for (R_iter = R.rowBegin (); R_iter != R.rowEnd (); ++R_iter)
		R_stream >> *R_iter;

	std::cout << "U = " << std::endl;
	U.write (std::cout, F, FORMAT_PRETTY);

	std::cout << "R = " << std::endl;
	MD.write (std::cout, R);

	MD.copy (Rorig, R);

	// Do triangular solve
	GJ.SparseTrSolve (R, U);

	std::cout << "U'⁻1 R = " << std::endl;
	MD.write (std::cout, R);

	// Copy rows of R to rows of L corresponding to pivots in U
	GaussJordan<Field>::DenseMatrix L (n, p);

	Matrix::RowIterator i_U;
	GaussJordan<Field>::DenseMatrix::RowIterator i_R = R.rowBegin ();

	for (i_U = U.rowBegin (); i_U != U.rowEnd () && !i_U->first.empty (); ++i_U, ++i_R)
		VD.copy (*(L.rowBegin () + i_U->first.front ()), *i_R);

	std::cout << "L = " << std::endl;
	MD.write (std::cout, L);

	GaussJordan<Field>::DenseMatrix UL (m, p);

	MD.gemm (one, U, L, zero, UL);

	std::cout << "U' U'⁻1 R = " << std::endl;
	MD.write (std::cout, UL);

	if (MD.areEqual (UL, Rorig))
		std::cout << "U'U'^-1 R == R okay" << std::endl;
	else
		std::cout << "U'U'^-1 R != R not okay" << std::endl;

	std::cout << "Finished testing SparseTrSolve" << std::endl;
}

#endif

int main (int argc, char **argv)
{
	int p = 2;
	int k = 10;

	Field F (p);

	testDenseGJ (F);
	testSparseGJ (F);
//	smallTestSparseTrSolve (F);
//	testSparseTrSolve (F);
}
