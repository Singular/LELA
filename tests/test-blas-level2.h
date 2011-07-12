/* tests/test-blas-level2.h
 * Copyright 2001, 2002, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Generic test suite for BLAS level 2 routines
 */

#ifndef __LINBOX_TESTS_TEST_BLAS_LEVEL2_H
#define __LINBOX_TESTS_TEST_BLAS_LEVEL2_H

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
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A1);

	report << "Input vector u: ";
	BLAS1::write (ctx, report, u) << endl;

	report << "Input vector v: ";
	BLAS1::write (ctx, report, v) << endl;

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	BLAS2::ger (ctx, a, u, v, A1);

	report << "Output matrix a * u * v^T + A: " << std::endl;
	BLAS3::write (ctx, report, A1);

	BLAS3::gemm (ctx, a, U, V, ctx.F.one (), A2);

	report << "Output matrix a * U * V + A: " << std::endl;
	BLAS3::write (ctx, report, A2);

	if (!BLAS3::equal (ctx, A1, A2)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported a * u * v^T + A != a * U * V + A" << endl;
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
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemv (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> ABgemm (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> ABgemv (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	typename DenseMatrix<typename Field::Element>::ColIterator i_AB = ABgemv.colBegin ();
	typename Matrix2::ConstColIterator i_B = B.colBegin ();

	for (; i_AB != ABgemv.colEnd (); ++i_B, ++i_AB)
		BLAS2::gemv (ctx, a, A, *i_B, ctx.F.zero (), *i_AB);

	report << "Output matrix AB (from gemv):" << endl;
	BLAS3::write (ctx, report, ABgemv);

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), ABgemm);

	report << "Output matrix AB (from gemm):" << endl;
	BLAS3::write (ctx, report, ABgemm);

	if (!BLAS3::equal (ctx, ABgemv, ABgemm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices from gemm and gemv are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Modules, class Matrix1, class Matrix2>
static bool testGemvGemmSpecialised (Context<Field, Modules> &ctx, const char *text, const Matrix1 &A, const Matrix2 &B, MatrixIteratorTypes::Row, MatrixIteratorTypes::Generic) 
{
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	str << "Testing " << text << " gemv (consistency with gemm)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	DenseMatrix<typename Field::Element> ABgemm (A.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> ABgemv (A.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B);

	typename Field::Element a;

	NonzeroRandIter<Field> r (ctx.F, typename Field::RandIter (ctx.F));

	r.random (a);

	report << "Coefficient a = ";
	ctx.F.write (report, a) << std::endl;

	typename DenseMatrix<typename Field::Element>::RowIterator i_AB = ABgemv.rowBegin ();
	typename Matrix1::ConstRowIterator i_A = A.rowBegin ();

	for (; i_AB != ABgemv.rowEnd (); ++i_A, ++i_AB)
		BLAS2::gemv (ctx, a, transpose (B), *i_A, ctx.F.zero (), *i_AB);

	report << "Output matrix AB (from gemv):" << endl;
	BLAS3::write (ctx, report, ABgemv);

	BLAS3::gemm (ctx, a, A, B, ctx.F.zero (), ABgemm);

	report << "Output matrix AB (from gemm):" << endl;
	BLAS3::write (ctx, report, ABgemm);

	if (!BLAS3::equal (ctx, ABgemv, ABgemm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices from gemm and gemv are not equal" << endl;
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

	typename LinBox::Vector<Field>::Dense aAv (A.rowdim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	BLAS3::write (ctx, report, A);

	report << "Input vector v: ";
	BLAS1::write (ctx, report, v) << endl;

	r.random (a);
	r.random (b);

	report << "Coefficients: a = ";
	ctx.F.write (report, a) << ", b = ";
	ctx.F.write (report, b) << std::endl;

	BLAS2::gemv (ctx, a, A, v, ctx.F.zero (), aAv);

	report << "Output vector a * A * v: ";
	BLAS1::write (ctx, report, aAv) << std::endl;

	ctx.F.inv (negainvb, a);
	ctx.F.mulin (negainvb, b);
	ctx.F.negin (negainvb);

	report << "Coefficient: -a^-1 b = ";
	ctx.F.write (report, negainvb) << std::endl;

	BLAS2::gemv (ctx, b, A, v, negainvb, aAv);

	report << "Output vector w := b * A * v - b * a^-1 * aAv: ";
	BLAS1::write (ctx, report, aAv) << std::endl;

	if (!BLAS1::is_zero (ctx, aAv)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported vector w is not zero" << endl;
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
	linbox_check (A.coldim () == B.rowdim ());

	ostringstream str;

	const char *type_str[] = { "upper", "lower" };

	str << "Testing " << text << " trsm (consistency with trsv, " << type_str[type] << " triangular matrix)" << ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	bool ret = true;

	Matrix1 U (A.rowdim (), A.coldim ());

	BLAS3::copy (ctx, A, U);
	makeNonsingDiag (ctx.F, U);

	DenseMatrix<typename Field::Element> UinvBtrsm (U.rowdim (), B.coldim ());
	DenseMatrix<typename Field::Element> UinvBtrsv (U.rowdim (), B.coldim ());

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix U:" << endl;
	BLAS3::write (ctx, report, U);

	report << "Input matrix B:" << endl;
	BLAS3::write (ctx, report, B);

	BLAS3::copy (ctx, B, UinvBtrsm);
	BLAS3::copy (ctx, B, UinvBtrsv);

	typename DenseMatrix<typename Field::Element>::ColIterator i_UinvB = UinvBtrsv.colBegin ();

	for (; i_UinvB != UinvBtrsv.colEnd (); ++i_UinvB)
		BLAS2::trsv (ctx, U, *i_UinvB, type, false);

	report << "Output matrix U^-1 * B (from trsv):" << endl;
	BLAS3::write (ctx, report, UinvBtrsv);

	BLAS3::trsm (ctx, ctx.F.one (), U, UinvBtrsm, type, false);

	report << "Output matrix U^-1 * B (from trsm):" << endl;
	BLAS3::write (ctx, report, UinvBtrsm);

	if (!BLAS3::equal (ctx, UinvBtrsv, UinvBtrsm)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: MatrixDomain reported matrices from trsm and trsv are not equal" << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, __FUNCTION__);

	return ret;
}

template <class Field, class Modules, class Matrix1, class Matrix2, class Vector1, class Vector2>
bool testBLAS2 (Context<Field, Modules> &ctx, const char *text,
		Matrix1 &M1, Matrix2 &M2,
		Vector1 &v1, Vector2 &v2,
		unsigned int iterations,
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
		unsigned int iterations,
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
		unsigned int iterations,
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
			 unsigned int iterations,
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
bool testgemvConsistency (LinBox::Context<Ring, Modules1> &ctx1,
			  LinBox::Context<Ring, Modules2> &ctx2,
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

	bool pass = true;
	
	typename VectorTraits<Ring,Vector3>::ContainerType w1;
	typename VectorTraits<Ring,Vector4>::ContainerType w2;
	typename VectorTraits<Ring,Vector2>::ContainerType w3;

	VectorUtils::ensureDim<Ring, Vector3> (w1, v1.size ());
	VectorUtils::ensureDim<Ring, Vector4> (w2, v2.size ());

	VectorUtils::ensureDim<Ring, Vector2> (w3, v2.size ());

	BLAS1::copy (ctx1, v1, w1);
	BLAS1::copy (ctx1, v2, w2);
	BLAS1::copy (ctx1, v2, w3);

	typename Matrix2::ContainerType A;

	BLAS3::copy(ctx1,M1,A);

	typename Ring::Element a,b;
	NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

	r.random (a);
	r.random (b);

	report << "Coefficient a: ";
	ctx1.F.write (report, a) << std::endl;

	report << "Coefficient b: ";
	ctx1.F.write (report, b) << std::endl;

	report << "Vector v_1: ";
	BLAS1::write (ctx1, report, v1) << std::endl;

	report << "Vector v_2: ";
	BLAS1::write (ctx1, report, v2) << std::endl;

	report << "Matrix A: "<< std::endl;
	BLAS3::write (ctx1, report, A);


	BLAS2::gemv(ctx1, a, M1, v1, b, w3);



	report << "Vector : a Av_1 + bv_2";
	BLAS1::write (ctx1, report, w3) << std::endl;

	report << "Vector w_1: ";
	BLAS1::write (ctx2, report, w1) << std::endl;

	report << "Vector w_2: ";
	BLAS1::write (ctx2, report, w2) << std::endl;

	BLAS2::gemv (ctx2, a, A, w1, b, w2);

	report << "Vector a Aw_1 + bw_2: ";
	BLAS1::write (ctx1, report, w2) << std::endl;

	if (!BLAS1::equal (ctx1, w3, w2))
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: a Av_1 + bv_2 != a Aw_1 + bw_2" << std::endl;
		pass = false;
	}


	commentator.stop (MSG_STATUS (pass));

	return pass;
}	

template <class Ring, class Modules>
bool testBLAS2Consistency (LinBox::Context<Ring, Modules> &ctx, const char *text, size_t m, size_t n, size_t k) 
{
	std::ostringstream str;
	str << "Testing BLAS1 consistency over <" << text << ">" << std::ends;
	LinBox::commentator.start (str.str ().c_str ());

	bool pass = true;

	typename LinBox::Vector<Ring>::Dense v1(n), v2(m);
	typename LinBox::Vector<Ring>::Sparse w1, w2;


	RandomDenseStream<Ring, typename DenseMatrix<typename Ring::Element>::Row> stream11 (ctx.F, n, m); 
	DenseMatrix<typename Ring::Element> M1 (stream11);

	RandomSparseStream<Ring, typename SparseMatrix<typename Ring::Element>::Row> stream21 (ctx.F, (double) k / (double) n, n, m);
	SparseMatrix<typename Ring::Element> M2 (stream21);

	TransposeMatrix<SparseMatrix<typename Ring::Element> > M3 (M2);


	LinBox::RandomDenseStream<Ring, typename LinBox::Vector<Ring>::Dense> stream1 (ctx.F, n), stream2 (ctx.F, m);
	LinBox::RandomSparseStream<Ring, typename LinBox::Vector<Ring>::Sparse> stream3 (ctx.F, (double) k / (double) n, n), stream4 (ctx.F, (double) k / (double) n, m);

	stream1 >>  v1;
	stream2 >>  v2;
	stream3 >>  w1;
	stream4 >>  w2;
	
	pass = testgemvConsistency (ctx, ctx, "sparse(row-wise)	/dense	/dense 		with dense/dense/dense", M2, v1, v2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(col-wise)	/dense	/dense 		with dense/dense/dense", M3, v1, v2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(row-wise)	/sparse	/sparse 	with dense/dense/dense", M2, w1, w2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(col-wise)	/sparse	/sparse 	with dense/dense/dense", M3, w1, w2, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(row-wise)	/sparse	/dense		with dense/dense/dense", M2, w1, v1, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "sparse(col-wise)	/sparse	/dense		with dense/dense/dense", M3, w1, v1, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "dense		/sparse	/dense		with dense/dense/dense", M1, w1, v1, M1, v1, v1) && pass; 
	pass = testgemvConsistency (ctx, ctx, "dense		/sparse	/sparse		with dense/dense/dense", M1, w1, w2, M1, v1, v1) && pass; 


	LinBox::commentator.stop (MSG_STATUS (pass));

	return pass;
}







#endif // __LINBOX_TESTS_TEST_BLAS_LEVEL2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
