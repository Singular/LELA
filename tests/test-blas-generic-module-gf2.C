/* tests/test-blas-generic-module-gf2.C
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test suite for BLAS-routines using GenericModule over GF(2)
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/util/commentator.h"
#include "lela/blas/context.h"
#include "lela/matrix/dense.h"
#include "lela/matrix/sparse.h"
#include "lela/vector/stream.h"
#include "lela/ring/gf2.h"

#include "test-common.h"
#include "test-blas-level1.h"
#include "test-blas-level2.h"
#include "test-blas-level3.h"

using namespace LELA;

typedef GF2 Field;

template <class Modules>
bool testBLAS1 (Context<GF2, Modules> &ctx, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing VectorDomain <" << text << ">" << std::ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<GF2, Vector<GF2>::Dense> stream1 (ctx.F, n, iterations), stream2 (ctx.F, n, iterations);
	RandomSparseStream<GF2, Vector<GF2>::Sparse> stream3 (ctx.F, 0.1, n, iterations), stream4 (ctx.F, 0.1, n, iterations);
	RandomHybridStream<GF2, Vector<GF2>::Hybrid> stream5 (ctx.F, 0.1, n, iterations), stream6 (ctx.F, 0.1, n, iterations);

	if (!testCopyEqual (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "dense/sparse", stream1, stream4)) pass = false;  stream1.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "dense/hybrid", stream1, stream6)) pass = false;  stream1.reset (); stream6.reset ();
	if (!testCopyEqual (ctx, "sparse/dense", stream3, stream2)) pass = false;  stream3.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "sparse/hybrid", stream3, stream6)) pass = false; stream3.reset (); stream6.reset ();
	if (!testCopyEqual (ctx, "hybrid/dense", stream5, stream2)) pass = false;  stream5.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "hybrid/sparse", stream5, stream4)) pass = false; stream5.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "hybrid/hybrid", stream5, stream6)) pass = false; stream5.reset (); stream6.reset ();

	if (!testInequality (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testInequality (ctx, "dense/sparse", stream1, stream4)) pass = false;  stream1.reset (); stream4.reset ();
	if (!testInequality (ctx, "dense/hybrid", stream1, stream6)) pass = false;  stream1.reset (); stream6.reset ();
	if (!testInequality (ctx, "sparse/dense", stream3, stream2)) pass = false;  stream3.reset (); stream2.reset ();
	if (!testInequality (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();
	if (!testInequality (ctx, "sparse/hybrid", stream3, stream6)) pass = false; stream3.reset (); stream6.reset ();
	if (!testInequality (ctx, "hybrid/dense", stream5, stream2)) pass = false;  stream5.reset (); stream2.reset ();
	if (!testInequality (ctx, "hybrid/sparse", stream5, stream4)) pass = false; stream5.reset (); stream4.reset ();
	if (!testInequality (ctx, "hybrid/hybrid", stream5, stream6)) pass = false; stream5.reset (); stream6.reset ();

	if (!testNonzero (ctx, "dense", stream1)) pass = false;   stream1.reset ();
	if (!testNonzero (ctx, "sparse", stream3)) pass = false;  stream3.reset ();
	if (!testNonzero (ctx, "hybrid", stream5)) pass = false;  stream5.reset ();

	if (!testDotProduct (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testDotProduct (ctx, "dense/sparse", stream1, stream3)) pass = false;  stream1.reset (); stream3.reset ();
	if (!testDotProduct (ctx, "dense/hybrid", stream1, stream5)) pass = false;  stream1.reset (); stream5.reset ();
	if (!testDotProduct (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();
	if (!testDotProduct (ctx, "sparse/hybrid", stream3, stream5)) pass = false; stream3.reset (); stream5.reset ();
	if (!testDotProduct (ctx, "hybrid/hybrid", stream5, stream6)) pass = false; stream5.reset (); stream6.reset ();

	if (!testScal (ctx, "dense", stream1)) pass = false;   stream1.reset ();
	if (!testScal (ctx, "sparse", stream3)) pass = false;  stream3.reset ();
	if (!testScal (ctx, "hybrid", stream5)) pass = false;  stream5.reset ();

	if (!testAXPY (ctx, "dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testAXPY (ctx, "sparse", stream3, stream4)) pass = false;  stream3.reset (); stream4.reset ();
	if (!testAXPY (ctx, "hybrid", stream5, stream6)) pass = false;  stream5.reset (); stream6.reset ();

	if (!testDotConsistency(ctx, ctx, "dense/sparse with dense/dense", stream1, stream3, stream1, stream2)) pass = false; stream1.reset (); stream3.reset();
        if (!testDotConsistency(ctx, ctx, "dense/hybrid with dense/dense", stream1, stream5, stream1, stream2)) pass = false; stream1.reset (); stream5.reset();       
	if (!testDotConsistency(ctx, ctx, "hybrid/sparse with dense/dense", stream5, stream3, stream1, stream2)) pass = false; stream3.reset(); stream5.reset();
	if (!testDotConsistency(ctx, ctx, "hybrid/dense with dense/dense", stream5, stream1, stream1, stream2)) pass = false; stream1.reset (); stream5.reset();
	if (!testDotConsistency(ctx, ctx, "sparse/dense with dense/dense", stream3, stream1, stream1, stream2)) pass = false; stream1.reset (); stream3.reset();
	if (!testDotConsistency(ctx, ctx, "sparse/hybrid with dense/dense", stream3, stream5, stream1, stream2)) pass = false; stream3.reset(); stream5.reset ();

	if (!testAxpyConsistency(ctx, ctx, "dense/sparse with dense/dense", stream1, stream3, stream1, stream2)) pass = false; stream1.reset ();  stream3.reset();
        if (!testAxpyConsistency(ctx, ctx, "dense/hybrid with dense/dense", stream1, stream5, stream1, stream2)) pass = false; stream1.reset ();  stream5.reset();
        if (!testAxpyConsistency(ctx, ctx, "hybrid/sparse with dense/dense", stream5, stream3, stream1, stream2)) pass = false; stream3.reset(); stream5.reset();
        if (!testAxpyConsistency(ctx, ctx, "hybrid/dense with dense/dense", stream5, stream1, stream1, stream2)) pass = false; stream1.reset (); stream5.reset();
        if (!testAxpyConsistency(ctx, ctx, "sparse/dense with dense/dense", stream3, stream1, stream1, stream2)) pass = false; stream1.reset (); stream3.reset();
        if (!testAxpyConsistency(ctx, ctx, "sparse/hybrid with dense/dense", stream3, stream5, stream1, stream2)) pass = false; stream3.reset(); stream5.reset ();

	if (!testScalConsistency(ctx, ctx, "sparse with dense", stream3, stream1)) pass = false; stream1. reset (); stream3.reset ();
        if (!testScalConsistency(ctx, ctx, "hybrid with dense", stream5, stream1)) pass = false; stream1. reset (); stream5.reset ();

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Modules>
bool testBLAS2Consistency (Context<GF2, Modules> &ctx, const char *text, size_t n, size_t m, size_t k)
{
	std::ostringstream str;
        str << "Testing BLAS2 Consistency <" << text << ">" << std::ends;
        commentator.start (str.str ().c_str ());

        bool pass = true;

        RandomDenseStream<GF2, Vector<GF2>::Dense> stream1 (ctx.F, n, 1), stream2 (ctx.F, m, 1);
	Vector<GF2>::Dense v1d (n);
	stream1 >> v1d;
	Vector<GF2>::Dense v2d (m);
	stream2 >> v2d;
        RandomSparseStream<GF2, Vector<GF2>::Sparse> stream3 (ctx.F, 0.1, n, 1), stream4 (ctx.F, 0.1, m, 1);
	Vector<GF2>::Sparse v1s;
	stream3 >> v1s;
	Vector<GF2>::Sparse v2s;
	stream4 >> v2s;
        RandomHybridStream<GF2, Vector<GF2>::Hybrid> stream5 (ctx.F, 0.1, n, 1), stream6 (ctx.F, 0.1, m, 1);
	Vector<GF2>::Hybrid v1h;
	stream5 >> v1h;
	Vector<GF2>::Hybrid v2h;
	stream6 >> v2h;

	RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream7 (ctx.F, n, m);
	DenseMatrix<bool> A1 (stream7);
	TransposeMatrix<DenseMatrix<bool> > A1t (A1); 
        RandomSparseStream<GF2, SparseMatrix<bool>::Row> stream9 (ctx.F, (double) k / (double) n, n, m);
	SparseMatrix<bool> A2 (stream9);
	TransposeMatrix<SparseMatrix<bool> > A2t (A2);
        RandomHybridStream<GF2, SparseMatrix<bool,Vector<GF2>::Hybrid>::Row> stream11 (ctx.F, (double) k / (double) n, n, m);
	SparseMatrix<bool, Vector<GF2>::Hybrid> A3 (stream11);
	TransposeMatrix<SparseMatrix<bool, Vector<GF2>::Hybrid> > A3t (A3);

        if( !testgemvConsistency (ctx, ctx, "sparse(row-wise)/dense/dense	with dense/dense/dense", A2, v1d, v2d, A1, v1d, v2d)) pass = false;  
        if( !testgemvConsistency (ctx, ctx, "sparse(col-wise)/dense/dense	with dense/dense/dense", A2t, v2d, v1d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "sparse(row-wise)/sparse/sparse	with dense/dense/dense", A2, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "sparse(col-wise)/sparse/sparse	with dense/dense/dense", A2t, v2s, v1s, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "sparse(row-wise)/sparse/dense	with dense/dense/dense", A2, v1s, v2d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "sparse(col-wise)/sparse/dense	with dense/dense/dense", A2t, v2s, v1d, A1, v1d, v2d)) pass = false;

        if( !testgemvConsistency (ctx, ctx, "hybrid(row-wise)/dense/dense       with dense/dense/dense", A3, v1d, v2d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "hybrid(col-wise)/dense/dense       with dense/dense/dense", A3t, v2d, v1d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "hybrid(row-wise)/sparse/sparse     with dense/dense/dense", A3, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "hybrid(col-wise)/sparse/sparse     with dense/dense/dense", A3t, v2s, v1s, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "hybrid(row-wise)/sparse/dense      with dense/dense/dense", A3, v1s, v2d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "hybrid(col-wise)/sparse/dense      with dense/dense/dense", A3t, v2s, v1d, A1, v1d, v2d)) pass = false;

        if( !testgemvConsistency (ctx, ctx, "dense/sparse/dense			with dense/dense/dense", A1, v1s, v2d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "dense/sparse/sparse		with dense/dense/dense", A1, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "dense/dense/sparse			with dense/dense/dense", A1, v1d, v2s, A1, v1d, v2d)) pass = false;

        if( !testgemvConsistency (ctx, ctx, "dense(col-wise)/sparse/dense       with dense/dense/dense", A1t, v2s, v1d, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "dense(col-wise)/sparse/sparse      with dense/dense/dense", A1t, v2s, v1s, A1, v1d, v2d)) pass = false;
        if( !testgemvConsistency (ctx, ctx, "dense(col-wise)/dense/sparse       with dense/dense/dense", A1t, v2d, v1s, A1, v1d, v2d)) pass = false;

        RandomDenseStream<GF2, DenseMatrix<bool>::Row> stream12 (ctx.F, n, n);
        DenseMatrix<bool> Aq1 (stream12);
        TransposeMatrix<DenseMatrix<bool> > Aq1t (Aq1);
        RandomSparseStream<GF2, SparseMatrix<bool>::Row> stream13 (ctx.F, (double) k / (double) n, n, n);
        SparseMatrix<bool> Aq2 (stream13);
        TransposeMatrix<SparseMatrix<bool> > Aq2t(Aq2);
        RandomHybridStream<GF2, SparseMatrix<bool,Vector<GF2>::Hybrid>::Row> stream14 (ctx.F, 0.1, n, n);
        SparseMatrix<bool, Vector<GF2>::Hybrid> Aq3 (stream14);
        TransposeMatrix<SparseMatrix<bool, Vector<GF2>::Hybrid> > Aq3t (Aq3);
	   
	if( !testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag = 1)", Aq2, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag = 1)", Aq2, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag != 1)", Aq2, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag != 1)", Aq2, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag = 1)", Aq2t, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag = 1)", Aq2t, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag != 1)", Aq2t, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag != 1)", Aq2t, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;

        if( !testtrmvConsistency (ctx, ctx, "hybrid(row-wise) with dense (LT, diag = 1)", Aq3, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(row-wise) with dense (UT, diag = 1)", Aq3, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(row-wise) with dense (LT, diag != 1)", Aq3, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(row-wise) with dense (UT, diag != 1)", Aq3, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(col-wise) with dense (LT, diag = 1)", Aq3t, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(col-wise) with dense (UT, diag = 1)", Aq3t, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(col-wise) with dense (LT, diag != 1)", Aq3t, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrmvConsistency (ctx, ctx, "hybrid(col-wise) with dense (UT, diag != 1)", Aq3t, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;

        SparseMatrix<bool> Aq2p (n,n);
	BLAS3::copy (ctx, Aq2, Aq2p);
	makeNonsingDiag(ctx.F, Aq2p, false);
        TransposeMatrix<SparseMatrix<bool> > Aq2tp (Aq2p);

        SparseMatrix<bool, Vector<GF2>::Hybrid> Aq3p (n,n);
	BLAS3::copy (ctx, Aq3, Aq3p);
        makeNonsingDiag(ctx.F, Aq3p, false);
        TransposeMatrix<SparseMatrix<bool, Vector<GF2>::Hybrid> > Aq3tp (Aq3p);

        if( !testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag = 1)", Aq2, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag = 1)", Aq2, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (LT, diag != 1)", Aq2p, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(row-wise) with dense (UT, diag != 1)", Aq2p, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag = 1)", Aq2t, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag = 1)", Aq2t, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (LT, diag != 1)", Aq2tp, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "sparse(col-wise) with dense (UT, diag != 1)", Aq2tp, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;

        if( !testtrsvConsistency (ctx, ctx, "hybrid(row-wise) with dense (LT, diag = 1)", Aq3, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(row-wise) with dense (UT, diag = 1)", Aq3, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(row-wise) with dense (LT, diag != 1)", Aq3p, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(row-wise) with dense (UT, diag != 1)", Aq3p, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(col-wise) with dense (LT, diag = 1)", Aq3t, v1d, Aq1, v1d, LowerTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(col-wise) with dense (UT, diag = 1)", Aq3t, v1d, Aq1, v1d, UpperTriangular, true)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(col-wise) with dense (LT, diag != 1)", Aq3tp, v1d, Aq1, v1d, LowerTriangular, false)) pass = false;
        if( !testtrsvConsistency (ctx, ctx, "hybrid(col-wise) with dense (UT, diag != 1)", Aq3tp, v1d, Aq1, v1d, UpperTriangular, false)) pass = false;

	if( !testgerConsistency (ctx, ctx, "sparse(row-wise)/dense/dense	with dense/dense/dense", A2, v2d, v1d, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "sparse(col-wise)/dense/dense	with dense/dense/dense", A2t, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "sparse(row-wise)/sparse/sparse	with dense/dense/dense", A2, v2s, v1s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "sparse(col-wise)/sparse/sparse	with dense/dense/dense", A2t, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "sparse(row-wise)/sparse/dense	with dense/dense/dense", A2, v2d, v1s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "sparse(col-wise)/sparse/dense	with dense/dense/dense", A2t, v1d, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "dense/sparse/dense			with dense/dense/dense", A1, v2d, v1s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "dense/sparse/sparse			with dense/dense/dense", A1, v2s, v1s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "dense/dense/sparse			with dense/dense/dense", A1, v2s, v1d, A1, v1d, v2d)) pass = false;

        if( !testgerConsistency (ctx, ctx, "hybrid(row-wise)/dense/dense        with dense/dense/dense", A3, v2d, v1d, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "hybrid(col-wise)/dense/dense        with dense/dense/dense", A3t, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "hybrid(row-wise)/sparse/sparse      with dense/dense/dense", A3, v2s, v1s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "hybrid(col-wise)/sparse/sparse      with dense/dense/dense", A3t, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "hybrid(row-wise)/sparse/dense       with dense/dense/dense", A3, v2d, v1s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "hybrid(col-wise)/sparse/dense       with dense/dense/dense", A3t, v1d, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "dense(col-wise)/sparse/dense        with dense/dense/dense", A1t, v1d, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "dense(col-wise)/sparse/sparse       with dense/dense/dense", A1t, v1s, v2s, A1, v1d, v2d)) pass = false;
        if( !testgerConsistency (ctx, ctx, "dense(col-wise)/dense/sparse        with dense/dense/dense", A1t, v1s, v2d, A1, v1d, v2d)) pass = false;


        commentator.stop (MSG_STATUS (pass));

        return pass;
}


int main (int argc, char **argv)
{
	bool pass = true;

	static long l = 374;
	static long m = 111;
	static long n = 412;
	static long p = 243;
	static long k = 100;
	static int iterations = 1;

	static Argument args[] = {
		{ 'l', "-l L", "Set row-dimension of matrix A to L.", TYPE_INT, &l },
		{ 'm', "-m M", "Set row-dimension of matrix B and column-dimension of A to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set row-dimension of matrix C and column-dimension of B to N.", TYPE_INT, &n },
		{ 'p', "-p P", "Set column-dimension of matrix C to P.", TYPE_INT, &p },
		{ 'k', "-k K", "K nonzero elements per row/column in sparse matrices.", TYPE_INT, &k },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Field F;
	Context<GF2, GenericModule<GF2> > ctx (F);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("BLAS GenericModule, GF2-specialisation, test suite", "BLASGenericModuleGF2");

	if (!testBLAS1 (ctx, "GF2", l, iterations)) pass = false;

	RandomDenseStream<Field, Vector<Field>::Dense> stream_v1 (F, l, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v2 (F, m, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v3 (F, n, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v4 (F, p, 1);

	Vector<Field>::Dense v1 (l), v2 (m), v3 (n), v4 (p);
	stream_v1 >> v1;
	stream_v2 >> v2;
	stream_v3 >> v3;
	stream_v4 >> v4;

	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream11 (F, m, l);
	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream12 (F, n, m);
	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream13 (F, p, n);
	RandomDenseStream<Field, DenseMatrix<Field::Element>::Row> stream14 (F, m, m);

	DenseMatrix<Field::Element> M1 (stream11);
	DenseMatrix<Field::Element> M2 (stream12);
	DenseMatrix<Field::Element> M3 (stream13);
	DenseMatrix<Field::Element> M4 (stream14);

	if (!testBLAS2 (ctx, "dense", M1, M2, v1, v2,
			DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "dense", M1, M2, M3, M4,
			DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;

	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream21 (F, (double) k / (double) m, m, l);
	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream22 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream23 (F, (double) k / (double) p, p, n);
	RandomSparseStream<Field, Vector<Field>::Sparse, Field::RandIter> stream24 (F, (double) k / (double) p, m, m);

	SparseMatrix<Field::Element, Vector<Field>::Sparse> M5 (stream21);
	SparseMatrix<Field::Element, Vector<Field>::Sparse> M6 (stream22);
	SparseMatrix<Field::Element, Vector<Field>::Sparse> M7 (stream23);
	SparseMatrix<Field::Element, Vector<Field>::Sparse> M8 (stream24);

	if (!testBLAS2 (ctx, "sparse", M5, M6, v1, v2,
			SparseMatrix<Field::Element, Vector<Field>::Sparse>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "sparse", M5, M6, M7, M8,
			SparseMatrix<Field::Element, Vector<Field>::Sparse>::IteratorType ()))
		pass = false;

	RandomHybridStream<Field, Vector<Field>::Hybrid> stream31 (F, (double) k / (double) m, m, l);
	RandomHybridStream<Field, Vector<Field>::Hybrid> stream32 (F, (double) k / (double) n, n, m);
	RandomHybridStream<Field, Vector<Field>::Hybrid> stream33 (F, (double) k / (double) p, p, n);
	RandomHybridStream<Field, Vector<Field>::Hybrid> stream34 (F, (double) k / (double) p, m, m);

	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M9 (stream31);
	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M10 (stream32);
	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M11 (stream33);
	SparseMatrix<Field::Element, Vector<Field>::Hybrid> M12 (stream34);

	if (!testBLAS2 (ctx, "hybrid", M9, M10, v1, v2,
			SparseMatrix<Field::Element, Vector<Field>::Hybrid>::IteratorType ()))
		pass = false;
	if (!testBLAS3 (ctx, "hybrid", M9, M10, M11, M12,
			SparseMatrix<Field::Element, Vector<Field>::Hybrid>::IteratorType ()))
		pass = false;

	RandomDenseStream<Field, Vector<Field>::Dense> stream_v1p (F, l - 8, 1);
	RandomDenseStream<Field, Vector<Field>::Dense> stream_v2p (F, m - 8, 1);

	Vector<Field>::Dense v1p (l - 8), v2p (m - 8);
	stream_v1p >> v1p;
	stream_v2p >> v2p;

	DenseMatrix<Field::Element>::SubmatrixType M13 (M1, 4, 4, l - 8, m - 8);
	DenseMatrix<Field::Element>::SubmatrixType M14 (M2, 4, 4, m - 8, n - 8);
	DenseMatrix<Field::Element>::SubmatrixType M15 (M3, 4, 4, n - 8, p - 8);
	DenseMatrix<Field::Element>::SubmatrixType M16 (M4, 4, 4, m - 8, m - 8);

	if (!testBLAS2Submatrix (ctx, "dense (submatrix)", M13, M14, v1p, v2p,
				 DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;
	if (!testBLAS3Submatrix (ctx, "dense (submatrix)", M13, M14, M15, M16,
				 DenseMatrix<Field::Element>::IteratorType ()))
		pass = false;

	if (!testBLAS2RepsConsistency(ctx, "BLAS2 Consistency", n, m, k))
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
