/* tests/test-vector-domain.C
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"
#include "linbox/ring/gf2.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/sparse.h"

#include "test-blas-level1.h"

using namespace LinBox;

template <class Field>
bool testVectorDomain (Field &F, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing VectorDomain <" << text << ">" << std::ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, typename Vector<Field>::Dense> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomSparseStream<Field, typename Vector<Field>::Sparse> stream3 (F, 0.1, n, iterations), stream4 (F, 0.1, n, iterations);

	Context<Field> ctx (F);

	if (!testCopyEqual (ctx, "dense/dense", stream1, stream1)) pass = false;
	if (!testCopyEqual (ctx, "dense/sparse", stream1, stream3)) pass = false;
	if (!testCopyEqual (ctx, "sparse/dense", stream3, stream1)) pass = false;
	if (!testCopyEqual (ctx, "sparse/sparse", stream3, stream3)) pass = false;

	if (!testDotProduct (ctx, "dense/dense", stream1, stream2)) pass = false;
	if (!testDotProduct (ctx, "sparse/dense", stream3, stream1)) pass = false;
	if (!testDotProduct (ctx, "sparse/sparse", stream3, stream4)) pass = false;

	if (!testScal (ctx, "dense", stream1)) pass = false;
	if (!testScal (ctx, "sparse", stream3)) pass = false;

	if (!testAXPY (ctx, "dense", stream1, stream2)) pass = false;
	if (!testAXPY (ctx, "sparse", stream3, stream4)) pass = false;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

bool testVectorDomain (GF2 &F, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing VectorDomain <" << text << ">" << std::ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<GF2, Vector<GF2>::Dense> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomSparseStream<GF2, Vector<GF2>::Sparse> stream3 (F, 0.1, n, iterations), stream4 (F, 0.1, n, iterations);
	RandomSparseStream<GF2, Vector<GF2>::Hybrid> stream5 (F, 0.1, n, iterations), stream6 (F, 0.1, n, iterations);

	Context<GF2> ctx (F);

#if 0
	if (!testCopyEqual (ctx, "dense/dense", stream1, stream2)) pass = false;
	if (!testCopyEqual (ctx, "dense/sparse", stream1, stream4)) pass = false;
	if (!testCopyEqual (ctx, "dense/hybrid", stream1, stream6)) pass = false;
	if (!testCopyEqual (ctx, "sparse/dense", stream3, stream2)) pass = false;
	if (!testCopyEqual (ctx, "sparse/sparse", stream3, stream4)) pass = false;
	if (!testCopyEqual (ctx, "sparse/hybrid", stream3, stream6)) pass = false;
	if (!testCopyEqual (ctx, "hybrid/dense", stream5, stream2)) pass = false;
	if (!testCopyEqual (ctx, "hybrid/sparse", stream5, stream4)) pass = false;
	if (!testCopyEqual (ctx, "hybrid/hybrid", stream5, stream6)) pass = false;

	if (!runDPTests (ctx, "dense/dense", stream1, stream2)) pass = false;
	if (!runDPTests (ctx, "dense/sparse", stream1, stream3)) pass = false;
	if (!runDPTests (ctx, "dense/hybrid", stream1, stream5)) pass = false;
	if (!runDPTests (ctx, "sparse/sparse", stream3, stream4)) pass = false;

	if (!testScal (ctx, "dense", stream1, stream2)) pass = false;
	if (!testScal (ctx, "sparse", stream3, stream4)) pass = false;
	if (!testScal (ctx, "hybrid", stream5, stream6)) pass = false;

	if (!testAXPY (ctx, "dense", stream1, stream2)) pass = false;
	if (!testAXPY (ctx, "sparse", stream3, stream4)) pass = false;
	if (!testAXPY (ctx, "hybrid", stream5, stream6)) pass = false;
#endif

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long n = 100;
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static int q4 = 101;
	static int iterations = 2;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to N.", TYPE_INT,     &n },
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus.", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16 modulus.", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"field\" GF(P) [1] for uint8 modulus.", TYPE_INTEGER, &q4 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Modular<integer> F_integer (q1);
	Modular<uint32> F_uint32 (q2.get_ui ());
	Modular<uint16> F_uint16 (q3.get_ui ());
	Modular<uint8> F_uint8 ((uint8) q4);
	GF2 F_2;

	commentator.start ("Vector domain test suite", "VectorDomain");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	if (!testVectorDomain (F_integer, "Modular <integer>", n, iterations)) pass = false;
	if (!testVectorDomain (F_uint32, "Modular <uint32>", n, iterations)) pass = false;
	if (!testVectorDomain (F_uint16, "Modular <uint16>", n, iterations)) pass = false;
	if (!testVectorDomain (F_uint8, "Modular <uint8>", n, iterations)) pass = false;
	if (!testVectorDomain (F_2, "GF2", n, iterations)) pass = false;

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
