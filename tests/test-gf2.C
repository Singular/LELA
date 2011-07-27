/* tests/test-gf2.C
 * Copyright 2003 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Tests for GF(2) field
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/lela-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "lela/ring/gf2.h"
#include "lela/ring/modular.h"
#include "lela/blas/context.h"

#include "test-ring.h"
#include "test-blas-level1.h"

using namespace LELA;
using LELA::uint16;
using LELA::uint32;

int main (int argc, char **argv)
{
	static unsigned int n = 1000;
	static int iterations = 2;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start ("GF2 field test suite", "GF2");
	bool pass = true;

	GF2 F;

	RandomDenseStream<GF2, Vector<GF2>::Dense> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomSparseStream<GF2, Vector<GF2>::Sparse, GF2RandIter, VectorRepresentationTypes::Sparse01>
		stream3 (F, 0.1, n, iterations),
		stream4 (F, 0.1, n, iterations);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	pass = runBasicRingTests (F, "GF(2)", iterations);

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
