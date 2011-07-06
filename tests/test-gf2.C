/* tests/test-gf2.C
 * Copyright (C) 2003 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to Modular interace.
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/ring/gf2.h"
#include "linbox/ring/modular.h"
#include "linbox/blas/context.h"

#include "test-generic-for-quad.h"
#include "test-blas-level1.h"

using namespace LinBox;
using LinBox::uint16;
using LinBox::uint32;

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

	commentator.start ("Testing GF2", "main", 10);
	
	if (!testFieldNegation         (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldInversion        (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldAxioms           (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldAssociativity    (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldCharacteristic   (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testGeometricSummation    (F, "GF2", iterations, 100)) pass = false; commentator.progress ();
	if (!testFreshmansDream        (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testArithmeticConsistency (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testAxpyConsistency       (F, "GF2", iterations))      pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "main");

	if (!testRandomIterator (F, "GF2", trials, categories, hist_level)) pass = false;

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
