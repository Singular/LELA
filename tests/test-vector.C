/* lela/tests/test-vector.C
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic tests for vectors
 * 
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include <vector>

#include "lela/integer.h"
#include "lela/ring/modular.h"
#include "lela/ring/gf2.h"
#include "lela/vector/sparse.h"

#include "test-common.h"
#include "test-vector.h"

using namespace LELA;

int main (int argc, char **argv)
{
	bool pass = true;

	static integer q = 101U;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the ring Z/Q for uint32 modulus.", TYPE_INTEGER, &q },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint32> Ring;

	Ring R (q);
	GF2 gf2;

	commentator.start ("Vector-test-suite", "Vector");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	pass = testVector<Vector<Ring>::Dense> (R) && pass;
	pass = testVector<Vector<Ring>::Sparse> (R) && pass;
	pass = testVector<Vector<GF2>::Dense> (gf2) && pass;
	pass = testVector<Vector<GF2>::Sparse> (gf2) && pass;

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
