/* tests/test-modular-double.C
 * Brazenly stolen by bds from 
 * tests/test-modular-short.C
 * Brazenly stolen by Zhendong Wan (Copyright 2003) from 
 * tests/test-modular.C
 * Copyright 2001, 2002 Bradford Hovinen,
 * Copyright 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>,
 *            Dave Saunders <saunders@cis.udel.edu>
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

#include "lela/ring/modular.h"

#include "test-common.h"
#include "test-ring.h"
#include "test-blas-level1.h"

using namespace LELA;

int main (int argc, char **argv)
{
	static int iterations = 1;

	static Argument args[] = {
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("Modular<double> ring test suite", "Modular<double>");
	bool pass = true;

	//Modular<double> F2 (2); 
	Modular<double> F3 (3); 
	Modular<double> F5 (5); 
	Modular<double> F7 (7); 
	Modular<double> F11 (11); 
	Modular<double> F (32749); 
	Modular<double> G (65521); 
	//Modular<double> H (1099511627689); 

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	std::ostream& report = commentator.report();
	report << "Ring F2" << std::endl;
	//if (!runRingTests (F2,  "Modular<double>",  iterations, false)) pass = false;
	report << "Ring F3" << std::endl;
	if (!runRingTests (F3,  "Modular<double>",  iterations, false)) pass = false;
	report << "Ring F5" << std::endl;
	if (!runRingTests (F5,  "Modular<double>",  iterations, false)) pass = false;
	report << "Ring F7" << std::endl;
	if (!runRingTests (F7,  "Modular<double>",  iterations, false)) pass = false;
	report << "Ring F11" << std::endl;
	if (!runRingTests (F11,  "Modular<double>",  iterations, false)) pass = false;
	report << "Ring F" << std::endl;
	if (!runRingTests (F,  "Modular<double>",  iterations, false)) pass = false;
	report << "Ring G" << std::endl;
	if (!runRingTests (G,  "Modular<double>",  iterations, false)) pass = false;

	commentator.stop(MSG_STATUS (pass));
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
