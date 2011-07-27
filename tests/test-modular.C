/* tests/test-modular.C
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
using namespace std;

/* Random number test
 *
 * Test that the random iterator over the given ring works
 */

template <class Element>
bool testRandomIteratorStep (const Modular<Element> &F,
			     const char *text,
			     unsigned int num_trials,
			     unsigned int num_categories,
			     unsigned int hist_len) 
{
	//std::ostringstream str;

	//str << "Testing " << text << "::RandIter" << std::ends;

	//LELA::commentator.start (str.str ().c_str (), "testRandomIteratorStep");
	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	bool ret = true;

	LELA::integer card;
	unsigned int i;
	std::vector<int> categories1 (num_categories, 0);
	std::vector<int> categories2 (num_categories, 0);
	std::list<std::vector<int> > diff_categories;
	std::list<Element> x_queue;

	F.cardinality (card);

	typename Modular<Element>::RandIter iter (F);
	Element x,  d;

	std::list<std::vector<int> >::iterator diff_cat_iter;

	for (i = 0; i < hist_len; ++i)
		diff_categories.push_back (std::vector<int> (num_categories, 0));

	// I make the simplifying assumption that ring elements are actually
	// C++ ints. Otherwise, I don't know how to place the numbers into
	// categories in any well-defined manner.
	for (i = 0; i < num_trials; ++i) {
		iter.random (x);
		integer ix (x), ixmodn;
		ixmodn = ix % num_categories;
		categories1[ixmodn.get_ui ()]++;
		categories2[(unsigned int) (ix.get_d () / card.get_d () * num_categories)]++;

		typename std::list<Element>::iterator x_queue_iter = x_queue.begin ();
		diff_cat_iter = diff_categories.begin ();

		for (; x_queue_iter != x_queue.end (); ++x_queue_iter, ++diff_cat_iter) {
			integer id (F.sub (d, *x_queue_iter, x));
			ixmodn = id % num_categories;
			(*diff_cat_iter)[ixmodn.get_ui ()]++;
		}

		x_queue.push_front (x);

		while (x_queue.size () > hist_len)
			x_queue.pop_back ();
	}

	double p, chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories1[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * (double)num_categories / (double)num_trials, (double)num_categories - 1.0);

	report << "Test of distribution uniformity (low-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (low-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) 
		reportError("Random iterator's values do not appear to be uniformly distributed", ret);

	chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories2[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

	report << "Test of distribution uniformity (high-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (high-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) 
		reportError("Consistency failure for addition", ret);

	diff_cat_iter = diff_categories.begin ();

	int idx = 0;

	for (; diff_cat_iter != diff_categories.end (); ++diff_cat_iter, ++idx) {
		chi_squared = 0.0;

		for (i = 0; i < num_categories; ++i)
			chi_squared += pow (double ((*diff_cat_iter)[i]) -
					    double (num_trials) / double (num_categories), 2);

		p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

		report << "Test of " << idx + 1 << " spacing: chi^2 = "
		       << chi_squared * num_categories / num_trials << std::endl;
		report << "Test of " << idx + 1 << " spacing:     p = " << p << std::endl;

		if (p < 0.05 || p > 0.95) 
			reportError("Difference values do not appear to be uniformly distributed", ret);
	}

	//LELA::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomIteratorStep");
	return ret;
}

/** Random number test
 *
 * Test that the random iterator over the given ring works.
 *
 * Test up to five times, accepting any one, to increase probability of 
 * passing statistical tests.
 */
template <class Element>
bool testRandomIterator (const Modular<Element> &F, const char *text,
			 unsigned int num_trials,
			 unsigned int num_categories,
			 unsigned int hist_len) 
{
	std::ostringstream str;

	str << "Testing " << text << "::RandIter" << std::ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());

	LELA::commentator.start (st, "testRandomIterator");

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	/* This test either passes or runs a lot of times */
	for (int i = 1; !testRandomIteratorStep (F, text, num_trials, num_categories, hist_len) && i < 20; ++i)
		if (0 == i % 10)  
			report << "Warning! Probable failure of uniformity" << std::endl;

	LELA::commentator.stop (MSG_STATUS (true), (const char *) 0, "testRandomIterator");

	delete[] st;
	return true;

}

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static int q4 = 101;
	static int iterations = 1;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"ring\" GF(Q) [1] for integer modulus.", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"ring\" GF(Q) [1] for uint32 modulus.", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"ring\" GF(Q) [1] for uint16 modulus.", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"ring\" GF(Q) [1] for uint8 modulus.", TYPE_INT, &q4 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("Modular test suite", "Modular ");
	bool pass = true;

	Modular<integer> F_integer (q1);
	Modular<uint32> F_uint32 (q2.get_ui ());
	Modular<uint16> F_uint16 (q3.get_ui ());
	Modular<uint8> F_uint8 ((uint8) q4);
	Modular<float> F_float ((float) q4);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (6);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runRingTests (F_integer, "Modular<integer>", iterations, false)) pass = false;
	if (!runRingTests (F_uint32,  "Modular<uint32>",  iterations, false)) pass = false;
	if (!runRingTests (F_uint16,  "Modular<uint16>",  iterations, false)) pass = false;
	if (!runRingTests (F_uint8,  "Modular<uint8>",  iterations, false)) pass = false;
	if (!runRingTests (F_float,  "Modular<float>",  iterations, false)) pass = false;

	//if (!testRandomIterator (F_integer, "Modular<integer>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint32,  "Modular<uint32>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint16,  "Modular<uint16>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint8,  "Modular<uint8>", trials, categories, hist_level)) pass = false;

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
