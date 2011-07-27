/* lela/tests/test-common.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_COMMON_H
#define __LELA_TESTS_TEST_COMMON_H

#include <iostream>
#include <fstream>
#include <vector>

#include "lela/ring/interface.h"
#include "lela/integer.h"
#include "lela/randiter/nonzero.h"

using namespace std;

enum ArgumentType {
	TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE, TYPE_STRING
};
#define TYPE_BOOL TYPE_NONE

struct Argument 
{
	char             c;
	const char      *example;
	const char      *helpString;
	ArgumentType     type;
	void            *data;
};
// example may be passed as null and will be generated intelligently
// eg "-b {YN+-}" for bools, "-v v" for all else

/* Force the input-matrix to have nonsingular entries on the diagonal */

template <class Field, class Matrix>
static Matrix &makeNonsingDiag (Field &F, Matrix &A, bool set_one)
{
	size_t i;

	typename Field::Element a;

	LELA::NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	for (i = 0; i < std::min (A.rowdim (), A.coldim ()); ++i) {
		if (set_one)
			A.setEntry (i, i, F.one ());
		else if (!A.getEntry (a, i, i) || F.isZero (a)) {
			r.random (a);
			A.setEntry (i, i, a);
		}
	}

	return A;
}

template <class Ring, class Polynomial>
void printPolynomial (Ring &F, ostream &output, const Polynomial &v) 
{
	int i;
	size_t val;

	for (val = 0; val < v.size () && F.isZero (v[val]); val++) ;

	if (v.size () == 0 || val == v.size ())
		output << "0";

	for (i = v.size () - 1; i >= 0; i--) {
		if (F.isZero (v[i]))
			continue;

		if (!F.isOne (v[i]) || i == 0)
			F.write (output, v[i]);

		if (i > 0)
			output << " x^" << i;

		if (i > (int) val)
			output << " + ";
	}

	output << endl;
}

void parseArguments (int argc, char **argv, Argument *args, bool printDefaults = true);
void printHelpMessage (const char *program, Argument *args, bool printDefaults = false);

/** writes the values of all arguments, preceded by the programName */
std::ostream& writeCommandString (std::ostream& os, Argument *args, char *programName);

bool isPower (LELA::integer n, LELA::integer m);

/* Give an approximation of the value of the incomplete gamma function at a, x,
 * to within the tolerance tol */

extern inline double incompleteGamma (double a, double x, double tol);

/* Give the value of the chi-squared cumulative density function for given
 * value of chi_sqr and the given degrees of freedom */

double chiSquaredCDF (double chi_sqr, double df);

#endif // __LELA_TESTS_TEST_COMMON_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
