/* lela/tests/test-ring.h
 * Copyright 2001, 2002 Bradford Hovinen
 * 
 * Generic tests for rings
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_RING_H
#define __LELA_TESTS_TEST_RING_H

#include <iostream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cmath>

#include "lela/util/commentator.h"
#include "lela/integer.h"

#include "test-common.h" 

/* Modular exponentiation */
using namespace std;
using namespace LELA;

template <class Ring>
typename Ring::Element &expt (const Ring &F, typename Ring::Element &res, const typename Ring::Element &a, LELA::integer &n) 
{
	if (n == 0) {
		F.init (res, 1);
	}
	else if (n == 1) {
		F.copy (res, a);
	}
	else if (n.get_ui () & 1) {
		n -= 1;
		expt (F, res, a, n);
		F.mulin (res, a);
	} else {
		n /= 2;
		expt (F, res, a, n);
		F.mulin (res, res);
	}

	return res;
}

bool reportError(string rep, bool& flag)
{
	ostream &report = commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
	report << "ERROR: " << rep << endl;
	return flag = false;
}

/** Check each ring or ring operation.
 *
 * Test various ring operations
 *
 * F - Ring over which to perform computations
 * title - String to use as the descriptive title of this test
 * ringp - use true if inv and div must work for all nonzero denominators
 *
 * Return true on success and false on failure
 */

template<class Ring>
bool testRing (Ring &F, const char *title, bool ringp = true) 
{
	commentator.start (title, "testRing", 5);
	ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	typename Ring::Element zero, one, two, three;
	F.init (zero, 0); F.init (one, 1); F.init (two, 2); F.init (three, 3);

	typename Ring::Element a, b, c, d, e, f;
	F.init (a, 0); F.init (b, 0); F.init (c, 0); F.init (d, 0); F.init (e, 0); F.init (f, 0);

	report << "Ring self description: ";
	F.write (report) << endl;

	LELA::integer n, m;
	bool pass = true, part_pass = true;

	commentator.start ("Testing characteristic/cardinality match");

	F.characteristic (n); 
	F.cardinality (m);

	if (n > 0 && !isPower (m, n)) part_pass = reportError("Characteristic, cardinality mismatch", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* tests for presence of members with minimal check of semantics */
	// these checks need improvement 

	commentator.start ("Testing correctness of 0, 1, and -1");
	part_pass = true;

	if (!F.isZero (zero))
		part_pass = reportError ("isZero (0) is false", pass);

	if (F.isZero (one))
		part_pass = reportError ("isZero (1) is true", pass);

	if (F.isOne (zero))
		part_pass = reportError ("isOne (0) is true", pass);

	if (!F.isOne (one))
		part_pass = reportError ("isOne (1) is false", pass);

	F.add (a, F.one (), F.minusOne ());

	if (!F.isZero (a))
		part_pass = reportError ("1 + -1 != 0", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing init");
	part_pass = true;

	int pos_cases[] = { 1, 10, 20, 97 };
	int neg_cases[] = { -1, -10, -20, -97 };

	// Test zero
	F.init (a, 0);

	report << "Result of F.init (a, 0): ";
	F.write (report, a);

	if (!F.isZero (a))
		part_pass = reportError ("Ring reports F.init (a, 0) != 0", pass);

	size_t j;
	int v;
	
	for (j = 0; j < sizeof (pos_cases) / sizeof (int); ++j) {
		F.init (a, pos_cases[j]);
		F.copy (b, F.zero ());

		for (v = 0; v < pos_cases[j]; ++v)
			F.addin (b, F.one ());

		if (!F.areEqual (a, b)) {
			std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "Ring reports F.init (a, " << pos_cases[j] << ") != 1 + ... + 1 (" << pos_cases[j] << " summands)" << std::endl
			      << "Reported value of a = ";
			F.write (error, a) << std::endl
					   << "Reported value of sum = ";
			F.write (error, b) << std::endl;
			part_pass = false;
		}
	}

	for (j = 0; j < sizeof (neg_cases) / sizeof (int); ++j) {
		F.init (a, neg_cases[j]);
		F.copy (b, F.zero ());

		for (v = 0; v < -neg_cases[j]; ++v)
			F.addin (b, F.minusOne ());

		if (!F.areEqual (a, b)) {
			std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "Ring reports F.init (a, " << neg_cases[j] << ") != -1 + ... + -1 (" << neg_cases[j] << " summands)" << std::endl
			      << "Reported value of a = ";
			F.write (error, a) << std::endl
					   << "Reported value of sum = ";
			F.write (error, b) << std::endl;
			part_pass = false;
		}
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing ring arithmetic");
	part_pass = true;

	F.init (b, n-2);
	F.init (d, n-2);
	F.init (e, 3);

	F.add (a, three, two);
	F.write (report << "Result of 2 + 3: ", a) << endl;
	F.copy (d, three);
	F.addin (d, two);
	F.write (report << "Result of 2 + 3 (inplace): ", d) << endl;

	if (!F.areEqual (a, F.init (f, 5)) || !F.areEqual (d, a)) 
		part_pass = reportError ("Results of add are incorrect", pass);

	F.neg (a, two);
	F.write (report << "Result of -2: ", a) << endl;
	F.copy (d, two);
	F.negin (d);
	F.write (report << "Result of -2 (inplace): ", d) << endl;
	F.init (f, -2);
	F.write (report << "-2 via init: ", f) << endl;

	if (!F.areEqual (a, f) || !F.areEqual (d, a)) 
		part_pass = reportError ("Results of neg are incorrect", pass);

	F.sub (a, three, two);
	F.write (report << "Result of 3 - 2: ", a) << endl;
	F.init (d, 3);
	F.subin (d, two);
	F.write (report << "Result of 3 - 2 (inplace): ", d) << endl;

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) 
		part_pass = reportError ("Results of neg sub incorrect", pass);

	F.mul (a, two, three);
	F.write (report << "Result of 2 * 3: ", a) << endl;
	F.copy (d, two);
	F.mulin (d, three);
	F.write (report << "Result of 2 * 3 (inplace): ", d) << endl;
	F.init (f, 6);

	if (!F.areEqual (a, f) || !F.areEqual (d, a)) 
		part_pass = reportError ("Results of mul incorrect", pass);

	F.inv (a, one);
	F.write (report << "Result of inverting 1: ", a) << endl;
	F.copy (d, one);
	F.invin (d);
	F.write (report << "Result of inverting 1 (inplace): ", d) << endl;

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) 
		part_pass = reportError ("Results of inv incorrect", pass);

	if (F.isZero (two))
		F.copy(f, three);
	else
		F.copy(f, two);

	F.div (a, f, f);
	F.write ( report << "Result of f/f: ", a) << endl;
	F.copy(d, f);
	F.divin (d, f);
	F.write ( report << "Result of f/f (inplace): ", d) << endl;

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of div incorrect", pass);

	F.axpy (a, two, three, two);
	F.write (report << "Result of axpy 2*3 + 2: ", a) << endl;
	F.copy (d, two);
	F.axpyin (d, two, three);
	F.write (report << "Result of axpy 2*3 + 2 (inplace): ", d) << endl;
	F.init (f, 8);

	if ( !F.areEqual (a, f) || !F.areEqual (d, a) ) 
		part_pass = reportError( "Results of axpy incorrect", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing summation of powers of 2");
	part_pass = true;

	//,..
	// 2^101 - 1 vs 1 + 2 + 4 + ... + 2^100

	F.init (a, 1);
	F.init (b, 1);
	F.init (c, 0);
	
	n = 101;
	expt (F, a, two, n);
	F.subin (a, one);
	F.write (report << "using expt, 2^101 - 1: ", a) << endl;

	for (int i = 1; i <= 101; ++i) {
		F.addin (c, b);
		F.mulin (b, two);
	}

	F.write (report << "using mul-add, 2^101 - 1: ", c) << endl;

	if (!F.areEqual (a, c)) 
		part_pass = reportError( "2^101 - 1 != 1 + 2 + .. + 2^100", pass);

	// 1 + 2*(1 + 2*( ... (1) ... )), 100 times.
	F.copy (d, one);
	for (int i = 1; i < 101; ++i) {
		F.axpy (b, two, d, one);
		F.copy (d, b);
	}

	F.write( report << "using axpy, 1 + 2(1 + ... + 2(1)...), with 100 '+'s: ", d) << endl;

	if (!F.areEqual (a, d)) 
		part_pass = reportError( "2^101 - 1 != 1 + 2(1 + ... + 2(1)...), with 100 '+'s: ", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* untested so far
	   ostream &write (ostream &os) const 
	   istream &read (istream &is)
	   ostream &write (ostream &os, const Element &x) const 
	   istream &read (istream &is, Element &x) const
	   RingArchetype (RingAbstract*, ElementAbstract*, RandIterAbstract* = 0)
	*/

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "testRing");
	
	return pass;
}


/** Tests of algebraic properties of rings */

/* Generic test 6: Negation of elements
 *
 * Negates random elements and checks that they are true negatives
 */

template <class Ring>
bool testRingNegation (const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " negation" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingNegation", iterations);

	typename Ring::Element a, neg_a, neg_a_a, zero;
	F.init (a,0); F.init (neg_a,0); F.init (neg_a_a,0); F.init (zero, 0);
	typename Ring::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);
		
		r.random (a);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.neg (neg_a, a);

		report << "-a = ";
		F.write (report, neg_a) << endl;

		F.add (neg_a_a, neg_a, a);

		report << "a + -a = ";
		F.write (report, neg_a_a) << endl;

		if (!F.areEqual (neg_a_a, zero)) reportError("a + -a != 0" , ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingNegation");
	delete[] st;
	return ret;
}

/** Generic test 5: Inversion of elements
 *
 * Inverts random elements and checks that they are true inverses
 */

template <class Ring>
bool testRingInversion (const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " inversion" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingInversion", iterations);

	typename Ring::Element a, ainv, aainv, one;
        F.init (a,0); F.init (ainv,0); F.init (aainv,0);
	F.init (one, 1);
	typename Ring::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		if (!F.inv (ainv, a)) {
			report << "Ring reports a^{-1} does not exist" << endl;
		} else {
			report << "a^{-1} = ";  F.write (report, ainv) << endl;

			F.mul (aainv, ainv, a);

			report << "a a^{-1} = ";  F.write (report, aainv) << endl;

			if (!F.areEqual (aainv, one)) reportError("a a^-1 != 1", ret);
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingInversion");
	delete[] st;
	return ret;
}

/** @brief Generic test 7a: Distributivity of multiplication over addition

 * Given random ring elements 'a', 'b', and 'c', checks that
 * (a + b) * c = a * c + b * c  and  c * (a + b) = c * a + c * b
 */


template <class Ring>
bool testRingDistributivity(const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " distributivity" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingDistributivity", iterations);

	typename Ring::Element a, b, c, a_b, a_bc, ac, bc, ac_bc, ca_b, ca, cb, ca_cb;
        F.init (a,0); F.init (b,0); F.init (c,0); 
        F.init (a_b,0); F.init (a_bc,0); F.init (ac,0); F.init (bc,0);
        F.init (ac_bc,0); 
		F.init (ca_b,0); F.init (ca,0); F.init (cb,0); F.init (ca_cb,0); 
        
	typename Ring::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);
		r.random (c);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << ", c = ";
		F.write (report, c) << endl;

		F.add (a_b, a, b);//a + b
		F.mul (a_bc, a_b, c);//(a+b)*c
		F.mul (ca_b, c, a_b);//c*(a+b)
		F.mul (ac, a, c); //a*c
		F.mul (bc, b, c); //b*c
		F.mul (ca, c, a); //c*a
		F.mul (cb, c, b); //c*b
		F.add (ac_bc, ac, bc); //a*c + b*c
		F.add (ca_cb, ca, cb); //c*a + c*b

		report << "(a + b) * c = ";
		F.write (report, a_bc) << endl;

		report << "a * c + b * c = ";
		F.write (report, ac_bc) << endl;

		report << "c * (a + b) = ";
		F.write (report, ca_b) << endl;

		report << "c * a + c * b = ";
		F.write (report, ca_cb) << endl;
		if (!F.areEqual (a_bc, ac_bc) || !F.areEqual (ca_b, ca_cb)) 
			reportError("Operations were not distributative", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingDistributivity");
	delete[] st;
	return ret;
}


/** @brief Generic test 7b: Commutativity of multiplication and addition

 * Given random ring elements 'a', 'b', checks that
 * a*b = b*a
 * a+b = b+a
 */


template <class Ring>
bool testRingCommutativity (const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " commutativity," << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingCommutativity", iterations);

	typename Ring::Element a, b, ab, ba, a_b, b_a;
        F.init (a,0); F.init (b,0);
        F.init (ab,0); F.init (ba,0);
        F.init (a_b,0); F.init (b_a,0);

        
	typename Ring::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;

		F.mul (ab, a, b);
		F.mul (ba, b, a);

		report << "a*b = ";
		F.write (report, ab) << endl;

		report << "b*a = ";
		F.write (report, ba) << endl;

		if (!F.areEqual (ab, ba)) 
			reportError("Multiplication was not commutative", ret);

		F.add (a_b, a, b);
		F.add (b_a, b, a);

		report << "a+b = ";
		F.write (report, a_b) << endl;

		report << "b+a = ";
		F.write (report, b_a) << endl;

		if (!F.areEqual (a_b, b_a)) 
			reportError("Addition was not commutative", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingCommutativity");
	delete[] st;
	return ret;
}


/** Generic test 7c: Associativity of addition and multiplication
 *
 * Given random ring elements 'a', 'b', and 'c', checks that
 * (a * b) * c = a * (b * c) and (a + b) + c = a + (b + c)
 */

template <class Ring>
bool testRingAssociativity (const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " associativity" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingAssociativity", iterations);

	typename Ring::Element a, b, c, a_b, b_c, a_bc, ab_c;
        F.init (a,0); F.init (b,0); F.init (c,0);
        F.init (a_b,0); F.init (b_c,0); F.init (a_bc,0); F.init (ab_c,0);
	typename Ring::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);
		r.random (c);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << ", c = ";
		F.write (report, c) << endl;

		F.add (a_b, a, b);
		F.add (ab_c, a_b, c);
		F.add (b_c, b, c);
		F.add (a_bc, a, b_c);

		report << "(a + b) + c = ";
		F.write (report, ab_c) << endl;

		report << "a + (b + c) = ";
		F.write (report, a_bc) << endl;

		if (!F.areEqual (ab_c, a_bc)) reportError( "Results are not equal", ret);

		F.mul (a_b, a, b);
		F.mul (ab_c, a_b, c);
		F.mul (b_c, b, c);
		F.mul (a_bc, a, b_c);

		report << "(a * b) * c = ";
		F.write (report, ab_c) << endl;

		report << "a * (b * c) = ";
		F.write (report, a_bc) << endl;

		if (!F.areEqual (ab_c, a_bc)) reportError( "Results are not equal", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingAssociativity");
	delete[] st;
	return ret;
}

/** Generic test 2: Geometric summation
 *
 * Generates a random ring element 'a' and raises it through repeated
 * exponentiation to the power n. Takes the sum k of all intermediate values and
 * checks that a^n = (k-1)/(a-1).
 */

template <class Ring>
bool testGeometricSummation (const Ring &F, const char *name, unsigned int iterations, unsigned int n) 
{
	std::ostringstream str;
	str << "Testing " << name << " geometric summation" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testGeometricSummation", iterations);

	typename Ring::Element a, a_n, k, zero, one;
	typename Ring::RandIter r (F);

	F.init (zero, 0);
	F.init (one, 1);
        F.init (a,0); F.init (a_n,0); F.init (k,0);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		do r.random (a); while (F.areEqual (a, one));

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.copy (k, one);
		F.copy (a_n, a);

		for (unsigned int j = 0; j < n; ++j) {
			F.addin (k, a_n);
			F.mulin (a_n, a);
		}

		report << "n = " << n << " a^n = ";
		F. write (report, a_n) << endl;
		F. write(report);
		report<<std::endl;

		report << "sum(a^i, i = 0..n-1) = ";
		F.write (report, k) << endl;

		F.subin (a_n, one);
		report << "(a^n - 1) = ";
		F.write (report, a_n) << endl;

		F.subin (a, one);
		report << "(a - 1) = ";
		F.write (report, a) << endl;

		F.divin (a_n, a);

		report << "(a^n - 1) / (a - 1) = ";
		F.write (report, a_n) << endl;

		if (!F.areEqual (k, a_n)) reportError("Ring elements are not equal", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testGeometricSummation");
	delete[] st;
	return ret;
}

/** Generic test 3: Test of ring characteristic
 *
 * Take random ring elements and add them p times, where p is the
 * characteristic of the ring. Checks that the sum is 0. The test is not too
 * useful when the characteristic of the ring is 0, but it should still work
 * correctly.
 */

template <class Ring>
bool testRingCharacteristic (const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " characteristic" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (string(str.str()).c_str(), "testRingCharacteristic", iterations);

	LELA::integer p, j;
	typename Ring::Element a, sigma, zero;
	typename Ring::RandIter r (F);

	F.characteristic (p);
	F.init (zero, 0);
        F.init (a,0); F.init (sigma,0);
        
	bool ret = true;

	ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Ring characteristic: " << p << endl;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.copy (sigma, zero);

		for (j = 0; j < p; j += 1)
			F.addin (sigma, a);

		report << "p a = ";
		F.write (report, sigma) << endl;

		if (!F.isZero (sigma)) reportError("p a != 0", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingCharacteristic");
	delete[] st;
	return ret;
}

/** Generic test 4: The Freshman's Dream
 *
 * Generates two random ring elements 'a' and 'b', and checks whether
 * (a + b)^p = a^p + b^p, where p is the characteristic of the ring. Bails out
 * (returning true) if the ring is of characteristic 0.
 */

template <class Ring>
bool testFreshmansDream (const Ring &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " Freshman's Dream" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFreshmansDream", iterations);

	LELA::integer c, j;

	F.characteristic (c);

	if (c == 0) {
		commentator.stop ("skipping", "Ring characteristic is 0, so this test makes no sense", "testFreshmansDream");
		delete[] st;
		return true;
	}

	bool ret = true;

	typename Ring::RandIter r (F);
	typename Ring::Element a, b, a_b, a_b_p, a_p, b_p, a_p_b_p;

        F.init (a,0); F.init (b,0); F.init (a_b,0);
        F.init (a_b_p,0); F.init (a_p,0); F.init (b_p,0); F.init (a_p_b_p,0);

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;

		F.add (a_b, a, b);

		j = c; expt (F, a_b_p, a_b, j);
		j = c; expt (F, a_p, a, j);
		j = c; expt (F, b_p, b, j);

		F.add (a_p_b_p, a_p, b_p);

		report << "(a + b)^p = ";
		F.write (report, a_b_p);
		report << endl;

		report << "      a^p = ";
		F.write (report, a_p);
		report << endl;

		report << "      b^p = ";
		F.write (report, b_p);
		report << endl;

		report << "a^p + b^p = ";
		F.write (report, a_p_b_p);
		report << endl;

		if (!F.areEqual (a_b_p, a_p_b_p)) reportError("(a + b)^p != a^p + b^p", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFreshmansDream");
	delete[]  st;
	return ret;
}


/* Tests of ring features */ 

/** Generic test 7: Consistency of in-place and out-of-place arithmetic
 *
 * Generates random elements 'a' and 'b' and performs all basic arithmetic
 * operations in-place and out-of-place, checking for consistency.
 *
 * Div and inv are checked in a separate function.
 */

template <class Ring>
bool testArithmeticConsistency (const Ring &F, const char *name, unsigned int iterations)
{
	bool pass = true;

	pass = testRingArithmeticConsistency (F, name, iterations) && pass;
	pass = testInvDivConsistency (F, name, iterations) && pass;

	return pass;
}

template <class Ring>
bool testRingArithmeticConsistency (const Ring &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "Testing " << name << " in-place/out-of-place arithmetic consistency" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingArithmeticConsistency", iterations);

	bool ret = true;

	typename Ring::RandIter r (F);
	typename Ring::Element a, b, c1, c2;
        F.init (a,0); F.init (b,0); F.init (c1,0); F.init (c2,0);

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);

		// This should be immaterial, since we have already "proven" commutativity
		if (F.isZero (a) && !F.isZero (b))
			std::swap (a, b);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;

		F.add (c1, a, b);
		F.copy (c2, a);
		F.addin (c2, b);

		report << "a + b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for addition", ret);

		F.sub (c1, a, b);
		F.copy (c2, a);
		F.subin (c2, b);

		report << "a - b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for subtraction", ret);

		F.neg (c1, a);
		F.copy (c2, a);
		F.negin (c2);

		report << "-a = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for negation", ret);

		F.mul (c1, a, b);
		F.copy (c2, a);
		F.mulin (c2, b);

		report << "a * b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for multiplication", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingArithmeticConsistency");
	delete[] st;
	return ret;
}

template <class Ring>
bool testInvDivConsistency (const Ring &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "Testing " << name << " in-place/out-of-place inv and div consistency" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testInvDivConsistency", iterations);

	bool ret = true;

	typename Ring::RandIter r (F);
	typename Ring::Element a, b, c1, c2;
	F.init (a,0); F.init (b,0); F.init (c1,0); F.init (c2,0);

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a); r.random (b);

		// This should be immaterial, since we have already "proven" commutativity
		if (F.isZero (a) && !F.isZero (b))
		std::swap (a, b);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;


		if (!F.isZero (a)) {
			F.div (c1, b, a);
			F.copy (c2, b);
			F.divin (c2, a);

			report << "b / a = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

			if (!F.areEqual (c1, c2)) reportError("Consistency failure for division", ret);

			F.inv (c1, a);
			F.copy (c2, a);
			F.invin (c2);

			report << "a^-1 = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

			if (!F.areEqual (c1, c2)) reportError("Consistency failure for inversion", ret);
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvDivConsistency");
	delete[] st;
	return ret;
}

/** Generic test 8: Consistency of axpy
 *
 * Generates random elements 'a', 'x', and 'y' and checks that a * x + y is the
 * same for axpy, axpyin, add/mul
 */

template <class Ring>
bool testAxpyConsistency (const Ring &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "Testing " << name << " axpy/add-mul consistency" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testAxpyConsistency", iterations);

	bool ret = true;

	typename Ring::RandIter r (F);
	typename Ring::Element a, x, y, c1, c2, c3;
        F.init (a,0); F.init (x,0); F.init (y,0); 
        F.init (c1,0); F.init (c2,0); F.init (c3,0);

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (x);
		r.random (y);

		ostream &report = commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", x = ";
		F.write (report, x) << ", y = ";
		F.write (report, y) << endl;

		F.mul (c1, a, x);
		F.addin (c1, y);
		F.axpy (c2, a, x, y);
		F.copy (c3, y);
		F.axpyin (c3, a, x);

		report << "a * x + y = (add-mul) ";
		F.write (report, c1) << ", (out-of-place) ";
		F.write (report, c2) << ", (in-place) ";
		F.write (report, c3) << endl;

		if (!F.areEqual (c1, c2) || !F.areEqual (c1, c3)) reportError("Consistency failure for axpy", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAxpyConsistency");
	delete[] st;
	return ret;
}

/** Generic test 9: Basic concept check of RandIter
 *
 * In a loop, generates random element 'a', and fails
 * if it is always zero.
 */
template <class Ring>
bool testRanditerBasic (const Ring &F, const char *name, unsigned int iterations)
{
	bool ret=false;
	std::ostringstream str;
	str << "Testing " << name << " randiter basic operation " << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRanditerBasic", iterations);

	 typename Ring::RandIter r (F);
	 typename Ring::Element a;
		 F.init (a,0);

	 if (iterations < 20) iterations = 20;
	 for (unsigned int i = 0; i < iterations; i++) {
		 r.random (a);
		 if ( ! F.isZero(a) ) {ret = true; break;}

	 }

	 commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRanditerBasic");
	 delete[] st;
	 return ret;
 }

template <class Ring>
bool runBasicRingTests (const Ring &F, const char *desc, unsigned int iterations, bool runCharacteristicTest = true) 
{
	bool pass = true;
	ostringstream str;

	str << "Testing " << desc << " ring" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());

	commentator.start (st, "runBasicRingTests", runCharacteristicTest ? 11 : 10);
	
	if (!testRing                 (F, string (str.str ()).c_str ())) pass = false; commentator.progress ();
	if (!testRingNegation         (F, desc, iterations))             pass = false; commentator.progress ();
	if (!testRingDistributivity   (F, desc, iterations))             pass = false; commentator.progress ();
	if (!testRingAssociativity    (F, desc, iterations))             pass = false; commentator.progress ();

	if (runCharacteristicTest) {
		if (!testRingCharacteristic (F, desc, iterations))
			pass = false;

		commentator.progress ();
	}

	if (!testGeometricSummation        (F, desc, iterations, 100))   pass = false; commentator.progress ();
	if (!testRingArithmeticConsistency (F, desc, iterations))        pass = false; commentator.progress ();
	if (!testAxpyConsistency           (F, desc, iterations))        pass = false; commentator.progress ();
	if (!testRanditerBasic             (F, desc, iterations))        pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runBasicRingTests");
	delete[] st;
	return pass;
}

 
 /* Convenience function to run all of the ring tests on a given ring */

template <class Ring>
bool runRingTests (const Ring &F, const char *desc, unsigned int iterations, bool runCharacteristicTest = true, bool runInversionTests = true) 
{
	ostringstream str;

	str << "Testing " << desc << " ring" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "runRingTests");
	bool ret = true;

	ret = runBasicRingTests (F, desc, iterations, runCharacteristicTest) && ret;

	if (runInversionTests) {
		ret = testInvDivConsistency (F, desc, iterations) && ret;
		ret = testRingInversion (F, desc, iterations) && ret;
	}

	ret = testRingCommutativity (F, desc, iterations) && ret;
	ret = testFreshmansDream (F, desc, iterations) && ret;

	commentator.stop (MSG_STATUS (ret));
	delete[] st;
	return ret;
}

#endif // __LELA_TESTS_TEST_RING_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
