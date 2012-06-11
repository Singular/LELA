/* tests/test-blas-level1.h
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic test suite for BLAS Level 1 routines
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_BLAS_LEVEL1_H
#define __LELA_TESTS_TEST_BLAS_LEVEL1_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>

#include "lela/util/commentator.h"
#include "lela/vector/stream.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"

#include "test-common.h" 

/** Test 1: Copy and areEqual
 *
 * Constructs a random vector and copies it to another vector. Then checks equality.
 *
 * F - Ring over which to perform computations
 * text - Text describing types of vectors
 * stream - Stream generating vectors
 * stream2 - Dummy stream of second vector type to trick the compiler
 *
 * Return true on success and false on failure
 */
template <class Ring, class Modules, class Vector1, class Vector2>
static bool testCopyEqual (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector1> &stream, LELA::VectorStream<Vector2> &stream2) 
{
	std::ostringstream str;

	str << "Testing " << text << " vector copy, equal" << std::ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__, stream.size ());

	bool ret = true;
	bool iter_passed;

	Vector1 v;
	Vector2 w;

	LELA::VectorUtils::ensureDim<Ring, Vector1> (v, stream.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector2> (w, stream.dim ());

	while (stream) {
		LELA::commentator.startIteration (stream.pos ());

		iter_passed = true;

		stream.get (v);
		LELA::BLAS1::copy (ctx, v, w);

		std::ostream &reportUI = LELA::commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		reportUI << "Input vector:   ";
		LELA::BLAS1::write (ctx, reportUI, v) << std::endl;

		reportUI << "Output vector:  ";
		LELA::BLAS1::write (ctx, reportUI, w) << std::endl;

		if (!LELA::BLAS1::equal (ctx, v, w))
			ret = iter_passed = false;

		if (!iter_passed)
			LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << std::endl;

		LELA::commentator.stop ("done");
		LELA::commentator.progress ();
	}

	LELA::commentator.stop (MSG_STATUS (ret));

	stream.reset ();
	stream2.reset ();

	return ret;
}

/** Test 1a: Inequality of vectors
 *
 * Take a vector and copy it. Change one copy in each position in
 * turn. Check that the two are not equal.
 *
 * ctx - Context-object
 * text - Text describing types of vectors
 * stream1 - Source for first vector
 * stream2 - Not used except as exemplar of type
 *
 * Return true on success and false on failure
 */

// Set the given entry of the input-vector (assumed zero) to the given value (assumed nonzero)
template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LELA::VectorRepresentationTypes::Dense)
	{ v[idx] = e; }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LELA::VectorRepresentationTypes::Sparse)
	{ v.push_back (typename Vector::value_type (idx, e)); }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LELA::VectorRepresentationTypes::Dense01)
	{ v[idx] = true; }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LELA::VectorRepresentationTypes::Sparse01)
	{ v.push_back (idx); }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LELA::VectorRepresentationTypes::Hybrid01)
	{ v.push_back (typename Vector::value_type (idx >> LELA::WordTraits<typename Vector::word_type>::logof_size,
						    Vector::Endianness::e_j (idx & LELA::WordTraits<typename Vector::word_type>::pos_mask))); }

template <class Ring, class Vector>
void set_entry (Vector &v, size_t idx, const typename Ring::Element &e)
	{ set_entry_spec (v, idx, e, typename LELA::VectorTraits<Ring, Vector>::RepresentationType ()); }

template <class Ring, class Modules, class Vector1, class Vector2>
bool testInequality (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector1> &stream1, LELA::VectorStream<Vector2> &stream2)
{
	bool pass = true;

	std::ostringstream str;

	str << "Testing " << text << " vector inequality" << std::ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	Vector1 v;
	typename LELA::VectorTraits<Ring, Vector2>::ContainerType v_c, e_i;

	LELA::VectorUtils::ensureDim<Ring, Vector1> (v, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector2> (v_c, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector2> (e_i, stream1.dim ());

	while (stream1) {
		stream1 >> v;

		std::ostream &reportUI = LELA::commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		reportUI << "Testing with vector: ";
		LELA::BLAS1::write (ctx, reportUI, v) << std::endl;

		for (size_t i = 0; i < stream1.dim (); ++i) {
			LELA::BLAS1::copy (ctx, v, v_c);
			LELA::BLAS1::scal (ctx, ctx.F.zero (), e_i);
			set_entry<Ring, typename LELA::VectorTraits<Ring, Vector2>::ContainerType> (e_i, i, ctx.F.one ());
			LELA::BLAS1::axpy (ctx, ctx.F.one (), e_i, v_c);

			if (LELA::BLAS1::equal (ctx, v, v_c)) {
				std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
				error << "ERROR: Vectors still reported as equal after change at position " << i << std::endl;
				error << "Original       : ";
				LELA::BLAS1::write (ctx, error, v) << std::endl;
				error << "Modified copy  : ";
				LELA::BLAS1::write (ctx, error, v_c) << std::endl;
				pass = false;
			}
		}
	}

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

/** Test 1b: Nonzero vectors
 *
 * Test whether the vector e_i is nonzero for each i from 0 to dim - 1
 *
 * ctx - Context-object
 * text - Text describing types of vectors
 * stream - Not used except as exemplar of type
 *
 * Return true on success and false on failure
 */

template <class Ring, class Modules, class Vector>
bool testNonzero (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector> &stream)
{
	bool pass = true;

	std::ostringstream str;

	str << "Testing " << text << " vector nonzero" << std::ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__, stream.size ());

	typename LELA::VectorTraits<Ring, Vector>::ContainerType e_i;

	LELA::VectorUtils::ensureDim<Ring, Vector> (e_i, stream.dim ());

	for (size_t i = 0; i < stream.dim (); ++i) {
		LELA::BLAS1::scal (ctx, ctx.F.zero (), e_i);
		set_entry<Ring, typename LELA::VectorTraits<Ring, Vector>::ContainerType> (e_i, i, ctx.F.one ());

		if (LELA::BLAS1::is_zero (ctx, e_i)) {
			std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "ERROR: Vector e_" << i << " reported zero" << std::endl;
			error << "Vector e_i: ";
			LELA::BLAS1::write (ctx, error, e_i) << std::endl;
			pass = false;
		}
	}

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

/** Test 2: Dot product of vectors
 *
 * Construct two random vectors and compute their dot product
 *
 * ctx - Context-object
 * text - Text describing types of vectors
 * stream1 - Stream for first family of vectors
 * stream2 - Stream for second family of vectors
 *
 * Return true on success and false on failure
 */

/* Construct the dot-product of two vectors manually. Depends on
 * proper functioning of VectorUtils::getEntry and of
 * Ring::axpyin.
 */

template <class Ring, class Vector1, class Vector2>
typename Ring::Element &manualDot (const Ring &F, typename Ring::Element &x, const Vector1 &v1, const Vector2 &v2, size_t dim)
{
	typename Ring::Element a, b;

	F.copy (x, F.zero ());

	for (size_t j = 0; j < dim; j++)
		if (LELA::VectorUtils::getEntry (v1, a, j) && LELA::VectorUtils::getEntry (v2, b, j))
			F.axpyin (x, a, b);

	return x;
}

template <class Ring, class Modules, class Vector1, class Vector2>
static bool testDotProduct (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector1> &stream1, LELA::VectorStream<Vector2> &stream2) 
{
	lela_check (stream1.dim () == stream2.dim ());

	std::ostringstream str;

	str << "Testing " << text << " dot product" << std::ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;

	Vector1 v1;
	Vector2 v2;
	typename Ring::Element sigma, rho;

	LELA::VectorUtils::ensureDim<Ring, Vector1> (v1, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector2> (v2, stream2.dim ());

	LELA::Timer timer;
	double totaltime = 0.0;

	while (stream1 && stream2) {
		LELA::commentator.startIteration (stream1.pos ());

		stream1.get (v1);
		stream2.get (v2);

		manualDot (ctx.F, sigma, v1, v2, stream1.dim ());

		std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		LELA::BLAS1::write (ctx, report, v1) << std::endl;

		report << "Input vector 2 of size " << v2.size() << ":  ";
		LELA::BLAS1::write (ctx, report, v2) << std::endl;

		timer.start ();
		LELA::BLAS1::dot (ctx, rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		report << "True dot product: ";
		ctx.F.write (report, sigma) << std::endl;

		report << "Dot product from vector domain: ";
		ctx.F.write (report, rho) << std::endl;

		if (!ctx.F.areEqual (sigma, rho)) {
			ret = false;
			std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "ERROR: Dot products are not equal" << std::endl;
		}

		LELA::commentator.stop ("done");
		LELA::commentator.progress ();
	}

	LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
		<< "Average time for dot product: " << totaltime / stream1.size () << std::endl;

	LELA::commentator.stop (MSG_STATUS (ret));

	stream1.reset ();
	stream2.reset ();

	return ret;
}

/** Test 2.1: swap
 *
 * Construct two random vector v1 and v2
 * Swap vector v1 and vector v2 
 * Check whether the results are equal.
 *
 * F - Ring over which to perform computations
 * text - Text describing types of vectors
 * stream1 - Source of random vectors v1
 *
 * Return true on success and false on failure
 */

template <class Ring, class Modules, class Vector1>
bool testswap (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector1> &stream1)
{
	ostringstream str;

	str << "Testing " << text << " test swap " << ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__);
	std::ostream &reportUI = LELA::commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;
	
	typename LELA::VectorTraits<Ring,Vector1>::ContainerType v1, v2, w1, w2;

	LELA::VectorUtils::ensureDim<Ring, Vector1> (v1, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector1> (v2, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector1> (w1, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector1> (w2, stream1.dim ());

	stream1 >> v1;

	stream1.reset ();

        stream1 >> v2;

        while(LELA::BLAS1::equal(ctx, v1 , v2)) {
		stream1.reset ();
		stream1 >> v2;
	}

	LELA::BLAS1::copy(ctx, v1, w1);
	LELA::BLAS1::copy(ctx, v2, w2);

	reportUI << "(Before swap) Vector v_1: ";
	LELA::BLAS1::write (ctx, reportUI, v1) << std::endl;

	reportUI << "(Before swap) Vector v_2: ";
	LELA::BLAS1::write (ctx, reportUI, v2) << std::endl;

	LELA::BLAS1::swap(ctx, w1, w2);

	reportUI << "(After swap) Vector v_1: ";
	LELA::BLAS1::write(ctx, reportUI, w1) << std::endl;

	reportUI << "(After swap) Vector v_2: ";
	LELA::BLAS1::write(ctx, reportUI, w2) << std::endl;

	if((LELA::BLAS1::equal(ctx, v1 , w1)) || (LELA::BLAS1::equal(ctx ,v2, w2)))
	{
		ret = false;
		std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		error << "ERROR: Swap unseccessfull! " << std::endl;
	}

	LELA::commentator.stop (MSG_STATUS (ret));

	return ret;
}

/** Test 3: Vector scal
 *
 * Construct a random vector x and a random element a and compute -a*x
 * + a*x using axpy and scal. Check whether the results are equal.
 *
 * F - Ring over which to perform computations
 * text - Text describing types of vectors
 * stream1 - Source of random vectors x
 *
 * Return true on success and false on failure
 */

template <class Ring, class Modules, class Vector>
static bool testScal (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector> &stream1) 
{
	std::ostringstream str;

	str << "Testing " << text << " vector scal" << std::ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2;
	typename Ring::Element a;
	typename Ring::Element nega;
	typename Ring::RandIter r (ctx.F);

	LELA::VectorUtils::ensureDim<Ring, Vector> (v1, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector> (v2, stream1.dim ());

	while (stream1) {
		LELA::commentator.startIteration (stream1.pos ());

		iter_passed = true;

		stream1.get (v1);

		do r.random (a); while (ctx.F.isZero (a));

		std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		std::ostream &reportUI = LELA::commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		report << "Input vector 1 of size " << v1.size() << ":  ";
		LELA::BLAS1::write (ctx, reportUI, v1) << std::endl;

		report << "Element a:  ";
		ctx.F.write (report, a) << std::endl;

		ctx.F.neg (nega, a);

		LELA::BLAS1::copy (ctx, v1, v2);

		LELA::BLAS1::scal (ctx, nega, v2);

		reportUI << "         -a * x = ";
		LELA::BLAS1::write (ctx, reportUI, v2) << std::endl;
		reportUI.flush ();

		LELA::BLAS1::axpy (ctx, a, v1, v2);

		reportUI << " a * x + -a * x = ";
		LELA::BLAS1::write (ctx, reportUI, v2) << std::endl;

		if (!LELA::BLAS1::is_zero (ctx, v2))
			ret = iter_passed = false;

		if (!iter_passed)
			LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a * x + -a * x != 0" << std::endl;

		LELA::commentator.stop ("done");
		LELA::commentator.progress ();
	}

	LELA::commentator.stop (MSG_STATUS (ret));

	stream1.reset ();

	return ret;
}

/** Test 4: Vector axpy
 *
 * Construct two random vectors x and y and a random element a and compute (x +
 * a*y) - a*(y + a^-1*x) using vector axpy. Check whether the result is 0.
 *
 * F - Ring over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Ring, class Modules, class Vector>
static bool testAXPY (LELA::Context<Ring, Modules> &ctx, const char *text, LELA::VectorStream<Vector> &stream1, LELA::VectorStream<Vector> &stream2) 
{
	std::ostringstream str;
	str << "Testing " << text << " vector axpy" << std::ends;
	LELA::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3;
	typename Ring::Element a, ainv, aneg;
	typename Ring::RandIter r (ctx.F);

	LELA::VectorUtils::ensureDim<Ring, Vector> (v1, stream1.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector> (v2, stream2.dim ());
	LELA::VectorUtils::ensureDim<Ring, Vector> (v3, stream1.dim ());

	while (stream1 && stream2) {
		LELA::commentator.startIteration (stream1.pos ());

		iter_passed = true;

		stream1.get (v1);
		stream2.get (v2);

		do r.random (a); while (ctx.F.isZero (a));

		std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		std::ostream &reportUI = LELA::commentator.report (LELA::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		reportUI << "Input vector 1 v_1 of size " << v1.size() << ":  ";
		LELA::BLAS1::write (ctx, reportUI, v1) << std::endl;

		reportUI << "Input vector 2 v_2 of size " << v2.size() << ":  ";
		LELA::BLAS1::write (ctx, reportUI, v2) << std::endl;

		report << "Element a:  ";
		ctx.F.write (report, a) << std::endl;

		ctx.F.inv (ainv, a);
		ctx.F.neg (aneg, a);

		LELA::BLAS1::copy (ctx, v1, v3);
		LELA::BLAS1::axpy (ctx, a, v2, v3);

		reportUI << "av_2 + v_1 = ";
		LELA::BLAS1::write (ctx, reportUI, v3) << std::endl;

		LELA::BLAS1::axpy (ctx, ainv, v1, v2);

		reportUI << "a^-1 v_1 + v_2 = ";
		LELA::BLAS1::write (ctx, reportUI, v3) << std::endl;

		LELA::BLAS1::axpy (ctx, aneg, v2, v3);

		reportUI << "(av_2 + v_1) + -a (a^-1 v_1 + v_2) = ";
		LELA::BLAS1::write (ctx, reportUI, v3) << std::endl;

		if (!LELA::BLAS1::is_zero (ctx, v3))
			ret = iter_passed = false;

		if (!iter_passed)
			LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (av_2 + v_1) + -a (a^-1 v_1 + v_2) != 0" << std::endl;

		LELA::commentator.stop ("done");
		LELA::commentator.progress ();
	}

	LELA::commentator.stop (MSG_STATUS (ret));

	stream1.reset ();
	stream2.reset ();

	return ret;
}

template <class Ring, class Modules>
bool testBLAS1 (LELA::Context<Ring, Modules> &ctx, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing BLAS1 over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;

	LELA::RandomDenseStream<Ring, typename LELA::Vector<Ring>::Dense> stream1 (ctx.F, n, iterations), stream2 (ctx.F, n, iterations);
	LELA::RandomSparseStream<Ring, typename LELA::Vector<Ring>::Sparse> stream3 (ctx.F, 0.1, n, iterations), stream4 (ctx.F, 0.1, n, iterations);

	if (!testCopyEqual (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "dense/sparse", stream1, stream4)) pass = false;  stream1.reset (); stream4.reset ();
	if (!testCopyEqual (ctx, "sparse/dense", stream3, stream2)) pass = false;  stream3.reset (); stream2.reset ();
	if (!testCopyEqual (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();

	if (!testInequality (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testInequality (ctx, "dense/sparse", stream1, stream4)) pass = false;  stream1.reset (); stream4.reset ();
	if (!testInequality (ctx, "sparse/dense", stream3, stream2)) pass = false;  stream3.reset (); stream2.reset ();
	if (!testInequality (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();

	if (!testNonzero (ctx, "dense", stream1)) pass = false;  stream1.reset ();
	if (!testNonzero (ctx, "sparse", stream3)) pass = false; stream3.reset ();

	if (!testDotProduct (ctx, "dense/dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testDotProduct (ctx, "sparse/dense", stream3, stream1)) pass = false;  stream3.reset (); stream1.reset ();
	if (!testDotProduct (ctx, "sparse/sparse", stream3, stream4)) pass = false; stream3.reset (); stream4.reset ();

	if (!testScal (ctx, "dense", stream1)) pass = false;  stream1.reset ();
	if (!testScal (ctx, "sparse", stream3)) pass = false; stream3.reset ();

	if (!testAXPY (ctx, "dense", stream1, stream2)) pass = false;   stream1.reset (); stream2.reset ();
	if (!testAXPY (ctx, "sparse", stream3, stream4)) pass = false;  stream3.reset (); stream4.reset ();

	if (!testswap(ctx, "dense", stream1)) pass = false; stream1.reset ();
        if (!testswap(ctx, "sparse", stream3)) pass = false; stream3.reset ();


	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

// Tests to check that the result is the same when using two modules or two representation-types

template <class Ring, class Modules1, class Modules2, class Vector1, class Vector2, class Vector3, class Vector4>
bool testDotConsistency (LELA::Context<Ring, Modules1> &ctx1,
			 LELA::Context<Ring, Modules2> &ctx2,
			 const char *text,
			 LELA::VectorStream<Vector1> &stream1, 
			 LELA::VectorStream<Vector2> &stream2,
			 LELA::VectorStream<Vector3> &stream3,
			 LELA::VectorStream<Vector4> &stream4)
{
	using namespace LELA;

	bool pass = true;

	std::ostringstream str;
	str << "Testing " << text << " vector dot consistency" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Vector1 v1;
	Vector2 v2;
	Vector3 w1;
	Vector4 w2;

	VectorUtils::ensureDim<Ring, Vector1> (v1, stream1.dim ());
	VectorUtils::ensureDim<Ring, Vector2> (v2, stream2.dim ());
	VectorUtils::ensureDim<Ring, Vector3> (w1, stream1.dim ());
	VectorUtils::ensureDim<Ring, Vector4> (w2, stream2.dim ());

	typename Ring::Element d1, d2;

	while (stream1 && stream2) {
		stream1 >> v1;
		stream2 >> v2;

		BLAS1::copy (ctx1, v1, w1);
		BLAS1::copy (ctx1, v2, w2);

		reportUI << "Vector v_1: ";
		BLAS1::write (ctx1, reportUI, v1) << std::endl;

		reportUI << "Vector v_2: ";
		BLAS1::write (ctx1, reportUI, v2) << std::endl;

		BLAS1::dot (ctx1, d1, v1, v2);

		reportUI << "Dot product <v_1,v_2>: ";
		ctx1.F.write (reportUI, d1) << std::endl;

		reportUI << "Vector w_1: ";
		BLAS1::write (ctx2, reportUI, w1) << std::endl;

		reportUI << "Vector w_2: ";
		BLAS1::write (ctx2, reportUI, w2) << std::endl;

		BLAS1::dot (ctx2, d2, w1, w2);

		reportUI << "Dot product <w_1,w_2>: ";
		ctx2.F.write (reportUI, d2) << std::endl;

		if (!ctx1.F.areEqual (d1, d2)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Dot products are not equal" << std::endl;
			pass = false;
		}
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Vector1, class Vector2, class Vector3, class Vector4>
bool testAxpyConsistency (LELA::Context<Ring, Modules1> &ctx1,
			  LELA::Context<Ring, Modules2> &ctx2,
			  const char *text,
			  LELA::VectorStream<Vector1> &stream1, 
			  LELA::VectorStream<Vector2> &stream2,
			  LELA::VectorStream<Vector3> &stream3,
			  LELA::VectorStream<Vector4> &stream4)
{
	using namespace LELA;

	bool pass = true;

	std::ostringstream str;
	str << "Testing " << text << " vector axpy consistency" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Vector1 v1;
	Vector2 v2;
	Vector3 w1;
	Vector4 w2;

	VectorUtils::ensureDim<Ring, Vector1> (v1, stream1.dim ());
	VectorUtils::ensureDim<Ring, Vector2> (v2, stream2.dim ());
	VectorUtils::ensureDim<Ring, Vector3> (w1, stream1.dim ());
	VectorUtils::ensureDim<Ring, Vector4> (w2, stream2.dim ());

	typename Ring::Element a;
	NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

	while (stream1 && stream2) {
		stream1 >> v1;
		stream2 >> v2;

		BLAS1::copy (ctx1, v1, w1);
		BLAS1::copy (ctx1, v2, w2);

		r.random (a);

		report << "Coefficient a: ";
		ctx1.F.write (report, a) << std::endl;

		reportUI << "Vector v_1: ";
		BLAS1::write (ctx1, reportUI, v1) << std::endl;

		reportUI << "Vector v_2: ";
		BLAS1::write (ctx1, reportUI, v2) << std::endl;

		BLAS1::axpy (ctx1, a, v1, v2);

		reportUI << "Vector a v_1 + v_2: ";
		BLAS1::write (ctx1, reportUI, v2) << std::endl;

		reportUI << "Vector w_1: ";
		BLAS1::write (ctx2, reportUI, w1) << std::endl;

		reportUI << "Vector w_2: ";
		BLAS1::write (ctx2, reportUI, w2) << std::endl;

		BLAS1::axpy (ctx2, a, w1, w2);

		reportUI << "Vector a w_1 + w_2: ";
		BLAS1::write (ctx1, reportUI, w2) << std::endl;

		if (!BLAS1::equal (ctx1, v2, w2)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a v_1 + v_2 != a w_1 + w_2" << std::endl;
			pass = false;
		}
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2, class Vector1, class Vector2>
bool testScalConsistency (LELA::Context<Ring, Modules1> &ctx1,
			  LELA::Context<Ring, Modules2> &ctx2,
			  const char *text,
			  LELA::VectorStream<Vector1> &stream1, 
			  LELA::VectorStream<Vector2> &stream2)
{
	using namespace LELA;

	bool pass = true;

	std::ostringstream str;
	str << "Testing " << text << " vector scal consistency" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Vector1 v;
	Vector2 w;

	VectorUtils::ensureDim<Ring, Vector1> (v, stream1.dim ());
	VectorUtils::ensureDim<Ring, Vector2> (w, stream1.dim ());

	typename Ring::Element a;
	NonzeroRandIter<Ring, typename Ring::RandIter> r (ctx1.F, typename Ring::RandIter (ctx1.F));

	while (stream1) {
		stream1 >> v;

		BLAS1::copy (ctx1, v, w);

		r.random (a);

		report << "Coefficient a: ";
		ctx1.F.write (report, a) << std::endl;

		reportUI << "Vector v: ";
		BLAS1::write (ctx1, reportUI, v) << std::endl;

		BLAS1::scal (ctx1, a, v);

		reportUI << "Vector a v: ";
		BLAS1::write (ctx1, reportUI, v) << std::endl;

		reportUI << "Vector w: ";
		BLAS1::write (ctx2, reportUI, w) << std::endl;

		BLAS1::scal (ctx2, a, w);

		reportUI << "Vector a w: ";
		BLAS1::write (ctx1, reportUI, w) << std::endl;

		if (!BLAS1::equal (ctx1, v, w)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a v != a w" << std::endl;
			pass = false;
		}
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules1, class Modules2>
bool testBLAS1ModulesConsistency (LELA::Context<Ring, Modules1> &ctx1, LELA::Context<Ring, Modules2> &ctx2, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing BLAS1 consistency over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;

	LELA::RandomDenseStream<Ring, typename LELA::Vector<Ring>::Dense> stream1 (ctx1.F, n, iterations), stream2 (ctx1.F, n, iterations);
	LELA::RandomSparseStream<Ring, typename LELA::Vector<Ring>::Sparse> stream3 (ctx1.F, 0.1, n, iterations), stream4 (ctx1.F, 0.1, n, iterations);

	pass = testDotConsistency (ctx1, ctx2, "dense/dense", stream1, stream2, stream1, stream1) && pass; stream1.reset (); stream2.reset ();
	pass = testDotConsistency (ctx1, ctx2, "dense/sparse", stream1, stream3, stream1, stream3) && pass; stream1.reset (); stream3.reset ();
	pass = testDotConsistency (ctx1, ctx2, "sparse/sparse", stream3, stream4, stream3, stream3) && pass; stream3.reset (); stream4.reset ();

	pass = testAxpyConsistency (ctx1, ctx2, "dense/dense", stream1, stream2, stream1, stream1) && pass; stream1.reset (); stream2.reset ();
	pass = testAxpyConsistency (ctx1, ctx2, "dense/sparse", stream1, stream3, stream1, stream3) && pass; stream1.reset (); stream3.reset ();
	pass = testAxpyConsistency (ctx1, ctx2, "sparse/sparse", stream3, stream4, stream3, stream3) && pass; stream3.reset (); stream4.reset ();

	pass = testScalConsistency (ctx1, ctx2, "dense", stream1, stream2) && pass; stream1.reset ();
	pass = testScalConsistency (ctx1, ctx2, "sparse", stream3, stream4) && pass; stream3.reset ();

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Ring, class Modules>
bool testBLAS1RepsConsistency (LELA::Context<Ring, Modules> &ctx, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing BLAS1 consistency over <" << text << ">" << std::ends;
	LELA::commentator.start (str.str ().c_str ());

	bool pass = true;

	LELA::RandomDenseStream<Ring, typename LELA::Vector<Ring>::Dense> stream1 (ctx.F, n, iterations), stream2 (ctx.F, n, iterations);
	LELA::RandomSparseStream<Ring, typename LELA::Vector<Ring>::Sparse> stream3 (ctx.F, 0.1, n, iterations), stream4 (ctx.F, 0.1, n, iterations);

	pass = testDotConsistency (ctx, ctx, "dense/sparse with dense/dense", stream1, stream3, stream1, stream1) && pass; stream1.reset (); stream3.reset ();
	pass = testDotConsistency (ctx, ctx, "sparse/sparse with dense/dense", stream3, stream4, stream1, stream1) && pass; stream3.reset (); stream4.reset ();

	pass = testAxpyConsistency (ctx, ctx, "dense/sparse with dense/dense", stream1, stream3, stream1, stream1) && pass; stream1.reset (); stream3.reset ();
	pass = testAxpyConsistency (ctx, ctx, "sparse/sparse with dense/dense", stream3, stream4, stream1, stream1) && pass; stream3.reset (); stream4.reset ();

	pass = testScalConsistency (ctx, ctx, "sparse with dense ", stream3, stream1) && pass; stream1.reset (); stream3.reset ();

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

#endif // __LELA_TESTS_TEST_BLAS_LEVEL1_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
