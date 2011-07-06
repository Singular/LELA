/* tests/test-blas-level1.h
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Generic test suite for BLAS Level 1 routines
 */

#ifndef __LINBOX_TEST_BLAS_LEVEL1_H
#define __LINBOX_TEST_BLAS_LEVEL1_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "linbox/blas/context.h"
#include "linbox/blas/level1.h"

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
static bool testCopyEqual (LinBox::Context<Ring, Modules> &ctx, const char *text, LinBox::VectorStream<Vector1> &stream, LinBox::VectorStream<Vector2> &stream2) 
{
	std::ostringstream str;

	str << "Testing " << text << " vector copy, equal" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream.size ());

	bool ret = true;
	bool iter_passed;

	Vector1 v;
	Vector2 w;

	LinBox::VectorUtils::ensureDim<Ring, Vector1> (v, stream.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector2> (w, stream.dim ());

	while (stream) {
		LinBox::commentator.startIteration (stream.pos ());

		iter_passed = true;

		stream.get (v);
		LinBox::BLAS1::copy (ctx, v, w);

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:   ";
		LinBox::BLAS1::write (ctx, report, v) << std::endl;

		report << "Output vector:  ";
		LinBox::BLAS1::write (ctx, report, w) << std::endl;

		if (!LinBox::BLAS1::equal (ctx, v, w))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret));

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
void set_entry_spec (Vector &v, size_t idx, Element &e, LinBox::VectorRepresentationTypes::Dense)
	{ v[idx] = e; }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LinBox::VectorRepresentationTypes::Sparse)
	{ v.push_back (typename Vector::value_type (idx, e)); }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LinBox::VectorRepresentationTypes::Dense01)
	{ v[idx] = true; }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LinBox::VectorRepresentationTypes::Sparse01)
	{ v.push_back (idx); }

template <class Element, class Vector>
void set_entry_spec (Vector &v, size_t idx, Element &e, LinBox::VectorRepresentationTypes::Hybrid01)
	{ v.push_back (typename Vector::value_type (idx >> LinBox::WordTraits<typename Vector::word_type>::logof_size,
						    Vector::Endianness::e_j (idx & LinBox::WordTraits<typename Vector::word_type>::pos_mask))); }

template <class Ring, class Vector>
void set_entry (Vector &v, size_t idx, const typename Ring::Element &e)
	{ set_entry_spec (v, idx, e, typename LinBox::VectorTraits<Ring, Vector>::RepresentationType ()); }

template <class Ring, class Modules, class Vector1, class Vector2>
bool testInequality (LinBox::Context<Ring, Modules> &ctx, const char *text, LinBox::VectorStream<Vector1> &stream1, LinBox::VectorStream<Vector2> &stream2)
{
	bool pass = true;

	std::ostringstream str;

	str << "Testing " << text << " vector inequality" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	Vector1 v;
	typename LinBox::VectorTraits<Ring, Vector2>::ContainerType v_c, e_i;

	LinBox::VectorUtils::ensureDim<Ring, Vector1> (v, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector2> (v_c, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector2> (e_i, stream1.dim ());

	while (stream1) {
		stream1 >> v;

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Testing with vector: ";
		LinBox::BLAS1::write (ctx, report, v) << std::endl;

		for (size_t i; i < stream1.dim (); ++i) {
			LinBox::BLAS1::copy (ctx, v, v_c);
			LinBox::BLAS1::scal (ctx, ctx.F.zero (), e_i);
			set_entry<Ring, typename LinBox::VectorTraits<Ring, Vector2>::ContainerType> (e_i, i, ctx.F.one ());
			LinBox::BLAS1::axpy (ctx, ctx.F.one (), e_i, v_c);

			if (LinBox::BLAS1::equal (ctx, v, v_c)) {
				std::ostream &error = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
				error << "ERROR: Vectors still reported as equal after change at position " << i << std::endl;
				error << "Modified copy  : ";
				LinBox::BLAS1::write (ctx, error, v_c) << std::endl;
				pass = false;
			}
		}
	}

	LinBox::commentator.stop (MSG_STATUS (pass));

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
bool testNonzero (LinBox::Context<Ring, Modules> &ctx, const char *text, LinBox::VectorStream<Vector> &stream)
{
	bool pass = true;

	std::ostringstream str;

	str << "Testing " << text << " vector nonzero" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream.size ());

	typename LinBox::VectorTraits<Ring, Vector>::ContainerType e_i;

	LinBox::VectorUtils::ensureDim<Ring, Vector> (e_i, stream.dim ());

	for (size_t i; i < stream.dim (); ++i) {
		LinBox::BLAS1::scal (ctx, ctx.F.zero (), e_i);
		set_entry<Ring, typename LinBox::VectorTraits<Ring, Vector>::ContainerType> (e_i, i, ctx.F.one ());

		if (LinBox::BLAS1::is_zero (ctx, e_i)) {
			std::ostream &error = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "ERROR: Vector e_" << i << " reported zero" << std::endl;
			error << "Vector e_i: ";
			LinBox::BLAS1::write (ctx, error, e_i) << std::endl;
			pass = false;
		}
	}

	LinBox::commentator.stop (MSG_STATUS (pass));

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

	F.assign (x, F.zero ());

	for (size_t j = 0; j < dim; j++)
		if (LinBox::VectorUtils::getEntry (v1, a, j) && LinBox::VectorUtils::getEntry (v2, b, j))
			F.axpyin (x, a, b);

	return x;
}

template <class Ring, class Modules, class Vector1, class Vector2>
static bool testDotProduct (LinBox::Context<Ring, Modules> &ctx, const char *text, LinBox::VectorStream<Vector1> &stream1, LinBox::VectorStream<Vector2> &stream2) 
{
	linbox_check (stream1.dim () == stream2.dim ());

	std::ostringstream str;

	str << "Testing " << text << " dot product)" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;

	Vector1 v1;
	Vector2 v2;
	typename Ring::Element sigma, rho;

	LinBox::VectorUtils::ensureDim<Ring, Vector1> (v1, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector2> (v2, stream2.dim ());

	LinBox::Timer timer;
	double totaltime = 0.0;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.pos ());

		stream1.get (v1);
		stream2.get (v2);

		manualDot (ctx.F, sigma, v1, v2, stream1.dim ());

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		LinBox::BLAS1::write (ctx, report, v1) << std::endl;

		report << "Input vector 2 of size " << v2.size() << ":  ";
		LinBox::BLAS1::write (ctx, report, v2) << std::endl;

		timer.start ();
		LinBox::BLAS1::dot (ctx, rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		report << "True dot product: ";
		ctx.F.write (report, sigma) << std::endl;

		report << "Dot product from vector domain: ";
		ctx.F.write (report, rho) << std::endl;

		if (!ctx.F.areEqual (sigma, rho)) {
			ret = false;
			std::ostream &error = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			error << "ERROR: Dot products are not equal" << std::endl;
		}

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
		<< "Average time for dot product: " << totaltime / stream1.size () << std::endl;

	LinBox::commentator.stop (MSG_STATUS (ret));

	stream1.reset ();
	stream2.reset ();

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
static bool testScal (LinBox::Context<Ring, Modules> &ctx, const char *text, LinBox::VectorStream<Vector> &stream1) 
{
	std::ostringstream str;

	str << "Testing " << text << " vector scal" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2;
	typename Ring::Element a;
	typename Ring::Element nega;
	typename Ring::RandIter r (ctx.F);

	LinBox::VectorUtils::ensureDim<Ring, Vector> (v1, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector> (v2, stream1.dim ());

	while (stream1) {
		LinBox::commentator.startIteration (stream1.pos ());

		iter_passed = true;

		stream1.get (v1);

		do r.random (a); while (ctx.F.isZero (a));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		LinBox::BLAS1::write (ctx, report, v1) << std::endl;

		report << "Element a:  ";
		ctx.F.write (report, a) << std::endl;

		ctx.F.neg (nega, a);

		LinBox::BLAS1::copy (ctx, v1, v2);

		LinBox::BLAS1::scal (ctx, nega, v2);

		report << "         -a * x = ";
		LinBox::BLAS1::write (ctx, report, v2) << std::endl;
		report.flush ();

		LinBox::BLAS1::axpy (ctx, a, v1, v2);

		report << " a * x + -a * x = ";
		LinBox::BLAS1::write (ctx, report, v2) << std::endl;

		if (!LinBox::BLAS1::is_zero (ctx, v2))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a * x + -a * x != 0" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret));

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
static bool testAXPY (LinBox::Context<Ring, Modules> &ctx, const char *text, LinBox::VectorStream<Vector> &stream1, LinBox::VectorStream<Vector> &stream2) 
{
	std::ostringstream str;
	str << "Testing " << text << " vector axpy" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3;
	typename Ring::Element a, ainv, aneg;
	typename Ring::RandIter r (ctx.F);

	LinBox::VectorUtils::ensureDim<Ring, Vector> (v1, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector> (v2, stream2.dim ());
	LinBox::VectorUtils::ensureDim<Ring, Vector> (v3, stream1.dim ());

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.pos ());

		iter_passed = true;

		stream1.get (v1);
		stream2.get (v2);

		do r.random (a); while (ctx.F.isZero (a));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 v_1 of size " << v1.size() << ":  ";
		LinBox::BLAS1::write (ctx, report, v1) << std::endl;

		report << "Input vector 2 v_2 of size " << v2.size() << ":  ";
		LinBox::BLAS1::write (ctx, report, v2) << std::endl;

		report << "Element a:  ";
		ctx.F.write (report, a) << std::endl;

		ctx.F.inv (ainv, a);
		ctx.F.neg (aneg, a);

		LinBox::BLAS1::copy (ctx, v1, v3);
		LinBox::BLAS1::axpy (ctx, a, v2, v3);

		report << "av_2 + v_1 = ";
		LinBox::BLAS1::write (ctx, report, v3) << std::endl;

		LinBox::BLAS1::axpy (ctx, ainv, v1, v2);

		report << "a^-1 v_1 + v_2 = ";
		LinBox::BLAS1::write (ctx, report, v3) << std::endl;

		LinBox::BLAS1::axpy (ctx, aneg, v2, v3);

		report << "(av_2 + v_1) + -a (a^-1 v_1 + v_2) = ";
		LinBox::BLAS1::write (ctx, report, v3) << std::endl;

		if (!LinBox::BLAS1::is_zero (ctx, v3))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (av_2 + v_1) + -a (a^-1 v_1 + v_2) != 0" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret));

	stream1.reset ();
	stream2.reset ();

	return ret;
}

template <class Ring, class Modules>
bool testBLAS1 (LinBox::Context<Ring, Modules> &ctx, const char *text, size_t n, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing VectorDomain <" << text << ">" << std::ends;
	LinBox::commentator.start (str.str ().c_str ());

	bool pass = true;

	LinBox::RandomDenseStream<Ring, typename LinBox::Vector<Ring>::Dense> stream1 (ctx.F, n, iterations), stream2 (ctx.F, n, iterations);
	LinBox::RandomSparseStream<Ring, typename LinBox::Vector<Ring>::Sparse> stream3 (ctx.F, 0.1, n, iterations), stream4 (ctx.F, 0.1, n, iterations);

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

	LinBox::commentator.stop (MSG_STATUS (pass));

	return pass;
}

#endif // __LINBOX_TEST_BLAS_LEVEL1_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
