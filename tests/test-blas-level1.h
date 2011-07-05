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
 * F - Field over which to perform computations
 * text - Text describing types of vectors
 * stream - Stream generating vectors
 * stream2 - Dummy stream of second vector type to trick the compiler
 *
 * Return true on success and false on failure
 */
template <class Field, class Modules, class Vector1, class Vector2>
static bool testCopyEqual (LinBox::Context<Field, Modules> &ctx, const char *text, LinBox::VectorStream<Vector1> &stream, LinBox::VectorStream<Vector2> &stream2) 
{
	std::ostringstream str;

	str << "Testing " << text << " vector copy, equal" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream.size ());

	bool ret = true;
	bool iter_passed;

	Vector1 v;
	Vector2 w;

	LinBox::VectorUtils::ensureDim<Field, Vector1> (v, stream.dim ());
	LinBox::VectorUtils::ensureDim<Field, Vector2> (w, stream.dim ());

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
 * Field::axpyin.
 */

template <class Field, class Vector1, class Vector2>
typename Field::Element &manualDot (const Field &F, typename Field::Element &x, const Vector1 &v1, const Vector2 &v2, size_t dim)
{
	typename Field::Element a, b;

	F.assign (x, F.zero ());

	for (size_t j = 0; j < dim; j++)
		if (LinBox::VectorUtils::getEntry (v1, a, j) && LinBox::VectorUtils::getEntry (v2, b, j))
			F.axpyin (x, a, b);

	return x;
}

template <class Field, class Modules, class Vector1, class Vector2>
static bool testDotProduct (LinBox::Context<Field, Modules> &ctx, const char *text, LinBox::VectorStream<Vector1> &stream1, LinBox::VectorStream<Vector2> &stream2) 
{
	linbox_check (stream1.dim () == stream2.dim ());

	std::ostringstream str;

	str << "Testing " << text << " dot product)" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;

	Vector1 v1;
	Vector2 v2;
	typename Field::Element sigma, rho;

	LinBox::VectorUtils::ensureDim<Field, Vector1> (v1, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Field, Vector2> (v2, stream2.dim ());

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
 * F - Field over which to perform computations
 * text - Text describing types of vectors
 * stream1 - Source of random vectors x
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Vector>
static bool testScal (LinBox::Context<Field, Modules> &ctx, const char *text, LinBox::VectorStream<Vector> &stream1) 
{
	std::ostringstream str;

	str << "Testing " << text << " vector scal" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2;
	typename Field::Element a;
	typename Field::Element nega;
	typename Field::RandIter r (ctx.F);

	LinBox::VectorUtils::ensureDim<Field, Vector> (v1, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Field, Vector> (v2, stream1.dim ());

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
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field, class Modules, class Vector>
static bool testAXPY (LinBox::Context<Field, Modules> &ctx, const char *text, LinBox::VectorStream<Vector> &stream1, LinBox::VectorStream<Vector> &stream2) 
{
	std::ostringstream str;
	str << "Testing " << text << " vector axpy" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3;
	typename Field::Element a, ainv, aneg;
	typename Field::RandIter r (ctx.F);

	LinBox::VectorUtils::ensureDim<Field, Vector> (v1, stream1.dim ());
	LinBox::VectorUtils::ensureDim<Field, Vector> (v2, stream2.dim ());
	LinBox::VectorUtils::ensureDim<Field, Vector> (v3, stream1.dim ());

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

#endif // __LINBOX_TEST_BLAS_LEVEL1_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
