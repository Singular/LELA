/* lela/tests/test-vector.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic tests for vectors
 *
 * See COPYING for license information.
 */

#ifndef __LELA_TESTS_TEST_VECTOR_H
#define __LELA_TESTS_TEST_VECTOR_H

#include <iostream>

#include "lela/util/commentator.h"
#include "lela/vector/traits.h"
#include "lela/blas/level1.h"

// Some parameters for the tests
#define DEFAULT_DIM 50;
#define DEFAULT_CYCLE_LEN 10;
#define DEFAULT_GAP_CYCLE_LEN 5;
#define DEFAULT_MID_SUBVECTOR_START_IDX 10;

template <class Vector, class Ring>
bool testDenseConstSubvector (const Ring &R, const Vector &v, typename Vector::ConstSubvectorType &v_sub, size_t offset)
{
	typename Vector::ConstSubvectorType::const_iterator i_sub;

	size_t idx;

	bool pass = true;

	for (j = v.begin (); j != v.end (); ++j)
		R.init (*j, (idx + offset) % DEFAULT_CYCLE_LEN);

	for (i_sub = v_sub.begin (), idx = 0; i_sub != v_sub.end (); ++i_sub, ++idx) {
		R.init (a, (idx + offset) % DEFAULT_CYCLE_LEN);
		if (!R.areEqual (a, *i_sub)) {
			error << "ERROR: Error at subvector index " << idx << " (subvector offset " << offset << "; const_iterator access): expected value ";
			R.write (error, a) << " but got ";
			R.write (error, *i_sub) << std::endl;
			pass = false;
		}

		if (!R.areEqual (a, v_sub[idx])) {
			error << "ERROR: Error at subvector index " << idx << " (subvector offset " << offset << "; array access): expected value ";
			R.write (error, a) << " but got ";
			R.write (error, v_sub[idx]) << std::endl;
			pass = false;
		}
	}

	if (idx != DEFAULT_CYCLE_LEN) {
		error << "ERROR: Subvector offset " << offset << " ended at wrong index " << idx << " instead of expected " << DEFAULT_CYCLE_LEN << std::endl;
		pass = false;
	}

	return pass;
}

template <class Vector, class Ring>
bool testDenseSubvector (const Ring &R, Vector &v, typename Vector::SubvectorType &v_sub, size_t offset)
{	
	typename Vector::SubvectorType::const_iterator i_sub;
	typename Vector::SubvectorType::iterator j_sub;

	size_t idx;

	bool pass = true;

	for (j_sub = v_sub.begin (); j_sub != v_sub.end (); ++j_sub)
		R.init (*j_sub, (idx + offset) % DEFAULT_CYCLE_LEN);

	for (i_sub = v_sub.begin (), idx = offset; i_sub != v_sub.end (); ++i_sub, ++offset)
		R.init (a, (idx + offset) % DEFAULT_CYCLE_LEN);
		if (!R.areEqual (a, *i_sub)) {
			error << "ERROR: Error at vector index " << idx << " (subvector offset " << offset << "; iterator access): expected value ";
			R.write (error, a) << " but got ";
			R.write (error, *i_sub) << std::endl;
			pass = false;
		}

		if (!R.areEqual (a, v[idx])) {
			error << "ERROR: Error at vector index " << idx << " (subvector offset " << offset << "; array access): expected value ";
			R.write (error, a) << " but got ";
			R.write (error, v[idx]) << std::endl;
			pass = false;
		}
	}

	if (idx != DEFAULT_CYCLE_LEN) {
		error << "ERROR: Subvector offset " << offset << " ended at wrong index " << idx << " instead of expected " << DEFAULT_CYCLE_LEN << std::endl;
		pass = false;
	}

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Dense)
{
	bool pass = true, part_pass;

	typename Vector::const_iterator i;
	typename Vector::iterator j;
	size_t idx;

	typename Ring::Element a;

	LELA::commentator.start ("Testing vector (dense)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	// Test 1: Array-access vs. const_iterator
	//
	// Initialise the vector with array-access and check the result with the const_iterator
	LELA::commentator.start ("Test 1: Array-access vs. const_iterator", __FUNCTION__);
	part_pass = true;

	Vector v1 (DEFAULT_DIM);

	for (idx = 0; idx < DEFAULT_DIM; ++idx)
		R.init (v1[idx], idx % DEFAULT_CYCLE_LEN);

	for (idx = 0, i = v1.begin (); i != v1.end (); ++i, ++idx) {
		R.init (a, idx % DEFAULT_CYCLE_LEN);
		if (!R.areEqual (a, *i)) {
			error << "ERROR: Error at position " << idx << ": should get ";
			R.write (error, a) << " but instead got ";
			R.write (error, *i) << std::endl;
			part_pass = false;
		}
	}

	if (idx != DEFAULT_DIM) {
		error << "ERROR: Vector const_iterator ended at wrong index: " << idx << " instead of " << DEFAULT_DIM << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 2: iterator vs. array-access
	//
	// Initialise the vector using an iterator and check the result with array-access
	LELA::commentator.start ("Test 2: iterator vs. array-access", __FUNCTION__);
	part_pass = true;

	Vector v2 (DEFAULT_DIM);

	for (j = v2.begin (); j != v2.end (); ++j)
		R.init (*j, idx % DEFAULT_CYCLE_LEN);

	for (idx = 0; idx < DEFAULT_DIM; ++idx) {
		R.init (a, idx % DEFAULT_CYCLE_LEN);
		if (!R.areEqual (a, v2[idx])) {
			error << "ERROR: Error at position " << idx << ": should get ";
			R.write (error, a) << " but instead got ";
			R.write (error, v2[idx]) << std::endl;
			pass = false;
		}
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 3: size and empty
	//
	// Check whether the value reported by size is correct
	LELA::commentator.start ("Test 3: size () and empty ()", __FUNCTION__);
	part_pass = true;

	Vector v3 (10);

	if (v3.size () != 10) {
		error << "ERROR: Vector::size () reports " << v3.size () << " when it should be " << 10 << std::endl;
		part_pass = false;
	}

	if (v3.empty ()) {
		error << "ERROR: Vector::empty () reports true when it should be false (subcase: size 10)" << std::endl;
		part_pass = false;
	}

	Vector v4 (1);

	if (v4.size () != 1) {
		error << "ERROR: Vector::size () reports " << v4.size () << " when it should be " << 1 << std::endl;
		part_pass = false;
	}

	if (v4.empty ()) {
		error << "ERROR: Vector::empty () reports true when it should be false (subcase: size 1)" << std::endl;
		part_pass = false;
	}

	Vector v5 (0);

	if (v5.size () != 0) {
		error << "ERROR: Vector::size () reports " << v5.size () << " when it should be " << 0 << std::endl;
		part_pass = false;
	}

	if (!v5.empty ()) {
		error << "ERROR: Vector::empty () reports false when it should be true (subcase: size 0)" << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 4: const subvector read
	//
	// Initialise elements of a vector and check that a given subvector is correct
	LELA::commentator.start ("Test 4: const subvector read", __FUNCTION__);
	part_pass = true;

	Vector v6 (DEFAULT_DIM);
	typename Vector::ConstSubvectorType v6_sub1 (v6, 0, DEFAULT_CYCLE_LEN);
	typename Vector::ConstSubvectorType v6_sub2 (v6, DEFAULT_MID_SUBVECTOR_START_IDX, DEFAULT_MID_SUBVECTOR_START_IDX + DEFAULT_CYCLE_LEN);
	typename Vector::ConstSubvectorType v6_sub3 (v6, DEFAULT_DIM - DEFAULT_CYCLE_LEN, DEFAULT_DIM);

	part_pass = testDenseConstSubvector (R, v6, v6_sub1, 0) && part_pass;
	part_pass = testDenseConstSubvector (R, v6, v6_sub2, DEFAULT_MID_SUBVECTOR_START_IDX) && part_pass;
	part_pass = testDenseConstSubvector (R, v6, v6_sub3, DEFAULT_DIM - DEFAULT_CYCLE_LEN) && part_pass;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 5: subvector write
	//
	// Initialise a subvector and check that the original vector is correct
	LELA::commentator.start ("Test 5: subvector write", __FUNCTION__);
	part_pass = true;

	Vector v7 (DEFAULT_DIM);
	typename Vector::ConstSubvectorType v7_sub1 (v7, 0, DEFAULT_CYCLE_LEN);
	typename Vector::ConstSubvectorType v7_sub2 (v7, DEFAULT_MID_SUBVECTOR_START_IDX, DEFAULT_MID_SUBVECTOR_START_IDX + DEFAULT_CYCLE_LEN);
	typename Vector::ConstSubvectorType v7_sub3 (v7, DEFAULT_DIM - DEFAULT_CYCLE_LEN, DEFAULT_DIM);

	part_pass = testDenseSubvector (R, v7, v7_sub1, 0) && part_pass;
	part_pass = testDenseSubvector (R, v7, v7_sub2, DEFAULT_MID_SUBVECTOR_START_IDX) && part_pass;
	part_pass = testDenseSubvector (R, v7, v7_sub3, DEFAULT_DIM - DEFAULT_CYCLE_LEN) && part_pass;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
bool buildSparseVector (const Ring &R, Vector &v, size_t &front_idx, typename Ring::Element &f, size_t &back_idx, typename Ring::Element &e) 
{
	size_t idx;
	typename Ring::Element a;

	first_idx = 0xffffffff;

	for (idx = 0, back_idx = 0; idx < DEFAULT_DIM; ++idx, back_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		R.init (e, idx % DEFAULT_CYCLE_LEN);

		if (!R.isZero (a)) {
			v1.push_back (typename Vector::value_type (back_idx, typename Ring::Element ()));
			R.assign (v1.back ().second, e);

			if (first_idx == 0xffffffff) {
				first_idx = back_idx;
				R.assign (f, e);
			}
		}
	}
}

template <class Vector, class Ring>
bool checkSparseVector (const Ring &R, const Vector &v, size_t front_idx, const typename Ring::Element &f, size_t back_idx, const typename Ring::Element &e) 
{
	bool pass = true;

	typename Vector::const_iterator i;
	size_t idx, curr_idx, curr_idx_1;

	typename Ring::Element a, b;

	if (v.front ().first != front_idx) {
		error << "ERROR: Incorrect index of front (): " << v.front ().first << ", expected " << front_idx << std::endl;
		pass = false;
	}

	if (!R.areEqual (v.front ().second, f)) {
		error << "ERROR: Incorrect element in front (): ";
		R.write (error, v.front ().second) << ", expected ";
		R.write (error, f) << std::endl;
		pass = false;
	}

	for (i = v.begin (), idx = curr_idx = curr_idx_1 = 0; i != v.end (); ++i, ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		R.init (a, idx % DEFAULT_CYCLE_LEN);

		if (!R.isZero (a)) {
			if (i->first != curr_idx) {
				error << "ERROR: Index of entry is " << i->first << ", expected " << curr_idx << std::endl;
				pass = false;
			}

			while (curr_idx_1 < i->first) {
				if (VectorUtils::getEntry (v1, b, curr_idx_1)) {
					error << "ERROR: VectorUtils::getEntry obtained an entry for index " << curr_idx_1
					      << " where there should not be one; entry is ";
					R.write (error, b) << std::endl;
					pass = false;
				}

				++curr_idx_1;
			}

			if (!VectorUtils::getEntry (v1, b, curr_idx_1)) {
				error << "ERROR: VectorUtils::getEntry failed to obtain an entry for index " << curr_idx_1
				      << " where there should be one; entry should be ";
				R.write (error, a) << std::endl;
				pass = false;
			}

			if (!R.areEqual (a, i->second)) {
				error << "ERROR: Entry at index " << i->first << " is ";
				R.write (error, i->second) << "; expected ";
				R.write (error, a) << std::endl;
				pass = false;
			}

			if (!R.areEqual (a, b)) {
				error << "ERROR: Entry at index " << curr_idx << " obtained by VectorUtils::getEntry is ";
				R.write (error, b) << "; expected ";
				R.write (error, a) << std::endl;
				pass = false;
			}
		}
	}

	if (back_idx != curr_idx) {
		error << "ERROR: Index at v.end () - 1 is " << curr_idx << " but expected " << back_idx;
	}

	if (v.back ().first != back_idx) {
		error << "ERROR: Incorrect index of back (): " << v.back ().first << ", expected " << back_idx << std::endl;
		pass = false;
	}

	if (!R.areEqual (v.back ().second, e)) {
		error << "ERROR: Incorrect element in back (): ";
		R.write (error, v.back ().second) << ", expected ";
		R.write (error, e) << std::endl;
		pass = false;
	}

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Sparse)
{
	bool pass = true, part_pass;

	typename Vector::const_iterator i;
	typename Vector::iterator j;
	size_t idx, end_idx, first_idx, exp_size;

	typename Ring::Element a, e;

	Context<Ring> ctx (R);

	LELA::commentator.start ("Testing vector (sparse)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	// Test 1: push_back, const_iterator, VectorUtils::getEntry
	//
	// Construct a sparse vector with push_back and check that the results are the same with const_iterator and VectorUtils::getEntry
	LELA::commentator.start ("Test 1: push_back vs. const_iterator", __FUNCTION__);
	part_pass = true;

	Vector v1;

	buildSparseVector (R, v1, first_idx, f, end_idx, e);

	report << "Vector as constructed: ";
	BLAS1::write (ctx, report, v1) << std::endl;

	part_pass = checkSparseVector (R, v1, first_idx, f, end_idx, e);

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 2: clear, empty, and size
	//
	// Construct a sparse vector, check its size, clear it, and check that it is empty
	LELA::commentator.start ("Test 2: clear, empty, size", __FUNCTION__);
	part_pass = true;

	Vector v2;

	for (idx = 0, curr_idx = 0, exp_size = 0; idx < DEFAULT_DIM; ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		R.init (a, idx % DEFAULT_CYCLE_LEN);

		if (!R.isZero (a)) {
			v1.push_back (typename Vector::value_type (curr_idx, typename Ring::Element ()));
			R.assign (v1.back ().second, a);
			++exp_size;
		}
	}

	report << "Vector as constructed: ";
	BLAS1::write (ctx, report, v2) << std::endl;

	if (v2.size () != exp_size) {
		error << "ERROR: v2.size () reports " << v2.size () << ", expected " << exp_size << std::endl;
		part_pass = false;
	}

	if (exp_size > 0 && v2.empty ()) {
		error << "ERROR: v2.empty () reports true expected false" << std::endl;
		part_pass = false;
	}
	else if (exp_size == 0)
		LELA::commentator.report (LELA::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "WARNING: Constructed empty vector for test 2; check parameters used" << std::endl;

	v2.clear ();

	if (!v2.empty ()) {
		error << "ERROR: v2.empty () is false after call to clear ()" << std::endl;
		part_pass = false;
	}
	else if (v2.size () != 0) {
		error << "ERROR: v2.size () != 0 after call to clear ()" << std::endl;
		part_pass = false;
	}
	else if (v2.begin () != v2.end ()) {
		error << "ERROR: v2.begin () != v2.end () after call to clear ()" << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 3: erase
	//
	// Construct a sparse vector then erase some entries
	LELA::commentator.start ("Test 2: clear, empty, size", __FUNCTION__);
	part_pass = true;

	Vector v3;

	v3.push_back (typename Vector::value_type (10, typename Ring::Element ()));
	R.init (v3.back ().second, 65521);

	v3.push_back (typename Vector::value_type (20, typename Ring::Element ()));
	R.init (v3.back ().second, 65521);

	v3.push_back (typename Vector::value_type (30, typename Ring::Element ()));
	R.init (v3.back ().second, 65521);

	v3.push_back (typename Vector::value_type (40, typename Ring::Element ()));
	R.init (v3.back ().second, 65521);

	report << "Vector as constructed: ";
	BLAS1::write (ctx, report, v2) << std::endl;

	// Erase the second entry
	j = v3.begin ();
	++j;

	j = v3.erase (j);

	if (v3.size () != 3) {
		error << "ERROR: v3.size () is " << v3.size () << " after call to erase (), expected 3" << std::endl;
		part_pass = false;
	}

	if (j->first != 30) {
		error << "ERROR: Return-value of erase () appears to point to wrong place; index is " << j->first << ", expected 30" << std::endl;
		part_pass = false;
	}

	// Erase the last entry
	j = v3.begin ();
	++j; ++j;

	v3.erase (j);

	if (v3.back ().first != 20) {
		error << "ERROR: v3.back ().first = " << v3.back ().first << " after erasure of last entry, expected 20" << std::endl;
		part_pass = false;
	}

	// Erase the first entry
	j = v3.begin ();
	v3.erase (j);

	if (v3.front ().first != 20) {
		error << "ERROR: v3.front ().first = " << v3.front ().first << " after erasure of first entry, expected 20" << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 3: assign
	//
	// Construct two sparse vectors, assign the second to the first
	LELA::commentator.start ("Test 3: assign", __FUNCTION__);
	part_pass = true;

	Vector v4;

	v4.push_back (typename Vector::value_type (10, typename Ring::Element ()));
	R.init (v4.back ().second, 65521);

	v4.push_back (typename Vector::value_type (20, typename Ring::Element ()));
	R.init (v4.back ().second, 65521);

	v4.push_back (typename Vector::value_type (30, typename Ring::Element ()));
	R.init (v4.back ().second, 65521);

	v4.push_back (typename Vector::value_type (40, typename Ring::Element ()));
	R.init (v4.back ().second, 65521);

	report << "Vector v4 as constructed: ";
	BLAS1::write (ctx, report, v4) << std::endl;

	report << "Vector v1 as constructed: ";
	BLAS1::write (ctx, report, v1) << std::endl;

	v4.assign (v1.begin (), v1.end ());

	if (v4.size () != v1.size ()) {
		error << "ERROR: after assign v4.size () = " << v4.size () << ", expected " << v1.size () << std::endl;
		part_pass = false;
	}

	part_pass = checkSparseVector (R, v4, first_idx, f, end_idx, e) && part_pass;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 3: Subvector read
	//
	// Construct a sparse vector and read from a subvector thereof
	LELA::commentator.start ("Test 3: assign", __FUNCTION__);
	part_pass = true;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 4: Subvector push_back
	//
	// Construct a subvector and call push_back to it
	LELA::commentator.start ("Test 3: assign", __FUNCTION__);
	part_pass = true;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 4: Subvector assign
	//
	// Construct a subvector and call assign to it
	LELA::commentator.start ("Test 3: assign", __FUNCTION__);
	part_pass = true;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Dense01)
{
	bool pass = true;

	LELA::commentator.start ("Testing vector (dense 0-1)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Sparse01)
{
	bool pass = true;

	LELA::commentator.start ("Testing vector (sparse 0-1)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Hybrid01)
{
	bool pass = true;

	LELA::commentator.start ("Testing vector (hybrid 0-1)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

/// Run a standard set of tests on the given vector-type
template <class Vector, class Ring>
bool testVector (const Ring &R)
	{ return testVectorSpec (typename VectorTraits<Vector, Ring>::RepresentationType (R)); }

#endif // __LELA_TESTS_TEST_VECTOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
