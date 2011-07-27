/* lela/tests/test-vector.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic tests for vectors
 * 
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_TESTS_TEST_VECTOR_H
#define __LELA_TESTS_TEST_VECTOR_H

#include <iostream>

#include "lela/util/commentator.h"
#include "lela/vector/traits.h"
#include "lela/blas/level1.h"

// Some parameters for the tests
#define DEFAULT_DIM 50
#define DEFAULT_CYCLE_LEN 10
#define DEFAULT_GAP_CYCLE_LEN 5
#define DEFAULT_MID_SUBVECTOR_START_IDX 10

template <class Vector, class Ring>
bool testDenseConstSubvector (const Ring &R, Vector &v, typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType &v_sub, size_t offset)
{
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType::const_iterator i_sub;
	typename Vector::iterator j;

	typename Ring::Element a;

	size_t idx;

	bool pass = true;

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	for (j = v.begin (), idx = 0; j != v.end (); ++j, ++idx)
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
bool testDenseSubvector (const Ring &R, Vector &v, typename LELA::VectorTraits<Ring, Vector>::SubvectorType &v_sub, size_t offset)
{	
	typename LELA::VectorTraits<Ring, Vector>::SubvectorType::const_iterator i_sub;
	typename LELA::VectorTraits<Ring, Vector>::SubvectorType::iterator j_sub;

	typename Ring::Element a;

	size_t idx;

	bool pass = true;

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	for (j_sub = v_sub.begin (), idx = 0; j_sub != v_sub.end (); ++j_sub, ++idx)
		R.init (*j_sub, (idx + offset) % DEFAULT_CYCLE_LEN);

	for (i_sub = v_sub.begin (), idx = 0; i_sub != v_sub.end (); ++i_sub, ++idx) {
		R.init (a, (idx + offset) % DEFAULT_CYCLE_LEN);
		if (!R.areEqual (a, *i_sub)) {
			error << "ERROR: Error at vector index " << idx << " (subvector offset " << offset << "; iterator access): expected value ";
			R.write (error, a) << " but got ";
			R.write (error, *i_sub) << std::endl;
			pass = false;
		}

		if (!R.areEqual (a, v[idx + offset])) {
			error << "ERROR: Error at vector index " << idx << " (subvector offset " << offset << "; array access): expected value ";
			R.write (error, a) << " but got ";
			R.write (error, v[idx + offset]) << std::endl;
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

	for (j = v2.begin (), idx = 0; j != v2.end (); ++j, ++idx)
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
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v6_sub1 (v6, 0, DEFAULT_CYCLE_LEN);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v6_sub2 (v6, DEFAULT_MID_SUBVECTOR_START_IDX, DEFAULT_MID_SUBVECTOR_START_IDX + DEFAULT_CYCLE_LEN);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v6_sub3 (v6, DEFAULT_DIM - DEFAULT_CYCLE_LEN, DEFAULT_DIM);

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
	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v7_sub1 (v7, 0, DEFAULT_CYCLE_LEN);
	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v7_sub2 (v7, DEFAULT_MID_SUBVECTOR_START_IDX, DEFAULT_MID_SUBVECTOR_START_IDX + DEFAULT_CYCLE_LEN);
	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v7_sub3 (v7, DEFAULT_DIM - DEFAULT_CYCLE_LEN, DEFAULT_DIM);

	part_pass = testDenseSubvector (R, v7, v7_sub1, 0) && part_pass;
	part_pass = testDenseSubvector (R, v7, v7_sub2, DEFAULT_MID_SUBVECTOR_START_IDX) && part_pass;
	part_pass = testDenseSubvector (R, v7, v7_sub3, DEFAULT_DIM - DEFAULT_CYCLE_LEN) && part_pass;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
void buildSparseVector (const Ring &R, Vector &v, size_t &front_idx, typename Ring::Element &f, size_t &back_idx, typename Ring::Element &e) 
{
	size_t idx, curr_idx;

	front_idx = 0xffffffff;

	for (idx = 0, curr_idx = 0; idx < DEFAULT_DIM; ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		R.init (e, idx % DEFAULT_CYCLE_LEN);

		if (!R.isZero (e)) {
			v.push_back (typename Vector::value_type (curr_idx, typename Ring::Element ()));
			R.copy (v.back ().second, e);

			if (front_idx == 0xffffffff) {
				front_idx = curr_idx;
				R.copy (f, e);
			}

			back_idx = curr_idx;
		}
	}
}

template <class Vector, class Ring>
bool checkSparseVector (const Ring &R, const Vector &v, size_t front_idx, const typename Ring::Element &f, size_t back_idx, const typename Ring::Element &e) 
{
	bool pass = true;

	typename Vector::const_iterator i;
	size_t idx, curr_idx, curr_idx_1, v_back_idx;

	typename Ring::Element a, b;

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

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

	for (i = v.begin (), idx = curr_idx = curr_idx_1 = 0; i != v.end (); ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		R.init (a, idx % DEFAULT_CYCLE_LEN);

		if (!R.isZero (a)) {
			if (i->first != curr_idx) {
				error << "ERROR: Index of entry is " << i->first << ", expected " << curr_idx << std::endl;
				pass = false;
			}

			while (curr_idx_1 < curr_idx) {
				if (LELA::VectorUtils::getEntry (v, b, curr_idx_1)) {
					error << "ERROR: LELA::VectorUtils::getEntry obtained an entry for index " << curr_idx_1
					      << " where there should not be one; entry is ";
					R.write (error, b) << std::endl;
					pass = false;
				}

				++curr_idx_1;
			}

			if (!LELA::VectorUtils::getEntry (v, b, curr_idx_1)) {
				error << "ERROR: LELA::VectorUtils::getEntry failed to obtain an entry for index " << curr_idx_1
				      << " where there should be one; entry should be ";
				R.write (error, a) << std::endl;
				pass = false;
			}

			++curr_idx_1;

			if (!R.areEqual (a, i->second)) {
				error << "ERROR: Entry at index " << i->first << " is ";
				R.write (error, i->second) << "; expected ";
				R.write (error, a) << std::endl;
				pass = false;
			}

			if (!R.areEqual (a, b)) {
				error << "ERROR: Entry at index " << curr_idx << " obtained by LELA::VectorUtils::getEntry is ";
				R.write (error, b) << "; expected ";
				R.write (error, a) << std::endl;
				pass = false;
			}

			v_back_idx = curr_idx;

			++i;
		}
	}

	if (back_idx != v_back_idx) {
		error << "ERROR: Index at v.end () - 1 is " << v_back_idx << " but expected " << back_idx;
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
bool checkSparseAfterSubConstruction (const Ring &R, const Vector &v, bool expect_first_elt)
{
	bool pass = true;

	typename Vector::const_iterator i;

	typename Ring::Element a;

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	i = v.begin ();

	if (expect_first_elt) {
		if (i->first != 10) {
			error << "ERROR: First index is " << i->first << ", expected 10" << std::endl;
			pass = false;
		}

		R.init (a, 1);

		if (!R.areEqual (a, i->second)) {
			error << "ERROR: First element is ";
			R.write (error, i->second) << ", expected ";
			R.write (error, a) << std::endl;
			pass = false;
		}

		++i;
	}

	if (i->first != 12) {
		error << "ERROR: Second index is " << i->first << ", expected 12" << std::endl;
		pass = false;
	}

	R.init (a, 2);

	if (!R.areEqual (a, i->second)) {
		error << "ERROR: Second element is ";
		R.write (error, i->second) << ", expected ";
		R.write (error, a) << std::endl;
		pass = false;
	}

	++i;

	if (i->first != 14) {
		error << "ERROR: Third index is " << i->first << ", expected 14" << std::endl;
		pass = false;
	}

	R.init (a, 3);

	if (!R.areEqual (a, i->second)) {
		error << "ERROR: Third element is ";
		R.write (error, i->second) << ", expected ";
		R.write (error, a) << std::endl;
		pass = false;
	}

	++i;

	if (i->first != 16) {
		error << "ERROR: Fourth index is " << i->first << ", expected 16" << std::endl;
		pass = false;
	}

	R.init (a, 4);

	if (!R.areEqual (a, i->second)) {
		error << "ERROR: Fourth element is ";
		R.write (error, i->second) << ", expected ";
		R.write (error, a) << std::endl;
		pass = false;
	}

	++i;

	if (i->first != 20) {
		error << "ERROR: Fifth index is " << i->first << ", expected 20" << std::endl;
		pass = false;
	}

	R.init (a, 5);

	if (!R.areEqual (a, i->second)) {
		error << "ERROR: Fifth element is ";
		R.write (error, i->second) << ", expected ";
		R.write (error, a) << std::endl;
		pass = false;
	}

	++i;

	if (i != v.end ()) {
		error << "ERROR: iterator does end at correct position" << std::endl;
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
	size_t idx, curr_idx, end_idx, first_idx, exp_size;

	typename Ring::Element a, f, e;

	LELA::Context<Ring> ctx (R);

	LELA::commentator.start ("Testing vector (sparse)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	// Test 1: push_back, const_iterator, LELA::VectorUtils::getEntry
	//
	// Construct a sparse vector with push_back and check that the results are the same with const_iterator and LELA::VectorUtils::getEntry
	LELA::commentator.start ("Test 1: push_back vs. const_iterator", __FUNCTION__);
	part_pass = true;

	Vector v1;

	buildSparseVector (R, v1, first_idx, f, end_idx, e);

	report << "Vector as constructed: ";
	LELA::BLAS1::write (ctx, report, v1) << std::endl;

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
			v2.push_back (typename Vector::value_type (curr_idx, typename Ring::Element ()));
			R.copy (v2.back ().second, a);
			++exp_size;
		}
	}

	report << "Vector as constructed: ";
	LELA::BLAS1::write (ctx, report, v2) << std::endl;

	if (v2.size () != exp_size) {
		error << "ERROR: v2.size () reports " << v2.size () << ", expected " << exp_size << std::endl;
		part_pass = false;
	}

	if (exp_size > 0 && v2.empty ()) {
		error << "ERROR: v2.empty () reports true expected false" << std::endl;
		part_pass = false;
	}
	else if (exp_size == 0)
		LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
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
	LELA::commentator.start ("Test 3: erase", __FUNCTION__);
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
	LELA::BLAS1::write (ctx, report, v3) << std::endl;

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

	if (v3.back ().first != 30) {
		error << "ERROR: v3.back ().first = " << v3.back ().first << " after erasure of last entry, expected 30" << std::endl;
		part_pass = false;
	}

	// Erase the first entry
	j = v3.begin ();
	v3.erase (j);

	if (v3.front ().first != 30) {
		error << "ERROR: v3.front ().first = " << v3.front ().first << " after erasure of first entry, expected 30" << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 4: assign
	//
	// Construct two sparse vectors, assign the second to the first
	LELA::commentator.start ("Test 4: assign", __FUNCTION__);
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
	LELA::BLAS1::write (ctx, report, v4) << std::endl;

	report << "Vector v1 as constructed: ";
	LELA::BLAS1::write (ctx, report, v1) << std::endl;

	v4.assign (v1.begin (), v1.end ());

	if (v4.size () != v1.size ()) {
		error << "ERROR: after assign v4.size () = " << v4.size () << ", expected " << v1.size () << std::endl;
		part_pass = false;
	}

	part_pass = checkSparseVector (R, v4, first_idx, f, end_idx, e) && part_pass;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 5: Subvector read
	//
	// Construct a sparse vector and read from a subvector thereof
	LELA::commentator.start ("Test 5: subvector read", __FUNCTION__);
	part_pass = true;

	Vector v5;

	v5.push_back (typename Vector::value_type (10, typename Ring::Element ()));
	R.init (v5.back ().second, 1);

	v5.push_back (typename Vector::value_type (20, typename Ring::Element ()));
	R.init (v5.back ().second, 2);

	v5.push_back (typename Vector::value_type (30, typename Ring::Element ()));
	R.init (v5.back ().second, 3);

	v5.push_back (typename Vector::value_type (40, typename Ring::Element ()));
	R.init (v5.back ().second, 4);

	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub1 (v5, 0, 10);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub2 (v5, 0, 20);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub3 (v5, 10, 20);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub4 (v5, 10, 25);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub5 (v5, 30, 50);

	typename LELA::VectorTraits<Ring, Vector>::SubvectorType::const_iterator i_sub;

	report << "Vector v5 as constructed: ";
	LELA::BLAS1::write (ctx, report, v5) << std::endl;

	if (!v5_sub1.empty ()) {
		error << "ERROR: Subvector indices 0..10 not empty, reported size " << v5_sub1.size () << ", first index "
		      << v5_sub1.front ().first << ", last index " << v5_sub1.back ().first << std::endl;
		part_pass = false;
	}
	else if (v5_sub1.begin () != v5_sub1.end ()) {
		error << "ERROR: Subvector indices 0..10: begin () != end ()" << std::endl;
		part_pass = false;
	}

	if (v5_sub2.size () != 1) {
		error << "ERROR: Subvector indices 0..20 has wrong size, reported size " << v5_sub2.size () << ", expected 1, last index is " << v5_sub2.size () << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub2.begin ();

	if (v5_sub2.front ().first != 10) {
		error << "ERROR: First index of subvector indices 0..20 is " << v5_sub2.front ().first << " (from front ()), expected 10" << std::endl;
		part_pass = false;
	}

	if (i_sub->first != 10) {
		error << "ERROR: First index of subvector indices 0..20 is " << i_sub->first << " (from *(begin ())), expected 10" << std::endl;
		part_pass = false;
	}

	R.init (a, 1);

	if (!R.areEqual (a, v5_sub2.front ().second)) {
		error << "ERROR: First element of subvector indices 0..20 is ";
		R.write (error, v5_sub2.front ().second) << " (from front ()), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	if (!R.areEqual (a, i_sub->second)) {
		error << "ERROR: First element of subvector indices 0..20 is ";
		R.write (error, v5_sub2.front ().second) << " (from *(begin ())), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub2.end ()) {
		error << "ERROR: Iterator on subvector indices 0..20 does not end at correct position." << std::endl;
		part_pass = false;
	}

	if (v5_sub3.size () != 1) {
		error << "ERROR: Subvector indices 10..20 has wrong size, reported size " << v5_sub3.size () << ", expected 1, last index is " << v5_sub3.size () << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub3.begin ();

	if (v5_sub3.front ().first != 0) {
		error << "ERROR: First index of subvector indices 10..20 is " << v5_sub3.front ().first << " (from front ()), expected 0" << std::endl;
		part_pass = false;
	}

	if (i_sub->first != 0) {
		error << "ERROR: First index of subvector indices 10..20 is " << i_sub->first << " (from *(begin ())), expected 0" << std::endl;
		part_pass = false;
	}

	R.init (a, 1);

	if (!R.areEqual (a, v5_sub3.front ().second)) {
		error << "ERROR: First element of subvector indices 10..20 is ";
		R.write (error, v5_sub3.front ().second) << " (from front ()), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	if (!R.areEqual (a, i_sub->second)) {
		error << "ERROR: First element of subvector indices 10..20 is ";
		R.write (error, v5_sub3.front ().second) << " (from *(begin ())), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub3.end ()) {
		error << "ERROR: Iterator on subvector indices 10..20 does not end at correct position." << std::endl;
		part_pass = false;
	}

	if (v5_sub4.size () != 2) {
		error << "ERROR: Subvector indices 10..25 has wrong size, reported size " << v5_sub4.size () << ", expected 2" << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub4.begin ();

	if (v5_sub4.front ().first != 0) {
		error << "ERROR: First index of subvector indices 10..25 is " << v5_sub4.front ().first << " (from front ()), expected 0" << std::endl;
		part_pass = false;
	}

	if (i_sub->first != 0) {
		error << "ERROR: First index of subvector indices 10..25 is " << i_sub->first << " (from *(begin ())), expected 0" << std::endl;
		part_pass = false;
	}

	R.init (a, 1);

	if (!R.areEqual (a, v5_sub4.front ().second)) {
		error << "ERROR: First element of subvector indices 10..25 is ";
		R.write (error, v5_sub4.front ().second) << " (from front ()), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	if (!R.areEqual (a, i_sub->second)) {
		error << "ERROR: First element of subvector indices 10..25 is ";
		R.write (error, i_sub->second) << " (from *(begin ())), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (v5_sub4.back ().first != 10) {
		error << "ERROR: Last index of subvector indices 10..25 is " << v5_sub4.back ().first << " (from back ()), expected 10" << std::endl;
		part_pass = false;
	}

	if (i_sub->first != 10) {
		error << "ERROR: Last index of subvector indices 10..25 is " << i_sub->first << " (from iterator), expected 10" << std::endl;
		part_pass = false;
	}

	R.init (a, 2);

	if (!R.areEqual (a, v5_sub4.back ().second)) {
		error << "ERROR: Last element of subvector indices 10..25 is ";
		R.write (error, v5_sub4.back ().second) << " (from back ()), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	if (!R.areEqual (a, i_sub->second)) {
		error << "ERROR: Last element of subvector indices 10..25 is ";
		R.write (error, i_sub->second) << " (from iterator), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub4.end ()) {
		error << "ERROR: Iterator on subvector indices 10..25 does not end at correct position." << std::endl;
		part_pass = false;
	}

	if (v5_sub5.size () != 2) {
		error << "ERROR: Subvector indices 30..50 has wrong size, reported size " << v5_sub5.size () << ", expected 2" << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub5.begin ();

	if (v5_sub5.front ().first != 0) {
		error << "ERROR: First index of subvector indices 30..50 is " << v5_sub5.front ().first << " (from front ()), expected 0" << std::endl;
		part_pass = false;
	}

	if (i_sub->first != 0) {
		error << "ERROR: First index of subvector indices 30..50 is " << i_sub->first << " (from iterator), expected 0" << std::endl;
		part_pass = false;
	}

	R.init (a, 3);

	if (!R.areEqual (a, v5_sub5.front ().second)) {
		error << "ERROR: First element of subvector indices 30..50 is ";
		R.write (error, v5_sub5.front ().second) << " (from front ()), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	if (!R.areEqual (a, i_sub->second)) {
		error << "ERROR: First element of subvector indices 30..50 is ";
		R.write (error, i_sub->second) << " (from iterator), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (v5_sub5.back ().first != 10) {
		error << "ERROR: Last index of subvector indices 30..50 is " << v5_sub5.back ().first << " (from back ()), expected 10" << std::endl;
		part_pass = false;
	}

	if (i_sub->first != 10) {
		error << "ERROR: Last index of subvector indices 30..50 is " << i_sub->first << " (from iterator), expected 10" << std::endl;
		part_pass = false;
	}

	R.init (a, 4);

	if (!R.areEqual (a, v5_sub5.back ().second)) {
		error << "ERROR: Last element of subvector indices 30..50 is ";
		R.write (error, v5_sub5.front ().second) << " (from back ()), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	if (!R.areEqual (a, i_sub->second)) {
		error << "ERROR: Last element of subvector indices 30..50 is ";
		R.write (error, i_sub->second) << " (from iterator), expected ";
		R.write (error, a) << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub5.end ()) {
		error << "ERROR: Iterator on subvector indices 30..50 does not end at correct position." << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 6: Subvector push_back
	//
	// Construct a subvector and call push_back to it
	LELA::commentator.start ("Test 6: subvector push_back", __FUNCTION__);
	part_pass = true;

	Vector v6;

	v6.push_back (typename Vector::value_type (10, typename Ring::Element ()));
	R.init (v6.back ().second, 1);

	v6.push_back (typename Vector::value_type (20, typename Ring::Element ()));
	R.init (v6.back ().second, 5);

	report << "Vector v6 as constructed: ";
	LELA::BLAS1::write (ctx, report, v6) << std::endl;

	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v6_sub (v6, 10, 20);

	v6_sub.push_back (typename Vector::value_type (2, typename Ring::Element ()));
	R.init (v6_sub.back ().second, 2);

	v6_sub.push_back (typename Vector::value_type (4, typename Ring::Element ()));
	R.init (v6_sub.back ().second, 3);

	v6_sub.push_back (typename Vector::value_type (6, typename Ring::Element ()));
	R.init (v6_sub.back ().second, 4);

	report << "Subvector v6_sub as constructed: ";
	LELA::BLAS1::write (ctx, report, v6_sub) << std::endl;

	report << "Vector v6 after construction of subvector: ";
	LELA::BLAS1::write (ctx, report, v6) << std::endl;

	part_pass = checkSparseAfterSubConstruction (R, v6, true);

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 6: Subvector assign
	//
	// Construct a subvector and call assign to it
	LELA::commentator.start ("Test 3: subvector assign", __FUNCTION__);
	part_pass = true;

	Vector v7, v8;

	v7.push_back (typename Vector::value_type (10, typename Ring::Element ()));
	R.init (v7.back ().second, 1);

	v7.push_back (typename Vector::value_type (20, typename Ring::Element ()));
	R.init (v7.back ().second, 5);

	report << "Vector v7 as constructed: ";
	LELA::BLAS1::write (ctx, report, v7) << std::endl;

	v8.push_back (typename Vector::value_type (2, typename Ring::Element ()));
	R.init (v8.back ().second, 2);

	v8.push_back (typename Vector::value_type (4, typename Ring::Element ()));
	R.init (v8.back ().second, 3);

	v8.push_back (typename Vector::value_type (6, typename Ring::Element ()));
	R.init (v8.back ().second, 4);

	report << "Vector v8 as constructed: ";
	LELA::BLAS1::write (ctx, report, v8) << std::endl;

	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v7_sub (v7, 10, 20);

	v7_sub.assign (v8.begin (), v8.end ());

	report << "Vector v7 after assign: ";
	LELA::BLAS1::write (ctx, report, v7) << std::endl;

	part_pass = checkSparseAfterSubConstruction (R, v7, false);

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Dense01)
{
	bool pass = true, part_pass;

	LELA::commentator.start ("Testing vector (dense 0-1)", __FUNCTION__);

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	// First we run all the normal tests on dense vectors
	testVectorSpec<Vector> (R, LELA::VectorRepresentationTypes::Dense ());

	// std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	// Test 1: word_iterator and back_word vs iterator
	//
	// Build the vector with word_iterator and back_word and check that the result with iterator is as expected
	LELA::commentator.start ("Test 1: word_iterator and back_word vs. iterator", __FUNCTION__);
	part_pass = true;

	Vector v1 (3 * LELA::WordTraits<typename Vector::word_type>::bits);

	typename Vector::word_iterator i_w;

	i_w = v1.word_begin ();
	*i_w = 0xffffffffffffffffULL;
	++i_w;
	*i_w = 0xaaaaaaaaaaaaaaaaULL;
	++i_w;

	if (i_w != v1.word_end ()) {
		error << "ERROR: word_iterator did not end where expected" << std::endl;
		part_pass = false;
	}

	v1.back_word () = 0x8888888888888888ULL;

	typename Vector::const_iterator i = v1.begin ();

	size_t idx;

	for (idx = 0; idx < LELA::WordTraits<typename Vector::word_type>::bits; ++idx, ++i) {
		if (!*i) {
			error << "ERROR: iterator reports 0 at position " << idx << "; expected 1" << std::endl;
			part_pass = false;
		}
	}

	for (idx = 0; idx < LELA::WordTraits<typename Vector::word_type>::bits; ++idx, ++i) {
		if (Vector::Endianness::e_0 & 1) {
			if ((idx & 1) && !*i) {
				error << "ERROR: iterator reports 0 at position " << idx << "; expected 1" << std::endl;
				part_pass = false;
			}
			else if (!(idx & 1) && *i) {
				error << "ERROR: iterator reports 1 at position " << idx << "; expected 0" << std::endl;
				part_pass = false;
			}
		} else {
			if (!(idx & 1) && !*i) {
				error << "ERROR: iterator reports 0 at position " << idx << "; expected 1" << std::endl;
				part_pass = false;
			}
			else if ((idx & 1) && *i) {
				error << "ERROR: iterator reports 1 at position " << idx << "; expected 0" << std::endl;
				part_pass = false;
			}
		}
	}

	for (idx = 0; idx < LELA::WordTraits<typename Vector::word_type>::bits; ++idx, ++i) {
		if (Vector::Endianness::e_0 & 1) {
			if ((idx % 4 == 3) && !*i) {
				error << "ERROR: iterator reports 0 at position " << idx << "; expected 1" << std::endl;
				part_pass = false;
			}
			else if ((idx % 4 != 3) && *i) {
				error << "ERROR: iterator reports 1 at position " << idx << "; expected 0" << std::endl;
				part_pass = false;
			}
		} else {
			if ((idx % 4 == 0) && !*i) {
				error << "ERROR: iterator reports 0 at position " << idx << "; expected 1" << std::endl;
				part_pass = false;
			}
			else if ((idx % 4 != 0) && *i) {
				error << "ERROR: iterator reports 1 at position " << idx << "; expected 0" << std::endl;
				part_pass = false;
			}
		}
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 2: iterator vs. word_iterator and back_word
	//
	// Reverse of test 1

	LELA::commentator.start ("Test 2: iterator vs. word_iterator and back_word", __FUNCTION__);
	part_pass = true;

	Vector v2 (3 * LELA::WordTraits<typename Vector::word_type>::bits);

	for (idx = 0; idx < LELA::WordTraits<typename Vector::word_type>::bits; ++idx)
		v2[idx] = true;

	for (; idx < 2 * LELA::WordTraits<typename Vector::word_type>::bits; ++idx) {
		if (Vector::Endianness::e_0 & 1)
			v2[idx] = ((idx & 1) != 0);
		else
			v2[idx] = ((idx & 1) == 0);
	}

	for (; idx < 3 * LELA::WordTraits<typename Vector::word_type>::bits; ++idx) {
		if (Vector::Endianness::e_0 & 1)
			v2[idx] = (idx % 4 == 3);
		else
			v2[idx] = (idx % 4 == 0);
	}

	typename Vector::const_word_iterator j_w = v2.word_begin ();

	typename Vector::word_type w = 0xffffffffffffffffULL;

	if (*j_w != w) {
		error << "ERROR: First word is " << std::hex << *j_w << ", expected " << w << std::endl;
		part_pass = false;
	}

	++j_w;
	w = 0xaaaaaaaaaaaaaaaaULL;

	if (*j_w != w) {
		error << "ERROR: Second word is " << std::hex << *j_w << ", expected " << w << std::endl;
		part_pass = false;
	}

	++j_w;

	if (j_w != v2.word_end ()) {
		error << "ERROR: const_word_iterator did not end at expected point" << std::endl;
		part_pass = false;
	}

	w = 0x8888888888888888ULL;

	if (v2.back_word () != w) {
		error << "ERROR: Last word is " << std::hex << v2.back_word () << ", expected " << w << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 3: word_size
	//
	// Check that word_size gives the right answer

	LELA::commentator.start ("Test 3: word_size", __FUNCTION__);
	part_pass = true;

	Vector v3 (4 * LELA::WordTraits<typename Vector::word_type>::bits);
	Vector v4 (4 * LELA::WordTraits<typename Vector::word_type>::bits + 10);

	if (v3.word_size () != 4) {
		error << "ERROR: word_size of v3 is " << v3.word_size () << ", expected 4" << std::endl;
		part_pass = false;
	}

	if (v4.word_size () != 5) {
		error << "ERROR: word_size of v4 is " << v4.word_size () << ", expected 5" << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Note: bit-subvectors are tested in test-bit-subvector.C

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
void buildSparse01Vector (const Ring &R, Vector &v, size_t &front_idx, size_t &back_idx) 
{
	size_t idx, curr_idx;

	front_idx = 0xffffffff;

	for (idx = 0, curr_idx = 0; idx < DEFAULT_DIM; ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		v.push_back (curr_idx);

		if (front_idx == 0xffffffff)
			front_idx = curr_idx;

		back_idx = curr_idx;
	}
}

template <class Vector, class Ring>
bool checkSparse01Vector (const Ring &R, const Vector &v, size_t front_idx, size_t back_idx) 
{
	bool pass = true;

	typename Vector::const_iterator i;
	size_t idx, curr_idx, curr_idx_1, v_back_idx;

	typename Ring::Element b;

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	if (v.front () != front_idx) {
		error << "ERROR: Incorrect index of front (): " << v.front () << ", expected " << front_idx << std::endl;
		pass = false;
	}

	for (i = v.begin (), idx = curr_idx = curr_idx_1 = 0; i != v.end (); ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		if (*i != curr_idx) {
			error << "ERROR: Index of entry is " << *i << ", expected " << curr_idx << std::endl;
			pass = false;
		}

		while (curr_idx_1 < curr_idx) {
			if (LELA::VectorUtils::getEntry (v, b, curr_idx_1)) {
				error << "ERROR: LELA::VectorUtils::getEntry obtained an entry for index " << curr_idx_1
				      << " where there should not be one" << std::endl;
				pass = false;
			}

			++curr_idx_1;
		}

		if (!LELA::VectorUtils::getEntry (v, b, curr_idx_1)) {
			error << "ERROR: LELA::VectorUtils::getEntry failed to obtain an entry for index " << curr_idx_1
			      << " where there should be one" << std::endl;
			pass = false;
		}

		++curr_idx_1;

		v_back_idx = curr_idx;

		++i;
	}

	if (back_idx != v_back_idx) {
		error << "ERROR: Index at v.end () - 1 is " << v_back_idx << " but expected " << back_idx;
	}

	if (v.back () != back_idx) {
		error << "ERROR: Incorrect index of back (): " << v.back () << ", expected " << back_idx << std::endl;
		pass = false;
	}

	return pass;
}

template <class Vector, class Ring>
bool checkSparse01AfterSubConstruction (const Ring &R, const Vector &v, bool expect_first_elt)
{
	bool pass = true;

	typename Vector::const_iterator i;

	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	i = v.begin ();

	if (expect_first_elt) {
		if (*i != 10) {
			error << "ERROR: First index is " << *i << ", expected 10" << std::endl;
			pass = false;
		}

		++i;
	}

	if (*i != 12) {
		error << "ERROR: Second index is " << *i << ", expected 12" << std::endl;
		pass = false;
	}

	++i;

	if (*i != 14) {
		error << "ERROR: Third index is " << *i << ", expected 14" << std::endl;
		pass = false;
	}

	++i;

	if (*i != 16) {
		error << "ERROR: Fourth index is " << *i << ", expected 16" << std::endl;
		pass = false;
	}

	++i;

	if (*i != 20) {
		error << "ERROR: Fifth index is " << *i << ", expected 20" << std::endl;
		pass = false;
	}

	++i;

	if (i != v.end ()) {
		error << "ERROR: iterator does end at correct position" << std::endl;
		pass = false;
	}

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Sparse01)
{
	bool pass = true, part_pass;

	typename Vector::const_iterator i;
	typename Vector::iterator j;
	size_t idx, curr_idx, end_idx, first_idx, exp_size;

	LELA::Context<Ring> ctx (R);

	LELA::commentator.start ("Testing vector (sparse 0-1)", __FUNCTION__);

	std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_ERROR);

	// Test 1: push_back, const_iterator, LELA::VectorUtils::getEntry
	//
	// Construct a sparse vector with push_back and check that the results are the same with const_iterator and LELA::VectorUtils::getEntry
	LELA::commentator.start ("Test 1: push_back vs. const_iterator", __FUNCTION__);
	part_pass = true;

	Vector v1;

	buildSparse01Vector (R, v1, first_idx, end_idx);

	report << "Vector as constructed: ";
	LELA::BLAS1::write (ctx, report, v1) << std::endl;

	part_pass = checkSparse01Vector (R, v1, first_idx, end_idx);

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 2: clear, empty, and size
	//
	// Construct a sparse vector, check its size, clear it, and check that it is empty
	LELA::commentator.start ("Test 2: clear, empty, size", __FUNCTION__);
	part_pass = true;

	Vector v2;

	for (idx = 0, curr_idx = 0, exp_size = 0; idx < DEFAULT_DIM; ++idx, curr_idx += idx % DEFAULT_GAP_CYCLE_LEN + 1) {
		v2.push_back (curr_idx);
		++exp_size;
	}

	report << "Vector as constructed: ";
	LELA::BLAS1::write (ctx, report, v2) << std::endl;

	if (v2.size () != exp_size) {
		error << "ERROR: v2.size () reports " << v2.size () << ", expected " << exp_size << std::endl;
		part_pass = false;
	}

	if (exp_size > 0 && v2.empty ()) {
		error << "ERROR: v2.empty () reports true expected false" << std::endl;
		part_pass = false;
	}
	else if (exp_size == 0)
		LELA::commentator.report (LELA::Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
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
	LELA::commentator.start ("Test 3: erase", __FUNCTION__);
	part_pass = true;

	Vector v3;

	v3.push_back (10);
	v3.push_back (20);
	v3.push_back (30);
	v3.push_back (40);

	report << "Vector as constructed: ";
	LELA::BLAS1::write (ctx, report, v3) << std::endl;

	// Erase the second entry
	j = v3.begin ();
	++j;

	j = v3.erase (j);

	if (v3.size () != 3) {
		error << "ERROR: v3.size () is " << v3.size () << " after call to erase (), expected 3" << std::endl;
		part_pass = false;
	}

	if (*j != 30) {
		error << "ERROR: Return-value of erase () appears to point to wrong place; index is " << *j << ", expected 30" << std::endl;
		part_pass = false;
	}

	// Erase the last entry
	j = v3.begin ();
	++j; ++j;

	v3.erase (j);

	if (v3.back () != 30) {
		error << "ERROR: v3.back () = " << v3.back () << " after erasure of last entry, expected 30" << std::endl;
		part_pass = false;
	}

	// Erase the first entry
	j = v3.begin ();
	v3.erase (j);

	if (v3.front () != 30) {
		error << "ERROR: v3.front () = " << v3.front () << " after erasure of first entry, expected 30" << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 4: assign
	//
	// Construct two sparse vectors, assign the second to the first
	LELA::commentator.start ("Test 4: assign", __FUNCTION__);
	part_pass = true;

	Vector v4;

	v4.push_back (10);
	v4.push_back (20);
	v4.push_back (30);
	v4.push_back (40);

	report << "Vector v4 as constructed: ";
	LELA::BLAS1::write (ctx, report, v4) << std::endl;

	report << "Vector v1 as constructed: ";
	LELA::BLAS1::write (ctx, report, v1) << std::endl;

	v4.assign (v1.begin (), v1.end ());

	if (v4.size () != v1.size ()) {
		error << "ERROR: after assign v4.size () = " << v4.size () << ", expected " << v1.size () << std::endl;
		part_pass = false;
	}

	part_pass = checkSparse01Vector (R, v4, first_idx, end_idx) && part_pass;

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 5: Subvector read
	//
	// Construct a sparse vector and read from a subvector thereof
	LELA::commentator.start ("Test 5: subvector read", __FUNCTION__);
	part_pass = true;

	Vector v5;

	v5.push_back (10);
	v5.push_back (20);
	v5.push_back (30);
	v5.push_back (40);

	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub1 (v5, 0, 10);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub2 (v5, 0, 20);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub3 (v5, 10, 20);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub4 (v5, 10, 25);
	typename LELA::VectorTraits<Ring, Vector>::ConstSubvectorType v5_sub5 (v5, 30, 50);

	typename LELA::VectorTraits<Ring, Vector>::SubvectorType::const_iterator i_sub;

	report << "Vector v5 as constructed: ";
	LELA::BLAS1::write (ctx, report, v5) << std::endl;

	if (!v5_sub1.empty ()) {
		error << "ERROR: Subvector indices 0..10 not empty, reported size " << v5_sub1.size () << ", first index "
		      << v5_sub1.front () << ", last index " << v5_sub1.back () << std::endl;
		part_pass = false;
	}
	else if (v5_sub1.begin () != v5_sub1.end ()) {
		error << "ERROR: Subvector indices 0..10: begin () != end ()" << std::endl;
		part_pass = false;
	}

	if (v5_sub2.size () != 1) {
		error << "ERROR: Subvector indices 0..20 has wrong size, reported size " << v5_sub2.size () << ", expected 1, last index is " << v5_sub2.size () << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub2.begin ();

	if (v5_sub2.front () != 10) {
		error << "ERROR: First index of subvector indices 0..20 is " << v5_sub2.front () << " (from front ()), expected 10" << std::endl;
		part_pass = false;
	}

	if (*i_sub != 10) {
		error << "ERROR: First index of subvector indices 0..20 is " << *i_sub << " (from *(begin ())), expected 10" << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub2.end ()) {
		error << "ERROR: Iterator on subvector indices 0..20 does not end at correct position." << std::endl;
		part_pass = false;
	}

	if (v5_sub3.size () != 1) {
		error << "ERROR: Subvector indices 10..20 has wrong size, reported size " << v5_sub3.size () << ", expected 1, last index is " << v5_sub3.size () << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub3.begin ();

	if (v5_sub3.front () != 0) {
		error << "ERROR: First index of subvector indices 10..20 is " << v5_sub3.front () << " (from front ()), expected 0" << std::endl;
		part_pass = false;
	}

	if (*i_sub != 0) {
		error << "ERROR: First index of subvector indices 10..20 is " << *i_sub << " (from *(begin ())), expected 0" << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub3.end ()) {
		error << "ERROR: Iterator on subvector indices 10..20 does not end at correct position." << std::endl;
		part_pass = false;
	}

	if (v5_sub4.size () != 2) {
		error << "ERROR: Subvector indices 10..25 has wrong size, reported size " << v5_sub4.size () << ", expected 2" << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub4.begin ();

	if (v5_sub4.front () != 0) {
		error << "ERROR: First index of subvector indices 10..25 is " << v5_sub4.front () << " (from front ()), expected 0" << std::endl;
		part_pass = false;
	}

	if (*i_sub != 0) {
		error << "ERROR: First index of subvector indices 10..25 is " << *i_sub << " (from *(begin ())), expected 0" << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (v5_sub4.back () != 10) {
		error << "ERROR: Last index of subvector indices 10..25 is " << v5_sub4.back () << " (from back ()), expected 10" << std::endl;
		part_pass = false;
	}

	if (*i_sub != 10) {
		error << "ERROR: Last index of subvector indices 10..25 is " << *i_sub << " (from iterator), expected 10" << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub4.end ()) {
		error << "ERROR: Iterator on subvector indices 10..25 does not end at correct position." << std::endl;
		part_pass = false;
	}

	if (v5_sub5.size () != 2) {
		error << "ERROR: Subvector indices 30..50 has wrong size, reported size " << v5_sub5.size () << ", expected 2" << std::endl;
		part_pass = false;
	}

	i_sub = v5_sub5.begin ();

	if (v5_sub5.front () != 0) {
		error << "ERROR: First index of subvector indices 30..50 is " << v5_sub5.front () << " (from front ()), expected 0" << std::endl;
		part_pass = false;
	}

	if (*i_sub != 0) {
		error << "ERROR: First index of subvector indices 30..50 is " << *i_sub << " (from iterator), expected 0" << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (v5_sub5.back () != 10) {
		error << "ERROR: Last index of subvector indices 30..50 is " << v5_sub5.back () << " (from back ()), expected 10" << std::endl;
		part_pass = false;
	}

	if (*i_sub != 10) {
		error << "ERROR: Last index of subvector indices 30..50 is " << *i_sub << " (from iterator), expected 10" << std::endl;
		part_pass = false;
	}

	++i_sub;

	if (i_sub != v5_sub5.end ()) {
		error << "ERROR: Iterator on subvector indices 30..50 does not end at correct position." << std::endl;
		part_pass = false;
	}

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 6: Subvector push_back
	//
	// Construct a subvector and call push_back to it
	LELA::commentator.start ("Test 6: subvector push_back", __FUNCTION__);
	part_pass = true;

	Vector v6;

	v6.push_back (10);
	v6.push_back (20);

	report << "Vector v6 as constructed: ";
	LELA::BLAS1::write (ctx, report, v6) << std::endl;

	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v6_sub (v6, 10, 20);

	v6_sub.push_back (2);
	v6_sub.push_back (4);
	v6_sub.push_back (6);

	report << "Subvector v6_sub as constructed: ";
	LELA::BLAS1::write (ctx, report, v6_sub) << std::endl;

	report << "Vector v6 after construction of subvector: ";
	LELA::BLAS1::write (ctx, report, v6) << std::endl;

	part_pass = checkSparse01AfterSubConstruction (R, v6, true);

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	// Test 6: Subvector assign
	//
	// Construct a subvector and call assign to it
	LELA::commentator.start ("Test 3: subvector assign", __FUNCTION__);
	part_pass = true;

	Vector v7, v8;

	v7.push_back (10);
	v7.push_back (20);

	report << "Vector v7 as constructed: ";
	LELA::BLAS1::write (ctx, report, v7) << std::endl;

	v8.push_back (2);
	v8.push_back (4);
	v8.push_back (6);

	report << "Vector v8 as constructed: ";
	LELA::BLAS1::write (ctx, report, v8) << std::endl;

	typename LELA::VectorTraits<Ring, Vector>::SubvectorType v7_sub (v7, 10, 20);

	v7_sub.assign (v8.begin (), v8.end ());

	report << "Vector v7 after assign: ";
	LELA::BLAS1::write (ctx, report, v7) << std::endl;

	part_pass = checkSparse01AfterSubConstruction (R, v7, false);

	LELA::commentator.stop (MSG_STATUS (part_pass));
	pass = part_pass && pass;

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Vector, class Ring>
bool testVectorSpec (const Ring &R, LELA::VectorRepresentationTypes::Hybrid01)
{
	bool pass = true;

	LELA::commentator.start ("Testing vector (hybrid 0-1)", __FUNCTION__);

	// std::ostream &report = LELA::commentator.report (LELA::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	// Note: In this version I haven't included anything, since
	// test-hybrid-vector.C seems to give hybrid vectors a good
	// workout already.

	LELA::commentator.stop (MSG_STATUS (pass));

	return pass;
}

/// Run a standard set of tests on the given vector-type
template <class Vector, class Ring>
bool testVector (const Ring &R)
	{ return testVectorSpec<Vector> (R, typename LELA::VectorTraits<Ring, Vector>::RepresentationType ()); }

#endif // __LELA_TESTS_TEST_VECTOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
