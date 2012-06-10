/* tests/test-hybrid-vector.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test for hybrid sparse-dense vector-format
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include <iostream>

#include "test-common.h"

#include "lela/ring/gf2.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"
#include "lela/blas/level3.h"
#include "lela/vector/sparse-subvector-hybrid.h"

using namespace LELA;

typedef GF2 Field;

typedef HybridVector<DefaultEndianness<uint8>, uint32, uint8> TestVector;

bool testAdd ()
{
	commentator.start ("Testing BLAS1::axpy", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Field F (2);
	Context<Field> ctx (F);

	TestVector u, v, w, t;

	// Test 1: Indices in same block
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (TestVector::value_type (0, 1));
	v.push_back (TestVector::value_type (0, 1 | (1 << 1)));
	t.push_back (TestVector::value_type (0, 1 << 1));

	BLAS1::copy (ctx, u, w);
	BLAS1::axpy (ctx, F.one (), v, w);

	report << "Test 1:" << std::endl;
	report << "    u = ";
	BLAS1::write (ctx, report, u) << std::endl;
	report << "    v = ";
	BLAS1::write (ctx, report, v) << std::endl;
	report << "u + v = ";
	BLAS1::write (ctx, report, w) << std::endl;
	report << "    t = ";
	BLAS1::write (ctx, report, t) << std::endl;

	if (!BLAS1::equal (ctx, w, t)) {
		error << "ERROR (Test 1): VectorDomain reports u + v != t" << std::endl;
		pass = false;
	}

	// Test 2: Indices in different blocks, u before v
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (TestVector::value_type (0, 1));
	v.push_back (TestVector::value_type (1, 1));
	t.push_back (TestVector::value_type (0, 1));
	t.push_back (TestVector::value_type (1, 1));

	BLAS1::copy (ctx, u, w);
	BLAS1::axpy (ctx, F.one (), v, w);

	report << "Test 2:" << std::endl;
	report << "    u = ";
	BLAS1::write (ctx, report, u) << std::endl;
	report << "    v = ";
	BLAS1::write (ctx, report, v) << std::endl;
	report << "u + v = ";
	BLAS1::write (ctx, report, w) << std::endl;
	report << "    t = ";
	BLAS1::write (ctx, report, t) << std::endl;

	if (!BLAS1::equal (ctx, w, t)) {
		error << "ERROR (Test 2): VectorDomain reports u + v != t" << std::endl;
		pass = false;
	}

	// Test 3: Indices in different blocks, v before u
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (TestVector::value_type (1, 1));
	v.push_back (TestVector::value_type (0, 1));
	t.push_back (TestVector::value_type (0, 1));
	t.push_back (TestVector::value_type (1, 1));

	BLAS1::copy (ctx, u, w);
	BLAS1::axpy (ctx, F.one (), v, w);

	report << "Test 3:" << std::endl;
	report << "    u = ";
	BLAS1::write (ctx, report, u) << std::endl;
	report << "    v = ";
	BLAS1::write (ctx, report, v) << std::endl;
	report << "u + v = ";
	BLAS1::write (ctx, report, w) << std::endl;
	report << "    t = ";
	BLAS1::write (ctx, report, t) << std::endl;

	if (!BLAS1::equal (ctx, w, t)) {
		error << "ERROR (Test 3): VectorDomain reports u + v != t" << std::endl;
		pass = false;
	}

	// Test 4: u = 0
	u.clear ();
	v.clear ();
	t.clear ();

	v.push_back (TestVector::value_type (0, 1));
	t.push_back (TestVector::value_type (0, 1));

	BLAS1::copy (ctx, u, w);
	BLAS1::axpy (ctx, F.one (), v, w);

	report << "Test 4:" << std::endl;
	report << "    u = ";
	BLAS1::write (ctx, report, u) << std::endl;
	report << "    v = ";
	BLAS1::write (ctx, report, v) << std::endl;
	report << "u + v = ";
	BLAS1::write (ctx, report, w) << std::endl;
	report << "    t = ";
	BLAS1::write (ctx, report, t) << std::endl;

	if (!BLAS1::equal (ctx, w, t)) {
		error << "ERROR (Test 4): VectorDomain reports u + v != t" << std::endl;
		pass = false;
	}

	// Test 5: v = 0
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (TestVector::value_type (0, 1));
	t.push_back (TestVector::value_type (0, 1));

	BLAS1::copy (ctx, u, w);
	BLAS1::axpy (ctx, F.one (), v, w);

	report << "Test 5:" << std::endl;
	report << "    u = ";
	BLAS1::write (ctx, report, u) << std::endl;
	report << "    v = ";
	BLAS1::write (ctx, report, v) << std::endl;
	report << "u + v = ";
	BLAS1::write (ctx, report, w) << std::endl;
	report << "    t = ";
	BLAS1::write (ctx, report, t) << std::endl;

	if (!BLAS1::equal (ctx, w, t)) {
		error << "ERROR (Test 5): VectorDomain reports u + v != t" << std::endl;
		pass = false;
	}

	// Test 6: u, v have lengths longer than 1
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (TestVector::value_type (0, 2ULL));
	u.push_back (TestVector::value_type (1, 1ULL));

	v.push_back (TestVector::value_type (0, 4ULL));
	v.push_back (TestVector::value_type (1, 1ULL));

	t.push_back (TestVector::value_type (0, 2 | 4));

	BLAS1::copy (ctx, u, w);
	BLAS1::axpy (ctx, F.one (), v, w);

	report << "Test 6:" << std::endl;
	report << "    u = ";
	BLAS1::write (ctx, report, u) << std::endl;
	report << "    v = ";
	BLAS1::write (ctx, report, v) << std::endl;
	report << "u + v = ";
	BLAS1::write (ctx, report, w) << std::endl;
	report << "    t = ";
	BLAS1::write (ctx, report, t) << std::endl;

	if (!BLAS1::equal (ctx, w, t)) {
		error << "ERROR (Test 6): VectorDomain reports u + v != t" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

bool testFirstNonzeroEntry ()
{
	commentator.start ("Testing BLAS1::head", __FUNCTION__);

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Field F (2);
	Context<Field> ctx (F);

	bool a;

	TestVector v;

	int idx;

	v.push_back (TestVector::value_type (2, TestVector::Endianness::e_j (0)));

	idx = BLAS1::head (ctx, a, v);

	if (idx != 2 * WordTraits<TestVector::word_type>::bits) {
		error << "ERROR (Test 1): BLAS1::head (ctx, a, v) = " << idx
		      << " (should be " << 2 * 8 * sizeof (TestVector::word_type) << ")" << std::endl;
		pass = false;
	}

	v.clear ();

	v.push_back (TestVector::value_type (2, TestVector::Endianness::e_j (WordTraits<TestVector::word_type>::bits / 2)));

	idx = BLAS1::head (ctx, a, v);

	if (idx != 5 * WordTraits<TestVector::word_type>::bits / 2) {
		error << "ERROR (Test 2): BLAS1::head (ctx, a, v) = " << idx
		      << " (should be " << 2 * 8 * sizeof (TestVector::word_type) + 16 << ")" << std::endl;
		pass = false;
	}

	v.clear ();

	v.push_back (TestVector::value_type (3, TestVector::Endianness::e_j (WordTraits<TestVector::word_type>::bits - 1)));

	idx = BLAS1::head (ctx, a, v);

	if (idx != 3 * WordTraits<TestVector::word_type>::bits + (WordTraits<TestVector::word_type>::bits - 1)) {
		error << "ERROR (Test 3): BLAS1::head (ctx, a, v) = " << idx
		      << " (should be " << 3 * 8 * sizeof (TestVector::word_type) + (8 * sizeof (TestVector::word_type) - 1) << ")" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

TestVector::word_type connect (TestVector::word_type word1, TestVector::word_type word2, int shift) 
{
	if (shift == 0)
		return word1;
	else
		return TestVector::Endianness::shift_left (word1, shift) | TestVector::Endianness::shift_right (word2, WordTraits<TestVector::word_type>::bits - shift);
}

bool testSparseSubvectorHybrid ()
{
	commentator.start ("Testing SparseSubvector<Hybrid>", __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	TestVector v;

	//	TestVector::word_type pattern[2] = { 0xf0f0f0f0f0f0f0f0ULL, 0xaaaaaaaaaaaaaaaaULL };
	TestVector::word_type pattern[2] = { 0xf0ULL, 0xaaULL };
	TestVector::word_type check;

	uint16 offset = 5;
	uint16 len = WordTraits<TestVector::word_type>::bits + 4;

	report << "Test 1" << std::endl;

	v.push_back (TestVector::value_type (0, pattern[0]));
	v.push_back (TestVector::value_type (1, pattern[1]));

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01> v1 (v, offset, offset + len);

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01>::const_iterator i_v1 = v1.begin ();

	if (i_v1->first != 0) {
		error << "ERROR: First index should be 0, is " << std::dec << i_v1->first << std::endl;
		pass = false;
	}

	check = connect (pattern[0], pattern[1], offset);

	if (i_v1->second != check) {
		error << "ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v1->second) << std::endl;
		pass = false;
	}

	++i_v1;

	if (i_v1->first != 1) {
		error << "ERROR: Second index should be 1, is " << std::dec << i_v1->first << std::endl;
		pass = false;
	}

	check = connect (pattern[1], 0ULL, offset) & TestVector::Endianness::mask_left (len % WordTraits<TestVector::word_type>::bits);

	if (i_v1->second != check) {
		error << "ERROR: Second word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v1->second) << std::endl;
		pass = false;
	}

	++i_v1;

	if (i_v1 != v1.end ()) {
		error << "ERROR: Did not hit end when expected." << std::endl;
		pass = false;
	}

	v.clear ();

	report << "Test 2" << std::endl;

	v.push_back (TestVector::value_type (1, pattern[0]));

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01> v2 (v, offset, offset + len);

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01>::const_iterator i_v2 = v2.begin ();

	if (i_v2->first != 0) {
		error << "ERROR: First index should be 0, is " << std::dec << i_v2->first << std::endl;
		pass = false;
	}

	check = connect (0ULL, pattern[0], offset);

	if (i_v2->second != check) {
		error << "ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v2->second) << std::endl;
		pass = false;
	}

	++i_v2;

	if (i_v2 == v2.end ()) {
		error << "ERROR: Iterator ended prematurely." << std::endl;
		pass = false;
	}

	if (i_v2->first != 1) {
		error << "ERROR: Second index should be 1, is " << std::dec << i_v2->first << std::endl;
		pass = false;
	}

	check = connect (pattern[0], 0ULL, offset) & TestVector::Endianness::mask_left (len % WordTraits<TestVector::word_type>::bits);

	if (i_v2->second != check) {
		error << "ERROR: Second word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v2->second) << std::endl;
		pass = false;
	}

	++i_v2;

	if (i_v2 != v2.end ()) {
		error << "ERROR: Did not hit end when expected." << std::endl;
		pass = false;
	}

	report << "Test 3" << std::endl;

	v.clear ();

	v.push_back (TestVector::value_type (0, pattern[0]));
	v.push_back (TestVector::value_type (2, pattern[1]));

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01> v3 (v, offset, WordTraits<TestVector::word_type>::bits + offset + len);

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01>::const_iterator i_v3 = v3.begin ();

	if (i_v3->first != 0) {
		error << "ERROR: First index should be 0, is " << std::dec << i_v3->first << std::endl;
		pass = false;
	}

	check = connect (pattern[0], 0ULL, offset);

	if (i_v3->second != check) {
		error << "ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v3->second) << std::endl;
		pass = false;
	}

	++i_v3;

	if (i_v3 == v3.end ()) {
		error << "ERROR: Iterator ended prematurely." << std::endl;
		pass = false;
	}

	if (i_v3->first != 1) {
		error << "ERROR: Second index should be 1, is " << std::dec << i_v3->first << std::endl;
		pass = false;
	}

	check = connect (0ULL, pattern[1], offset);

	if (i_v3->second != check) {
		error << "ERROR: Second word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v3->second) << std::endl;
		pass = false;
	}

	++i_v3;

	if (i_v3 == v3.end ()) {
		error << "ERROR: Iterator ended prematurely." << std::endl;
		pass = false;
	}

	if (i_v3->first != 2) {
		error << "ERROR: Third index should be 2, is " << std::dec << i_v3->first << std::endl;
		pass = false;
	}

	check = connect (pattern[1], 0ULL, offset) & TestVector::Endianness::mask_left (len % WordTraits<TestVector::word_type>::bits);

	if (i_v3->second != check) {
		error << "ERROR: Third word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v3->second) << std::endl;
		pass = false;
	}

	++i_v3;

	if (i_v3 != v3.end ()) {
		error << "ERROR: Did not hit end when expected." << std::endl;
		pass = false;
	}

	report << "Test 4" << std::endl;

	v.clear ();

	v.push_back (TestVector::value_type (0, pattern[0]));

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01> v4 (v, offset, WordTraits<TestVector::word_type>::bits + offset + len);

	SparseSubvector<const TestVector, VectorRepresentationTypes::Hybrid01>::const_iterator i_v4 = v4.begin ();

	if (i_v4->first != 0) {
		error << "ERROR: First index should be 0, is " << std::dec << i_v4->first << std::endl;
		pass = false;
	}

	check = connect (pattern[0], 0ULL, offset);

	if (i_v4->second != check) {
		error << "ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v4->second) << std::endl;
		pass = false;
	}

	++i_v4;

	if (i_v4 != v4.end ()) {
		error << "ERROR: Did not hit end when expected." << std::endl;
		pass = false;
	}
	
	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("Hybrid 0-1 vector test-suite", "HybridVector");

	pass = testAdd () && pass;
	pass = testFirstNonzeroEntry () && pass;
	pass = testSparseSubvectorHybrid () && pass;

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
