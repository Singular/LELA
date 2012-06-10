/* tests/test-faugere-lachartre.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Test for implementation of algorithm of Faugère and Lachartre
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include <iostream>
#include <ctime>
#include <cmath>

#include "test-common.h"

#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/ring/gf2.h"
#include "lela/randiter/mersenne-twister.h"
#include "lela/algorithms/faugere-lachartre.h"
#include "lela/algorithms/elimination.h"

using namespace LELA;

typedef GF2 Ring;

static const double nonzero_density = 0.1;

typedef Vector<GF2>::Hybrid::Endianness Endianness;

// Generate a random sparse vector in which all nonzero entries occur after column col and the entry at col is guaranteed to be one

template <class Ring, class Vector>
void randomVectorStartingAtSpec (const Ring &R, Vector &v, size_t col, size_t coldim, MersenneTwister &MT, VectorRepresentationTypes::Sparse)
{
	double val;
	int skip;
	typename Vector::value_type::first_type idx = col;

	NonzeroRandIter<Ring> ri (R, typename Ring::RandIter (R));

	v.clear ();

	while (idx < coldim) {
		v.push_back (typename Vector::value_type (idx, typename Ring::Element ()));
		ri.random (v.back ().second);
		val = MT.randomDouble ();
		skip = (int) (ceil (log (val) / log (1 - nonzero_density)));
		idx += std::max (skip, 1);
	}
}

template <class Ring, class Vector>
void randomVectorStartingAtSpec (const Ring &R, Vector &v, size_t col, size_t coldim, MersenneTwister &MT, VectorRepresentationTypes::Dense01)
{
	size_t colend = (coldim + WordTraits<typename Vector::word_type>::bits - 1) & ~WordTraits<typename Vector::word_type>::pos_mask;
	size_t t = 0, idx;

	typename Vector::word_type w;
	typename Vector::word_type mask;

	for (idx = col & ~WordTraits<typename Vector::word_type>::pos_mask; idx < colend; idx += t) {
		w = 0ULL;

		if (idx <= col) {
			t = col & WordTraits<typename Vector::word_type>::pos_mask;
			mask = Endianness::e_j (t);
			w |= mask;  // Force leading entry to be one
		} else {
			t = 0;
			mask = Endianness::e_0;
		}

		for (; t < WordTraits<typename Vector::word_type>::bits && idx + t < coldim; mask = Endianness::shift_right (mask, 1), ++t)
			if (MT.randomDoubleRange (0.0, 1.0) < nonzero_density)
				w |= mask;

		*(v.word_begin () + (idx >> WordTraits<typename Vector::word_type>::logof_size)) = w;
	}
}

template <class Ring, class Vector>
void randomVectorStartingAtSpec (const Ring &R, Vector &v, size_t col, size_t coldim, MersenneTwister &MT, VectorRepresentationTypes::Sparse01)
{
	double val;
	int skip;
	typename Vector::value_type idx = col;

	v.clear ();

	while (idx < coldim) {
		v.push_back (idx);
		val = MT.randomDouble ();
		skip = (int) (ceil (log (val) / log (1 - nonzero_density)));
		idx += std::max (skip, 1);
	}
}

template <class Ring, class Vector>
void randomVectorStartingAtSpec (const Ring &R, Vector &v, size_t col, size_t coldim, MersenneTwister &MT, VectorRepresentationTypes::Hybrid01)
{
	size_t t;
	typename Vector::index_type idx;

	typename Vector::word_type w;
	typename Vector::word_type mask;

	typename Vector::index_type colstart = col >> WordTraits<typename Vector::word_type>::logof_size;
	typename Vector::index_type colend = coldim >> WordTraits<typename Vector::word_type>::logof_size;

	v.clear ();

	for (idx = colstart; idx <= colend; ++idx) {
		w = 0ULL;

		if (static_cast<size_t> (idx) << WordTraits<typename Vector::word_type>::logof_size <= col) {
			t = col & WordTraits<typename Vector::word_type>::pos_mask;
			mask = Endianness::e_j (t);
			w |= mask;  // Force leading entry to be one
		} else {
			t = 0;
			mask = Endianness::e_0;
		}

		for (; mask != 0 && (static_cast<size_t> (idx) << WordTraits<typename Vector::word_type>::logof_size) + t < coldim; mask = Endianness::shift_right (mask, 1), ++t)
			if (MT.randomDoubleRange (0.0, 1.0) < nonzero_density)
				w |= mask;

		if (w != 0) {
			v.push_back (typename Vector::value_type (idx, w));
		}
	}
}

template <class Ring, class Vector>
void randomVectorStartingAt (const Ring &R, Vector &v, size_t col, size_t coldim, MersenneTwister &MT)
	{ randomVectorStartingAtSpec (R, v, col, coldim, MT, typename VectorTraits<Ring, Vector>::RepresentationType ()); }

template <class Ring, class Matrix>
void createRandomF4Matrix (const Ring &R, Matrix &A)
{
	MersenneTwister MT;
	typename Matrix::RowIterator i_A;

	size_t col = 0;

	for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A) {
		LELA::uint32 v = MT.randomIntRange (0, 7);

		switch (v) {
		case 0:
			break;

		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
			++col;
			break;

		case 6:
			col = MT.randomIntRange (col, A.coldim () - (A.rowEnd () - i_A + 1));
			break;
		}

		if (col >= A.coldim ())
			break;

		randomVectorStartingAt (R, *i_A, col, A.coldim (), MT);
	}
}

// Small version of the test, for debugging

template <class Ring, class Matrix>
bool testFaugereLachartre (const Ring &R, const char *text, size_t m, size_t n, bool reduced)
{
	bool pass = true;

	std::ostringstream str;
	str << "Testing Faugère-Lachartre implementation over " << text << " (" << (reduced ? "reduced" : "non-reduced") << " variant)" << std::ends;

	commentator.start (str.str ().c_str (), __FUNCTION__);

	Matrix A (m, n), C (m, n);
	DenseMatrix<typename Ring::Element> L (m, m);
	typename GaussJordan<Ring>::Permutation P;

	createRandomF4Matrix (R, A);

	Context<Ring> ctx (R);
	FaugereLachartre<Ring> Solver (ctx);
	Elimination<Ring> elim (ctx);

	size_t rank;
	typename Ring::Element det;

	BLAS3::copy (ctx, A, C);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << std::endl;
	BLAS3::write (ctx, report, A);

	Solver.echelonize (A, A, rank, det, reduced);

	report << "Output matrix:" << std::endl;
	BLAS3::write (ctx, report, A);

	report << "Computed rank: " << rank << std::endl
	       << "Computed determinant: ";
	R.write (report, det) << std::endl;

	if (!reduced) {
		elim.echelonize_reduced (A, L, P, rank, det);

		report << "Output matrix after reduction:" << std::endl;
		BLAS3::write (ctx, report, A);
	}

	size_t rank1;
	typename Ring::Element det1;

	elim.echelonize_reduced (C, L, P, rank1, det1);

	report << "True row-echelon form:" << std::endl;
	BLAS3::write (ctx, report, C);

	report << "True rank: " << rank1 << std::endl;

	if (!BLAS3::equal (ctx, A, C)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << "ERROR: Output-matrices are not equal!" << std::endl;
		pass = false;
	}

	if (rank != rank1) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << "ERROR: Computed ranks are not equal!" << std::endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	static long m = 96;
	static long n = 128;

	bool pass = true;

	static Argument args[] = {
		{ 'm', "-m M", "Set row-dimension of matrix to M.", TYPE_INT, &m },
		{ 'n', "-n N", "Set column-dimension of matrix to N.", TYPE_INT, &n },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	commentator.start ("Faugère-Lachartre test suite", "FaugereLachartre");

	Modular<float> R (101);

	pass = testFaugereLachartre<Modular<float>, DefaultSparseMatrix<Modular<float> >::Type> (R, "GF(5)", m, n, false);
	pass = testFaugereLachartre<Modular<float>, DefaultSparseMatrix<Modular<float> >::Type> (R, "GF(5)", m, n, true) && pass;

	GF2 gf2;

	typedef SparseMatrix<bool, HybridVector<DefaultEndianness<uint64>, uint16, uint64> > GF2Matrix;
	// typedef DenseMatrix<bool> GF2Matrix;
	// typedef SparseMatrix<bool, Vector<GF2>::Sparse> GF2Matrix;

	pass = testFaugereLachartre<GF2, GF2Matrix> (gf2, "GF(2)", m, n, false) && pass;
	pass = testFaugereLachartre<GF2, GF2Matrix> (gf2, "GF(2)", m, n, true) && pass;

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
