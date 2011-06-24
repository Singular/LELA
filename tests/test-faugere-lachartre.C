/* tests/test-faugere-lachartre.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Test for implementation of algorithm of Faugère and Lachartre
 */

#include <iostream>
#include <ctime>

#include "test-common.h"

#include "linbox/util/commentator.h"
#include "linbox/ring/gf2.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/algorithms/faugere-lachartre.h"

using namespace LinBox;

typedef GF2 Ring;

static const double nonzero_density = 0.1;

typedef Adaptor<Ring>::Endianness Endianness;

// Generate a random sparse vector in which all nonzero entries occur after column col and the entry at col is guaranteed to be one

template <class Vector>
void randomVectorStartingAtSpec (Vector &v, size_t col, size_t coldim, MersenneTwister &MT, VectorCategories::HybridZeroOneVectorTag)
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

template <class Vector>
void randomVectorStartingAt (Vector &v, size_t col, size_t coldim, MersenneTwister &MT)
	{ randomVectorStartingAtSpec (v, col, coldim, MT, typename VectorTraits<Ring, Vector>::VectorCategory ()); }

template <class Matrix>
void createRandomF4Matrix (Matrix &A)
{
	MersenneTwister MT (time (NULL));
	typename Matrix::RowIterator i_A;

	size_t col = 0;

	for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A) {
		LinBox::uint32 v = MT.randomIntRange (0, 7);

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

		randomVectorStartingAt (*i_A, col, A.coldim (), MT);
	}
}

// Small version of the test, for debugging

template <class Ring>
bool testFaugereLachartre (const Ring &R, size_t m, size_t n)
{
	bool pass = true;

	commentator.start ("Testing Faugère-Lachartre implementation", __FUNCTION__);

	typename DefaultSparseMatrix<Ring>::Type A (m, n), C (m, n);
	DenseMatrix<typename Ring::Element> L (m, m);
	typename GaussJordan<Ring>::Permutation P;

	createRandomF4Matrix (A);

	Context<Ring> ctx (R);
	FaugereLachartre<Ring> Solver (ctx);
	GaussJordan<Ring> GJ (ctx);

	size_t rank;
	typename Ring::Element det;

	BLAS3::copy (ctx, A, C);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << std::endl;
	BLAS3::write (ctx, report, A);

	Solver.RowEchelonForm (A, A, rank, det);

	report << "Output matrix:" << std::endl;
	BLAS3::write (ctx, report, A);

	report << "Computed rank: " << rank << std::endl
	       << "Computed determinant: ";
	R.write (report, det) << std::endl;

	size_t rank1;
	typename Ring::Element det1;

	GJ.StandardRowEchelonForm (C, L, P, rank1, det1, true, false);

	report << "True reduced row-echelon form:" << std::endl;
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
	static long m = 128;
	static long n = 96;

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

	GF2 gf2;

	pass = testFaugereLachartre (gf2, m, n);

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}
