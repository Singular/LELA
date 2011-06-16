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

void randomVectorStartingAt (GaussJordan<Ring>::SparseMatrix::Row &v, size_t col, size_t coldim, MersenneTwister &MT)
{
	size_t t;
	GaussJordan<Ring>::SparseMatrix::Row::index_type idx;

	GaussJordan<Ring>::SparseMatrix::Row::word_type w;
	GaussJordan<Ring>::SparseMatrix::Row::word_type mask;

	GaussJordan<Ring>::SparseMatrix::Row::index_type colstart = col >> WordTraits<GaussJordan<Ring>::SparseMatrix::Row::word_type>::logof_size;
	GaussJordan<Ring>::SparseMatrix::Row::index_type colend = coldim >> WordTraits<GaussJordan<Ring>::SparseMatrix::Row::word_type>::logof_size;

	v.clear ();

	for (idx = colstart; idx <= colend; ++idx) {
		w = 0ULL;

		if (static_cast<size_t> (idx) << WordTraits<GaussJordan<Ring>::SparseMatrix::Row::word_type>::logof_size <= col) {
			t = col & WordTraits<GaussJordan<Ring>::SparseMatrix::Row::word_type>::pos_mask;
			mask = Endianness::e_j (t);
			w |= mask;  // Force leading entry to be one
		} else {
			t = 0;
			mask = Endianness::e_0;
		}

		for (; mask != 0 && (static_cast<size_t> (idx) << WordTraits<GaussJordan<Ring>::SparseMatrix::Row::word_type>::logof_size) + t < coldim; mask = Endianness::shift_right (mask, 1), ++t)
			if (MT.randomDoubleRange (0.0, 1.0) < nonzero_density)
				w |= mask;

		if (w != 0) {
			v.push_back (GaussJordan<Ring>::SparseMatrix::Row::value_type (idx, w));
		}
	}
}

void createRandomF4Matrix (GaussJordan<Ring>::SparseMatrix &A)
{
	MersenneTwister MT (time (NULL));
	GaussJordan<Ring>::SparseMatrix::RowIterator i_A;

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

bool smallTest () {
	bool pass = true;

	commentator.start ("Testing Faugère-Lachartre implementation", __FUNCTION__);

	size_t m = 64;
	size_t n = 96;

	GaussJordan<Ring>::SparseMatrix A (m, n), C (m, n);
	GaussJordan<Ring>::DenseMatrix L (m, m);
	GaussJordan<Ring>::Permutation P;

	createRandomF4Matrix (A);

	Ring F (2);
	Context<Ring> ctx (F);
	FaugereLachartre<Ring> Solver (ctx);
	GaussJordan<Ring> GJ (ctx);

	size_t rank;
	Ring::Element det;

	BLAS3::copy (ctx, A, C);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << std::endl;
	BLAS3::write (ctx, report, A);

	Solver.RowEchelonForm (A, A, rank, det);

	report << "Output matrix:" << std::endl;
	BLAS3::write (ctx, report, A);

	report << "Computed rank: " << rank << std::endl
	       << "Computed determinant: ";
	F.write (report, det) << std::endl;

	size_t rank1;
	Ring::Element det1;

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

// Get and print memory-usage of a sparse matrix

size_t ComputeMemoryUsage (const GaussJordan<Ring>::SparseMatrix &A) 
{
	GaussJordan<Ring>::SparseMatrix::ConstRowIterator i;

	size_t total_len = 0;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i)
		total_len += i->size ();

	return total_len;
}

// Check hybrid-vector to make sure it is valid

void CheckHybridMatrix (const GaussJordan<Ring>::SparseMatrix &A)
{
	GaussJordan<Ring>::SparseMatrix::ConstRowIterator i;

	commentator.start ("Checking that rows of matrix are valid", __FUNCTION__);

	for (i = A.rowBegin (); i != A.rowEnd (); ++i) {
		if (!VectorWrapper::isValid<GF2> (*i)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Error at row " << i - A.rowBegin () << std::endl;
			commentator.stop ("FAILED");
		}
	}

	commentator.stop ("passed");
}

// Real test using data from a PNG-file

void fileTest (char *filename, char *output) {
	Ring F (2);
	Context<Ring> ctx (F);

	GaussJordan<Ring>::SparseMatrix A;

	commentator.start ("Testing Faugère-Lachartre with PNG-file", __FUNCTION__);

	std::ifstream ifile (filename);

	if (!ifile.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open file " << filename << std::endl;
		return;
	}

	BLAS3::read (ctx, ifile, A, FORMAT_PNG);

	size_t total = ComputeMemoryUsage (A);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Total length of row-vectors: " << total << std::endl;

	CheckHybridMatrix (A);

	FaugereLachartre<Ring> Solver (ctx);

	size_t rank;
	Ring::Element det;

	Solver.RowEchelonForm (A, A, rank, det);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Computed rank: " << rank << std::endl << "Computed determinant: ";
	F.write (commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION), det) << std::endl;

	CheckHybridMatrix (A);

	std::ofstream f (output);
	BLAS3::write (ctx, f, A, FORMAT_PNG);

	commentator.stop ("done");
}

int main (int argc, char **argv)
{
	bool pass = true;

	static bool runFileTest = false;
	static char *inputFile = NULL;
	static char *outputFile = NULL;

	static Argument args[] = {
		{ 'f', "-f", "Run file-test", TYPE_NONE, &runFileTest },
		{ 'i', "-i", "Input-file for file-test", TYPE_STRING, &inputFile },
		{ 'o', "-o", "Output-file for file-test", TYPE_STRING, &outputFile },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (0);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	commentator.start ("Faugère-Lachartre test suite", "FaugereLachartre");

	if (runFileTest)
		fileTest (inputFile, outputFile);
	else
		pass = smallTest ();

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}
