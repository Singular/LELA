/* -*- mode: c++;-*-
 * test-solver.C
 * Test code to reduce F4-matrix to row-echelon form
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include <iostream>
#include <ctime>

#include "linbox/util/commentator.h"
#include "linbox/field/gf2.h"
#include "linbox/randiter/mersenne-twister.h"

#include "solver.h"

using namespace LinBox;
using namespace F4;

typedef GF2 Field;

static const double nonzero_density = 0.1;

namespace F4Tests {

typedef Adaptor<Field>::Endianness Endianness;
typedef GaussJordan<Field>::SparseMatrix SparseMatrix;
typedef GaussJordan<Field>::SparseMatrix::Row SparseVector;

// Generate a random sparse vector in which all nonzero entries occur after column col and the entry at col is guaranteed to be one

void randomVectorStartingAt (SparseVector &v, size_t col, size_t coldim, MersenneTwister &MT)
{
	size_t t;
	SparseVector::index_type idx;

	SparseVector::word_type w;
	SparseVector::word_type mask;

	SparseVector::index_type colstart = col >> WordTraits<SparseVector::word_type>::logof_size;
	SparseVector::index_type colend = coldim >> WordTraits<SparseVector::word_type>::logof_size;

	v.clear ();

	for (idx = colstart; idx <= colend; ++idx) {
		w = 0ULL;

		if (static_cast<size_t> (idx) << WordTraits<SparseVector::word_type>::logof_size <= col) {
			t = col & WordTraits<SparseVector::word_type>::pos_mask;
			mask = Endianness::e_j (t);
			w |= mask;  // Force leading entry to be one
		} else {
			t = 0;
			mask = Endianness::e_0;
		}

		for (; mask != 0 && (static_cast<size_t> (idx) << WordTraits<SparseVector::word_type>::logof_size) + t < coldim; mask = Endianness::shift_right (mask, 1), ++t)
			if (MT.randomDoubleRange (0.0, 1.0) < nonzero_density)
				w |= mask;

		if (w != 0) {
			v.push_back (SparseVector::value_type (idx, w));
		}
	}
}

void createRandomF4Matrix (SparseMatrix &A)
{
	MersenneTwister MT (time (NULL));
	SparseMatrix::RowIterator i_A;

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

void smallTest () {
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	size_t m = 64;
	size_t n = 96;

	SparseMatrix A (m, n), C (m, n);
	GaussJordan<Field>::DenseMatrix L (m, m);
	GaussJordan<Field>::Permutation P;

	createRandomF4Matrix (A);

	Field F (2);
	MatrixDomain<Field> MD (F);
	F4Solver<Field> Solver (F);
	GaussJordan<Field> GJ (F);

	size_t rank;
	Field::Element det;

	MD.copy (C, A);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << std::endl;
	MD.write (report, A);

	Solver.RowEchelonForm (A, A, rank, det);

	report << "Output matrix:" << std::endl;
	MD.write (report, A);

	report << "Computed rank: " << rank << std::endl
	       << "Computed determinant: ";
	F.write (report, det) << std::endl;

	size_t rank1;
	Field::Element det1;

	GJ.StandardRowEchelonForm (C, L, P, rank1, det1, true, false);

	report << "True reduced row-echelon form:" << std::endl;
	MD.write (report, C);

	report << "True rank: " << rank1 << std::endl;

	if (!MD.areEqual (A, C))
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << "ERROR: Output-matrices are not equal!" << std::endl;
	if (rank != rank1)
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << "ERROR: Computed ranks are not equal!" << std::endl;
}

// Get and print memory-usage of a sparse matrix

size_t ComputeMemoryUsage (const SparseMatrix &A) 
{
	SparseMatrix::ConstRowIterator i;

	size_t total_len = 0;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i)
		total_len += i->size ();

	return total_len;
}

// Check hybrid-vector to make sure it is valid

void CheckHybridMatrix (const SparseMatrix &A)
{
	SparseMatrix::ConstRowIterator i;

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
	Field F (2);
	MatrixDomain<Field> MD (F);

	SparseMatrix A;

	commentator.start ("Testing F4-solver with PNG-file", __FUNCTION__);

	std::ifstream ifile (filename);

	if (!ifile.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open file " << filename << std::endl;
		return;
	}

	MD.read (ifile, A, FORMAT_PNG);

	size_t total = ComputeMemoryUsage (A);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Total length of row-vectors: " << total << std::endl;

	CheckHybridMatrix (A);

	F4Solver<Field> Solver (F);

	size_t rank;
	Field::Element det;

	Solver.RowEchelonForm (A, A, rank, det);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Computed rank: " << rank << std::endl << "Computed determinant: ";
	F.write (commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION), det) << std::endl;

	CheckHybridMatrix (A);

	std::ofstream f (output);
	MD.write (f, A, FORMAT_PNG);

	commentator.stop ("done");
}

} // namespace F4Tests

int main (int argc, char **argv)
{
	commentator.setBriefReportParameters (Commentator::OUTPUT_PIPE, false, false, false);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (0);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);
	commentator.setBriefReportStream (commentator.cnull);
	commentator.setReportStream (std::cout);

	if (argc == 2 && !strcmp (argv[1], "-s"))
		F4Tests::smallTest ();
	else if (argc == 4 && !strcmp (argv[1], "-f"))
		F4Tests::fileTest (argv[2], argv[3]);
	else {
		std::cout << "Usage: test-solver -s|{-f <png-input> <matrix-output>}" << std::endl;
		return 0;
	}
}
