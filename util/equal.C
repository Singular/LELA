/* util/equal.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Utility to check whether two matrices are equal
 */

#include "linbox/util/commentator.h"
#include "linbox/blas/context.h"
#include "linbox/ring/gf2.h"
#include "linbox/ring/modular.h"
#include "linbox/blas/level1.h"
#include "linbox/blas/level3.h"

#include "support.h"

using namespace LinBox;

template <class Ring>
int check_equal (const Ring &R, const char *input1, FileFormatTag input1_format, const char *input2, FileFormatTag input2_format)
{
	bool res;

	Context<Ring> ctx (R);

	commentator.start ("Checking equality", __FUNCTION__);

	DenseMatrix<typename Ring::Element> A, B;

	commentator.start ("Reading input-matrix 1");

	std::ifstream ifile1 (input1);

	if (!ifile1.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open first input-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	try {
		BLAS3::read (ctx, ifile1, A, input1_format);
	}
	catch (LinboxError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	commentator.stop (MSG_DONE);

	commentator.start ("Reading input-matrix 2");

	std::ifstream ifile2 (input2);

	if (!ifile2.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open second input-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	try {
		BLAS3::read (ctx, ifile2, B, input2_format);
	}
	catch (LinboxError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	commentator.stop (MSG_DONE);

	res = BLAS3::equal (ctx, A, B);

	commentator.stop (res ? "equal" : "not equal");

	return res ? 0 : 1;
}

int main (int argc, char **argv)
{
	static const char *ringString = "modular";
	static integer p = 101;
	static bool floatingPoint = false;
	static const char *input1FileFormat = "guess";
	static const char *input2FileFormat = "guess";
	static char *input1 = NULL;
	static char *input2 = NULL;

	static Argument args[] = {
		{ 'k', "-k", "Ring over which to compute ('gf2', 'modular')", TYPE_STRING, &ringString },
		{ 'p', "-p", "Modulus of ring, when ring is 'modular'", TYPE_INT, &p },
		{ 'f', "-f", "Compute using floating point, when ring is 'modular'", TYPE_NONE, &floatingPoint },
		{ '1', "-1", "File format of first input-matrix ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png')", TYPE_STRING, &input1FileFormat },
		{ '2', "-2", "File format of second input-matrix ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png')", TYPE_STRING, &input2FileFormat },
		{ '\0' }
	};

	parseArguments (argc, argv, args, 2, &input1, &input2);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	FileFormatTag input1_format = get_format_tag (input1FileFormat);
	FileFormatTag input2_format = get_format_tag (input2FileFormat);

	if (input1_format == FORMAT_UNKNOWN || input1_format == FORMAT_PRETTY) {
		std::cerr << "Invalid first input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (input2_format == FORMAT_UNKNOWN || input2_format == FORMAT_PRETTY) {
		std::cerr << "Invalid second input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	RingType ring_type = get_ring_type (ringString);

	if (ring_type == RING_GF2)
		return check_equal (GF2 (), input1, input1_format, input2, input2_format);
	else if (ring_type == RING_MODULAR) {
		if (floatingPoint) {
			if (p <= Modular<float>::getMaxModulus ())
				return check_equal (Modular<float> (p), input1, input1_format, input2, input2_format);
			else if (p <= Modular<double>::getMaxModulus ())
				return check_equal (Modular<double> (p), input1, input1_format, input2, input2_format);
			else
				return check_equal (Modular<integer> (p), input1, input1_format, input2, input2_format);
		} else {
			if (p <= Modular<uint8>::getMaxModulus ())
				return check_equal (Modular<uint8> (p), input1, input1_format, input2, input2_format);
			else if (p <= Modular<uint16>::getMaxModulus ())
				return check_equal (Modular<uint16> (p), input1, input1_format, input2, input2_format);
			else if (p <= Modular<uint32>::getMaxModulus ())
				return check_equal (Modular<uint32> (p), input1, input1_format, input2, input2_format);
			else
				return check_equal (Modular<integer> (p), input1, input1_format, input2, input2_format);
		}
	}
	else if (ring_type == RING_UNKNOWN) {
		std::cerr << "Invalid ring (use 'gf2' or 'modular')" << std::endl;
		return -1;
	}
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
