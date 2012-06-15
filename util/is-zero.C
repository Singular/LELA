/* util/is-zero.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Utility to check whether a matrix is zero
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/util/commentator.h"
#include "lela/blas/context.h"
#include "lela/ring/gf2.h"
#include "lela/ring/modular.h"
#include "lela/blas/level1.h"
#include "lela/blas/level3.h"
#include "lela/matrix/sparse.h"

#include "support.h"

using namespace LELA;

template <class Ring>
struct DefaultSparseMatrix
{
	typedef SparseMatrix<typename Ring::Element> Type;
};

template <>
struct DefaultSparseMatrix<GF2>
{
	typedef SparseMatrix<bool, Vector<GF2>::Hybrid> Type;
};

template <class Ring>
int check_is_zero (const Ring &R, const char *input, FileFormatTag input_format)
{
	bool res;

	Context<Ring> ctx (R);

	commentator.start ("Checking whether matrix is zero", __FUNCTION__);

	typename DefaultSparseMatrix<Ring>::Type A;

	commentator.start ("Reading input-matrix");

	std::ifstream ifile (input);

	if (!ifile.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open first input-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	try {
		BLAS3::read (ctx, ifile, A, input_format);
	}
	catch (LELAError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}
	catch (UnrecognisedFormat) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Unable to determine format of input-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	commentator.stop (MSG_DONE);

	res = BLAS3::is_zero (ctx, A);

	commentator.stop (res ? "zero" : "not zero");

	return res ? 0 : 1;
}

int main (int argc, char **argv)
{
	static const char *ringString = "guess";
	static integer p = 65521;
	static bool floatingPoint = false;
	static const char *inputFileFormat = "guess";
	static char *input = NULL;

	static Argument args[] = {
		{ 'k', "-k", "Ring over which to compute ('guess', 'gf2', 'modular')", TYPE_STRING, &ringString },
		{ 'p', "-p", "Modulus of ring, when ring is 'modular'", TYPE_INT, &p },
		{ 'f', "-f", "Compute using floating point, when ring is 'modular'", TYPE_NONE, &floatingPoint },
		{ 'i', "-i", "File format of input-matrix ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png')", TYPE_STRING, &inputFileFormat },
		{ '\0' }
	};

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);

	parseArguments (argc, argv, args, "<first matrix filename> <second matrix filename>", 1, &input);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	if (input == NULL) {
		printHelpMessage (argv[0], args, "<matrix filename>", true);
		return -1;
	}

	FileFormatTag input_format = get_format_tag (inputFileFormat);

	if (input_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (input_format == FORMAT_DETECT)
		input_format = guess_format (input);

	RingType ring_type = get_ring_type (ringString);

	if (ring_type == RING_GUESS) {
#ifdef __LELA_HAVE_LIBPNG
		if (input_format == FORMAT_PNG)
			ring_type = RING_GF2;
		else
			ring_type = RING_MODULAR;
#else // !__LELA_HAVE_LIBPNG
		ring_type = RING_MODULAR;
#endif // __LELA_HAVE_LIBPNG
	}

	if (ring_type == RING_GF2)
		return check_is_zero (GF2 (), input, input_format);
	else if (ring_type == RING_MODULAR) {
		if (floatingPoint) {
			if (ModularTraits<float>::valid_modulus (p))
				return check_is_zero (Modular<float> (p), input, input_format);
			else if (ModularTraits<double>::valid_modulus (p))
				return check_is_zero (Modular<double> (p), input, input_format);
			else
				return check_is_zero (Modular<integer> (p), input, input_format);
		} else {
			if (ModularTraits<uint8>::valid_modulus (p))
				return check_is_zero (Modular<uint8> (p), input, input_format);
			else if (ModularTraits<uint16>::valid_modulus (p))
				return check_is_zero (Modular<uint16> (p), input, input_format);
			else if (ModularTraits<uint32>::valid_modulus (p))
				return check_is_zero (Modular<uint32> (p), input, input_format);
			else
				return check_is_zero (Modular<integer> (p), input, input_format);
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
