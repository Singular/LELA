/* util/diff.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Utility to subtract two matrices and output the result
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
#include "lela/matrix/dense.h"

#include "support.h"

using namespace LELA;

template <class Ring>
int run_diff (const Ring &R, const char *input1, FileFormatTag input1_format, const char *input2, FileFormatTag input2_format,
	      const char *output, FileFormatTag output_format)
{
	Context<Ring> ctx (R);

	commentator.start ("Computing difference", __FUNCTION__);

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
	catch (LELAError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}
	catch (UnrecognisedFormat) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Unable to determine format of first input-file" << std::endl;
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
	catch (LELAError e) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}
	catch (UnrecognisedFormat) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Unable to determine format of second input-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	commentator.stop (MSG_DONE);

	BLAS3::axpy (ctx, ctx.F.minusOne (), B, A);

	std::ostringstream str;
	str << "Writing output-matrix (format " << format_names[output_format] << ")" << std::ends;
	commentator.start (str.str ().c_str ());

	std::ofstream ofile (output);

	if (!ofile.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open output-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	BLAS3::write (ctx, ofile, A, output_format);

	commentator.stop (MSG_DONE);

	commentator.stop (MSG_DONE);

	return 0;
}

int main (int argc, char **argv)
{
	static const char *ringString = "guess";
	static integer p = 65521;
	static bool floatingPoint = false;
	static const char *input1FileFormat = "guess";
	static const char *input2FileFormat = "guess";
	static const char *outputFileFormat = "guess";
	static char *input1 = NULL;
	static char *input2 = NULL;
	static char *output = NULL;

	static Argument args[] = {
		{ 'k', "-k", "Ring over which to compute ('guess', 'gf2', 'modular')", TYPE_STRING, &ringString },
		{ 'p', "-p", "Modulus of ring, when ring is 'modular'", TYPE_INT, &p },
		{ 'f', "-f", "Compute using floating point, when ring is 'modular'", TYPE_NONE, &floatingPoint },
		{ '1', "-1", "File format of first input-matrix ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png')", TYPE_STRING, &input1FileFormat },
		{ '2', "-2", "File format of second input-matrix ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png')", TYPE_STRING, &input2FileFormat },
		{ 'o', "-o", "Output file format ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png', 'pretty')", TYPE_STRING, &outputFileFormat },
		{ '\0' }
	};

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);

	parseArguments (argc, argv, args, "<first matrix filename> <second matrix filename> <output matrix filename>", 3, &input1, &input2, &output);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	if (input1 == NULL || input2 == NULL) {
		printHelpMessage (argv[0], args, "<first matrix filename> <second matrix filename> <output matrix filename>", true);
		return -1;
	}

	FileFormatTag input1_format = get_format_tag (input1FileFormat);
	FileFormatTag input2_format = get_format_tag (input2FileFormat);
	FileFormatTag output_format = get_format_tag (outputFileFormat);

	if (input1_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid first input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (input2_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid second input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (output_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid output-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (input1_format == FORMAT_DETECT)
		input1_format = guess_format (input1);

	if (input2_format == FORMAT_DETECT)
		input2_format = guess_format (input2);

	if (output_format == FORMAT_DETECT) {
		output_format = guess_format_tag (output);

		if (output_format == FORMAT_UNKNOWN) {
			std::cerr << "Could not guess output file-format" << std::endl;
			return -1;
		}
	}

	RingType ring_type = get_ring_type (ringString);

	if (ring_type == RING_GUESS) {
#ifdef __LELA_HAVE_LIBPNG
		if (input1_format == FORMAT_PNG || input2_format == FORMAT_PNG)
			ring_type = RING_GF2;
		else
			ring_type = RING_MODULAR;
#else // !__LELA_HAVE_LIBPNG
		ring_type = RING_MODULAR;
#endif // __LELA_HAVE_LIBPNG
	}

	if (ring_type == RING_GF2)
		return run_diff (GF2 (), input1, input1_format, input2, input2_format, output, output_format);
	else if (ring_type == RING_MODULAR) {
		if (floatingPoint) {
			if (ModularTraits<float>::valid_modulus (p))
				return run_diff (Modular<float> (p), input1, input1_format, input2, input2_format, output, output_format);
			else if (ModularTraits<double>::valid_modulus (p))
				return run_diff (Modular<double> (p), input1, input1_format, input2, input2_format, output, output_format);
			else
				return run_diff (Modular<integer> (p), input1, input1_format, input2, input2_format, output, output_format);
		} else {
			if (ModularTraits<uint8>::valid_modulus (p))
				return run_diff (Modular<uint8> (p), input1, input1_format, input2, input2_format, output, output_format);
			else if (ModularTraits<uint16>::valid_modulus (p))
				return run_diff (Modular<uint16> (p), input1, input1_format, input2, input2_format, output, output_format);
			else if (ModularTraits<uint32>::valid_modulus (p))
				return run_diff (Modular<uint32> (p), input1, input1_format, input2, input2_format, output, output_format);
			else
				return run_diff (Modular<integer> (p), input1, input1_format, input2, input2_format, output, output_format);
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
