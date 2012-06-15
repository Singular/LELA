/* util/submatrix.C
 * Copyright 2012 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Utility to extract a submatrix from a matrix
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
int get_submatrix (const Ring &R, const char *input, FileFormatTag input_format, size_t start_row, size_t start_col, size_t rows, size_t cols, const char *output, FileFormatTag output_format)
{
	Context<Ring> ctx (R);

	commentator.start ("Constructing submatrix", __FUNCTION__);

	typename DefaultSparseMatrix<Ring>::Type A;

	commentator.start ("Reading input-matrix");

	std::ifstream ifile1 (input);

	if (!ifile1.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open first input-file" << std::endl;
		commentator.stop ("error");
		commentator.stop ("error");
		return -1;
	}

	try {
		BLAS3::read (ctx, ifile1, A, input_format);
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

	typename DefaultSparseMatrix<Ring>::Type::ConstSubmatrixType B (A, start_row, start_col, rows, cols);

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

	BLAS3::write (ctx, ofile, B, output_format);

	commentator.stop (MSG_DONE);

	commentator.stop (MSG_DONE);

	return 0;
}

int main (int argc, char **argv)
{
	static const char *ringString = "guess";
	static integer p = 65521;
	static bool floatingPoint = false;
	static const char *inputFileFormat = "guess";
	static const char *outputFileFormat = "guess";
	static char *input = NULL;
	static char *start_row_str = NULL;
	static char *start_col_str = NULL;
	static char *rows_str = NULL;
	static char *cols_str = NULL;
	static char *output = NULL;

	static Argument args[] = {
		{ 'k', "-k", "Ring over which to compute ('guess', 'gf2', 'modular')", TYPE_STRING, &ringString },
		{ 'p', "-p", "Modulus of ring, when ring is 'modular'", TYPE_INT, &p },
		{ 'f', "-f", "Compute using floating point, when ring is 'modular'", TYPE_NONE, &floatingPoint },
		{ 'i', "-i", "File format of input-matrix ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png')", TYPE_STRING, &inputFileFormat },
		{ 'o', "-o", "Output file format ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png', 'pretty')", TYPE_STRING, &outputFileFormat },
		{ '\0' }
	};

	size_t start_row, start_col, rows, cols;

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);

	parseArguments (argc, argv, args, "<input matrix filename> <start row> <start col> <rows> <cols> <output matrix filename>", 6,
			&input, &start_row_str, &start_col_str, &rows_str, &cols_str, &output);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	if (input == NULL || start_row_str == NULL || start_col_str == NULL || rows_str == NULL || cols_str == NULL || output == NULL) {
		printHelpMessage (argv[0], args, "<input matrix filename> <start row> <start col> <rows> <cols> <output matrix filename>", true);
		return -1;
	}

	start_row = atoi (start_row_str);
	start_col = atoi (start_col_str);
	rows = atoi (rows_str);
	cols = atoi (cols_str);

	FileFormatTag input_format = get_format_tag (inputFileFormat);
	FileFormatTag output_format = get_format_tag (outputFileFormat);

	if (input_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (output_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid output-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', or 'png')" << std::endl;
		return -1;
	}

	if (input_format == FORMAT_DETECT)
		input_format = guess_format (input);

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
		if (input_format == FORMAT_PNG)
			ring_type = RING_GF2;
		else
			ring_type = RING_MODULAR;
#else // !__LELA_HAVE_LIBPNG
		ring_type = RING_MODULAR;
#endif // __LELA_HAVE_LIBPNG
	}

	if (ring_type == RING_GF2)
		return get_submatrix (GF2 (), input, input_format, start_row, start_col, rows, cols, output, output_format);
	else if (ring_type == RING_MODULAR) {
		if (floatingPoint) {
			if (ModularTraits<float>::valid_modulus (p))
				return get_submatrix (Modular<float> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
			else if (ModularTraits<double>::valid_modulus (p))
				return get_submatrix (Modular<double> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
			else
				return get_submatrix (Modular<integer> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
		} else {
			if (ModularTraits<uint8>::valid_modulus (p))
				return get_submatrix (Modular<uint8> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
			else if (ModularTraits<uint16>::valid_modulus (p))
				return get_submatrix (Modular<uint16> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
			else if (ModularTraits<uint32>::valid_modulus (p))
				return get_submatrix (Modular<uint32> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
			else
				return get_submatrix (Modular<integer> (p), input, input_format, start_row, start_col, rows, cols, output, output_format);
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
