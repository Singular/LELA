/* util/row-echelon-form.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Utility to convert a matrix to row-echelon form
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
#include "lela/solutions/echelon-form.h"
#include "lela/solutions/echelon-form-gf2.h"

#include "support.h"

using namespace LELA;

template <class Ring, class Matrix>
int row_echelon_form (const Ring &R, const char *input, FileFormatTag input_format, const char *output, FileFormatTag output_format,
		      typename EchelonForm<Ring>::Method method, bool reduced)
{
	Context<Ring> ctx (R);

	commentator.start ("Converting matrix to row-echelon-form", __FUNCTION__);

	Matrix A;

	commentator.start ("Reading input-matrix");

	std::ifstream ifile (input);

	if (!ifile.good ()) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Could not open input-file" << std::endl;
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

	EchelonForm<Ring> EF (ctx);

	if (reduced)
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Computing reduced form" << std::endl;
	else
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Computing non-reduced form" << std::endl;

	EF.echelonize (A, reduced, method);

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

template <class Ring>
typename EchelonForm<Ring>::Method get_method (const char *str)
{
	if (!strcmp (str, "standard"))
		return EchelonForm<Ring>::METHOD_STANDARD_GJ;
	if (!strcmp (str, "afast"))
		return EchelonForm<Ring>::METHOD_ASYMPTOTICALLY_FAST_GJ;
	if (!strcmp (str, "f4"))
		return EchelonForm<Ring>::METHOD_FAUGERE_LACHARTRE;

	return EchelonForm<Ring>::METHOD_UNKNOWN;
}

template <>
EchelonForm<GF2>::Method get_method<GF2> (const char *str)
{
	if (!strcmp (str, "standard"))
		return EchelonForm<GF2>::METHOD_STANDARD_GJ;
	if (!strcmp (str, "afast"))
		return EchelonForm<GF2>::METHOD_ASYMPTOTICALLY_FAST_GJ;
#ifdef __LELA_HAVE_M4RI
	if (!strcmp (str, "m4ri"))
		return EchelonForm<GF2>::METHOD_M4RI;
#endif // __LELA_HAVE_M4RI
	if (!strcmp (str, "f4"))
		return EchelonForm<GF2>::METHOD_FAUGERE_LACHARTRE;

	return EchelonForm<GF2>::METHOD_UNKNOWN;
}

template <class Ring>
struct ErrorText
{
	static const char *method;
	static const char *type;
};

template <class Ring>
const char *ErrorText<Ring>::method = "Invalid method (use 'standard', 'afast', or 'f4')";

template <class Ring>
const char *ErrorText<Ring>::type = "Invalid matrix-type (use 'dense' or 'sparse')";

template <>
struct ErrorText<GF2>
{
	static const char *method;
	static const char *type;
};

#ifdef __LELA_HAVE_M4RI
const char *ErrorText<GF2>::method = "Invalid method (use 'standard', 'afast', 'm4ri', or 'f4')";
#else
const char *ErrorText<GF2>::method = "Invalid method (use 'standard', 'afast', or 'f4')";
#endif // __LELA_HAVE_M4RI

const char *ErrorText<GF2>::type = "Invalid matrix-type (use 'dense', 'sparse', or 'hybrid')";

template <class Ring>
int inner_run_row_echelon_form (const Ring &R, typename EchelonForm<Ring>::Method method, MatrixType type, bool reduced,
				const char *input, FileFormatTag input_format, const char *output, FileFormatTag output_format)
{
	if (type == MATRIX_DENSE)
		return row_echelon_form<Ring, DenseMatrix<typename Ring::Element> > (R, input, input_format, output, output_format, method, reduced);
	else if (type == MATRIX_SPARSE)
		return row_echelon_form<Ring, SparseMatrix<typename Ring::Element> > (R, input, input_format, output, output_format, method, reduced);
	else {
		std::cerr << ErrorText<Ring>::type << std::endl;
		return -1;
	}
}

int inner_run_row_echelon_form (const GF2 &R, EchelonForm<GF2>::Method method, MatrixType type, bool reduced,
				const char *input, FileFormatTag input_format, const char *output, FileFormatTag output_format)
{
	if (type == MATRIX_DENSE)
		return row_echelon_form<GF2, DenseMatrix<GF2::Element> > (R, input, input_format, output, output_format, method, reduced);
	else if (type == MATRIX_SPARSE)
		return row_echelon_form<GF2, SparseMatrix<GF2::Element, Vector<GF2>::Sparse> > (R, input, input_format, output, output_format, method, reduced);
	else if (type == MATRIX_HYBRID)
		return row_echelon_form<GF2, SparseMatrix<GF2::Element, Vector<GF2>::Hybrid> > (R, input, input_format, output, output_format, method, reduced);
	else {
		std::cerr << ErrorText<GF2>::type << std::endl;
		return -1;
	}
}

template <class Ring>
int run_row_echelon_form (const Ring &R, const char *method_string, const char *matrix_type_string, bool reduced,
			  const char *input, FileFormatTag input_format, const char *output, FileFormatTag output_format)
{
	typename EchelonForm<Ring>::Method method = get_method<Ring> (method_string);

	if (method == EchelonForm<Ring>::METHOD_UNKNOWN) {
		std::cerr << ErrorText<Ring>::method << std::endl;
		return -1;
	}

	MatrixType type = get_matrix_type<Ring> (matrix_type_string);

	return inner_run_row_echelon_form (R, method, type, reduced, input, input_format, output, output_format);
}

int main (int argc, char **argv)
{
	static bool reduced = false;
	static const char *ringString = "guess";
	static integer p = 65521;
	static bool floatingPoint = false;
	static const char *methodString = "f4";
	static const char *inputFileFormat = "guess";
	static const char *outputFileFormat = "guess";
	static const char *matrixType = "dense";
	static char *input = NULL;
	static char *output = NULL;

	static Argument args[] = {
		{ 'r', "-r", "Compute the reduced row-echelon form", TYPE_NONE, &reduced },
		{ 'k', "-k", "Ring over which to compute ('guess', 'gf2', 'modular')", TYPE_STRING, &ringString },
		{ 'p', "-p", "Modulus of ring, when ring is 'modular'", TYPE_INT, &p },
		{ 'f', "-f", "Compute using floating point, when ring is 'modular'", TYPE_NONE, &floatingPoint },
#ifdef __LELA_HAVE_M4RI
		{ 'm', "-m", "Method to be used ('standard', 'afast', 'm4ri', or 'f4')", TYPE_STRING, &methodString },
#else
		{ 'm', "-m", "Method to be used ('standard', 'afast', or 'f4')", TYPE_STRING, &methodString },
#endif // __LELA_HAVE_M4RI
		{ 'i', "-i", "Input file format ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png', 'pretty')", TYPE_STRING, &inputFileFormat },
		{ 'o', "-o", "Output file format ('guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png', 'pretty')", TYPE_STRING, &outputFileFormat },
		{ 't', "-t", "Type to use for matrix ('dense', 'sparse', 'hybrid')", TYPE_STRING, &matrixType },
		{ '\0' }
	};

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);

	parseArguments (argc, argv, args, "<input-filename> <output-filename>", 2, &input, &output);

	if (input == NULL || output == NULL) {
		printHelpMessage (argv[0], args, "<first matrix filename> <second matrix filename>", true);
		return -1;
	}

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	FileFormatTag input_format = get_format_tag (inputFileFormat);
	FileFormatTag output_format = get_format_tag (outputFileFormat);

	if (input_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid input-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png', or 'pretty')" << std::endl;
		return -1;
	}

	if (output_format == FORMAT_UNKNOWN) {
		std::cerr << "Invalid output-file-format (use 'guess', 'dumas', 'turner', 'maple', 'matlab', 'sage', 'png', or 'pretty')" << std::endl;
		return -1;
	}

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
		if (input_format == FORMAT_PNG || output_format == FORMAT_PNG)
			ring_type = RING_GF2;
		else
			ring_type = RING_MODULAR;
#else // !__LELA_HAVE_LIBPNG
		ring_type = RING_MODULAR;
#endif // __LELA_HAVE_LIBPNG
	}

	if (ring_type == RING_GF2)
		return run_row_echelon_form (GF2 (), methodString, matrixType, reduced, input, input_format, output, output_format);
	else if (ring_type == RING_MODULAR) {
		if (floatingPoint) {
			if (p < 1U << 7)
				return run_row_echelon_form (Modular<float> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
			else if (p < 1U << 21)
				return run_row_echelon_form (Modular<double> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
			else
				return run_row_echelon_form (Modular<integer> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
		} else {
			if (ModularTraits<uint8>::valid_modulus (p))
				return run_row_echelon_form (Modular<uint8> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
			else if (ModularTraits<uint16>::valid_modulus (p))
				return run_row_echelon_form (Modular<uint16> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
			else if (ModularTraits<uint32>::valid_modulus (p))
				return run_row_echelon_form (Modular<uint32> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
			else
				return run_row_echelon_form (Modular<integer> (p), methodString, matrixType, reduced, input, input_format, output, output_format);
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
