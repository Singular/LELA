/* util/support.h
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Support-routines for utilities: command-line-processing etc.
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_UTIL_SUPPORT_H
#define __LELA_UTIL_SUPPORT_H

#ifndef __SUPPORT_C
#  include "lela/ring/gf2.h"
#  include "lela/matrix/io.h"
#endif

enum ArgumentType {
	TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE, TYPE_STRING
};
#define TYPE_BOOL TYPE_NONE

struct Argument 
{
	char             c;
	const char      *example;
	const char      *helpString;
	ArgumentType     type;
	void            *data;
};

extern const char *format_names[];

#ifndef __SUPPORT_C

enum RingType
{
	RING_UNKNOWN, RING_GUESS, RING_GF2, RING_MODULAR
};

enum MatrixType
{
	MATRIX_UNKNOWN, MATRIX_DENSE, MATRIX_SPARSE, MATRIX_HYBRID
};

RingType get_ring_type (const char *str)
{
	if (!strcmp (str, "guess"))
		return RING_GUESS;
	if (!strcmp (str, "gf2"))
		return RING_GF2;
	if (!strcmp (str, "modular"))
		return RING_MODULAR;

	return RING_UNKNOWN;
}

LELA::FileFormatTag get_format_tag (const char *str)
{
	if (!strcmp (str, "guess"))
		return LELA::FORMAT_DETECT;
	if (!strcmp (str, "dumas"))
		return LELA::FORMAT_DUMAS;
	if (!strcmp (str, "turner"))
		return LELA::FORMAT_TURNER;
	if (!strcmp (str, "maple"))
		return LELA::FORMAT_MAPLE;
	if (!strcmp (str, "matlab"))
		return LELA::FORMAT_MATLAB;
	if (!strcmp (str, "sage"))
		return LELA::FORMAT_SAGE;
#ifdef __LELA_HAVE_LIBPNG
	if (!strcmp (str, "png"))
		return LELA::FORMAT_PNG;
#endif // __LELA_HAVE_LIBPNG
	if (!strcmp (str, "pretty"))
		return LELA::FORMAT_PRETTY;

	return LELA::FORMAT_UNKNOWN;
}

// Try to guess the format-tag by looking at the filename
LELA::FileFormatTag guess_format_tag (const char *filename)
{
	const char *filename_ext = strrchr (filename, '.');

	if (filename_ext == NULL || filename_ext[1] == '\0')
		return LELA::FORMAT_UNKNOWN;

	++filename_ext;

#ifdef __LELA_HAVE_LIBPNG
	if (!strcmp (filename_ext, "png"))
		return LELA::FORMAT_PNG;
#endif // __LELA_HAVE_LIBPNG
	if (!strcmp (filename_ext, "m"))
		return LELA::FORMAT_MATLAB;
	if (!strcmp (filename_ext, "sage"))
		return LELA::FORMAT_SAGE;

	return LELA::FORMAT_UNKNOWN;
}

template <class Ring>
MatrixType get_matrix_type (const char *str)
{
	if (!strcmp (str, "dense"))
		return MATRIX_DENSE;
	if (!strcmp (str, "sparse"))
		return MATRIX_SPARSE;

	return MATRIX_UNKNOWN;
}

template <>
MatrixType get_matrix_type<LELA::GF2> (const char *str)
{
	if (!strcmp (str, "dense"))
		return MATRIX_DENSE;
	if (!strcmp (str, "sparse"))
		return MATRIX_SPARSE;
	if (!strcmp (str, "hybrid"))
		return MATRIX_HYBRID;

	return MATRIX_UNKNOWN;
}

LELA::FileFormatTag guess_format (char *filename)
{
	std::ifstream file (filename);

	if (!file.good ())
		return LELA::FORMAT_UNKNOWN;

	return LELA::MatrixReader<LELA::GF2>::detectFormat (file);
}

#endif // __SUPPORT_C

/* template <class Ring> */
/* MatrixType get_matrix_type (const char *str); */

void parseArguments (int argc, char **argv, Argument *args, const char *freeArgsText, int freeArgs, ...);
void printHelpMessage (const char *program, Argument *args, const char *freeArgsText, bool printDefaults = false);

#endif // __LELA_UTIL_SUPPORT_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
