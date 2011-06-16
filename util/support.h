/* util/support.h
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Support-routines for utilities: command-line-processing etc.
 */

#ifndef __LINBOX_UTIL_SUPPORT_H
#define __LINBOX_UTIL_SUPPORT_H

#ifndef __SUPPORT_C
#  include "linbox/ring/gf2.h"
#  include "linbox/matrix/io.h"
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

#ifndef __SUPPORT_C

enum RingType
{
	RING_UNKNOWN, RING_GF2, RING_MODULAR
};

enum MatrixType
{
	MATRIX_UNKNOWN, MATRIX_DENSE, MATRIX_SPARSE, MATRIX_HYBRID
};

RingType get_ring_type (const char *str)
{
	if (!strcmp (str, "gf2"))
		return RING_GF2;
	if (!strcmp (str, "modular"))
		return RING_MODULAR;

	return RING_UNKNOWN;
}

LinBox::FileFormatTag get_format_tag (const char *str)
{
	if (!strcmp (str, "guess"))
		return LinBox::FORMAT_DETECT;
	if (!strcmp (str, "dumas"))
		return LinBox::FORMAT_GUILLAUME;
	if (!strcmp (str, "turner"))
		return LinBox::FORMAT_TURNER;
	if (!strcmp (str, "maple"))
		return LinBox::FORMAT_MAPLE;
	if (!strcmp (str, "matlab"))
		return LinBox::FORMAT_MATLAB;
	if (!strcmp (str, "sage"))
		return LinBox::FORMAT_SAGE;
	if (!strcmp (str, "png"))
		return LinBox::FORMAT_PNG;
	if (!strcmp (str, "pretty"))
		return LinBox::FORMAT_PRETTY;

	return LinBox::FORMAT_UNKNOWN;
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
MatrixType get_matrix_type<LinBox::GF2> (const char *str)
{
	if (!strcmp (str, "dense"))
		return MATRIX_DENSE;
	if (!strcmp (str, "sparse"))
		return MATRIX_SPARSE;
	if (!strcmp (str, "hybrid"))
		return MATRIX_HYBRID;

	return MATRIX_UNKNOWN;
}

#endif // __SUPPORT_C

/* template <class Ring> */
/* MatrixType get_matrix_type (const char *str); */

void parseArguments (int argc, char **argv, Argument *args, int freeArgs, ...);
void printHelpMessage (const char *program, Argument *args, bool printDefaults = false);

#endif // __LINBOX_UTIL_SUPPORT_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
