/* lela/matrix/io.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Support for reading and writing matrices to and from streams
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_IO_TCC
#define __LELA_MATRIX_IO_TCC

#include <sstream>
#include <cmath>
#include <cctype>

#include <regex.h>

#include "lela/lela-config.h"
#include "lela/util/error.h"
#include "lela/util/commentator.h"
#include "lela/matrix/io.h"
#include "lela/blas/level3.h"

#define BUF_SIZE 32

namespace LELA
{

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::read (std::istream &is, Matrix &A, FileFormatTag format) const
{
	if (format == FORMAT_DETECT)
		format = detectFormat (is);

	switch (format) {
	case FORMAT_UNKNOWN:
		throw UnrecognisedFormat ();

	case FORMAT_TURNER:
		return readTurner (is, A);

	case FORMAT_ONE_BASED:
		return readOneBased (is, A);

	case FORMAT_DUMAS:
		return readDumas (is, A);

	case FORMAT_MAPLE:
		return readMaple (is, A);

	case FORMAT_MATLAB:
		return readMatlab (is, A);

	case FORMAT_SAGE:
		return readSage (is, A);

	case FORMAT_PRETTY:
		return readPretty (is, A);

#ifdef __LELA_HAVE_LIBPNG
	case FORMAT_PNG:
		return readPNG (is, A);
#endif // __LELA_HAVE_LIBPNG

	default:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "bad input-format");
	}
}

template <class Ring>
bool MatrixReader<Ring>::isDumas (char *buf, std::streamsize n)
{
	regex_t re;

	if (regcomp (&re, "^[:digit:]+ [:digit:]+ M$", REG_EXTENDED) != 0)
		throw LELAError ("regcomp failure (isDumas)");

	return regexec (&re, buf, 0, NULL, 0) == 0;
}

template <class Ring>
bool MatrixReader<Ring>::isTurner (char *buf, std::streamsize n)
{
	regex_t re;

	if (regcomp (&re, "^[:digit:]+ [:digit:]+ [:digit:]+$", REG_EXTENDED) != 0)
		throw LELAError ("regcomp failure (isTurner)");

	return regexec (&re, buf, 0, NULL, 0) == 0;
}

template <class Ring>
bool MatrixReader<Ring>::isMaple (char *buf, std::streamsize n)
{
	regex_t re;

	if (regcomp (&re, "(\\[[:space:]+\\])|(\\[[:space:]+\\[[:space:]+([:digit:]+[:space:]+,[:space:]+)*[:digit:]+)", REG_EXTENDED) != 0)
		throw LELAError ("regcomp failure (isMaple)");

	return regexec (&re, buf, 0, NULL, 0) == 0;
}

template <class Ring>
bool MatrixReader<Ring>::isMatlab (char *buf, std::streamsize n)
{
	regex_t re;

	if (regcomp (&re, "(\\[[:space:]+\\])|(\\[[:space:]+([:digit:]+[:space:]+,[:space:]+)*[:digit:]+)", REG_EXTENDED) != 0)
		throw LELAError ("regcomp failure (isMatlab)");

	return regexec (&re, buf, 0, NULL, 0) == 0;
}

template <class Ring>
bool MatrixReader<Ring>::isSage (char *buf, std::streamsize n)
{
	regex_t re;

	if (regcomp (&re, "matrix\\((GF\\([:digit:]+\\),|ZZ|QQ)?[ \t\n]+\\[[ \t\n]+\\[[:space:]+([:digit:]+[:space:]+,[:space:]+)*[:digit:]+", REG_EXTENDED) != 0)
		throw LELAError ("regcomp failure (isSage)");

	return regexec (&re, buf, 0, NULL, 0) == 0;
}

template <class Ring>
bool MatrixReader<Ring>::isPretty (char *buf, std::streamsize n)
{
	regex_t re;

	if (regcomp (&re, "((\\[[:space:]+\\])|(\\[[:space:]+((\\.|[:digit:]+)[:space:]+)*(\\.|[:digit:]+)))", REG_EXTENDED) != 0)
		throw LELAError ("regcomp failure (isPretty)");

	return regexec (&re, buf, 0, NULL, 0) == 0;
}

template <class Ring>
FileFormatTag MatrixReader<Ring>::detectFormat (std::istream &is)
{
	FileFormatTag format = FORMAT_UNKNOWN;

	char line[BUF_SIZE];

	is.get (line, BUF_SIZE);

	if (isDumas (line, BUF_SIZE))
		format = FORMAT_DUMAS;
	else if (isTurner (line, BUF_SIZE))
		format = FORMAT_TURNER;
	else if (isMaple (line, BUF_SIZE))
		format = FORMAT_MAPLE;
	else if (isMatlab (line, BUF_SIZE))
		format = FORMAT_MATLAB;
	else if (isSage (line, BUF_SIZE))
		format = FORMAT_SAGE;
	else if (isPretty (line, BUF_SIZE))
		format = FORMAT_PRETTY;
#ifdef __LELA_HAVE_LIBPNG
	else {
		std::streamsize s = is.gcount ();
		while (s-- != 0)
			is.unget ();

		if (isPNG (is))
			format = FORMAT_PNG;
	}
#endif // __LELA_HAVE_LIBPNG

	std::streamsize s = is.gcount ();
	while (s-- != 0)
		is.unget ();

	return format;
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readTurner (std::istream &is, Matrix &A) const
{
	size_t i, j;

	Context<Ring> ctx (_F);
	BLAS3::scal (ctx, _F.zero (), A);

	typename Ring::Element x;

	char buf[BUF_SIZE];

	is.getline (buf, BUF_SIZE);

	do {
		std::istringstream str (buf);

		str >> i;

		if (i == (size_t) -1) break; // return also if row index is -1

		str >> j;

		_F.read (str, x);

		if (! _F.isZero(x))
			A.setEntry (i, j, x);

		is.getline (buf, 80);
	} while (is);

	return is;
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readDumas (std::istream &is, Matrix &A) const
{
	size_t m, n, i, j;
	char c;

	is >> m >> n >> c;

	if (c != 'M')
		throw InvalidMatrixInput ();

	A.resize (m, n);

	Context<Ring> ctx (_F);
	BLAS3::scal (ctx, _F.zero (), A);

	typename Ring::Element x;

	while (is >> i) {
		if (i == 0 || i == (size_t) -1) {
			is >> j;
			_F.read(is, x);
			break;
		}

		is >> j;

		if (i > A.rowdim () || j > A.coldim ())
			throw InvalidMatrixInput ();

		_F.read (is, x);

		if (! _F.isZero (x))
			A.setEntry (i - 1, j - 1, x);
	}

	return is;
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readOneBased (std::istream &is, Matrix &A) const
{
	throw NotImplemented ();
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readMaple (std::istream &is, Matrix &A) const
{
	throw NotImplemented ();
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readMatlab (std::istream &is, Matrix &A) const
{
	size_t i = 0, j = 0;
	char c;
	typename Ring::Element a_ij;

	Context<Ring> ctx (_F);
	BLAS3::scal (ctx, _F.zero (), A);

	do
		is >> c;
	while (isspace (c));

	if (c != '[')
		throw InvalidMatrixInput ();

	while (1) {
		do
			is >> c;
		while (isspace (c));

		if (!is || !(isdigit (c) || c == '-'))
			throw InvalidMatrixInput ();

		is.putback (c);

		_F.read (is, a_ij);

		if (!_F.isZero (a_ij))
			A.setEntry (i, j, a_ij);

		do
			is >> c;
		while (isspace (c));

		if (c == ';') {
			++i;
			j = 0;
		}
		else if (c == ']')
			break;
		else if (c == ',')
			++j;
		else
			throw InvalidMatrixInput ();
	}

	return is;
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readSage (std::istream &is, Matrix &A) const
{
	throw NotImplemented ();
}

template <class Ring>
template <class Matrix>
std::istream &MatrixReader<Ring>::readPretty (std::istream &is, Matrix &A) const
{
	size_t i, j;
	typename Ring::Element a_ij;
	char c;

	i = 0;

	Context<Ring> ctx (_F);
	BLAS3::scal (ctx, _F.zero (), A);

	while (!is.eof ()) {
		while (isspace (is.peek ()))
			is.ignore (1);

		if (is.eof ())
			break;

		is >> c;

		if (c != '[')
			throw InvalidMatrixInput ();

		j = 0;

		while (!is.eof ()) {
			while (isspace (is.peek ()))
				is.ignore (1);

			c = is.peek ();

			switch (c) {
			case ']':
				is >> c;
				goto end_of_line;

			case '.':
				is >> c;
				j++;
				continue;

			default:
				_F.read (is, a_ij);

				if (!_F.isZero (a_ij))
					A.setEntry (i, j, a_ij);

				j++;
			}
		}

	end_of_line:

		i++;
	}

	return is;
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Dense) const
{
	if (v.size () <= index)
		v.resize (index + 1);

	v[index] = a;
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Dense01) const
{
	if (v.size () <= index)
		v.resize (index + 1);

	v[index] = a;
}

template <class Ring>
template <class Vector>
void MatrixReader<Ring>::appendEntrySpecialised (Vector &v, size_t index, const typename Ring::Element &a, VectorRepresentationTypes::Hybrid01) const
{
	throw NotImplemented ();
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::write (std::ostream &os, const Matrix &A, FileFormatTag format) const
{
	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	switch (format) {
	case FORMAT_TURNER:
		return writeTurner (os, A);
		break;

	case FORMAT_ONE_BASED:
		return writeOneBased (os, A);
		break;

	case FORMAT_DUMAS:
		return writeDumas (os, A);
		break;

	case FORMAT_MAPLE:
		return writeMaple (os, A);
		break;

	case FORMAT_MATLAB:
		return writeMatlab (os, A);
		break;

	case FORMAT_SAGE:
		return writeSage (os, A);
		break;

	case FORMAT_PRETTY:
		return writePretty (os, A);
		break;

#ifdef __LELA_HAVE_LIBPNG
	case FORMAT_PNG:
		return writePNG (os, A);
		break;
#endif // __LELA_HAVE_LIBPNG

	default:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "bad output-format");
	}

	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writeTurner (std::ostream &os, const Matrix &A) const
{
	typename Matrix::ConstRawIterator i_elt;
	typename Matrix::ConstRawIndexedIterator i_idx;

	for (i_idx = A.rawIndexedBegin (), i_elt = A.rawBegin (); i_idx != A.rawIndexedEnd (); ++i_idx, ++i_elt) {
		if (!_F.isZero (*i_elt)) {
			os << i_idx->first << ' ' << i_idx->second << ' ';
			_F.write (os, *i_elt) << std::endl;
		}
	}

	os << "-1" << std::endl;

	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writeOneBased (std::ostream &os, const Matrix &A) const
{
	typename Matrix::ConstRawIterator i_elt;
	typename Matrix::ConstRawIndexedIterator i_idx;

	for (i_idx = A.rawIndexedBegin (), i_elt = A.rawBegin (); i_idx != A.rawIndexedEnd (); ++i_idx, ++i_elt) {
		if (!_F.isZero (*i_elt)) {
			os << i_idx->first + 1 << ' ' << i_idx->second + 1 << ' ';
			_F.write (os, *i_elt) << std::endl;
		}
	}

	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writeDumas (std::ostream &os, const Matrix &A) const
{
	os << A.rowdim () << ' ' << A.coldim () << " M" << std::endl;
	writeOneBased (os, A);
	os << "0 0 0" << std::endl;
	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writeMaple (std::ostream &os, const Matrix &A) const
{
	size_t i, j;

	typename Ring::Element a;

	if (A.rowdim () == 0) {
		os << "[]";
		return os;
	}

	os << "[ ";

	for (i = 0; i < A.rowdim (); ++i) {
		os << "[";

		for (j = 0; j < A.coldim (); ++j) {
			if (A.getEntry (a, i, j))
				_F.write (os, a);
			else
				os << "0";

			if (j < A.coldim () - 1)
				os << ", ";
		}

		if (i < A.rowdim () - 1)
			os << "], ";
		else
			os << "] ]" << std::endl;
	}

	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writeMatlab (std::ostream &os, const Matrix &A) const
{
	size_t i, j;

	typename Ring::Element a;

	if (A.rowdim () == 0) {
		os << "[]";
		return os;
	}

	os << "[ ";

	for (i = 0; i < A.rowdim (); ++i) {
		for (j = 0; j < A.coldim (); ++j) {
			if (A.getEntry (a, i, j))
				_F.write (os, a);
			else
				os << "0";

			if (j < A.coldim () - 1)
				os << ", ";
		}

		if (i < A.rowdim () - 1)
			os << "; ";
	}

	os << "]" << std::endl;

	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writeSage (std::ostream &os, const Matrix &A) const
{
	size_t i, j;

	if (A.rowdim () == 0) {
		os << "matrix([])" << std::endl;
		return os;
	}

	os << "matrix([ ";

	typename Ring::Element a;

	for (i = 0; i < A.rowdim (); ++i) {
		if (i == 0)
			os << "[ ";
		else
			os << "         [ ";

		for (j = 0; j < A.coldim (); ++j) {
			if (A.getEntry (a, i, j))
				_F.write (os, a);
			else
				os << "0";

			if (j < A.coldim () - 1)
				os << ", ";
		}

		if (i < A.rowdim () - 1)
			os << " ]," << std::endl;
		else
			os << " ] ])" << std::endl;
	}

	return os;
}

template <class Ring>
template <class Matrix>
std::ostream &MatrixWriter<Ring>::writePretty (std::ostream &os, const Matrix &A) const
{
	integer c;

	unsigned col_width = _F.elementWidth ();

	typename Ring::Element a;

	size_t i, j;

	if (A.rowdim () == 0 || A.coldim () == 0) {
		os << "(empty " << A.rowdim () << "x" << A.coldim () << " matrix)" << std::endl;
		return os;
	}

	for (i = 0; i < A.rowdim (); ++i) {
		os << "  [ ";

		for (j = 0; j < A.coldim (); ++j) {
			if (!A.getEntry (a, i, j) || _F.isZero (a)) {
				for (unsigned t = 0; t < col_width - 1; ++t)
					os << ' ';

				os << ". ";
			} else {
				os.width (col_width);
				_F.write (os, a) << " ";
			}
		}

		os << "]" << std::endl;
	}

	return os;
}

} // namespace LELA

#ifdef __LELA_HAVE_LIBPNG
#  include "lela/matrix/io-png.tcc"
#endif // __LELA_HAVE_LIBPNG

#endif // __LELA_MATRIX_IO_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
