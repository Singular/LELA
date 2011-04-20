/* linbox/matrix/io.inl
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_MATRIX_IO_INL
#define __LINBOX_MATRIX_IO_INL

#include <sstream>
#include <cmath>

#include "linbox/util/commentator.h"
#include "linbox/matrix/io.h"

#define BUF_SIZE 80

namespace LinBox
{

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::read (std::istream &is, Matrix &A, FileFormatTag format) const
{
	if (format == FORMAT_DETECT)
		format = detectFormat (is);

	switch (format) {
	case FORMAT_UNKNOWN:
		throw UnrecognisedFormat ();

	case FORMAT_TURNER:
		return readTurner (is, A);
		break;

	case FORMAT_ONE_BASED:
		return readOneBased (is, A);
		break;

	case FORMAT_GUILLAUME:
		return readGuillaume (is, A);
		break;

	case FORMAT_MAPLE:
		return readMaple (is, A);
		break;

	case FORMAT_MATLAB:
		return readMatlab (is, A);
		break;

	case FORMAT_SAGE:
		return readSage (is, A);
		break;

	case FORMAT_PRETTY:
		return readPretty (is, A);
		break;

#ifdef __LINBOX_HAVE_LIBPNG
	case FORMAT_PNG:
		return readPNG (is, A);
		break;
#endif // __LINBOX_HAVE_LIBPNG

	default:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "bad input-format");
	}
}

template <class Field>
FileFormatTag MatrixReader<Field>::detectFormat (std::istream &is) const
{
	throw NotImplemented ();
}

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readTurner (std::istream &is, Matrix &A) const
{
	size_t i, j;

	typename Field::Element x;

	char buf[BUF_SIZE];

	is.getline (buf, 80);

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

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readGuillaume (std::istream &is, Matrix &A) const
{
	size_t m, n, i, j;
	char c;

	is >> m >> n >> c;

	if (c != 'M')
		throw InvalidMatrixInput ();

	A.resize (m, n);

	typename Field::Element x;

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

		if (! _F.isZero(x))
			A.setEntry (i - 1, j - 1, x);
	}

	return is;
}

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readOneBased (std::istream &is, Matrix &A) const
{
	throw NotImplemented ();
}

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readMaple (std::istream &is, Matrix &A) const
{
	throw NotImplemented ();
}

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readMatlab (std::istream &is, Matrix &A) const
{
	size_t i = 0, j = 0;
	char c;
	typename Field::Element a_ij;

	while (1) {
		do
			is >> c;
		while (is && !isdigit (c));

		if (!is)
			break;

		is.putback (c);

		_F.read (is, a_ij);
		A.setEntry (i, j++, a_ij);

		do
			is >> c;
		while (is && c != ',' && c != ';' && c != ']');

		if (!is)
			break;

		if (c == ';') {
			++i;
			j = 0;
		}
		else if (c == ']')
			break;
	}

	return is;
}

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readSage (std::istream &is, Matrix &A) const
{
	throw NotImplemented ();
}

template <class Field>
template <class Matrix>
std::istream &MatrixReader<Field>::readPretty (std::istream &is, Matrix &A) const
{
	size_t i, j;
	typename Field::Element a_ij;
	char c;

	i = 0;

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

template <class Field>
template <class Vector>
void MatrixReader<Field>::appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::DenseVectorTag) const
{
	if (v.size () <= index)
		v.resize (index + 1);

	v[index] = a;
}

template <class Field>
template <class Vector>
void MatrixReader<Field>::appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::DenseZeroOneVectorTag) const
{
	if (v.size () <= index)
		v.resize (index + 1);

	v[index] = a;
}

template <class Field>
template <class Vector>
void MatrixReader<Field>::appendEntrySpecialised (Vector &v, size_t index, const typename Field::Element &a, VectorCategories::HybridZeroOneVectorTag) const
{
	throw NotImplemented ();
}

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::write (std::ostream &os, const Matrix &A, FileFormatTag format) const
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

	case FORMAT_GUILLAUME:
		return writeGuillaume (os, A);
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

#ifdef __LINBOX_HAVE_LIBPNG
	case FORMAT_PNG:
		return writePNG (os, A);
		break;
#endif // __LINBOX_HAVE_LIBPNG

	default:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "bad output-format");
	}

	return os;
}

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writeTurner (std::ostream &os, const Matrix &A) const
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

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writeOneBased (std::ostream &os, const Matrix &A) const
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

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writeGuillaume (std::ostream &os, const Matrix &A) const
{
	os << A.rowdim () << ' ' << A.coldim () << " M" << std::endl;
	writeOneBased (os, A);
	os << "0 0 0" << std::endl;
	return os;
}

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writeMaple (std::ostream &os, const Matrix &A) const
{
	size_t i, j;

	typename Field::Element a;

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

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writeMatlab (std::ostream &os, const Matrix &A) const
{
	size_t i, j;

	typename Field::Element a;

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

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writeSage (std::ostream &os, const Matrix &A) const
{
	size_t i, j;

	if (A.rowdim () == 0) {
		os << "matrix([])" << std::endl;
		return os;
	}

	os << "matrix([ ";

	typename Field::Element a;

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

template <class Field>
template <class Matrix>
std::ostream &MatrixWriter<Field>::writePretty (std::ostream &os, const Matrix &A) const
{
	integer c;

	_F.characteristic (c);
	unsigned col_width = (int) ceil (log (c.get_d ()) / M_LN10);

	typename Field::Element a;

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

} // namespace LinBox

#ifdef __LINBOX_HAVE_LIBPNG
#  include "linbox/matrix/io-png.inl"
#endif // __LINBOX_HAVE_LIBPNG

#endif // __LINBOX_MATRIX_IO_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
