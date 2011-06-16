/* linbox/blas/level3-ll.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-3 interface, to be used by implementation-functions
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_LL_H
#define __BLAS_LEVEL3_LL_H

#include "linbox/blas/context.h"
#include "linbox/matrix/io.h"

namespace LinBox
{

/** This namespace contains the level 3 BLAS interface */
namespace BLAS3
{

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_copy (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B);

template <class Ring, class Modules, class Matrix>
Matrix &_scal (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A);

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_axpy (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B);

template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C);

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trmm (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne);

template <class Ring, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trsm (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &_permute_rows (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A);

template <class Ring, class Modules, class Iterator, class Matrix>
Matrix &_permute_cols (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A);

template <class Ring, class Modules, class Matrix1, class Matrix2>
bool _equal (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B);

template <class Ring, class Modules, class Matrix>
bool _is_zero (const Ring &F, Modules &M, const Matrix &A);

template <class Ring, class Modules, class Matrix>
std::istream &_read (const Ring &F, Modules &M, std::istream &is, Matrix &A, FileFormatTag format = FORMAT_DETECT);

template <class Ring, class Modules, class Matrix>
std::ostream &_write (const Ring &F, Modules &M, std::ostream &os, const Matrix &A, FileFormatTag format);

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_LL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
