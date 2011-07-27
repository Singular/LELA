/* lela/blas/level3-ll.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-3 interface, to be used by implementation-functions
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_LL_H
#define __BLAS_LEVEL3_LL_H

#include "lela/blas/context.h"
#include "lela/matrix/io.h"

namespace LELA
{

/** This namespace contains the level 3 BLAS interface */
namespace BLAS3
{

template <class Ring, class ModulesTag>
class _copy
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B)
		{ return _copy<Ring, typename ModulesTag::Parent>::op (F, M, A, B); }
};

template <class Ring, class ModulesTag>
class _scal
{
public:
	template <class Modules, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A)
		{ return _scal<Ring, typename ModulesTag::Parent>::op (F, M, a, A); }
};

template <class Ring, class ModulesTag>
class _axpy
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B)
		{ return _axpy<Ring, typename ModulesTag::Parent>::op (F, M, a, A, B); }
};

template <class Ring, class ModulesTag>
class _gemm
{
public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
		{ return _gemm<Ring, typename ModulesTag::Parent>::op (F, M, a, A, B, b, C); }
};

template <class Ring, class ModulesTag>
class _trmm
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trmm<Ring, typename ModulesTag::Parent>::op (F, M, a, A, B, type, diagIsOne); }
};

template <class Ring, class ModulesTag>
class _trsm
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trsm<Ring, typename ModulesTag::Parent>::op (F, M, a, A, B, type, diagIsOne); }
};

template <class Ring, class ModulesTag>
class _permute_rows
{
public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return _permute_rows<Ring, typename ModulesTag::Parent>::op (F, M, P_begin, P_end, A); }
};

template <class Ring, class ModulesTag>
class _permute_cols
{
public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return _permute_cols<Ring, typename ModulesTag::Parent>::op (F, M, P_begin, P_end, A); }
};

template <class Ring, class ModulesTag>
class _equal
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static bool op (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B)
		{ return _equal<Ring, typename ModulesTag::Parent>::op (F, M, A, B); }
};

template <class Ring, class ModulesTag>
class _is_zero
{
public:
	template <class Modules, class Matrix>
	static bool op (const Ring &F, Modules &M, const Matrix &A)
		{ return _is_zero<Ring, typename ModulesTag::Parent>::op (F, M, A); }
};

template <class Ring, class ModulesTag>
class _read
{
public:
	template <class Modules, class Matrix>
	static std::istream &op (const Ring &F, Modules &M, std::istream &is, Matrix &A, FileFormatTag format = FORMAT_DETECT)
		{ return _read<Ring, typename ModulesTag::Parent>::op (F, M, is, A, format); }
};

template <class Ring, class ModulesTag>
class _write
{
public:
	template <class Modules, class Matrix>
	static std::ostream &op (const Ring &F, Modules &M, std::ostream &os, const Matrix &A, FileFormatTag format)
		{ return _write<Ring, typename ModulesTag::Parent>::op (F, M, os, A, format); }
};

} // namespace BLAS3

} // namespace LELA

#include "lela/blas/level3-generic.h"

#endif // __BLAS_LEVEL3_LL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
