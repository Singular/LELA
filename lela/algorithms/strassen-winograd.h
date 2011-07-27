/* lela/algorithms/strassen-winograd.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementation of Straßen-Winograd fast matrix multiplication based on
 *
 * Boyer, B., Dumas, J.-G., Pernet, C., & Zhou, W. (2007). Memory efficient
 * scheduling of Strassen-Winogradʼs matrix multiplication algorithm.
 * Proceedings of the 2009 international symposium on Symbolic and algebraic
 * computation ISSAC 09, 55. ACM Press. Retrieved from
 * http://arxiv.org/abs/0707.2347
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ALGORITHMS_STRASSEN_WINOGRAD_H
#define __LELA_ALGORITHMS_STRASSEN_WINOGRAD_H

#include <iostream>
#include <iomanip>
#include <cassert>

#include "lela/util/commentator.h"
#include "lela/util/timer.h"
#include "lela/blas/context.h"

namespace LELA
{

/** Implementation of Strassen-Winograd fast matrix-multiplication
 *
 * Based on
 *
 * Boyer, B., Dumas, J.-G., Pernet, C., & Zhou, W. (2007). Memory efficient
 * scheduling of Strassen-Winogradʼs matrix multiplication algorithm.
 * Proceedings of the 2009 international symposium on Symbolic and algebraic
 * computation ISSAC 09, 55. ACM Press. Retrieved from
 * http://arxiv.org/abs/0707.2347
 *
 * \ingroup algorithms
 */
template <class ParentTag>
class StrassenWinograd
{
	size_t _cutoff;
	bool _use_ip;

	// FIXME
	static const size_t default_cutoff = 2048;

	// Calculate the product of the residual part of the input classically
	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &gemm_res (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
			   size_t m, size_t k, size_t n);

	// Ordinary variants: require extra storage and are fast
	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &mul (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, Matrix3 &C);

	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &addmul (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C);

	// In-place variants: require no extra storage and are not quite as fast
	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &mul_ip (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, Matrix3 &C);

	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &addmul_ip (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C);

	// Overwriting variants: require no extra storage and are fast, but overwrite the inputs A and B
	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &mul_ow (const Ring &R, Modules &M, const typename Ring::Element &a, Matrix1 &A, Matrix2 &B, Matrix3 &C);

	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	Matrix3 &addmul_ow (const Ring &R, Modules &M, const typename Ring::Element &a, Matrix1 &A, Matrix2 &B, const typename Ring::Element &b, Matrix3 &C);

public:
	/** Main constructor
	 *
	 * @param cutoff Size below which to switch to classical multiplication
	 *
	 * @param use_ip Set to true to use the in-place variant,
	 * which requires no additional memory but is slower than the
	 * ordinary variant
	 */
	StrassenWinograd (size_t cutoff = default_cutoff, bool use_ip = false)
		: _cutoff (cutoff), _use_ip (use_ip) {}

	StrassenWinograd (const StrassenWinograd &sw)
		: _cutoff (sw._cutoff), _use_ip (sw._use_ip) {}

	StrassenWinograd &operator = (const StrassenWinograd &sw)
		{ _cutoff = sw._cutoff; _use_ip = sw._use_ip; return *this; }

	/** C <- a A * B + b * C using Strassen-Winograd
	 */
	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	inline Matrix3 &gemm (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
	{
		lela_check (A.coldim () == B.rowdim ());
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());

		if (_use_ip) {
			if (R.isZero (b))
				return mul_ip (R, M, a, A, B, C);
			else
				return addmul_ip (R, M, a, A, B, b, C);
		} else {
			if (R.isZero (b))
				return mul (R, M, a, A, B, C);
			else
				return addmul (R, M, a, A, B, b, C);
		}
	}

	/** C <- a A * B + b * C using Strassen-Winograd
	 *
	 * This variant overwrites the inputs A and B. It requires no
	 * additional memory and is faster than the in-place variant.
	 */
	template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
	inline Matrix3 &gemm_overwrite (const Ring &R, Modules &M, const typename Ring::Element &a, Matrix1 &A, Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
	{
		lela_check (A.coldim () == B.rowdim ());
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());

		if (R.isZero (b))
			return mul_ow (R, M, a, A, B, C);
		else
			return addmul_ow (R, M, a, A, B, b, C);
	}
};

template <class Ring, class ParentModule>
struct StrassenModuleTag { typedef typename ParentModule::Tag Parent; };

template <class Ring, class ParentModule>
struct StrassenModule : public ParentModule
{
	typedef StrassenModuleTag<Ring, ParentModule> Tag;

	StrassenWinograd<typename ParentModule::Tag> sw;

	StrassenModule (const Ring &R) : ParentModule (R) {}
};

} // namespace LELA

#include "lela/algorithms/strassen-winograd.tcc"

#endif // __LELA_ALGORITHMS_STRASSEN_WINOGRAD_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
