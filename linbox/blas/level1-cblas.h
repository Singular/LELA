/* linbox/blas/level1-cblas.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * CBLAS-wrapper for level 1 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_CBLAS_H
#define __BLAS_LEVEL1_CBLAS_H

#include <algorithm>
#include <iostream>

#ifndef __LINBOX_BLAS_AVAILABLE
#  error "A working installation of BLAS is required to use this header-file."
#endif // __LINBOX_BLAS_AVAILABLE

#include <cblas.h>

#include "linbox/blas/context.h"
#include "linbox/vector/traits.h"
#include "linbox/ring/unparametric.h"

namespace LinBox
{

namespace BLAS1 
{

template <class Vector1, class Vector2>
float dot_impl (const UnparametricRing<float> &F, BLASModule &M, float &res, const Vector1 &x, const Vector2 &y,
		size_t start_idx, size_t end_idx,
		VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());
	return cblas_sdot (((end_idx == (size_t) -1) ? x.size () : end_idx) - start_idx, &x[start_idx], &x[1] - &x[0], &y[start_idx], &y[1] - &y[0]);
}

template <class Vector1, class Vector2>
double dot_impl (const UnparametricRing<double> &F, BLASModule &M, double &res, const Vector1 &x, const Vector2 &y,
		 size_t start_idx, size_t end_idx,
		 VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());
	return cblas_ddot (((end_idx == (size_t) -1) ? x.size () : end_idx) - start_idx, &x[start_idx], &x[1] - &x[0], &y[start_idx], &y[1] - &y[0]);
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const UnparametricRing<float> &F, BLASModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());
	cblas_scopy (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
	return y;
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const UnparametricRing<double> &F, BLASModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());
	cblas_dcopy (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
	return y;
}

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const UnparametricRing<float> &F, BLASModule &M, float a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());
	cblas_saxpy (x.size (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
	return y;
}

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const UnparametricRing<double> &F, BLASModule &M, double a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());
	cblas_daxpy (x.size (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
	return y;
}

template <class Vector>
Vector &scal_impl (const UnparametricRing<double> &F, BLASModule &M, float a, Vector &x, VectorCategories::DenseVectorTag)
{
	cblas_sscal (x.size (), a, &x[0], &x[1] - &x[0]);
	return x;
}

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_CBLAS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
