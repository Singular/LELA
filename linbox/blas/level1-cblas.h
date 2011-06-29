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
#include "linbox/blas/level1-ll.h"

namespace LinBox
{

namespace BLAS1 
{

template <>
class _dot<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static float dot_impl (const UnparametricRing<float> &F, Modules &M, float &res, const Vector1 &x, const Vector2 &y,
			       size_t start_idx, size_t end_idx,
			       VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _dot<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, res, x, y, start_idx, end_idx); }

	template <class Modules, class Vector1, class Vector2>
	static float dot_impl (const UnparametricRing<float> &F, Modules &M, float &res, const Vector1 &x, const Vector2 &y,
			       size_t start_idx, size_t end_idx,
			       VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		linbox_check (x.size () == y.size ());
		return cblas_sdot (((end_idx == (size_t) -1) ? x.size () : end_idx) - start_idx, &x[start_idx], &x[1] - &x[0], &y[start_idx], &y[1] - &y[0]);
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static float op (const UnparametricRing<float> &F, Modules &M, float &res, const Vector1 &x, const Vector2 &y,
			 size_t start_idx = 0, size_t end_idx = (size_t) -1)
		{ return dot_impl (F, M, res, x, y, start_idx, end_idx,
				   typename VectorTraits<UnparametricRing<float>, Vector1>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<float>, Vector1>::StorageType (),
				   typename VectorTraits<UnparametricRing<float>, Vector2>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<float>, Vector2>::StorageType ()); }
};

template <>
class _dot<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static double dot_impl (const UnparametricRing<double> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y,
				size_t start_idx, size_t end_idx,
				VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _dot<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, res, x, y, start_idx, end_idx); }

	template <class Modules, class Vector1, class Vector2>
	static double dot_impl (const UnparametricRing<double> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y,
				size_t start_idx, size_t end_idx,
				VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		linbox_check (x.size () == y.size ());
		return cblas_ddot (((end_idx == (size_t) -1) ? x.size () : end_idx) - start_idx, &x[start_idx], &x[1] - &x[0], &y[start_idx], &y[1] - &y[0]);
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static double op (const UnparametricRing<double> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y,
			  size_t start_idx = 0, size_t end_idx = (size_t) -1)
		{ return dot_impl (F, M, res, x, y, start_idx, end_idx,
				   typename VectorTraits<UnparametricRing<double>, Vector1>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<double>, Vector1>::StorageType (),
				   typename VectorTraits<UnparametricRing<double>, Vector2>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<double>, Vector2>::StorageType ()); }
};

template <>
class _copy<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const UnparametricRing<float> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _copy<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const UnparametricRing<float> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		linbox_check (x.size () == y.size ());
		cblas_scopy (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const UnparametricRing<float> &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return copy_impl (F, M, x, y,
				    typename VectorTraits<UnparametricRing<float>, Vector1>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector1>::StorageType (),
				    typename VectorTraits<UnparametricRing<float>, Vector2>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector2>::StorageType ()); }
};

template <>
class _copy<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const UnparametricRing<double> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _copy<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const UnparametricRing<double> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		linbox_check (x.size () == y.size ());
		cblas_dcopy (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const UnparametricRing<double> &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return copy_impl (F, M, x, y,
				    typename VectorTraits<UnparametricRing<double>, Vector1>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector1>::StorageType (),
				    typename VectorTraits<UnparametricRing<double>, Vector2>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector2>::StorageType ()); }
};

template <>
class _axpy<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const UnparametricRing<float> &F, Modules &M, float a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _axpy<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const UnparametricRing<float> &F, Modules &M, float a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		linbox_check (x.size () == y.size ());
		cblas_saxpy (x.size (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const UnparametricRing<float> &F, Modules &M, float a, const Vector1 &x, Vector2 &y)
		{ return axpy_impl (F, M, a, x, y,
				    typename VectorTraits<UnparametricRing<float>, Vector1>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector1>::StorageType (),
				    typename VectorTraits<UnparametricRing<float>, Vector2>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector2>::StorageType ()); }
};

template <>
class _axpy<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const UnparametricRing<double> &F, Modules &M, double a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _axpy<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const UnparametricRing<double> &F, Modules &M, double a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		linbox_check (x.size () == y.size ());
		cblas_daxpy (x.size (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const UnparametricRing<double> &F, Modules &M, double a, const Vector1 &x, Vector2 &y)
		{ return axpy_impl (F, M, a, x, y,
				    typename VectorTraits<UnparametricRing<double>, Vector1>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector1>::StorageType (),
				    typename VectorTraits<UnparametricRing<double>, Vector2>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector2>::StorageType ()); }
};

template <>
class _scal<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Vector>
	static Vector &scal_impl (const UnparametricRing<float> &F, Modules &M, float a, Vector &x,
				  VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _scal<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, x); }

	template <class Modules, class Vector>
	static Vector &scal_impl (const UnparametricRing<float> &F, Modules &M, float a, Vector &x,
				  VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		cblas_sscal (x.size (), a, &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Vector>
	static Vector &op (const UnparametricRing<float> &F, Modules &M, float a, Vector &x)
		{ return scal_impl (F, M, a, x,
				    typename VectorTraits<UnparametricRing<float>, Vector>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector>::StorageType ()); }
};

template <>
class _scal<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Vector>
	static Vector &scal_impl (const UnparametricRing<double> &F, Modules &M, double a, Vector &x,
				  VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _scal<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, x); }

	template <class Modules, class Vector>
	static Vector &scal_impl (const UnparametricRing<double> &F, Modules &M, double a, Vector &x,
				  VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		cblas_sscal (x.size (), a, &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Vector>
	static Vector &op (const UnparametricRing<double> &F, Modules &M, double a, Vector &x)
		{ return scal_impl (F, M, a, x,
				    typename VectorTraits<UnparametricRing<double>, Vector>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector>::StorageType ()); }
};

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
