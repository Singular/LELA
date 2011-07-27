/* lela/blas/level1-cblas.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * CBLAS-wrapper for level 1 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_CBLAS_H
#define __BLAS_LEVEL1_CBLAS_H

#include <algorithm>
#include <iostream>

#ifndef __LELA_BLAS_AVAILABLE
#  error "A working installation of BLAS is required to use this header-file."
#endif // __LELA_BLAS_AVAILABLE

#include "lela/cblas.h"

#include "lela/blas/context.h"
#include "lela/vector/traits.h"
#include "lela/ring/type-wrapper.h"
#include "lela/blas/level1-ll.h"

namespace LELA
{

namespace BLAS1 
{

template <>
class _dot<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static float &dot_impl (const TypeWrapperRing<float> &F, Modules &M, float &res, const Vector1 &x, const Vector2 &y,
			       VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _dot<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, res, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static float &dot_impl (const TypeWrapperRing<float> &F, Modules &M, float &res, const Vector1 &x, const Vector2 &y,
			       VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		lela_check (x.size () == y.size ());
		res = cblas_sdot (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return res;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static float &op (const TypeWrapperRing<float> &F, Modules &M, float &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<TypeWrapperRing<float>, Vector1>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<float>, Vector1>::StorageType (),
				   typename VectorTraits<TypeWrapperRing<float>, Vector2>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<float>, Vector2>::StorageType ()); }

	template <class Modules, class Vector1, class Vector2>
	static double &op (const TypeWrapperRing<float> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y)
		{ return _dot<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, res, x, y); }
};

template <>
class _dot<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static double &dot_impl (const TypeWrapperRing<double> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _dot<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, res, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static double &dot_impl (const TypeWrapperRing<double> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		lela_check (x.size () == y.size ());
		res = cblas_ddot (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return res;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static double &op (const TypeWrapperRing<double> &F, Modules &M, double &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<TypeWrapperRing<double>, Vector1>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<double>, Vector1>::StorageType (),
				   typename VectorTraits<TypeWrapperRing<double>, Vector2>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<double>, Vector2>::StorageType ()); }
};

#if 0 // Disabled, since cblas_scopy and cblas_dcopy don't exist

template <>
class _copy<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const TypeWrapperRing<float> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _copy<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const TypeWrapperRing<float> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		lela_check (x.size () == y.size ());
		cblas_scopy (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const TypeWrapperRing<float> &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return copy_impl (F, M, x, y,
				    typename VectorTraits<TypeWrapperRing<float>, Vector1>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector1>::StorageType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector2>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector2>::StorageType ()); }
};

template <>
class _copy<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const TypeWrapperRing<double> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _copy<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const TypeWrapperRing<double> &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		lela_check (x.size () == y.size ());
		cblas_dcopy (x.size (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const TypeWrapperRing<double> &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return copy_impl (F, M, x, y,
				    typename VectorTraits<TypeWrapperRing<double>, Vector1>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector1>::StorageType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector2>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector2>::StorageType ()); }
};

#endif // Disabled

template <>
class _axpy<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _axpy<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		lela_check (x.size () == y.size ());
		cblas_saxpy (x.size (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, Vector2 &y)
		{ return axpy_impl (F, M, a, x, y,
				    typename VectorTraits<TypeWrapperRing<float>, Vector1>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector1>::StorageType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector2>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector2>::StorageType ()); }
};

template <>
class _axpy<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Generic, VectorStorageTypes::Generic, VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _axpy<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, x, y); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorStorageTypes::Real, VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		lela_check (x.size () == y.size ());
		cblas_daxpy (x.size (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, Vector2 &y)
		{ return axpy_impl (F, M, a, x, y,
				    typename VectorTraits<TypeWrapperRing<double>, Vector1>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector1>::StorageType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector2>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector2>::StorageType ()); }
};

template <>
class _scal<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Vector>
	static Vector &scal_impl (const TypeWrapperRing<float> &F, Modules &M, float a, Vector &x,
				  VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _scal<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, x); }

	template <class Modules, class Vector>
	static Vector &scal_impl (const TypeWrapperRing<float> &F, Modules &M, float a, Vector &x,
				  VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		cblas_sscal (x.size (), a, &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Vector>
	static Vector &op (const TypeWrapperRing<float> &F, Modules &M, float a, Vector &x)
		{ return scal_impl (F, M, a, x,
				    typename VectorTraits<TypeWrapperRing<float>, Vector>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector>::StorageType ()); }
};

template <>
class _scal<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Vector>
	static Vector &scal_impl (const TypeWrapperRing<double> &F, Modules &M, double a, Vector &x,
				  VectorRepresentationTypes::Generic, VectorStorageTypes::Generic)
		{ return _scal<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, x); }

	template <class Modules, class Vector>
	static Vector &scal_impl (const TypeWrapperRing<double> &F, Modules &M, double a, Vector &x,
				  VectorRepresentationTypes::Dense, VectorStorageTypes::Real)
	{
		cblas_dscal (x.size (), a, &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Vector>
	static Vector &op (const TypeWrapperRing<double> &F, Modules &M, double a, Vector &x)
		{ return scal_impl (F, M, a, x,
				    typename VectorTraits<TypeWrapperRing<double>, Vector>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector>::StorageType ()); }
};

} // namespace BLAS1

} // namespace LELA

#endif // __BLAS_LEVEL1_CBLAS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
