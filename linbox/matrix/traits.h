/* linbox/matrix/traits.h
 * Copyright 2011 Bradford Hovinen
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_MATRIX_TRAITS_H
#define __LINBOX_MATRIX_TRAITS_H

namespace LinBox {

/** \brief For specializing matrix arithmetic
 *
 * This class defines matrix categories that allow us to specialize the matrix
 * arithmetic in \ref MatrixDomain for different matrix representations. For
 * example, a sparse matrix may have an efficient iterator over row vectors but
 * not over column vectors. Therefore, an algorithm that tries to iterate over
 * column vectors will run very slowly. Hence a specialization that avoids using
 * column vectors is used instead.
 */

struct MatrixCategories 
{
	struct BlackboxTag { };
	struct RowMatrixTag : public virtual BlackboxTag { };
	struct ColMatrixTag : public virtual BlackboxTag { };
	struct RowColMatrixTag : public RowMatrixTag, public ColMatrixTag { };
	struct ZeroOneRowMatrixTag : public RowMatrixTag { };
	struct ZeroOneColMatrixTag : public ColMatrixTag { };
};

template <class Matrix> struct MatrixTraits
{
	typedef Matrix MatrixType;
	typedef typename Matrix::MatrixCategory MatrixCategory;
};

/** Trait describing what iterators the matrix provides
 *
 * This gets around the fact that the C++-compiler only instantiates a
 * partial specialisation when the type exactly matches that in the
 * specification, not when the type inherits it.
 */

template <class Trait> struct MatrixIteratorTypes
{
	typedef Trait MatrixCategory;
};

template <> struct MatrixIteratorTypes<MatrixCategories::ZeroOneRowMatrixTag>
{
	typedef MatrixCategories::RowMatrixTag MatrixCategory;
};

template <> struct MatrixIteratorTypes<MatrixCategories::ZeroOneColMatrixTag>
{
	typedef MatrixCategories::ColMatrixTag MatrixCategory;
};

} // namespace LinBox

#endif // __LINBOX_MATRIX_TRAITS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
