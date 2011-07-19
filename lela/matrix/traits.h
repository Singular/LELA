/* lela/matrix/traits.h
 * Copyright 2011 Bradford Hovinen
 *
 * See COPYING for license information.
 */

#ifndef __LELA_MATRIX_TRAITS_H
#define __LELA_MATRIX_TRAITS_H

namespace LELA {

/** Matrix-iterator-types
 *
 * These tags indicate which iterators a matrix supports: generic
 * (unspecified), row-iterators only, column-iterators only, or both
 * row- and column-iterators.
 */

namespace MatrixIteratorTypes
{
	struct Generic {};
	struct Row : public virtual Generic {};
	struct Col : public virtual Generic {};
	struct RowCol : public Row, public Col {};
}

/** Matrix-storage-types
 *
 * These tags indicate how a matrix is stored.
 *
 * Generic means that no assumptions are made about storage. This can
 * be used for "virtual" matrices which do not physically exist in
 * memory.
 *
 * Rows means that the matrix maintains a vector of row-vectors.
 *
 * Dense means that the matrix is stored as an array of elements. The
 * matrix should provide the method disp () which indicates the
 * displacement from one row to the next in the array.
 *
 * DenseTranspose is similar to dense, but the matrix is stored in
 * column-major order rather than row-major order, i.e. it is the
 * transpose of an ordinary dense matrix.
 *
 * M4RI means that the matrix is a wrapper for a matrix in
 * libm4ri. This is only meaningful if libm4ri is enabled.
 */
namespace MatrixStorageTypes
{
	struct Generic {};
	struct Rows : public Generic {};
	struct Dense : public Generic {};
	struct DenseTranspose : public Generic {};
	struct M4RI : public Generic {};
}

} // namespace LELA

#endif // __LELA_MATRIX_TRAITS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
