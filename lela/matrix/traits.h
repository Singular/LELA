/* lela/matrix/traits.h
 * Copyright 2011 Bradford Hovinen
 *
 * Properties of matrices on which algorithms are specialised
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_TRAITS_H
#define __LELA_MATRIX_TRAITS_H

namespace LELA
{

/** Matrix-iterator-types
 *
 * These tags indicate which iterators a matrix supports: generic
 * (unspecified), row-iterators only, column-iterators only, or both
 * row- and column-iterators.
 *
 * \ingroup matrix
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
 * \ingroup matrix
 */
namespace MatrixStorageTypes
{
	/** Unspecified storage-type
	 *
	 * Generic means that no assumptions are made about storage. This can
	 * be used for "virtual" matrices which do not physically exist in
	 * memory.
	 */
	struct Generic {};

	/** Storage by rows
	 *
	 * Rows means that the matrix maintains a vector of row-vectors.
	 */
	struct Rows : public Generic {};

	/** Dense storage
	 *
	 * Dense means that the matrix is stored as an array of elements. The
	 * matrix should provide the method disp () which indicates the
	 * displacement from one row to the next in the array.
	 */
	struct Dense : public Generic {};

	/** Transposed dense storage
	 *
	 * DenseTranspose is similar to dense, but the matrix is stored in
	 * column-major order rather than row-major order, i.e. it is the
	 * transpose of an ordinary dense matrix.
	 */
	struct DenseTranspose : public Generic {};

	/** M4RI matrix
	 *
	 * M4RI means that the matrix is a wrapper for a matrix in
	 * libm4ri. This is only meaningful if libm4ri is enabled.
	 */
	struct M4RI : public Generic {};

	/** Transposed M4RI matrix
	 *
	 * M4RITranspose refers to the transpose of a M4RI-matrix. This is
	 * only meaningful if libm4ri is enabled.
	 */
	struct M4RITranspose : public Generic {};
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
