/* lela/blas/context.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_CONTEXT_H
#define __BLAS_CONTEXT_H

namespace LELA
{

/** Generic module
 *
 * This module enables the naive algorithms used when nothing else is
 * available.
 *
 * Included in every module-structure should be an empty structure
 * called Tag. Inside Tag should be a reference to the tag of the
 * parent-module. With this mechanism the BLAS-methods are specialised
 * to different modules.
 *
 * All module-structures must have a constructor which takes a const
 * Ring.
 */
template <class Ring>
struct GenericModule
{
	struct Tag {};

	GenericModule (const Ring &R) {}
	GenericModule () {}
};

/** All modules
 *
 * This module enables all modules available for the given ring. It
 * should be specialised for each ring based on what is available.
 */
template <class Ring>
struct AllModules : public GenericModule<Ring>
{
	struct Tag { typedef typename GenericModule<Ring>::Tag Parent; };

	AllModules (const Ring &R) : GenericModule<Ring> (R) {}
};

/** Context-object
 *
 * This object encapsulates parameters required to do calculations on
 * a matrix. It stores the ring of computation, settings (such as the
 * cutoff for fast matrix-multiplication), and temporary storage.
 *
 * This object is not designed to be reentrant, therefore there must
 * be at least one Context for each thread. However, under that
 * condition, all arithmetic in this library is reentrant
 *
 * Contexts use a @ref Module object, which allow the selection of
 * different implementations of arithmetic. To create a Context object
 * supporting multiple Modules, just create a single structure which
 * inherits each of the modules and instantiate the Context on that
 * structure. The default is @ref AllModules, which uses all modules
 * available.
 */
template <class Ring, class Modules = AllModules<Ring> >
class Context
{
public:
	const Ring &F;
	Modules M;

	/// Construct a Context from a ring
	Context (const Ring &_F) : F (_F), M (_F) {}

	/// Copy-constructor
	Context (const Context &ctx) : F (ctx.F), M (ctx.M) {}
};

/// @name Enumerations used in arithmetic operations
//@{

/** Upper versus lower triangular matrices, for trmv, trsv, trmm, trsm */
enum TriangularMatrixType {
	UpperTriangular, LowerTriangular
};

/** Whether a matrix is on the left or the right side, for trmv, trsv, trmm, trsm, and permute */
enum MatrixSide {
	LeftSide, RightSide
};

//@} Enumerations used in arithmetic operations

} // namespace LELA

#endif // __BLAS_CONTEXT_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
