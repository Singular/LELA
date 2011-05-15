/* linbox/blas/context.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_H
#define __BLAS_LEVEL3_H

namespace LinBox
{

/** Generic module
 *
 * This module enables the naive algorithms used when nothing else is
 * available.
 *
 * Other computation-modules should inherit this as a virtual public
 * base so that this is always available as a fallback. It should not
 * be inherited when collecting modules together, since that would
 * lead to ambiguous instantiation.
 */
struct GenericModule {};

/** All modules
 *
 * This module enables all modules available for the given field. It
 * should be specialised for each field based on what is available.
 */
template <class Field>
struct AllModules : public GenericModule {};

/** Context-object
 *
 * This object encapsulates parameters required to do calculations on
 * a matrix. It stores the field of computation, settings (such as the
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
template <class Field, class Modules = AllModules<Field> >
class Context
{
public:
	const Field &F;
	Modules M;

	/// Construct a Context from a field
	Context (const Field &_F) : F (_F) {}

	/// Copy-constructor
	Context (const Context &ctx) : F (ctx.F), M (ctx.M) {}
};

/// @name Enumerations used in arithmetic operations
//@{

#if 0
/** Upper versus lower triangular matrices, for trmv, trsv, trmm, trsm */
enum TriangularMatrixType {
	UpperTriangular, LowerTriangular
};
#endif

/** Whether a triangular matrix is on the left or the right side, for trmv, trsv, trmm, trsm */
enum TriangularMatrixSide {
	LeftSide, RightSide
};

//@} Enumerations used in arithmetic operations

} // namespace LinBox

#endif // __BLAS_LEVEL3_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
