/* lela/blas/context.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
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
 *
 * \ingroup blas
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
 *
 * \ingroup blas
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
 * Contexts are parametrised by Modules, which indicate which
 * arithmetic-operations should be used. The default AllModules<Ring>
 * is sufficient for nearly all purposes; it may be overridden for
 * purposes of testing or comparative benchmarking.
 *
 * \ingroup blas
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
///
/// \ingroup blas
//@{

/** Upper versus lower triangular matrices, for trmv, trsv, trmm, trsm
 *
 * \ingroup blas
 */
enum TriangularMatrixType {
	UpperTriangular, LowerTriangular
};

/** Whether a matrix is on the left or the right side, for trmv, trsv, trmm, trsm, and permute
 *
 * \ingroup blas
 */
enum MatrixSide {
	LeftSide, RightSide
};

//@} Enumerations used in arithmetic operations

/** Notes on structure of low-level interface
 *
 * Each BLAS-routine has its own class whose name is the name of the
 * routine preceeded by an underscore. Within the routine is a single
 * static method called op. The classes are parametrised by ring and
 * an empty structure called ModulesTag, which can be retrieved from a
 * Modules-class with Modules::Tag. The method is then parametrised
 * normally by Modules (whose tag need not be ModulesTag, in case
 * ModulesTag refers to a less specialised implementation) and the
 * vector- and matrix-types.
 *
 * To create a specialised implementation or a higher level algorithm,
 * just create the corresponding Modules class with an empty Tag and
 * partially specialise the corresponding class to that tag. If the
 * implementation only works with a given ring, then the class may
 * also be specialised to the corresponding ring.
 *
 * The default implementation of the method op simply retrieves the
 * parent of the tag (accessible via ModulesTag::Parent) and invokes
 * the same method with the same operations on the corresponding class
 * instanted with that tag. This permits that a given module,
 * identified by a certain ModulesTag, need not implement all
 * BLAS-routines -- when it fails to implement a routine, a call to
 * that routine is just passed to the parent.
 *
 * The reason to put these methods into classes, as opposed to using
 * naked functions, is to ensure that the correct specialisation is
 * always invoked. Otherwise the compiler may in some cases (which are
 * notoriously tricky and hard to predict) actually prefer a more
 * generic parametrised method to a more specialised one, perhaps even
 * invoking the most generic method in a never ending recursion.
 *
 * When a method must be further specified by, say, vector-type or
 * matrix-iterator-type, these specialisations should appear as
 * private static methods invoked by the method op.
 *
 * \ingroup blas
 */

} // namespace LELA

#endif // __BLAS_CONTEXT_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
