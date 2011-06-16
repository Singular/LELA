/* linbox/vector/stream.h
 * Copyright 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *  
 * See COPYING for license information.
 */

#ifndef __LINBOX_vector_stream_H
#define __LINBOX_vector_stream_H

#include <vector>
#include <cmath>

#include "linbox/vector/traits.h"
#include "linbox/util/debug.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/vector/bit-iterator.h"

namespace LinBox 
{

/** \brief Vector factory

 * This is an abstract base class that generates a sequence of vectors
 * in a generic way. Typical uses would be in tests, where the same test
 * might be run on a sequence of random vectors or on e_1, ..., e_n.
 \ingroup vector
 */
template <class _Vector>
class VectorStream 
{
    public:
	typedef _Vector Vector;
        typedef VectorStream<Vector> Self_t;
     

	virtual ~VectorStream () {}

	/** Get the next vector from the factory and store it in v
	 */

	virtual Vector &get (Vector &v) = 0;

	/** Extraction operator form
	 */
	Self_t &operator >> (Vector &v)
		{ get (v); return *this; }

	/** Get the number of vectors to be constructed in this stream
	 */
	virtual size_t size () const = 0;

	/** Get the number of vectors constructed so far
	 */
	virtual size_t pos () const = 0;

	/** Get the dimension of each vector
	 */
	virtual size_t dim () const = 0;

	/** Return true if and only if the vector stream still has more vectors
	 * to construct
	 */
	virtual operator bool () const = 0;

	/** Reset the vector stream to the beginning.
	 */
	virtual void reset () = 0;

	/** Alias for reset
	 */
	void rewind () { reset (); }
};

/** Constant vector factory
 * Returns the same vector repeatedly
 */
template <class _Vector>
class ConstantVectorStream : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef ConstantVectorStream<Vector> Self_t;

	/** Constructor
	 * Construct a new factory with the given ring and vector size.
	 * @param v Vector to return on next
	 * @param m Number of vectors to return (0 for unlimited)
	 */
	ConstantVectorStream (Vector &v, size_t m) : _v (v), _m (m), _j (0) {}

	/** Retrieve vector
	 * @param v Vector to use
	 */
	Vector &get (Vector &v) 
		{ if (_m == 0 || _j < _m) copy (_v.begin (), _v.end (), v.begin ()); return v; }

	/** Extraction operator form
	 */
	Self_t &operator >> (Vector &v)
		{ get (v); return *this; }

	/** Number of vectors to be created
	 */
	size_t size () const { return _m; }

	/** Number of vectors created so far
	 */
	size_t pos () const { return _j; }

	/** Dimension of the space
	 */
	size_t dim () const { return _v.size (); }

	/** Check whether we have reached the end
	 */
	operator bool () const { return _m == 0 || _j < _m; }

	/** Reset the factory to start at the beginning
	 */
	void reset () { _j = 0; }

    private:
	Vector &_v;
	size_t  _m;
	size_t  _j;
};

/** Random dense vector stream
 * Generates a sequence of random dense vectors over a given ring
 */
template <class Ring, class _Vector = typename LinBox::Vector<Ring>::Dense, class RandIter = typename Ring::RandIter, class Trait = typename VectorTraits<Ring, _Vector>::VectorCategory>
class RandomDenseStream : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomDenseStream<Ring, Vector, RandIter, Trait> Self_t;

	/** Constructor
	 * Construct a new stream with the given ring and vector size.
	 * @param F Ring over which to create random vectors
	 * @param n Size of vectors
	 * @param m Number of vectors to return (0 for unlimited)
	 */
	RandomDenseStream (const Ring &F, size_t n, size_t m = 0);

	/** Constructor
	 * Construct a new stream with the given ring and vector size.
	 * @param F Ring over which to create random vectors
	 * @param n Size of vectors
	 * @param m Number of vectors to return (0 for unlimited)
	 */
	RandomDenseStream (const Ring &F, const RandIter &r, size_t n, size_t m = 0);

	/** Get next element
	 * @param v Vector into which to generate random vector
	 * @return reference to new random vector
	 */
	Vector &get (Vector &v);

	/** Extraction operator form
	 */
	Self_t &operator >> (Vector &v)
		{ get (v); return *this; }

	/** Number of vectors to be created
	 */
	size_t size () const;

	/** Number of vectors created so far
	 */
	size_t pos () const;

	/** Dimension of the space
	 */
	size_t dim () const;

	/** Check whether we have reached the end
	 */
	operator bool () const;

	/** Reset the stream to start at the beginning
	 */
	void reset ();
};

// Specialization of random dense stream for dense vectors

template <class Ring, class _Vector, class RandIter>
class RandomDenseStream<Ring, _Vector, RandIter, VectorCategories::DenseVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomDenseStream<Ring, Vector, RandIter, VectorCategories::DenseVectorTag > Self_t;

	RandomDenseStream (const Ring &F, size_t n, size_t m = 0)
		: _F (F), _r (F), _n (n), _m (m), _j (0)
	{}

	RandomDenseStream (const Ring &F, const RandIter &r, size_t n, size_t m = 0)
		: _F (F), _r (r), _n (n), _m (m), _j (0)
	{}

	Vector &get (Vector &v);

	/** Extraction operator form
	 */
	Self_t &operator >> (Vector &v)
		{ get (v); return *this; }
	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

    private:
	const Ring &_F;
	RandIter     _r;
	size_t       _n;
	size_t       _m;
	size_t       _j;
};

/** Random sparse vector stream
 * Generates a sequence of random sparse vectors over a given ring
 */
template <class Ring, class _Vector = typename LinBox::Vector<Ring>::Sparse, class RandIter = typename Ring::RandIter, class Trait = typename VectorTraits<Ring, _Vector>::VectorCategory>
class RandomSparseStream : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomSparseStream<Ring, Vector, RandIter, Trait > Self_t;

	/** Constructor
	 * Construct a new stream with the given ring and vector size.
	 * @param F Ring over which to create random vectors
	 * @param p Proportion of nonzero entries
	 * @param n Size of vectors
	 * @param m Number of vectors to return (0 for unlimited)
	 */
	RandomSparseStream (const Ring &F, double p, size_t n, size_t m = 0);

	/** Constructor
	 * Construct a new stream with the given ring and vector size.
	 * @param F Ring over which to create random vectors
	 * @param n Size of vectors
	 * @param p Proportion of nonzero entries
	 * @param m Number of vectors to return (0 for unlimited)
	 */
	RandomSparseStream (const Ring &F, const RandIter &r, double p, size_t n, size_t m = 0);

	/** Get next element
	 * @param v Vector into which to generate random vector
	 * @return reference to new random vector
	 */
	Vector &get (Vector &v);

	/** Extraction operator form
	 */
	Self_t &operator >> (Vector &v)
		{ get (v); return *this; }

	/** Number of vectors to be created
	 */
	size_t size () const;

	/** Number of vectors created so far
	 */
	size_t pos () const;

	/** Dimension of the space
	 */
	size_t dim () const;

	/** Check whether we have reached the end
	 */
	operator bool () const;

	/** Reset the stream to start at the beginning
	 */
	void reset ();

	/** Set the probability of a nonzero entry
	 */
	void setP (double p);
};

// Specialization of RandomSparseStream for dense vectors

template <class Ring, class _Vector, class RandIter>
class RandomSparseStream<Ring, _Vector, RandIter, VectorCategories::DenseVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomSparseStream<Ring, Vector, RandIter, VectorCategories::DenseVectorTag> Self_t;

	RandomSparseStream (const Ring &F, double p, size_t n, size_t m = 0)
		: _F (F), _r1 (F), _r (F, _r1),
		  _n (n), _p (p), _m (m), _j (0),
		  _MT (time (NULL))
		{ linbox_check ((p >= 0.0) && (p <= 1.0)); _F.init (_zero, 0); }

	RandomSparseStream (const Ring &F, const RandIter &r, double p, size_t n, size_t m = 0, int seed = 0)
		: _F (F), _r1 (r), _r (F, _r1), _n (n), _p (p), _m (m), _j (0),
		  _MT (seed)
		{ linbox_check ((p >= 0.0) && (p <= 1.0)); _F.init (_zero, 0); }

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }
	void setP (double p) { linbox_check ((p >= 0.0) && (p <= 1.0)); _p = p; }

    private:
	const Ring                      &_F;
	typename Ring::Element           _zero;
	RandIter                          _r1;
	NonzeroRandIter<Ring, RandIter>  _r;
	size_t                            _n;
	double                            _p;
	size_t                            _m;
	size_t                            _j;
	MersenneTwister                   _MT;
};

// Specialization of RandomSparseStream for sparse sequence vectors

template <class Ring, class _Vector, class RandIter>
class RandomSparseStream<Ring, _Vector, RandIter, VectorCategories::SparseVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomSparseStream<Ring, Vector, RandIter, VectorCategories::SparseVectorTag> Self_t;

	RandomSparseStream (const Ring &F, double p, size_t n, size_t m = 0)
		: _F (F), _r1 (F), _r (F, _r1), _n (n), _m (m), _j (0),
		  _MT (time (NULL))
		{ setP (p); }

	RandomSparseStream (const Ring &F, const RandIter &r, double p, size_t n, size_t m = 0, int seed = 0)
		: _F (F), _r1 (r), _r (F, _r1), _n (n), _p (p), _m (m), _j (0),
		  _MT (seed)
		{ setP (p); }

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

	void setP (double p)
	{
		linbox_check ((p >= 0.0) && (p <= 1.0)); 
		_p = p;
		_1_log_1mp   = 1 / log (1 - _p);
	}

    private:
	const Ring                      &_F;
	RandIter                          _r1;
	NonzeroRandIter<Ring, RandIter>  _r;
	size_t                            _n;
	double                            _p;
	double                            _1_log_1mp;
	size_t                            _m;
	size_t                            _j;
	MersenneTwister                   _MT;
};

// Specialisation for dense zero-one vectors

template <class Ring, class _Vector, class RandIter>
class RandomDenseStream<Ring, _Vector, RandIter, VectorCategories::DenseZeroOneVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomDenseStream<Ring, Vector, RandIter, VectorCategories::DenseZeroOneVectorTag> Self_t;

	RandomDenseStream (const Ring &F, size_t n, size_t m = 0)
		: _F (F), _r (F), _n (n), _m (m), _j (0)
	{}

	RandomDenseStream (const Ring &F, const RandIter &r, size_t n, size_t m = 0)
		: _F (F), _r (r), _n (n), _m (m), _j (0)
	{}

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

    private:
	const Ring &_F;
	RandIter     _r;
	size_t       _n;
	size_t       _m;
	size_t       _j;
};

// Specialization of RandomSparseStream for sparse zero-one vectors

template <class Ring, class _Vector, class RandIter>
class RandomSparseStream<Ring, _Vector, RandIter, VectorCategories::SparseZeroOneVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomSparseStream<Ring, Vector, RandIter, VectorCategories::SparseZeroOneVectorTag> Self_t;

	RandomSparseStream (const Ring &F, double p, size_t n, size_t m = 0)
		: _F (F), _n (n), _m (m), _j (0), _MT (0)
		{ setP (p); }

	RandomSparseStream (const Ring &F, const RandIter &r, double p, size_t n, size_t m = 0)
		: _F (F), _n (n), _m (m), _j (0), _MT (0)
		{ setP (p); }

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

	void setP (double p)
	{
		linbox_check ((p >= 0.0) && (p <= 1.0)); 
		_p = p;
		_1_log_1mp   = 1 / log (1 - _p);
	}

    private:
	const Ring                      &_F;
	size_t                            _n;
	double                            _p;
	double                            _1_log_1mp;
	size_t                            _m;
	size_t                            _j;
	MersenneTwister                   _MT;
};

// Specialization of RandomSparseStream for hybrid zero-one vectors

template <class Ring, class _Vector, class RandIter>
class RandomSparseStream<Ring, _Vector, RandIter, VectorCategories::HybridZeroOneVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef RandomSparseStream<Ring, Vector, RandIter, VectorCategories::HybridZeroOneVectorTag> Self_t;

	RandomSparseStream (const Ring &F, double p, size_t n, size_t m = 0)
		: _F (F), _n (n), _m (m), _j (0), _MT (0)
		{ setP (p); }

	RandomSparseStream (const Ring &F, const RandIter &r, double p, size_t n, size_t m = 0)
		: _F (F), _n (n), _m (m), _j (0), _MT (0)
		{ setP (p); }

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

	void setP (double p)
	{
		linbox_check ((p >= 0.0) && (p <= 1.0)); 
		_p = p;
		_1_log_1mp   = 1 / log (1 - _p);
	}

    private:
	const Ring                      &_F;
	size_t                            _n;
	double                            _p;
	double                            _1_log_1mp;
	size_t                            _m;
	size_t                            _j;
	MersenneTwister                   _MT;
};

/** Stream for e_1,...,e_n
 * Generates the sequence e_1,...,e_n over a given ring
 * 
 * This class is generic with respect to the underlying vector
 * representation.
 */

template <class Ring, class _Vector, class Trait = typename VectorTraits<Ring, _Vector>::VectorCategory>
class StandardBasisStream : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef StandardBasisStream<Ring, Vector, Trait> Self_t;

	/** Constructor
	 * Construct a new stream with the given ring and vector size.
	 * @param F Ring over which to create vectors
	 * @param n Size of vectors
	 */
	StandardBasisStream (Ring &F, size_t n);

	/** Get next element
	 * @param v Vector into which to generate vector
	 * @return reference to new vector
	 */
	Vector &get (Vector &v);

	/** Extraction operator form
	 */
	Self_t &operator >> (Vector &v)
		{ get (v); return *this; }
	/** Number of vectors to be created
	 */
	size_t size () const;

	/** Number of vectors created so far
	 */
	size_t pos () const;

	/** Dimension of the space
	 */
	size_t dim () const;

	/** Check whether we have reached the end
	 */
	operator bool () const;

	/** Advance the stream k positions
	 */
	void advance (int k);

	/** Reset the stream to start at the beginning
	 */
	void reset ();

    private:
	const Ring              &_F;
	size_t                    _n;
	size_t                    _j;
};

// Specialization of standard basis stream for dense vectors

template <class Ring, class _Vector>
class StandardBasisStream<Ring, _Vector, VectorCategories::DenseVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef StandardBasisStream<Ring, Vector, VectorCategories::DenseVectorTag> Self_t;

	StandardBasisStream (const Ring &F, size_t n)
		: _F (F), _n (n), _j (0)
	{}

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _n; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _j < _n; }
	void advance (int k) { _j += k; }
	void reset () { _j = 0; }

    private:
	const Ring              &_F;
	size_t                    _n;
	size_t                    _j;
};

// Specialization of standard basis stream for sparse sequence vectors

template <class Ring, class _Vector>
class StandardBasisStream<Ring, _Vector, VectorCategories::SparseVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef StandardBasisStream<Ring, Vector, VectorCategories::SparseVectorTag> Self_t;

	StandardBasisStream (Ring &F, size_t n)
		: _F (F), _n (n), _j (0)
		{ _F.init (_one, 1); }

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _n; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _j < _n; }
	void advance (int k) { _j += k; }
	void reset () { _j = 0; }

    private:
	const Ring              &_F;
	size_t                    _n;
	size_t                    _j;
	typename Ring::Element   _one;
};

template <class Ring, class _Vector>
class StandardBasisStream<Ring, _Vector, VectorCategories::DenseZeroOneVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef StandardBasisStream<Ring, Vector, VectorCategories::DenseZeroOneVectorTag> Self_t;

	StandardBasisStream (const Ring &F, size_t n)
		: _F (F), _n (n), _j (0)
	{}

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _n; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _j < _n; }
	void advance (int k) { _j += k; }
	void reset () { _j = 0; }

    private:
	const Ring              &_F;
	size_t                    _n;
	size_t                    _j;
};

template <class Ring, class _Vector>
class StandardBasisStream<Ring, _Vector, VectorCategories::HybridZeroOneVectorTag> : public VectorStream<_Vector>
{
    public:
	typedef _Vector Vector;
        typedef StandardBasisStream<Ring, Vector, VectorCategories::HybridZeroOneVectorTag> Self_t;

	StandardBasisStream (const Ring &F, size_t n)
		: _F (F), _n (n), _j (0)
	{}

	Vector &get (Vector &v);
	Self_t &operator >> (Vector &v) { get (v); return *this; }
	size_t size () const { return _n; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _j < _n; }
	void advance (int k) { _j += k; }
	void reset () { _j = 0; }

    private:
	const Ring              &_F;
	size_t                    _n;
	size_t                    _j;
};

} // namespace LinBox

#include "linbox/vector/stream.tcc"

#endif // __LINBOX_vector_stream_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
