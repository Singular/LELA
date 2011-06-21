/* linbox/util/splicer.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#ifndef __LINBOX_UTIL_SPLICER_H
#define __LINBOX_UTIL_SPLICER_H

#include <iostream>

#include "linbox/vector/traits.h"
#include "linbox/matrix/traits.h"

namespace LinBox
{

/** Class representing a horizontal or vertical block for the splicer. */
class Block {
private:

	unsigned int _source;
	unsigned int _dest;
	size_t _src_idx;
	size_t _dest_idx;
	size_t _size;

public:
	/** Construct a new block
	 *
	 * @param source Index of the source-matrix
	 * @param src_idx Row- resp. column-index of block in source-matrix
	 * @param dest_idx Row- resp. column-index of block in destination-matrix
	 * @param size Size of block
	 */
	Block (unsigned int source, unsigned int dest, size_t src_idx, size_t dest_idx, size_t size)
		: _source (source), _dest (dest), _src_idx (src_idx), _dest_idx (dest_idx), _size (size) {}

	/** Get index of source-matrix */
	inline unsigned int source () const { return _source; }

	/** Get index of destination-matrix */
	inline unsigned int dest () const { return _dest; }

	/** Get column- resp. row-index of block in source-matrix */
	inline size_t sourceIndex () const { return _src_idx; }

	/** Get column- resp. row-index of block in destination-matrix */
	inline size_t destIndex () const { return _dest_idx; }

	/** Get column- resp. row-index of block after this one in destination-matrix */
	inline size_t destIndexNextBlock () const { return _dest_idx + _size; }

	/** Map a column- resp. row-index in the source-matrix to one in the destination-matrix */
	inline size_t sourceToDestIndex (size_t idx) const { return idx - _src_idx + _dest_idx; }

	/** Map a column- resp. row-index in the destination-matrix to one in the source-matrix */
	inline size_t destToSourceIndex (size_t idx) const { return idx - _dest_idx + _src_idx; }

	/** Get the size of the block */
	inline size_t size () const { return _size; }

	/** Determine whether the given column- resp. row-index in the source-matrix is in this block */
	inline bool isSourceIndexInBlock (size_t idx) const { return idx >= _src_idx && idx < _src_idx + _size; }

	/** Determine whether the given column- resp. row-index in the destination-matrix is in this block */
	inline bool isDestIndexInBlock (size_t idx) const { return idx >= _dest_idx && idx < _dest_idx + _size; }
};

/** Class representing a part of a matrix to be spliced by the Splicer */
template <class Ring>
class MatrixPart
{
public:
	/// Copy from another matrix to this
	virtual void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart *source) = 0;
};

/** Zero-matrix */
template <class Ring>
class MatrixPartZero : public MatrixPart<Ring>
{
public:
	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source) {}
};

/** Identity-matrix */
template <class Ring, class Matrix, class Trait = typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory>
class MatrixPartIdentity : public MatrixPart<Ring>
{
public:
	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source) {}
};

template <class Ring, class Matrix>
class MatrixPartIdentity<Ring, Matrix, MatrixCategories::RowMatrixTag> : public MatrixPart<Ring>
{
public:
	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source) {}
};

template <class Ring, class Matrix>
class MatrixPartIdentity<Ring, Matrix, MatrixCategories::ColMatrixTag> : public MatrixPart<Ring>
{
public:
	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source) {}
};

template <class Ring, class Matrix>
class MatrixPartIdentity<Ring, Matrix, MatrixCategories::RowColMatrixTag> : public MatrixPartIdentity<Ring, Matrix, MatrixCategories::RowMatrixTag> {};

/** Source-matrix */
template <class Ring, class Matrix, class OtherMatrix, class Trait = typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory>
class MatrixPartMatrix : public MatrixPart<Ring>
{
	friend class MatrixPartMatrix<Ring, OtherMatrix, Matrix, Trait>;

protected:
	Matrix *_A;

public:
	/** Construct a MatrixPart-object from a matrix
	 *
	 * @param A Source-matrix
	 */
	MatrixPartMatrix (Matrix &A) : _A (&A) {}

	/** Copy-constructor */
	MatrixPartMatrix (const MatrixPartMatrix &S) : _A (S._A) {}

	/** Assignment-operator */
	MatrixPartMatrix &operator = (const MatrixPartMatrix &S) { _A = S._A; return *this; }

	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source);
};

template <class Ring, class Matrix, class OtherMatrix>
class MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::RowMatrixTag> : public MatrixPart<Ring>
{
	friend class MatrixPartMatrix<Ring, OtherMatrix, Matrix, MatrixCategories::RowMatrixTag>;

protected:
	Matrix *_A;

public:
	MatrixPartMatrix (Matrix &A) : _A (&A) {}
	MatrixPartMatrix (const MatrixPartMatrix &S) : _A (S._A) {}
	MatrixPartMatrix &operator = (const MatrixPartMatrix &S) { _A = S._A; return *this; }

	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source);
};

template <class Ring, class Matrix, class OtherMatrix>
class MatrixPartMatrix<Ring, const Matrix, OtherMatrix, MatrixCategories::RowMatrixTag> : public MatrixPart<Ring>
{
	friend class MatrixPartMatrix<Ring, OtherMatrix, const Matrix, MatrixCategories::RowMatrixTag>;

protected:
	const Matrix *_A;

public:
	MatrixPartMatrix (const Matrix &A) : _A (&A) {}
	MatrixPartMatrix (const MatrixPartMatrix &S) : _A (S._A) {}
	MatrixPartMatrix &operator = (const MatrixPartMatrix &S) { _A = S._A; return *this; }

	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source) {}
};

template <class Ring, class Matrix, class OtherMatrix>
class MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::ColMatrixTag> : public MatrixPart<Ring>
{
	friend class MatrixPartMatrix<Ring, OtherMatrix, Matrix, MatrixCategories::ColMatrixTag>;

protected:
	Matrix *_A;

public:
	MatrixPartMatrix (Matrix &A) : _A (&A) {}
	MatrixPartMatrix (const MatrixPartMatrix &S) : _A (S._A) {}
	MatrixPartMatrix &operator = (const MatrixPartMatrix &S) { _A = S._A; return *this; }

	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source);
};

template <class Ring, class Matrix, class OtherMatrix>
class MatrixPartMatrix<Ring, const Matrix, OtherMatrix, MatrixCategories::ColMatrixTag> : public MatrixPart<Ring>
{
	friend class MatrixPartMatrix<Ring, OtherMatrix, const Matrix, MatrixCategories::ColMatrixTag>;

protected:
	const Matrix *_A;

public:
	MatrixPartMatrix (const Matrix &A) : _A (&A) {}
	MatrixPartMatrix (const MatrixPartMatrix &S) : _A (S._A) {}
	MatrixPartMatrix &operator = (const MatrixPartMatrix &S) { _A = S._A; return *this; }

	void copy (const Ring &R, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source) {}
};

template <class Ring, class Matrix, class OtherMatrix>
class MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::RowColMatrixTag> : public MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::RowMatrixTag>
{
	typedef MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::RowMatrixTag> parent_type;

public:
	MatrixPartMatrix (Matrix &A) : parent_type (A) {}
	MatrixPartMatrix (const MatrixPartMatrix &S) : parent_type (S) {}
	MatrixPartMatrix &operator = (const MatrixPartMatrix &S) { parent_type::operator = (S); return *this; }
};

/** Class to chop matrices up and splice them together.
 *
 * This class maintains a set of data-structures which describe how a
 * large virtual matrix is divided in two ways (a source and a
 * destination) into smaller matrices. It then provides a facility,
 * @ref splice, which, given appropriate sets of matrices (actually
 * arrays of type @ref MatrixPart), converts from the
 * source-representation to the destination-representation.
 *
 * The virtual matrix is divided horizontally and vertically into
 * blocks represented by the class @ref Block. Each block indicates a
 * source-id, a destination-id, a source-index, a destination-index,
 * and a size. A pair consisting of a horizontal block and a vertical
 * block represents an atomic division of the virtual matrix and
 * contains the information required to transfer the data from the
 * source-representation to the destination-representation.
 *
 * The pair (horizontal id, vertical id) identifies which matrix from
 * the source respectively destination is used. It gives the indices
 * into two-dimensional arrays of type @ref MatrixPart which are
 * passed to @ref splice. The indices are then of the row respectively
 * column where the block begins.
 *
 * Horizontal and vertical blocks must be arranged from left to right
 * respectively top to bottom in both the source- and
 * destination-matrix with no gaps in either and no overlaps. The
 * method @ref check checks whether this is so and reports
 * problems. The method @ref removeGaps searches for gaps in the
 * horizontal and vertical blocks and removes them.
 */
class Splicer {
	std::vector<Block> _vert_blocks, _horiz_blocks;

	std::vector<Block> &map_blocks (std::vector<Block> &output,
					const std::vector<Block> &outer_blocks,
					const std::vector<Block> &inner_blocks,
					unsigned int inner_source,
					unsigned int only_source) const;

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorCategories::DenseVectorTag)
		{ row[idx] = F.one (); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorCategories::SparseVectorTag)
		{ row.push_back (typename Vector::value_type (idx, F.one ())); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorCategories::DenseZeroOneVectorTag)
		{ row[idx] = F.one (); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorCategories::SparseZeroOneVectorTag)
		{ row.push_back (typename Vector::value_type (idx)); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorCategories::HybridZeroOneVectorTag)
		{ row.push_back (typename Vector::value_type
				 (idx >> WordTraits<typename Vector::word_type>::logof_size,
				  Vector::Endianness::e_j (idx & WordTraits<typename Vector::word_type>::pos_mask))); }

	static size_t round_up (size_t v, size_t x)
		{ return ((v + x - 1) / x) * x; }

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
		{ std::copy (in.begin () + src_idx, in.begin () + (src_idx + size), out.begin () + dest_idx); }

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

	template <class Vector>
	static void append_word (Vector &v, size_t index, typename Vector::word_type word);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

	bool check_blocks (const std::vector<Block> &blocks, const char *type) const;
	void consolidate_blocks (std::vector<Block> &blocks);
	void remove_gaps_from_blocks (std::vector<Block> &blocks);
	void reverse_blocks (std::vector<Block> &out_blocks, const std::vector<Block> &in_blocks) const;
	void fill_blocks (std::vector<Block> &blocks, unsigned int sid, unsigned int did, size_t dim);

	friend std::ostream &operator << (std::ostream &os, const Splicer &splicer);

public:
	/// @name Management of the blocks
	//@{

	/** Add a new block which divides the matrix horizontally */
	void addHorizontalBlock (const Block &block)
		{ _horiz_blocks.push_back (block); }

	/** Add a new block which divides the matrix vertically */
	void addVerticalBlock (const Block &block)
		{ _vert_blocks.push_back (block); }

	/** Clear all horizontal blocks */
	void clearHorizontalBlocks ()
		{ _horiz_blocks.clear (); }

	/** Clear all vertical blocks */
	void clearVerticalBlocks ()
		{ _vert_blocks.clear (); }

	/** Add a horizontal block with the given source so that the
	 * total row-dimension matches what is requested.
	 *
	 * @param sid Source-id of block to be added
	 * @param did Destination-id of block to be added
	 * @param rowdim Desired row-dimension
	 */
	void fillHorizontal (unsigned int sid, unsigned int did, size_t rowdim)
		{ fill_blocks (_horiz_blocks, sid, did, rowdim); }

	/** Add a vertical block with the given source so that the
	 * total column-dimension matches what is requested.
	 *
	 * @param sid Source-id of block to be added
	 * @param did Destination-id of block to be added
	 * @param rowdim Desired column-dimension
	 */
	void fillVertical (unsigned int sid, unsigned int did, size_t coldim)
		{ fill_blocks (_vert_blocks, sid, did, coldim); }

	/** Find successive blocks with identical source- and
	 * destination-ids and consolidate them into single blocks
	 */
	void consolidate ();

	/** Find and remove gaps in the blocks, updating indices
	 * appropriately.
	 */
	void removeGaps ();

	/** Compose this Splicer with another and produce a new Splicer
	 *
	 * To make a direct copy of either horizontal or vertical
	 * blocks, just pass a value for horiz_inner_source
	 * resp. vert_inner_source which doesn't correspond to any
	 * source in this splicer.
	 *
	 * @param output Output-splicer
	 * @param inner Splicer with which to compose this one
	 * @param horiz_inner_source Horizontal index of source to which inner blocks correspond
	 * @param vert_inner_source Vertical index of source to which inner blocks correspond
	 * @param horiz_only_source If not -1, only include horizontal blocks with this source-id and ignore all others
	 * @param vert_only_source If not -1, only include vertical blocks with this source-id and ignore all others
	 * @returns Reference to output
	 */
	Splicer &compose (Splicer &output,
			  const Splicer &inner,
			  unsigned int horiz_inner_source,
			  unsigned int vert_inner_source,
			  unsigned int horiz_only_source = -1,
			  unsigned int vert_only_source = -1) const;

	/** Construct the reverse of this splicer and place it in output
	 *
	 * @param output Splicer into which to place the reverse
	 * @returns Reference to output
	 */
	Splicer &reverse (Splicer &output) const;

	/** Check that the blocks are valid. Useful for debugging.
	 *
	 * If an error is found, prints a report using the commentator.
	 *
	 * @returns true if everything is okay, false if error found
	 */
	bool check () const;

	//@} Management of the blocks

	/** Splice the given set of matrices input into the set of
	 * matrices output
	 *
	 * @param F Ring over which matrices are defined. Needed to
	 * build identity-matrix when requested and to copy vectors.
	 *
	 * @param input Array of arrays of @ref MatrixPart objects
	 * from which data are copied.
	 *
	 * @param output Array of arrays of @ref MatrixPart objects to
	 * which data are copied. If the type of an output-matrix is
	 * MatrixTypeSource, then its corresponding matrix should
	 * already be allocated to the correct size and set to
	 * zero. Otherwise it is ignored by this operation.
	 */
	template <class Ring>
	void splice (const Ring &F, MatrixPart<Ring> ***input, MatrixPart<Ring> ***output) const;

	/// Functions for attaching pieces to vectors
	template <class Ring, class Vector>
	static void attach_e_i (const Ring &F, Vector &row, size_t idx)
		{ attach_e_i_specialised (F, row, idx, typename VectorTraits<Ring, Vector>::VectorCategory ()); }

	template <class Ring, class Vector1, class Vector2>
	static void attach_block (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size)
		{ attach_block_specialised (F, out, in, src_idx, dest_idx, size,
					    typename VectorTraits<Ring, Vector1>::VectorCategory (),
					    typename VectorTraits<Ring, Vector2>::VectorCategory ()); }
};

} // namespace LinBox

#endif // __LINBOX_UTIL_SPLICER_H

#include "linbox/util/splicer.tcc"

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
