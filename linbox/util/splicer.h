/* linbox/util/splicer.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#ifndef __LINBOX_UTIL_SPLICER_H
#define __LINBOX_UTIL_SPLICER_H

#include <iostream>

#include "linbox/vector/vector-traits.h"

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

/** Class representing a matrix to be spliced by the Splicer */
template<class Matrix>
class MatrixPart {
public:
	enum Type { TYPE_ZERO, TYPE_IDENTITY, TYPE_MATRIX };

private:
	Matrix *_A;
	const Type _type;

public:
	/** Construct a MatrixPart-object from a matrix
	 *
	 * @param A Source-matrix
	 */
	MatrixPart (Matrix &A) : _A (&A), _type (TYPE_MATRIX) {}

	/** Construct a MatrixPart-object with no matrix
	 *
	 * @param type Type of source. Must be TYPE_ZERO or
	 * TYPE_IDENTITY.
	 */
	MatrixPart (Type type) : _A (NULL), _type (type) {}

	/** Copy-constructor */
	MatrixPart (const MatrixPart &S) : _A (S._A), _type (S._type) {}

	/** Assignment-operator */
	MatrixPart &operator = (const MatrixPart &S) { _A = S._A; _type = S._type; }

	Type type () const { return _type; }
	Matrix &A () { return *_A; }
};

/** Class to chop a matrix up into smaller matrices which may be
 * distributed across non-continuous rows and columns and to assemble
 * smaller matrices into a large matrix in the same manner.
 *
 * Horizontal respectively vertical blocks must be arranged from left
 * to right respectively top to bottom in both the source- and
 * destination-matrix with no gaps in either. The method check ()
 * checks whether this is so and reports problems.
 */
class Splicer {
	std::vector<Block> _vert_blocks, _horiz_blocks;

	std::vector<Block> &map_blocks (std::vector<Block> &output,
					const std::vector<Block> &outer_blocks,
					const std::vector<Block> &inner_blocks,
					unsigned int inner_source,
					unsigned int only_source) const;

	template <class Field, class Vector>
	void attach_e_i_specialised (const Field &F, Vector &row, size_t idx, VectorCategories::DenseVectorTag) const
		{ row[idx] = F.one (); }

	template <class Field, class Vector>
	void attach_e_i_specialised (const Field &F, Vector &row, size_t idx, VectorCategories::SparseVectorTag) const
		{ row.push_back (typename Vector::value_type (idx, F.one ())); }

	template <class Field, class Vector>
	void attach_e_i_specialised (const Field &F, Vector &row, size_t idx, VectorCategories::DenseZeroOneVectorTag) const
		{ row[idx] = F.one (); }

	template <class Field, class Vector>
	void attach_e_i_specialised (const Field &F, Vector &row, size_t idx, VectorCategories::SparseZeroOneVectorTag) const
		{ row.push_back (typename Vector::value_type (idx)); }

	template <class Field, class Vector>
	void attach_e_i_specialised (const Field &F, Vector &row, size_t idx, VectorCategories::HybridZeroOneVectorTag) const
		{ row.push_back (typename Vector::value_type
				 (idx >> WordTraits<typename Vector::word_type>::logof_size,
				  Vector::Endianness::e_j (idx & WordTraits<typename Vector::word_type>::pos_mask))); }

	template <class Field, class Vector>
	void attach_e_i (const Field &F, Vector &row, size_t idx) const
		{ attach_e_i_specialised (F, row, idx, typename VectorTraits<Field, Vector>::VectorCategory ()); }

	size_t round_up (size_t v, size_t x) const
		{ return ((v + x - 1) / x) * x; }

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag) const
		{ std::copy (in.begin () + src_idx, in.begin () + (src_idx + size), out.begin () + dest_idx); }

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
				       VectorCategories::GenericVectorTag, VectorCategories::GenericVectorTag) const;

	template <class Field, class Vector1, class Vector2>
	void attach_block (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size) const
		{ attach_block_specialised (F, out, in, src_idx, dest_idx, size,
					    typename VectorTraits<Field, Vector1>::VectorCategory (),
					    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

	template <class Field, class Matrix>
	void attach_identity (const Field &F, Matrix &A, const Block &horiz_block, const Block &vert_block) const;

	template <class Field, class Matrix1, class Matrix2>
	void attach_source (const Field &F, Matrix1 &A, const Matrix2 &S, const Block &horiz_block, const Block &vert_block) const;

	bool check_blocks (const std::vector<Block> &blocks, const char *type) const;

	void consolidate_blocks (std::vector<Block> &blocks);

	void remove_gaps_from_blocks (std::vector<Block> &blocks);

	void reverse_blocks (std::vector<Block> &out_blocks, const std::vector<Block> &in_blocks) const;

	friend std::ostream &operator << (std::ostream &os, const Splicer &splicer);

public:
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
	void fillHorizontal (unsigned int sid, unsigned int did, size_t rowdim);

	/** Add a vertical block with the given source so that the
	 * total column-dimension matches what is requested.
	 *
	 * @param sid Source-id of block to be added
	 * @param did Destination-id of block to be added
	 * @param rowdim Desired column-dimension
	 */
	void fillVertical (unsigned int sid, unsigned int did, size_t coldim);

	/** Find successive blocks with an identical source-id and
	 * consolidate them into single blocks
	 */
	void consolidate ();

	/** Find and remove gaps in the blocks, updating
	 * destination-ids appropriately.
	 */
	void removeGaps ();

	/** Splice the given set of matrices into a single matrix
	 * based on the given block-decomposition
	 *
	 * @param F Field over which matrices are defined. Needed to
	 * build identity-matrix when requested and to copy vectors.
	 *
	 * @param input Array of arrays of MatrixPart-objects. The pointer
	 * at index i,j should point to the source-matrix for blocks
	 * with that source-index. If the source-type for a block is
	 * SOURCE_MATRIX, then the pointer at the corresponding index
	 * in this vector must be non-null.
	 *
	 * @param output Array of arrays of MatrixPart-objects. Matrices
	 * should already be allocated and be of the right sizes. The
	 * pointer at index i,j in the vector is whither data from
	 * blocks with destination-index i,j are to be copied. If a
	 * pointer at a given index is null, then the matrix is
	 * ignored and the data are not copied.
	 */
	template <class Field, class Matrix1, class Matrix2>
	void splice (const Field &F, MatrixPart<Matrix1> **input, MatrixPart<Matrix2> **output) const;

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
