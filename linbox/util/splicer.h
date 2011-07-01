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
 * the source respectively destination is used. This may be read from
 * the corresponding horizontal and vertical @ref Block objects which
 * are passed to the closure grid when the method splice runs.
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
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorRepresentationTypes::Dense)
		{ row[idx] = F.one (); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorRepresentationTypes::Sparse)
		{ row.push_back (typename Vector::value_type (idx, F.one ())); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorRepresentationTypes::Dense01)
		{ row[idx] = F.one (); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorRepresentationTypes::Sparse01)
		{ row.push_back (typename Vector::value_type (idx)); }

	template <class Ring, class Vector>
	static void attach_e_i_specialised (const Ring &F, Vector &row, size_t idx, VectorRepresentationTypes::Hybrid01)
		{ row.push_back (typename Vector::value_type
				 (idx >> WordTraits<typename Vector::word_type>::logof_size,
				  Vector::Endianness::e_j (idx & WordTraits<typename Vector::word_type>::pos_mask))); }

	static size_t round_up (size_t v, size_t x)
		{ return ((v + x - 1) / x) * x; }

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
		{ std::copy (in.begin () + src_idx, in.begin () + (src_idx + size), out.begin () + dest_idx); }

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01);

	template <class Vector>
	static void append_word (Vector &v, size_t index, typename Vector::word_type word);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01);

	template <class Ring, class Vector1, class Vector2>
	static void attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					      VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01);

	/// Functions for attaching pieces to vectors
	template <class Ring, class Vector>
	static void attach_e_i (const Ring &F, Vector &row, size_t idx)
		{ attach_e_i_specialised (F, row, idx, typename VectorTraits<Ring, Vector>::RepresentationType ()); }

	template <class Ring, class Vector1, class Vector2>
	static void attach_block (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size)
		{ attach_block_specialised (F, out, in, src_idx, dest_idx, size,
					    typename VectorTraits<Ring, Vector1>::RepresentationType (),
					    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }

	bool check_blocks (const std::vector<Block> &blocks, const char *type) const;
	void consolidate_blocks (std::vector<Block> &blocks);
	void remove_gaps_from_blocks (std::vector<Block> &blocks);
	void reverse_blocks (std::vector<Block> &out_blocks, const std::vector<Block> &in_blocks) const;
	void fill_blocks (std::vector<Block> &blocks, unsigned int sid, unsigned int did, size_t dim);

	friend std::ostream &operator << (std::ostream &os, const Splicer &splicer);

	template <class Ring, class Matrix1, class Matrix2>
	static void copyBlockSpecialised (const Ring &R, const Matrix1 &source, Matrix2 &dest, const Block &horiz_block, const Block &vert_block,
					  MatrixIteratorTypes::Row, MatrixIteratorTypes::Row);

	template <class Ring, class Matrix1, class Matrix2>
	static void copyBlockSpecialised (const Ring &R, const Matrix1 &source, Matrix2 &dest, const Block &horiz_block, const Block &vert_block,
					  MatrixIteratorTypes::Col, MatrixIteratorTypes::Col);

	template <class Ring, class Matrix1, class Matrix2>
	static void copyBlockSpecialised (const Ring &R, const Matrix1 &source, Matrix2 &dest, const Block &horiz_block, const Block &vert_block,
					  MatrixIteratorTypes::RowCol, MatrixIteratorTypes::RowCol)
		{ copyBlockSpecialised (R, source, dest, horiz_block, vert_block, MatrixIteratorTypes::Row (), MatrixIteratorTypes::Row ()); }

	template <class Ring, class Matrix>
	static void copyIdentitySpecialised (const Ring &R, Matrix &dest, const Block &horiz_block, const Block &vert_block,
					     MatrixIteratorTypes::Row);

	template <class Ring, class Matrix>
	static void copyIdentitySpecialised (const Ring &R, Matrix &dest, const Block &horiz_block, const Block &vert_block,
					     MatrixIteratorTypes::Col);

	template <class Ring, class Matrix>
	static void copyIdentitySpecialised (const Ring &R, Matrix &dest, const Block &horiz_block, const Block &vert_block,
					     MatrixIteratorTypes::RowCol)
		{ copyIdentitySpecialised (R, dest, horiz_block, vert_block, MatrixIteratorTypes::Row ()); }

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
	 * @param grid Closure which performs copying-operations. It
	 * should have the signature grid (horizontal_block,
	 * vertical_block) (horizontal_block and vertical_block
	 * instances of type @ref Block) and should perform the actual
	 * copying according to the information in those blocks. The
	 * static methods @ref copyBlock and @ref copyIdentity are
	 * available to assist.
	 */
	template <class Grid>
	void splice (Grid grid) const;

	/** Copy from the source-matrix to the destination-matrix
	 * according to the information in horiz_block and
	 * vert_block. To be used inside a Grid-object in @ref splice
	 *
	 * @param R Ring over which to operate
	 * @param source Source-matrix
	 * @param dest Destination-matrix
	 * @param horiz_block Horizontal block
	 * @param vert_block Vertical block
	 */
	template <class Ring, class Matrix1, class Matrix2>
	static void copyBlock (const Ring &R, const Matrix1 &source, Matrix2 &dest, const Block &horiz_block, const Block &vert_block)
		{ copyBlockSpecialised (R, source, dest, horiz_block, vert_block,
					typename Matrix1::IteratorType (),
					typename Matrix2::IteratorType ()); }

	/** Copy a part of the identity-matrix to the
	 * destination-matrix according to the information in
	 * horiz_block and vert_block. To be used inside a Grid-object
	 * in @ref splice
	 *
	 * @param R Ring over which to operate
	 * @param dest Destination-matrix
	 * @param horiz_block Horizontal block
	 * @param vert_block Vertical block
	 */
	template <class Ring, class Matrix>
	static void copyIdentity (const Ring &R, Matrix &dest, const Block &horiz_block, const Block &vert_block)
		{ copyIdentitySpecialised (R, dest, horiz_block, vert_block,
					   typename Matrix::IteratorType ()); }
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
