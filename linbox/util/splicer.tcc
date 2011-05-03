/* linbox/util/splicer.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#ifndef __LINBOX_UTIL_SPLICER_TCC
#define __LINBOX_UTIL_SPLICER_TCC

#include "linbox/util/splicer.h"
#include "linbox/util/commentator.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector.h"
#include "linbox/matrix/submatrix.h"
#include "linbox/matrix/matrix-domain.h"

namespace LinBox
{

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size, VectorWrapper::CompareSparseEntries ());

	for (; i != i_end; ++i)
		out.push_back (typename Vector1::value_type (i->first - src_idx + dest_idx, i->second));
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag) const
{
	typename Vector2::const_iterator i = in.begin () + src_idx;
	typename Vector2::const_iterator i_end = in.begin () + (src_idx + size);

	for (; i != i_end; ++i, ++dest_idx)
		out.push_back (typename Vector1::value_type (dest_idx, *i->second));
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag) const
{
	VectorDomain<Field> VD (F);

	typedef BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> DestSubvector;
	typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> SourceSubvector;

	SourceSubvector v1 (in.begin () + src_idx, in.begin () + (src_idx + size));
	DestSubvector v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	VD.copy (v2, v1);
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx);
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size);

	for (; i != i_end; ++i)
		out.push_back (typename Vector1::value_type (*i - src_idx + dest_idx));
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag) const
	{}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag) const
{
	VectorDomain<Field> VD (F);

	typedef BitVector<typename Vector1::Endianness> TmpVector;
	typedef BitSubvector<typename TmpVector::iterator, typename TmpVector::const_iterator> TmpSubvector;
	typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> ConstTmpSubvector;

	size_t dest_idx_offset = dest_idx & WordTraits<typename Vector2::word_type>::pos_mask;

	TmpVector tmp (round_up (size + dest_idx_offset, WordTraits<typename Vector2::word_type>::bits));
	ConstTmpSubvector v1 (in.begin () + src_idx, in.begin () + (src_idx + size));
	TmpSubvector v2 (tmp.begin () + dest_idx_offset, tmp.begin () + (dest_idx_offset + size));

	VD.copy (v2, v1);

	typename TmpVector::word_iterator i = tmp.wordBegin ();
	size_t idx = dest_idx >> WordTraits<typename Vector1::word_type>::logof_size;

	if (!out.empty () && out.back ().first == idx) {
		out.back ().second |= *i;
		++i;
		++idx;
	}

	for (; i != tmp.wordEnd (); ++i, ++idx)
		if (*i)
			out.push_back (typename Vector1::value_type (idx, *i));
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::GenericVectorTag, VectorCategories::GenericVectorTag) const
	{}

template <class Field, class Matrix>
void Splicer::attach_identity (const Field &F, Matrix &A, const Block &horiz_block, const Block &vert_block) const
{
	typename Matrix::RowIterator i_A = A.rowBegin () + horiz_block.destIndex ();
	typename Matrix::RowIterator i_A_end = A.rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
	size_t idx = horiz_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++idx)
		if (vert_block.isSourceIndexInBlock (idx))
			attach_e_i (F, *i_A, vert_block.sourceToDestIndex (idx));
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::attach_source (const Field &F, Matrix1 &A, const Matrix2 &S, const Block &horiz_block, const Block &vert_block) const
{
	typename Matrix1::RowIterator i_A = A.rowBegin () + horiz_block.destIndex ();
	typename Matrix1::RowIterator i_A_end = A.rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
	typename Matrix2::ConstRowIterator i_S = S.rowBegin () + horiz_block.sourceIndex ();
	typename Matrix2::ConstRowIterator i_S_end = S.rowBegin () + (horiz_block.sourceIndex () + horiz_block.size ());

	for (; i_A != i_A_end; ++i_A, ++i_S)
		attach_block (F, *i_A, *i_S, vert_block.sourceIndex (), vert_block.destIndex (), vert_block.size ());
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::chop (const Field &F, SourceMatrix<Matrix1> **output, const Matrix2 &A) const
{
	linbox_check (!_horiz_blocks.empty ());
	linbox_check (!_vert_blocks.empty ());

	linbox_check (check ());

	linbox_check (_horiz_blocks.back ().destIndex () + _horiz_blocks.back ().size () == A.rowdim ());
	linbox_check (_vert_blocks.back ().destIndex () + _vert_blocks.back ().size () == A.coldim ());

	commentator.start ("Chopping matrix", __FUNCTION__);

	MatrixDomain<Field> MD (F);

	typename std::vector<Block>::const_iterator horiz_block, vert_block;

	for (horiz_block = _horiz_blocks.begin (); horiz_block != _horiz_blocks.end (); ++horiz_block) {
		for (vert_block = _vert_blocks.begin (); vert_block != _vert_blocks.end (); ++vert_block) {
			if (output[horiz_block->source ()][vert_block->source ()].type () == SourceMatrix<Matrix1>::TYPE_MATRIX) {
				Submatrix<const Matrix2> dest_part (A, horiz_block->destIndex (), vert_block->destIndex (), horiz_block->size (), vert_block->size ());
				Submatrix<Matrix1> source_part (output[horiz_block->source ()][vert_block->source ()].A (),
								horiz_block->sourceIndex (), vert_block->sourceIndex (), horiz_block->size (), vert_block->size ());

				MD.copy (source_part, dest_part);
			}
		}
	}

	commentator.stop (MSG_DONE, NULL, __FUNCTION__);
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::splice (const Field &F, Matrix1 &A, SourceMatrix<Matrix2> **input) const
{
	linbox_check (!_horiz_blocks.empty ());
	linbox_check (!_vert_blocks.empty ());

	linbox_check (_horiz_blocks.back ().destIndex () + _horiz_blocks.back ().size () == A.rowdim ());
	linbox_check (_vert_blocks.back ().destIndex () + _vert_blocks.back ().size () == A.coldim ());

	linbox_check (check ());

	commentator.start ("Splicing matrices", __FUNCTION__);

	typename std::vector<Block>::const_iterator horiz_block, vert_block;

	MatrixDomain<Field> MD (F);

	MD.scal (A, F.zero ());

	for (horiz_block = _horiz_blocks.begin (); horiz_block != _horiz_blocks.end (); ++horiz_block) {
		for (vert_block = _vert_blocks.begin (); vert_block != _vert_blocks.end (); ++vert_block) {
			switch (input[horiz_block->source ()][vert_block->source ()].type ()) {
			case SourceMatrix<Matrix2>::TYPE_ZERO:
				break;

			case SourceMatrix<Matrix2>::TYPE_IDENTITY:
				attach_identity (F, A, *horiz_block, *vert_block);
				break;

			case SourceMatrix<Matrix2>::TYPE_MATRIX:
				attach_source (F, A, input[horiz_block->source ()][vert_block->source ()].A (), *horiz_block, *vert_block);
				break;
			}
		}
	}

	commentator.stop (MSG_DONE);
}

} // namespace LinBox

#endif // __LINBOX_UTIL_SPLICER_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
