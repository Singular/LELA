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
#include "linbox/vector/sparse-subvector.h"
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
					VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size, VectorWrapper::CompareSparseEntries ());

	for (; i != i_end; ++i)
		out[i->first - src_idx + dest_idx] = i->second;
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

template <class Vector>
void Splicer::append_word (Vector &v, size_t index, typename Vector::word_type word) const
{
	typedef WordTraits<typename Vector::word_type> WT;

	if (v.empty () || v.back ().second != index >> WT::logof_size) {
		if (index & WT::pos_mask == 0)
			v.push_back (typename Vector::value_type (index >> WT::logof_size, word));
		else {
			v.push_back (typename Vector::value_type (index >> WT::logof_size,
								  Vector::Endianness::shift_right (word, index & WT::pos_mask)));
			v.push_back (typename Vector::value_type ((index >> WT::logof_size) + 1,
								  Vector::Endianness::shift_left (word, WT::bits - (index & WT::pos_mask))));
		}
	} else {
		v.back ().second |= Vector::Endianness::shift_right (word, index & WT::pos_mask);
		v.push_back (typename Vector::value_type ((index >> WT::logof_size) + 1,
							  Vector::Endianness::shift_left (word, WT::bits - (index & WT::pos_mask))));
	}
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag) const
{
	SparseSubvector<Vector2, VectorCategories::HybridZeroOneVectorTag> v1 (in, src_idx, src_idx + size);
	typename SparseSubvector<Vector2, VectorCategories::HybridZeroOneVectorTag>::iterator i_v1;

	for (i_v1 = v1.begin (); i_v1 != v1.end (); ++i_v1)
		append_word (out, dest_idx + (i_v1->first << WordTraits<typename Vector2::word_type>::logof_size), i_v1->second);
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag) const
{
	VectorDomain<Field> VD (F);

	SparseSubvector<Vector2, typename VectorTraits<Field, Vector2>::VectorCategory> v1 (in, src_idx, src_idx + size);
	BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	VD.copy (v2, v1);
}

template <class Field, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Field &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i = in.begin () + src_idx;
	typename Vector2::const_iterator i_end = in.begin () + (src_idx + size);
	size_t idx = dest_idx;

	for (; i != i_end; ++i, ++idx)
		if (*i)
			out.push_back (idx);
}

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
					VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag) const
{
	VectorDomain<Field> VD (F);

	SparseSubvector<Vector2, typename VectorTraits<Field, Vector2>::VectorCategory> v1 (in, src_idx, src_idx + size);
	BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	VD.copy (v2, v1);
}

template <class Field, class Matrix>
void Splicer::attach_identity (const Field &F, Matrix &A, const Block &horiz_block, const Block &vert_block, MatrixCategories::RowMatrixTag) const
{
	typename Matrix::RowIterator i_A = A.rowBegin () + horiz_block.destIndex ();
	typename Matrix::RowIterator i_A_end = A.rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
	size_t idx = horiz_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++idx)
		if (vert_block.isSourceIndexInBlock (idx))
			attach_e_i (F, *i_A, vert_block.sourceToDestIndex (idx));
}

template <class Field, class Matrix>
void Splicer::attach_identity (const Field &F, Matrix &A, const Block &horiz_block, const Block &vert_block, MatrixCategories::ColMatrixTag) const
{
	typename Matrix::ColIterator i_A = A.colBegin () + vert_block.destIndex ();
	typename Matrix::ColIterator i_A_end = A.colBegin () + (vert_block.destIndex () + vert_block.size ());
	size_t idx = vert_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++idx)
		if (horiz_block.isSourceIndexInBlock (idx))
			attach_e_i (F, *i_A, horiz_block.sourceToDestIndex (idx));
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::attach_source (const Field &F, Matrix1 &A, const Matrix2 &S, const Block &horiz_block, const Block &vert_block, MatrixCategories::RowMatrixTag) const
{
	typename Matrix1::RowIterator i_A = A.rowBegin () + horiz_block.destIndex ();
	typename Matrix1::RowIterator i_A_end = A.rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
	typename Matrix2::ConstRowIterator i_S = S.rowBegin () + horiz_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++i_S)
		attach_block (F, *i_A, *i_S, vert_block.sourceIndex (), vert_block.destIndex (), vert_block.size ());
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::attach_source (const Field &F, Matrix1 &A, const Matrix2 &S, const Block &horiz_block, const Block &vert_block, MatrixCategories::ColMatrixTag) const
{
	typename Matrix1::ColIterator i_A = A.colBegin () + vert_block.destIndex ();
	typename Matrix1::ColIterator i_A_end = A.colBegin () + (vert_block.destIndex () + vert_block.size ());
	typename Matrix2::ConstColIterator i_S = S.colBegin () + vert_block.sourceIndex ();

	for (; i_A != i_A_end; ++i_A, ++i_S)
		attach_block (F, *i_A, *i_S, horiz_block.sourceIndex (), horiz_block.destIndex (), horiz_block.size ());
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::splice (const Field &F, MatrixPart<Matrix1, MatrixCategories::RowMatrixTag> **input, MatrixPart<Matrix2, MatrixCategories::RowMatrixTag> **output) const
{
	linbox_check (!_horiz_blocks.empty ());
	linbox_check (!_vert_blocks.empty ());

	linbox_check (check ());

	commentator.start ("Splicing matrices", __FUNCTION__);

	typename std::vector<Block>::const_iterator horiz_block, vert_block;

	for (horiz_block = _horiz_blocks.begin (); horiz_block != _horiz_blocks.end (); ++horiz_block) {
		for (vert_block = _vert_blocks.begin (); vert_block != _vert_blocks.end (); ++vert_block) {
			if (output[horiz_block->dest ()][vert_block->dest ()].type () == MatrixPart<Matrix2>::TYPE_MATRIX) {
				switch (input[horiz_block->source ()][vert_block->source ()].type ()) {
				case MatrixPart<Matrix1>::TYPE_ZERO:
					break;

				case MatrixPart<Matrix1>::TYPE_IDENTITY:
					attach_identity (F, output[horiz_block->dest ()][vert_block->dest ()].A (),
							 *horiz_block, *vert_block, MatrixCategories::RowMatrixTag ());
					break;

				case MatrixPart<Matrix1>::TYPE_MATRIX:
					attach_source (F, output[horiz_block->dest ()][vert_block->dest ()].A (),
						       input[horiz_block->source ()][vert_block->source ()].A (),
						       *horiz_block, *vert_block, MatrixCategories::RowMatrixTag ());
					break;
				}
			}
		}
	}

	commentator.stop (MSG_DONE);
}

template <class Field, class Matrix1, class Matrix2>
void Splicer::splice (const Field &F, MatrixPart<Matrix1, MatrixCategories::ColMatrixTag> **input, MatrixPart<Matrix2, MatrixCategories::ColMatrixTag> **output) const
{
	linbox_check (!_horiz_blocks.empty ());
	linbox_check (!_vert_blocks.empty ());

	linbox_check (check ());

	commentator.start ("Splicing matrices", __FUNCTION__);

	typename std::vector<Block>::const_iterator horiz_block, vert_block;

	for (horiz_block = _horiz_blocks.begin (); horiz_block != _horiz_blocks.end (); ++horiz_block) {
		for (vert_block = _vert_blocks.begin (); vert_block != _vert_blocks.end (); ++vert_block) {
			if (output[horiz_block->dest ()][vert_block->dest ()].type () == MatrixPart<Matrix2>::TYPE_MATRIX) {
				switch (input[horiz_block->source ()][vert_block->source ()].type ()) {
				case MatrixPart<Matrix1>::TYPE_ZERO:
					break;

				case MatrixPart<Matrix1>::TYPE_IDENTITY:
					attach_identity (F, output[horiz_block->dest ()][vert_block->dest ()].A (),
							 *horiz_block, *vert_block, MatrixCategories::ColMatrixTag ());
					break;

				case MatrixPart<Matrix1>::TYPE_MATRIX:
					attach_source (F, output[horiz_block->dest ()][vert_block->dest ()].A (),
						       input[horiz_block->source ()][vert_block->source ()].A (),
						       *horiz_block, *vert_block, MatrixCategories::ColMatrixTag ());
					break;
				}
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
