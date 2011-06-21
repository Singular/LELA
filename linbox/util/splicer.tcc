/* linbox/util/splicer.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#ifndef __LINBOX_UTIL_SPLICER_TCC
#define __LINBOX_UTIL_SPLICER_TCC

#include "linbox/util/splicer.h"
#include "linbox/util/commentator.h"
#include "linbox/util/debug.h"
#include "linbox/blas/context.h"
#include "linbox/blas/level1.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector.h"
#include "linbox/vector/sparse-subvector.h"
#include "linbox/matrix/submatrix.h"

namespace LinBox
{

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag)
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size, VectorWrapper::CompareSparseEntries ());

	for (; i != i_end; ++i)
		out.push_back (typename Vector1::value_type (i->first - src_idx + dest_idx, i->second));
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag)
{
	typename Vector2::const_iterator i = in.begin () + src_idx;
	typename Vector2::const_iterator i_end = in.begin () + (src_idx + size);

	for (; i != i_end; ++i, ++dest_idx)
		out.push_back (typename Vector1::value_type (dest_idx, *i));
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag)
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size, VectorWrapper::CompareSparseEntries ());

	for (; i != i_end; ++i)
		out[i->first - src_idx + dest_idx] = i->second;
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	typedef BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> DestSubvector;
	typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> SourceSubvector;

	SourceSubvector v1 (in.begin () + src_idx, in.begin () + (src_idx + size));
	DestSubvector v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
{
	typename Vector2::const_iterator i = std::lower_bound (in.begin (), in.end (), src_idx);
	typename Vector2::const_iterator i_end = std::lower_bound (in.begin (), in.end (), src_idx + size);

	for (; i != i_end; ++i)
		out.push_back (typename Vector1::value_type (*i - src_idx + dest_idx));
}

template <class Vector>
void Splicer::append_word (Vector &v, size_t index, typename Vector::word_type word)
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

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
{
	SparseSubvector<const Vector2, VectorCategories::HybridZeroOneVectorTag> v1 (in, src_idx, src_idx + size);
	typename SparseSubvector<const Vector2, VectorCategories::HybridZeroOneVectorTag>::iterator i_v1;

	for (i_v1 = v1.begin (); i_v1 != v1.end (); ++i_v1)
		append_word (out, dest_idx + (i_v1->first << WordTraits<typename Vector2::word_type>::logof_size), i_v1->second);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
{
	SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::VectorCategory> v1 (in, src_idx, src_idx + size);
	BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	typename Vector2::const_iterator i = in.begin () + src_idx;
	typename Vector2::const_iterator i_end = in.begin () + (src_idx + size);
	size_t idx = dest_idx;

	for (; i != i_end; ++i, ++idx)
		if (*i)
			out.push_back (idx);
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	typedef BitVector<typename Vector1::Endianness> TmpVector;
	typedef BitSubvector<typename TmpVector::iterator, typename TmpVector::const_iterator> TmpSubvector;
	typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> ConstTmpSubvector;

	size_t dest_idx_offset = dest_idx & WordTraits<typename Vector2::word_type>::pos_mask;

	TmpVector tmp (round_up (size + dest_idx_offset, WordTraits<typename Vector2::word_type>::bits));
	ConstTmpSubvector v1 (in.begin () + src_idx, in.begin () + (src_idx + size));
	TmpSubvector v2 (tmp.begin () + dest_idx_offset, tmp.begin () + (dest_idx_offset + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);

	typename TmpVector::word_iterator i = tmp.word_begin ();
	size_t idx = dest_idx >> WordTraits<typename Vector1::word_type>::logof_size;

	if (i == tmp.word_end ()) {
		if (!out.empty () && out.back ().first == idx)
			out.back ().second |= tmp.back_word ();
		else if (tmp.back_word ())
			out.push_back (typename Vector1::value_type (idx, tmp.back_word ()));
	} else {
		if (!out.empty () && out.back ().first == idx) {
			out.back ().second |= *i;
			++i;
			++idx;
		}

		for (; i != tmp.word_end (); ++i, ++idx)
			if (*i)
				out.push_back (typename Vector1::value_type (idx, *i));

		if (tmp.back_word ())
			out.push_back (typename Vector1::value_type (idx, tmp.back_word ()));
	}
}

template <class Ring, class Vector1, class Vector2>
void Splicer::attach_block_specialised (const Ring &F, Vector1 &out, const Vector2 &in, size_t src_idx, size_t dest_idx, size_t size,
					VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
{
	SparseSubvector<const Vector2, typename VectorTraits<Ring, Vector2>::VectorCategory> v1 (in, src_idx, src_idx + size);
	BitSubvector<typename Vector1::iterator, typename Vector1::const_iterator> v2 (out.begin () + dest_idx, out.begin () + (dest_idx + size));

	Context<Ring> ctx (F);

	BLAS1::copy (ctx, v1, v2);
}

template <class Ring, class Matrix, class OtherMatrix>
void MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::RowMatrixTag>::copy (const Ring &F, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source)
{
	MatrixPartZero<Ring> *zero_part = dynamic_cast<MatrixPartZero<Ring> *> (source);

	if (zero_part != NULL)
		return;

	MatrixPartIdentity<Ring, Matrix, MatrixCategories::RowMatrixTag> *ident_part = dynamic_cast <MatrixPartIdentity<Ring, Matrix, MatrixCategories::RowMatrixTag> *> (source);

	if (ident_part != NULL) {
		typename Matrix::RowIterator i_A = _A->rowBegin () + horiz_block.destIndex ();
		typename Matrix::RowIterator i_A_end = _A->rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
		size_t idx = horiz_block.sourceIndex ();

		for (; i_A != i_A_end; ++i_A, ++idx)
			if (vert_block.isSourceIndexInBlock (idx))
				Splicer::attach_e_i (F, *i_A, vert_block.sourceToDestIndex (idx));

		return;
	}

	MatrixPartMatrix<Ring, OtherMatrix, Matrix, MatrixCategories::RowMatrixTag> *source_part =
		dynamic_cast<MatrixPartMatrix<Ring, OtherMatrix, Matrix, MatrixCategories::RowMatrixTag> *> (source);

	if (source_part != NULL) {
		typename Matrix::RowIterator i_A = _A->rowBegin () + horiz_block.destIndex ();
		typename Matrix::RowIterator i_A_end = _A->rowBegin () + (horiz_block.destIndex () + horiz_block.size ());
		typename OtherMatrix::ConstRowIterator i_S = source_part->_A->rowBegin () + horiz_block.sourceIndex ();

		for (; i_A != i_A_end; ++i_A, ++i_S)
			Splicer::attach_block (F, *i_A, *i_S, vert_block.sourceIndex (), vert_block.destIndex (), vert_block.size ());

		return;
	}

	throw PreconditionFailed (__FUNCTION__, __FILE__, __LINE__, "Unrecognised source-type (incompatible matrix)");
}

template <class Ring, class Matrix, class OtherMatrix>
void MatrixPartMatrix<Ring, Matrix, OtherMatrix, MatrixCategories::ColMatrixTag>::copy (const Ring &F, const Block &horiz_block, const Block &vert_block, MatrixPart<Ring> *source)
{
	MatrixPartZero<Ring> *zero_part = dynamic_cast<MatrixPartZero<Ring> *> (source);

	if (zero_part != NULL)
		return;

	MatrixPartIdentity<Ring, Matrix, MatrixCategories::ColMatrixTag> *ident_part = dynamic_cast <MatrixPartIdentity<Ring, Matrix, MatrixCategories::ColMatrixTag> *> (source);

	if (ident_part != NULL) {
		typename Matrix::ColIterator i_A = _A->colBegin () + vert_block.destIndex ();
		typename Matrix::ColIterator i_A_end = _A->colBegin () + (vert_block.destIndex () + vert_block.size ());
		size_t idx = vert_block.sourceIndex ();

		for (; i_A != i_A_end; ++i_A, ++idx)
			if (vert_block.isSourceIndexInBlock (idx))
				Splicer::attach_e_i (F, *i_A, horiz_block.sourceToDestIndex (idx));

		return;
	}

	MatrixPartMatrix<Ring, OtherMatrix, Matrix, MatrixCategories::ColMatrixTag> *source_part =
		dynamic_cast<MatrixPartMatrix<Ring, OtherMatrix, Matrix, MatrixCategories::ColMatrixTag> *> (source);

	if (source_part != NULL) {
		typename Matrix::ColIterator i_A = _A->colBegin () + vert_block.destIndex ();
		typename Matrix::ColIterator i_A_end = _A->colBegin () + (vert_block.destIndex () + vert_block.size ());
		typename OtherMatrix::ConstColIterator i_S = source_part->_A->colBegin () + vert_block.sourceIndex ();

		for (; i_A != i_A_end; ++i_A, ++i_S)
			Splicer::attach_block (F, *i_A, *i_S, horiz_block.sourceIndex (), horiz_block.destIndex (), horiz_block.size ());

		return;
	}

	throw PreconditionFailed (__FUNCTION__, __FILE__, __LINE__, "Unrecognised source-type (incompatible matrix)");
}

template <class Ring>
void Splicer::splice (const Ring &F, MatrixPart<Ring> ***input, MatrixPart<Ring> ***output) const
{
	linbox_check (!_horiz_blocks.empty ());
	linbox_check (!_vert_blocks.empty ());

	linbox_check (check ());

	commentator.start ("Splicing matrices", __FUNCTION__);

	typename std::vector<Block>::const_iterator horiz_block, vert_block;

	for (horiz_block = _horiz_blocks.begin (); horiz_block != _horiz_blocks.end (); ++horiz_block)
		for (vert_block = _vert_blocks.begin (); vert_block != _vert_blocks.end (); ++vert_block)
			output[horiz_block->dest ()][vert_block->dest ()]->copy (F, *horiz_block, *vert_block, input[horiz_block->source ()][vert_block->source ()]);

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
