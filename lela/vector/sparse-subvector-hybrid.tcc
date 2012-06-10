/* lela/vector/sparse-subvector-hybrid.tcc
 * Copyright 2011 Bradford Hovinen
 *
 * Evolved from sparse-subvector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SPARSE_SUBVECTOR_HYBRID_TCC
#define __VECTOR_SPARSE_SUBVECTOR_HYBRID_TCC

#include <utility>

#include "lela/vector/sparse-subvector-hybrid.h"
#include "lela/vector/bit-subvector.h"

namespace LELA {

template <class Vector>
void SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01>::set_start_end (typename Vector::const_iterator begin,
										  typename Vector::const_iterator end)
{
	_begin = std::lower_bound (begin, end, _start >> WordTraits<word_type>::logof_size, VectorUtils::FindSparseEntryLB ());
	_end = end;
}

template <class Vector>
typename SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01>::const_iterator &SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01>::const_iterator::operator ++ ()
{
	size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _start_idx;
	size_t col_index = static_cast<size_t> (_pos->first) << WordTraits<word_type>::logof_size;

	if (virtual_index < col_index) {
		++_ref.first;
		_ref_second_valid = false;
	} else {
		++_pos;
		_ref_first_valid = false;
		_ref_second_valid = false;
	}

	return *this;
}

template <class Vector>
bool SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01>::const_iterator::operator == (const const_iterator &c)
{
	if (c._end) {
		if (_pos == c._pos)
			return true;
		else {
			update_ref ();
			size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _start_idx;
			return virtual_index >= _end_idx;
		}
	}
	else
		return _pos == c._pos;
}

template <class Vector>
void SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01>::const_iterator::update_ref ()
{
	size_t shift, virtual_index, col_index;

	if (!_ref_first_valid || !_ref_second_valid) {
		shift = _start_idx & WordTraits<word_type>::pos_mask;

		virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _start_idx;
		col_index = static_cast<size_t> (_pos->first) << WordTraits<word_type>::logof_size;
	}

	if (!_ref_first_valid) {
		if (_start) {
			if (col_index < _start_idx)
				_ref.first = 0;
			else
				_ref.first = (col_index - _start_idx) >> WordTraits<word_type>::logof_size;
			
			_start = false;
		}
		else if (virtual_index + WordTraits<word_type>::bits < col_index)
			_ref.first = (col_index - _start_idx) >> WordTraits<word_type>::logof_size;
		else
			++_ref.first;

		virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _start_idx;

		_ref_first_valid = true;
	}

	if (!_ref_second_valid) {
		typename Endianness::word_pair w;

		if (virtual_index < col_index) {
			w.parts.low = 0ULL;
			w.parts.high = _pos->second;
		}
		else if ((_pos + 1)->first > _pos->first + 1) {
			w.parts.low = _pos->second;
			w.parts.high = 0ULL;
		}
		else {
			w.parts.low = _pos->second;
			w.parts.high = (_pos + 1)->second;
		}

		w.full = Endianness::shift_left (w.full, shift);
		_ref.second = w.parts.low;

		if ((static_cast<size_t> (_ref.first) + 1) << WordTraits<word_type>::logof_size > _end_idx - _start_idx)
			_ref.second &= Endianness::mask_left ((_end_idx - _start_idx) & WordTraits<word_type>::pos_mask);

		_ref_second_valid = true;
	}
}

} // namespace LELA

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
