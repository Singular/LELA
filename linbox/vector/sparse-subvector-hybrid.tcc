/* linbox/vector/sparse-subvector-hybrid.tcc
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

#include "linbox/vector/sparse-subvector-hybrid.h"
#include "linbox/vector/bit-subvector.h"

namespace LinBox {

template <class Vector>
void SparseSubvector<const Vector, VectorCategories::HybridZeroOneVectorTag>::set_start_end (typename Vector::const_iterator begin,
											     typename Vector::const_iterator end)
{
	_begin = std::lower_bound (begin, end, _start >> WordTraits<word_type>::logof_size, VectorWrapper::CompareSparseEntries ());
	_end = std::upper_bound (begin, end, _finish >> WordTraits<word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

	if (_end != _begin && (static_cast<size_t> ((_end - 1)->first) << WordTraits<word_type>::logof_size) + (_start & WordTraits<word_type>::pos_mask) >= _finish)
		_end_marker = _end - 1;
	else
		_end_marker = _end;
}

template <class Vector>
typename SparseSubvector<const Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator &SparseSubvector<const Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator::operator ++ ()
{
	size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
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
bool SparseSubvector<const Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator::operator != (const const_iterator &c)
{
	if (c._pos == _v->_end_marker) {
		if (_pos == _v->_end)
			return false;
		else {
			update_ref ();
			size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
			return virtual_index < _v->_finish;
		}
	}
	else
		return _pos != c._pos;
}

template <class Vector>
void SparseSubvector<const Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator::update_ref ()
{
	size_t shift, virtual_index, col_index;

	if (!_ref_first_valid || !_ref_second_valid) {
		shift = _v->_start & WordTraits<word_type>::pos_mask;

		virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
		col_index = static_cast<size_t> (_pos->first) << WordTraits<word_type>::logof_size;
	}

	if (!_ref_first_valid) {
		if (_start) {
			if (col_index < _v->_start)
				_ref.first = 0;
			else
				_ref.first = (col_index - _v->_start) >> WordTraits<word_type>::logof_size;
			
			_start = false;
		}
		else if (virtual_index + WordTraits<word_type>::bits < col_index)
			_ref.first = (col_index - _v->_start) >> WordTraits<word_type>::logof_size;
		else
			++_ref.first;

		virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;

		_ref_first_valid = true;
	}

	if (!_ref_second_valid) {
		if (shift == 0)
			_ref.second = _pos->second;
		else if (virtual_index < col_index)
			_ref.second = Endianness::shift_right (_pos->second, WordTraits<word_type>::bits - shift);
		else if ((_pos + 1 == _v->_end && _v->_end_is_end) || (_pos + 1)->first > _pos->first + 1)
			_ref.second = Endianness::shift_left (_pos->second, shift);
		else
			_ref.second = Endianness::shift_left (_pos->second, shift) | Endianness::shift_right ((_pos + 1)->second, WordTraits<word_type>::bits - shift);

		if ((static_cast<size_t> (_ref.first) + 1) << WordTraits<word_type>::logof_size > _v->_finish - _v->_start)
			_ref.second &= Endianness::mask_left ((_v->_finish - _v->_start) & WordTraits<word_type>::pos_mask);

		_ref_second_valid = true;
	}
}

} // namespace LinBox

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
