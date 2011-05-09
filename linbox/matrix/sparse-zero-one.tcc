/* linbox/matrix/sparse-zero-one.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of SparseMatrix and helpers for 0-1 matrices
 * 
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_sparse_zero_one_TCC
#define __LINBOX_matrix_sparse_zero_one_TCC

namespace LinBox
{

template <class Element, class Row>
bool SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag>::getEntry (Element &x, size_t i, size_t j) const
{
	const Row &v = (*this)[i];

	typename Row::const_iterator idx = std::lower_bound (v.begin (), v.end (), j);

	if (idx != v.end () && *idx == j) {
		x = true;
		return true;
	} else
		return false;
}

template <class Element, class Row>
bool SparseMatrix<Element, Row, VectorCategories::HybridZeroOneVectorTag>::getEntry (Element &x, size_t i, size_t j) const
{
	const Row &v = (*this)[i];

	typename Row::const_iterator idx;

	idx = std::lower_bound (v.begin (), v.end (), j >> WordTraits<typename Row::word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

	if (idx != v.end () && idx->first == j >> WordTraits<typename Row::word_type>::logof_size) {
		x = idx->second & Row::Endianness::e_j (j & WordTraits<typename Row::word_type>::pos_mask);
		return true;
	} else
		return false;
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag>::eraseEntry (size_t i, size_t j)
{
	Row &v = (*this)[i];

	typename Row::iterator idx = std::lower_bound (v.begin (), v.end (), j);

	if (idx != v.end () && *idx == j)
		v.erase (idx);
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	Row &v = _A[i];
	typename Row::iterator iter;

	if (value && v.size () == 0) {
		v.push_back (j);
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j);

		if (value && (iter == v.end () || *iter != j))
			iter = v.insert (iter, j);
		else if (!value && iter != v.end () && *iter == j)
			v.erase (iter);
	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::HybridZeroOneVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	Row &v = _A[i];
	typename Row::iterator it;

	typename Row::word_type m = Row::Endianness::e_j (j & WordTraits<typename Row::word_type>::pos_mask);

	if (value && v.empty ()) {
		v.push_back (typename Row::value_type (j >> WordTraits<typename Row::word_type>::logof_size, m));
	} else {
		it = std::lower_bound (v.begin (), v.end (), j >> WordTraits<typename Row::word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

		if (it == v.end () || it->first != (j >> WordTraits<typename Row::word_type>::logof_size)) {
			if (value)
				it = v.insert (it, typename Row::value_type (j >> WordTraits<typename Row::word_type>::logof_size, m));
		}
		else {
			if (value)
				it->second |= m;
			else {
				it->second &= ~m;

				if (!it->second)
					v.erase (it);
			}
		}
	}
}

} // namespace LinBox

#endif // __LINBOX_matrix_sparse_zero_one_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
