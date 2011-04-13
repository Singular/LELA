/* linbox/algorithms/gauss-jordan.inl
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Gauss-Jordan elimination
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ALGORITHMS_GAUSS_JORDAN_INL
#define __LINBOX_ALGORITHMS_GAUSS_JORDAN_INL

#include "linbox/algorithms/gauss-jordan.h"

#ifdef DETAILED_PROFILE
#  define TIMER_DECLARE(part) LinBox::UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_START(part) part##_timer.start ()
#  define TIMER_STOP(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define TIMER_REPORT(part) \
	commentator.report (Commentator::LEVEL_NORMAL, TIMING_MEASURE) \
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
#else
#  define TIMER_DECLARE(part)
#  define TIMER_START(part)
#  define TIMER_STOP(part)
#  define TIMER_REPORT(part)
#endif

#ifndef PROGRESS_STEP
#  define PROGRESS_STEP 1024
#endif // PROGRESS_STEP

namespace LinBox
{

template <class Field>
template <class Matrix>
int GaussJordan<Field>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
					     VectorCategories::DenseVectorTag) const
{
	typename Matrix::ConstRowIterator i;

	for (; col < A.coldim (); ++col) {
		int k = start_row;

		for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k)
			if (!F.isZero ((*i)[col]))
				return k;
	}

	return -1;
}

template <class Field>
template <class Matrix>
int GaussJordan<Field>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
					     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Matrix::ConstRowIterator i;

	for (; col < A.coldim (); ++col) {
		int k = start_row;

		for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k)
			if ((*i)[col])
				return k;
	}

	return -1;
}

template <class Field>
template <class Matrix>
int GaussJordan<Field>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
					     VectorCategories::SparseVectorTag) const
{
	typename SparseMatrix::ConstRowIterator i;

	size_t min_nonzero = (size_t) -1, pivot = -1, k = start_row, s;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty ()) {
			if (i->front ().first < col) {
				col = i->front ().first;
				min_nonzero = i->size ();
				pivot = k;
			}
			else if (i->front ().first == col && (s = i->size ()) < min_nonzero) {
				min_nonzero = s;
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Field>
template <class Matrix>
int GaussJordan<Field>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
					     VectorCategories::SparseZeroOneVectorTag) const
{
	typename SparseMatrix::RowIterator i;

	size_t min_nonzero = (size_t) -1, pivot = -1, k = start_row, s;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty ()) {
			if (i->front () < col) {
				col = i->front ();
				min_nonzero = i->size ();
				pivot = k;
			}
			else if (i->front () == col && (s = i->size ()) < min_nonzero) {
				min_nonzero = s;
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Field>
template <class Matrix>
int GaussJordan<Field>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
					     VectorCategories::HybridZeroOneVectorTag) const
{
	typename SparseMatrix::ConstRowIterator i;
	typename Field::Element a;

	size_t min_blocks = (size_t) -1, pivot = -1, k = start_row, s;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty () && i->front ().first << WordTraits<typename Matrix::Row::word_type>::logof_size <= (int) col) {
			int idx = VD.firstNonzeroEntry (a, *i);
			if ((size_t) idx < col) {
				col = idx;
				min_blocks = i->size ();
				pivot = k;
			}
			else if ((size_t) idx == col && (s = i->size ()) < min_blocks) {
				min_blocks = s;
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Field>
template <class Matrix>
void GaussJordan<Field>::SetIdentity (Matrix &U, size_t start_row) const
{
	StandardBasisStream<Field, typename Matrix::Row> stream (F, U.coldim ());
	typename Matrix::RowIterator i = U.rowBegin () + start_row;

	for (; i != U.rowEnd (); ++i)
		stream >> *i;
}

template <class Field>
void GaussJordan<Field>::GaussJordanTransform (DenseMatrix  &A,
					       int           k,
					       Element       d_0,
					       DenseMatrix  &U,
					       Permutation  &P,
					       size_t       &r,
					       int          &h,
					       Element      &d,
					       DenseMatrix  &S,
					       DenseMatrix  &T) const
{
	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	// DEBUG
	// report << "enter" << std::endl;
	// report << "A =" << std::endl;
	// MD.write (report, A);
	// report << "U =" << std::endl;
	// MD.write (report, U);
	// report << "k = " << k << ", d_0 = " << d_0 << std::endl;

	DenseMatrix Aw (A, k, 0, A.rowdim () - k, A.coldim ());

	if (MD.isZero (Aw)) {
		// DEBUG
		// report << "A is 0" << std::endl;

		r = 0;
		h = A.rowdim () - k;
		F.assign (d, d_0);
	}
	else if (A.coldim () <= _cutoff) {
		// DEBUG
		// report << "A.coldim <= " << _cutoff << ", using standard elimination" << std::endl;

		size_t l = P.size ();

		// DEBUG
		// report << "U before elimination:" << std::endl;
		// MD.write (report, U);

		Submatrix<DenseMatrix> Up (U, 0, k, U.rowdim (), U.coldim () - k);

		StandardRowEchelonForm (A, Up, P, r, d, true, true, k);

		// DEBUG
		// report << "U after elimination:" << std::endl;
		// MD.write (report, U);

		if (!P.empty ())
			h = std::min_element (P.begin () + l, P.end (), CompareSecond ())->second;
		else
			h = k;
	}
	else {
		// DEBUG
		// report << "A.coldim () > " << _cutoff << std::endl;

		int m_1 = round_up (A.coldim () / 2, _cutoff), m_2 = A.coldim () - m_1;
		DenseMatrix A_1 (A, 0, 0,   A.rowdim (), m_1);
		DenseMatrix A_2 (A, 0, m_1, A.rowdim (), m_2);

		size_t r_1, r_2;
		int h_1, h_2;
		Element d_1;

		GaussJordanTransform (A_1, k, d_0, U, P, r_1, h_1, d_1, S, T);

		// DEBUG
		// report << "Status after first recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// MD.write (report, U);
		// report << "P = ";
		// MD.writePermutation (report, P) << std::endl;
		// report << "R =" << std::endl;
		// MD.write (report, R);
		// report << "r_1 = " << r_1 << std::endl;
		// report << "d_1 = " << d_1 << std::endl;

		Submatrix<DenseMatrix> U_2  (U, 0,       k, U.rowdim (),           r_1);
		Submatrix<DenseMatrix> U_21 (U, 0,       k, k,                     r_1);
		Submatrix<DenseMatrix> U_22 (U, k,       k, r_1,                   r_1);
		Submatrix<DenseMatrix> U_23 (U, r_1 + k, k, U.rowdim () - r_1 - k, r_1);

		T.resize (r_1, m_2);

		DenseMatrix A_21    (A, 0,       m_1,     k,                     m_2);
		DenseMatrix A_22    (A, k,       m_1,     r_1,                   m_2);
		DenseMatrix A_23    (A, k + r_1, m_1,     A.rowdim () - r_1 - k, m_2);

		MD.permuteRows (A_2, P.begin (), P.end ());

		MD.copy (T, A_22);
		MD.gemm (F.one (), U_21, T, F.one (),  A_21);
		MD.gemm (F.one (), U_22, T, F.zero (), A_22);
		MD.gemm (F.one (), U_23, T, F.one (),  A_23);

		Permutation P_2;

		// DEBUG
		// report << "Status before second recursive call:" << std::endl;
		// report << "B =" << std::endl;
		// MD.write (report, R_2);

		GaussJordanTransform (A_2, k + r_1, d_1, U, P_2, r_2, h_2, d, S, T);

		// DEBUG
		// report << "Status after second recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// MD.write (report, U);
		// report << "R =" << std::endl;
		// MD.write (report, R);
		// report << "r_2 = " << r_2 << std::endl;
		// report << "d = " << d << std::endl;

		Submatrix<DenseMatrix> P_2U_2       (U,  0,             k,       U.rowdim (),                  r_1);
		Submatrix<DenseMatrix> U_212        (U,  0,             k,       k + r_1,                      r_1);
		Submatrix<DenseMatrix> P_2U_23      (U,  k + r_1,       k,       r_2,                          r_1);
		Submatrix<DenseMatrix> P_2U_24      (U,  k + r_1 + r_2, k,       U.rowdim () - r_1 - r_2 - k,  r_1);
		Submatrix<DenseMatrix> U_312        (U,  0,             k + r_1, k + r_1,                      r_2);
		Submatrix<DenseMatrix> U_33         (U,  k + r_1,       k + r_1, r_2,                          r_2);
		Submatrix<DenseMatrix> U_34         (U,  k + r_1 + r_2, k + r_1, U.rowdim () - r_1 - r_2 - k,  r_2);

		// DEBUG
		// report << "U_2 =" << std::endl;
		// MD.write (report, U);
		// report << "P_2 = ";
		// MD.writePermutation (report, P_2) << std::endl;

		MD.permuteRows (P_2U_2, P_2.begin (), P_2.end ());
				
		// DEBUG
		// report << "P_2U_2 =" << std::endl;
		// MD.write (report, P_2U_2);
		// report << "P_2U_23 =" << std::endl;
		// MD.write (report, P_2U_23);
		// report << "U_3 =" << std::endl;;
		// MD.write (report, U_3);

		S.resize (r_2, r_1);
		MD.copy (S, P_2U_23);

		MD.gemm (F.one (), U_312, S, F.one (), U_212);
		MD.gemm (F.one (), U_33, S, F.zero (), P_2U_23);
		MD.gemm (F.one (), U_34, S, F.one (), P_2U_24);

		// DEBUG
		// report << "U_3P_2U_23 =" << std::endl;
		// MD.write (report, U_3P_2U_23);

		P.insert (P.end (), P_2.begin (), P_2.end ());
		r = r_1 + r_2;
		h = std::min (h_1, h_2);
	}

	// DEBUG
	// report << "Status at end:" << std::endl;
	// report << "U =" << std::endl;
	// MD.write (report, U);
	// report << "P = ";
	// MD.writePermutation (report, P) << std::endl;
	// report << "R =" << std::endl;
	// MD.write (report, R);
	// report << "k = " << k << ", r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
}

template <class Field>
void GaussJordan<Field>::GaussTransform (DenseMatrix             &A,
					 Element                  d_0,
					 Submatrix<DenseMatrix>  &U,
					 Permutation             &P,
					 size_t                  &r,
					 int                     &h,
					 Element                 &d,
					 DenseMatrix             &T) const
{
	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	// DEBUG
	// report << "enter" << std::endl;
	// report << "A =" << std::endl;
	// MD.write (report, A);
	// report << "U =" << std::endl;
	// MD.write (report, U);
	// report << "k = " << k << ", d_0 = " << d_0 << std::endl;

	if (MD.isZero (A)) {
		// DEBUG
		// report << "A is 0" << std::endl;

		r = 0;
		h = A.rowdim ();
		F.assign (d, d_0);
	}
	else if (A.coldim () <= _cutoff) {
		// DEBUG
		// report << "A.coldim <= " << _cutoff << ", using standard elimination" << std::endl;

		size_t l = P.size ();

		// DEBUG
		// report << "U before elimination:" << std::endl;
		// MD.write (report, U);

		StandardRowEchelonForm (A, U, P, r, d, false, true);

		// DEBUG
		// report << "U after elimination:" << std::endl;
		// MD.write (report, U);

		if (!P.empty ())
			h = std::min_element (P.begin () + l, P.end (), CompareSecond ())->second;
		else
			h = 0;
	}
	else {
		// DEBUG
		// report << "A.coldim () > " << _cutoff << std::endl;

		int m_1 = round_up (A.coldim () / 2, _cutoff), m_2 = A.coldim () - m_1;
		DenseMatrix A_1 (A, 0, 0, A.rowdim (), m_1);

		size_t r_1, r_2;
		int h_1, h_2;
		Element d_1;

		GaussTransform (A_1, d_0, U, P, r_1, h_1, d_1, T);

		// DEBUG
		// report << "Status after first recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// MD.write (report, U);
		// report << "P = ";
		// MD.writePermutation (report, P) << std::endl;
		// report << "R =" << std::endl;
		// MD.write (report, R);
		// report << "r_1 = " << r_1 << std::endl;
		// report << "d_1 = " << d_1 << std::endl;

		Submatrix<DenseMatrix> U_11 (U, 0,   0,   r_1,               r_1);
		Submatrix<DenseMatrix> U_12 (U, r_1, 0,   U.rowdim () - r_1, r_1);
		Submatrix<DenseMatrix> U_22 (U, r_1, r_1, U.rowdim () - r_1, U.coldim () - r_1);

		DenseMatrix A_2  (A, 0,   m_1, A.rowdim (),       m_2);
		DenseMatrix A_21 (A, 0,   m_1, r_1,               m_2);
		DenseMatrix A_22 (A, r_1, m_1, A.rowdim () - r_1, m_2);

		MD.permuteRows (A_2, P.begin (), P.end ());

		MD.gemm (F.one (), U_12, A_21, F.one (), A_22);

		// FIXME: This should be replaced by trmm when it becomes available
		T.resize (r_1, m_2);
		MD.copy (T, A_21);
		MD.gemm (F.one (), U_11, T, F.zero (), A_21);

		Permutation P_2;

		// DEBUG
		// report << "Status before second recursive call:" << std::endl;
		// report << "B =" << std::endl;
		// MD.write (report, R_2);

		GaussTransform (A_22, d_1, U_22, P_2, r_2, h_2, d, T);

		// DEBUG
		// report << "Status after second recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// MD.write (report, U);
		// report << "R =" << std::endl;
		// MD.write (report, R);
		// report << "r_2 = " << r_2 << std::endl;
		// report << "d = " << d << std::endl;

		// DEBUG
		// report << "U_2 =" << std::endl;
		// MD.write (report, U);
		// report << "P_2 = ";
		// MD.writePermutation (report, P_2) << std::endl;

		MD.permuteRows (U_12, P_2.begin (), P_2.end ());

		// FIXME: This should be replaced by trmm when it becomes available
		T.resize (U_12.rowdim (), U_12.coldim ());
		MD.copy (T, U_12);
		MD.gemm (F.one (), U_22, T, F.zero (), U_12);
				
		// DEBUG
		// report << "P_2U_2 =" << std::endl;
		// MD.write (report, P_2U_2);

		size_t P_len = P.size ();
		P.insert (P.end (), P_2.begin (), P_2.end ());

		// Update indices in permutation to refer to whole matrix
		typename Permutation::iterator i;

		for (i = P.begin () + P_len; i != P.end (); ++i) {
			i->first += r_1;
			i->second += r_1;
		}

		r = r_1 + r_2;
		h = std::min (h_1, h_2);
	}

	// DEBUG
	// report << "Status at end:" << std::endl;
	// report << "U =" << std::endl;
	// MD.write (report, U);
	// report << "P = ";
	// MD.writePermutation (report, P) << std::endl;
	// report << "R =" << std::endl;
	// MD.write (report, R);
	// report << "k = " << k << ", r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
}

template <class Field>
template <class Vector>
Vector &GaussJordan<Field>::FastAddinSpecialised (Vector &v, const Vector &w, size_t idx,
						  VectorCategories::SparseZeroOneVectorTag) const
{
	static Vector tmp;

	VD.add (tmp, Subvector<typename Vector::iterator> (v.begin () + idx, v.end ()), w);
	v.resize (idx + tmp.size ());
	std::copy (tmp.begin (), tmp.end (), v.begin () + idx);

	return v;
}

template <class Field>
template <class Vector>
Vector &GaussJordan<Field>::FastAddinSpecialised (Vector &v, const Vector &w, size_t idx,
						  VectorCategories::HybridZeroOneVectorTag) const
{
	static Vector tmp;

	Subvector<typename Vector::iterator> t (v.begin () + idx, v.end ());

	VD.add (tmp, t, w);
	v.resize (idx + tmp.size ());
	std::copy (tmp.begin (), tmp.end (), v.begin () + idx);

	return v;
}

template <class Field>
bool GaussJordan<Field>::testFastAddinHybridVector () const
{
	commentator.start ("Testing FastAddin", __FUNCTION__);

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Vector<GF2>::Hybrid v, w;

	error << std::hex << std::setfill ('0');

	v.push_back (Vector<GF2>::Hybrid::value_type (0, 0xffff0000ffff0000ULL));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, 0xffff0000ffff0000ULL));

	w.push_back (Vector<GF2>::Hybrid::value_type (1, 0x00ffff0000ffff00ULL));

	FastAddin (v, w, 1);

	if (v.front ().second != 0xffff0000ffff0000ULL) {
		error << "Test 1 not okay: first word is " << std::setw (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 4) << v.front ().second << " but should be ffff0000ffff0000" << std::endl;
		pass = false;
	}

	if ((v.begin () + 1)->second != 0xff00ff00ff00ff00ULL) {
		error << "Test 2 not okay: second word is " << std::setw (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 4) << (v.begin () + 1)->second << " but should be ff00ff00ff00ff00" << std::endl;
		pass = false;
	}

	v.clear ();
	w.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (0, 0xffff0000ffff0000ULL));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, 0xffff0000ffff0000ULL));

	w.push_back (Vector<GF2>::Hybrid::value_type (0, 0x00ffff0000ffff00ULL));

	FastAddin (v, w, 0);

	if (v.front ().second != 0xff00ff00ff00ff00ULL) {
		error << "Test 3 not okay: first word is " << std::setw (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 4) << v.front ().second << " but should be ff00ff00ff00ff00" << std::endl;
		pass = false;
	}

	if ((v.begin () + 1)->second != 0xffff0000ffff0000ULL) {
		error << "Test 4 not okay: second word is " << std::setw (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 4) << (v.begin () + 1)->second << " but should be ffff0000ffff0000" << std::endl;
		pass = false;
	}

	v.clear ();
	w.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (0, 0xffff0000ffff0000ULL));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, 0xffff0000ffff0000ULL));

	w.push_back (Vector<GF2>::Hybrid::value_type (0, 0x00ffff0000ffff00ULL));
	w.push_back (Vector<GF2>::Hybrid::value_type (1, 0x00ffff0000ffff00ULL));

	FastAddin (v, w, 0);

	if (v.front ().second != 0xff00ff00ff00ff00ULL) {
		error << "Test 5 not okay: first word is " << std::setw (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 4) << v.front ().second << " but should be ff00ff00ff00ff00" << std::endl;
		pass = false;
	}

	if ((v.begin () + 1)->second != 0xff00ff00ff00ff00ULL) {
		error << "Test 6 not okay: second word is " << std::setw (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 4) << (v.begin () + 1)->second << " but should be ff00ff00ff00ff00" << std::endl;
		pass = false;
	}

	error << std::dec << std::setfill (' ');

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
							  VectorCategories::DenseVectorTag) const
{
	commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	if (rank == 0)
		return A;

	typename Matrix1::RowIterator i_A, j_A;
	typename Matrix2::RowIterator i_L;

	if (compute_L)
		i_L = L.rowBegin () + (rank - 1);

	int current_row = rank - 1, elim_row;

	typename Field::Element a;

	i_A = A.rowBegin () + current_row;

	do {
		// DEBUG
		// report << "Row " << current_row << ", current A:" << std::endl;
		// MD.write (report, A);

		for (elim_row = rank - 1, j_A = A.rowBegin () + elim_row; elim_row > std::max (current_row, (int) start_row - 1); --elim_row, --j_A) {
			size_t col = VD.firstNonzeroEntry (a, *j_A);

			if (!F.isZero ((*i_A)[col])) {
				// DEBUG
				// report << "Eliminating " << current_row << " from " << elim_row << std::endl;

				VD.addin (*i_A, *j_A);

				if (compute_L)
					VD.addin (*i_L, *(L.rowBegin () + elim_row));
			}
		}

		if (compute_L)
			--i_L;

		if ((rank - current_row) % PROGRESS_STEP == 0)
			commentator.progress ();

		--current_row;
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "GaussJordan::ReduceRowEchelonSpecialised");

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
							  VectorCategories::SparseVectorTag) const
{
	commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

	typename Matrix1::RowIterator i_A, j_A, i_Ae;
	typename Matrix2::RowIterator i_L, j_L;

	Element negx, xinv;

	if (compute_L)
		i_L = L.rowBegin () + (rank - 1);

	i_A = A.rowBegin () + (rank - 1);

	do {
		if (i_A - A.rowBegin () >= (long) start_row)
			i_Ae = i_A;

		for (j_A = A.rowBegin () + (rank - 1), j_L = L.rowBegin () + (rank - 1); j_A != i_Ae; --j_A, --j_L) {
			// We must start over each time because the operations below invalidate iterators...
			std::reverse_iterator<typename SparseMatrix::Row::iterator> i (i_A->end ());
			std::reverse_iterator<typename SparseMatrix::Row::iterator> i_stop (i_A->begin ());

			while (i != i_stop && i->first > j_A->front ().first)
				++i;

			if (i != i_stop && i->first == j_A->front ().first) {
				F.neg (negx, i->second);

				VD.axpyin (*i_A, negx, *j_A);

				if (compute_L)
					VD.axpyin (*i_L, negx, *j_L);
			}
		}

		F.inv (xinv, i_A->front ().second);
		VD.mulin (*i_A, xinv);

		if (compute_L) {
			VD.mulin (*i_L, xinv);
			--i_L;
		}

		if ((A.rowEnd () - i_A) % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "GaussJordan::ReduceRowEchelonSpecialised");

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
							  VectorCategories::SparseZeroOneVectorTag) const
{
	commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

	typename Matrix1::RowIterator i_A, j_A, i_Ae;
	typename Matrix2::RowIterator i_L, j_L;

	if (compute_L)
		i_L = L.rowBegin () + (rank - 1);

	i_A = A.rowBegin () + (rank - 1);

	do {
		if (i_A - A.rowBegin () >= (long) start_row)
			i_Ae = i_A;

		size_t prev_idx = i_A->size ();

		if (compute_L)
			j_L = L.rowBegin () + (rank - 1);

		for (j_A = A.rowBegin () + (rank - 1); j_A != i_Ae; --j_A) {
			if (prev_idx > i_A->size ())
				prev_idx = i_A->size ();

			// We must start over each time because the operations below invalidate iterators...
			std::reverse_iterator<typename SparseMatrix::Row::iterator> i_idx (i_A->begin () + prev_idx);
			std::reverse_iterator<typename SparseMatrix::Row::iterator> i_stop (i_A->begin ());

			while (i_idx != i_stop && *i_idx > j_A->front ()) {
				++i_idx;
				--prev_idx;
			}

			if (i_idx != i_stop && *i_idx == j_A->front ()) {
				FastAddin (*i_A, *j_A, prev_idx - 1);

				if (compute_L)
					VD.addin (*i_L, *j_L);
			}

			if (compute_L)
				--j_L;
		}

		if (compute_L)
			--i_L;

		if ((A.rowEnd () - i_A) % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "GaussJordan::ReduceRowEchelonSpecialised");

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
							  VectorCategories::HybridZeroOneVectorTag) const
{
	commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	if (rank == 0)
		return A;

	typename Matrix1::RowIterator i_A, i_Ae, j_A;
	typename Matrix2::RowIterator i_L;

	typename Matrix1::Row::iterator i;

	typename Matrix1::Row::word_type v, v1, mask;

	if (compute_L)
		i_L = L.rowBegin () + (rank - 1);

	i_A = A.rowBegin () + (rank - 1);

	int current_row = rank - 1, elim_row;

	do {
		// DEBUG
		// report << "Row " << current_row << ", state of A:" << std::endl;
		// MD.write (report, A);

		if (i_A - A.rowBegin () >= (int) start_row)
			i_Ae = i_A;

		size_t prev_idx = i_A->size () - 1;
		elim_row = rank - 1;

		j_A = A.rowBegin () + elim_row;

		while (elim_row > std::max (current_row, (int) start_row - 1)) {
			if (prev_idx >= i_A->size ())
				prev_idx = i_A->size () - 1;

			typename Matrix1::Row::index_type j_idx = j_A->front ().first;

			// Note: We don't need to test for validity because, if the input
			// is valid, i_idx will *never* go past the beginning
			for (i = i_A->begin () + prev_idx; i->first > j_idx; --i, --prev_idx);

			if (i->first == j_A->front ().first) {
				v = j_A->front ().second;
				mask = Adaptor<Field>::Endianness::first_position (v);

				if ((i_A->begin () + prev_idx)->second & mask) {
					// DEBUG
					// report << "ReduceRowEchelonSpecialised: eliminating row " << current_row << " from row " << elim_row << std::endl;

					FastAddin (*i_A, *j_A, prev_idx);

					if (compute_L)
						VD.addin (*i_L, *(L.rowBegin () + elim_row));
				}
			}

			v1 = (i_A->begin () + prev_idx)->second;
			mask = 0;

			while (j_A->front ().first >= (*i_A)[prev_idx].first && !(v1 & mask)) {
				--j_A;
				--elim_row;

				if (j_A->front ().first > (*i_A)[prev_idx].first) {
					j_A = std::upper_bound (A.rowBegin () + current_row + 1, A.rowBegin () + elim_row, (*i_A)[prev_idx].first,
								PivotRowCompare<typename SparseMatrix::Row> ());
					--j_A;
					elim_row = j_A - A.rowBegin ();

					if (elim_row == current_row)
						break;
				}

				v = j_A->front ().second;
				mask = Adaptor<Field>::Endianness::first_position (v);
			}
		}

		if (compute_L)
			--i_L;

		if ((rank - current_row) % PROGRESS_STEP == 0)
			commentator.progress ();

		--current_row;
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "GaussJordan::ReduceRowEchelonSpecialised");

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
							  VectorCategories::DenseZeroOneVectorTag) const
{
	commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	if (rank == 0)
		return A;

	typename Matrix1::RowIterator i_A, j_A;
	typename Matrix2::RowIterator i_L;

	if (compute_L)
		i_L = L.rowBegin () + (rank - 1);

	int current_row = rank - 1, elim_row;

	typename Field::Element a;

	i_A = A.rowBegin () + current_row;

	do {
		// DEBUG
		// report << "Row " << current_row << ", current A:" << std::endl;
		// MD.write (report, A);

		for (elim_row = rank - 1, j_A = A.rowBegin () + elim_row; elim_row > std::max (current_row, (int) start_row - 1); --elim_row, --j_A) {
			size_t col = VD.firstNonzeroEntry (a, *j_A);

			if ((*i_A)[col]) {
				// DEBUG
				// report << "Eliminating " << current_row << " from " << elim_row << std::endl;

				VD.addin (*i_A, *j_A);

				if (compute_L)
					VD.addin (*i_L, *(L.rowBegin () + elim_row));
			}
		}

		if (compute_L)
			--i_L;

		if ((rank - current_row) % PROGRESS_STEP == 0)
			commentator.progress ();

		--current_row;
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "GaussJordan::ReduceRowEchelonSpecialised");

	return A;
}

template <class Field>
void GaussJordan<Field>::DenseRowEchelonForm (DenseMatrix &A,
					      DenseMatrix &U,
					      Permutation &P,
					      size_t      &rank,
					      Element     &det,
					      bool         reduced)
{
	linbox_check (U.rowdim () == A.rowdim ());
	linbox_check (U.rowdim () == U.coldim ());

	commentator.start ("Dense row-echelon form", "GaussJordan::DenseRowEchelonForm");

	int h;

	SetIdentity (U);  // Necessary in case where A.rowdim () > A.coldim ()
	P.clear ();

	if (reduced) {
		DenseMatrix S, T;

		GaussJordanTransform (A, 0, F.one (), U, P, rank, h, det, S, T);
	} else {
		DenseMatrix T;
		Submatrix<DenseMatrix> Up (U, 0, 0, U.rowdim (), U.coldim ());

		GaussTransform (A, F.one (), Up, P, rank, h, det, T);
	}

	commentator.stop (MSG_DONE, NULL, "GaussJordan::DenseRowEchelonForm");
}

template <class Field>
template <class Matrix1, class Matrix2>
void GaussJordan<Field>::StandardRowEchelonForm (Matrix1       &A,
						 Matrix2       &L,
						 Permutation   &P,
						 size_t        &rank,
						 Element       &det,
						 bool           reduced,
						 bool           compute_L,
						 size_t         start_row) const
{
	commentator.start ("Standard row-echelon form", __FUNCTION__, A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	TIMER_DECLARE(GetPivot);
	TIMER_DECLARE(Permute);
	TIMER_DECLARE(ElimBelow);

	typename Matrix1::RowIterator i_A, j_A;

	size_t col = 0;
	int k = start_row;
	Element a, x, xinv, negxinv, negaxinv;

	F.init (a, 0);
	F.init (x, 0);

	typename Matrix2::RowIterator i_L, j_L;

	if (compute_L)
		SetIdentity (L, start_row);

	P.clear ();
	rank = 0;
	F.init (det, 1);

	for (i_A = A.rowBegin () + k, i_L = L.rowBegin () + k; i_A != A.rowEnd (); ++k, ++i_A, ++i_L) {
		TIMER_START(GetPivot);
		int pivot = GetPivot (A, k, col);
		TIMER_STOP(GetPivot);

		if (pivot == -1)
			break;

		TIMER_START(Permute);
		if (k != pivot) {
			// DEBUG
			// report << "Permuting " << k << " and " << pivot << std::endl;
			// report << "L before permutation:" << std::endl;
			// MD.write (report, L);

			Transposition t (k, pivot);
			P.push_back (t);
			MD.permuteRows (A, &t, &t + 1);

			if (compute_L) {
				typename SubmatrixTypename<Matrix2>::Type Lp (L, 0, 0, L.rowdim (), k - start_row);
				MD.permuteRows (Lp, &t, &t + 1);
			}

			// DEBUG
			// report << "L after permutation:" << std::endl;
			// MD.write (report, L);
		}
		TIMER_STOP(Permute);

		// DEBUG
		// report << "Row " << k << ", pivot is " << pivot << std::endl;
		// report << "Current A:" << std::endl;
		// MD.write (report, A);
		// report << "Current L:" << std::endl;
		// MD.write (report, L);

		VD.firstNonzeroEntry (x, *i_A);

		F.mulin (det, x);

		F.inv (xinv, x);
		F.neg (negxinv, xinv);

		TIMER_START(ElimBelow);
		j_L = i_L;
		++j_L;

		for (j_A = i_A; ++j_A != A.rowEnd (); ++j_L) {
			if (VD.firstNonzeroEntry (a, *j_A) == (int) col) {
				// DEBUG
				// report << "Eliminating row " << j_A - A.rowBegin () << " from row " << k << std::endl;

				F.mul (negaxinv, negxinv, a);
				VD.axpyin (*j_A, negaxinv, *i_A);

				if (compute_L)
					VD.axpyin (*j_L, negaxinv, *i_L);
			}
		}
		TIMER_STOP(ElimBelow);

		++rank;

		if ((i_A - A.rowBegin ()) % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	}

	if (reduced)
		ReduceRowEchelon (A, L, compute_L, rank + start_row, start_row);

	TIMER_REPORT(GetPivot);
	TIMER_REPORT(Permute);
	TIMER_REPORT(ElimBelow);

	commentator.stop (MSG_DONE);
}

template <class Field>
void GaussJordan<Field>::RunTests () const
{
	commentator.start ("GaussJordan: Running internal tests", __FUNCTION__);
	bool pass = testFastAddinHybridVector ();
	commentator.stop (MSG_STATUS (pass));
}

} // namespace LinBox

#endif // __LINBOX_ALGORITHMS_GAUSS_JORDAN_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
