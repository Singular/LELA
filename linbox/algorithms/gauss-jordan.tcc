/* linbox/algorithms/gauss-jordan.tcc
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

#ifndef __LINBOX_ALGORITHMS_GAUSS_JORDAN_TCC
#define __LINBOX_ALGORITHMS_GAUSS_JORDAN_TCC

#include "linbox/algorithms/gauss-jordan.h"

#include "linbox/vector/subvector.h"
#include "linbox/blas/level1.h"
#include "linbox/blas/level2.h"
#include "linbox/blas/level3.h"

#include "linbox/matrix/sparse-zero-one.h"

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

template <class Field, class Modules>
template <class Matrix>
int GaussJordan<Field, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						      VectorCategories::DenseVectorTag) const
{
	typename Matrix::ConstRowIterator i;

	for (; col < A.coldim (); ++col) {
		int k = start_row;

		for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k)
			if (!ctx.F.isZero ((*i)[col]))
				return k;
	}

	return -1;
}

template <class Field, class Modules>
template <class Matrix>
int GaussJordan<Field, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
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

template <class Field, class Modules>
template <class Matrix>
int GaussJordan<Field, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						      VectorCategories::SparseVectorTag) const
{
	typename SparseMatrix::ConstRowIterator i;

	size_t min_nonzero = (size_t) -1, pivot = -1, k = start_row;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty ()) {
			if (i->front ().first < col) {
				col = i->front ().first;
				min_nonzero = i->size ();
				pivot = k;
			}
			else if (i->front ().first == col && i->size () < min_nonzero) {
				min_nonzero = i->size ();
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Field, class Modules>
template <class Matrix>
int GaussJordan<Field, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
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

template <class Field, class Modules>
template <class Matrix>
int GaussJordan<Field, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						      VectorCategories::HybridZeroOneVectorTag) const
{
	typename SparseMatrix::ConstRowIterator i;
	typename Field::Element a;

	size_t min_blocks = (size_t) -1, pivot = -1, k = start_row, s;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty () && i->front ().first << WordTraits<typename Matrix::Row::word_type>::logof_size <= (int) col) {
			int idx = BLAS1::head (ctx, a, *i);
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

template <class Field, class Modules>
template <class Matrix>
void GaussJordan<Field, Modules>::SetIdentity (Matrix &U, size_t start_row) const
{
	StandardBasisStream<Field, typename Matrix::Row> stream (ctx.F, U.coldim ());
	typename Matrix::RowIterator i = U.rowBegin () + start_row;

	for (; i != U.rowEnd (); ++i)
		stream >> *i;
}

template <class Field, class Modules>
void GaussJordan<Field, Modules>::GaussJordanTransform (DenseMatrix  &A,
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
	// BLAS3::write (ctx, report, A);
	// report << "U =" << std::endl;
	// BLAS3::write (ctx, report, U);
	// report << "k = " << k << ", d_0 = " << d_0 << std::endl;

	DenseMatrix Aw (A, k, 0, A.rowdim () - k, A.coldim ());

	if (BLAS3::is_zero (ctx, Aw)) {
		// DEBUG
		// report << "A is 0" << std::endl;

		r = 0;
		h = A.rowdim () - k;
		ctx.F.assign (d, d_0);
	}
	else if (A.coldim () <= _cutoff) {
		// DEBUG
		// report << "A.coldim <= " << _cutoff << ", using standard elimination" << std::endl;

		size_t l = P.size ();

		// DEBUG
		// report << "U before elimination:" << std::endl;
		// BLAS3::write (ctx, report, U);

		Submatrix<DenseMatrix> Up (U, 0, k, U.rowdim (), U.coldim () - k);

		StandardRowEchelonForm (A, Up, P, r, d, true, true, k);

		// DEBUG
		// report << "U after elimination:" << std::endl;
		// BLAS3::write (ctx, report, U);

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
		// BLAS3::write (ctx, report, U);
		// report << "P = ";
		// BLAS3::write_permutation (report, P) << std::endl;
		// report << "R =" << std::endl;
		// BLAS3::write (ctx, report, R);
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

		BLAS3::permute_rows (ctx, P.begin (), P.end (), A_2);

		BLAS3::copy (ctx, A_22, T);
		BLAS3::gemm (ctx, ctx.F.one (), U_21, T, ctx.F.one (),  A_21);
		BLAS3::gemm (ctx, ctx.F.one (), U_22, T, ctx.F.zero (), A_22);
		BLAS3::gemm (ctx, ctx.F.one (), U_23, T, ctx.F.one (),  A_23);

		Permutation P_2;

		// DEBUG
		// report << "Status before second recursive call:" << std::endl;
		// report << "B =" << std::endl;
		// BLAS3::write (ctx, report, R_2);

		GaussJordanTransform (A_2, k + r_1, d_1, U, P_2, r_2, h_2, d, S, T);

		// DEBUG
		// report << "Status after second recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "R =" << std::endl;
		// BLAS3::write (ctx, report, R);
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
		// BLAS3::write (ctx, report, U);
		// report << "P_2 = ";
		// BLAS3::write_permutation (report, P_2) << std::endl;

		BLAS3::permute_rows (ctx, P_2.begin (), P_2.end (), P_2U_2);
				
		// DEBUG
		// report << "P_2U_2 =" << std::endl;
		// BLAS3::write (ctx, report, P_2U_2);
		// report << "P_2U_23 =" << std::endl;
		// BLAS3::write (ctx, report, P_2U_23);
		// report << "U_3 =" << std::endl;;
		// BLAS3::write (ctx, report, U_3);

		S.resize (r_2, r_1);
		BLAS3::copy (ctx, P_2U_23, S);

		BLAS3::gemm (ctx, ctx.F.one (), U_312, S, ctx.F.one (), U_212);
		BLAS3::gemm (ctx, ctx.F.one (), U_33, S, ctx.F.zero (), P_2U_23);
		BLAS3::gemm (ctx, ctx.F.one (), U_34, S, ctx.F.one (), P_2U_24);

		// DEBUG
		// report << "U_3P_2U_23 =" << std::endl;
		// BLAS3::write (ctx, report, U_3P_2U_23);

		P.insert (P.end (), P_2.begin (), P_2.end ());
		r = r_1 + r_2;
		h = std::min (h_1, h_2);
	}

	// DEBUG
	// report << "Status at end:" << std::endl;
	// report << "U =" << std::endl;
	// BLAS3::write (ctx, report, U);
	// report << "P = ";
	// BLAS3::write_permutation (report, P) << std::endl;
	// report << "R =" << std::endl;
	// BLAS3::write (ctx, report, R);
	// report << "k = " << k << ", r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
}

template <class Field, class Modules>
void GaussJordan<Field, Modules>::GaussTransform (DenseMatrix             &A,
						  Element                  d_0,
						  Submatrix<DenseMatrix>  &U,
						  Permutation             &P,
						  size_t                  &r,
						  int                     &h,
						  Element                 &d) const
{
	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	// DEBUG
	// report << "enter" << std::endl;
	// report << "A =" << std::endl;
	// BLAS3::write (ctx, report, A);
	// report << "U =" << std::endl;
	// BLAS3::write (ctx, report, U);

	if (BLAS3::is_zero (ctx, A)) {
		// DEBUG
		// report << "A is 0" << std::endl;

		r = 0;
		h = A.rowdim ();
		ctx.F.assign (d, d_0);
	}
	else if (A.coldim () <= _cutoff) {
		// DEBUG
		// report << "A.coldim <= " << _cutoff << ", using standard elimination" << std::endl;

		size_t l = P.size ();

		// DEBUG
		// report << "U before elimination:" << std::endl;
		// BLAS3::write (ctx, report, U);

		StandardRowEchelonForm (A, U, P, r, d, false, true);

		// DEBUG
		// report << "U after elimination:" << std::endl;
		// BLAS3::write (ctx, report, U);

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

		GaussTransform (A_1, d_0, U, P, r_1, h_1, d_1);

		// DEBUG
		// report << "Status after first recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "P = ";
		// BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;
		// report << "A =" << std::endl;
		// BLAS3::write (ctx, report, A);
		// report << "r_1 = " << r_1 << std::endl;
		// report << "d_1 = " << d_1 << std::endl;

		Submatrix<DenseMatrix> U_11 (U, 0,   0,   r_1,               r_1);
		Submatrix<DenseMatrix> U_12 (U, r_1, 0,   U.rowdim () - r_1, r_1);
		Submatrix<DenseMatrix> U_22 (U, r_1, r_1, U.rowdim () - r_1, U.coldim () - r_1);

		DenseMatrix A_2  (A, 0,   m_1, A.rowdim (),       m_2);
		DenseMatrix A_21 (A, 0,   m_1, r_1,               m_2);
		DenseMatrix A_22 (A, r_1, m_1, A.rowdim () - r_1, m_2);

		BLAS3::permute_rows (ctx, P.begin (), P.end (), A_2);
		BLAS3::gemm (ctx, ctx.F.one (), U_12, A_21, ctx.F.one (), A_22);
		BLAS3::trmm (ctx, ctx.F.one (), U_11, A_21, LowerTriangular, true);

		Permutation P_2;

		// DEBUG
		// report << "Status before second recursive call:" << std::endl;
		// report << "A_22 =" << std::endl;
		// BLAS3::write (ctx, report, A_22);

		GaussTransform (A_22, d_1, U_22, P_2, r_2, h_2, d);

		// DEBUG
		// report << "Status after second recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "R =" << std::endl;
		// BLAS3::write (ctx, report, R);
		// report << "r_2 = " << r_2 << std::endl;
		// report << "d = " << d << std::endl;

		// DEBUG
		// report << "U_2 =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "P_2 = ";
		// BLAS1::write_permutation (report, P_2.begin (), P_2.end ()) << std::endl;

		BLAS3::permute_rows (ctx, P_2.begin (), P_2.end (), U_12);
		BLAS3::trmm (ctx, ctx.F.one (), U_22, U_12, LowerTriangular, true);
				
		// DEBUG
		// report << "P_2U_2 =" << std::endl;
		// BLAS3::write (ctx, report, P_2U_2);

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
	// BLAS3::write (ctx, report, U);
	// report << "P = ";
	// BLAS3::write_permutation (report, P) << std::endl;
	// report << "R =" << std::endl;
	// BLAS3::write (ctx, report, R);
	// report << "k = " << k << ", r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
}

template <>
class MutableSubvector<Vector<GF2>::Hybrid> : public MutableSubvector<SparseVector<Vector<GF2>::Hybrid::word_type, std::vector<Vector<GF2>::Hybrid::index_type>, std::vector<Vector<GF2>::Hybrid::word_type> > >
{
public:
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
	typedef Vector<GF2>::Hybrid::Endianness Endianness;
	typedef Vector<GF2>::Hybrid::index_type index_type;
	typedef Vector<GF2>::Hybrid::word_type word_type;

	MutableSubvector (Vector<GF2>::Hybrid &v, Vector<GF2>::Hybrid::iterator begin, Vector<GF2>::Hybrid::iterator end)
		: MutableSubvector<SparseVector<Vector<GF2>::Hybrid::word_type, std::vector<Vector<GF2>::Hybrid::index_type>, std::vector<Vector<GF2>::Hybrid::word_type> > > (v, begin, end) {}
};

template <class Field, class Modules>
template <class Vector>
Vector &GaussJordan<Field, Modules>::FastAxpy (Vector &v, const typename Field::Element &a, const Vector &w, size_t idx) const
{
	MutableSubvector<Vector> t (v, v.begin () + idx, v.end ());
	BLAS1::axpy (ctx, a, w, t);
	return v;
}

template <class Field, class Modules>
bool GaussJordan<Field, Modules>::testFastAxpyHybridVector () const
{
	commentator.start ("Testing FastAxpy", __FUNCTION__);

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	bool pass = true;

	Vector<GF2>::Hybrid v, w;

	error << std::hex << std::setfill ('0');

	v.push_back (Vector<GF2>::Hybrid::value_type (0, 0xffff0000ffff0000ULL));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, 0xffff0000ffff0000ULL));

	w.push_back (Vector<GF2>::Hybrid::value_type (1, 0x00ffff0000ffff00ULL));

	FastAxpy (v, ctx.F.one (), w, 1);

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

	FastAxpy (v, ctx.F.one (), w, 0);

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

	FastAxpy (v, ctx.F.one (), w, 0);

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

template <class Field, class Modules>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field, Modules>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
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
		// BLAS3::write (ctx, report, A);

		for (elim_row = rank - 1, j_A = A.rowBegin () + elim_row; elim_row > std::max (current_row, (int) start_row - 1); --elim_row, --j_A) {
			size_t col = BLAS1::head (ctx, a, *j_A);

			if (!ctx.F.isZero ((*i_A)[col])) {
				// DEBUG
				// report << "Eliminating " << current_row << " from " << elim_row << std::endl;

				BLAS1::axpy (ctx, ctx.F.one (), *j_A, *i_A);

				if (compute_L)
					BLAS1::axpy (ctx, ctx.F.one (), *(L.rowBegin () + elim_row), *i_L);
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

template <class Field, class Modules>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field, Modules>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
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
				ctx.F.neg (negx, i->second);

				BLAS1::axpy (ctx, negx, *j_A, *i_A);

				if (compute_L)
					BLAS1::axpy (ctx, negx, *j_L, *i_L);
			}
		}

		ctx.F.inv (xinv, i_A->front ().second);
		BLAS1::scal (ctx, xinv, *i_A);

		if (compute_L) {
			BLAS1::scal (ctx, xinv, *i_L);
			--i_L;
		}

		if ((A.rowEnd () - i_A) % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "GaussJordan::ReduceRowEchelonSpecialised");

	return A;
}

template <class Field, class Modules>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field, Modules>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
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
				FastAxpy (*i_A, ctx.F.one (), *j_A, prev_idx - 1);

				if (compute_L)
					BLAS1::axpy (ctx, ctx.F.one (), *j_L, *i_L);
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

template <class Field, class Modules>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field, Modules>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
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
		// BLAS3::write (ctx, report, A);

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

					FastAxpy (*i_A, ctx.F.one (), *j_A, prev_idx);

					if (compute_L)
						BLAS1::axpy (ctx, ctx.F.one (), *(L.rowBegin () + elim_row), *i_L);
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

template <class Field, class Modules>
template <class Matrix1, class Matrix2>
Matrix1 &GaussJordan<Field, Modules>::ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
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
		// BLAS3::write (ctx, report, A);

		for (elim_row = rank - 1, j_A = A.rowBegin () + elim_row; elim_row > std::max (current_row, (int) start_row - 1); --elim_row, --j_A) {
			size_t col = BLAS1::head (ctx, a, *j_A);

			if ((*i_A)[col]) {
				// DEBUG
				// report << "Eliminating " << current_row << " from " << elim_row << std::endl;

				BLAS1::axpy (ctx, ctx.F.one (), *j_A, *i_A);

				if (compute_L)
					BLAS1::axpy (ctx, ctx.F.one (), *(L.rowBegin () + elim_row), *i_L);
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

template <class Field, class Modules>
void GaussJordan<Field, Modules>::DenseRowEchelonForm (DenseMatrix &A,
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

		GaussJordanTransform (A, 0, ctx.F.one (), U, P, rank, h, det, S, T);
	} else {
		Submatrix<DenseMatrix> Up (U, 0, 0, U.rowdim (), U.coldim ());

		GaussTransform (A, ctx.F.one (), Up, P, rank, h, det);
	}

	commentator.stop (MSG_DONE, NULL, "GaussJordan::DenseRowEchelonForm");
}

template <class Field, class Modules>
template <class Matrix1, class Matrix2>
void GaussJordan<Field, Modules>::StandardRowEchelonForm (Matrix1       &A,
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

	ctx.F.assign (a, ctx.F.zero ());
	ctx.F.assign (x, ctx.F.zero ());

	typename Matrix2::RowIterator i_L, j_L;

	if (compute_L)
		SetIdentity (L, start_row);

	P.clear ();
	rank = 0;
	ctx.F.init (det, 1);

	for (i_A = A.rowBegin () + k, i_L = L.rowBegin () + k; i_A != A.rowEnd (); ++k, ++i_A, ++i_L) {
		TIMER_START(GetPivot);
		int pivot = GetPivot (A, k, col);
		TIMER_STOP(GetPivot);

		if (pivot == -1)
			break;

		TIMER_START(Permute);
		if (k != pivot) {
			// DEBUG
			// report << "Permuting " << k << " and " << pivot << " (first " << k - start_row << " columns)" << std::endl;
			// report << "L before permutation:" << std::endl;
			// BLAS3::write (ctx, report, L);

			Transposition t (k, pivot);
			P.push_back (t);
			BLAS3::permute_rows (ctx, &t, &t + 1, A);

			if (compute_L) {
				Submatrix<typename RealMatrixType<Matrix2>::Type> Lp (L, 0, 0, L.rowdim (), k - start_row);
				BLAS3::permute_rows (ctx, &t, &t + 1, Lp);
			}

			// DEBUG
			// report << "L after permutation:" << std::endl;
			// BLAS3::write (ctx, report, L);
		}
		TIMER_STOP(Permute);

		// DEBUG
		// report << "Row " << k << ", pivot is " << pivot << std::endl;
		// report << "Current A:" << std::endl;
		// BLAS3::write (ctx, report, A);
		// report << "Current L:" << std::endl;
		// BLAS3::write (ctx, report, L);

		BLAS1::head (ctx, x, *i_A);

		ctx.F.mulin (det, x);

		ctx.F.inv (xinv, x);
		ctx.F.neg (negxinv, xinv);

		TIMER_START(ElimBelow);
		j_L = i_L;
		++j_L;

		for (j_A = i_A; ++j_A != A.rowEnd (); ++j_L) {
			if (BLAS1::head (ctx, a, *j_A) == (int) col) {
				// DEBUG
				// report << "Eliminating row " << j_A - A.rowBegin () << " from row " << k << std::endl;

				ctx.F.mul (negaxinv, negxinv, a);
				BLAS1::axpy (ctx, negaxinv, *i_A, *j_A);

				if (compute_L)
					BLAS1::axpy (ctx, negaxinv, *i_L, *j_L);
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

template <class Field, class Modules>
void GaussJordan<Field, Modules>::RunTests () const
{
	commentator.start ("GaussJordan: Running internal tests", __FUNCTION__);
	bool pass = testFastAxpyHybridVector ();
	commentator.stop (MSG_STATUS (pass));
}

} // namespace LinBox

#endif // __LINBOX_ALGORITHMS_GAUSS_JORDAN_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
