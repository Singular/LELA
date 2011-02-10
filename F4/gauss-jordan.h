/* -*- mode: c++;-*-
 * gauss-jordan.h
 * Gauss-Jordan elimination
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#ifndef __F4_GAUSS_JORDAN_H
#define __F4_GAUSS_JORDAN_H

#include <iostream>
#include <iomanip>
#include <cassert>

#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/field/gf2.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/dense-submatrix.h>
#include <linbox/matrix/sparse.h>
#include <linbox/matrix/dense-zero-one.h>
#include <linbox/matrix/sparse-zero-one.h>
#include <linbox/matrix/submatrix.h>
#include <linbox/vector/bit-subvector-word-aligned.h>
#include <linbox/vector/sparse-subvector-hybrid.h>

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

namespace LinBox {
	template <class Iterator, class ConstIterator, class Endianness>
	struct VectorTraits<std::pair<Subvector<std::vector<size_t>::iterator>, BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> > >
	{ 
		typedef std::pair<Subvector<std::vector<size_t>::iterator>, BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> > VectorType;
		typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
	};
}

namespace F4 {
	using namespace LinBox;

	template <class Field>
	class Adaptor {
	public:
		typedef typename RawVector<typename Field::Element>::Sparse SparseVector;
		typedef SparseMatrixBase<typename Field::Element, typename RawVector<typename Field::Element>::Sparse> SparseMatrix;
		typedef DenseMatrixBase<typename Field::Element> DenseMatrix;
		static const size_t cutoff = 1;
	};

	template <>
	class Adaptor<GF2> {
	public:
		typedef RawVector<bool>::Hybrid SparseVector;
		typedef SparseMatrixBase<bool, RawVector<bool>::Hybrid, VectorCategories::HybridZeroOneVectorTag> SparseMatrix;
		typedef DenseZeroOneMatrix<> DenseMatrix;
		static const size_t cutoff = __LINBOX_BITSOF_LONG;
	};

	/**
	 * \brief Implementation of Gauss-Jordan elimination which
	 * does not permute columns in the input-matrix
	 *
	 * The pivoting-strategy which it uses is to find the row with
	 * the least number of nonzero entries and a nonzero entry in
	 * the pivot-column.
	 *
	 * This class depends on the sparse-vector-facilities of LinBox, see http://www.linalg.org/
	 *
	 * This class depends on the following template-parameters:
	 *
	 * @param Field The field in which arithmetic takes
	 * place. Must satisfy the field-archetype in LinBox.
	 */
	template <class Field>
	class GaussJordan {
	public:
		typedef typename Field::Element Element;
		typedef typename MatrixDomain<Field>::Permutation Permutation;
		typedef typename MatrixDomain<Field>::Transposition Transposition;
		typedef typename Adaptor<Field>::SparseMatrix SparseMatrix;
		typedef typename Adaptor<Field>::DenseMatrix DenseMatrix;

	private:
		const Field &F;   // Field over which operations take place
		const VectorDomain<Field> VD;
		const MatrixDomain<Field> MD;

		Element zero, one;

		const size_t _cutoff;

		// Find a suitable pivot from the matrix A starting at
		// start_row. If no pivot can be found (i.e. all rows
		// from start_row onwards are already 0) return
		// -1. Otherwise fill in col with the pivot-column.
		template <class Matrix>
		int GetPivot (Matrix &A, int start_row, int &col) const
			{ return GetPivotSpecialised (A, start_row, col, typename VectorTraits<typename Matrix::Row>::VectorCategory ()); }

		// Find the first nonzero element in the given column
		// starting at the row of the same index. Return -1 if
		// none found.
		template <class Matrix>
		int GetPivotSpecialised (Matrix &A, int start_row, int &col, VectorCategories::DenseZeroOneVectorTag) const
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

		template <class Matrix>
		int GetPivotSpecialised (Matrix &A, int start_row, int &col, VectorCategories::SparseParallelVectorTag) const
		{
			typename SparseMatrix::RowIterator i;

			size_t min_nonzero = (size_t) -1, pivot = -1, k = start_row, s;
			col = A.coldim ();

		        for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
				if (!i->first.empty ()) {
					if (i->first.front () < col) {
						col = i->first.front ();
						min_nonzero = i->first.size ();
						pivot = k;
					}
					else if (i->first.front () == col && (s = i->first.size ()) < min_nonzero) {
						min_nonzero = s;
						pivot = k;
					}
				}
			}

			return pivot;
		}

		template <class Matrix>
		int GetPivotSpecialised (Matrix &A, int start_row, int &col, VectorCategories::SparseZeroOneVectorTag) const
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

		template <class Matrix>
		int GetPivotSpecialised (Matrix &A, int start_row, int &col, VectorCategories::HybridZeroOneVectorTag) const
		{
			typename SparseMatrix::RowIterator i;
			typename Field::Element a;

			size_t min_blocks = (size_t) -1, pivot = -1, k = start_row, s;
			col = A.coldim ();

		        for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
				if (!i->first.empty () && i->first.front () <= col) {
					int idx = VD.firstNonzeroEntry (a, *i);
					if (idx < col) {
						col = idx;
						min_blocks = i->first.size ();
						pivot = k;
					}
					else if (idx == col && (s = i->first.size ()) < min_blocks) {
						min_blocks = s;
						pivot = k;
					}
				}
			}

			return pivot;
		}

		// Set the given matrix to the identity-matrix
		template <class Matrix>
		void SetIdentity (Matrix &U, size_t start_row = 0) const
		{
#if 0
			StandardBasisStream<Field, typename DenseSubmatrix<typename Matrix1::Element>::Row> stream (MD.field (), res.coldim ());
			typename DenseSubmatrix<typename Matrix1::Element>::RowIterator ip = M2.rowBegin ();

			for (; ip != U.rowEnd (); ++ip)
				stream >> *ip;
#endif

			int k;

			// DEBUG
			// std::cout << __FUNCTION__ << ": U at start is " << std::endl;
			// MD.write (std::cout, U);

			MD.scal (U, false);

			// DEBUG
			// std::cout << __FUNCTION__ << ": U after clear is " << std::endl;
			// MD.write (std::cout, U);

			for (k = 0; k < U.coldim (); ++k)
				U.setEntry (k + start_row, k, one);

			// DEBUG
			// std::cout << __FUNCTION__ << ": U at end is " << std::endl;
			// MD.write (std::cout, U);
		}

		std::ostream &writePermutation (std::ostream &os, const Permutation P) const
		{
			typename Permutation::const_iterator i;

			for (i = P.begin (); i != P.end (); ++i)
				os << "(" << i->first << " " << i->second << ")";

			return os;
		}

		class CompareSecond {
		public:
			bool operator () (const Transposition &t1, const Transposition &t2) const { return t1.second < t2.second; }
		};

		// Internal recursive procedure for the kth indexed
		// Gauss-Jordan transform. Uses a divide and conquer
		// method to maximise use of fast
		// matrix-multiplication. See Chapter 2 of "Algorithms
		// for Matrix Canonical Forms", Ph.D thesis by Arne
		// Storjohann.
		void GaussTransform (SubmatrixBase<DenseMatrix> &A,
				     int                         k,
				     Element                    &d_0,
				     DenseMatrix                &U,
				     Permutation                &P,
				     SubmatrixBase<DenseMatrix> &R,
				     size_t                     &r,
				     int                        &h,
				     Element                    &d,
			             DenseMatrix                &T,
			             DenseMatrix                &Up) const
		{
			// DEBUG
			// std::cout << __FUNCTION__ << ": enter" << std::endl;
			// std::cout << __FUNCTION__ << ": A =" << std::endl;
			// MD.write (std::cout, A);
			// std::cout << __FUNCTION__ << ": U =" << std::endl;
			// MD.write (std::cout, U);
			// std::cout << __FUNCTION__ << ": k = " << k << ", d_0 = " << d_0 << std::endl;
			
			SubmatrixBase<DenseMatrix> Aw (A, k, 0, A.rowdim () - k, A.coldim ());
			
			if (MD.isZero (Aw)) {
				// DEBUG
				// std::cout << __FUNCTION__ << ": A is 0" << std::endl;

				r = 0;
				h = A.rowdim () - k;
				F.assign (d, d_0);
			}
			else if (A.coldim () <= _cutoff) {
				// DEBUG
				// std::cout << __FUNCTION__ << ": A.coldim <= " << _cutoff << ", using standard elimination" << std::endl;

				size_t l = P.size ();

				MD.copy (R, A);

				// DEBUG
				// std::cout << __FUNCTION__ << ": U before elimination:" << std::endl;
				// MD.write (std::cout, U);

				SubmatrixBase<DenseMatrix> Up (U, 0, k, U.rowdim (), U.coldim () - k);

				StandardRowEchelonForm (R, Up, P, r, d, true, true, k);

				// DEBUG
				// std::cout << __FUNCTION__ << ": U after elimination:" << std::endl;
				// MD.write (std::cout, U);

				if (!P.empty ())
					h = std::min_element (P.begin () + l, P.end (), CompareSecond ())->second;
				else
					h = k;
			}
			else {
				// DEBUG
				// std::cout << __FUNCTION__ << ": A.coldim () > " << _cutoff << std::endl;

				int m_1 = A.coldim () / 2, m_2 = A.coldim () - m_1;
				SubmatrixBase<DenseMatrix> A_1 (A, 0, 0,   A.rowdim (), m_1);
				SubmatrixBase<DenseMatrix> A_2 (A, 0, m_1, A.rowdim (), m_2);
				SubmatrixBase<DenseMatrix> R_1 (R, 0, 0,   A.rowdim (), m_1);
				
				size_t r_1, r_2;
				int h_1, h_2;
				Element d_1;
				
				GaussTransform (A_1, k, d_0, U, P, R_1, r_1, h_1, d_1, T, Up);

				// DEBUG
				// std::cout << __FUNCTION__ << ": Status after first recursive call:" << std::endl;
				// std::cout << "U =" << std::endl;
				// MD.write (std::cout, U);
				// std::cout << "P = ";
				// writePermutation (std::cout, P) << std::endl;
				// std::cout << "R =" << std::endl;
				// MD.write (std::cout, R);
				// std::cout << "r_1 = " << r_1 << std::endl;
				// std::cout << "d_1 = " << d_1 << std::endl;
				
				SubmatrixBase<DenseMatrix> U_2     (U, 0,       k,       U.rowdim (),           r_1);
				
				SubmatrixBase<DenseMatrix> P_1A_2  (T, 0,       k + m_1, T.rowdim (),           m_2);
				SubmatrixBase<DenseMatrix> P_1A_21 (T, 0,       k + m_1, k,                     m_2);
				SubmatrixBase<DenseMatrix> P_1A_22 (T, k,       k + m_1, r_1,                   m_2);
				SubmatrixBase<DenseMatrix> P_1A_23 (T, r_1 + k, k + m_1, T.rowdim () - r_1 - k, m_2);
				
				SubmatrixBase<DenseMatrix> R_2     (R, 0,       m_1,     R.rowdim (),           m_2);
				SubmatrixBase<DenseMatrix> R_21    (R, 0,       m_1,     k,                     m_2);
				SubmatrixBase<DenseMatrix> R_23    (R, k + r_1, m_1,     R.rowdim () - r_1 - k, m_2);

				MD.copy (P_1A_2, A_2);
				MD.permuteRows (P_1A_2, P.begin (), P.end ());
				MD.mul (R_2, U_2, P_1A_22);
				MD.addin (R_21, P_1A_21);
				MD.addin (R_23, P_1A_23);
				
				Permutation P_2;

				// DEBUG
				// std::cout << __FUNCTION__ << ": Status before second recursive call:" << std::endl;
				// std::cout << "B =" << std::endl;
				// MD.write (std::cout, R_2);
				
				GaussTransform (R_2, k + r_1, d_1, U, P_2, R_2, r_2, h_2, d, T, Up);

				// DEBUG
				// std::cout << __FUNCTION__ << ": Status after second recursive call:" << std::endl;
				// std::cout << "U =" << std::endl;
				// MD.write (std::cout, U);
				// std::cout << "R =" << std::endl;
				// MD.write (std::cout, R);
				// std::cout << "r_2 = " << r_2 << std::endl;
				// std::cout << "d = " << d << std::endl;
				
				SubmatrixBase<DenseMatrix> U_212        (U,  0,             k,       k + r_1,                      r_1);
				SubmatrixBase<DenseMatrix> U_3          (U,  0,             k + r_1, U.rowdim (),                  r_2);
				SubmatrixBase<DenseMatrix> P_2U_2       (U,  0,             k,       U.rowdim (),                  r_1);
				SubmatrixBase<DenseMatrix> P_2U_23      (U,  k + r_1,       k,       r_2,                          r_1);
				SubmatrixBase<DenseMatrix> P_2U_24      (U,  k + r_1 + r_2, k,       U.rowdim () - r_1 - r_2 - k,  r_1);
				SubmatrixBase<DenseMatrix> U_3P_2U_23   (Up, 0,             k,       Up.rowdim (),                 r_1);
				SubmatrixBase<DenseMatrix> U_312P_2U_23 (Up, 0,             k,       k + r_1,                      r_1);
				SubmatrixBase<DenseMatrix> U_34P_2U_23  (Up, k + r_1 + r_2, k,       Up.rowdim () - k - r_1 - r_2, r_1);

				// DEBUG
				// std::cout << "U_2 =" << std::endl;
				// MD.write (std::cout, U);
				// std::cout << "P_2 = ";
				// writePermutation (std::cout, P_2) << std::endl;

				MD.permuteRows (P_2U_2, P_2.begin (), P_2.end ());
				
				// DEBUG
				// std::cout << "P_2U_2 =" << std::endl;
				// MD.write (std::cout, P_2U_2);
				// std::cout << "P_2U_23 =" << std::endl;
				// MD.write (std::cout, P_2U_23);
				// std::cout << "U_3 =" << std::endl;;
				// MD.write (std::cout, U_3);
				
				MD.mul (U_3P_2U_23, U_3, P_2U_23);

				// DEBUG
				// std::cout << "U_3P_2U_23 =" << std::endl;
				// MD.write (std::cout, U_3P_2U_23);
				
				MD.addin (U_312P_2U_23, U_212);
				MD.addin (U_34P_2U_23, P_2U_24);
				
				MD.copy (U_2, U_3P_2U_23);
				
				P.insert (P.end (), P_2.begin (), P_2.end ());
				r = r_1 + r_2;
				h = (h_1 < h_2) ? h_1 : h_2;
			}

			// DEBUG
			// std::cout << __FUNCTION__ << ": Status at end:" << std::endl;
			// std::cout << "U =" << std::endl;
			// MD.write (std::cout, U);
			// std::cout << "P = ";
			// writePermutation (std::cout, P) << std::endl;
			// std::cout << "R =" << std::endl;
			// MD.write (std::cout, R);
			// std::cout << "k = " << k << ", r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
		}

		// Optimised version of VD.addin which can take
		// advantage of knowledge of where in v the entries of
		// w start
		template <class Vector>
		Vector &FastAddin (Vector &v, const Vector &w, size_t idx) const
			{ return FastAddinSpecialised (v, w, idx, typename VectorTraits<Vector>::VectorCategory ()); }

		template <class Vector>
		Vector &FastAddinSpecialised (Vector &v, const Vector &w, size_t idx,
					      VectorCategories::SparseZeroOneVectorTag) const
		{
			static Vector tmp;

			VD.add (tmp, Subvector<typename Vector::iterator> (v.begin () + idx, v.end ()), w);
			v.resize (idx + tmp.size ());
			std::copy (tmp.begin (), tmp.end (), v.begin () + idx);
		}

		template <class Vector>
		Vector &FastAddinSpecialised (Vector &v, const Vector &w, size_t idx,
					      VectorCategories::HybridZeroOneVectorTag) const
		{
			static Vector tmp;
			typedef BitSubvectorWordAligned<typename Vector::second_type::word_iterator, typename Vector::second_type::const_word_iterator,
				LittleEndian<typename Vector::second_type::word_iterator::value_type> > elt_type;

			std::pair<Subvector<typename Vector::first_type::iterator>, elt_type>
				t (Subvector<typename Vector::first_type::iterator> (v.first.begin () + idx, v.first.end ()),
				   elt_type (v.second.wordBegin () + idx, v.second.wordEnd ()));

			VD.add (tmp, t, w);
			v.first.resize (idx + tmp.first.size ());
			v.second.resize ((idx + tmp.second.word_size ()) << __LINBOX_LOGOF_SIZE);
			std::copy (tmp.first.begin (), tmp.first.end (), v.first.begin () + idx);
			std::copy (tmp.second.wordBegin (), tmp.second.wordEnd (), v.second.wordBegin () + idx);
		}

		void testFastAddinHybridVector () const
		{
			std::cout << "GaussJordan: Testing FastAddin" << std::endl;

			RawVector<bool>::Hybrid v, w;

			std::cout << std::hex << std::setfill ('0');

			v.first.push_back (0);
			v.second.push_word_back (0xffff0000ffff0000ULL);
			v.first.push_back (__LINBOX_BITSOF_LONG);
			v.second.push_word_back (0xffff0000ffff0000ULL);

			w.first.push_back (__LINBOX_BITSOF_LONG);
			w.second.push_word_back (0x00ffff0000ffff00ULL);

			FastAddin (v, w, 1);

			if (v.second.front_word () != 0xffff0000ffff0000ULL)
				std::cout << "Test 1 not okay: first word is " << std::setw (__LINBOX_BITSOF_LONG / 4) << v.second.front_word () << " but should be ffff0000ffff0000" << std::endl;
			else
				std::cout << "Test 1 okay" << std::endl;

			if (*(v.second.wordBegin () + 1) != 0xff00ff00ff00ff00ULL)
				std::cout << "Test 2 not okay: second word is " << std::setw (__LINBOX_BITSOF_LONG / 4) << *(v.second.wordBegin () + 1) << " but should be ff00ff00ff00ff00" << std::endl;
			else
				std::cout << "Test 2 okay" << std::endl;

			v.first.clear ();
			v.second.clear ();
			w.first.clear ();
			w.second.clear ();

			v.first.push_back (0);
			v.second.push_word_back (0xffff0000ffff0000ULL);
			v.first.push_back (__LINBOX_BITSOF_LONG);
			v.second.push_word_back (0xffff0000ffff0000ULL);

			w.first.push_back (0);
			w.second.push_word_back (0x00ffff0000ffff00ULL);

			FastAddin (v, w, 0);

			if (v.second.front_word () != 0xff00ff00ff00ff00ULL)
				std::cout << "Test 3 not okay: first word is " << std::setw (__LINBOX_BITSOF_LONG / 4) << v.second.front_word () << " but should be ff00ff00ff00ff00" << std::endl;
			else
				std::cout << "Test 3 okay" << std::endl;

			if (*(v.second.wordBegin () + 1) != 0xffff0000ffff0000ULL)
				std::cout << "Test 4 not okay: second word is " << std::setw (__LINBOX_BITSOF_LONG / 4) << *(v.second.wordBegin () + 1) << " but should be ffff0000ffff0000" << std::endl;
			else
				std::cout << "Test 4 okay" << std::endl;

			v.first.clear ();
			v.second.clear ();
			w.first.clear ();
			w.second.clear ();

			v.first.push_back (0);
			v.second.push_word_back (0xffff0000ffff0000ULL);
			v.first.push_back (__LINBOX_BITSOF_LONG);
			v.second.push_word_back (0xffff0000ffff0000ULL);

			w.first.push_back (0);
			w.second.push_word_back (0x00ffff0000ffff00ULL);
			w.first.push_back (__LINBOX_BITSOF_LONG);
			w.second.push_word_back (0x00ffff0000ffff00ULL);

			FastAddin (v, w, 0);

			if (v.second.front_word () != 0xff00ff00ff00ff00ULL)
				std::cout << "Test 5 not okay: first word is " << std::setw (__LINBOX_BITSOF_LONG / 4) << v.second.front_word () << " but should be ff00ff00ff00ff00" << std::endl;
			else
				std::cout << "Test 5 okay" << std::endl;

			if (*(v.second.wordBegin () + 1) != 0xff00ff00ff00ff00ULL)
				std::cout << "Test 6 not okay: second word is " << std::setw (__LINBOX_BITSOF_LONG / 4) << *(v.second.wordBegin () + 1) << " but should be ff00ff00ff00ff00" << std::endl;
			else
				std::cout << "Test 6 okay" << std::endl;

			std::cout << std::dec << std::setfill (' ');

			std::cout << "GaussJordan: Finished testing FastAddin" << std::endl;
		}

		template <class Matrix1, class Matrix2>
		void StandardRowEchelonFormSpecialised (Matrix1       &A,
							Matrix2       &L,
							Permutation   &P,
							size_t        &rank,
							Element       &det,
							bool           reduced,
							bool           compute_L,
							size_t         start_row,
							VectorCategories::SparseParallelVectorTag) const
		{
			commentator.start ("Sparse row-echelon form", "GaussJordan::StandardRowEchelonForm", A.rowdim () / PROGRESS_STEP);

			TIMER_DECLARE(GetPivot);
			TIMER_DECLARE(Permute);
			TIMER_DECLARE(ElimBelow);

			typename Matrix1::RowIterator i_A, j_A;

			int col, k = start_row;
			Element a, xinv, negxinv, negaxinv;

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
					Transposition t (k, pivot);
					P.push_back (t);
					MD.permuteRows (A, &t, &t + 1);

					if (compute_L) {
						SubmatrixBase<Matrix2> Lp (L, 0, 0, L.rowdim (), k - start_row);
						MD.permuteRows (Lp, &t, &t + 1);
					}
				}
				TIMER_STOP(Permute);

				F.inv (xinv, i_A->second.front ());
				F.neg (negxinv, xinv);

				F.mulin (det, i_A->second.front ());

				TIMER_START(ElimBelow);
				j_L = i_L;
				++j_L;

				for (j_A = i_A; ++j_A != A.rowEnd (); ++j_L) {
					if (VD.firstNonzeroEntry (a, *j_A) == col) {
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

			commentator.stop (MSG_DONE, NULL, "GaussJordan::StandardRowEchelonForm");
		}

		template <class Matrix1, class Matrix2>
		void StandardRowEchelonFormSpecialised (Matrix1       &A,
							Matrix2       &L,
							Permutation   &P,
							size_t        &rank,
							Element       &det,
							bool           reduced,
							bool           compute_L,
							size_t         start_row,
							VectorCategories::SparseZeroOneVectorTag) const
		{
			commentator.start ("Sparse row-echelon form", "GaussJordan::StandardRowEchelonForm", A.rowdim () / PROGRESS_STEP);

			TIMER_DECLARE(GetPivot);
			TIMER_DECLARE(Permute);
			TIMER_DECLARE(ElimBelow);

			typename Matrix1::RowIterator i_A, j_A;

			int col, k = start_row;
			Element a;

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
					Transposition t (k, pivot);
					P.push_back (t);
					MD.permuteRows (A, &t, &t + 1);

					if (compute_L) {
						SubmatrixBase<Matrix2> Lp (L, 0, 0, L.rowdim (), k - start_row);
						MD.permuteRows (Lp, &t, &t + 1);
					}
				}
				TIMER_STOP(Permute);

				TIMER_START(ElimBelow);
				j_L = i_L;
				++j_L;

				for (j_A = i_A; ++j_A != A.rowEnd (); ++j_L) {
					if (VD.firstNonzeroEntry (a, *j_A) == col) {
						VD.addin (*j_A, *i_A);

						if (compute_L)
							VD.addin (*j_L, *i_L);
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

			commentator.stop (MSG_DONE, NULL, "GaussJordan::StandardRowEchelonForm");
		}

		template <class Matrix1, class Matrix2>
		void StandardRowEchelonFormSpecialised (Matrix1       &A,
							Matrix2       &L,
							Permutation   &P,
							size_t        &rank,
							Element       &det,
							bool           reduced,
							bool           compute_L,
							size_t         start_row,
							VectorCategories::HybridZeroOneVectorTag) const
		{
			commentator.start ("Sparse row-echelon form", "GaussJordan::StandardRowEchelonForm", A.rowdim () / PROGRESS_STEP);

			TIMER_DECLARE(GetPivot);
			TIMER_DECLARE(Permute);
			TIMER_DECLARE(ElimBelow);

			typename Matrix1::RowIterator i_A, j_A;

			int col, k = start_row;
			Element a;

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
					Transposition t (k, pivot);
					P.push_back (t);
					MD.permuteRows (A, &t, &t + 1);

					if (compute_L) {
						SubmatrixBase<Matrix2> Lp (L, 0, 0, L.rowdim (), k - start_row);
						MD.permuteRows (Lp, &t, &t + 1);
					}
				}
				TIMER_STOP(Permute);

				TIMER_START(ElimBelow);
				j_L = i_L;
				++j_L;

				for (j_A = i_A; ++j_A != A.rowEnd (); ++j_L) {
					if (VD.firstNonzeroEntry (a, *j_A) == col) {
						VD.addin (*j_A, *i_A);

						if (compute_L)
							VD.addin (*j_L, *i_L);
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

			commentator.stop (MSG_DONE, NULL, "GaussJordan::StandardRowEchelonForm");
		}

		template <class Matrix1, class Matrix2>
		void StandardRowEchelonFormSpecialised (Matrix1       &A,
							Matrix2       &L,
							Permutation   &P,
							size_t        &rank,
							Element       &det,
							bool           reduced,
							bool           compute_L,
							size_t         start_row,
							VectorCategories::DenseZeroOneVectorTag) const
		{
			commentator.start ("Sparse row-echelon form", "GaussJordan::StandardRowEchelonForm", A.rowdim () / PROGRESS_STEP);

			TIMER_DECLARE(GetPivot);
			TIMER_DECLARE(Permute);
			TIMER_DECLARE(ElimBelow);

			typename Matrix1::RowIterator i_A, j_A;

			int col = 0, k = start_row;
			Element a;

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
					Transposition t (k, pivot);
					P.push_back (t);
					MD.permuteRows (A, &t, &t + 1);

					if (compute_L) {
						// DEBUG
						// std::cout << __FUNCTION__ << ": Permuting " << k << " and " << pivot << std::endl;
						// std::cout << __FUNCTION__ << ": L before permutation:" << std::endl;
						// MD.write (std::cout, L);

						SubmatrixBase<DenseMatrix> Lp (L, 0, 0, L.rowdim (), k - start_row);
						MD.permuteRows (Lp, &t, &t + 1);

						// DEBUG
						// std::cout << __FUNCTION__ << ": L after permutation:" << std::endl;
						// MD.write (std::cout, L);
					}
				}
				TIMER_STOP(Permute);

				// DEBUG
				// std::cout << __FUNCTION__ << ": Row " << k << ", pivot is " << pivot << std::endl;
				// std::cout << __FUNCTION__ << ": Current A:" << std::endl;
				// MD.write (std::cout, A);
				// std::cout << __FUNCTION__ << ": Current L:" << std::endl;
				// MD.write (std::cout, L);

				TIMER_START(ElimBelow);
				j_L = i_L;
				++j_L;

				int j_idx = k + 1;

				for (j_A = i_A; ++j_A != A.rowEnd (); ++j_L, ++j_idx) {
					if (VD.firstNonzeroEntry (a, *j_A) == col) {
						// DEBUG
						// std::cout << __FUNCTION__ << ": Eliminating row " << j_idx << " from row " << k << std::endl;

						VD.addin (*j_A, *i_A);

						if (compute_L)
							VD.addin (*j_L, *i_L);
					}
				}
				TIMER_STOP(ElimBelow);

				++rank;

				if (rank % PROGRESS_STEP == PROGRESS_STEP - 1)
					commentator.progress ();
			}

			if (reduced)
				ReduceRowEchelon (A, L, compute_L, rank + start_row, start_row);

			TIMER_REPORT(GetPivot);
			TIMER_REPORT(Permute);
			TIMER_REPORT(ElimBelow);

			commentator.stop (MSG_DONE, NULL, "GaussJordan::StandardRowEchelonForm");
		}

		template <class Matrix1, class Matrix2>
		Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
						     VectorCategories::SparseParallelVectorTag) const
		{
			commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

			typename Matrix1::RowIterator i_A, j_A, i_Ae;
			typename Matrix2::RowIterator i_L, j_L;

			Element negx, xinv;

			if (compute_L)
				i_L = L.rowBegin () + (rank - 1);

			i_A = A.rowBegin () + (rank - 1);

			do {
				if (i_A - A.rowBegin () >= start_row)
					i_Ae = i_A;

				for (j_A = A.rowBegin () + (rank - 1), j_L = L.rowBegin () + (rank - 1); j_A != i_Ae; --j_A, --j_L) {
					// We must start over each time because the operations below invalidate iterators...
					std::reverse_iterator<typename SparseMatrix::Row::first_type::iterator> i_idx (i_A->first.end ());
					std::reverse_iterator<typename SparseMatrix::Row::second_type::iterator> i_elt (i_A->second.end ());

					std::reverse_iterator<typename SparseMatrix::Row::first_type::iterator> i_stop (i_A->first.begin ());

					while (i_idx != i_stop && *i_idx > j_A->first.front ()) {
						++i_idx;
						++i_elt;
					}

					if (i_idx != i_stop && *i_idx == j_A->first.front ()) {
						F.neg (negx, *i_elt);

						VD.axpyin (*i_A, negx, *j_A);

						if (compute_L)
							VD.axpyin (*i_L, negx, *j_L);
					}
				}

				F.inv (xinv, i_A->second.front ());
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

		template <class Matrix1, class Matrix2>
		Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
						     VectorCategories::SparseZeroOneVectorTag) const
		{
			commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

			typename Matrix1::RowIterator i_A, j_A, i_Ae;
			typename Matrix2::RowIterator i_L, j_L;

			if (compute_L)
				i_L = L.rowBegin () + (rank - 1);

			i_A = A.rowBegin () + (rank - 1);

			do {
				if (i_A - A.rowBegin () >= start_row)
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

		template <class Vector>
		class PivotRowCompare {
		public:
			inline bool operator () (const typename Vector::first_type::value_type &x, const Vector &v1) const {
				return x < v1.first.front ();
			}
		};

		template <class Matrix1, class Matrix2>
		Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
						      VectorCategories::HybridZeroOneVectorTag) const
		{
			return ReduceRowEchelonHybridSpecialised (A, L, compute_L, rank, start_row,
								  LittleEndian<typename std::iterator_traits<typename Matrix1::Row::second_type::const_word_iterator>::value_type> ());
		}

		template <class Matrix1, class Matrix2, class Endianness>
		Matrix1 &ReduceRowEchelonHybridSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row, Endianness) const
		{
			commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

			if (rank == 0)
				return A;

			typename Matrix1::RowIterator i_A, i_Ae, j_A;
			typename Matrix2::RowIterator i_L;

			typename Matrix1::Row::first_type::iterator i_idx;

			typename SparseMatrix::Row::second_type::word_iterator::value_type v, v1, mask;

			if (compute_L)
				i_L = L.rowBegin () + (rank - 1);

			i_A = A.rowBegin () + (rank - 1);

			int current_row = rank - 1, elim_row;

			do {
				// DEBUG
				// std::cout << __FUNCTION__ << ": Row " << current_row << ", state of A:" << std::endl;
				// MD.write (std::cout, A);

				if (i_A - A.rowBegin () >= start_row)
					i_Ae = i_A;

				size_t prev_idx = i_A->first.size () - 1;
				elim_row = rank - 1;

				j_A = A.rowBegin () + elim_row;

				while (elim_row > std::max (current_row, (int) start_row - 1)) {
					if (prev_idx >= i_A->first.size ())
						prev_idx = i_A->first.size () - 1;

					typename SparseMatrix::Row::first_type::value_type j_idx = j_A->first.front ();

					// Note: We don't need to test for validity because, if the input
					// is valid, i_idx will *never* go past the beginning
					for (i_idx = i_A->first.begin () + prev_idx; *i_idx > j_idx; --i_idx, --prev_idx);

					if (*i_idx == j_A->first.front ()) {
						v = j_A->second.front_word ();
						mask = Endianness::first_position (v);

						if (*(i_A->second.wordBegin () + prev_idx) & mask) {
							// DEBUG
							// std::cout << "ReduceRowEchelonSpecialised: eliminating row " << current_row << " from row " << elim_row << std::endl;

							FastAddin (*i_A, *j_A, prev_idx);

							if (compute_L)
								VD.addin (*i_L, *(L.rowBegin () + elim_row));
						}
					}

					v1 = *(i_A->second.wordBegin () + prev_idx);
					mask = 0;

					while (j_A->first.front () >= i_A->first[prev_idx] && !(v1 & mask)) {
						--j_A;
						--elim_row;

						if (j_A->first.front () > i_A->first[prev_idx]) {
							j_A = std::upper_bound (A.rowBegin () + current_row + 1, A.rowBegin () + elim_row, i_A->first[prev_idx],
										PivotRowCompare<typename SparseMatrix::Row> ());
							--j_A;
							elim_row = j_A - A.rowBegin ();

							if (elim_row == current_row)
								break;
						}

						v = j_A->second.front_word ();
						mask = Endianness::first_position (v);
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

		template <class Matrix1, class Matrix2>
		Matrix1 &ReduceRowEchelonSpecialised (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row,
						     VectorCategories::DenseZeroOneVectorTag) const
		{
			commentator.start ("Reducing row-echelon form", "GaussJordan::ReduceRowEchelonSpecialised", A.rowdim () / PROGRESS_STEP);

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
				// std::cout << __FUNCTION__ << ": Row " << current_row << ", current A:" << std::endl;
				// MD.write (std::cout, A);

				for (elim_row = rank - 1, j_A = A.rowBegin () + elim_row; elim_row > std::max (current_row, (int) start_row - 1); --elim_row, --j_A) {
					size_t col = VD.firstNonzeroEntry (a, *j_A);

					if ((*i_A)[col]) {
						// DEBUG
						// std::cout << __FUNCTION__ << ": Eliminating " << current_row << " from " << elim_row << std::endl;

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

	public:
		/**
		 * \brief Constructor
		 *
		 * @param _F Field over which operations take place
		 */
		GaussJordan (const Field &_F) : F (_F), VD (_F), MD (_F), _cutoff (Adaptor<Field>::cutoff) {
			F.init (zero, 0);
			F.init (one, 1);
		}

		/**
		 * \brief Convert the matrix A into reduced
		 * row-echelon form, preserving as much sparsity as
		 * possible. At the end the matrices satisfy the
		 * equation R=UPA, with R in row-echelon form.
		 *
		 * If A is invertible, then U will be the inverse of
		 * A.
		 *
		 * R and A may be the same matrix, in which case A
		 * will be replaced by its row-echelon form.
		 *
		 * @param R Sparse matrix-object into which to store
		 * row-echelon form of A; should have same dimensions
		 * as A.
		 *
		 * @param U Dense matrix-object into which to store
		 * the matrix U. Should be n x n with n equal to the
		 * row-dimension of A.
		 *
		 * @param P Permutation in which to store the
		 * row-permutations of A made by the choice of pivots.
		 *
		 * @param A Matrix to be converted to row-echelon
		 * form. Not altered.
		 *
		 * @param rank Integer into which to store the
		 * computed rank
		 *
		 * @param det Field-element into which to store the
		 * computed determinant
		 */
		void DenseRowEchelonForm (DenseMatrix &R,
					  DenseMatrix &U,
					  Permutation &P,
					  DenseMatrix &A,
					  size_t      &rank,
					  Element     &det)
		{
			commentator.start ("Dense row-echelon form", "GaussJordan::DenseRowEchelonForm");

			int h;
			
			SetIdentity (U);  // Necessary in case where A.rowdim () > A.coldim ()
			P.clear ();
			
			DenseMatrix T (A.rowdim (), A.coldim ());
			DenseMatrix Up (A.rowdim (), A.rowdim ());

			SubmatrixBase<DenseMatrix> Rp (R, 0, 0, R.rowdim (), R.coldim ());
			SubmatrixBase<DenseMatrix> Ap (A, 0, 0, A.rowdim (), A.coldim ());
			
			GaussTransform (Ap, 0, one, U, P, Rp, rank, h, det, T, Up);

			commentator.stop (MSG_DONE, NULL, "GaussJordan::DenseRowEchelonForm");
		}

		/**
		 * \brief Compute the (reduced or non-reduced)
		 * row-echelon form of a matrix using a sparse
		 * algorithm
		 *
		 * At conclusion, the parameters will have the
		 * property that A_out=LPA_in, where A_out is the
		 * matrix A at output and A_in is the matrix A at
		 * input. R is in reduced row-echelon form, L is lower
		 * triangular, and P is a permutation.
		 *
		 * In comparison with @see RowEchelonForm, this
		 * version does not take advantage of fast
		 * matrix-multiplication and does not use a
		 * divide-and-conquer method. It also modifies the
		 * input-matrix A.
		 *
		 * The pivot-strategy is to find the row with the
		 * fewest elements. This seems to be the only sensible
		 * approach, given that we are not allowed to permute
		 * columns.
		 *
		 * @param A The sparse matrix whose reduced
		 * row-echelon form to compute. Will be replaced by
		 * its reduced row-echelon form during computation.
		 *
		 * @param L The dense matrix into which to store the
		 * matrix L as defined above. Should be of size n x n,
		 * with n the row-dimension of A.
		 *
		 * @param P The permutation into which to store the
		 * permutation P as defined above.
		 *
		 * @param rank An integer into which to store the
		 * computed rank of A.
		 *
		 * @param det A field-element into which to store the
		 * computed determinant of the submatrix of A formed
		 * by taking pivot-rows and -columns.
		 *
		 * @param reduced True if the routine should compute
		 * the reduced row-echelon form and false if it should
		 * only compute the (non-reduced) row-echelon form
		 *
		 * @param compute_L True if the matrix L should be
		 * computed. If false, then L is ignored.
		 *
		 * @param start_row Start at this row. Intended for
		 * internal use.
		 */
		template <class Matrix1, class Matrix2>
		void StandardRowEchelonForm (Matrix1       &A,
					     Matrix2       &L,
					     Permutation   &P,
					     size_t        &rank,
					     Element       &det,
					     bool           reduced = false,
					     bool           compute_L = true,
					     size_t         start_row = 0) const
			{ StandardRowEchelonFormSpecialised (A, L, P, rank, det, reduced, compute_L, start_row,
							     typename VectorTraits<typename Matrix1::Row>::VectorCategory ()); }

		/** \brief Take a matrix of known rank in row-echelon
		 * form and convert to reduced row-echelon form
		 *
		 * @param A Input matrix A in row-echelon form;
		 * replaced by its reduced row-echelon form
		 *
		 * @param L Dense matrix L into which to store
		 * conversion-matrix
		 *
		 * @param compute_L bool, true if L should be
		 * computed; if false then L is left unchanged
		 *
		 * @param rank Rank of A; must be known a priori
		 * (though can be easily computed by scanning rows of
		 * A)
		 *
		 * @param start_row Pivot-rows begin at this
		 * row. Intended for internal use only.
		 *
		 * @return Reference to A
		 */
		template <class Matrix1, class Matrix2>
		Matrix1 &ReduceRowEchelon (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row = 0) const
			{ return ReduceRowEchelonSpecialised (A, L, compute_L, rank, start_row,
							      typename VectorTraits<typename Matrix1::Row>::VectorCategory ()); }

		/** Run internal tests
		 */
		void RunTests () const
		{
			std::cout << "GaussJordan: Running internal tests..." << std::endl;
			testFastAddinHybridVector ();
			std::cout << "GaussJordan: Finished running internal tests." << std::endl;
		}
	};
}

#endif // __F4_GAUSS_JORDAN_H
