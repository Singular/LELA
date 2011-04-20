/* -*- mode: c++;-*-
 * solver.h
 * Linear system-solver for F4-implementation
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#ifndef __F4_SOLVER_H
#define __F4_SOLVER_H

#include <linbox/vector/vector-domain.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/algorithms/gauss-jordan.h>

#ifndef PROGRESS_STEP
#  define PROGRESS_STEP 1024
#endif // PROGRESS_STEP

namespace F4 {
	using namespace LinBox;

	/**
	 * \brief Implementation of algorithm for solving linear systems generated in F4
	 *
	 * The algorithm used is based loosely on the algorithm
	 * described in J.C. Faugère's paper "Parallel Gaussian
	 * Elimination for Gröbner bases computations in finite
	 * fields", PASCO 2010.
	 *
	 * This implementation uses LinBox for underlying matrix-arithmetic; see http://www.linalg.org/
	 */
        template <class Field>
	class F4Solver {
	public:
		typedef typename Field::Element Element;
		typedef typename GaussJordan<Field>::SparseMatrix SparseMatrix;
		typedef typename GaussJordan<Field>::DenseMatrix DenseMatrix;
		typedef typename Adaptor<Field>::Endianness Endianness;

	private:
		const Field &F;
		const VectorDomain<Field> VD;
		const MatrixDomain<Field> MD;
		GaussJordan<Field> GJ;

		Element one, neg_one;

		const double threshold;

		struct Block {
			size_t col_idx;
			size_t size;
			size_t D_size;

			Block (size_t _col_idx, size_t _size, size_t _D_size) : col_idx (_col_idx), size (_size), D_size (_D_size) {}
		};

		size_t GetNextBlockColumn (size_t coldim, typename std::vector<Block>::const_iterator &b, typename std::vector<Block>::const_iterator &b_end) const
		{
			return b->col_idx + b->size + b->D_size;
		}

		template <class Matrix>
		void FindPivotRows (const Matrix &A, std::vector<Block> &blocks, std::vector<size_t> &pivot_rows) const
		{
			commentator.start ("Finding pivot-rows", __FUNCTION__);

			typename Matrix::ConstRowIterator i_A;
			int last_col = -1, col, first_col_in_block = 0, height = 0;
			typename Field::Element a;

			for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A) {
				col = VD.firstNonzeroEntry (a, *i_A);

				if (col == -1)
					break;
				else if (col == last_col + 1) {
					pivot_rows.push_back (i_A - A.rowBegin ());
					last_col = col;
					++height;
				}
				else if (col > last_col + 1) {
					blocks.push_back (Block (first_col_in_block, height, col - last_col - 1));
					pivot_rows.push_back (i_A - A.rowBegin ());
					height = 1;
					first_col_in_block = last_col = col;
				}
			}

			if (col != -1)
				blocks.push_back (Block (first_col_in_block, height, A.coldim () - first_col_in_block - height));

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}

		void AssemblePivotRows (const SparseMatrix &A, SparseMatrix &Agj, SparseMatrix &B, const std::vector<size_t> &pivot_rows) const
		{
			commentator.start ("Assembly of pivot-rows", __FUNCTION__);

			typename SparseMatrix::ConstRowIterator i_A;
			typename SparseMatrix::RowIterator i_Agj = Agj.rowBegin (), i_B = B.rowBegin ();
			std::vector<size_t>::const_iterator i = pivot_rows.begin ();

			for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A) {
				if (i != pivot_rows.end () && i_A - A.rowBegin () == *i) {
					VD.copy (*i_Agj, *i_A);
					++i_Agj;
					++i;
				} else {
					VD.copy (*i_B, *i_A);
					++i_B;
				}
			}

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}

		template <class Matrix1, class Matrix2>
		void SplitMatrix (const Matrix1 &X, Matrix2 &A, Matrix2 &B, const std::vector<Block> &blocks) const
		{
			commentator.start ("Splitting input-matrix into parts", __FUNCTION__);

			typename std::vector<Block>::const_iterator i;

			size_t prev_col_A = 0, prev_col_B = 0;

			for (i = blocks.begin (); i != blocks.end (); ++i) {
				Submatrix<const Matrix1> A_part_source (X, 0, i->col_idx, X.rowdim (), i->size);
				Submatrix<Matrix2> A_part_dest (A, 0, prev_col_A, A.rowdim (), i->size);

				Submatrix<const Matrix1> B_part_source (X, 0, i->col_idx + i->size, X.rowdim (), i->D_size);
				Submatrix<Matrix2> B_part_dest (B, 0, prev_col_B, B.rowdim (), i->D_size);

				MD.copy (A_part_dest, A_part_source);
				MD.copy (B_part_dest, B_part_source);

				prev_col_A += i->size;
				prev_col_B += i->D_size;
			}

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}

		std::vector<Block> &CombineBlocks (std::vector<Block> &output, const std::vector<Block> &B_blocks, const std::vector<Block> &D_blocks) const
		{
			typename std::vector<Block>::const_iterator i_B, i_D = D_blocks.begin ();
			size_t D_col = 0;

			for (i_B = B_blocks.begin (); i_B != B_blocks.end (); ++i_B) {
				while (i_D != D_blocks.end () && i_D->col_idx + i_D->size < D_col + i_B->D_size) {
					output.push_back (Block (i_B->col_idx + i_D->col_idx - D_col, i_D->size, i_D->D_size));
					++i_D;
				}

				D_col += i_B->D_size;
			}
		}

		size_t round_up (size_t v, size_t x) const
			{ return ((v + x - 1) / x) * x; }

		template <class Vector1, class Vector2>
		Vector1 &CopyBlockSpecialised (Vector1 &dest, const Vector2 &source, size_t source_idx, size_t len, size_t dest_idx,
					       VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag) const
		{
			typedef BitVector<typename Vector1::Endianness> TmpVector;
			typedef BitSubvector<typename TmpVector::iterator, typename TmpVector::const_iterator> TmpSubvector;
			typedef BitSubvector<typename Vector2::const_iterator, typename Vector2::const_iterator> ConstTmpSubvector;

			size_t dest_idx_offset = dest_idx & WordTraits<typename Vector2::word_type>::pos_mask;

			TmpVector tmp (round_up (len + dest_idx_offset, WordTraits<typename Vector2::word_type>::bits));
			ConstTmpSubvector v1 (source.begin () + source_idx, source.begin () + (source_idx + len));
			TmpSubvector v2 (tmp.begin () + dest_idx_offset, tmp.begin () + (dest_idx_offset + len));

			VD.copy (v2, v1);

			typename TmpVector::word_iterator i = tmp.wordBegin ();
			size_t idx = dest_idx >> WordTraits<typename Vector1::word_type>::logof_size;

			if (!dest.empty () && dest.back ().first == idx) {
				dest.back ().second |= *i;
				++i;
				++idx;
			}

			for (; i != tmp.wordEnd (); ++i, ++idx)
				dest.push_back (typename Vector1::value_type (idx, *i));

			return dest;
		}

		template <class Vector1, class Vector2>
		Vector1 &CopyBlock (Vector1 &dest, const Vector2 &source, size_t source_idx, size_t len, size_t dest_idx) const
			{ return CopyBlockSpecialised (dest, source, source_idx, len, dest_idx,
						       typename VectorTraits<Field, Vector1>::VectorCategory (),
						       typename VectorTraits<Field, Vector2>::VectorCategory ()); }

		void AssembleOutput (SparseMatrix &R, const DenseMatrix &B, const DenseMatrix &D, const std::vector<Block> &blocks, const std::vector<size_t> &pivot_rows) const
		{
			commentator.start ("Assembly of final output", __FUNCTION__, R.rowdim () / 1024);

			typename SparseMatrix::RowIterator i_R;
			typename DenseMatrix::ConstRowIterator i_B = B.rowBegin (), i_D = D.rowBegin ();
			typename std::vector<Block>::const_iterator block, block_begin = blocks.begin ();
			typename Field::Element a;

			std::vector<size_t>::const_iterator pr = pivot_rows.begin ();

			size_t idx = 0, col_idx;

			for (i_R = R.rowBegin (); i_R != R.rowEnd () && block_begin != blocks.end (); ++i_R) {
				i_R->push_back (typename SparseMatrix::Row::value_type
						(idx >> WordTraits<typename SparseMatrix::Row::word_type>::logof_size,
						 SparseMatrix::Row::Endianness::e_j (idx & WordTraits<typename SparseMatrix::Row::word_type>::pos_mask)));

				for (block = block_begin, col_idx = 0; block != blocks.end (); ++block) {
					if (pr != pivot_rows.end () && *pr == idx) {
						CopyBlock (*i_R, *i_B, col_idx, block->D_size, block->col_idx + block->size);
						++i_B; ++pr;
					} else {
						CopyBlock (*i_R, *i_D, col_idx, block->D_size, block->col_idx + block->size);
						++i_D;
					}

					col_idx += block->D_size;
				}

				if (idx >= block_begin->col_idx + block_begin->size) {
					++block_begin;

					if (block_begin != blocks.end ())
						idx = block_begin->col_idx;
				}
				else
					++idx;
			}

			commentator.stop (MSG_DONE);
		}

	public:
		/**
		 * \brief Construct a new F4Solver
		 *
		 * @param _F Field over which to do computations
		 *
		 * @param _threshold Real number between 0 and 1. If
		 * the portion of nonzero entries in a row is at least
		 * this value, then the row is considered to be
		 * dense. Default 0.2.
		 */
		F4Solver (const Field &_F, double _threshold = 0.2) : F (_F), VD (_F), MD (_F), GJ (_F), threshold (_threshold) {
			F.init (one, 1);
			F.init (neg_one, -1);
		}

		/** 
		 * \brief Convert the matrix A into reduced
		 * row-echelon form
		 *
		 * @param R Matrix into which to store the reduced
		 * row-echelon form. Should have the same dimensions
		 * as A. May be the same matrix as A, in which case A
		 * is replaced by its reduced row-echelon form.
		 *
		 * @param A Matrix to be converted to row-echelon
		 * form. Not altered, unless R is the same matrix.
		 *
		 * @param rank Integer-reference into which to store
		 * computed rank
		 *
		 * @param det Field-element-reference into which to
		 * store computed determinant of pivot-submatrix
		 */
		void RowEchelonForm (SparseMatrix &R, const SparseMatrix &X, size_t &rank, Element &det) {
			commentator.start ("Reduction of F4-matrix to reduced row-echelon form", __FUNCTION__);

			std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

			std::vector<Block> blocks1, blocks2, blocks_combined;
			std::vector<size_t> pivot_rows, pivot_rows2;

			FindPivotRows (X, blocks1, pivot_rows);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Found " << pivot_rows.size () << " pivots in " << blocks1.size () << " blocks" << std::endl;

			SparseMatrix X1 (pivot_rows.size (), X.coldim ()), X2 (X.rowdim () - pivot_rows.size (), X.coldim ());

			AssemblePivotRows (X, X1, X2, pivot_rows);

			DenseMatrix A (pivot_rows.size (), pivot_rows.size ());
			DenseMatrix B (pivot_rows.size (), X.coldim () - pivot_rows.size ());
			DenseMatrix C (X.rowdim () - pivot_rows.size (), pivot_rows.size ());
			DenseMatrix D (X.rowdim () - pivot_rows.size (), X.coldim () - pivot_rows.size ());

			SplitMatrix (X1, A, B, blocks1);
			SplitMatrix (X2, C, D, blocks1);

			reportUI << "Matrix A:" << std::endl;
			MD.write (reportUI, A);
			reportUI << "Matrix B:" << std::endl;
			MD.write (reportUI, B);
			reportUI << "Matrix C:" << std::endl;
			MD.write (reportUI, C);
			reportUI << "Matrix D:" << std::endl;
			MD.write (reportUI, D);

			commentator.start ("Constructing D - C A^-1 B");

			MD.trsm (F.one (), A, B, UpperTriangular, true);
			MD.gemm (F.minusOne (), C, B, F.one (), D);

			commentator.stop (MSG_DONE);

			reportUI << "A^-1 B:" << std::endl;
			MD.write (reportUI, B);

			reportUI << "D - C A^-1 B:" << std::endl;
			MD.write (reportUI, D);

			DenseMatrix U (D.rowdim (), D.rowdim ());
			typename GaussJordan<Field>::Permutation P;
			size_t r_D;
			typename Field::Element d_D;

			GJ.DenseRowEchelonForm (D, U, P, r_D, d_D);

			reportUI << "Row-echelon form of D - C A^-1 B:" << std::endl;
			MD.write (reportUI, D);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Rank of dense part is " << r_D << std::endl;

			FindPivotRows (D, blocks2, pivot_rows2);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "(In D) found " << pivot_rows2.size () << " pivots in " << blocks2.size () << " blocks" << std::endl;

			DenseMatrix B1 (B.rowdim (), pivot_rows2.size ());
			DenseMatrix B2 (B.rowdim (), D.coldim () - pivot_rows2.size ());
			DenseMatrix D1 (D.rowdim (), pivot_rows2.size ());
			DenseMatrix D2 (D.rowdim (), D.coldim () - pivot_rows2.size ());

			SplitMatrix (B, B1, B2, blocks2);
			SplitMatrix (D, D1, D2, blocks2);

			reportUI << "Matrix B1:" << std::endl;
			MD.write (reportUI, B1);
			reportUI << "Matrix B2:" << std::endl;
			MD.write (reportUI, B2);
			reportUI << "Matrix D1:" << std::endl;
			MD.write (reportUI, D1);
			reportUI << "Matrix D2:" << std::endl;
			MD.write (reportUI, D2);

			DenseMatrix D1p (D1, 0, 0, pivot_rows2.size (), pivot_rows2.size ());
			DenseMatrix D2p (D2, 0, 0, pivot_rows2.size (), D2.coldim ());

			commentator.start ("Constructing B2 - B1 D1^-1 D2");

			MD.trsm (F.one (), D1p, D2p, UpperTriangular, true);			
			MD.gemm (F.minusOne (), B1, D2p, F.one (), B2);

			commentator.stop (MSG_DONE);

			CombineBlocks (blocks_combined, blocks1, blocks2);

			AssembleOutput (R, B2, D2, blocks_combined, pivot_rows);

			rank = pivot_rows.size () + r_D;
			det = d_D;

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}
	};
}

#endif // __F4_SOLVER_H
