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

	class Block {
		size_t orig_idx;
		size_t dest_idx;
		size_t tr_size;
		size_t res_size;

		friend std::ostream &operator << (std::ostream &os, const Block &b);

	public:
		Block (size_t _orig_idx, size_t _dest_idx, size_t _tr_size, size_t _res_size)
			: dest_idx (_dest_idx), orig_idx (_orig_idx), tr_size (_tr_size), res_size (_res_size) {}

		inline size_t dest_idx_residual_part () const { return dest_idx; }
		inline size_t orig_idx_triangular_part () const { return orig_idx; }
		inline size_t orig_idx_residual_part () const { return orig_idx + tr_size; }
		inline size_t orig_idx_next_block () const { return orig_idx + tr_size + res_size; }
		inline size_t map_dest_to_orig_idx (size_t idx) const { return idx - dest_idx + orig_idx; }

		inline size_t get_tr_size () const { return tr_size; }
		inline size_t get_res_size () const { return res_size; }

		inline bool is_orig_idx_in_tr_block (size_t idx) const { return idx >= orig_idx && idx < orig_idx + tr_size; }
		inline bool is_orig_idx_in_res_block (size_t idx) const { return idx >= orig_idx + tr_size && idx < orig_idx + tr_size + res_size; }
		inline bool is_dest_idx_in_block (size_t idx) const { return idx >= dest_idx && idx < dest_idx + res_size; }
	};

	std::ostream &operator << (std::ostream &os, const Block &b)
	{
		os << "Block (" << b.orig_idx << "," << b.dest_idx << "," << b.tr_size << "," << b.res_size << ")";
		return os;
	}

	std::ostream &operator << (std::ostream &os, const std::vector<Block> &v)
	{
		std::vector<Block>::const_iterator i;

		for (i = v.begin (); i != v.end (); ++i)
			os << *i << std::endl;

		return os;
	}

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

		template <class Matrix>
		void FindPivotRows (const Matrix &A, std::vector<Block> &blocks, std::vector<size_t> &pivot_rows) const
		{
			commentator.start ("Finding pivot-rows", __FUNCTION__);

			typename Matrix::ConstRowIterator i_A;
			int last_col = -1, col, first_col_in_block = 0, height = 0, dest_col = 0;
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
					blocks.push_back (Block (first_col_in_block, dest_col, height, col - last_col - 1));
					pivot_rows.push_back (i_A - A.rowBegin ());
					dest_col += col - last_col - 1;
					height = 1;
					first_col_in_block = last_col = col;
				}
			}

			blocks.push_back (Block (first_col_in_block, dest_col, height, A.coldim () - first_col_in_block - height));

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
				Submatrix<const Matrix1> A_part_source (X, 0, i->orig_idx_triangular_part (), X.rowdim (), i->get_tr_size ());
				Submatrix<Matrix2> A_part_dest (A, 0, prev_col_A, A.rowdim (), i->get_tr_size ());

				Submatrix<const Matrix1> B_part_source (X, 0, i->orig_idx_residual_part (), X.rowdim (), i->get_res_size ());
				Submatrix<Matrix2> B_part_dest (B, 0, prev_col_B, B.rowdim (), i->get_res_size ());

				MD.copy (A_part_dest, A_part_source);
				MD.copy (B_part_dest, B_part_source);

				prev_col_A += i->get_tr_size ();
				prev_col_B += i->get_res_size ();
			}

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
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
				if (*i)
					dest.push_back (typename Vector1::value_type (idx, *i));

			return dest;
		}

		std::vector<Block> &MapBlocks (std::vector<Block> &output, const std::vector<Block> &B_blocks, const std::vector<Block> &D_blocks)
		{
			std::vector<Block>::const_iterator B_block = B_blocks.begin ();			
			std::vector<Block>::const_iterator D_block;

			for (D_block = D_blocks.begin (); D_block != D_blocks.end (); ++D_block) {
				while (B_block != B_blocks.end () && !B_block->is_dest_idx_in_block (D_block->orig_idx_residual_part ()))
					++B_block;

				linbox_check (B_block != B_blocks.end ());

				output.push_back (Block (B_block->map_dest_to_orig_idx (D_block->orig_idx_triangular_part ()),
							 D_block->dest_idx_residual_part (),
							 D_block->get_tr_size (),
							 D_block->get_res_size ()));
			}

			return output;
		}

		template <class Vector1, class Vector2>
		Vector1 &CopyBlock (Vector1 &dest, const Vector2 &source, size_t source_idx, size_t len, size_t dest_idx) const
			{ return CopyBlockSpecialised (dest, source, source_idx, len, dest_idx,
						       typename VectorTraits<Field, Vector1>::VectorCategory (),
						       typename VectorTraits<Field, Vector2>::VectorCategory ()); }

		void AssembleOutput (SparseMatrix &R, const DenseMatrix &B, const DenseMatrix &D, const std::vector<Block> &B_blocks, const std::vector<Block> &D_blocks) const
		{
			commentator.start ("Assembly of final output", __FUNCTION__, R.rowdim () / 1024);

			typename SparseMatrix::RowIterator i_R = R.rowBegin (), i_R_end;
			typename DenseMatrix::ConstRowIterator i_B = B.rowBegin (), i_D = D.rowBegin ();
			typename std::vector<Block>::const_iterator B_block = B_blocks.begin (), D_block = D_blocks.begin (), block;
			typename Field::Element a;

			size_t pivot_idx, col_idx, B_col = 0;

			while (i_R != R.rowEnd () && B_block != B_blocks.end ()) {
				i_R_end = i_R + B_block->get_tr_size ();
				pivot_idx = B_block->orig_idx_triangular_part ();

				for (; i_R != i_R_end; ++i_R, ++i_B, ++pivot_idx) {
					linbox_check (i_R != R.rowEnd ());
					linbox_check (i_B != B.rowEnd ());

					i_R->clear ();

					i_R->push_back (typename SparseMatrix::Row::value_type
							(pivot_idx >> WordTraits<typename SparseMatrix::Row::word_type>::logof_size,
							 SparseMatrix::Row::Endianness::e_j (pivot_idx & WordTraits<typename SparseMatrix::Row::word_type>::pos_mask)));

					for (block = B_block, col_idx = 0; block != B_blocks.end (); ++block) {
						CopyBlock (*i_R, *i_B, col_idx, block->get_res_size (), block->orig_idx_residual_part ());
						col_idx += block->get_res_size ();
					}
				}

				B_col += B_block->get_res_size ();

				for (; D_block != D_blocks.end () && B_block->is_orig_idx_in_res_block (D_block->orig_idx_residual_part ()); ++D_block) {
					i_R_end = i_R + D_block->get_tr_size ();
					pivot_idx = D_block->orig_idx_triangular_part ();

					for (; i_R != i_R_end; ++i_R, ++i_D, ++pivot_idx) {
						linbox_check (i_R != R.rowEnd ());
						linbox_check (i_D != D.rowEnd ());

						i_R->clear ();

						i_R->push_back (typename SparseMatrix::Row::value_type
								(pivot_idx >> WordTraits<typename SparseMatrix::Row::word_type>::logof_size,
								 SparseMatrix::Row::Endianness::e_j (pivot_idx & WordTraits<typename SparseMatrix::Row::word_type>::pos_mask)));

						for (block = D_block, col_idx = 0; block != D_blocks.end (); ++block) {
							CopyBlock (*i_R, *i_D, col_idx, block->get_res_size (), block->orig_idx_residual_part ());
							col_idx += block->get_res_size ();
						}
					}
				}

				++B_block;
			}

			for (; i_R != R.rowEnd (); ++i_R)
				i_R->clear ();

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

			std::vector<Block> B_blocks, D_blocks, D_blocks_mapped;
			std::vector<size_t> pivot_rows, pivot_rows2;

			FindPivotRows (X, B_blocks, pivot_rows);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Found " << pivot_rows.size () << " pivots in " << B_blocks.size () << " blocks" << std::endl;
			reportUI << "Blocks:" << std::endl << B_blocks << std::endl;

			SparseMatrix X1 (pivot_rows.size (), X.coldim ()), X2 (X.rowdim () - pivot_rows.size (), X.coldim ());

			AssemblePivotRows (X, X1, X2, pivot_rows);

			DenseMatrix A (pivot_rows.size (), pivot_rows.size ());
			DenseMatrix B (pivot_rows.size (), X.coldim () - pivot_rows.size ());
			DenseMatrix C (X.rowdim () - pivot_rows.size (), pivot_rows.size ());
			DenseMatrix D (X.rowdim () - pivot_rows.size (), X.coldim () - pivot_rows.size ());

			SplitMatrix (X1, A, B, B_blocks);
			SplitMatrix (X2, C, D, B_blocks);

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

			FindPivotRows (D, D_blocks, pivot_rows2);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "(In D) found " << pivot_rows2.size () << " pivots in " << D_blocks.size () << " blocks" << std::endl;
			reportUI << "Blocks:" << std::endl << D_blocks << std::endl;

			DenseMatrix B1 (B.rowdim (), pivot_rows2.size ());
			DenseMatrix B2 (B.rowdim (), D.coldim () - pivot_rows2.size ());
			DenseMatrix D1 (D.rowdim (), pivot_rows2.size ());
			DenseMatrix D2 (D.rowdim (), D.coldim () - pivot_rows2.size ());

			SplitMatrix (B, B1, B2, D_blocks);
			SplitMatrix (D, D1, D2, D_blocks);

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

			reportUI << "B2 - B1 D1^-1 D2:" << std::endl;
			MD.write (reportUI, B2);

			MapBlocks (D_blocks_mapped, B_blocks, D_blocks);

			reportUI << "Mapped blocks:" << std::endl << D_blocks_mapped << std::endl;

			AssembleOutput (R, B2, D2, B_blocks, D_blocks_mapped);

			rank = pivot_rows.size () + r_D;
			det = d_D;

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}
	};
}

#endif // __F4_SOLVER_H
