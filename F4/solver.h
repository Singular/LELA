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
#include <linbox/solutions/echelon-form.h>
#include <linbox/solutions/echelon-form-gf2.h>

#ifndef PROGRESS_STEP
#  define PROGRESS_STEP 1024
#endif // PROGRESS_STEP

namespace F4 {
	using namespace LinBox;

	class Block {
	public:
		enum BlockType { TYPE_TRIANGULAR, TYPE_RESIDUAL };

	private:

		BlockType type;
		size_t orig_idx;
		size_t dest_idx;
		size_t size;

		friend std::ostream &operator << (std::ostream &os, const Block &b);

	public:
		Block (BlockType _type, size_t _orig_idx, size_t _dest_idx, size_t _size)
			: type (_type), dest_idx (_dest_idx), orig_idx (_orig_idx), size (_size) {}

		inline BlockType get_type () const { return type; }

		inline size_t get_dest_idx () const { return dest_idx; }
		inline size_t get_orig_idx () const { return orig_idx; }
		inline size_t get_orig_idx_next_block () const { return orig_idx + size; }
		inline size_t map_dest_to_orig_idx (size_t idx) const { return idx - dest_idx + orig_idx; }

		inline size_t get_size () const { return size; }

		inline bool is_orig_idx_in_block (size_t idx) const { return idx >= orig_idx && idx < orig_idx + size; }
		inline bool is_dest_idx_in_block (size_t idx) const { return idx >= dest_idx && idx < dest_idx + size; }
	};

	std::ostream &operator << (std::ostream &os, const Block &b)
	{
		os << "Block (" << (b.type == Block::TYPE_TRIANGULAR ? "triangular," : "residual,") << b.orig_idx << "," << b.dest_idx << "," << b.size << ")";
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
		EchelonForm<Field> EF;

		Element one, neg_one;

		const double threshold;

		template <class Matrix>
		void FindPivotRows (const Matrix &A, std::vector<Block> &blocks, std::vector<size_t> &pivot_rows) const
		{
			commentator.start ("Finding pivot-rows", __FUNCTION__);

			typename Matrix::ConstRowIterator i_A;
			int last_col = -1, col, first_col_in_block = 0, height = 0, dest_col_tr = 0, dest_col_res = 0;
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
					if (height > 0)
						blocks.push_back (Block (Block::TYPE_TRIANGULAR, first_col_in_block, dest_col_tr, height));

					if (col - last_col > 1)
						blocks.push_back (Block (Block::TYPE_RESIDUAL, first_col_in_block + height, dest_col_res, col - last_col - 1));

					pivot_rows.push_back (i_A - A.rowBegin ());
					dest_col_tr += height;
					dest_col_res += col - last_col - 1;
					height = 1;
					first_col_in_block = last_col = col;
				}
			}

			if (height > 0)
				blocks.push_back (Block (Block::TYPE_TRIANGULAR, first_col_in_block, dest_col_tr, height));

			if (A.coldim () - first_col_in_block > height)
				blocks.push_back (Block (Block::TYPE_RESIDUAL, first_col_in_block + height, dest_col_res, A.coldim () - first_col_in_block - height));

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

			for (i = blocks.begin (); i != blocks.end (); ++i) {
				if (i->get_type () == Block::TYPE_TRIANGULAR) {
					Submatrix<const Matrix1> A_part_source (X, 0, i->get_orig_idx (), X.rowdim (), i->get_size ());
					Submatrix<Matrix2> A_part_dest (A, 0, i->get_dest_idx (), A.rowdim (), i->get_size ());

					MD.copy (A_part_dest, A_part_source);
				}
				else if (i->get_type () == Block::TYPE_RESIDUAL) {
					Submatrix<const Matrix1> B_part_source (X, 0, i->get_orig_idx (), X.rowdim (), i->get_size ());
					Submatrix<Matrix2> B_part_dest (B, 0, i->get_dest_idx (), B.rowdim (), i->get_size ());

					MD.copy (B_part_dest, B_part_source);
				}
			}

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}

		std::vector<Block> &MapBlocks (std::vector<Block> &output, const std::vector<Block> &X_blocks, const std::vector<Block> &D_blocks)
		{
			std::vector<Block>::const_iterator X_block;			
			std::vector<Block>::const_iterator D_block = D_blocks.begin ();

			size_t curr_orig_idx = 0, curr_dest_idx = 0, rest_size = 0, curr_size = 0;

			commentator.start ("Mapping blocks", __FUNCTION__);

			for (X_block = X_blocks.begin (); X_block != X_blocks.end (); ++X_block) {
				if (X_block->get_type () == Block::TYPE_TRIANGULAR)
					continue;

				curr_orig_idx = X_block->get_orig_idx ();

				if (rest_size > 0) {
					curr_size = std::min (rest_size, X_block->get_size ());

					output.push_back (Block (D_block->get_type (), curr_orig_idx, curr_dest_idx, curr_size));

					commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
						<< "New block created (leftover): " << output.back () << std::endl
						<< "Source-block from D: " << *D_block << std::endl
						<< "Source-block from X: " << *X_block << std::endl;

					curr_orig_idx += curr_size;
					curr_dest_idx += curr_size;
					rest_size -= curr_size;

					if (rest_size == 0)
						++D_block;
				}

				while (D_block != D_blocks.end () && rest_size == 0 && X_block->is_dest_idx_in_block (D_block->get_orig_idx ())) {
					curr_size = std::min (D_block->get_size (), X_block->get_orig_idx () + X_block->get_size () - curr_orig_idx);

					output.push_back (Block (D_block->get_type (), curr_orig_idx, D_block->get_dest_idx (), curr_size));

					commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
						<< "New block created: " << output.back () << std::endl
						<< "Source-block from D: " << *D_block << std::endl
						<< "Source-block from X: " << *X_block << std::endl;

					curr_orig_idx += curr_size;
					curr_dest_idx = D_block->get_dest_idx () + curr_size;
					rest_size = D_block->get_size () - curr_size;

					if (rest_size == 0)
						++D_block;
				}
			}

			commentator.stop (MSG_DONE);

			return output;
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

		template <class Vector1, class Vector2>
		Vector1 &CopyBlock (Vector1 &dest, const Vector2 &source, size_t source_idx, size_t len, size_t dest_idx) const
			{ return CopyBlockSpecialised (dest, source, source_idx, len, dest_idx,
						       typename VectorTraits<Field, Vector1>::VectorCategory (),
						       typename VectorTraits<Field, Vector2>::VectorCategory ()); }

		template <class Iterator1, class Iterator2>
		void AssembleRows (Iterator2 i_B,
				   Iterator1 i_R,
				   Iterator1 i_R_end,
				   typename std::vector<Block>::const_iterator start_block,
				   typename std::vector<Block>::const_iterator end_block,
				   size_t pivot_idx) const
		{
			typename std::vector<Block>::const_iterator block;

			for (; i_R != i_R_end; ++i_R, ++i_B, ++pivot_idx) {
				i_R->clear ();

				i_R->push_back (typename SparseMatrix::Row::value_type
						(pivot_idx >> WordTraits<typename SparseMatrix::Row::word_type>::logof_size,
						 SparseMatrix::Row::Endianness::e_j (pivot_idx & WordTraits<typename SparseMatrix::Row::word_type>::pos_mask)));

				for (block = start_block; block != end_block; ++block)
					if (block->get_type () == Block::TYPE_RESIDUAL)
						CopyBlock (*i_R, *i_B, block->get_dest_idx (), block->get_size (), block->get_orig_idx ());
			}
		}

		void AssembleOutput (SparseMatrix &R, const DenseMatrix &B, const DenseMatrix &D, const std::vector<Block> &X_blocks, const std::vector<Block> &D_blocks) const
		{
			commentator.start ("Assembly of final output", __FUNCTION__, R.rowdim () / 1024);

			typename SparseMatrix::RowIterator i_R = R.rowBegin (), i_R_end;
			typename DenseMatrix::ConstRowIterator i_B = B.rowBegin (), i_D = D.rowBegin ();
			typename std::vector<Block>::const_iterator X_block, D_block = D_blocks.begin ();

			for (X_block = X_blocks.begin (); X_block != X_blocks.end (); ++X_block) {
				if (X_block->get_type () == Block::TYPE_TRIANGULAR) {
					linbox_check (R.rowEnd () - i_R >= X_block->get_size ());
					linbox_check (B.rowEnd () - i_B >= X_block->get_size ());

					i_R_end = i_R + X_block->get_size ();
					AssembleRows (i_B, i_R, i_R_end, D_block, D_blocks.end (), X_block->get_orig_idx ());

					i_R = i_R_end;
					i_B += X_block->get_size ();
				}
				else if (X_block->get_type () == Block::TYPE_RESIDUAL) {
					for (; D_block != D_blocks.end () && X_block->is_orig_idx_in_block (D_block->get_orig_idx ()); ++D_block) {
						if (D_block->get_type () != Block::TYPE_TRIANGULAR)
							continue;

						linbox_check (R.rowEnd () - i_R >= D_block->get_size ());
						linbox_check (D.rowEnd () - i_D >= D_block->get_size ());

						i_R_end = i_R + D_block->get_size ();
						AssembleRows (i_D, i_R, i_R_end, D_block, D_blocks.end (), D_block->get_orig_idx ());

						i_R = i_R_end;
						i_D += D_block->get_size ();
					}
				}
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
		F4Solver (const Field &_F, double _threshold = 0.2) : F (_F), VD (_F), MD (_F), EF (_F), threshold (_threshold) {
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

			std::vector<Block> X_blocks, D_blocks, D_blocks_mapped;
			std::vector<size_t> pivot_rows, pivot_rows2;

			FindPivotRows (X, X_blocks, pivot_rows);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Found " << pivot_rows.size () << " pivots in " << X_blocks.size () << " blocks" << std::endl;
			reportUI << "Blocks:" << std::endl << X_blocks << std::endl;

			SparseMatrix X1 (pivot_rows.size (), X.coldim ()), X2 (X.rowdim () - pivot_rows.size (), X.coldim ());

			AssemblePivotRows (X, X1, X2, pivot_rows);

			DenseMatrix A (pivot_rows.size (), pivot_rows.size ());
			DenseMatrix B (pivot_rows.size (), X.coldim () - pivot_rows.size ());
			DenseMatrix C (X.rowdim () - pivot_rows.size (), pivot_rows.size ());
			DenseMatrix D (X.rowdim () - pivot_rows.size (), X.coldim () - pivot_rows.size ());

			SplitMatrix (X1, A, B, X_blocks);
			SplitMatrix (X2, C, D, X_blocks);

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

			// size_t r_D;

			EF.RowEchelonForm (D);

			reportUI << "Row-echelon form of D - C A^-1 B:" << std::endl;
			MD.write (reportUI, D);

			// commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			// 	<< "Rank of dense part is " << r_D << std::endl;

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

			MapBlocks (D_blocks_mapped, X_blocks, D_blocks);

			reportUI << "Mapped blocks:" << std::endl << D_blocks_mapped << std::endl;

			AssembleOutput (R, B2, D2, X_blocks, D_blocks_mapped);

			// rank = pivot_rows.size () + r_D;
			// det = d_D;

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}
	};
}

#endif // __F4_SOLVER_H
