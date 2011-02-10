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

#include "gauss-jordan.h"

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

		size_t DenseSize (const SparseMatrix &A, const std::vector<size_t> &pivot_rows, size_t &rowdim, size_t &coldim) const
		{
			rowdim = A.rowdim () - pivot_rows.size ();
			coldim = A.coldim () - pivot_rows.size ();
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

		void ConstructDensePart (const SparseMatrix &Agj, const SparseMatrix &B, DenseMatrix &D, const std::vector<Block> &blocks) const
		{
			commentator.start ("Construction of dense matrix for reduction", __FUNCTION__, blocks.size ());

			typename std::vector<Block>::const_iterator b, c, b_end = blocks.end ();
			size_t D_col = 0, A_row = 0, b_width;

			for (b = blocks.begin (); b != b_end; ++b) {
				b_width = GetNextBlockColumn (B.coldim (), b, b_end) - (b->col_idx + b->size);

				SubmatrixBase<DenseMatrix> Dp (D, 0, D_col, D.rowdim (), b_width);
				SubmatrixBase<const SparseMatrix> Ds (B, 0, b->col_idx + b->size, D.rowdim (), b_width);

				MD.copy (Dp, Ds);

				// DEBUG
				// std::cout << __FUNCTION__ << ": block start " << b->col_idx << ", size " << b->size << ", width " << b_width << ", source-matrix is" << std::endl;
				// MD.write (std::cout, Ds);
				// std::cout << __FUNCTION__ << ": Initial dest-matrix is" << std::endl;
				// MD.write (std::cout, Dp);

				A_row = 0;

				for (c = blocks.begin (); c <= b; ++c) {
					SubmatrixBase<const SparseMatrix> B1 (B, 0, c->col_idx, D.rowdim (), c->size);
					SubmatrixBase<const SparseMatrix> A1 (Agj, A_row, b->col_idx + b->size, c->size, b_width);

					MD.gemm (neg_one, B1, A1, one, Dp);

					A_row += c->size;
				}

				D_col += b_width;

				commentator.progress (b - blocks.begin ());
			}

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}

		// FIXME: This is quite inefficient, since it ends up copying zero much of the time
		template <class Vector1, class Vector2>
		void ExpandVector (Vector1 &x, const Vector2 &y, const std::vector<Block> &blocks, size_t coldim) const
		{
			static BitVector v;

			v.resize (coldim);

			std::fill (v.wordBegin (), v.wordEnd (), 0ULL);

			typename std::vector<Block>::const_iterator b, b_end = blocks.end ();
			size_t next_block_col, start_block_y = 0, next_block_y;

			for (b = blocks.begin (); b != b_end; ++b) {
				next_block_col = GetNextBlockColumn (coldim, b, b_end);
				next_block_y = start_block_y + (next_block_col - (b->col_idx + b->size));

				BitSubvector<BitVector::iterator> vs (v.begin () + (b->col_idx + b->size), v.begin () + next_block_col);
				BitSubvector<typename Vector2::const_iterator> ys (y.begin () + start_block_y, y.begin () + next_block_y);

				VD.copy (vs, ys);

				start_block_y = next_block_y;
			}

			VD.copy (x, v);
		}

		size_t EffectiveColumn (size_t idx, const std::vector<Block> &blocks) const
		{
			typename std::vector<Block>::const_iterator b, b_end;
			size_t next_column;

			for (b = blocks.begin (); b != b_end; ++b) {
				if (idx < b->D_size)
					return idx + b->col_idx + b->size;

				idx -= b->D_size;
			}

			return next_column;
		}

		void AssembleOutput (const SparseMatrix &Agj, const DenseMatrix &D, SparseMatrix &R, const std::vector<Block> &blocks) const
		{
			commentator.start ("Assembly of final output", __FUNCTION__, R.rowdim () / 1024);

			typename SparseMatrix::RowIterator i_R;
			typename SparseMatrix::ConstRowIterator i_Agj = Agj.rowBegin ();
			typename DenseMatrix::ConstRowIterator i_D = D.rowBegin ();

			typename Field::Element a;

			size_t first_col_A, first_col_D;
			int idx;

			for (i_R = R.rowBegin (); i_R != R.rowEnd (); ++i_R) {
				if (i_Agj != Agj.rowEnd ())
					first_col_A = VD.firstNonzeroEntry (a, *i_Agj);
				else
					first_col_A = Agj.coldim ();

				if (i_D != D.rowEnd ()) {
					idx = VD.firstNonzeroEntry (a, *i_D);

					if (idx == -1)
						first_col_D = Agj.coldim ();
					else
						first_col_D = EffectiveColumn ((size_t) idx, blocks);
				}
				else
					first_col_D = Agj.coldim ();

				if (first_col_A < first_col_D) {
					VD.copy (*i_R, *i_Agj);
					++i_Agj;
				}
				else if (first_col_A > first_col_D) {
					ExpandVector (*i_R, *i_D, blocks, R.coldim ());
					++i_D;
				}
				else if (i_Agj == Agj.rowEnd ())
					VD.mulin (*i_R, false);
				else
					throw LinboxError ("Got two pivots in the same column -- this shouldn't ever happen!");

				if (!(i_R - R.rowBegin () + 1) % 1024)
					commentator.progress ((i_R - R.rowBegin () + 1) / 1024);
			}

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
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
		void RowEchelonForm (SparseMatrix &R, const SparseMatrix &A, size_t &rank, Element &det) {
			commentator.start ("Reduction of F4-matrix to reduced row-echelon form", __FUNCTION__);

			std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

			std::vector<Block> blocks;
			std::vector<size_t> pivot_rows;

			size_t D_rowdim, D_coldim;

			FindPivotRows (A, blocks, pivot_rows);
			DenseSize (A, pivot_rows, D_rowdim, D_coldim);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Found " << pivot_rows.size () << " pivots in " << blocks.size () << " blocks" << std::endl;

			SparseMatrix Agj (A.rowdim () - D_rowdim, A.coldim ()), B (D_rowdim, A.coldim ());
			DenseMatrix D (D_rowdim, D_coldim);

			AssemblePivotRows (A, Agj, B, pivot_rows);

			reportUI << "Matrix of pivots:" << std::endl;
			MD.write (reportUI, Agj);
			reportUI << "Residual block:" << std::endl;
			MD.write (reportUI, B);

			GJ.ReduceRowEchelon (Agj, D, false, Agj.rowdim ());

			reportUI << "Reduced echelon-form of pivot-matrix:" << std::endl;
			MD.write (reportUI, Agj);

			ConstructDensePart (Agj, B, D, blocks);

			reportUI << "Dense block:" << std::endl;
			MD.write (reportUI, D);

			DenseMatrix U (D_rowdim, D_rowdim);
			typename GaussJordan<Field>::Permutation P;
			size_t r_D;
			typename Field::Element d_D;

			GJ.DenseRowEchelonForm (D, U, P, D, r_D, d_D);

			reportUI << "Reduced row-echelon form of dense block:" << std::endl;
			MD.write (reportUI, D);

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Rank of dense part is " << r_D << std::endl;

			AssembleOutput (Agj, D, R, blocks);

			rank = pivot_rows.size () + r_D;
			det = true;

			GJ.ReduceRowEchelon (R, D, false, rank);

			commentator.stop (MSG_DONE, NULL, __FUNCTION__);
		}
	};
}

#endif // __F4_SOLVER_H
