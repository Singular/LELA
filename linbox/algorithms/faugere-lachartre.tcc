/* linbox/algorithms/faugere-lachartre.tcc
 * Copyright 2010 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 *
 * Implementation of algorithm for computing reduced row-echelon form
 * of a matrix coming from the F4-algorithm
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ALGORITHMS_FAUGERE_LACHARTRE_TCC
#define __LINBOX_ALGORITHMS_FAUGERE_LACHARTRE_TCC

#include "linbox/algorithms/faugere-lachartre.h"
#include "linbox/blas/level1.h"
#include "linbox/blas/level3.h"
#include "linbox/solutions/echelon-form.h"
#include "linbox/solutions/echelon-form-gf2.h"

#ifndef PROGRESS_STEP
#  define PROGRESS_STEP 1024
#endif // PROGRESS_STEP

namespace LinBox
{

template <class Ring, class Modules>
FaugereLachartre<Ring, Modules>::FaugereLachartre (Context<Ring, Modules> &_ctx)
	: ctx (_ctx), EF (_ctx) {}

template <class Ring, class Modules>
template <class Matrix>
void FaugereLachartre<Ring, Modules>::setup_splicer (Splicer &splicer, Splicer &reconst_splicer, const Matrix &A, size_t &num_pivot_rows, typename Ring::Element &det) const
{
	commentator.start ("Finding pivot-rows", __FUNCTION__);

	typename Matrix::ConstRowIterator i_A;
	int last_col = -1, col, first_col_in_block = 0, height = 0, dest_col_tr = 0, dest_col_res = 0;
	int row = 0, first_row_in_block = 0, dest_row_tr = 0, dest_row_res = 0;
	bool last_was_same_col = false;
	typename Ring::Element a;

	num_pivot_rows = 0;

	for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A, ++row) {
		col = BLAS1::head (ctx, a, *i_A);

		if (col == -1) {
			if (!last_was_same_col) {
				splicer.addHorizontalBlock (Block (0, 0, first_row_in_block, dest_row_tr, row - first_row_in_block));
				dest_row_tr += row - first_row_in_block;
				first_row_in_block = row;
			}
			break;
		}
		else if (col == last_col) {
			if (!last_was_same_col) {
				splicer.addHorizontalBlock (Block (0, 0, first_row_in_block, dest_row_tr, row - first_row_in_block));
				dest_row_tr += row - first_row_in_block;
				first_row_in_block = row;
				last_was_same_col = true;
			}
		}
		else {
			if (last_was_same_col) {
				splicer.addHorizontalBlock (Block (0, 1, first_row_in_block, dest_row_res, row - first_row_in_block));
				dest_row_res += row - first_row_in_block;
				first_row_in_block = row;
				last_was_same_col = false;
			}

			if (col == last_col + 1) {
				last_col = col;
				++height;
				++num_pivot_rows;
			} else {
				if (height > 0) {
					Block newblock (0, 0, dest_col_tr, first_col_in_block, height);
							
					splicer.addVerticalBlock (Block (0, 0, first_col_in_block, dest_col_tr, height));
					reconst_splicer.addVerticalBlock (newblock);
					reconst_splicer.addHorizontalBlock (newblock);
				}

				if (col - last_col > 1) {
					Block newblock (1, 0, dest_col_res, first_col_in_block + height, col - last_col - 1);
							
					splicer.addVerticalBlock (Block (0, 1, first_col_in_block + height, dest_col_res, col - last_col - 1));
					reconst_splicer.addVerticalBlock (newblock);
					reconst_splicer.addHorizontalBlock (newblock);
				}

				dest_col_tr += height;
				dest_col_res += col - last_col - 1;
				height = 1;
				first_col_in_block = last_col = col;
				++num_pivot_rows;
			}
		}

		ctx.F.mulin (det, a);
	}

	if (height > 0) {
		Block newblock (0, 0, dest_col_tr, first_col_in_block, height);
				
		splicer.addVerticalBlock (Block (0, 0, first_col_in_block, dest_col_tr, height));
		reconst_splicer.addVerticalBlock (newblock);
		reconst_splicer.addHorizontalBlock (newblock);
	}

	if (first_col_in_block + height < (int) A.coldim ()) {
		Block newblock (1, 0, dest_col_res, first_col_in_block + height, A.coldim () - first_col_in_block - height);
				
		splicer.addVerticalBlock (Block (0, 1, first_col_in_block + height, dest_col_res, A.coldim () - first_col_in_block - height));
		reconst_splicer.addVerticalBlock (newblock);
		reconst_splicer.addHorizontalBlock (newblock);
	}

	if (first_row_in_block < (int) A.rowdim ()) {
		if (row < (int) A.rowdim () || last_was_same_col)
			splicer.addHorizontalBlock (Block (0, 1, first_row_in_block, dest_row_res, A.rowdim () - first_row_in_block));
		else
			splicer.addHorizontalBlock (Block (0, 0, first_row_in_block, dest_row_tr, A.rowdim () - first_row_in_block));
	}

	commentator.stop (MSG_DONE, NULL, __FUNCTION__);
}

template <class Ring, class Modules>
template <class Matrix>
void FaugereLachartre<Ring, Modules>::RowEchelonForm (Matrix &R, const Matrix &X, size_t &rank, typename Ring::Element &det)
{
	commentator.start ("Reduction of F4-matrix to reduced row-echelon form", __FUNCTION__);

	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Splicer X_splicer, X_reconst_splicer;

	size_t num_pivot_rows;

	setup_splicer (X_splicer, X_reconst_splicer, X, num_pivot_rows, det);
	rank = num_pivot_rows;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Found " << num_pivot_rows << " pivots" << std::endl;
	reportUI << "Splicer:" << std::endl << X_splicer << std::endl;

	DenseMatrix A (num_pivot_rows, num_pivot_rows);
	DenseMatrix B (num_pivot_rows, X.coldim () - num_pivot_rows);
	DenseMatrix C (X.rowdim () - num_pivot_rows, num_pivot_rows);
	DenseMatrix D (X.rowdim () - num_pivot_rows, X.coldim () - num_pivot_rows);

	MatrixPart<DenseMatrix> X_targets_inner_1[] = { MatrixPart<DenseMatrix> (A), MatrixPart<DenseMatrix> (B) };
	MatrixPart<DenseMatrix> X_targets_inner_2[] = { MatrixPart<DenseMatrix> (C), MatrixPart<DenseMatrix> (D) };

	MatrixPart<DenseMatrix> *X_targets[] = { X_targets_inner_1, X_targets_inner_2 };

	MatrixPart<const Matrix> X_sources_1[] = { MatrixPart<const Matrix> (X) };
	MatrixPart<const Matrix> *X_sources[] = { X_sources_1 };

	X_splicer.splice<Ring, const Matrix, DenseMatrix> (ctx.F, X_sources, X_targets);

	reportUI << "Matrix A:" << std::endl;
	BLAS3::write (ctx, reportUI, A);
	reportUI << "Matrix B:" << std::endl;
	BLAS3::write (ctx, reportUI, B);
	reportUI << "Matrix C:" << std::endl;
	BLAS3::write (ctx, reportUI, C);
	reportUI << "Matrix D:" << std::endl;
	BLAS3::write (ctx, reportUI, D);

	// std::ofstream Aout ("A.png");
	// BLAS3::write (ctx, Aout, A, FORMAT_PNG);
	// std::ofstream Bout ("B.png");
	// BLAS3::write (ctx, Bout, B, FORMAT_PNG);
	// std::ofstream Cout ("C.png");
	// BLAS3::write (ctx, Cout, C, FORMAT_PNG);
	// std::ofstream Dout ("D.png");
	// BLAS3::write (ctx, Dout, D, FORMAT_PNG);

	commentator.start ("Constructing D - C A^-1 B");

	BLAS3::trsm (ctx, ctx.F.one (), A, B, UpperTriangular, true);
	BLAS3::gemm (ctx, ctx.F.minusOne (), C, B, ctx.F.one (), D);

	commentator.stop (MSG_DONE);

	// std::ofstream ABout ("AB.png");
	// BLAS3::write (ctx, ABout, B, FORMAT_PNG);
	// std::ofstream DCABout ("D-CAB.png");
	// BLAS3::write (ctx, DCABout, D, FORMAT_PNG);

	reportUI << "A^-1 B:" << std::endl;
	BLAS3::write (ctx, reportUI, B);

	reportUI << "D - C A^-1 B:" << std::endl;
	BLAS3::write (ctx, reportUI, D);

	// size_t r_D;

	EF.RowEchelonForm (D);

	reportUI << "Row-echelon form of D - C A^-1 B:" << std::endl;
	BLAS3::write (ctx, reportUI, D);

	// commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	// 	<< "Rank of dense part is " << r_D << std::endl;

	Splicer D_splicer, D_reconst_splicer;

	setup_splicer (D_splicer, D_reconst_splicer, D, num_pivot_rows, det);
	rank += num_pivot_rows;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "(In D) found " << num_pivot_rows << " pivots" << std::endl;
	reportUI << "Splicer:" << std::endl << D_splicer << std::endl;

	DenseMatrix B1 (B.rowdim (), num_pivot_rows);
	DenseMatrix B2 (B.rowdim (), D.coldim () - num_pivot_rows);
	DenseMatrix D1 (num_pivot_rows, num_pivot_rows);
	DenseMatrix D2 (num_pivot_rows, D.coldim () - num_pivot_rows);

	MatrixPart<DenseMatrix> B_targets_inner[] = { MatrixPart<DenseMatrix> (B1), MatrixPart<DenseMatrix> (B2) };
	MatrixPart<DenseMatrix> D_targets_inner_1[] = { MatrixPart<DenseMatrix> (D1), MatrixPart<DenseMatrix> (D2) };
	MatrixPart<DenseMatrix> D_targets_inner_2[] = { MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO),
							MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO) };

	MatrixPart<DenseMatrix> *B_targets[] = { B_targets_inner };
	MatrixPart<DenseMatrix> *D_targets[] = { D_targets_inner_1, D_targets_inner_2 };

	MatrixPart<DenseMatrix> B_sources_1[] = { MatrixPart<DenseMatrix> (B) };
	MatrixPart<DenseMatrix> D_sources_1[] = { MatrixPart<DenseMatrix> (D) };

	MatrixPart<DenseMatrix> *B_sources[] = { B_sources_1 };
	MatrixPart<DenseMatrix> *D_sources[] = { D_sources_1 };

	Splicer B_splicer (D_splicer);
	B_splicer.clearHorizontalBlocks ();
	B_splicer.addHorizontalBlock (Block (0, 0, 0, 0, B.rowdim ()));

	B_splicer.splice<Ring, DenseMatrix, DenseMatrix> (ctx.F, B_sources, B_targets);
	D_splicer.splice<Ring, DenseMatrix, DenseMatrix> (ctx.F, D_sources, D_targets);

	reportUI << "Matrix B1:" << std::endl;
	BLAS3::write (ctx, reportUI, B1);
	reportUI << "Matrix B2:" << std::endl;
	BLAS3::write (ctx, reportUI, B2);
	reportUI << "Matrix D1:" << std::endl;
	BLAS3::write (ctx, reportUI, D1);
	reportUI << "Matrix D2:" << std::endl;
	BLAS3::write (ctx, reportUI, D2);

	commentator.start ("Constructing B2 - B1 D1^-1 D2");

	BLAS3::trsm (ctx, ctx.F.one (), D1, D2, UpperTriangular, true);			
	BLAS3::gemm (ctx, ctx.F.minusOne (), B1, D2, ctx.F.one (), B2);

	commentator.stop (MSG_DONE);

	reportUI << "B2 - B1 D1^-1 D2:" << std::endl;
	BLAS3::write (ctx, reportUI, B2);

	Splicer composed_splicer, subst_splicer, D_splicer_rev;

	D_splicer.reverse (D_splicer_rev);
	X_reconst_splicer.compose (subst_splicer, D_reconst_splicer, 1, (unsigned int) -1, 0, (unsigned int) -1);
	subst_splicer.removeGaps ();
	subst_splicer.consolidate ();
	subst_splicer.compose (composed_splicer, D_splicer_rev, 1, 1);
	composed_splicer.fillHorizontal (2, 0, X.rowdim ());

	reportUI << "Splicer after substitution:" << std::endl << subst_splicer << std::endl;
	reportUI << "Composed splicer:" << std::endl << composed_splicer << std::endl;

	MatrixPart<DenseMatrix> R_sources_inner_1[] = { MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_IDENTITY),
							MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO),
							MatrixPart<DenseMatrix> (B2) };
	MatrixPart<DenseMatrix> R_sources_inner_2[] = { MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO),
							MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_IDENTITY),
							MatrixPart<DenseMatrix> (D2) };
	MatrixPart<DenseMatrix> R_sources_inner_3[] = { MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO),
							MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO),
							MatrixPart<DenseMatrix> (MatrixPart<DenseMatrix>::TYPE_ZERO) };

	MatrixPart<DenseMatrix> *R_sources[] = { R_sources_inner_1, R_sources_inner_2, R_sources_inner_3 };

	MatrixPart<Matrix> R_targets_1[] = { MatrixPart<Matrix> (R) };
	MatrixPart<Matrix> *R_targets[] = { R_targets_1 };

	BLAS3::scal (ctx, ctx.F.zero (), R);

	composed_splicer.splice<Ring, DenseMatrix, Matrix> (ctx.F, R_sources, R_targets);

	commentator.stop (MSG_DONE, NULL, __FUNCTION__);
}

} // namespace LinBox

#endif // __LINBOX_ALGORITHMS_FAUGERE_LACHARTRE_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
