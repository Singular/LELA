/* lela/algorithms/faugere-lachartre.tcc
 * Copyright 2010 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 *
 * Implementation of algorithm for computing reduced row-echelon form
 * of a matrix coming from the F4-algorithm
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ALGORITHMS_FAUGERE_LACHARTRE_TCC
#define __LELA_ALGORITHMS_FAUGERE_LACHARTRE_TCC

#include "lela/algorithms/faugere-lachartre.h"
#include "lela/blas/level1.h"
#include "lela/blas/level3.h"
#include "lela/solutions/echelon-form.h"
#include "lela/solutions/echelon-form-gf2.h"

#ifndef PROGRESS_STEP
#  define PROGRESS_STEP 1024
#endif // PROGRESS_STEP

namespace LELA
{

// Exception thrown when the matrix is not in the right form
class WrongMatrixForm {
	size_t _row;

public:
	WrongMatrixForm (size_t row) { _row = row; }
	friend std::ostream &operator << (std::ostream &out, const WrongMatrixForm &e);
};

std::ostream &operator << (std::ostream &out, const WrongMatrixForm &e)
	{ out << "Bad input-matrix for Faugère-Lachartre at row-index " << e._row << std::endl; return out; }

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
		else if (col < last_col)
			throw WrongMatrixForm (row);
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

template <class Ring, class Matrix1, class Matrix2, class Matrix3>
class MatrixGrid1
{
	const Ring &R;
	Matrix1 &X;
	Matrix2 &A;
	Matrix3 &B, &C, &D;

	template <class Container>
	void moveRowSpecialised (const Block &horiz_block, int row, const Container &vert_blocks, VectorRepresentationTypes::Generic)
	{
		typename Container::const_iterator vert_block;
		typename Matrix1::ConstRowIterator v_X = X.rowBegin () + (horiz_block.sourceIndex () + row);
		typename Matrix1::ConstRow::const_iterator i = v_X->begin ();

		if (horiz_block.dest () == 0) {
			typename Matrix2::RowIterator v_A = A.rowBegin () + (horiz_block.destIndex () + row);
			typename Matrix3::RowIterator v_B = B.rowBegin () + (horiz_block.destIndex () + row);

			for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
				if (vert_block->dest () == 0)
					i = Splicer::moveBlock (R, *v_A, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
								typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
				else
					i = Splicer::moveBlock (R, *v_B, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
								typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
			}
		} else {
			typename Matrix3::RowIterator v_C = C.rowBegin () + (horiz_block.destIndex () + row);
			typename Matrix3::RowIterator v_D = D.rowBegin () + (horiz_block.destIndex () + row);

			for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
				if (vert_block->dest () == 0)
					i = Splicer::moveBlock (R, *v_C, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
								typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
				else
					i = Splicer::moveBlock (R, *v_D, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
								typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
			}
		}
	}

	template <class Container>
	void moveRowSpecialised (const Block &horiz_block, int row, const Container &vert_blocks, VectorRepresentationTypes::Dense01)
	{
		typename Container::const_iterator vert_block;
		typename Matrix1::ConstRowIterator v_X = X.rowBegin () + (horiz_block.sourceIndex () + row);

		if (horiz_block.dest () == 0) {
			typename Matrix2::RowIterator v_A = A.rowBegin () + (horiz_block.destIndex () + row);
			typename Matrix3::RowIterator v_B = B.rowBegin () + (horiz_block.destIndex () + row);

			for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
				typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType v (*v_X, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

				if (vert_block->dest () == 0)
					Splicer::moveBitBlockDense (R, *v_A, v, vert_block->sourceIndex (), vert_block->destIndex ());
				else
					Splicer::moveBitBlockDense (R, *v_B, v, vert_block->sourceIndex (), vert_block->destIndex ());
			}
		} else {
			typename Matrix3::RowIterator v_C = C.rowBegin () + (horiz_block.destIndex () + row);
			typename Matrix3::RowIterator v_D = D.rowBegin () + (horiz_block.destIndex () + row);

			for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
				typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType v (*v_X, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

				if (vert_block->dest () == 0)
					Splicer::moveBitBlockDense (R, *v_C, v, vert_block->sourceIndex (), vert_block->destIndex ());
				else
					Splicer::moveBitBlockDense (R, *v_D, v, vert_block->sourceIndex (), vert_block->destIndex ());
			}
		}
	}

	template <class Container>
	void moveRowSpecialised (const Block &horiz_block, int row, const Container &vert_blocks, VectorRepresentationTypes::Hybrid01)
	{
		typename Container::const_iterator vert_block;
		typename Matrix1::ConstRowIterator v_X = X.rowBegin () + (horiz_block.sourceIndex () + row);
		typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType::const_iterator i;

		typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType vp (*v_X, 0, vert_blocks.back ().sourceIndex () + vert_blocks.back ().size ());
		i = vp.begin ();

		if (horiz_block.dest () == 0) {
			typename Matrix2::RowIterator v_A = A.rowBegin () + (horiz_block.destIndex () + row);
			typename Matrix3::RowIterator v_B = B.rowBegin () + (horiz_block.destIndex () + row);

			for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
				typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType w (*v_X, i, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

				if (vert_block->dest () == 0)
					i = Splicer::moveBitBlockHybrid (R, *v_A, w, vert_block->destIndex ());
				else
					i = Splicer::moveBitBlockHybrid (R, *v_B, w, vert_block->destIndex ());
			}
		} else {
			typename Matrix3::RowIterator v_C = C.rowBegin () + (horiz_block.destIndex () + row);
			typename Matrix3::RowIterator v_D = D.rowBegin () + (horiz_block.destIndex () + row);

			for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
				typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType w (*v_X, i, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

				if (vert_block->dest () == 0)
					i = Splicer::moveBitBlockHybrid (R, *v_C, w, vert_block->destIndex ());
				else
					i = Splicer::moveBitBlockHybrid (R, *v_D, w, vert_block->destIndex ());
			}
		}
	}

public:
	typedef GridTypeRowOptimised GridType;

	MatrixGrid1 (const Ring &__R, Matrix1 &__X, Matrix2 &__A, Matrix3 &__B, Matrix3 &__C, Matrix3 &__D)
		: R (__R), X (__X), A (__A), B (__B), C (__C), D (__D)
		{}

	template <class Container>
	void moveRow (const Block &horiz_block, int row, const Container &vert_blocks)
		{ moveRowSpecialised (horiz_block, row, vert_blocks, typename VectorTraits<Ring, typename Matrix1::Row>::RepresentationType ()); }
};

template <class Ring, class Matrix>
class MatrixGrid2
{
	const Ring &R;
	Matrix &B, &B1, &B2;

public:
	typedef GridTypeNormal GridType;

	MatrixGrid2 (const Ring &__R, Matrix &__B, Matrix &__B1, Matrix &__B2)
		: R (__R), B (__B), B1 (__B1), B2 (__B2)
		{}

	void operator () (const Block &horiz_block, const Block &vert_block)
	{
		if (horiz_block.dest () == 0) {
			if (vert_block.dest () == 0)
				Splicer::copyBlock (R, B, B1, horiz_block, vert_block);
			else
				Splicer::copyBlock (R, B, B2, horiz_block, vert_block);
		}
	}
};

template <class Ring, class Matrix1, class Matrix2>
class MatrixGrid3
{
	const Ring &R;
	Matrix1 &B2, &D2;
	Matrix2 &X;

public:
	typedef GridTypeNormal GridType;

	MatrixGrid3 (const Ring &__R, Matrix1 &__B2, Matrix1 &__D2, Matrix2 &__X)
		: R (__R), B2 (__B2), D2 (__D2), X (__X)
		{}

	void operator () (const Block &horiz_block, const Block &vert_block)
	{
		if (horiz_block.source () == 0) {
			if (vert_block.source () == 0)
				Splicer::copyIdentity (R, X, horiz_block, vert_block);
			else if (vert_block.source () == 2)
				Splicer::copyBlock (R, B2, X, horiz_block, vert_block);
		}
		else if (horiz_block.source () == 1) {
			if (vert_block.source () == 1)
				Splicer::copyIdentity (R, X, horiz_block, vert_block);
			else if (vert_block.source () == 2)
				Splicer::copyBlock (R, D2, X, horiz_block, vert_block);
		}
	}
};

template <class Ring, class Matrix1, class Matrix2>
class MatrixGrid4
{
	const Ring &R;
	Matrix1 &D2;
	Matrix2 &X;

public:
	typedef GridTypeNormal GridType;

	MatrixGrid4 (const Ring &__R, Matrix1 &__D2, Matrix2 &__X)
		: R (__R), D2 (__D2), X (__X)
		{}

	void operator () (const Block &horiz_block, const Block &vert_block)
	{
		if (horiz_block.source () == 0) {
			if (vert_block.source () == 1)
				Splicer::copyIdentity (R, X, horiz_block, vert_block);
			else if (vert_block.source () == 2)
				Splicer::copyBlock (R, D2, X, horiz_block, vert_block);
		}
	}
};

template <class Ring, class Matrix1, class Matrix2>
class MatrixGrid5
{
	const Ring &R;
	Matrix1 &B, &D;
	Matrix2 &X;

public:
	typedef GridTypeNormal GridType;

	MatrixGrid5 (const Ring &__R, Matrix1 &__B, Matrix1 &__D, Matrix2 &__X)
		: R (__R), B (__B), D (__D), X (__X)
		{}

	void operator () (const Block &horiz_block, const Block &vert_block)
	{
		if (horiz_block.source () == 0) {
			if (vert_block.source () == 0)
				Splicer::copyIdentity (R, X, horiz_block, vert_block);
			else
				Splicer::copyBlock (R, B, X, horiz_block, vert_block);
		} else {
			if (vert_block.source () == 1)
				Splicer::copyBlock (R, D, X, horiz_block, vert_block);
		}
	}
};

template <class Ring>
struct DefaultSparseMatrix
{
	typedef SparseMatrix<typename Ring::Element> Type;
};

template <>
struct DefaultSparseMatrix<GF2>
{
	typedef SparseMatrix<bool, Vector<GF2>::Sparse> Type;
};

template <class Ring, class Modules>
template <class Matrix>
void FaugereLachartre<Ring, Modules>::echelonize (Matrix &R, const Matrix &X, size_t &rank, typename Ring::Element &det, bool reduced, bool only_D)
{
	commentator.start ("Reduction of F4-matrix to reduced row-echelon form", __FUNCTION__);

	std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Splicer X_splicer, X_reconst_splicer;

	size_t num_pivot_rows;

	ctx.F.copy (det, ctx.F.one ());

	setup_splicer (X_splicer, X_reconst_splicer, X, num_pivot_rows, det);
	rank = num_pivot_rows;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Found " << num_pivot_rows << " pivots" << std::endl;
	reportUI << "Splicer:" << std::endl << X_splicer << std::endl;

	typename DefaultSparseMatrix<Ring>::Type A (num_pivot_rows, num_pivot_rows);
	// DenseMatrix<typename Ring::Element> A (num_pivot_rows, num_pivot_rows);
	DenseMatrix<typename Ring::Element> B (num_pivot_rows, X.coldim () - num_pivot_rows);
	DenseMatrix<typename Ring::Element> C (X.rowdim () - num_pivot_rows, num_pivot_rows);
	DenseMatrix<typename Ring::Element> D (X.rowdim () - num_pivot_rows, X.coldim () - num_pivot_rows);
	// typename DefaultSparseMatrix<Ring>::Type B (num_pivot_rows, X.coldim () - num_pivot_rows);
	// typename DefaultSparseMatrix<Ring>::Type C (X.rowdim () - num_pivot_rows, num_pivot_rows);
	// typename DefaultSparseMatrix<Ring>::Type D (X.rowdim () - num_pivot_rows, X.coldim () - num_pivot_rows);

	X_splicer.splice (MatrixGrid1<Ring, const Matrix, typename DefaultSparseMatrix<Ring>::Type, DenseMatrix<typename Ring::Element> > (ctx.F, X, A, B, C, D));
	// X_splicer.splice (MatrixGrid1<Ring, const Matrix, DenseMatrix<typename Ring::Element>, DenseMatrix<typename Ring::Element> > (ctx.F, X, A, B, C, D));
	// X_splicer.splice (MatrixGrid1<Ring, const Matrix, typename DefaultSparseMatrix<Ring>::Type, typename DefaultSparseMatrix<Ring>::Type > (ctx.F, X, A, B, C, D));

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

	commentator.start ("Constructing A^-1 B");

	BLAS3::trsm (ctx, ctx.F.one (), A, B, UpperTriangular, false);

	commentator.stop (MSG_DONE);

	commentator.start ("Constructing D - C A^-1 B");

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

	EF.echelonize (D);

	reportUI << "Row-echelon form of D - C A^-1 B:" << std::endl;
	BLAS3::write (ctx, reportUI, D);

	// commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	// 	<< "Rank of dense part is " << r_D << std::endl;

	Splicer D_splicer, D_reconst_splicer;
	Splicer composed_splicer, subst_splicer, D_splicer_rev;

	setup_splicer (D_splicer, D_reconst_splicer, D, num_pivot_rows, det);
	rank += num_pivot_rows;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "(In D) found " << num_pivot_rows << " pivots" << std::endl;
	reportUI << "Splicer:" << std::endl << D_splicer << std::endl;

	X_reconst_splicer.compose (subst_splicer, D_reconst_splicer, 1, Splicer::noSource, 0, Splicer::noSource);
	subst_splicer.removeGaps ();
	subst_splicer.consolidate ();

	reportUI << "Splicer after substitution:" << std::endl << subst_splicer << std::endl;

	BLAS3::scal (ctx, ctx.F.zero (), R);

	if (reduced) {
		DenseMatrix<typename Ring::Element> D1 (num_pivot_rows, num_pivot_rows);
		DenseMatrix<typename Ring::Element> D2 (num_pivot_rows, D.coldim () - num_pivot_rows);

		D_splicer.splice (MatrixGrid2<Ring, DenseMatrix<typename Ring::Element> > (ctx.F, D, D1, D2));

		reportUI << "Matrix D1:" << std::endl;
		BLAS3::write (ctx, reportUI, D1);
		reportUI << "Matrix D2:" << std::endl;
		BLAS3::write (ctx, reportUI, D2);

		commentator.start ("Constructing D1^-1 D2");

		BLAS3::trsm (ctx, ctx.F.one (), D1, D2, UpperTriangular, false);

		commentator.stop (MSG_DONE);

		D_splicer.reverse (D_splicer_rev);

		if (only_D) {
			subst_splicer.compose (composed_splicer, D_splicer_rev, 1, 1, 1, Splicer::noSource);
			composed_splicer.clearHorizontalBlocks ();
			composed_splicer.addHorizontalBlock (Block (0, 0, 0, 0, D1.rowdim ()));
			R.resize (D1.rowdim (), R.coldim ());
			composed_splicer.splice (MatrixGrid4<Ring, DenseMatrix<typename Ring::Element>, Matrix> (ctx.F, D2, R));

			reportUI << "Composed splicer:" << std::endl << composed_splicer << std::endl;
		} else {
			DenseMatrix<typename Ring::Element> B1 (B.rowdim (), num_pivot_rows);
			DenseMatrix<typename Ring::Element> B2 (B.rowdim (), D.coldim () - num_pivot_rows);

			Splicer B_splicer (D_splicer);
			B_splicer.clearHorizontalBlocks ();
			B_splicer.addHorizontalBlock (Block (0, 0, 0, 0, B.rowdim ()));

			B_splicer.splice (MatrixGrid2<Ring, DenseMatrix<typename Ring::Element> > (ctx.F, B, B1, B2));

			reportUI << "Matrix B1:" << std::endl;
			BLAS3::write (ctx, reportUI, B1);
			reportUI << "Matrix B2:" << std::endl;
			BLAS3::write (ctx, reportUI, B2);

			commentator.start ("Constructing B2 - B1 D1^-1 D2");

			BLAS3::gemm (ctx, ctx.F.minusOne (), B1, D2, ctx.F.one (), B2);

			commentator.stop (MSG_DONE);

			reportUI << "B2 - B1 D1^-1 D2:" << std::endl;
			BLAS3::write (ctx, reportUI, B2);

			subst_splicer.compose (composed_splicer, D_splicer_rev, 1, 1);
			composed_splicer.fillHorizontal (2, 0, X.rowdim ());
			composed_splicer.splice (MatrixGrid3<Ring, DenseMatrix<typename Ring::Element>, Matrix> (ctx.F, B2, D2, R));

			reportUI << "Composed splicer:" << std::endl << composed_splicer << std::endl;
		}
	} else {
		if (only_D) {
			subst_splicer.clearHorizontalBlocks ();
			subst_splicer.addHorizontalBlock (Block (0, 0, 0, 0, num_pivot_rows));
			R.resize (num_pivot_rows, R.coldim ());
			subst_splicer.splice (MatrixGrid4<Ring, DenseMatrix<typename Ring::Element>, Matrix> (ctx.F, D, R));
		} else {
			subst_splicer.splice (MatrixGrid5<Ring, DenseMatrix<typename Ring::Element>, Matrix> (ctx.F, B, D, R));
		}
	}

	commentator.stop (MSG_DONE, NULL, __FUNCTION__);
}

} // namespace LELA

#endif // __LELA_ALGORITHMS_FAUGERE_LACHARTRE_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
