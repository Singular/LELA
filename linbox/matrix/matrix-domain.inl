/* linbox/matrix/matrix-domain.inl
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_domain_INL
#define __LINBOX_matrix_domain_INL

#include "linbox/matrix/transpose-matrix.h"

namespace LinBox
{

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvRowSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::DenseVectorTag) const
{
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector2::iterator j = y.begin ();

	typename Field::Element d;

            // JGD 02.09.2008 : when sizes differ
            // A must decide if dot is possible, not w 
// 	for (; j != w.end (); ++j, ++i)
// 		_VD.dot (*j, v, *i);
	for (; i != A.rowEnd (); ++j, ++i) {
		_VD.dot (d, x, *i);
		_F.mulin (*j, beta);
		_F.axpyin (*j, alpha, d);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvRowSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseSequenceVectorTag) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	std::vector<std::pair<size_t, typename Field::Element> > yp;

	if (_F.isZero (beta))
		y.clear ();
	else
		_VD.mulin (y, beta);

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, x, *i);
		_F.mulin (t, alpha);

		if (!_F.isZero (t))
			yp.push_back (std::pair<size_t, typename Field::Element> (idx, t));
	}

	return _VD.addin (y, yp);
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvRowSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseAssociativeVectorTag) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	if (_F.isZero (beta))
		y.clear ();
	else
		_VD.mulin (y, beta);

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, x, *i);
		_F.mulin (t, alpha);

		if (!_F.isZero (t))
			y[idx] += t;
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvRowSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseParallelVectorTag) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	std::pair<std::vector<size_t>, std::vector<typename Field::Element> > yp;

	if (_F.isZero (beta)) {
		y.first.clear ();
		y.second.clear ();
	}
	else
		_VD.mulin (y, beta);

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, x, *i);
		_F.mulin (t, alpha);

		if (!_F.isZero (t)) {
			yp.first.push_back (idx);
			yp.second.push_back (t);
		}
	}

	return _VD.addin (y, yp);
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Field>::gemvColDense (const VectorDomain<Field>     &VD,
					       const typename Field::Element &alpha,
					       const Matrix                  &A,
					       const Vector1                 &x,
					       const typename Field::Element &beta,
					       Vector2                       &y) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector1::const_iterator j = x.begin ();
	typename Field::Element d;

	VD.mulin (y, beta);

	for (; j != x.end (); ++j, ++i) {
		VD.field ().mul (d, alpha, *j);
		VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvColSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseSequenceVectorTag,
								VectorCategories::DenseVectorTag) const
{
	linbox_check (A.rowdim () == y.size ());

	typename Vector1::const_iterator j = x.begin ();
	typename Field::Element d;

	_VD.mulin (y, beta);

	for (; j != x.end (); ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		_F.mul (d, alpha, j->second);
		_VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvColSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta, 
								Vector2                       &y,
								VectorCategories::SparseAssociativeVectorTag,
								VectorCategories::DenseVectorTag) const
{
	linbox_check (A.rowdim () == y.size ());

	typename Vector1::const_iterator j = x.begin ();
	typename Field::Element d;

	_VD.mulin (y, beta);

	for (; j != x.end (); ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		_F.mul (d, alpha, j->second);
		_VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvColSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseParallelVectorTag,
								VectorCategories::DenseVectorTag) const
{
	linbox_check (A.rowdim () == y.size ());

	typename Vector1::first_type::const_iterator j_idx = x.first.begin ();
	typename Vector1::second_type::const_iterator j_elt = x.second.begin ();
	typename Field::Element d;

	_VD.mulin (y, beta);

	for (; j_idx != x.first.end (); ++j_idx, ++j_elt) {
		typename Matrix::ConstColIterator i = A.colBegin () + *j_idx;
		_F.mul (d, alpha, *j_elt);
		_VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvColSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseParallelVectorTag,
								VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.rowdim () == y.size ());

	typename Vector1::first_type::const_iterator j_idx = x.first.begin ();
	typename Vector1::second_type::const_iterator j_elt = x.second.begin ();
	typename Field::Element d;

	_VD.mulin (y, beta);

	for (; j_idx != x.first.end (); ++j_idx, ++j_elt) {
		typename Matrix::ConstColIterator i = A.colBegin () + *j_idx;
		_F.mul (d, alpha, *j_elt);
		_VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Matrix, class Vector>
Vector &MatrixDomainSupportGeneric<Field>::trsvSpecialized (const Matrix &A, Vector &x,
							    MatrixCategories::RowMatrixTag,
							    VectorCategories::DenseVectorTag)
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (A.rowdim () == x.size ());

	typename Field::Element ai, ai_p_1, neg_ai_inv, d;
	int i = A.rowdim () - 1;

	while (--i >= 0) {
		if (_VD.firstNonzeroEntry (ai, *(A.rowBegin () + i)) == -1)
			continue;

		_VD.dot (d, *(A.rowBegin () + i), x);

		_F.add (ai_p_1, ai, _one);
		_F.mulin (x[i], ai_p_1);
		_F.inv (neg_ai_inv, ai);
		_F.negin (neg_ai_inv);
		_F.axpyin (x[i], neg_ai_inv, d);
	}

	return x;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomainSupportGeneric<Field>::copyRow (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.copy (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomainSupportGeneric<Field>::copyCol (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.copy (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
bool MatrixDomainSupportGeneric<Field>::areEqualRow (const Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		if (!_VD.areEqual (*i, *j))
			return false;

	return true;
}

template <class Field>
template <class Matrix1, class Matrix2>
bool MatrixDomainSupportGeneric<Field>::areEqualCol (const Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		if (!_VD.areEqual (*i, *j))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
bool MatrixDomainSupportGeneric<Field>::isZeroRow (const Matrix &A) const
{
	typename Matrix::ConstRowIterator i;

	i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
bool MatrixDomainSupportGeneric<Field>::isZeroCol (const Matrix &A) const
{
	typename Matrix::ConstColIterator i;

	i = A.colBegin ();

	for (; i != A.colEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomainSupportGeneric<Field>::scalRow (Matrix &B, const typename Field::Element &a) const
{
	typename Matrix::RowIterator i;

	for (i = B.rowBegin (); i != B.rowEnd (); ++i)
		_VD.mulin (*i, a);

	return B;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomainSupportGeneric<Field>::scalCol (Matrix &B, const typename Field::Element &a) const
{
	typename Matrix::ColIterator i;

	for (i = B.colBegin (); i != B.colEnd (); ++i)
		_VD.mulin (*i, a);

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::axpyRow (const typename Field::Element &a,
						     const Matrix1                 &A,
						     Matrix2                       &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::RowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.axpyin (*j, a, *i);

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::axpyCol (const typename Field::Element &a,
						     const Matrix1                 &A,
						     Matrix2                       &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.axpyin (*j, a, *i);

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmRowRowCol (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::RowIterator l1;
	typename Matrix3::Row::iterator l2;
	typename Field::Element d;

	for (i = A.rowBegin (), l1 = C.rowBegin (); i != A.rowEnd (); ++i, ++l1) {
		for (j = B.colBegin (), l2 = l1->begin (); j != B.colEnd (); ++j, ++l2) {
			_VD.dot (d, *i, *j);
			_F.mulin (*l2, beta);
			_F.axpyin (*l2, alpha, d);
		}
	}

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmColRowCol (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::ColIterator l1;
	typename Matrix3::Col::iterator l2;
	typename Field::Element d;

	for (j = B.colBegin (), l1 = C.colBegin (); j != B.colEnd (); ++j, ++l1) {
		for (i = A.rowBegin (), l2 = l1->begin (); i != A.rowEnd (); ++i, ++l2) {
			_VD.dot (d, *i, *j);
			_F.mulin (*l2, beta);
			_F.axpyin (*l2, alpha, d);
		}
	}

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmRowRowRow (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i = A.rowBegin ();
	typename Matrix3::RowIterator j = C.rowBegin ();

	TransposeMatrix<const Matrix2> BT (B);

	for (; i != A.rowEnd (); ++i, ++j)
		gemv (alpha, BT, *i, beta, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmColColCol (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstColIterator i = B.colBegin ();
	typename Matrix3::ColIterator j = C.colBegin ();

	for (; i != B.colEnd (); ++i, ++j)
		gemv (alpha, A, *i, beta, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::trsmSpecialized (const typename Field::Element &alpha, const Matrix1 &A, Matrix2 &B,
							     MatrixCategories::RowMatrixTag,
							     MatrixCategories::RowMatrixTag)
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	typename Field::Element ai, ai_p_1, neg_ai_inv, d;
	int i = A.rowdim () - 1;

	TransposeMatrix<const Matrix1> AT (A);

	while (--i >= 0) {
		if (_VD.firstNonzeroEntry (ai, *(A.rowBegin () + i)) == -1)
			continue;

		_F.add (ai_p_1, ai, _one);
		_F.inv (neg_ai_inv, ai);
		_F.negin (neg_ai_inv);

		gemv (neg_ai_inv, AT, *(A.rowBegin () + i), ai_p_1, *(B.rowBegin () + i));
	}

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::trsmSpecialized (const typename Field::Element &alpha, const Matrix1 &A, Matrix2 &B,
							     MatrixCategories::RowMatrixTag,
							     MatrixCategories::ColMatrixTag)
{
	typename Matrix2::ColIterator i_B;

	for (i_B = B.colBegin (); i_B != B.colEnd (); ++i_B) {
		trsv (A, *i_B);
		_VD.mulin (*i_B, alpha);
	}

	return B;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteRowsByRow (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	Iterator i;
	typename Matrix::RowIterator j, k;

	for (i = P_start; i != P_end; ++i) {
		j = A.rowBegin () + i->first;
		k = A.rowBegin () + i->second;

		_VD.swap (*j, *k);
	}

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteRowsByCol (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	typename Matrix::ColIterator j;

	for (j = A.colBegin (); j != A.colEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteColsByRow (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	typename Matrix::RowIterator j;

	for (j = A.rowBegin (); j != A.rowEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteColsByCol (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	Iterator i;
	typename Matrix::ColIterator j, k;

	for (i = P_start; i != P_end; ++i) {
		j = A.colBegin () + i->first;
		k = A.colBegin () + i->second;

		_VD.swap (*j, *k);
	}

	return A;
}

// FIXME: Possibly useless methods which should be dropped, I think

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomain<Field>::leftMulin (const Matrix1 &A, Matrix2 &B) const
{
	linbox_check (A.rowdim () == A.coldim ());
	linbox_check (A.coldim () == B.rowdim ());

	typename LinBox::Vector<Field>::Dense t (A.rowdim ());

	typename Matrix2::ColIterator i;
	typename Matrix1::ConstRowIterator j;

	typename LinBox::Vector<Field>::Dense::iterator k;

	for (i = B.colBegin (); i != B.colEnd (); ++i) {
		for (j = A.rowBegin (), k = t.begin (); j != A.rowEnd (); ++j, ++k)
			MatrixDomainSupportGeneric<Field>::_VD.dot (*k, *j, *i);

		MatrixDomainSupportGeneric<Field>::_VD.copy (*i, t);
	}

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomain<Field>::rightMulin (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (B.rowdim () == B.coldim ());

	typename LinBox::Vector<Field>::Dense t (B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstColIterator j;

	typename LinBox::Vector<Field>::Dense::iterator k;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i) {
		for (j = B.colBegin (), k = t.begin (); j != B.colEnd (); ++j, ++k)
			MatrixDomainSupportGeneric<Field>::_VD.dot (*k, *i, *j);

		MatrixDomainSupportGeneric<Field>::_VD.copy (*i, t);
	}

	return A;
}

/* FIXME: These methods are undocumented, and I'm unclear what they are supposed
 * to do 
 */

#if 0

/*M1<-M2**k;
 */
template<class Matrix1, class Matrix2>
Matrix1& MatrixDomain::pow_apply(Matrix1& M1, const Matrix2& M2, unsigned long int k) const
{
	linbox_check((M1.rowdim()==M1.coldim())&&
		     (M2.rowdim()==M2.coldim())&&
		     (M1.rowdim()==M2.rowdim()));

  
	typename Matrix1::RawIterator p=M1.rawBegin();
	for(;p!=M1.rawEnd();++p)
		M1.field().init(*p,0);
	for(p=M1.rawBegin();p<M1.rawEnd();)
	{
		M1.field().init(*p,1);
		p=p+M1.rowdim()+1;
	}
    
    
	for(int i=0;i<k;++i)
		mulin_R(M1,M2);
    
	return M1;
}
    
  
template<class Matrix1, class Matrix2>
Matrix1& MatrixDomain::pow_horn(Matrix1& M1, const Matrix2& M2, unsigned long int k) const
{
	linbox_check((M1.rowdim()==M1.coldim())&&
		     (M2.rowdim()==M2.coldim())&&
		     (M1.rowdim()==M2.rowdim()));
    
	if(k==0)
	{
		typename Matrix1::RawIterator p=M1.rawBegin();
		for(;p!=M1.rawEnd();++p)
			M1.field().init(*p,0);
		for(p=M1.rawBegin();p<M1.rawEnd();)
		{
			M1.field().init(*p,1);
			p+=M1.rowdim()+1;
		}
		return M1;
	}
    
	typename Matrix1::RawIterator p1;
	typename Matrix2::ConstRawIterator p2;
	for(p1=M1.rawBegin(),p2=M2.rawBegin();p1!=M1.rawEnd();++p1,++p2)
		M1.field().assign(*p1,*p2);
  
	std::vector<bool> bit;
	bit.reserve(sizeof(unsigned long)*4);
	while(k>0)
	{
		bit.push_back(k%2);
		k/=2;
	};

    
	std::vector<bool>::reverse_iterator p=bit.rbegin();
	++p;
	Matrix1 temp(M1);
	for(;p!=bit.rend();++p)
	{
		temp=M1;
		mulin_L(M1,temp);
		if(*p)
			mulin_L(M1,M2);

	}
      
	return M1;     
}

#endif

} // namespace LinBox

#endif // __LINBOX_matrix_domain_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
