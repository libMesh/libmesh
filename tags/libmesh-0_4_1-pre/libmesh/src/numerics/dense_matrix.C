// $Id: dense_matrix.C,v 1.13 2003-09-02 18:02:43 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Local Includes
#include "dense_matrix.h"




// // ------------------------------------------------------------
// // Dense Matrix member functions
template<typename T>
void DenseMatrix<T>::left_multiply (const DenseMatrixBase<T>& M2)
{
  // (*this) <- M2 * M3
  // Where:
  // (*this) = (m x n),
  // M2      = (m x p),
  // M3      = (p x n)

  // M3 is a copy of *this before it gets resize()d
  DenseMatrix<T> M3(*this);

  // Resize *this so that the result can fit
  this->resize (M2.m(), M3.n());

  // Call the multiply function in the base class
  this->multiply(*this, M2, M3);
}




template<typename T>
void DenseMatrix<T>::left_multiply_transpose(const DenseMatrix<T>& A)
{

  DenseMatrix<T> B(*this);
  
  this->resize (A.n(), B.n());
      
  assert (A.m() == B.m());
  assert (this->m() == A.n());
  assert (this->n() == B.n());
      
  const unsigned int m_s = A.n();
  const unsigned int p_s = A.m(); 
  const unsigned int n_s = this->n();
  
  // Do it this way because there is a
  // decent chance (at least for constraint matrices)
  // that A.transpose(i,k) = 0.
  for (unsigned int i=0; i<m_s; i++)
    for (unsigned int k=0; k<p_s; k++)
      if (A.transpose(i,k) != 0.)
	for (unsigned int j=0; j<n_s; j++)
	  (*this)(i,j) += A.transpose(i,k)*B(k,j);  

}






template<typename T>
void DenseMatrix<T>::right_multiply (const DenseMatrixBase<T>& M3)
{
  // (*this) <- M2 * M3
  // Where:
  // (*this) = (m x n),
  // M2      = (m x p),
  // M3      = (p x n)

  // M2 is a copy of *this before it gets resize()d
  DenseMatrix<T> M2(*this);

  // Resize *this so that the result can fit
  this->resize (M2.m(), M3.n());
  
  this->multiply(*this, M2, M3);
}




template<typename T>
void DenseMatrix<T>::right_multiply_transpose (const DenseMatrix<T>& B)
{
  DenseMatrix<T> A(*this);
  
  this->resize (A.m(), B.m());
      
  assert (A.n() == B.n());
  assert (this->m() == A.m());
  assert (this->n() == B.m());
      
  const unsigned int m_s = A.m();
  const unsigned int p_s = A.n(); 
  const unsigned int n_s = this->n();

  // Do it this way because there is a
  // decent chance (at least for constraint matrices)
  // that B.transpose(k,j) = 0.
  for (unsigned int j=0; j<n_s; j++)
    for (unsigned int k=0; k<p_s; k++)
      if (B.transpose(k,j) != 0.)
	for (unsigned int i=0; i<m_s; i++)
	  (*this)(i,j) += A(i,k)*B.transpose(k,j);	            
}




//--------------------------------------------------------------
// Explicit instantiations
template class DenseMatrix<Real>;

#ifdef USE_COMPLEX_NUMBERS
template class DenseMatrix<Complex>;
#endif
