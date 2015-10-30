// $Id: dense_matrix.C,v 1.15 2004-01-03 15:37:43 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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

#include "libmesh.h"



// ------------------------------------------------------------
// Dense Matrix member functions
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



template<typename T>
void DenseMatrix<T>::lu_solve (DenseVector<T>& b,
			       DenseVector<T>& x,
			       const bool partial_pivot)
{
  // Perform LU decomposition
  this->lu_decompose (partial_pivot);

  // Perform back substitution
  this->lu_back_substitute (b, x, partial_pivot);
}

  

template<typename T>
void DenseMatrix<T>::lu_back_substitute (DenseVector<T>& b,
					 DenseVector<T>& x,
					 const bool ) const
{
  const unsigned int
    m = this->m(),
    n = this->n();

  assert (m == n);
  assert (b.size() == m);
  
  x.resize (n);

  
  // Transform the RHS vector
  for (unsigned int i=0; i<(n-1); i++)
    {
      // Get the diagonal entry and take its inverse
      const T diag = (*this)(i,i);
      
      assert (diag != libMesh::zero);
  
      const T diag_inv = 1./diag;

      // Get the entry b(i) and store it
      const T bi = b(i);
      
      for (unsigned int j=(i+1); j<n; j++)
	b(j) -= bi*(*this)(j,i)*diag_inv;
    }


  // Perform back-substitution
  {
    x(n-1) = b(n-1)/(*this)(n-1,n-1);

    for (unsigned int i=0; i<=(n-1); i++)
      {
	const unsigned int ib = (n-1)-i;

	// Get the diagonal and take its inverse
	const T diag = (*this)(ib,ib);

	assert (diag != libMesh::zero);

	const T diag_inv = 1./diag;
	
	for (unsigned int j=(ib+1); j<n; j++)
	  {
	    b(ib) -= (*this)(ib,j)*x(j);
	    x(ib)  = b(ib)*diag_inv;
	  }
      }
  }
}


template<typename T>
void DenseMatrix<T>::lu_decompose (const bool partial_pivot)
{

  // Get the matrix size and make sure it is square
  const unsigned int
    m = this->m(),
    n = this->n();

  assert (m == n);


  
  // Straight, vanilla LU factorization without pivoting
  if (!partial_pivot)
    {
      
      // For each row in the matrix
      for (unsigned int i=0; i<m; i++)
	{
	  // Get the diagonal entry and take its inverse
	  const T diag = (*this)(i,i);

	  assert (diag != libMesh::zero);

	  const T diag_inv = 1./diag;
	  
	  // For each row in the submatrix
	  for (unsigned int j=i+1; j<m; j++)
	    {
	      // Get the scale factor for this row
	      const T fact = (*this)(j,i)*diag_inv;
	      
	      // For each column in the subrow scale it
	      // by the factor
	      for (unsigned int k=i+1; k<m; k++)
		(*this)(j,k) -= fact*(*this)(i,k);	      
	    }
	}      
    }

  // Do partial pivoting.
  else
    {
      error();
    }
}



//--------------------------------------------------------------
// Explicit instantiations
template class DenseMatrix<Real>;

#ifdef USE_COMPLEX_NUMBERS
template class DenseMatrix<Complex>;
#endif
