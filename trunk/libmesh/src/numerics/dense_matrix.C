// $Id: dense_matrix.C,v 1.24 2005-06-07 12:52:21 spetersen Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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


// C++ Includes
#include <cmath> // for sqrt

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
  // Check for a previous decomposition
  switch(this->_decomposition_type)
    {
    case NONE:
      {
	this->_lu_decompose (partial_pivot);
	break;
      }

    case LU:
      {
	// Already factored, just need to call back_substitute.
	break;
      }

    default:
      {
	std::cerr << "Error! This matrix already has a "
		  << "different decomposition..."
		  << std::endl;
	error();
      }
    }

  // Perform back substitution
  this->_lu_back_substitute (b, x, partial_pivot);
}

  

template<typename T>
void DenseMatrix<T>::_lu_back_substitute (DenseVector<T>& b,
					  DenseVector<T>& x,
					  const bool ) const
{
  const unsigned int
    m = this->m(),
    n = this->n();

  assert (m == n);
  assert (b.size() == m);
  
  x.resize (n);

  // A convenient reference to *this
  const DenseMatrix<T>& A = *this;

  // Transform the RHS vector
  for (unsigned int i=0; i<(n-1); i++)
    {
      // Get the diagonal entry and take its inverse
      const T diag = A(i,i);
      
      assert (diag != libMesh::zero);
  
      const T diag_inv = 1./diag;

      // Get the entry b(i) and store it
      const T bi = b(i);
      
      for (unsigned int j=(i+1); j<n; j++)
	b(j) -= bi*A(j,i)*diag_inv;
    }


  // Perform back-substitution
  {
    x(n-1) = b(n-1)/A(n-1,n-1);

    for (unsigned int i=0; i<=(n-1); i++)
      {
	const unsigned int ib = (n-1)-i;

	// Get the diagonal and take its inverse
	const T diag = A(ib,ib);

	assert (diag != libMesh::zero);

	const T diag_inv = 1./diag;
	
	for (unsigned int j=(ib+1); j<n; j++)
	  {
	    b(ib) -= A(ib,j)*x(j);
	    x(ib)  = b(ib)*diag_inv;
	  }
      }
  }
}


template<typename T>
void DenseMatrix<T>::_lu_decompose (const bool partial_pivot)
{
  // If this function was called, there better not be any
  // previous decomposition of the matrix.
  assert(this->_decomposition_type == NONE);
  
  // Get the matrix size and make sure it is square
  const unsigned int
    m = this->m(),
    n = this->n();

  assert (m == n);

  // A convenient reference to *this
  DenseMatrix<T>& A = *this;
  
  // Straight, vanilla LU factorization without pivoting
  if (!partial_pivot)
    {
      
      // For each row in the matrix
      for (unsigned int i=0; i<m; i++)
	{
	  // Get the diagonal entry and take its inverse
	  const T diag = A(i,i);

	  assert (diag != libMesh::zero);

	  const T diag_inv = 1./diag;
	  
	  // For each row in the submatrix
	  for (unsigned int j=i+1; j<m; j++)
	    {
	      // Get the scale factor for this row
	      const T fact = A(j,i)*diag_inv;
	      
	      // For each column in the subrow scale it
	      // by the factor
	      for (unsigned int k=i+1; k<m; k++)
		A(j,k) -= fact*A(i,k);	      
	    }
	}
    }
  
  // Do partial pivoting.
  else
    {
      error();
    }
  
  // Set the flag for LU decomposition
  this->_decomposition_type = LU;
}



template<typename T>
T DenseMatrix<T>::det ()
{
  // First LU decompose the matrix (without partial pivoting).
  // Note that the lu_decompose routine will check to see if the
  // matrix is square so we don't worry about it.
  if (this->_decomposition_type == NONE)
    this->_lu_decompose(false);
  else if (this->_decomposition_type != LU)
    {
      std::cerr << "Error! Can't compute the determinant under "
		<< "the current decomposition."
		<< std::endl;
      error();
    }
  
  // A variable to keep track of the running product of diagonal terms.
  T determinant = 1.;
  
  // Loop over diagonal terms, computing the product.  In practice,
  // be careful because this value could easily become too large to
  // fit in a double or float.  To be safe, one should keep track of
  // the power (of 10) of the determinant in a separate variable
  // and maintain an order 1 value for the determinant itself.
  for (unsigned int i=0; i<this->m(); i++)
    determinant *= (*this)(i,i);

  // Return the determinant
  return determinant;
}



// The cholesky solve function first decomposes the matrix
// with cholesky_decompose and then uses the cholesky_back_substitute
// routine to find the solution x.
template <typename T>
template <typename T2>
void DenseMatrix<T>::cholesky_solve (DenseVector<T2>& b,
				     DenseVector<T2>& x)
{
  // Check for a previous decomposition
  switch(this->_decomposition_type)
    {
    case NONE:
      {
	this->_cholesky_decompose ();
	break;
      }

    case CHOLESKY:
      {
	// Already factored, just need to call back_substitute.
	break;
      }
      
    default:
      {
	std::cerr << "Error! This matrix already has a "
		  << "different decomposition..."
		  << std::endl;
	error();
      }
    }

  // Perform back substitution
  this->_cholesky_back_substitute (b, x);
}




// This algorithm is based on the Cholesky decomposition in
// the Numerical Recipes in C book.
template<typename T>
void DenseMatrix<T>::_cholesky_decompose ()
{
  // If we called this function, there better not be any
  // previous decomposition of the matrix.
  assert(this->_decomposition_type == NONE);
  
  // Shorthand notation for number of rows and columns.
  const unsigned int
    m = this->m(),
    n = this->n();

  // Just to be really sure...
  assert(m==n);

  // A convenient reference to *this
  DenseMatrix<T>& A = *this;
  
  for (unsigned int i=0; i<m; ++i)
    {
      for (unsigned int j=i; j<n; ++j)
	{
	  for (unsigned int k=0; k<i; ++k)
	    A(i,j) -= A(i,k) * A(j,k);

	  if (i == j)
	    {
#ifndef USE_COMPLEX_NUMBERS
	      if (A(i,j) <= 0.0)
		{
		  std::cerr << "Error! Can only use Cholesky decomposition "
			    << "with symmetric positive definite matrices."
			    << std::endl;
		  error();
		}
#endif

	      A(i,i) = std::sqrt(A(i,j));
	    }
	  else
	    A(j,i) = A(i,j) / A(i,i);
	}
    }
  
  // Set the flag for CHOLESKY decomposition
  this->_decomposition_type = CHOLESKY;
}



template <typename T>
template <typename T2>
void DenseMatrix<T>::_cholesky_back_substitute (DenseVector<T2>& b,
						DenseVector<T2>& x) const
{
  // Shorthand notation for number of rows and columns.
  const unsigned int
    m = this->m(),
    n = this->n();

  // Just to be really sure...
  assert(m==n);

  // A convenient reference to *this
  const DenseMatrix<T>& A = *this;
   
  // Now compute the solution to Ax =b using the factorization.
  x.resize(m);
  
  // Solve for Ly=b
  for (unsigned int i=0; i<n; ++i)
    {
      for (unsigned int k=0; k<i; ++k)
	b(i) -= A(i,k)*x(k);
      
      x(i) = b(i) / A(i,i);
    }
  
  // Solve for L^T x = y
  for (unsigned int i=0; i<n; ++i)
    {
      const unsigned int ib = (n-1)-i;

      for (unsigned int k=(ib+1); k<n; ++k)
	x(ib) -= A(k,ib) * x(k);
      
      x(ib) /= A(ib,ib);
    }
}



// This routine is commented out since it is not really a memory
// efficient implementation.  Also, you don't *need* the inverse
// for anything, instead just use lu_solve to solve Ax=b.
// template<typename T>
// void DenseMatrix<T>::inverse ()
// {
//   // First LU decompose the matrix (without partial pivoting).
//   // Note that the lu_decompose routine will check to see if the
//   // matrix is square so we don't worry about it.
//   if (!this->_lu_decomposed)
//     this->_lu_decompose(false);

//   // A unit vector which will be used as a rhs
//   // to pick off a single value each time.
//   DenseVector<T> e;
//   e.resize(this->m());

//   // An empty vector which will be used to hold the solution
//   // to the back substitutions.
//   DenseVector<T> x;
//   x.resize(this->m());

//   // An empty dense matrix to store the resulting inverse
//   // temporarily until we can overwrite A.
//   DenseMatrix<T> inv;
//   inv.resize(this->m(), this->n());
  
//   // Resize the passed in matrix to hold the inverse
//   inv.resize(this->m(), this->n());
  
//   for (unsigned int j=0; j<this->n(); ++j)
//     {
//       e.zero();
//       e(j) = 1.;
//       this->_lu_back_substitute(e, x, false);
//       for (unsigned int i=0; i<this->n(); ++i)
// 	inv(i,j) = x(i);
//     }

//   // Now overwrite all the entries
//   *this = inv;
// }


//--------------------------------------------------------------
// Explicit instantiations
template class DenseMatrix<Real>;
template void DenseMatrix<Real>::cholesky_solve(DenseVector<Real>&, DenseVector<Real>&);
template void DenseMatrix<Real>::_cholesky_back_substitute(DenseVector<Real>&, DenseVector<Real>&);

#ifdef USE_COMPLEX_NUMBERS
template class DenseMatrix<Complex>;
template void DenseMatrix<Complex>::cholesky_solve(DenseVector<Complex>&,DenseVector<Complex>&);
template void DenseMatrix<Complex>::_cholesky_back_substitute(DenseVector<Complex>&, DenseVector<Complex>&);
template void DenseMatrix<Real>::cholesky_solve(DenseVector<Complex>&, DenseVector<Complex>&);
template void DenseMatrix<Real>::_cholesky_back_substitute(DenseVector<Complex>&, DenseVector<Complex>&);
#endif
