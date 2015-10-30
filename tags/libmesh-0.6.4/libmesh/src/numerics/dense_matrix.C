// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "dense_vector.h"
#include "libmesh.h"


#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)
#include "petsc_macro.h"

EXTERN_C_FOR_PETSC_BEGIN
#include <petscblaslapack.h>
EXTERN_C_FOR_PETSC_END
#endif

// ------------------------------------------------------------
// Dense Matrix member functions



template<typename T>
void DenseMatrix<T>::left_multiply (const DenseMatrixBase<T>& M2)
{
  if (this->use_blas)
    this->_multiply_blas(M2, LEFT_MULTIPLY);
  else
    {
      // (*this) <- M2 * (*this)
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
}




template<typename T>
void DenseMatrix<T>::left_multiply_transpose(const DenseMatrix<T>& A)
{
  if (this->use_blas)
    this->_multiply_blas(A, LEFT_MULTIPLY_TRANSPOSE);
  else
    {
      //Check to see if we are doing (A^T)*A
      if (this == &A)
	{
	  //libmesh_here();
	  DenseMatrix<T> B(*this);
      
	  // Simple but inefficient way
	  // return this->left_multiply_transpose(B);

	  // More efficient, but more code way
	  // If A is mxn, the result will be a square matrix of Size n x n.
	  const unsigned int m = A.m();
	  const unsigned int n = A.n();

	  // resize() *this and also zero out all entries.
	  this->resize(n,n);

	  // Compute the lower-triangular part
	  for (unsigned int i=0; i<n; ++i)
	    for (unsigned int j=0; j<=i; ++j)
	      for (unsigned int k=0; k<m; ++k) // inner products are over m
		(*this)(i,j) += B(k,i)*B(k,j);

	  // Copy lower-triangular part into upper-triangular part
	  for (unsigned int i=0; i<n; ++i)
	    for (unsigned int j=i+1; j<n; ++j)
	      (*this)(i,j) = (*this)(j,i);
	}

      else
	{
	  DenseMatrix<T> B(*this);
  
	  this->resize (A.n(), B.n());
      
	  libmesh_assert (A.m() == B.m());
	  libmesh_assert (this->m() == A.n());
	  libmesh_assert (this->n() == B.n());
      
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
    }
  
}






template<typename T>
void DenseMatrix<T>::right_multiply (const DenseMatrixBase<T>& M3)
{
  if (this->use_blas)
    this->_multiply_blas(M3, RIGHT_MULTIPLY);
  else
    {
      // (*this) <- M3 * (*this)
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
}




template<typename T>
void DenseMatrix<T>::right_multiply_transpose (const DenseMatrix<T>& B)
{
  if (this->use_blas)
    this->_multiply_blas(B, RIGHT_MULTIPLY_TRANSPOSE);
  else
    {
      //Check to see if we are doing B*(B^T)
      if (this == &B)
	{
	  //libmesh_here();
	  DenseMatrix<T> A(*this);
      
	  // Simple but inefficient way
	  // return this->right_multiply_transpose(A);

	  // More efficient, more code
	  // If B is mxn, the result will be a square matrix of Size m x m.
	  const unsigned int m = B.m();
	  const unsigned int n = B.n();

	  // resize() *this and also zero out all entries.
	  this->resize(m,m);

	  // Compute the lower-triangular part
	  for (unsigned int i=0; i<m; ++i)
	    for (unsigned int j=0; j<=i; ++j)
	      for (unsigned int k=0; k<n; ++k) // inner products are over n
		(*this)(i,j) += A(i,k)*A(j,k);

	  // Copy lower-triangular part into upper-triangular part
	  for (unsigned int i=0; i<m; ++i)
	    for (unsigned int j=i+1; j<m; ++j)
	      (*this)(i,j) = (*this)(j,i);
	}

      else
	{
	  DenseMatrix<T> A(*this);
  
	  this->resize (A.m(), B.m());
      
	  libmesh_assert (A.n() == B.n());
	  libmesh_assert (this->m() == A.m());
	  libmesh_assert (this->n() == B.m());
      
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
    }
}

template<typename T>
void DenseMatrix<T>::vector_mult (DenseVector<T>& dest, 
                                  const DenseVector<T>& arg) const
{
  const unsigned int n_rows = this->m();
  const unsigned int n_cols = this->n();

  // Make sure the sizes are compatible
  libmesh_assert(n_cols == arg.size());
  libmesh_assert(n_rows == dest.size());

  dest.zero();
  DenseMatrix<T> A(*this);

  for(unsigned int i=0; i<n_rows; i++)
    for(unsigned int j=0; j<n_cols; j++)
      dest(i) += A(i,j)*arg(j);
}

template<typename T>
void DenseMatrix<T>::vector_mult_add (DenseVector<T>& dest, 
                                      const T factor,
                                      const DenseVector<T>& arg) const
{
  DenseVector<T> temp(arg.size());
  this->vector_mult(temp, arg);
  dest.add(factor, temp);
}

template<typename T>
void DenseMatrix<T>::get_principal_submatrix (unsigned int sub_m,
                                              unsigned int sub_n,
                                              DenseMatrix<T>& dest) const
{
  libmesh_assert( (sub_m <= this->m()) && (sub_n <= this->n()) );

  dest.resize(sub_m, sub_n);
  for(unsigned int i=0; i<sub_m; i++)
    for(unsigned int j=0; j<sub_n; j++)
      dest(i,j) = (*this)(i,j);
}

template<typename T>
void DenseMatrix<T>::get_principal_submatrix (unsigned int sub_m, DenseMatrix<T>& dest) const
{
  get_principal_submatrix(sub_m, sub_m, dest);
}

template<typename T>
void DenseMatrix<T>::get_transpose (DenseMatrix<T>& dest) const
{
  dest.resize(this->n(), this->m());

  for (unsigned int i=0; i<dest.m(); i++)
    for (unsigned int j=0; j<dest.n(); j++)
      dest(i,j) = (*this)(j,i);
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
	libmesh_error();
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
    n = this->n();

  libmesh_assert (this->m() == n);
  libmesh_assert (this->m() == b.size());
  
  x.resize (n);

  // A convenient reference to *this
  const DenseMatrix<T>& A = *this;

  // Transform the RHS vector
  for (unsigned int i=0; i<(n-1); i++)
    {
      // Get the diagonal entry and take its inverse
      const T diag = A(i,i);
      
      libmesh_assert (diag != libMesh::zero);
  
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

	libmesh_assert (diag != libMesh::zero);

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
  libmesh_assert(this->_decomposition_type == NONE);
  
  // Get the matrix size and make sure it is square
  const unsigned int
    m = this->m();

  libmesh_assert (m == this->n());

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

	  libmesh_assert (diag != libMesh::zero);

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
      libmesh_error();
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
      libmesh_error();
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
	libmesh_error();
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
  libmesh_assert(this->_decomposition_type == NONE);
  
  // Shorthand notation for number of rows and columns.
  const unsigned int
    m = this->m(),
    n = this->n();

  // Just to be really sure...
  libmesh_assert(m==n);

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
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
	      if (A(i,j) <= 0.0)
		{
		  std::cerr << "Error! Can only use Cholesky decomposition "
			    << "with symmetric positive definite matrices."
			    << std::endl;
		  libmesh_error();
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
  libmesh_assert(m==n);

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





#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_multiply_blas(const DenseMatrixBase<T>& other,
				    _BLAS_Multiply_Flag flag)
{
  int result_size = 0;
  
  // For each case, determine the size of the final result make sure
  // that the inner dimensions match
  switch (flag)
    {
    case LEFT_MULTIPLY:
      {
	result_size = other.m() * this->n();
	if (other.n() == this->m())
	  break;
      }
    case RIGHT_MULTIPLY:
      {
	result_size = other.n() * this->m();
	if (other.m() == this->n())
	  break;
      }
    case LEFT_MULTIPLY_TRANSPOSE:
      {
	result_size = other.n() * this->n();
	if (other.m() == this->m())
	  break;
      }
    case RIGHT_MULTIPLY_TRANSPOSE:
      {
	result_size = other.m() * this->m();
	if (other.n() == this->n())
	  break;
      }
    default:
      {
	std::cout << "Unknown flag selected or matrices are ";
	std::cout << "incompatible for multiplication." << std::endl;
	libmesh_error();
      }
    }

  // For this to work, the passed arg. must actually be a DenseMatrix<T>
  const DenseMatrix<T>* const_that = dynamic_cast< const DenseMatrix<T>* >(&other);
  if (!const_that)
    {
      std::cerr << "Unable to cast input matrix to usable type." << std::endl;
      libmesh_error();
    }

  // Also, although 'that' is logically const in this BLAS routine,
  // the PETSc BLAS interface does not specify that any of the inputs are
  // const.  To use it, I must cast away const-ness.
  DenseMatrix<T>* that = const_cast< DenseMatrix<T>* > (const_that);

  // Initialize A, B pointers for LEFT_MULTIPLY* cases
  DenseMatrix<T>
    *A = this,
    *B = that;

  // For RIGHT_MULTIPLY* cases, swap the meaning of A and B.
  // Here is a full table of combinations we can pass to BLASgemm, and what the answer is when finished:
  // pass A B   -> (Fortran) -> A^T B^T -> (C++) -> (A^T B^T)^T -> (identity) -> B A   "lt multiply"
  // pass B A   -> (Fortran) -> B^T A^T -> (C++) -> (B^T A^T)^T -> (identity) -> A B   "rt multiply"
  // pass A B^T -> (Fortran) -> A^T B   -> (C++) -> (A^T B)^T   -> (identity) -> B^T A "lt multiply t"
  // pass B^T A -> (Fortran) -> B A^T   -> (C++) -> (B A^T)^T   -> (identity) -> A B^T "rt multiply t"
  if (flag==RIGHT_MULTIPLY || flag==RIGHT_MULTIPLY_TRANSPOSE)
    std::swap(A,B);

  // transa, transb values to pass to blas
  char
    transa[] = "n",
    transb[] = "n";

  // Integer values to pass to BLAS:
  //
  // M  
  // In Fortran, the number of rows of op(A),
  // In the BLAS documentation, typically known as 'M'.
  //
  // In C/C++, we set:
  // M = n_cols(A) if (transa='n')
  //     n_rows(A) if (transa='t')
  int M = static_cast<int>( A->n() );

  // N
  // In Fortran, the number of cols of op(B), and also the number of cols of C.
  // In the BLAS documentation, typically known as 'N'.
  //	    
  // In C/C++, we set:
  // N = n_rows(B) if (transb='n')
  //     n_cols(B) if (transb='t')
  int N = static_cast<int>( B->m() );

  // K
  // In Fortran, the number of cols of op(A), and also
  // the number of rows of op(B). In the BLAS documentation,
  // typically known as 'K'.
  //
  // In C/C++, we set:
  // K = n_rows(A) if (transa='n')
  //     n_cols(A) if (transa='t')
  int K = static_cast<int>( A->m() );

  // LDA (leading dimension of A). In our cases,
  // LDA is always the number of columns of A.
  int LDA = static_cast<int>( A->n() );

  // LDB (leading dimension of B).  In our cases,
  // LDB is always the number of columns of B.
  int LDB = static_cast<int>( B->n() );

  if (flag == LEFT_MULTIPLY_TRANSPOSE)
    {
      transb[0] = 't';
      N = static_cast<int>( B->n() );
    }

  else if (flag == RIGHT_MULTIPLY_TRANSPOSE)
    {
      transa[0] = 't';
      std::swap(M,K);
    }

  // LDC (leading dimension of C).  LDC is the
  // number of columns in the solution matrix.
  int LDC = M;
  
  // Scalar values to pass to BLAS
  //
  // scalar multiplying the whole product AB 
  T alpha = 1.;
  
  // scalar multiplying C, which is the original matrix.
  T beta  = 0.;

  // Storage for the result
  std::vector<T> result (result_size);

  // Finally ready to call the BLAS
  BLASgemm_(transa, transb, &M, &N, &K, &alpha, &(A->_val[0]), &LDA, &(B->_val[0]), &LDB, &beta, &result[0], &LDC);

  // Update the relevant dimension for this matrix.
  switch (flag)
    {
    case LEFT_MULTIPLY:            { this->_m = other.m(); break; }
    case RIGHT_MULTIPLY:           { this->_n = other.n(); break; }
    case LEFT_MULTIPLY_TRANSPOSE:  { this->_m = other.n(); break; }
    case RIGHT_MULTIPLY_TRANSPOSE: { this->_n = other.m(); break; }
    default:
      {
	std::cout << "Unknown flag selected." << std::endl;
	libmesh_error();
      }
    }
  
  // Swap my data vector with the result
  this->_val.swap(result);
}

#else

template<typename T>
void DenseMatrix<T>::_multiply_blas(const DenseMatrixBase<T>& ,
				    _BLAS_Multiply_Flag )
{
  std::cerr << "No PETSc-provided BLAS available!" << std::endl;
  libmesh_error();
}

#endif


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
template void DenseMatrix<Real>::_cholesky_back_substitute(DenseVector<Real>&, DenseVector<Real>&) const;
template void DenseMatrix<Real>::cholesky_solve(DenseVector<Complex>&, DenseVector<Complex>&);
template void DenseMatrix<Real>::_cholesky_back_substitute(DenseVector<Complex>&, DenseVector<Complex>&) const;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class DenseMatrix<Complex>;
template void DenseMatrix<Complex>::cholesky_solve(DenseVector<Complex>&,DenseVector<Complex>&);
template void DenseMatrix<Complex>::_cholesky_back_substitute(DenseVector<Complex>&, DenseVector<Complex>&) const;
#endif
