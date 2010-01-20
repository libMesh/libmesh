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


// Local Includes
#include "dense_matrix.h"
#include "dense_vector.h"


#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)
#include "petsc_macro.h"

EXTERN_C_FOR_PETSC_BEGIN
#include <petscblaslapack.h>
EXTERN_C_FOR_PETSC_END
#endif



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
  std::cerr << "No PETSc-provided BLAS/LAPACK available!" << std::endl;
  libmesh_error();
}

#endif







#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_lu_decompose_lapack ()
{
  // If this function was called, there better not be any
  // previous decomposition of the matrix.
  libmesh_assert(this->_decomposition_type == NONE);

  // The calling sequence for dgetrf is:
  // dgetrf(M, N, A, lda, ipiv, info)
  
  //    M       (input) int*
  //            The number of rows of the matrix A.  M >= 0.
  // In C/C++, pass the number of *cols* of A
  int M = this->n();
    
  //    N       (input) int*
  //            The number of columns of the matrix A.  N >= 0.
  // In C/C++, pass the number of *rows* of A
  int N = this->m();
  
  //    A (input/output) double precision array, dimension (LDA,N)
  //      On entry, the M-by-N matrix to be factored.
  //      On exit, the factors L and U from the factorization
  //      A = P*L*U; the unit diagonal elements of L are not stored.
  // Here, we pass &(_val[0]).
  
  //    LDA     (input) int*
  //            The leading dimension of the array A.  LDA >= max(1,M).
  int LDA = M;
  
  //    ipiv    (output) integer array, dimension (min(m,n))
  //            The pivot indices; for 1 <= i <= min(m,n), row i of the
  //            matrix was interchanged with row IPIV(i).
  // Here, we pass &(_pivots[0]), a private class member used to store pivots
  this->_pivots.resize( std::min(M,N) );
  
  //    info    (output) int*
  //            = 0:  successful exit
  //            < 0:  if INFO = -i, the i-th argument had an illegal value
  //            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
  //                  has been completed, but the factor U is exactly
  //                  singular, and division by zero will occur if it is used
  //                  to solve a system of equations.
  int INFO = 0;

  // Ready to call the actual factorization routine through PETSc's interface
  LAPACKgetrf_(&M, &N, &(this->_val[0]), &LDA, &(_pivots[0]), &INFO);

  // Check return value for errors
  if (INFO != 0)
    {
      std::cout << "INFO="
		<< INFO
		<< ", Error during Lapack LU factorization!" << std::endl;
      libmesh_error();
    }
  
  // Set the flag for LU decomposition
  this->_decomposition_type = LU_BLAS_LAPACK;
}

#else

template<typename T>
void DenseMatrix<T>::_lu_decompose_lapack ()
{
  std::cerr << "No PETSc-provided BLAS/LAPACK available!" << std::endl;
  libmesh_error();
}

#endif





#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_lu_back_substitute_lapack (const DenseVector<T>& b,
						 DenseVector<T>& x) 
{
  // The calling sequence for getrs is:
  // dgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

  //    trans   (input) char*
  //            'n' for no tranpose, 't' for transpose
  char TRANS[] = "t";
  
  //    N       (input) int*
  //            The order of the matrix A.  N >= 0.
  int N = this->m();
  
 
  //    NRHS    (input) int*
  //            The number of right hand sides, i.e., the number of columns
  //            of the matrix B.  NRHS >= 0.
  int NRHS = 1;
 
  //    A       (input) DOUBLE PRECISION array, dimension (LDA,N)
  //            The factors L and U from the factorization A = P*L*U
  //            as computed by dgetrf.
  // Here, we pass &(_val[0])
  
  //    LDA     (input) int*
  //            The leading dimension of the array A.  LDA >= max(1,N).
  int LDA = N;
    
  //    ipiv    (input) int array, dimension (N)
  //            The pivot indices from DGETRF; for 1<=i<=N, row i of the
  //            matrix was interchanged with row IPIV(i).
  // Here, we pass &(_pivots[0]) which was computed in _lu_decompose_lapack
 
  //    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
  //            On entry, the right hand side matrix B.
  //            On exit, the solution matrix X.
  // Here, we pass a copy of the rhs vector's data array in x, so that the
  // passed right-hand side b is unmodified.  I don't see a way around this
  // copy if we want to maintain an unmodified rhs in LibMesh.
  x = b;
  std::vector<T>& x_vec = x.get_values();

  // We can avoid the copy if we don't care about overwriting the RHS: just
  // pass b to the Lapack routine and then swap with x before exiting
  // std::vector<T>& x_vec = b.get_values();
 
  //    LDB     (input) int*
  //            The leading dimension of the array B.  LDB >= max(1,N).
  int LDB = N;
 
  //    INFO    (output) int*
  //            = 0:  successful exit
  //            < 0:  if INFO = -i, the i-th argument had an illegal value
  int INFO = 0;
  
  // Finally, ready to call the Lapack getrs function
  LAPACKgetrs_(TRANS, &N, &NRHS, &(_val[0]), &LDA, &(_pivots[0]), &(x_vec[0]), &LDB, &INFO);

  // Check return value for errors
  if (INFO != 0)
    {
      std::cout << "INFO="
		<< INFO
		<< ", Error during Lapack LU solve!" << std::endl;
      libmesh_error();
    }

  // Don't do this if you already made a copy of b above
  // Swap b and x.  The solution will then be in x, and whatever was originally
  // in x, maybe garbage, maybe nothing, will be in b.
  // FIXME: Rewrite the LU and Cholesky solves to just take one input, and overwrite
  // the input.  This *should* make user code simpler, as they don't have to create
  // an extra vector just to pass it in to the solve function!
  // b.swap(x);
}

#else

template<typename T>
void DenseMatrix<T>::_lu_back_substitute_lapack (const DenseVector<T>& ,
						 DenseVector<T>& )
{
  std::cerr << "No PETSc-provided BLAS/LAPACK available!" << std::endl;
  libmesh_error();
}

#endif





#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_matvec_blas(T alpha, T beta,
				  DenseVector<T>& dest,
				  const DenseVector<T>& arg) const
{
  // Ensure that dest and arg are proper size
  // dest  ~ A     * arg
  // (mx1)   (mxn) * (nx1)
  if ((dest.size() != this->m()) || (arg.size() != this->n()))
    {
      std::cout << "Improper input argument sizes!" << std::endl;
      libmesh_error();
    }
  
  // Calling sequence for dgemv:
  //
  // dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  
  //   TRANS  - CHARACTER*1, 't' for transpose, 'n' for non-transpose multiply
  char TRANS[] = "t";
  
  //   M      - INTEGER.
  //            On entry, M specifies the number of rows of the matrix A.
  // In C/C++, pass the number of *cols* of A
  int M = this->n();
    
  //   N      - INTEGER.
  //            On entry, N specifies the number of columns of the matrix A.
  // In C/C++, pass the number of *rows* of A
  int N = this->m();
  
  //   ALPHA  - DOUBLE PRECISION.
  // The scalar constant passed to this function

  //   A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  //            Before entry, the leading m by n part of the array A must
  //            contain the matrix of coefficients.
  // The matrix, *this.  Note that _matvec_blas is called from
  // a const function, vector_mult(), and so we have made this function const
  // as well.  Since BLAS knows nothing about const, we have to cast it away
  // now.
  DenseMatrix<T>& a_ref = const_cast< DenseMatrix<T>& > ( *this );
  std::vector<T>& a = a_ref.get_values();

  //   LDA    - INTEGER.
  //            On entry, LDA specifies the first dimension of A as declared
  //            in the calling (sub) program. LDA must be at least
  //            max( 1, m ).
  int LDA = M;
  
  //   X      - DOUBLE PRECISION array of DIMENSION at least
  //            ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  //            and at least
  //            ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  //            Before entry, the incremented array X must contain the
  //            vector x.
  // Here, we must cast away the const-ness of "arg" since BLAS knows
  // nothing about const
  DenseVector<T>& x_ref = const_cast< DenseVector<T>& > ( arg );
  std::vector<T>& x = x_ref.get_values();
  
  //   INCX   - INTEGER.
  //            On entry, INCX specifies the increment for the elements of
  //            X. INCX must not be zero.
  int INCX = 1;
  
  //   BETA   - DOUBLE PRECISION.
  //            On entry, BETA specifies the scalar beta. When BETA is
  //            supplied as zero then Y need not be set on input.
  // The second scalar constant passed to this function

  //   Y      - DOUBLE PRECISION array of DIMENSION at least
  //            ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  //            and at least
  //            ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  //            Before entry with BETA non-zero, the incremented array Y
  //            must contain the vector y. On exit, Y is overwritten by the
  //            updated vector y.
  // The input vector "dest"
  std::vector<T>& y = dest.get_values();

  //   INCY   - INTEGER.
  //            On entry, INCY specifies the increment for the elements of
  //            Y. INCY must not be zero.
  int INCY = 1;
  
  // Finally, ready to call the BLAS function
  BLASgemv_(TRANS, &M, &N, &alpha, &(a[0]), &LDA, &(x[0]), &INCX, &beta, &(y[0]), &INCY);
}


#else


template<typename T>
void DenseMatrix<T>::_matvec_blas(T , T,
				  DenseVector<T>& ,
				  const DenseVector<T>& ) const
{
  std::cerr << "No PETSc-provided BLAS/LAPACK available!" << std::endl;
  libmesh_error();
}


#endif


//--------------------------------------------------------------
// Explicit instantiations
template void DenseMatrix<Real>::_multiply_blas(const DenseMatrixBase<Real>&, _BLAS_Multiply_Flag);
template void DenseMatrix<Real>::_lu_decompose_lapack();
template void DenseMatrix<Real>::_lu_back_substitute_lapack(const DenseVector<Real>& ,
							    DenseVector<Real>&);
template void DenseMatrix<Real>::_matvec_blas(Real, Real,
					      DenseVector<Real>& ,
					      const DenseVector<Real>& ) const;

#if !(LIBMESH_USE_REAL_NUMBERS)
template void DenseMatrix<Number>::_multiply_blas(const DenseMatrixBase<Number>&, _BLAS_Multiply_Flag);
template void DenseMatrix<Number>::_lu_decompose_lapack();
template void DenseMatrix<Number>::_lu_back_substitute_lapack(const DenseVector<Number>& ,
							      DenseVector<Number>&);
template void DenseMatrix<Number>::_matvec_blas(Number, Number,
					        DenseVector<Number>& ,
					        const DenseVector<Number>& ) const;
#endif

