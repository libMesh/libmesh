// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"


#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)
#include "libmesh/petsc_macro.h"

EXTERN_C_FOR_PETSC_BEGIN
#include <petscblaslapack.h>
EXTERN_C_FOR_PETSC_END
#endif

#if (LIBMESH_HAVE_SLEPC && LIBMESH_USE_REAL_NUMBERS)
#include <slepcblaslapack.h>
#endif

namespace libMesh
{



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
      libmesh_error_msg("Unknown flag selected or matrices are incompatible for multiplication.");
    }

  // For this to work, the passed arg. must actually be a DenseMatrix<T>
  const DenseMatrix<T>* const_that = cast_ptr< const DenseMatrix<T>* >(&other);

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
      libmesh_error_msg("Unknown flag selected.");
    }

  // Swap my data vector with the result
  this->_val.swap(result);
}

#else

template<typename T>
void DenseMatrix<T>::_multiply_blas(const DenseMatrixBase<T>& ,
                                    _BLAS_Multiply_Flag )
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}

#endif







#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_lu_decompose_lapack ()
{
  // If this function was called, there better not be any
  // previous decomposition of the matrix.
  libmesh_assert_equal_to (this->_decomposition_type, NONE);

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
    libmesh_error_msg("INFO=" << INFO << ", Error during Lapack LU factorization!");

  // Set the flag for LU decomposition
  this->_decomposition_type = LU_BLAS_LAPACK;
}

#else

template<typename T>
void DenseMatrix<T>::_lu_decompose_lapack ()
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}

#endif



template<typename T>
void DenseMatrix<T>::_svd_lapack (DenseVector<T>& sigma)
{
  // The calling sequence for dgetrf is:
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )


  //  JOBU    (input) CHARACTER*1
  //          Specifies options for computing all or part of the matrix U:
  //          = 'A':  all M columns of U are returned in array U:
  //          = 'S':  the first min(m,n) columns of U (the left singular
  //                  vectors) are returned in the array U;
  //          = 'O':  the first min(m,n) columns of U (the left singular
  //                  vectors) are overwritten on the array A;
  //          = 'N':  no columns of U (no left singular vectors) are
  //                  computed.
  char JOBU = 'N';

  //  JOBVT   (input) CHARACTER*1
  //          Specifies options for computing all or part of the matrix
  //          V**T:
  //          = 'A':  all N rows of V**T are returned in the array VT;
  //          = 'S':  the first min(m,n) rows of V**T (the right singular
  //                  vectors) are returned in the array VT;
  //          = 'O':  the first min(m,n) rows of V**T (the right singular
  //                  vectors) are overwritten on the array A;
  //          = 'N':  no rows of V**T (no right singular vectors) are
  //                  computed.
  char JOBVT = 'N';

  std::vector<T> sigma_val;
  std::vector<T> U_val;
  std::vector<T> VT_val;

  _svd_helper(JOBU, JOBVT, sigma_val, U_val, VT_val);

  // Load the singular values into sigma, ignore U_val and VT_val
  const unsigned int n_sigma_vals =
    cast_int<unsigned int>(sigma_val.size());
  sigma.resize(n_sigma_vals);
  for(unsigned int i=0; i<n_sigma_vals; i++)
    sigma(i) = sigma_val[i];

}

template<typename T>
void DenseMatrix<T>::_svd_lapack (DenseVector<T>& sigma, DenseMatrix<T>& U, DenseMatrix<T>& VT)
{
  // The calling sequence for dgetrf is:
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )


  //  JOBU    (input) CHARACTER*1
  //          Specifies options for computing all or part of the matrix U:
  //          = 'A':  all M columns of U are returned in array U:
  //          = 'S':  the first min(m,n) columns of U (the left singular
  //                  vectors) are returned in the array U;
  //          = 'O':  the first min(m,n) columns of U (the left singular
  //                  vectors) are overwritten on the array A;
  //          = 'N':  no columns of U (no left singular vectors) are
  //                  computed.
  char JOBU = 'S';

  //  JOBVT   (input) CHARACTER*1
  //          Specifies options for computing all or part of the matrix
  //          V**T:
  //          = 'A':  all N rows of V**T are returned in the array VT;
  //          = 'S':  the first min(m,n) rows of V**T (the right singular
  //                  vectors) are returned in the array VT;
  //          = 'O':  the first min(m,n) rows of V**T (the right singular
  //                  vectors) are overwritten on the array A;
  //          = 'N':  no rows of V**T (no right singular vectors) are
  //                  computed.
  char JOBVT = 'S';

  std::vector<T> sigma_val;
  std::vector<T> U_val;
  std::vector<T> VT_val;

  _svd_helper(JOBU, JOBVT, sigma_val, U_val, VT_val);

  // Load the singular values into sigma, ignore U_val and VT_val
  const unsigned int n_sigma_vals =
    cast_int<unsigned int>(sigma_val.size());
  sigma.resize(n_sigma_vals);
  for(unsigned int i=0; i<n_sigma_vals; i++)
    sigma(i) = sigma_val[i];

  int M = this->n();
  int N = this->m();
  int min_MN = (M < N) ? M : N;
  U.resize(M,min_MN);
  for(unsigned int i=0; i<U.m(); i++)
    for(unsigned int j=0; j<U.n(); j++)
      {
        unsigned int index = i + j*U.n();  // Column major storage
        U(i,j) = U_val[index];
      }

  VT.resize(min_MN,N);
  for(unsigned int i=0; i<VT.m(); i++)
    for(unsigned int j=0; j<VT.n(); j++)
      {
        unsigned int index = i + j*U.n(); // Column major storage
        VT(i,j) = VT_val[index];
      }

}

#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_svd_helper (char JOBU,
                                  char JOBVT,
                                  std::vector<T>& sigma_val,
                                  std::vector<T>& U_val,
                                  std::vector<T>& VT_val)
{

  //    M       (input) int*
  //            The number of rows of the matrix A.  M >= 0.
  // In C/C++, pass the number of *cols* of A
  int M = this->n();

  //    N       (input) int*
  //            The number of columns of the matrix A.  N >= 0.
  // In C/C++, pass the number of *rows* of A
  int N = this->m();

  int min_MN = (M < N) ? M : N;
  int max_MN = (M > N) ? M : N;

  //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //          On entry, the M-by-N matrix A.
  //          On exit,
  //          if JOBU = 'O',  A is overwritten with the first min(m,n)
  //                          columns of U (the left singular vectors,
  //                          stored columnwise);
  //          if JOBVT = 'O', A is overwritten with the first min(m,n)
  //                          rows of V**T (the right singular vectors,
  //                          stored rowwise);
  //          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
  //                          are destroyed.
  // Here, we pass &(_val[0]).

  //    LDA     (input) int*
  //            The leading dimension of the array A.  LDA >= max(1,M).
  int LDA = M;

  //  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
  //          The singular values of A, sorted so that S(i) >= S(i+1).
  sigma_val.resize( min_MN );

  //  LDU     (input) INTEGER
  //          The leading dimension of the array U.  LDU >= 1; if
  //          JOBU = 'S' or 'A', LDU >= M.
  int LDU = M;

  //  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
  //          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
  //          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
  //          if JOBU = 'S', U contains the first min(m,n) columns of U
  //          (the left singular vectors, stored columnwise);
  //          if JOBU = 'N' or 'O', U is not referenced.
  U_val.resize( LDU*M );

  //  LDVT    (input) INTEGER
  //          The leading dimension of the array VT.  LDVT >= 1; if
  //          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
  int LDVT = N;

  //  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
  //          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
  //          V**T;
  //          if JOBVT = 'S', VT contains the first min(m,n) rows of
  //          V**T (the right singular vectors, stored rowwise);
  //          if JOBVT = 'N' or 'O', VT is not referenced.
  VT_val.resize( LDVT*N );

  //  LWORK   (input) INTEGER
  //          The dimension of the array WORK.
  //          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
  //          For good performance, LWORK should generally be larger.
  //
  //          If LWORK = -1, then a workspace query is assumed; the routine
  //          only calculates the optimal size of the WORK array, returns
  //          this value as the first entry of the WORK array, and no error
  //          message related to LWORK is issued by XERBLA.
  int larger = (3*min_MN+max_MN > 5*min_MN) ? 3*min_MN+max_MN : 5*min_MN;
  int LWORK  = (larger > 1) ? larger : 1;


  //  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
  //          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
  //          superdiagonal elements of an upper bidiagonal matrix B
  //          whose diagonal is in S (not necessarily sorted). B
  //          satisfies A = U * B * VT, so it has the same singular values
  //          as A, and singular vectors related by U and VT.
  std::vector<T> WORK( LWORK );

  //  INFO    (output) INTEGER
  //          = 0:  successful exit.
  //          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //          > 0:  if DBDSQR did not converge, INFO specifies how many
  //                superdiagonals of an intermediate bidiagonal form B
  //                did not converge to zero. See the description of WORK
  //                above for details.
  int INFO = 0;

  // Ready to call the actual factorization routine through PETSc's interface
  LAPACKgesvd_(&JOBU, &JOBVT, &M, &N, &(_val[0]), &LDA, &(sigma_val[0]), &(U_val[0]),
               &LDU, &(VT_val[0]), &LDVT, &(WORK[0]), &LWORK, &INFO);

  // Check return value for errors
  if (INFO != 0)
    libmesh_error_msg("INFO=" << INFO << ", Error during Lapack SVD calculation!");
}


#else

template<typename T>
void DenseMatrix<T>::_svd_helper (char,
                                  char,
                                  std::vector<T>&,
                                  std::vector<T>&,
                                  std::vector<T>&)
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}

#endif




#if (LIBMESH_HAVE_SLEPC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_evd_lapack (DenseVector<T>& lambda_real, DenseVector<T>& lambda_imag)
{
  // The calling sequence for dgeev is:
  // DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, VL, LDVL, VR,
  //         LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )


  //  BALANC  (input) CHARACTER*1
  //          Indicates how the input matrix should be diagonally scaled
  //          and/or permuted to improve the conditioning of its
  //          eigenvalues.
  //          = 'N': Do not diagonally scale or permute;
  char BALANC = 'N';

  //  JOBVL   (input) CHARACTER*1
  //          = 'N': left eigenvectors of A are not computed;
  //          = 'V': left eigenvectors of A are computed.
  char JOBVL = 'N';

  //  JOBVR   (input) CHARACTER*1
  //          = 'N': right eigenvectors of A are not computed;
  //          = 'V': right eigenvectors of A are computed.
  char JOBVR = 'N';

  //  SENSE   (input) CHARACTER*1
  //          Determines which reciprocal condition numbers are computed.
  //          = 'N': None are computed;
  //          = 'E': Computed for eigenvalues only;
  //          = 'V': Computed for right eigenvectors only;
  //          = 'B': Computed for eigenvalues and right eigenvectors.
  char SENSE = 'N';

  //    N       (input) int*
  //            The number of rows/cols of the matrix A.  N >= 0.
  libmesh_assert( this->m() == this->n() );
  int N = this->m();

  //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //          On entry, the N-by-N matrix A.
  //          On exit, A has been overwritten.
  // Here, we pass &(_val[0]).

  //    LDA     (input) int*
  //            The leading dimension of the array A.  LDA >= max(1,N).
  int LDA = N;

  //  WR      (output) DOUBLE PRECISION array, dimension (N)
  //  WI      (output) DOUBLE PRECISION array, dimension (N)
  //          WR and WI contain the real and imaginary parts,
  //          respectively, of the computed eigenvalues.  Complex
  //          conjugate pairs of eigenvalues appear consecutively
  //          with the eigenvalue having the positive imaginary part
  //          first.
  lambda_real.resize(N);
  lambda_imag.resize(N);

  //  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
  //          If JOBVL = 'V', the left eigenvectors u(j) are stored one
  //          after another in the columns of VL, in the same order
  //          as their eigenvalues.
  //          If JOBVL = 'N', VL is not referenced.
  //          If the j-th eigenvalue is real, then u(j) = VL(:,j),
  //          the j-th column of VL.
  //          If the j-th and (j+1)-st eigenvalues form a complex
  //          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
  //          u(j+1) = VL(:,j) - i*VL(:,j+1).
  // Just set to NULL here.

  //  LDVL    (input) INTEGER
  //          The leading dimension of the array VL.  LDVL >= 1; if
  //          JOBVL = 'V', LDVL >= N.
  int LDVL = 1;

  //  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
  //          If JOBVR = 'V', the right eigenvectors v(j) are stored one
  //          after another in the columns of VR, in the same order
  //          as their eigenvalues.
  //          If JOBVR = 'N', VR is not referenced.
  //          If the j-th eigenvalue is real, then v(j) = VR(:,j),
  //          the j-th column of VR.
  //          If the j-th and (j+1)-st eigenvalues form a complex
  //          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
  //          v(j+1) = VR(:,j) - i*VR(:,j+1).
  // Just set to NULL here.

  //  LDVR    (input) INTEGER
  //          The leading dimension of the array VR.  LDVR >= 1; if
  //          JOBVR = 'V', LDVR >= N.
  int LDVR = 1;

  // Outputs (unused)
  int ILO = 0;
  int IHI = 0;
  std::vector<T> SCALE(N);
  T ABNRM;
  std::vector<T> RCONDE(N);
  std::vector<T> RCONDV(N);

  //  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //
  //  LWORK   (input) INTEGER
  //          The dimension of the array WORK.
  int LWORK = 3*N;
  std::vector<T> WORK( LWORK );

  //  IWORK   (workspace) INTEGER array, dimension (2*N-2)
  //          If SENSE = 'N' or 'E', not referenced.
  // Just set to NULL


  //  INFO    (output) INTEGER
  //          = 0:  successful exit
  //          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //          > 0:  if INFO = i, the QR algorithm failed to compute all the
  //                eigenvalues, and no eigenvectors or condition numbers
  //                have been computed; elements 1:ILO-1 and i+1:N of WR
  //                and WI contain eigenvalues which have converged.
  int INFO = 0;

  // Get references to raw data
  std::vector<T>& lambda_real_val = lambda_real.get_values();
  std::vector<T>& lambda_imag_val = lambda_imag.get_values();

  // Ready to call the actual factorization routine through SLEPc's interface
  LAPACKgeevx_( &BALANC, &JOBVL, &JOBVR, &SENSE, &N, &(_val[0]), &LDA, &lambda_real_val[0],
                &lambda_imag_val[0], NULL, &LDVL, NULL, &LDVR, &ILO, &IHI, &SCALE[0], &ABNRM,
                &RCONDE[0], &RCONDV[0], &WORK[0], &LWORK, NULL, &INFO );

  // Check return value for errors
  if (INFO != 0)
    libmesh_error_msg("INFO=" << INFO << ", Error during Lapack eigenvalue calculation!");
}

#else

template<typename T>
void DenseMatrix<T>::_evd_lapack (DenseVector<T>& , DenseVector<T>& )
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
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
    libmesh_error_msg("INFO=" << INFO << ", Error during Lapack LU solve!");

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
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}

#endif





#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_matvec_blas(T alpha, T beta,
                                  DenseVector<T>& dest,
                                  const DenseVector<T>& arg,
                                  bool trans) const
{
  // Ensure that dest and arg sizes are compatible
  if (!trans)
    {
      // dest  ~ A     * arg
      // (mx1)   (mxn) * (nx1)
      if ((dest.size() != this->m()) || (arg.size() != this->n()))
        libmesh_error_msg("Improper input argument sizes!");
    }

  else // trans == true
    {
      // Ensure that dest and arg are proper size
      // dest  ~ A^T   * arg
      // (nx1)   (nxm) * (mx1)
      if ((dest.size() != this->n()) || (arg.size() != this->m()))
        libmesh_error_msg("Improper input argument sizes!");
    }

  // Calling sequence for dgemv:
  //
  // dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

  //   TRANS  - CHARACTER*1, 't' for transpose, 'n' for non-transpose multiply
  // We store everything in row-major order, so pass the transpose flag for
  // non-transposed matvecs and the 'n' flag for transposed matvecs
  char TRANS[] = "t";
  if (trans)
    TRANS[0] = 'n';

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
                                  const DenseVector<T>&,
                                  bool ) const
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
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
                                              const DenseVector<Real>&,
                                              bool ) const;
template void DenseMatrix<Real>::_svd_lapack(DenseVector<Real>&);
template void DenseMatrix<Real>::_svd_lapack(DenseVector<Real>&, DenseMatrix<Real>&, DenseMatrix<Real>&);
template void DenseMatrix<Real>::_svd_helper (char, char, std::vector<Real>&,
                                              std::vector<Real>&,
                                              std::vector<Real>& );
template void DenseMatrix<Real>::_evd_lapack(DenseVector<Real>&, DenseVector<Real>&);

#if !(LIBMESH_USE_REAL_NUMBERS)
template void DenseMatrix<Number>::_multiply_blas(const DenseMatrixBase<Number>&, _BLAS_Multiply_Flag);
template void DenseMatrix<Number>::_lu_decompose_lapack();
template void DenseMatrix<Number>::_lu_back_substitute_lapack(const DenseVector<Number>& ,
                                                              DenseVector<Number>&);
template void DenseMatrix<Number>::_matvec_blas(Number, Number,
                                                DenseVector<Number>& ,
                                                const DenseVector<Number>&,
                                                bool ) const;
template void DenseMatrix<Number>::_svd_lapack(DenseVector<Number>&);
template void DenseMatrix<Number>::_svd_lapack(DenseVector<Number>&, DenseMatrix<Number>&, DenseMatrix<Number>&);
template void DenseMatrix<Number>::_svd_helper (char, char, std::vector<Number>&,
                                                std::vector<Number>&,
                                                std::vector<Number>& );
template void DenseMatrix<Number>::_evd_lapack(DenseVector<Number>&, DenseVector<Number>&);

#endif

} // namespace libMesh
