// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#if (LIBMESH_HAVE_PETSC)
# include "libmesh/petsc_macro.h"
# include <petscblaslapack.h>
#endif

namespace libMesh
{

#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_multiply_blas(const DenseMatrixBase<T> & other,
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
  const DenseMatrix<T> * const_that = cast_ptr<const DenseMatrix<T> *>(&other);

  // Also, although 'that' is logically const in this BLAS routine,
  // the PETSc BLAS interface does not specify that any of the inputs are
  // const.  To use it, I must cast away const-ness.
  DenseMatrix<T> * that = const_cast< DenseMatrix<T> * > (const_that);

  // Initialize A, B pointers for LEFT_MULTIPLY* cases
  DenseMatrix<T> * A = this;
  DenseMatrix<T> * B = that;

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
  PetscBLASInt M = static_cast<PetscBLASInt>( A->n() );

  // N
  // In Fortran, the number of cols of op(B), and also the number of cols of C.
  // In the BLAS documentation, typically known as 'N'.
  //
  // In C/C++, we set:
  // N = n_rows(B) if (transb='n')
  //     n_cols(B) if (transb='t')
  PetscBLASInt N = static_cast<PetscBLASInt>( B->m() );

  // K
  // In Fortran, the number of cols of op(A), and also
  // the number of rows of op(B). In the BLAS documentation,
  // typically known as 'K'.
  //
  // In C/C++, we set:
  // K = n_rows(A) if (transa='n')
  //     n_cols(A) if (transa='t')
  PetscBLASInt K = static_cast<PetscBLASInt>( A->m() );

  // LDA (leading dimension of A). In our cases,
  // LDA is always the number of columns of A.
  PetscBLASInt LDA = static_cast<PetscBLASInt>( A->n() );

  // LDB (leading dimension of B).  In our cases,
  // LDB is always the number of columns of B.
  PetscBLASInt LDB = static_cast<PetscBLASInt>( B->n() );

  if (flag == LEFT_MULTIPLY_TRANSPOSE)
    {
      transb[0] = 't';
      N = static_cast<PetscBLASInt>( B->n() );
    }

  else if (flag == RIGHT_MULTIPLY_TRANSPOSE)
    {
      transa[0] = 't';
      std::swap(M,K);
    }

  // LDC (leading dimension of C).  LDC is the
  // number of columns in the solution matrix.
  PetscBLASInt LDC = M;

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
void DenseMatrix<T>::_multiply_blas(const DenseMatrixBase<T> &,
                                    _BLAS_Multiply_Flag)
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

  // M (input)
  //   The number of rows of the matrix A.  M >= 0.
  // In C/C++, pass the number of *cols* of A
  PetscBLASInt M = this->n();

  // N (input)
  //   The number of columns of the matrix A.  N >= 0.
  // In C/C++, pass the number of *rows* of A
  PetscBLASInt N = this->m();

  // A (input/output) double precision array, dimension (LDA,N)
  //   On entry, the M-by-N matrix to be factored.
  //   On exit, the factors L and U from the factorization
  //   A = P*L*U; the unit diagonal elements of L are not stored.
  // Here, we pass &(_val[0]).

  // LDA (input)
  //     The leading dimension of the array A.  LDA >= max(1,M).
  PetscBLASInt LDA = M;

  // ipiv (output) integer array, dimension (min(m,n))
  //      The pivot indices; for 1 <= i <= min(m,n), row i of the
  //      matrix was interchanged with row IPIV(i).
  // Here, we pass &(_pivots[0]), a private class member used to store pivots
  this->_pivots.resize( std::min(M,N) );

  // info (output)
  //      = 0:  successful exit
  //      < 0:  if INFO = -i, the i-th argument had an illegal value
  //      > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
  //            has been completed, but the factor U is exactly
  //            singular, and division by zero will occur if it is used
  //            to solve a system of equations.
  PetscBLASInt INFO = 0;

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
void DenseMatrix<T>::_svd_lapack (DenseVector<Real> & sigma)
{
  // The calling sequence for dgetrf is:
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

  // JOBU (input)
  //      Specifies options for computing all or part of the matrix U:
  //      = 'A':  all M columns of U are returned in array U:
  //      = 'S':  the first min(m,n) columns of U (the left singular
  //              vectors) are returned in the array U;
  //      = 'O':  the first min(m,n) columns of U (the left singular
  //              vectors) are overwritten on the array A;
  //      = 'N':  no columns of U (no left singular vectors) are
  //              computed.
  char JOBU = 'N';

  // JOBVT (input)
  //       Specifies options for computing all or part of the matrix
  //       V**T:
  //       = 'A':  all N rows of V**T are returned in the array VT;
  //       = 'S':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are returned in the array VT;
  //       = 'O':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are overwritten on the array A;
  //       = 'N':  no rows of V**T (no right singular vectors) are
  //               computed.
  char JOBVT = 'N';

  std::vector<Real> sigma_val;
  std::vector<Number> U_val;
  std::vector<Number> VT_val;

  _svd_helper(JOBU, JOBVT, sigma_val, U_val, VT_val);

  // Copy the singular values into sigma, ignore U_val and VT_val
  sigma.resize(cast_int<unsigned int>(sigma_val.size()));
  for (unsigned int i=0; i<sigma.size(); i++)
    sigma(i) = sigma_val[i];

}

template<typename T>
void DenseMatrix<T>::_svd_lapack (DenseVector<Real> & sigma,
                                  DenseMatrix<Number> & U,
                                  DenseMatrix<Number> & VT)
{
  // The calling sequence for dgetrf is:
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

  // JOBU (input)
  //      Specifies options for computing all or part of the matrix U:
  //      = 'A':  all M columns of U are returned in array U:
  //      = 'S':  the first min(m,n) columns of U (the left singular
  //              vectors) are returned in the array U;
  //      = 'O':  the first min(m,n) columns of U (the left singular
  //              vectors) are overwritten on the array A;
  //      = 'N':  no columns of U (no left singular vectors) are
  //              computed.
  char JOBU = 'S';

  // JOBVT (input)
  //       Specifies options for computing all or part of the matrix
  //       V**T:
  //       = 'A':  all N rows of V**T are returned in the array VT;
  //       = 'S':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are returned in the array VT;
  //       = 'O':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are overwritten on the array A;
  //       = 'N':  no rows of V**T (no right singular vectors) are
  //               computed.
  char JOBVT = 'S';

  // Note: Lapack is going to compute the singular values of A^T.  If
  // A=U * S * V^T, then A^T = V * S * U^T, which means that the
  // values returned in the "U_val" array actually correspond to the
  // entries of the V matrix, and the values returned in the VT_val
  // array actually correspond to the entries of U^T.  Therefore, we
  // pass VT in the place of U and U in the place of VT below!
  std::vector<Real> sigma_val;
  int M = this->n();
  int N = this->m();
  int min_MN = (M < N) ? M : N;

  // Size user-provided storage appropriately. Inside svd_helper:
  // U_val is sized to (M x min_MN)
  // VT_val is sized to (min_MN x N)
  // So, we set up U to have the shape of "VT_val^T", and VT to
  // have the shape of "U_val^T".
  //
  // Finally, since the results are stored in column-major order by
  // Lapack, but we actually want the transpose of what Lapack
  // returns, this means (conveniently) that we don't even have to
  // copy anything after the call to _svd_helper, it should already be
  // in the correct order!
  U.resize(N, min_MN);
  VT.resize(min_MN, M);

  _svd_helper(JOBU, JOBVT, sigma_val, VT.get_values(), U.get_values());

  // Copy the singular values into sigma.
  sigma.resize(cast_int<unsigned int>(sigma_val.size()));
  for (unsigned int i=0; i<sigma.size(); i++)
    sigma(i) = sigma_val[i];
}

#if (LIBMESH_HAVE_PETSC)

template<typename T>
void DenseMatrix<T>::_svd_helper (char JOBU,
                                  char JOBVT,
                                  std::vector<Real> & sigma_val,
                                  std::vector<Number> & U_val,
                                  std::vector<Number> & VT_val)
{

  // M (input)
  //   The number of rows of the matrix A.  M >= 0.
  // In C/C++, pass the number of *cols* of A
  PetscBLASInt M = this->n();

  // N (input)
  //   The number of columns of the matrix A.  N >= 0.
  // In C/C++, pass the number of *rows* of A
  PetscBLASInt N = this->m();

  PetscBLASInt min_MN = (M < N) ? M : N;
  PetscBLASInt max_MN = (M > N) ? M : N;

  // A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //   On entry, the M-by-N matrix A.
  //   On exit,
  //   if JOBU  = 'O', A is overwritten with the first min(m,n)
  //                   columns of U (the left singular vectors,
  //                   stored columnwise);
  //   if JOBVT = 'O', A is overwritten with the first min(m,n)
  //                   rows of V**T (the right singular vectors,
  //                   stored rowwise);
  //   if JOBU != 'O' and JOBVT != 'O', the contents of A are destroyed.
  // Here, we pass &(_val[0]).

  // LDA (input)
  //     The leading dimension of the array A.  LDA >= max(1,M).
  PetscBLASInt LDA = M;

  // S (output) DOUBLE PRECISION array, dimension (min(M,N))
  //   The singular values of A, sorted so that S(i) >= S(i+1).
  sigma_val.resize( min_MN );

  // LDU (input)
  //     The leading dimension of the array U.  LDU >= 1; if
  //     JOBU = 'S' or 'A', LDU >= M.
  PetscBLASInt LDU = M;

  // U (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
  //   (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
  //   If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
  //   if JOBU = 'S', U contains the first min(m,n) columns of U
  //   (the left singular vectors, stored columnwise);
  //   if JOBU = 'N' or 'O', U is not referenced.
  if (JOBU == 'S')
    U_val.resize( LDU*min_MN );
  else
    U_val.resize( LDU*M );

  // LDVT (input)
  //      The leading dimension of the array VT.  LDVT >= 1; if
  //      JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
  PetscBLASInt LDVT = N;
  if (JOBVT == 'S')
    LDVT = min_MN;

  // VT (output) DOUBLE PRECISION array, dimension (LDVT,N)
  //    If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
  //    V**T;
  //    if JOBVT = 'S', VT contains the first min(m,n) rows of
  //    V**T (the right singular vectors, stored rowwise);
  //    if JOBVT = 'N' or 'O', VT is not referenced.
  VT_val.resize( LDVT*N );

  // LWORK (input)
  //       The dimension of the array WORK.
  //       LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
  //       For good performance, LWORK should generally be larger.
  //
  //       If LWORK = -1, then a workspace query is assumed; the routine
  //       only calculates the optimal size of the WORK array, returns
  //       this value as the first entry of the WORK array, and no error
  //       message related to LWORK is issued by XERBLA.
  PetscBLASInt larger = (3*min_MN+max_MN > 5*min_MN) ? 3*min_MN+max_MN : 5*min_MN;
  PetscBLASInt LWORK  = (larger > 1) ? larger : 1;


  // WORK (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //      On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
  //      if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
  //      superdiagonal elements of an upper bidiagonal matrix B
  //      whose diagonal is in S (not necessarily sorted). B
  //      satisfies A = U * B * VT, so it has the same singular values
  //      as A, and singular vectors related by U and VT.
  std::vector<Number> WORK( LWORK );

  // INFO (output)
  //      = 0:  successful exit.
  //      < 0:  if INFO = -i, the i-th argument had an illegal value.
  //      > 0:  if DBDSQR did not converge, INFO specifies how many
  //            superdiagonals of an intermediate bidiagonal form B
  //            did not converge to zero. See the description of WORK
  //            above for details.
  PetscBLASInt INFO = 0;

  // Ready to call the actual factorization routine through PETSc's interface.
#ifdef LIBMESH_USE_REAL_NUMBERS
  // Note that the call to LAPACKgesvd_ may modify _val
  LAPACKgesvd_(&JOBU, &JOBVT, &M, &N, &(_val[0]), &LDA, &(sigma_val[0]), &(U_val[0]),
               &LDU, &(VT_val[0]), &LDVT, &(WORK[0]), &LWORK, &INFO);
#else
  // When we have LIBMESH_USE_COMPLEX_NUMBERS then we must pass an array of Complex
  // numbers to LAPACKgesvd_, but _val may contain Reals so we copy to Number below to
  // handle both the real-valued and complex-valued cases.
  std::vector<Number> val_copy(_val.size());
  for (std::size_t i=0; i<_val.size(); i++)
      val_copy[i] = _val[i];

  std::vector<Real> RWORK(5 * min_MN);
  LAPACKgesvd_(&JOBU, &JOBVT, &M, &N, &(val_copy[0]), &LDA, &(sigma_val[0]), &(U_val[0]),
               &LDU, &(VT_val[0]), &LDVT, &(WORK[0]), &LWORK, &(RWORK[0]), &INFO);
#endif

  // Check return value for errors
  if (INFO != 0)
    libmesh_error_msg("INFO=" << INFO << ", Error during Lapack SVD calculation!");
}


#else

template<typename T>
void DenseMatrix<T>::_svd_helper (char,
                                  char,
                                  std::vector<Real> &,
                                  std::vector<Number> &,
                                  std::vector<Number> &)
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}

#endif



#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)
#if !PETSC_VERSION_LESS_THAN(3,1,0)

template<typename T>
void DenseMatrix<T>::_svd_solve_lapack(const DenseVector<T> & rhs,
                                       DenseVector<T> & x,
                                       Real rcond) const
{
  // Since BLAS is expecting column-major storage, we first need to
  // make a transposed copy of *this, then pass it to the gelss
  // routine instead of the original.  This extra copy is kind of a
  // bummer, it might be better if we could use the full SVD to
  // compute the least-squares solution instead...  Note that it isn't
  // completely terrible either, since A_trans gets overwritten by
  // Lapack, and we usually would end up making a copy of A outside
  // the function call anyway.
  DenseMatrix<T> A_trans;
  this->get_transpose(A_trans);

  // M
  // The number of rows of the input matrix. M >= 0.
  // This is actually the number of *columns* of A_trans.
  PetscBLASInt M = A_trans.n();

  // N
  // The number of columns of the matrix A. N >= 0.
  // This is actually the number of *rows* of A_trans.
  PetscBLASInt N = A_trans.m();

  // We'll use the min and max of (M,N) several times below.
  PetscBLASInt max_MN = std::max(M,N);
  PetscBLASInt min_MN = std::min(M,N);

  // NRHS
  // The number of right hand sides, i.e., the number of columns
  // of the matrices B and X. NRHS >= 0.
  // This could later be generalized to solve for multiple right-hand
  // sides...
  PetscBLASInt NRHS = 1;

  // A is double precision array, dimension (LDA,N)
  // On entry, the M-by-N matrix A.
  // On exit, the first min(m,n) rows of A are overwritten with
  // its right singular vectors, stored rowwise.
  //
  // The data vector that will be passed to Lapack.
  std::vector<T> & A_trans_vals = A_trans.get_values();

  // LDA
  // The leading dimension of the array A.  LDA >= max(1,M).
  PetscBLASInt LDA = M;

  // B is double precision array, dimension (LDB,NRHS)
  // On entry, the M-by-NRHS right hand side matrix B.
  // On exit, B is overwritten by the N-by-NRHS solution
  // matrix X.  If m >= n and RANK = n, the residual
  // sum-of-squares for the solution in the i-th column is given
  // by the sum of squares of elements n+1:m in that column.
  //
  // Since we don't want the user's rhs vector to be overwritten by
  // the solution, we copy the rhs values into the solution vector "x"
  // now.  x needs to be long enough to hold both the (Nx1) solution
  // vector or the (Mx1) rhs, so size it to the max of those.
  x.resize(max_MN);
  for (unsigned i=0; i<rhs.size(); ++i)
    x(i) = rhs(i);

  // Make the syntax below simpler by grabbing a reference to this array.
  std::vector<T> & B = x.get_values();

  // LDB
  // The leading dimension of the array B. LDB >= max(1,max(M,N)).
  PetscBLASInt LDB = x.size();

  // S is double precision array, dimension (min(M,N))
  // The singular values of A in decreasing order.
  // The condition number of A in the 2-norm = S(1)/S(min(m,n)).
  std::vector<T> S(min_MN);

  // RCOND
  // Used to determine the effective rank of A.  Singular values
  // S(i) <= RCOND*S(1) are treated as zero.  If RCOND < 0, machine
  // precision is used instead.
  Real RCOND = rcond;

  // RANK
  // The effective rank of A, i.e., the number of singular values
  // which are greater than RCOND*S(1).
  PetscBLASInt RANK = 0;

  // LWORK
  // The dimension of the array WORK. LWORK >= 1, and also:
  // LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
  // For good performance, LWORK should generally be larger.
  //
  // If LWORK = -1, then a workspace query is assumed; the routine
  // only calculates the optimal size of the WORK array, returns
  // this value as the first entry of the WORK array, and no error
  // message related to LWORK is issued by XERBLA.
  //
  // The factor of 1.5 is arbitrary and is used to satisfy the "should
  // generally be larger" clause.
  PetscBLASInt LWORK = 1.5 * (3*min_MN + std::max(2*min_MN, std::max(max_MN, NRHS)));

  // WORK is double precision array, dimension (MAX(1,LWORK))
  // On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  std::vector<T> WORK(LWORK);

  // INFO
  // = 0:  successful exit
  // < 0:  if INFO = -i, the i-th argument had an illegal value.
  // > 0:  the algorithm for computing the SVD failed to converge;
  //       if INFO = i, i off-diagonal elements of an intermediate
  //       bidiagonal form did not converge to zero.
  PetscBLASInt INFO = 0;

  // LAPACKgelss_(const PetscBLASInt *, // M
  //              const PetscBLASInt *, // N
  //              const PetscBLASInt *, // NRHS
  //              PetscScalar *,        // A
  //              const PetscBLASInt *, // LDA
  //              PetscScalar *,        // B
  //              const PetscBLASInt *, // LDB
  //              PetscReal *,          // S(out) = singular values of A in increasing order
  //              const PetscReal *,    // RCOND = tolerance for singular values
  //              PetscBLASInt *,       // RANK(out) = number of "non-zero" singular values
  //              PetscScalar *,        // WORK
  //              const PetscBLASInt *, // LWORK
  //              PetscBLASInt *);      // INFO
  LAPACKgelss_(&M, &N, &NRHS, &A_trans_vals[0], &LDA, &B[0], &LDB, &S[0], &RCOND, &RANK, &WORK[0], &LWORK, &INFO);

  // Check for errors in the Lapack call
  if (INFO < 0)
    libmesh_error_msg("Error, argument " << -INFO << " to LAPACKgelss_ had an illegal value.");
  if (INFO > 0)
    libmesh_error_msg("The algorithm for computing the SVD failed to converge!");

  // Debugging: print singular values and information about condition number:
  // libMesh::err << "RCOND=" << RCOND << std::endl;
  // libMesh::err << "Singular values: " << std::endl;
  // for (std::size_t i=0; i<S.size(); ++i)
  //   libMesh::err << S[i] << std::endl;
  // libMesh::err << "The condition number of A is approximately: " << S[0]/S.back() << std::endl;

  // Lapack has already written the solution into B, but it will be
  // the wrong size for non-square problems, so we need to resize it
  // correctly.  The size of the solution vector should be the number
  // of columns of the original A matrix.  Unfortunately, resizing a
  // DenseVector currently also zeros it out (unlike a std::vector) so
  // we'll resize the underlying storage directly (the size is not
  // stored independently elsewhere).
  x.get_values().resize(this->n());
}

#else

template<typename T>
void DenseMatrix<T>::_svd_solve_lapack(const DenseVector<T> & /*rhs*/,
                                       DenseVector<T> & /*x*/,
                                       Real /*rcond*/) const
{
  libmesh_error_msg("svd_solve() requires PETSc >= 3.1!");
}

#endif // !PETSC_VERSION_LESS_THAN(3,1,0)

#else
template<typename T>
void DenseMatrix<T>::_svd_solve_lapack(const DenseVector<T>& /*rhs*/,
                                       DenseVector<T> & /*x*/,
                                       Real /*rcond*/) const
{
  libmesh_not_implemented();
}
#endif // (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)



#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_evd_lapack (DenseVector<T> & lambda_real,
                                  DenseVector<T> & lambda_imag,
                                  DenseMatrix<T> * VL,
                                  DenseMatrix<T> * VR)
{
  // This algorithm only works for square matrices, so verify this up front.
  if (this->m() != this->n())
    libmesh_error_msg("Can only compute eigen-decompositions for square matrices.");

  // If the user requests left or right eigenvectors, we have to make
  // sure and pass the transpose of this matrix to Lapack, otherwise
  // it will compute the inverse transpose of what we are
  // after... since we know the matrix is square, we can just swap
  // entries in place.  If the user does not request eigenvectors, we
  // can skip this extra step, since the eigenvalues for A and A^T are
  // the same.
  if (VL || VR)
    {
      for (unsigned int i=0; i<this->_m; ++i)
        for (unsigned int j=0; j<i; ++j)
          std::swap((*this)(i,j), (*this)(j,i));
    }

  // The calling sequence for dgeev is:
  // DGEEV (character  JOBVL,
  //        character  JOBVR,
  //        integer N,
  //        double precision, dimension( lda, * )  A,
  //        integer  LDA,
  //        double precision, dimension( * )  WR,
  //        double precision, dimension( * )  WI,
  //        double precision, dimension( ldvl, * )  VL,
  //        integer  LDVL,
  //        double precision, dimension( ldvr, * )  VR,
  //        integer  LDVR,
  //        double precision, dimension( * )  WORK,
  //        integer  LWORK,
  //        integer  INFO)

  // JOBVL (input)
  //       = 'N': left eigenvectors of A are not computed;
  //       = 'V': left eigenvectors of A are computed.
  char JOBVL = VL ? 'V' : 'N';

  // JOBVR (input)
  //       = 'N': right eigenvectors of A are not computed;
  //       = 'V': right eigenvectors of A are computed.
  char JOBVR = VR ? 'V' : 'N';

  // N (input)
  //   The number of rows/cols of the matrix A.  N >= 0.
  PetscBLASInt N = this->m();

  // A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //   On entry, the N-by-N matrix A.
  //   On exit, A has been overwritten.
  // Here, we pass &(_val[0]).

  // LDA (input)
  //     The leading dimension of the array A.  LDA >= max(1,N).
  PetscBLASInt LDA = N;

  // WR (output) double precision array, dimension (N)
  // WI (output) double precision array, dimension (N)
  //    WR and WI contain the real and imaginary parts,
  //    respectively, of the computed eigenvalues.  Complex
  //    conjugate pairs of eigenvalues appear consecutively
  //    with the eigenvalue having the positive imaginary part
  //    first.
  lambda_real.resize(N);
  lambda_imag.resize(N);

  // VL (output) double precision array, dimension (LDVL,N)
  //    If JOBVL = 'V', the left eigenvectors u(j) are stored one
  //    after another in the columns of VL, in the same order
  //    as their eigenvalues.
  //    If JOBVL = 'N', VL is not referenced.
  //    If the j-th eigenvalue is real, then u(j) = VL(:,j),
  //    the j-th column of VL.
  //    If the j-th and (j+1)-st eigenvalues form a complex
  //    conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
  //    u(j+1) = VL(:,j) - i*VL(:,j+1).
  // Will be set below if needed.

  // LDVL (input)
  //      The leading dimension of the array VL.  LDVL >= 1; if
  //      JOBVL = 'V', LDVL >= N.
  PetscBLASInt LDVL = VL ? N : 1;

  // VR (output) DOUBLE PRECISION array, dimension (LDVR,N)
  //    If JOBVR = 'V', the right eigenvectors v(j) are stored one
  //    after another in the columns of VR, in the same order
  //    as their eigenvalues.
  //    If JOBVR = 'N', VR is not referenced.
  //    If the j-th eigenvalue is real, then v(j) = VR(:,j),
  //    the j-th column of VR.
  //    If the j-th and (j+1)-st eigenvalues form a complex
  //    conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
  //    v(j+1) = VR(:,j) - i*VR(:,j+1).
  // Will be set below if needed.

  // LDVR (input)
  //      The leading dimension of the array VR.  LDVR >= 1; if
  //      JOBVR = 'V', LDVR >= N.
  PetscBLASInt LDVR = VR ? N : 1;

  // WORK (workspace/output) double precision array, dimension (MAX(1,LWORK))
  //      On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //
  // LWORK (input)
  //       The dimension of the array WORK.  LWORK >= max(1,3*N), and
  //       if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
  //       performance, LWORK must generally be larger.
  //
  //       If LWORK = -1, then a workspace query is assumed; the routine
  //       only calculates the optimal size of the WORK array, returns
  //       this value as the first entry of the WORK array, and no error
  //       message related to LWORK is issued by XERBLA.
  PetscBLASInt LWORK = (VR || VL) ? 4*N : 3*N;
  std::vector<T> WORK(LWORK);

  // INFO (output)
  //      = 0:  successful exit
  //      < 0:  if INFO = -i, the i-th argument had an illegal value.
  //      > 0:  if INFO = i, the QR algorithm failed to compute all the
  //            eigenvalues, and no eigenvectors or condition numbers
  //            have been computed; elements 1:ILO-1 and i+1:N of WR
  //            and WI contain eigenvalues which have converged.
  PetscBLASInt INFO = 0;

  // Get references to raw data
  std::vector<T> & lambda_real_val = lambda_real.get_values();
  std::vector<T> & lambda_imag_val = lambda_imag.get_values();

  // Set up eigenvector storage if necessary.
  T * VR_ptr = libmesh_nullptr;
  if (VR)
    {
      VR->resize(N, N);
      VR_ptr = &(VR->get_values()[0]);
    }

  T * VL_ptr = libmesh_nullptr;
  if (VL)
    {
      VL->resize(N, N);
      VL_ptr = &(VL->get_values()[0]);
    }

  // Ready to call the Lapack routine through PETSc's interface
  LAPACKgeev_(&JOBVL,
              &JOBVR,
              &N,
              &(_val[0]),
              &LDA,
              &lambda_real_val[0],
              &lambda_imag_val[0],
              VL_ptr,
              &LDVL,
              VR_ptr,
              &LDVR,
              &WORK[0],
              &LWORK,
              &INFO);

  // Check return value for errors
  if (INFO != 0)
    libmesh_error_msg("INFO=" << INFO << ", Error during Lapack eigenvalue calculation!");

  // If the user requested either right or left eigenvectors, LAPACK
  // has now computed the transpose of the desired matrix, i.e. V^T
  // instead of V.  We could leave this up to user code to handle, but
  // rather than risking getting very unexpected results, we'll just
  // transpose it in place before handing it back.
  if (VR)
    {
      for (unsigned int i=0; i<static_cast<unsigned int>(N); ++i)
        for (unsigned int j=0; j<i; ++j)
          std::swap((*VR)(i,j), (*VR)(j,i));
    }

 if (VL)
   {
     for (unsigned int i=0; i<static_cast<unsigned int>(N); ++i)
       for (unsigned int j=0; j<i; ++j)
         std::swap((*VL)(i,j), (*VL)(j,i));
   }
}

#else

template<typename T>
void DenseMatrix<T>::_evd_lapack (DenseVector<T> &,
                                  DenseVector<T> &,
                                  DenseMatrix<T> *,
                                  DenseMatrix<T> *)
{
  libmesh_error_msg("_evd_lapack is currently only available when LIBMESH_USE_REAL_NUMBERS is defined!");
}

#endif





#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_lu_back_substitute_lapack (const DenseVector<T> & b,
                                                 DenseVector<T> & x)
{
  // The calling sequence for getrs is:
  // dgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

  // trans (input)
  //       'n' for no tranpose, 't' for transpose
  char TRANS[] = "t";

  // N (input)
  //   The order of the matrix A.  N >= 0.
  PetscBLASInt N = this->m();


  // NRHS (input)
  //      The number of right hand sides, i.e., the number of columns
  //      of the matrix B.  NRHS >= 0.
  PetscBLASInt NRHS = 1;

  // A (input) double precision array, dimension (LDA,N)
  //   The factors L and U from the factorization A = P*L*U
  //   as computed by dgetrf.
  // Here, we pass &(_val[0])

  // LDA (input)
  //     The leading dimension of the array A.  LDA >= max(1,N).
  PetscBLASInt LDA = N;

  // ipiv (input) int array, dimension (N)
  //      The pivot indices from DGETRF; for 1<=i<=N, row i of the
  //      matrix was interchanged with row IPIV(i).
  // Here, we pass &(_pivots[0]) which was computed in _lu_decompose_lapack

  // B (input/output) double precision array, dimension (LDB,NRHS)
  //   On entry, the right hand side matrix B.
  //   On exit, the solution matrix X.
  // Here, we pass a copy of the rhs vector's data array in x, so that the
  // passed right-hand side b is unmodified.  I don't see a way around this
  // copy if we want to maintain an unmodified rhs in LibMesh.
  x = b;
  std::vector<T> & x_vec = x.get_values();

  // We can avoid the copy if we don't care about overwriting the RHS: just
  // pass b to the Lapack routine and then swap with x before exiting
  // std::vector<T> & x_vec = b.get_values();

  // LDB (input)
  //     The leading dimension of the array B.  LDB >= max(1,N).
  PetscBLASInt LDB = N;

  // INFO (output)
  //      = 0:  successful exit
  //      < 0:  if INFO = -i, the i-th argument had an illegal value
  PetscBLASInt INFO = 0;

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
void DenseMatrix<T>::_lu_back_substitute_lapack (const DenseVector<T> &,
                                                 DenseVector<T> &)
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}

#endif





#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)

template<typename T>
void DenseMatrix<T>::_matvec_blas(T alpha,
                                  T beta,
                                  DenseVector<T> & dest,
                                  const DenseVector<T> & arg,
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

  // TRANS (input)
  //       't' for transpose, 'n' for non-transpose multiply
  // We store everything in row-major order, so pass the transpose flag for
  // non-transposed matvecs and the 'n' flag for transposed matvecs
  char TRANS[] = "t";
  if (trans)
    TRANS[0] = 'n';

  // M (input)
  //   On entry, M specifies the number of rows of the matrix A.
  // In C/C++, pass the number of *cols* of A
  PetscBLASInt M = this->n();

  // N (input)
  //   On entry, N specifies the number of columns of the matrix A.
  // In C/C++, pass the number of *rows* of A
  PetscBLASInt N = this->m();

  // ALPHA (input)
  // The scalar constant passed to this function

  // A (input) double precision array of DIMENSION ( LDA, n ).
  //   Before entry, the leading m by n part of the array A must
  //   contain the matrix of coefficients.
  // The matrix, *this.  Note that _matvec_blas is called from
  // a const function, vector_mult(), and so we have made this function const
  // as well.  Since BLAS knows nothing about const, we have to cast it away
  // now.
  DenseMatrix<T> & a_ref = const_cast< DenseMatrix<T> &> ( *this );
  std::vector<T> & a = a_ref.get_values();

  // LDA (input)
  //     On entry, LDA specifies the first dimension of A as declared
  //     in the calling (sub) program. LDA must be at least
  //     max( 1, m ).
  PetscBLASInt LDA = M;

  // X (input) double precision array of DIMENSION at least
  //   ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  //   and at least
  //   ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  //   Before entry, the incremented array X must contain the
  //   vector x.
  // Here, we must cast away the const-ness of "arg" since BLAS knows
  // nothing about const
  DenseVector<T> & x_ref = const_cast< DenseVector<T> &> ( arg );
  std::vector<T> & x = x_ref.get_values();

  // INCX (input)
  //      On entry, INCX specifies the increment for the elements of
  //      X. INCX must not be zero.
  PetscBLASInt INCX = 1;

  // BETA (input)
  //      On entry, BETA specifies the scalar beta. When BETA is
  //      supplied as zero then Y need not be set on input.
  // The second scalar constant passed to this function

  // Y (input) double precision array of DIMENSION at least
  //   ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  //   and at least
  //   ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  //   Before entry with BETA non-zero, the incremented array Y
  //   must contain the vector y. On exit, Y is overwritten by the
  //   updated vector y.
  // The input vector "dest"
  std::vector<T> & y = dest.get_values();

  // INCY (input)
  //      On entry, INCY specifies the increment for the elements of
  //      Y. INCY must not be zero.
  PetscBLASInt INCY = 1;

  // Finally, ready to call the BLAS function
  BLASgemv_(TRANS, &M, &N, &alpha, &(a[0]), &LDA, &(x[0]), &INCX, &beta, &(y[0]), &INCY);
}


#else


template<typename T>
void DenseMatrix<T>::_matvec_blas(T,
                                  T,
                                  DenseVector<T> &,
                                  const DenseVector<T> &,
                                  bool) const
{
  libmesh_error_msg("No PETSc-provided BLAS/LAPACK available!");
}


#endif


//--------------------------------------------------------------
// Explicit instantiations
template void DenseMatrix<Real>::_multiply_blas(const DenseMatrixBase<Real> &, _BLAS_Multiply_Flag);
template void DenseMatrix<Real>::_lu_decompose_lapack();
template void DenseMatrix<Real>::_lu_back_substitute_lapack(const DenseVector<Real> &,
                                                            DenseVector<Real> &);
template void DenseMatrix<Real>::_matvec_blas(Real,
                                              Real,
                                              DenseVector<Real> &,
                                              const DenseVector<Real> &,
                                              bool) const;
template void DenseMatrix<Real>::_svd_lapack(DenseVector<Real> &);
template void DenseMatrix<Real>::_svd_lapack(DenseVector<Real> &,
                                             DenseMatrix<Number> &,
                                             DenseMatrix<Number> &);
template void DenseMatrix<Real>::_svd_helper (char,
                                              char,
                                              std::vector<Real> &,
                                              std::vector<Number> &,
                                              std::vector<Number> &);
template void DenseMatrix<Real>::_svd_solve_lapack (const DenseVector<Real> &, DenseVector<Real> &, Real) const;
template void DenseMatrix<Real>::_evd_lapack(DenseVector<Real> &,
                                             DenseVector<Real> &,
                                             DenseMatrix<Real> *,
                                             DenseMatrix<Real> *);

#if !(LIBMESH_USE_REAL_NUMBERS)
template void DenseMatrix<Number>::_multiply_blas(const DenseMatrixBase<Number> &, _BLAS_Multiply_Flag);
template void DenseMatrix<Number>::_lu_decompose_lapack();
template void DenseMatrix<Number>::_lu_back_substitute_lapack(const DenseVector<Number> &,
                                                              DenseVector<Number> &);
template void DenseMatrix<Number>::_matvec_blas(Number,
                                                Number,
                                                DenseVector<Number> &,
                                                const DenseVector<Number> &,
                                                bool) const;
template void DenseMatrix<Number>::_svd_lapack(DenseVector<Real> &);
template void DenseMatrix<Number>::_svd_lapack(DenseVector<Real> &,
                                               DenseMatrix<Number> &,
                                               DenseMatrix<Number> &);
template void DenseMatrix<Number>::_svd_helper (char,
                                                char,
                                                std::vector<Real> &,
                                                std::vector<Number> &,
                                                std::vector<Number> &);
template void DenseMatrix<Number>::_svd_solve_lapack (const DenseVector<Number> &, DenseVector<Number> &, Real) const;
template void DenseMatrix<Number>::_evd_lapack(DenseVector<Number> &,
                                               DenseVector<Number> &,
                                               DenseMatrix<Number> *,
                                               DenseMatrix<Number> *);

#endif

} // namespace libMesh
