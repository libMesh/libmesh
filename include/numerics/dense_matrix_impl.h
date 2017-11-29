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


#ifndef LIBMESH_DENSE_MATRIX_IMPL_H
#define LIBMESH_DENSE_MATRIX_IMPL_H

// C++ Includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for sqrt

// Local Includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/libmesh.h"

namespace libMesh
{



// ------------------------------------------------------------
// Dense Matrix member functions

template<typename T>
void DenseMatrix<T>::left_multiply (const DenseMatrixBase<T> & M2)
{
  if (this->use_blas_lapack)
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
template<typename T2>
void DenseMatrix<T>::left_multiply (const DenseMatrixBase<T2> & M2)
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



template<typename T>
void DenseMatrix<T>::left_multiply_transpose(const DenseMatrix<T> & A)
{
  if (this->use_blas_lapack)
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
          const unsigned int n_rows = A.m();
          const unsigned int n_cols = A.n();

          // resize() *this and also zero out all entries.
          this->resize(n_cols,n_cols);

          // Compute the lower-triangular part
          for (unsigned int i=0; i<n_cols; ++i)
            for (unsigned int j=0; j<=i; ++j)
              for (unsigned int k=0; k<n_rows; ++k) // inner products are over n_rows
                (*this)(i,j) += B(k,i)*B(k,j);

          // Copy lower-triangular part into upper-triangular part
          for (unsigned int i=0; i<n_cols; ++i)
            for (unsigned int j=i+1; j<n_cols; ++j)
              (*this)(i,j) = (*this)(j,i);
        }

      else
        {
          DenseMatrix<T> B(*this);

          this->resize (A.n(), B.n());

          libmesh_assert_equal_to (A.m(), B.m());
          libmesh_assert_equal_to (this->m(), A.n());
          libmesh_assert_equal_to (this->n(), B.n());

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
template<typename T2>
void DenseMatrix<T>::left_multiply_transpose(const DenseMatrix<T2> & A)
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
      const unsigned int n_rows = A.m();
      const unsigned int n_cols = A.n();

      // resize() *this and also zero out all entries.
      this->resize(n_cols,n_cols);

      // Compute the lower-triangular part
      for (unsigned int i=0; i<n_cols; ++i)
        for (unsigned int j=0; j<=i; ++j)
          for (unsigned int k=0; k<n_rows; ++k) // inner products are over n_rows
            (*this)(i,j) += B(k,i)*B(k,j);

      // Copy lower-triangular part into upper-triangular part
      for (unsigned int i=0; i<n_cols; ++i)
        for (unsigned int j=i+1; j<n_cols; ++j)
          (*this)(i,j) = (*this)(j,i);
    }

  else
    {
      DenseMatrix<T> B(*this);

      this->resize (A.n(), B.n());

      libmesh_assert_equal_to (A.m(), B.m());
      libmesh_assert_equal_to (this->m(), A.n());
      libmesh_assert_equal_to (this->n(), B.n());

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



template<typename T>
void DenseMatrix<T>::right_multiply (const DenseMatrixBase<T> & M3)
{
  if (this->use_blas_lapack)
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
template<typename T2>
void DenseMatrix<T>::right_multiply (const DenseMatrixBase<T2> & M3)
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




template<typename T>
void DenseMatrix<T>::right_multiply_transpose (const DenseMatrix<T> & B)
{
  if (this->use_blas_lapack)
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
          const unsigned int n_rows = B.m();
          const unsigned int n_cols = B.n();

          // resize() *this and also zero out all entries.
          this->resize(n_rows,n_rows);

          // Compute the lower-triangular part
          for (unsigned int i=0; i<n_rows; ++i)
            for (unsigned int j=0; j<=i; ++j)
              for (unsigned int k=0; k<n_cols; ++k) // inner products are over n_cols
                (*this)(i,j) += A(i,k)*A(j,k);

          // Copy lower-triangular part into upper-triangular part
          for (unsigned int i=0; i<n_rows; ++i)
            for (unsigned int j=i+1; j<n_rows; ++j)
              (*this)(i,j) = (*this)(j,i);
        }

      else
        {
          DenseMatrix<T> A(*this);

          this->resize (A.m(), B.m());

          libmesh_assert_equal_to (A.n(), B.n());
          libmesh_assert_equal_to (this->m(), A.m());
          libmesh_assert_equal_to (this->n(), B.m());

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
template<typename T2>
void DenseMatrix<T>::right_multiply_transpose (const DenseMatrix<T2> & B)
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
      const unsigned int n_rows = B.m();
      const unsigned int n_cols = B.n();

      // resize() *this and also zero out all entries.
      this->resize(n_rows,n_rows);

      // Compute the lower-triangular part
      for (unsigned int i=0; i<n_rows; ++i)
        for (unsigned int j=0; j<=i; ++j)
          for (unsigned int k=0; k<n_cols; ++k) // inner products are over n_cols
            (*this)(i,j) += A(i,k)*A(j,k);

      // Copy lower-triangular part into upper-triangular part
      for (unsigned int i=0; i<n_rows; ++i)
        for (unsigned int j=i+1; j<n_rows; ++j)
          (*this)(i,j) = (*this)(j,i);
    }

  else
    {
      DenseMatrix<T> A(*this);

      this->resize (A.m(), B.m());

      libmesh_assert_equal_to (A.n(), B.n());
      libmesh_assert_equal_to (this->m(), A.m());
      libmesh_assert_equal_to (this->n(), B.m());

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




template<typename T>
void DenseMatrix<T>::vector_mult (DenseVector<T> & dest,
                                  const DenseVector<T> & arg) const
{
  // Make sure the input sizes are compatible
  libmesh_assert_equal_to (this->n(), arg.size());

  // Resize and clear dest.
  // Note: DenseVector::resize() also zeros the vector.
  dest.resize(this->m());

  // Short-circuit if the matrix is empty
  if(this->m() == 0 || this->n() == 0)
    return;

  if (this->use_blas_lapack)
    this->_matvec_blas(1., 0., dest, arg);
  else
    {
      const unsigned int n_rows = this->m();
      const unsigned int n_cols = this->n();

      for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
          dest(i) += (*this)(i,j)*arg(j);
    }
}



template<typename T>
template<typename T2>
void DenseMatrix<T>::vector_mult (DenseVector<typename CompareTypes<T,T2>::supertype> & dest,
                                  const DenseVector<T2> & arg) const
{
  // Make sure the input sizes are compatible
  libmesh_assert_equal_to (this->n(), arg.size());

  // Resize and clear dest.
  // Note: DenseVector::resize() also zeros the vector.
  dest.resize(this->m());

  // Short-circuit if the matrix is empty
  if (this->m() == 0 || this->n() == 0)
    return;

  const unsigned int n_rows = this->m();
  const unsigned int n_cols = this->n();

  for (unsigned int i=0; i<n_rows; i++)
    for (unsigned int j=0; j<n_cols; j++)
      dest(i) += (*this)(i,j)*arg(j);
}



template<typename T>
void DenseMatrix<T>::vector_mult_transpose (DenseVector<T> & dest,
                                            const DenseVector<T> & arg) const
{
  // Make sure the input sizes are compatible
  libmesh_assert_equal_to (this->m(), arg.size());

  // Resize and clear dest.
  // Note: DenseVector::resize() also zeros the vector.
  dest.resize(this->n());

  // Short-circuit if the matrix is empty
  if (this->m() == 0)
    return;

  if (this->use_blas_lapack)
    {
      this->_matvec_blas(1., 0., dest, arg, /*trans=*/true);
    }
  else
    {
      const unsigned int n_rows = this->m();
      const unsigned int n_cols = this->n();

      // WORKS
      // for (unsigned int j=0; j<n_cols; j++)
      //   for (unsigned int i=0; i<n_rows; i++)
      //     dest(j) += (*this)(i,j)*arg(i);

      // ALSO WORKS, (i,j) just swapped
      for (unsigned int i=0; i<n_cols; i++)
        for (unsigned int j=0; j<n_rows; j++)
          dest(i) += (*this)(j,i)*arg(j);
    }
}



template<typename T>
template<typename T2>
void DenseMatrix<T>::vector_mult_transpose (DenseVector<typename CompareTypes<T,T2>::supertype> & dest,
                                            const DenseVector<T2> & arg) const
{
  // Make sure the input sizes are compatible
  libmesh_assert_equal_to (this->m(), arg.size());

  // Resize and clear dest.
  // Note: DenseVector::resize() also zeros the vector.
  dest.resize(this->n());

  // Short-circuit if the matrix is empty
  if (this->m() == 0)
    return;

  const unsigned int n_rows = this->m();
  const unsigned int n_cols = this->n();

  // WORKS
  // for (unsigned int j=0; j<n_cols; j++)
  //   for (unsigned int i=0; i<n_rows; i++)
  //     dest(j) += (*this)(i,j)*arg(i);

  // ALSO WORKS, (i,j) just swapped
  for (unsigned int i=0; i<n_cols; i++)
    for (unsigned int j=0; j<n_rows; j++)
      dest(i) += (*this)(j,i)*arg(j);
}



template<typename T>
void DenseMatrix<T>::vector_mult_add (DenseVector<T> & dest,
                                      const T factor,
                                      const DenseVector<T> & arg) const
{
  // Short-circuit if the matrix is empty
  if (this->m() == 0)
    {
      dest.resize(0);
      return;
    }

  if (this->use_blas_lapack)
    this->_matvec_blas(factor, 1., dest, arg);
  else
    {
      DenseVector<T> temp(arg.size());
      this->vector_mult(temp, arg);
      dest.add(factor, temp);
    }
}



template<typename T>
template<typename T2, typename T3>
void DenseMatrix<T>::vector_mult_add (DenseVector<typename CompareTypes<T, typename CompareTypes<T2,T3>::supertype>::supertype> & dest,
                                      const T2 factor,
                                      const DenseVector<T3> & arg) const
{
  // Short-circuit if the matrix is empty
  if (this->m() == 0)
    {
      dest.resize(0);
      return;
    }

  DenseVector<typename CompareTypes<T,T3>::supertype>
    temp(arg.size());
  this->vector_mult(temp, arg);
  dest.add(factor, temp);
}



template<typename T>
void DenseMatrix<T>::get_principal_submatrix (unsigned int sub_m,
                                              unsigned int sub_n,
                                              DenseMatrix<T> & dest) const
{
  libmesh_assert( (sub_m <= this->m()) && (sub_n <= this->n()) );

  dest.resize(sub_m, sub_n);
  for (unsigned int i=0; i<sub_m; i++)
    for (unsigned int j=0; j<sub_n; j++)
      dest(i,j) = (*this)(i,j);
}



template<typename T>
void DenseMatrix<T>::get_principal_submatrix (unsigned int sub_m, DenseMatrix<T> & dest) const
{
  get_principal_submatrix(sub_m, sub_m, dest);
}



template<typename T>
void DenseMatrix<T>::get_transpose (DenseMatrix<T> & dest) const
{
  dest.resize(this->n(), this->m());

  for (unsigned int i=0; i<dest.m(); i++)
    for (unsigned int j=0; j<dest.n(); j++)
      dest(i,j) = (*this)(j,i);
}




template<typename T>
void DenseMatrix<T>::lu_solve (const DenseVector<T> & b,
                               DenseVector<T> & x)
{
  // Check to be sure that the matrix is square before attempting
  // an LU-solve.  In general, one can compute the LU factorization of
  // a non-square matrix, but:
  //
  // Overdetermined systems (m>n) have a solution only if enough of
  // the equations are linearly-dependent.
  //
  // Underdetermined systems (m<n) typically have infinitely many
  // solutions.
  //
  // We don't want to deal with either of these ambiguous cases here...
  libmesh_assert_equal_to (this->m(), this->n());

  switch(this->_decomposition_type)
    {
    case NONE:
      {
        if (this->use_blas_lapack)
          this->_lu_decompose_lapack();
        else
          this->_lu_decompose ();
        break;
      }

    case LU_BLAS_LAPACK:
      {
        // Already factored, just need to call back_substitute.
        if (this->use_blas_lapack)
          break;
      }
      libmesh_fallthrough();

    case LU:
      {
        // Already factored, just need to call back_substitute.
        if (!(this->use_blas_lapack))
          break;
      }
      libmesh_fallthrough();

    default:
      libmesh_error_msg("Error! This matrix already has a different decomposition...");
    }

  if (this->use_blas_lapack)
    this->_lu_back_substitute_lapack (b, x);
  else
    this->_lu_back_substitute (b, x);
}






template<typename T>
void DenseMatrix<T>::_lu_back_substitute (const DenseVector<T> & b,
                                          DenseVector<T> & x ) const
{
  const unsigned int
    n_cols = this->n();

  libmesh_assert_equal_to (this->m(), n_cols);
  libmesh_assert_equal_to (this->m(), b.size());

  x.resize (n_cols);

  // A convenient reference to *this
  const DenseMatrix<T> & A = *this;

  // Temporary vector storage.  We use this instead of
  // modifying the RHS.
  DenseVector<T> z = b;

  // Lower-triangular "top to bottom" solve step, taking into account pivots
  for (unsigned int i=0; i<n_cols; ++i)
    {
      // Swap
      if (_pivots[i] != static_cast<pivot_index_t>(i))
        std::swap( z(i), z(_pivots[i]) );

      x(i) = z(i);

      for (unsigned int j=0; j<i; ++j)
        x(i) -= A(i,j)*x(j);

      x(i) /= A(i,i);
    }

  // Upper-triangular "bottom to top" solve step
  const unsigned int last_row = n_cols-1;

  for (int i=last_row; i>=0; --i)
    {
      for (int j=i+1; j<static_cast<int>(n_cols); ++j)
        x(i) -= A(i,j)*x(j);
    }
}








template<typename T>
void DenseMatrix<T>::_lu_decompose ()
{
  // If this function was called, there better not be any
  // previous decomposition of the matrix.
  libmesh_assert_equal_to (this->_decomposition_type, NONE);

  // Get the matrix size and make sure it is square
  const unsigned int
    n_rows = this->m();

  // A convenient reference to *this
  DenseMatrix<T> & A = *this;

  _pivots.resize(n_rows);

  for (unsigned int i=0; i<n_rows; ++i)
    {
      // Find the pivot row by searching down the i'th column
      _pivots[i] = i;

      // std::abs(complex) must return a Real!
      Real the_max = std::abs( A(i,i) );
      for (unsigned int j=i+1; j<n_rows; ++j)
        {
          Real candidate_max = std::abs( A(j,i) );
          if (the_max < candidate_max)
            {
              the_max = candidate_max;
              _pivots[i] = j;
            }
        }

      // libMesh::out << "the_max=" << the_max << " found at row " << _pivots[i] << std::endl;

      // If the max was found in a different row, interchange rows.
      // Here we interchange the *entire* row, in Gaussian elimination
      // you would only interchange the subrows A(i,j) and A(p(i),j), for j>i
      if (_pivots[i] != static_cast<pivot_index_t>(i))
        {
          for (unsigned int j=0; j<n_rows; ++j)
            std::swap( A(i,j), A(_pivots[i], j) );
        }


      // If the max abs entry found is zero, the matrix is singular
      if (A(i,i) == libMesh::zero)
        libmesh_error_msg("Matrix A is singular!");

      // Scale upper triangle entries of row i by the diagonal entry
      // Note: don't scale the diagonal entry itself!
      const T diag_inv = 1. / A(i,i);
      for (unsigned int j=i+1; j<n_rows; ++j)
        A(i,j) *= diag_inv;

      // Update the remaining sub-matrix A[i+1:m][i+1:m]
      // by subtracting off (the diagonal-scaled)
      // upper-triangular part of row i, scaled by the
      // i'th column entry of each row.  In terms of
      // row operations, this is:
      // for each r > i
      //   SubRow(r) = SubRow(r) - A(r,i)*SubRow(i)
      //
      // If we were scaling the i'th column as well, like
      // in Gaussian elimination, this would 'zero' the
      // entry in the i'th column.
      for (unsigned int row=i+1; row<n_rows; ++row)
        for (unsigned int col=i+1; col<n_rows; ++col)
          A(row,col) -= A(row,i) * A(i,col);

    } // end i loop

  // Set the flag for LU decomposition
  this->_decomposition_type = LU;
}



template<typename T>
void DenseMatrix<T>::svd (DenseVector<Real> & sigma)
{
  // We use the LAPACK svd implementation
  _svd_lapack(sigma);
}


template<typename T>
void DenseMatrix<T>::svd (DenseVector<Real> & sigma,
                          DenseMatrix<Number> & U,
                          DenseMatrix<Number> & VT)
{
  // We use the LAPACK svd implementation
  _svd_lapack(sigma, U, VT);
}



template<typename T>
void DenseMatrix<T>::svd_solve(const DenseVector<T> & rhs,
                               DenseVector<T> & x,
                               Real rcond) const
{
  _svd_solve_lapack(rhs, x, rcond);
}



template<typename T>
void DenseMatrix<T>::evd (DenseVector<T> & lambda_real,
                          DenseVector<T> & lambda_imag)
{
  // We use the LAPACK eigenvalue problem implementation
  _evd_lapack(lambda_real, lambda_imag);
}



template<typename T>
void DenseMatrix<T>::evd_left(DenseVector<T> & lambda_real,
                              DenseVector<T> & lambda_imag,
                              DenseMatrix<T> & VL)
{
  // We use the LAPACK eigenvalue problem implementation
  _evd_lapack(lambda_real, lambda_imag, &VL, libmesh_nullptr);
}



template<typename T>
void DenseMatrix<T>::evd_right(DenseVector<T> & lambda_real,
                               DenseVector<T> & lambda_imag,
                               DenseMatrix<T> & VR)
{
  // We use the LAPACK eigenvalue problem implementation
  _evd_lapack(lambda_real, lambda_imag, libmesh_nullptr, &VR);
}



template<typename T>
void DenseMatrix<T>::evd_left_and_right(DenseVector<T> & lambda_real,
                                        DenseVector<T> & lambda_imag,
                                        DenseMatrix<T> & VL,
                                        DenseMatrix<T> & VR)
{
  // We use the LAPACK eigenvalue problem implementation
  _evd_lapack(lambda_real, lambda_imag, &VL, &VR);
}



template<typename T>
T DenseMatrix<T>::det ()
{
  switch(this->_decomposition_type)
    {
    case NONE:
      {
        // First LU decompose the matrix.
        // Note that the lu_decompose routine will check to see if the
        // matrix is square so we don't worry about it.
        if (this->use_blas_lapack)
          this->_lu_decompose_lapack();
        else
          this->_lu_decompose ();
      }
    case LU:
    case LU_BLAS_LAPACK:
      {
        // Already decomposed, don't do anything
        break;
      }
    default:
      libmesh_error_msg("Error! Can't compute the determinant under the current decomposition.");
    }

  // A variable to keep track of the running product of diagonal terms.
  T determinant = 1.;

  // Loop over diagonal terms, computing the product.  In practice,
  // be careful because this value could easily become too large to
  // fit in a double or float.  To be safe, one should keep track of
  // the power (of 10) of the determinant in a separate variable
  // and maintain an order 1 value for the determinant itself.
  unsigned int n_interchanges = 0;
  for (unsigned int i=0; i<this->m(); i++)
    {
      if (this->_decomposition_type==LU)
        if (_pivots[i] != static_cast<pivot_index_t>(i))
          n_interchanges++;

      // Lapack pivots are 1-based!
      if (this->_decomposition_type==LU_BLAS_LAPACK)
        if (_pivots[i] != static_cast<pivot_index_t>(i+1))
          n_interchanges++;

      determinant *= (*this)(i,i);
    }

  // Compute sign of determinant, depends on number of row interchanges!
  // The sign should be (-1)^{n}, where n is the number of interchanges.
  Real sign = n_interchanges % 2 == 0 ? 1. : -1.;

  return sign*determinant;
}



// The cholesky solve function first decomposes the matrix
// with cholesky_decompose and then uses the cholesky_back_substitute
// routine to find the solution x.
template <typename T>
template <typename T2>
void DenseMatrix<T>::cholesky_solve (const DenseVector<T2> & b,
                                     DenseVector<T2> & x)
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
      libmesh_error_msg("Error! This matrix already has a different decomposition...");
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
  libmesh_assert_equal_to (this->_decomposition_type, NONE);

  // Shorthand notation for number of rows and columns.
  const unsigned int
    n_rows = this->m(),
    n_cols = this->n();

  // Just to be really sure...
  libmesh_assert_equal_to (n_rows, n_cols);

  // A convenient reference to *this
  DenseMatrix<T> & A = *this;

  for (unsigned int i=0; i<n_rows; ++i)
    {
      for (unsigned int j=i; j<n_cols; ++j)
        {
          for (unsigned int k=0; k<i; ++k)
            A(i,j) -= A(i,k) * A(j,k);

          if (i == j)
            {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
              if (A(i,j) <= 0.0)
                libmesh_error_msg("Error! Can only use Cholesky decomposition with symmetric positive definite matrices.");
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
void DenseMatrix<T>::_cholesky_back_substitute (const DenseVector<T2> & b,
                                                DenseVector<T2> & x) const
{
  // Shorthand notation for number of rows and columns.
  const unsigned int
    n_rows = this->m(),
    n_cols = this->n();

  // Just to be really sure...
  libmesh_assert_equal_to (n_rows, n_cols);

  // A convenient reference to *this
  const DenseMatrix<T> & A = *this;

  // Now compute the solution to Ax =b using the factorization.
  x.resize(n_rows);

  // Solve for Ly=b
  for (unsigned int i=0; i<n_cols; ++i)
    {
      T2 temp = b(i);

      for (unsigned int k=0; k<i; ++k)
        temp -= A(i,k)*x(k);

      x(i) = temp / A(i,i);
    }

  // Solve for L^T x = y
  for (unsigned int i=0; i<n_cols; ++i)
    {
      const unsigned int ib = (n_cols-1)-i;

      for (unsigned int k=(ib+1); k<n_cols; ++k)
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
//   // First LU decompose the matrix
//   // Note that the lu_decompose routine will check to see if the
//   // matrix is square so we don't worry about it.
//   if (!this->_lu_decomposed)
//     this->_lu_decompose();

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
// inv(i,j) = x(i);
//     }

//   // Now overwrite all the entries
//   *this = inv;
// }


} // namespace libMesh

#endif // LIBMESH_DENSE_MATRIX_IMPL_H
