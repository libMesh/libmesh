// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DENSE_MATRIX_H
#define LIBMESH_DENSE_MATRIX_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_matrix_base.h"
#include "libmesh/int_range.h"

// For the definition of PetscBLASInt.
#if (LIBMESH_HAVE_PETSC)
# include "libmesh/petsc_macro.h"
# ifdef I
#  define LIBMESH_SAW_I
# endif

#include "libmesh/ignore_warnings.h"
# include <petscsys.h>
#include "libmesh/restore_warnings.h"

# ifndef LIBMESH_SAW_I
#  undef I // Avoid complex.h contamination
# endif
#endif

// C++ includes
#include <vector>
#include <algorithm>

#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/dualnumber_forward.h"
#include "metaphysicl/raw_type.h"

namespace std
{
// These declarations must be visible to the DenseMatrix method declarations that use
// a std::abs trailing return type in order to instantiate a DenseMatrix<DualNumber>
template <typename T, typename D, bool asd>
MetaPhysicL::DualNumber<T, D, asd> abs(const MetaPhysicL::DualNumber<T, D, asd> & in);
template <typename T, typename D, bool asd>
MetaPhysicL::DualNumber<T, D, asd> abs(MetaPhysicL::DualNumber<T, D, asd> && in);
}
#endif

namespace libMesh
{

// Forward Declarations
template <typename T> class DenseVector;

/**
 * Defines a dense matrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation into
 * a global matrix.  All overridden virtual functions are documented
 * in dense_matrix_base.h.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A matrix object used for finite element assembly and numerics.
 */
template<typename T>
class DenseMatrix : public DenseMatrixBase<T>
{
public:

  /**
   * Constructor.  Creates a dense matrix of dimension \p m by \p n.
   */
  DenseMatrix(const unsigned int new_m=0,
              const unsigned int new_n=0);

  /**
   * The 5 special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  DenseMatrix (DenseMatrix &&) = default;
  DenseMatrix (const DenseMatrix &) = default;
  DenseMatrix & operator= (const DenseMatrix &) = default;
  DenseMatrix & operator= (DenseMatrix &&) = default;
  virtual ~DenseMatrix() = default;

  /**
   * Sets all elements of the matrix to 0 and resets any decomposition
   * flag which may have been previously set.  This allows e.g. a new
   * LU decomposition to be computed while reusing the same storage.
   */
  virtual void zero() override;

  /**
   * Get submatrix with the smallest row and column indices and the submatrix size.
   */
  DenseMatrix sub_matrix(unsigned int row_id, unsigned int row_size,
                         unsigned int col_id, unsigned int col_size) const;

  /**
   * \returns The \p (i,j) element of the matrix.
   */
  T operator() (const unsigned int i,
                const unsigned int j) const;

  /**
   * \returns The \p (i,j) element of the matrix as a writable reference.
   */
  T & operator() (const unsigned int i,
                  const unsigned int j);

  virtual T el(const unsigned int i,
               const unsigned int j) const override
  { return (*this)(i,j); }

  virtual T & el(const unsigned int i,
                 const unsigned int j) override
  { return (*this)(i,j); }

  virtual void left_multiply (const DenseMatrixBase<T> & M2) override;

  /**
   * Left multiplies by the matrix \p M2 of different type
   */
  template <typename T2>
  void left_multiply (const DenseMatrixBase<T2> & M2);

  virtual void right_multiply (const DenseMatrixBase<T> & M2) override;

  /**
   * Right multiplies by the matrix \p M2 of different type
   */
  template <typename T2>
  void right_multiply (const DenseMatrixBase<T2> & M2);

  /**
   * Performs the matrix-vector multiplication,
   * \p dest := (*this) * \p arg.
   */
  void vector_mult (DenseVector<T> & dest,
                    const DenseVector<T> & arg) const;

  /**
   * Performs the matrix-vector multiplication,
   * \p dest := (*this) * \p arg
   * on mixed types
   */
  template <typename T2>
  void vector_mult (DenseVector<typename CompareTypes<T,T2>::supertype> & dest,
                    const DenseVector<T2> & arg) const;

  /**
   * Performs the matrix-vector multiplication,
   * \p dest := (*this)^T * \p arg.
   */
  void vector_mult_transpose (DenseVector<T> & dest,
                              const DenseVector<T> & arg) const;

  /**
   * Performs the matrix-vector multiplication,
   * \p dest := (*this)^T * \p arg.
   * on mixed types
   */
  template <typename T2>
  void vector_mult_transpose (DenseVector<typename CompareTypes<T,T2>::supertype> & dest,
                              const DenseVector<T2> & arg) const;

  /**
   * Performs the scaled matrix-vector multiplication,
   * \p dest += \p factor * (*this) * \p arg.
   */
  void vector_mult_add (DenseVector<T> & dest,
                        const T factor,
                        const DenseVector<T> & arg) const;

  /**
   * Performs the scaled matrix-vector multiplication,
   * \p dest += \p factor * (*this) * \p arg.
   * on mixed types
   */
  template <typename T2, typename T3>
  void vector_mult_add (DenseVector<typename CompareTypes<T, typename CompareTypes<T2,T3>::supertype>::supertype> & dest,
                        const T2 factor,
                        const DenseVector<T3> & arg) const;

  /**
   * Put the \p sub_m x \p sub_n principal submatrix into \p dest.
   */
  void get_principal_submatrix (unsigned int sub_m, unsigned int sub_n, DenseMatrix<T> & dest) const;

  /**
   * Put the \p sub_m x \p sub_m principal submatrix into \p dest.
   */
  void get_principal_submatrix (unsigned int sub_m, DenseMatrix<T> & dest) const;

  /**
   * Computes the outer (dyadic) product of two vectors and stores in (*this).
   *
   * The outer product of two real-valued vectors \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$ is
   * \f[
   *   (\mathbf{a}\mathbf{b}^T)_{i,j} = \mathbf{a}_i \mathbf{b}_j .
   * \f]
   * The outer product of two complex-valued vectors \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$ is
   * \f[
   *   (\mathbf{a}\mathbf{b}^H)_{i,j} = \mathbf{a}_i \mathbf{b}^*_j ,
   * \f]
   * where \f$H\f$ denotes the conjugate transpose of the vector and \f$*\f$
   * denotes the complex conjugate.
   *
   * \param[in] a   Vector whose entries correspond to rows in the product matrix.
   * \param[in] b   Vector whose entries correspond to columns in the product matrix.
   */
  void outer_product(const DenseVector<T> & a, const DenseVector<T> & b);

  /**
   * Assignment-from-other-matrix-type operator.
   *
   * Copies the dense matrix of type T2 into the present matrix.  This
   * is useful for copying real matrices into complex ones for further
   * operations.
   *
   * \returns A reference to *this.
   */
  template <typename T2>
  DenseMatrix<T> & operator = (const DenseMatrix<T2> & other_matrix);

  /**
   * STL-like swap method
   */
  void swap(DenseMatrix<T> & other_matrix);

  /**
   * Resizes the matrix to the specified size and calls zero().  Will
   * never free memory, but may allocate more. Note: when the matrix
   * is zero()'d, any decomposition (LU, Cholesky, etc.) is also
   * cleared, forcing a new decomposition to be computed the next time
   * e.g. lu_solve() is called.
   */
  void resize(const unsigned int new_m,
              const unsigned int new_n);

  /**
   * Multiplies every element in the matrix by \p factor.
   */
  void scale (const T factor);

  /**
   * Multiplies every element in the column \p col matrix by \p factor.
   */
  void scale_column (const unsigned int col, const T factor);

  /**
   * Multiplies every element in the matrix by \p factor.
   *
   * \returns A reference to *this.
   */
  DenseMatrix<T> & operator *= (const T factor);

  /**
   * Adds \p factor times \p mat to this matrix.
   *
   * \returns A reference to *this.
   */
  template<typename T2, typename T3>
  typename boostcopy::enable_if_c<
    ScalarTraits<T2>::value, void >::type add (const T2 factor,
                                               const DenseMatrix<T3> & mat);

  /**
   * \returns \p true if \p mat is exactly equal to this matrix, \p false otherwise.
   */
  bool operator== (const DenseMatrix<T> & mat) const;

  /**
   * \returns \p true if \p mat is not exactly equal to this matrix, false otherwise.
   */
  bool operator!= (const DenseMatrix<T> & mat) const;

  /**
   * Adds \p mat to this matrix.
   *
   * \returns A reference to *this.
   */
  DenseMatrix<T> & operator+= (const DenseMatrix<T> & mat);

  /**
   * Subtracts \p mat from this matrix.
   *
   * \returns A reference to *this.
   */
  DenseMatrix<T> & operator-= (const DenseMatrix<T> & mat);

  /**
   * \returns The minimum entry in the matrix, or the minimum real
   * part in the case of complex numbers.
   */
  auto min () const -> decltype(libmesh_real(T(0)));

  /**
   * \returns The maximum entry in the matrix, or the maximum real
   * part in the case of complex numbers.
   */
  auto max () const -> decltype(libmesh_real(T(0)));

  /**
   * \returns The l1-norm of the matrix, that is, the max column sum:
   *
   * \f$ |M|_1 = max_{all columns j} \sum_{all rows i} |M_ij| \f$,
   *
   * This is the natural matrix norm that is compatible to the l1-norm
   * for vectors, i.e. \f$ |Mv|_1 \leq |M|_1 |v|_1 \f$.
   */
  auto l1_norm () const -> decltype(std::abs(T(0)));

  /**
   * \returns The linfty-norm of the matrix, that is, the max row sum:
   *
   * \f$ |M|_\infty = max_{all rows i} \sum_{all columns j} |M_ij| \f$,
   *
   * This is the natural matrix norm that is compatible to the
   * linfty-norm of vectors, i.e. \f$ |Mv|_\infty \leq |M|_\infty |v|_\infty \f$.
   */
  auto linfty_norm () const -> decltype(std::abs(T(0)));

  /**
   * Left multiplies by the transpose of the matrix \p A.
   */
  void left_multiply_transpose (const DenseMatrix<T> & A);

  /**
   * Left multiplies by the transpose of the matrix \p A which
   * contains a different numerical type.
   */
  template <typename T2>
  void left_multiply_transpose (const DenseMatrix<T2> & A);


  /**
   * Right multiplies by the transpose of the matrix \p A
   */
  void right_multiply_transpose (const DenseMatrix<T> & A);

  /**
   * Right multiplies by the transpose of the matrix \p A which
   * contains a different numerical type.
   */
  template <typename T2>
  void right_multiply_transpose (const DenseMatrix<T2> & A);

  /**
   * \returns The \p (i,j) element of the transposed matrix.
   */
  T transpose (const unsigned int i,
               const unsigned int j) const;

  /**
   * Put the tranposed matrix into \p dest.
   */
  void get_transpose(DenseMatrix<T> & dest) const;

  /**
   * \returns A reference to the underlying data storage vector.
   *
   * This should be used with caution (i.e. one should not change the
   * size of the vector, etc.) but is useful for interoperating with
   * low level BLAS routines which expect a simple array.
   */
  std::vector<T> & get_values() { return _val; }

  /**
   * \returns A constant reference to the underlying data storage vector.
   */
  const std::vector<T> & get_values() const { return _val; }

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,
                const unsigned int j,
                const T val,
                DenseVector<T> & rhs)
  { DenseMatrixBase<T>::condense (i, j, val, rhs); }

  /**
   * Solve the system Ax=b given the input vector b.  Partial pivoting
   * is performed by default in order to keep the algorithm stable to
   * the effects of round-off error.
   *
   * Important note: once you call lu_solve(), you must _not_ modify
   * the entries of the matrix via calls to operator(i,j) and call
   * lu_solve() again without first calling either zero() or resize(),
   * otherwise the code will skip computing the decomposition of the
   * matrix and go directly to the back substitution step. This is
   * done on purpose for efficiency, so that the same LU decomposition
   * can be used with multiple right-hand sides, but it does also make
   * it possible to "shoot yourself in the foot", so be careful!
   */
  void lu_solve (const DenseVector<T> & b,
                 DenseVector<T> & x);

  /**
   * For symmetric positive definite (SPD) matrices. A Cholesky factorization
   * of A such that A = L L^T is about twice as fast as a standard LU
   * factorization.  Therefore you can use this method if you know a-priori
   * that the matrix is SPD.  If the matrix is not SPD, an error is generated.
   * One nice property of Cholesky decompositions is that they do not require
   * pivoting for stability.
   *
   * Important note: once you call cholesky_solve(), you must _not_
   * modify the entries of the matrix via calls to operator(i,j) and
   * call cholesky_solve() again without first calling either zero()
   * or resize(), otherwise the code will skip computing the
   * decomposition of the matrix and go directly to the back
   * substitution step. This is done on purpose for efficiency, so
   * that the same decomposition can be used with multiple right-hand
   * sides, but it does also make it possible to "shoot yourself in
   * the foot", so be careful!
   *
   * \note This method may also be used when A is real-valued and x
   * and b are complex-valued.
   */
  template <typename T2>
  void cholesky_solve(const DenseVector<T2> & b,
                      DenseVector<T2> & x);

  /**
   * Compute the singular value decomposition of the matrix.
   * On exit, sigma holds all of the singular values (in
   * descending order).
   *
   * The implementation uses PETSc's interface to BLAS/LAPACK.
   * If this is not available, this function throws an error.
   */
  void svd(DenseVector<Real> & sigma);

  /**
   * Compute the "reduced" singular value decomposition of the matrix.
   * On exit, sigma holds all of the singular values (in
   * descending order), U holds the left singular vectors,
   * and VT holds the transpose of the right singular vectors.
   * In the reduced SVD, U has min(m,n) columns and VT has
   * min(m,n) rows. (In the "full" SVD, U and VT would be square.)
   *
   * The implementation uses PETSc's interface to BLAS/LAPACK.
   * If this is not available, this function throws an error.
   */
  void svd(DenseVector<Real> & sigma,
           DenseMatrix<Number> & U,
           DenseMatrix<Number> & VT);

  /**
   * Solve the system of equations \f$ A x = rhs \f$ for \f$ x \f$ in the
   * least-squares sense. \f$ A \f$ may be non-square and/or rank-deficient.
   * You can control which singular values are treated as zero by
   * changing the "rcond" parameter.  Singular values S(i) for which
   * S(i) <= rcond*S(1) are treated as zero for purposes of the solve.
   * Passing a negative number for rcond forces a "machine precision"
   * value to be used instead.
   *
   * This function is marked const, since due to various
   * implementation details, we do not need to modify the contents of
   * A in order to compute the SVD (a copy is made internally
   * instead).
   *
   * Requires PETSc >= 3.1 since this was the first version to provide
   * the LAPACKgelss_ wrapper.
   */
  void svd_solve(const DenseVector<T> & rhs,
                 DenseVector<T> & x,
                 Real rcond=std::numeric_limits<Real>::epsilon()) const;

  /**
   * Compute the eigenvalues (both real and imaginary parts) of a general matrix.
   *
   * Warning: the contents of \p *this are overwritten by this function!
   *
   * The implementation requires the LAPACKgeev_ function which is wrapped by PETSc.
   */
  void evd(DenseVector<T> & lambda_real,
           DenseVector<T> & lambda_imag);

  /**
   * Compute the eigenvalues (both real and imaginary parts) and left
   * eigenvectors of a general matrix, \f$ A \f$.
   *
   * Warning: the contents of \p *this are overwritten by this function!
   *
   * The left eigenvector \f$ u_j \f$ of \f$ A \f$ satisfies:
   * \f$ u_j^H A = lambda_j u_j^H \f$
   * where \f$ u_j^H \f$ denotes the conjugate-transpose of \f$ u_j \f$.
   *
   * If the j-th and (j+1)-st eigenvalues form a complex conjugate
   * pair, then the j-th and (j+1)-st columns of VL "share" their
   * real-valued storage in the following way:
   * u_j     = VL(:,j) + i*VL(:,j+1) and
   * u_{j+1} = VL(:,j) - i*VL(:,j+1).
   *
   * The implementation requires the LAPACKgeev_ routine which is provided by PETSc.
   */
  void evd_left(DenseVector<T> & lambda_real,
                DenseVector<T> & lambda_imag,
                DenseMatrix<T> & VL);

  /**
   * Compute the eigenvalues (both real and imaginary parts) and right
   * eigenvectors of a general matrix, \f$ A \f$.
   *
   * Warning: the contents of \p *this are overwritten by this function!
   *
   * The right eigenvector \f$ v_j \f$ of \f$ A \f$ satisfies:
   * \f$ A v_j = lambda_j v_j \f$
   * where \f$ lambda_j \f$ is its corresponding eigenvalue.
   *
   * \note If the j-th and (j+1)-st eigenvalues form a complex
   * conjugate pair, then the j-th and (j+1)-st columns of VR "share"
   * their real-valued storage in the following way:
   * v_j     = VR(:,j) + i*VR(:,j+1) and
   * v_{j+1} = VR(:,j) - i*VR(:,j+1).
   *
   * The implementation requires the LAPACKgeev_ routine which is provided by PETSc.
   */
  void evd_right(DenseVector<T> & lambda_real,
                 DenseVector<T> & lambda_imag,
                 DenseMatrix<T> & VR);

  /**
   * Compute the eigenvalues (both real and imaginary parts) as well as the left
   * and right eigenvectors of a general matrix.
   *
   * Warning: the contents of \p *this are overwritten by this function!
   *
   * See the documentation of the \p evd_left() and \p evd_right()
   * functions for more information.  The implementation requires the
   * LAPACKgeev_ routine which is provided by PETSc.
   */
  void evd_left_and_right(DenseVector<T> & lambda_real,
                          DenseVector<T> & lambda_imag,
                          DenseMatrix<T> & VL,
                          DenseMatrix<T> & VR);

  /**
   * \returns The determinant of the matrix.
   *
   * \note Implemented by computing an LU decomposition and then
   * taking the product of the diagonal terms.  Therefore this is a
   * non-const method which modifies the entries of the matrix.
   */
  T det();

  /**
   * Computes the inverse of the dense matrix (assuming it is invertible)
   * by first computing the LU decomposition and then performing multiple
   * back substitution steps.  Follows the algorithm from Numerical Recipes
   * in C that is available on the web.
   *
   * This routine is commented out since it is not really a memory- or
   * computationally- efficient implementation.  Also, you typically
   * don't need the actual inverse for anything, and can use something
   * like lu_solve() instead.
   */
  // void inverse();

  /**
   * Run-time selectable option to turn on/off BLAS support.
   * This was primarily used for testing purposes, and could be
   * removed...
   */
  bool use_blas_lapack;

  /**
   * Helper structure for determining whether to use blas_lapack
   */
  struct UseBlasLapack
  {
    static const bool value = false;
  };

private:

  /**
   * The actual data values, stored as a 1D array.
   */
  std::vector<T> _val;

  /**
   * Form the LU decomposition of the matrix.  This function
   * is private since it is only called as part of the implementation
   * of the lu_solve(...) function.
   */
  void _lu_decompose ();

  /**
   * Solves the system Ax=b through back substitution.  This function
   * is private since it is only called as part of the implementation
   * of the lu_solve(...) function.
   */
  void _lu_back_substitute (const DenseVector<T> & b,
                            DenseVector<T> & x) const;

  /**
   * Decomposes a symmetric positive definite matrix into a
   * product of two lower triangular matrices according to
   * A = LL^T.
   *
   * \note This program generates an error if the matrix is not SPD.
   */
  void _cholesky_decompose();

  /**
   * Solves the equation Ax=b for the unknown value x and rhs
   * b based on the Cholesky factorization of A.
   *
   * \note This method may be used when A is real-valued and b and x
   * are complex-valued.
   */
  template <typename T2>
  void _cholesky_back_substitute(const DenseVector<T2> & b,
                                 DenseVector<T2> & x) const;

  /**
   * The decomposition schemes above change the entries of the matrix
   * A.  It is therefore an error to call A.lu_solve() and subsequently
   * call A.cholesky_solve() since the result will probably not match
   * any desired outcome.  This typedef keeps track of which decomposition
   * has been called for this matrix.
   */
  enum DecompositionType {LU=0, CHOLESKY=1, LU_BLAS_LAPACK, NONE};

  /**
   * This flag keeps track of which type of decomposition has been
   * performed on the matrix.
   */
  DecompositionType _decomposition_type;

  /**
   * Enumeration used to determine the behavior of the _multiply_blas
   * function.
   */
  enum _BLAS_Multiply_Flag {
    LEFT_MULTIPLY = 0,
    RIGHT_MULTIPLY,
    LEFT_MULTIPLY_TRANSPOSE,
    RIGHT_MULTIPLY_TRANSPOSE
  };

  /**
   * The _multiply_blas function computes A <- op(A) * op(B) using
   * BLAS gemm function.  Used in the right_multiply(),
   * left_multiply(), right_multiply_transpose(), and
   * left_multiply_transpose() routines.
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _multiply_blas(const DenseMatrixBase<T> & other,
                      _BLAS_Multiply_Flag flag);

  /**
   * Computes an LU factorization of the matrix using the
   * Lapack routine "getrf".  This routine should only be
   * used by the "use_blas_lapack" branch of the lu_solve()
   * function.  After the call to this function, the matrix
   * is replaced by its factorized version, and the
   * DecompositionType is set to LU_BLAS_LAPACK.
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _lu_decompose_lapack();

  /**
   * Computes an SVD of the matrix using the
   * Lapack routine "getsvd".
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _svd_lapack(DenseVector<Real> & sigma);

  /**
   * Computes a "reduced" SVD of the matrix using the
   * Lapack routine "getsvd".
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _svd_lapack(DenseVector<Real> & sigma,
                   DenseMatrix<Number> & U,
                   DenseMatrix<Number> & VT);

  /**
   * Called by svd_solve(rhs).
   */
  void _svd_solve_lapack(const DenseVector<T> & rhs,
                         DenseVector<T> & x,
                         Real rcond) const;

  /**
   * Helper function that actually performs the SVD.
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _svd_helper (char JOBU,
                    char JOBVT,
                    std::vector<Real> & sigma_val,
                    std::vector<Number> & U_val,
                    std::vector<Number> & VT_val);

  /**
   * Computes the eigenvalues of the matrix using the Lapack routine
   * "DGEEV".  If VR and/or VL are not nullptr, then the matrix of right
   * and/or left eigenvectors is also computed and returned by this
   * function.
   *
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _evd_lapack(DenseVector<T> & lambda_real,
                   DenseVector<T> & lambda_imag,
                   DenseMatrix<T> * VL = nullptr,
                   DenseMatrix<T> * VR = nullptr);

  /**
   * Array used to store pivot indices.  May be used by whatever
   * factorization is currently active, clients of the class should
   * not rely on it for any reason.
   */
#if (LIBMESH_HAVE_PETSC && LIBMESH_USE_REAL_NUMBERS)
  typedef PetscBLASInt pivot_index_t;
#else
  typedef int pivot_index_t;
#endif
  std::vector<pivot_index_t> _pivots;

  /**
   * Companion function to _lu_decompose_lapack().  Do not use
   * directly, called through the public lu_solve() interface.
   * This function is logically const in that it does not modify
   * the matrix, but since we are just calling LAPACK routines,
   * it's less const_cast hassle to just declare the function
   * non-const.
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _lu_back_substitute_lapack (const DenseVector<T> & b,
                                   DenseVector<T> & x);

  /**
   * Uses the BLAS GEMV function (through PETSc) to compute
   *
   * dest := alpha*A*arg + beta*dest
   *
   * where alpha and beta are scalars, A is this matrix, and
   * arg and dest are input vectors of appropriate size.  If
   * trans is true, the transpose matvec is computed instead.
   * By default, trans==false.
   *
   * [ Implementation in dense_matrix_blas_lapack.C ]
   */
  void _matvec_blas(T alpha, T beta,
                    DenseVector<T> & dest,
                    const DenseVector<T> & arg,
                    bool trans=false) const;
};





// ------------------------------------------------------------
/**
 * Provide Typedefs for dense matrices
 */
namespace DenseMatrices
{

/**
 * Convenient definition of a real-only
 * dense matrix.
 */
typedef DenseMatrix<Real> RealDenseMatrix;

/**
 * This typedef may be either a real-only matrix, or a truly complex
 * matrix, depending on how \p Number was defined in \p
 * libmesh_common.h.  Also, be aware of the fact that \p
 * DenseMatrix<T> is likely to be more efficient for real than for
 * complex data.
 */
typedef DenseMatrix<Complex> ComplexDenseMatrix;

}



using namespace DenseMatrices;

// The PETSc Lapack wrappers are only for PetscScalar, therefore we
// can't e.g. get a Lapack version of DenseMatrix<Real>::lu_solve()
// when libmesh/PETSc are compiled with complex numbers.
#if defined(LIBMESH_HAVE_PETSC) && \
  defined(LIBMESH_USE_REAL_NUMBERS) && \
  defined(LIBMESH_DEFAULT_DOUBLE_PRECISION)
template <>
struct DenseMatrix<double>::UseBlasLapack
{
  static const bool value = true;
};
#endif

// ------------------------------------------------------------
// Dense Matrix member functions
template<typename T>
inline
DenseMatrix<T>::DenseMatrix(const unsigned int new_m,
                            const unsigned int new_n) :
  DenseMatrixBase<T>(new_m,new_n),
  use_blas_lapack(DenseMatrix<T>::UseBlasLapack::value),
  _val(),
  _decomposition_type(NONE)
{
  this->resize(new_m,new_n);
}



template<typename T>
inline
void DenseMatrix<T>::swap(DenseMatrix<T> & other_matrix)
{
  std::swap(this->_m, other_matrix._m);
  std::swap(this->_n, other_matrix._n);
  _val.swap(other_matrix._val);
  DecompositionType _temp = _decomposition_type;
  _decomposition_type = other_matrix._decomposition_type;
  other_matrix._decomposition_type = _temp;
}


template <typename T>
template <typename T2>
inline
DenseMatrix<T> &
DenseMatrix<T>::operator=(const DenseMatrix<T2> & mat)
{
  unsigned int mat_m = mat.m(), mat_n = mat.n();
  this->resize(mat_m, mat_n);
  for (unsigned int i=0; i<mat_m; i++)
    for (unsigned int j=0; j<mat_n; j++)
      (*this)(i,j) = mat(i,j);

  return *this;
}



template<typename T>
inline
void DenseMatrix<T>::resize(const unsigned int new_m,
                            const unsigned int new_n)
{
  _val.resize(new_m*new_n);

  this->_m = new_m;
  this->_n = new_n;

  // zero and set decomposition_type to NONE
  this->zero();
}



template<typename T>
inline
void DenseMatrix<T>::zero()
{
  _decomposition_type = NONE;

  std::fill (_val.begin(), _val.end(), static_cast<T>(0));
}



template<typename T>
inline
DenseMatrix<T> DenseMatrix<T>::sub_matrix(unsigned int row_id, unsigned int row_size,
                                          unsigned int col_id, unsigned int col_size) const
{
  libmesh_assert_less (row_id + row_size - 1, this->_m);
  libmesh_assert_less (col_id + col_size - 1, this->_n);

  DenseMatrix<T> sub;
  sub._m = row_size;
  sub._n = col_size;
  sub._val.resize(row_size * col_size);

  unsigned int end_col = this->_n - col_size - col_id;
  unsigned int p = row_id * this->_n;
  unsigned int q = 0;
  for (unsigned int i=0; i<row_size; i++)
  {
    // skip the beginning columns
    p += col_id;
    for (unsigned int j=0; j<col_size; j++)
      sub._val[q++] = _val[p++];
    // skip the rest columns
    p += end_col;
  }

  return sub;
}



template<typename T>
inline
T DenseMatrix<T>::operator () (const unsigned int i,
                               const unsigned int j) const
{
  libmesh_assert_less (i*j, _val.size());
  libmesh_assert_less (i, this->_m);
  libmesh_assert_less (j, this->_n);


  //  return _val[(i) + (this->_m)*(j)]; // col-major
  return _val[(i)*(this->_n) + (j)]; // row-major
}



template<typename T>
inline
T & DenseMatrix<T>::operator () (const unsigned int i,
                                 const unsigned int j)
{
  libmesh_assert_less (i*j, _val.size());
  libmesh_assert_less (i, this->_m);
  libmesh_assert_less (j, this->_n);

  //return _val[(i) + (this->_m)*(j)]; // col-major
  return _val[(i)*(this->_n) + (j)]; // row-major
}





template<typename T>
inline
void DenseMatrix<T>::scale (const T factor)
{
  for (auto & v : _val)
    v *= factor;
}


template<typename T>
inline
void DenseMatrix<T>::scale_column (const unsigned int col, const T factor)
{
  for (auto i : make_range(this->m()))
    (*this)(i, col) *= factor;
}



template<typename T>
inline
DenseMatrix<T> & DenseMatrix<T>::operator *= (const T factor)
{
  this->scale(factor);
  return *this;
}



template<typename T>
template<typename T2, typename T3>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<T2>::value, void >::type
DenseMatrix<T>::add (const T2 factor,
                     const DenseMatrix<T3> & mat)
{
  libmesh_assert_equal_to (this->m(), mat.m());
  libmesh_assert_equal_to (this->n(), mat.n());

  for (auto i : make_range(this->m()))
    for (auto j : make_range(this->n()))
      (*this)(i,j) += factor * mat(i,j);
}



template<typename T>
inline
bool DenseMatrix<T>::operator == (const DenseMatrix<T> & mat) const
{
  for (auto i : index_range(_val))
    if (_val[i] != mat._val[i])
      return false;

  return true;
}



template<typename T>
inline
bool DenseMatrix<T>::operator != (const DenseMatrix<T> & mat) const
{
  for (auto i : index_range(_val))
    if (_val[i] != mat._val[i])
      return true;

  return false;
}



template<typename T>
inline
DenseMatrix<T> & DenseMatrix<T>::operator += (const DenseMatrix<T> & mat)
{
  for (auto i : index_range(_val))
    _val[i] += mat._val[i];

  return *this;
}



template<typename T>
inline
DenseMatrix<T> & DenseMatrix<T>::operator -= (const DenseMatrix<T> & mat)
{
  for (auto i : index_range(_val))
    _val[i] -= mat._val[i];

  return *this;
}



template<typename T>
inline
auto DenseMatrix<T>::min () const -> decltype(libmesh_real(T(0)))
{
  libmesh_assert (this->_m);
  libmesh_assert (this->_n);
  auto my_min = libmesh_real((*this)(0,0));

  for (unsigned int i=0; i!=this->_m; i++)
    {
      for (unsigned int j=0; j!=this->_n; j++)
        {
          auto current = libmesh_real((*this)(i,j));
          my_min = (my_min < current? my_min : current);
        }
    }
  return my_min;
}



template<typename T>
inline
auto DenseMatrix<T>::max () const -> decltype(libmesh_real(T(0)))
{
  libmesh_assert (this->_m);
  libmesh_assert (this->_n);
  auto my_max = libmesh_real((*this)(0,0));

  for (unsigned int i=0; i!=this->_m; i++)
    {
      for (unsigned int j=0; j!=this->_n; j++)
        {
          auto current = libmesh_real((*this)(i,j));
          my_max = (my_max > current? my_max : current);
        }
    }
  return my_max;
}



template<typename T>
inline
auto DenseMatrix<T>::l1_norm () const -> decltype(std::abs(T(0)))
{
  libmesh_assert (this->_m);
  libmesh_assert (this->_n);

  auto columnsum = std::abs(T(0));
  for (unsigned int i=0; i!=this->_m; i++)
    {
      columnsum += std::abs((*this)(i,0));
    }
  auto my_max = columnsum;
  for (unsigned int j=1; j!=this->_n; j++)
    {
      columnsum = 0.;
      for (unsigned int i=0; i!=this->_m; i++)
        {
          columnsum += std::abs((*this)(i,j));
        }
      my_max = (my_max > columnsum? my_max : columnsum);
    }
  return my_max;
}



template<typename T>
inline
auto DenseMatrix<T>::linfty_norm () const -> decltype(std::abs(T(0)))
{
  libmesh_assert (this->_m);
  libmesh_assert (this->_n);

  auto rowsum = std::abs(T(0));
  for (unsigned int j=0; j!=this->_n; j++)
    {
      rowsum += std::abs((*this)(0,j));
    }
  auto my_max = rowsum;
  for (unsigned int i=1; i!=this->_m; i++)
    {
      rowsum = 0.;
      for (unsigned int j=0; j!=this->_n; j++)
        {
          rowsum += std::abs((*this)(i,j));
        }
      my_max = (my_max > rowsum? my_max : rowsum);
    }
  return my_max;
}



template<typename T>
inline
T DenseMatrix<T>::transpose (const unsigned int i,
                             const unsigned int j) const
{
  // Implement in terms of operator()
  return (*this)(j,i);
}





// template<typename T>
// inline
// void DenseMatrix<T>::condense(const unsigned int iv,
//       const unsigned int jv,
//       const T val,
//       DenseVector<T> & rhs)
// {
//   libmesh_assert_equal_to (this->_m, rhs.size());
//   libmesh_assert_equal_to (iv, jv);


//   // move the known value into the RHS
//   // and zero the column
//   for (auto i : make_range(this->m()))
//     {
//       rhs(i) -= ((*this)(i,jv))*val;
//       (*this)(i,jv) = 0.;
//     }

//   // zero the row
//   for (auto j : make_range(this->n()))
//     (*this)(iv,j) = 0.;

//   (*this)(iv,jv) = 1.;
//   rhs(iv) = val;

// }




} // namespace libMesh

#ifdef LIBMESH_HAVE_METAPHYSICL
namespace MetaPhysicL
{
template <typename T>
struct RawType<libMesh::DenseMatrix<T>>
{
  typedef libMesh::DenseMatrix<typename RawType<T>::value_type> value_type;

  static value_type value (const libMesh::DenseMatrix<T> & in)
    {
      value_type ret(in.m(), in.n());
      for (unsigned int i = 0; i < in.m(); ++i)
        for (unsigned int j = 0; j < in.n(); ++j)
          ret(i,j) = raw_value(in(i,j));

      return ret;
    }
};
}
#endif


#endif // LIBMESH_DENSE_MATRIX_H
