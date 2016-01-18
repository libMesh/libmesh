// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DENSE_MATRIX_BASE_H
#define LIBMESH_DENSE_MATRIX_BASE_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/compare_types.h"

// C++ includes

namespace libMesh
{

// Forward Delcarations
template <typename T> class DenseVectorBase;

/**
 * Defines an abstract dense matrix base class for use in Finite Element-type
 * computations.  Specialized dense matrices, for example DenseSubMatrices,
 * can be derived from this class.
 *
 * \author John W. Peterson
 * \date 2003
 */
template<typename T>
class DenseMatrixBase
{

protected:
  /**
   * Constructor.  Creates a dense matrix of dimension \p m by \p n.
   * Protected so that there is no way the user can create one.
   */
  DenseMatrixBase(const unsigned int new_m=0,
                  const unsigned int new_n=0) : _m(new_m), _n(new_n) {}

public:
  /**
   * Destructor. Empty.
   */
  virtual ~DenseMatrixBase() {}

  /**
   * Set every element in the matrix to 0.  You must redefine
   * what you mean by zeroing the matrix since it depends on
   * how your values are stored.
   */
  virtual void zero() = 0;

  /**
   * @returns the \p (i,j) element of the matrix.
   * Since internal data representations may differ, you
   * must redefine this function.
   */
  virtual T el(const unsigned int i,
               const unsigned int j) const = 0;

  /**
   * @returns the \p (i,j) element of the matrix as a writeable reference.
   * Since internal data representations may differ, you
   * must redefine this function.
   */
  virtual T & el(const unsigned int i,
                 const unsigned int j) = 0;

  /**
   * Performs the operation: (*this) <- M2 * (*this)
   */
  virtual void left_multiply (const DenseMatrixBase<T> & M2) = 0;

  /**
   * Performs the operation: (*this) <- (*this) * M3
   */
  virtual void right_multiply (const DenseMatrixBase<T> & M3) = 0;

  /**
   * @returns the row-dimension of the matrix.
   */
  unsigned int m() const { return _m; }

  /**
   * @returns the column-dimension of the matrix.
   */
  unsigned int n() const { return _n; }

  /**
   * Pretty-print the matrix, by default to \p libMesh::out.
   */
  void print(std::ostream & os = libMesh::out) const;

  /**
   * Formatted print as above but allows you to do
   * DenseMatrix K;
   * libMesh::out << K << std::endl;
   */
  friend std::ostream & operator << (std::ostream & os, const DenseMatrixBase<T> & m)
  {
    m.print(os);
    return os;
  }

  /**
   * Prints the matrix entries with more decimal places in
   * scientific notation.
   */
  void print_scientific(std::ostream & os, unsigned precision=8) const;

  /**
   * Adds \p factor to every element in the matrix.
   * This should only work if T += T2 * T3 is valid C++ and
   * if T2 is scalar.  Return type is void
   */
  template <typename T2, typename T3>
  typename boostcopy::enable_if_c<
    ScalarTraits<T2>::value, void >::type
  add (const T2 factor,
       const DenseMatrixBase<T3> & mat);

protected:

  /**
   * Helper function -
   * Performs the computation M1 = M2 * M3 where:
   * M1 = (m x n)
   * M2 = (m x p)
   * M3 = (p x n)
   */
  static void multiply (DenseMatrixBase<T> & M1,
                        const DenseMatrixBase<T> & M2,
                        const DenseMatrixBase<T> & M3);

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,
                const unsigned int j,
                const T val,
                DenseVectorBase<T> & rhs);

  /**
   * The row dimension.
   */
  unsigned int _m;

  /**
   * The column dimension.
   */
  unsigned int _n;
};







template<typename T>
template<typename T2, typename T3>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<T2>::value, void >::type
DenseMatrixBase<T>::add (const T2 factor,
                         const DenseMatrixBase<T3> & mat)
{
  libmesh_assert_equal_to (this->m(), mat.m());
  libmesh_assert_equal_to (this->n(), mat.n());

  for (unsigned int j=0; j<this->n(); j++)
    for (unsigned int i=0; i<this->m(); i++)
      this->el(i,j) += factor*mat.el(i,j);
}


} // namespace libMesh

#endif // LIBMESH_DENSE_MATRIX_BASE_H
