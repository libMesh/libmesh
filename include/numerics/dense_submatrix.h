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



#ifndef LIBMESH_DENSE_SUBMATRIX_H
#define LIBMESH_DENSE_SUBMATRIX_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_subvector.h"

// C++ includes

namespace libMesh
{

/**
 * Defines a dense submatrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation
 * into a global matrix, particularly when you have systems of equations.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
template<typename T>
class DenseSubMatrix : public DenseMatrixBase<T>
{
public:

  /**
   * Constructor.  Creates a dense submatrix of the matrix
   * \p parent.  The submatrix has dimensions \f$(m \times n)\f$,
   * and the \f$(0,0)\f$ entry of the submatrix is located
   * at the \f$(ioff,joff)\f$ location in the parent matrix.
   */
  DenseSubMatrix(DenseMatrix<T> & new_parent,
                 const unsigned int ioff=0,
                 const unsigned int joff=0,
                 const unsigned int m=0,
                 const unsigned int n=0);

  /**
   * Copy Constructor.
   */
  DenseSubMatrix (const DenseSubMatrix<T> & other_matrix);

  /**
   * Destructor.  Empty.
   */
  virtual ~DenseSubMatrix() {}


  /**
   * @returns a reference to the parent matrix.
   */
  DenseMatrix<T> & parent () { return _parent_matrix; }

  /**
   * Set every element in the submatrix to 0.
   */
  virtual void zero() libmesh_override;

  /**
   * @returns the \p (i,j) element of the submatrix.
   */
  T operator() (const unsigned int i,
                const unsigned int j) const;

  /**
   * @returns the \p (i,j) element of the submatrix as a writeable reference.
   */
  T & operator() (const unsigned int i,
                  const unsigned int j);

  /**
   * @returns the \p (i,j) element of the matrix as a writeable reference.
   */
  virtual T el(const unsigned int i,
               const unsigned int j) const libmesh_override
  { return (*this)(i,j); }

  /**
   * @returns the \p (i,j) element of the matrix as a writeable reference.
   */
  virtual T & el(const unsigned int i,
                 const unsigned int j) libmesh_override
  { return (*this)(i,j); }

  /**
   * Performs the operation: (*this) <- M2 * (*this)
   */
  virtual void left_multiply (const DenseMatrixBase<T> & M2) libmesh_override;

  /**
   * Performs the operation: (*this) <- (*this) * M3
   */
  virtual void right_multiply (const DenseMatrixBase<T> & M3) libmesh_override;

  /**
   * Changes the location of the submatrix in the parent matrix.
   */
  void reposition(const unsigned int ioff,
                  const unsigned int joff,
                  const unsigned int new_m,
                  const unsigned int new_n);

  /**
   * @returns the row offset into the parent matrix.
   */
  unsigned int i_off() const { return _i_off; }

  /**
   * @returns the column offset into the parent matrix.
   */
  unsigned int j_off() const { return _j_off; }

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,
                const unsigned int j,
                const T val,
                DenseSubVector<T> & rhs)
  {
    this->parent().condense(this->i_off()+i,
                            this->j_off()+j,
                            val, rhs.parent());
  }

private:

  /**
   * The parent matrix that contains this submatrix.
   */
  DenseMatrix<T> & _parent_matrix;

  /**
   * The row offset into the parent matrix.
   */
  unsigned int _i_off;

  /**
   * The column offset into the parent matrix.
   */
  unsigned int _j_off;
};


// --------------------------------------------------
// Constructor
template<typename T>
inline
DenseSubMatrix<T>::DenseSubMatrix(DenseMatrix<T> & new_parent,
                                  const unsigned int ioff,
                                  const unsigned int joff,
                                  const unsigned int new_m,
                                  const unsigned int new_n) :
  DenseMatrixBase<T>(new_m,new_n),
  _parent_matrix(new_parent)
{
  this->reposition (ioff, joff, new_m, new_n);
}


// Copy Constructor
template<typename T>
inline
DenseSubMatrix<T>::DenseSubMatrix(const DenseSubMatrix<T> & other_matrix) :
  DenseMatrixBase<T>(other_matrix._m, other_matrix._n),
  _parent_matrix(other_matrix._parent_matrix)
{
  _i_off = other_matrix._i_off;
  _j_off = other_matrix._j_off;
}


template<typename T>
inline
void DenseSubMatrix<T>::reposition(const unsigned int ioff,
                                   const unsigned int joff,
                                   const unsigned int new_m,
                                   const unsigned int new_n)
{
  _i_off = ioff;
  _j_off = joff;
  this->_m = new_m;
  this->_n = new_n;

  // Make sure we still fit in the parent matrix.
  libmesh_assert_less_equal ((this->i_off() + this->m()), _parent_matrix.m());
  libmesh_assert_less_equal ((this->j_off() + this->n()), _parent_matrix.n());
}



template<typename T>
inline
void DenseSubMatrix<T>::zero()
{
  for (unsigned int i=0; i<this->m(); i++)
    for (unsigned int j=0; j<this->n(); j++)
      _parent_matrix(i + this->i_off(),
                     j + this->j_off()) = 0.;
}



template<typename T>
inline
T DenseSubMatrix<T>::operator () (const unsigned int i,
                                  const unsigned int j) const
{
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());
  libmesh_assert_less (i + this->i_off(), _parent_matrix.m());
  libmesh_assert_less (j + this->j_off(), _parent_matrix.n());

  return _parent_matrix (i + this->i_off(),
                         j + this->j_off());
}


template<typename T>
inline
T & DenseSubMatrix<T>::operator () (const unsigned int i,
                                    const unsigned int j)
{
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());
  libmesh_assert_less (i + this->i_off(), _parent_matrix.m());
  libmesh_assert_less (j + this->j_off(), _parent_matrix.n());

  return _parent_matrix (i + this->i_off(),
                         j + this->j_off());
}


} // namespace libMesh


#endif // LIBMESH_DENSE_SUBMATRIX_H
