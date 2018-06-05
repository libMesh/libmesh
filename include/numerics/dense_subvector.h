// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DENSE_SUBVECTOR_H
#define LIBMESH_DENSE_SUBVECTOR_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_vector.h"

// C++ includes

namespace libMesh
{

/**
 * Defines a dense subvector for use in finite element computations.
 * Useful for storing element load vectors before summation into a
 * global vector, particularly when you have systems of equations.
 * All overridden virtual functions are documented in dense_vector_base.h.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
template<typename T>
class DenseSubVector : public DenseVectorBase<T>
{
public:
  /**
   * Constructor.  Creates a dense subvector of the vector
   * \p parent.  The subvector has dimensions \f$(m \times n)\f$,
   * and the \f$(0,0)\f$ entry of the subvector is located
   * at the \f$(ioff,joff)\f$ location in the parent vector.
   */
  DenseSubVector(DenseVector<T> & new_parent,
                 const unsigned int ioff=0,
                 const unsigned int n=0);

  /**
   * Destructor.  Does nothing.
   */
  virtual ~DenseSubVector() {}

  /**
   * \returns A reference to the parent vector.
   */
  DenseVector<T> & parent () { return _parent_vector; }

  virtual void zero() override;

  /**
   * \returns The \p (i,j) element of the subvector as a const
   * reference.
   */
  const T & operator() (const unsigned int i) const;

  /**
   * \returns The \p (i,j) element of the subvector as a writable reference.
   */
  T & operator() (const unsigned int i);

  virtual T el(const unsigned int i) const override
  { return (*this)(i); }

  virtual T & el(const unsigned int i) override
  { return (*this)(i); }

  virtual unsigned int size() const override
  { return _n; }

  virtual bool empty() const override
  { return (_n == 0); }

  /**
   * \returns The row offset into the parent vector.
   */
  unsigned int i_off() const { return _i_off; }

  /**
   * Changes the location of the subvector in the parent vector.
   */
  void reposition(const unsigned int ioff,
                  const unsigned int n);

  /**
   * \returns The minimum element in the vector, or the minimum real
   * part in the case of complex numbers.
   */
  Real min () const;

  /**
   * \returns The maximum element in the vector, or the maximum real
   * part in the case of complex numbers.
   */
  Real max () const;

  /**
   * \returns The \f$l_1\f$-norm of the vector, i.e. the sum of the
   * absolute values of the entries.
   */
  Real l1_norm () const;

  /**
   * \returns The \f$l_2\f$-norm of the vector, i.e. the square root
   * of the sum of the squares of the entries.
   */
  Real l2_norm () const;

  /**
   * \returns The \f$l_\infty\f$-norm of the vector, i.e. the maximum
   * absolute value of the entries.
   */
  Real linfty_norm () const;

private:


  /**
   * The parent vector that contains this subvector.
   */
  DenseVector<T> & _parent_vector;

  /**
   * The length of this subvector.
   */
  unsigned int _n;

  /**
   * The offset into the parent vector.
   */
  unsigned int _i_off;
};



// ------------------------------------------------------------
// Dense Vector member functions
template<typename T>
inline
DenseSubVector<T>::DenseSubVector(DenseVector<T> & new_parent,
                                  const unsigned int ioff,
                                  const unsigned int n) :
  _parent_vector(new_parent)
{
  reposition (ioff, n);
}



template<typename T>
inline
void DenseSubVector<T>::reposition(const unsigned int ioff,
                                   const unsigned int n)
{
  _i_off = ioff;
  _n = n;

  // Make sure we still fit in the parent vector.
  libmesh_assert_less_equal ((this->i_off() + this->size()), _parent_vector.size());
}



template<typename T>
inline
void DenseSubVector<T>::zero()
{
  for (unsigned int i=0; i<this->size(); i++)
    _parent_vector (i + this->i_off()) = 0.;
}



template<typename T>
inline
const T & DenseSubVector<T>::operator () (const unsigned int i) const
{
  libmesh_assert_less (i, this->size());
  libmesh_assert_less (i + this->i_off(), _parent_vector.size());

  return _parent_vector (i + this->i_off());
}


template<typename T>
inline
T & DenseSubVector<T>::operator () (const unsigned int i)
{
  libmesh_assert_less (i, this->size());
  libmesh_assert_less (i + this->i_off(), _parent_vector.size());

  return _parent_vector (i + this->i_off());
}

template<typename T>
inline
Real DenseSubVector<T>::min () const
{
  libmesh_assert (this->size());
  Real my_min = libmesh_real(_parent_vector (this->i_off()));

  for (unsigned int i=1; i!=this->size(); i++)
    {
      Real current = libmesh_real(_parent_vector (i + this->i_off()));
      my_min = (my_min < current? my_min : current);
    }
  return my_min;
}



template<typename T>
inline
Real DenseSubVector<T>::max () const
{
  libmesh_assert (this->size());
  Real my_max = libmesh_real(_parent_vector (this->i_off()));

  for (unsigned int i=1; i!=this->size(); i++)
    {
      Real current = libmesh_real(_parent_vector (i + this->i_off()));
      my_max = (my_max > current? my_max : current);
    }
  return my_max;
}



template<typename T>
inline
Real DenseSubVector<T>::l1_norm () const
{
  Real my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)
    {
      my_norm += std::abs(_parent_vector (i + this->i_off()));
    }
  return my_norm;
}



template<typename T>
inline
Real DenseSubVector<T>::l2_norm () const
{
  Real my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)
    {
      my_norm += TensorTools::norm_sq(_parent_vector (i + this->i_off()));
    }
  return sqrt(my_norm);
}



template<typename T>
inline
Real DenseSubVector<T>::linfty_norm () const
{
  if (!this->size())
    return 0.;
  Real my_norm = TensorTools::norm_sq(_parent_vector (this->i_off()));

  for (unsigned int i=1; i!=this->size(); i++)
    {
      Real current = TensorTools::norm_sq(_parent_vector (i + this->i_off()));
      my_norm = (my_norm > current? my_norm : current);
    }
  return sqrt(my_norm);
}

} // namespace libMesh


#endif // LIBMESH_DENSE_SUBVECTOR_H
