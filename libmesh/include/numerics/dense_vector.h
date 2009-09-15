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



#ifndef __dense_vector_h__
#define __dense_vector_h__

// C++ includes
#include <vector>

// Local Includes
#include "libmesh_common.h"
#include "compare_types.h"
#include "dense_vector_base.h"

// Forward Declarations



/**
 * Defines a dense vector for use in Finite Element-type computations.
 * This class is to basically compliment the \p DenseMatix class.  It
 * has additional capabilities over the \p std::vector that make it
 * useful for finite elements, particulary for systems of equations.
 *
 * @author Benjamin S. Kirk, 2003
 */ 

// ------------------------------------------------------------
// DenseVector class definition
template<typename T>
class DenseVector : public DenseVectorBase<T>
{
public:

  /**
   * Constructor.  Creates a dense vector of dimension \p n.
   */
  explicit
  DenseVector(const unsigned int n=0);

  /**
   * Copy-constructor.
   */
  template <typename T2>
  DenseVector (const DenseVector<T2>& other_vector);

  /**
   * Copy-constructor, from a \p std::vector.
   */
  template <typename T2>
  DenseVector (const std::vector<T2>& other_vector);
  
  /**
   * Destructor.  Does nothing.
   */     
  ~DenseVector() {}

  /**
   * @returns the size of the vector.
   */
  virtual unsigned int size() const { return _val.size(); }

  /**
   * Set every element in the vector to 0.
   */
  virtual void zero();
  
  /**
   * @returns the \p (i) element of the vector.
   */
  T operator() (const unsigned int i) const;

  /**
   * @returns the \p (i,j) element of the vector as a writeable reference.
   */
  T & operator() (const unsigned int i);

  /**
   * @returns the \p (i) element of the vector.
   */
  virtual T el(const unsigned int i) const { return (*this)(i); }

  /**
   * @returns the \p (i) element of the vector as a writeable reference.
   */
  virtual T & el(const unsigned int i)     { return (*this)(i); }
  
  /**
   * Assignment operator.
   */
  template <typename T2>
  DenseVector<T>& operator = (const DenseVector<T2>& other_vector);
  
  /**
   * STL-like swap method
   */
  void swap(DenseVector<T>& other_vector);
  
  /**
   * Resize the vector. Sets all elements to 0.
   */
  void resize (const unsigned int n);

  /**
   * Multiplies every element in the vector by \p factor.
   */
  void scale (const T factor);
  
  /**
   * Multiplies every element in the vector by \p factor.
   */
  DenseVector<T>& operator*= (const T factor);
  
  /**
   * Adds \p factor times \p vec to this vector.
   * This should only work if T += T2 * T3 is valid C++ and
   * if T2 is scalar.  Return type is void
   */
  template <typename T2, typename T3>
  typename boostcopy::enable_if_c<
    ScalarTraits<T2>::value, void >::type
  add (const T2 factor,
       const DenseVector<T3>& vec);

  /**
   * Evaluate dot product with \p vec.
   */
  template <typename T2>
  Number dot (const DenseVector<T2> &vec) const;

  /**
   * Tests if \p vec is exactly equal to this vector.
   */
  template <typename T2>
  bool operator== (const DenseVector<T2> &vec) const;
  
  /**
   * Tests if \p vec is not exactly equal to this vector.
   */
  template <typename T2>
  bool operator!= (const DenseVector<T2> &vec) const;
  
  /**
   * Adds \p vec to this vector.
   */
  template <typename T2>
  DenseVector<T>& operator+= (const DenseVector<T2> &vec);
  
  /**
   * Subtracts \p vec from this vector.
   */
  template <typename T2>
  DenseVector<T>& operator-= (const DenseVector<T2> &vec);
  
  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  Real min () const;

  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  Real max () const;

  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  Real l1_norm () const;

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  Real l2_norm () const;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
   */
  Real linfty_norm () const;

  /**
   * Access to the values array. This should be used with
   * caution but can  be used to speed up code compilation
   * significantly.
   */
  std::vector<T>& get_values() { return _val; }

  /**
   * Access to the values array. This should be used with
   * caution but can  be used to speed up code compilation
   * significantly.
   */
  const std::vector<T>& get_values() const { return _val; }

private:

  /**
   * The actual data values, stored as a 1D array.
   */
  std::vector<T> _val;

};



// ------------------------------------------------------------
// DenseVector member functions
template<typename T>
inline
DenseVector<T>::DenseVector(const unsigned int n) :
  _val (n, 0.)
{
}



template<typename T>
template<typename T2>
inline
DenseVector<T>::DenseVector (const DenseVector<T2>& other_vector) :
  DenseVectorBase<T>()
{
  const std::vector<T2> &other_vals = other_vector.get_values();

  _val.clear();
  _val.reserve(other_vals.size());

  for (unsigned int i=0; i<other_vals.size(); i++)
    _val.push_back(other_vals[i]);
}



template<typename T>
template<typename T2>
inline
DenseVector<T>::DenseVector (const std::vector<T2>& other_vector) :
  _val(other_vector)
{  
}





template<typename T>
template<typename T2>
inline
DenseVector<T>& DenseVector<T>::operator = (const DenseVector<T2>& other_vector)
{
  //  _val = other_vector._val;

  const std::vector<T2> &other_vals = other_vector.get_values();

  _val.clear();
  _val.reserve(other_vals.size());

  for (unsigned int i=0; i<other_vals.size(); i++)
    _val.push_back(other_vals[i]);

  return *this;
}



template<typename T>
inline
void DenseVector<T>::swap(DenseVector<T>& other_vector)
{
  _val.swap(other_vector._val);
}



template<typename T>
inline
void DenseVector<T>::resize(const unsigned int n)
{
  _val.resize(n);

  zero();
}



template<typename T>
inline
void DenseVector<T>::zero()
{
  std::fill (_val.begin(),
	     _val.end(),
	     0.);
}



template<typename T>
inline
T DenseVector<T>::operator () (const unsigned int i) const
{
  libmesh_assert (i < _val.size());
  
  return _val[i]; 
}



template<typename T>
inline
T & DenseVector<T>::operator () (const unsigned int i)
{
  libmesh_assert (i < _val.size());
  
  return _val[i]; 
}
      


template<typename T>
inline
void DenseVector<T>::scale (const T factor)
{
  for (unsigned int i=0; i<_val.size(); i++)
    _val[i] *= factor;
}



template<typename T>
inline
DenseVector<T>& DenseVector<T>::operator*= (const T factor)
{
  this->scale(factor);
  return *this;
}



template<typename T>
template<typename T2, typename T3>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<T2>::value, void >::type
DenseVector<T>::add (const T2 factor,
		     const DenseVector<T3>& vec)
{
  libmesh_assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    (*this)(i) += factor*vec(i);
}

template<typename T>
template<typename T2>
inline
Number DenseVector<T>::dot (const DenseVector<T2>& vec) const
{
  libmesh_assert(this->size() == vec.size());

  Number val = 0.;

  for (unsigned int i=0; i<this->size(); i++)
    val += (*this)(i)*vec(i);

  return val;
}

template<typename T>
template<typename T2>
inline
bool DenseVector<T>::operator== (const DenseVector<T2>& vec) const
{
  libmesh_assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    if ((*this)(i) != vec(i))
      return false;

  return true;
}



template<typename T>
template<typename T2>
inline
bool DenseVector<T>::operator!= (const DenseVector<T2>& vec) const
{
  libmesh_assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    if ((*this)(i) != vec(i))
      return true;

  return false;
}



template<typename T>
template<typename T2>
inline
DenseVector<T>& DenseVector<T>::operator+= (const DenseVector<T2>& vec)
{
  libmesh_assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    (*this)(i) += vec(i);

  return *this;
}



template<typename T>
template<typename T2>
inline
DenseVector<T>& DenseVector<T>::operator-= (const DenseVector<T2>& vec)
{
  libmesh_assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    (*this)(i) -= vec(i);

  return *this;
}



template<typename T>
inline
Real DenseVector<T>::min () const
{
  libmesh_assert (this->size());
  Real my_min = libmesh_real((*this)(0));

  for (unsigned int i=1; i!=this->size(); i++)
    {
      Real current = libmesh_real((*this)(i));
      my_min = (my_min < current? my_min : current);
    }
  return my_min;
}



template<typename T>
inline
Real DenseVector<T>::max () const
{
  libmesh_assert (this->size());
  Real my_max = libmesh_real((*this)(0));

  for (unsigned int i=1; i!=this->size(); i++)
    {
      Real current = libmesh_real((*this)(i));
      my_max = (my_max > current? my_max : current);
    }
  return my_max;
}



template<typename T>
inline
Real DenseVector<T>::l1_norm () const
{
  Real my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)
    {
      my_norm += std::abs((*this)(i));
    }
  return my_norm;
}



template<typename T>
inline
Real DenseVector<T>::l2_norm () const
{
  Real my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)
    {
      my_norm += libmesh_norm((*this)(i));
    }
  return sqrt(my_norm);
}



template<typename T>
inline
Real DenseVector<T>::linfty_norm () const
{
  if (!this->size())
    return 0.;
  Real my_norm = libmesh_norm((*this)(0));

  for (unsigned int i=1; i!=this->size(); i++)
    {
      Real current = libmesh_norm((*this)(i));
      my_norm = (my_norm > current? my_norm : current);
    }
  return sqrt(my_norm);
}


#endif // #ifndef __dense_vector_h__

