// $Id: dense_vector.h,v 1.6 2003-09-02 18:02:37 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_common.h"
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
  DenseVector(const unsigned int n=0);

  /**
   * Copy-constructor.
   */
  DenseVector (const DenseVector<T>& other_vector);

  /**
   * Copy-constructor, from a \p std::vector.
   */
  DenseVector (const std::vector<T>& other_vector);
  
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
  DenseVector<T>& operator = (const DenseVector<T>& other_vector);
  
  /**
   * Resize the vector. Sets all elements to 0.
   */
  void resize (const unsigned int n);

  /**
   * Multiplies every element in the vector by \p factor.
   */
  void scale (const T factor);
  
  /**
   * Adds \p factor to every element in the vector.
   */
  void add (const T factor,
	    const DenseVector<T>& vec);


#ifdef USE_COMPLEX_NUMBERS

  /**
   * For a complex-valued vector, let a real-valued
   * vector being added to us (Note that the other
   * way around would be wrong!). 
   */
  void add (const Complex factor,
	    const DenseVector<Real>& vec);

#endif


  /**
   * Access to the values array. This should be used with
   * caution but can  be used to speed up code compilation
   * significantly.
   */
  std::vector<T>& get_values() { return _val; }

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
inline
DenseVector<T>::DenseVector (const DenseVector<T>& other_vector) :
  DenseVectorBase<T>(),
  _val(other_vector._val)
{  
}



template<typename T>
inline
DenseVector<T>::DenseVector (const std::vector<T>& other_vector) :
  _val(other_vector)
{  
}





template<typename T>
inline
DenseVector<T>& DenseVector<T>::operator = (const DenseVector<T>& other_vector)
{
  _val = other_vector._val;
  
  return *this;
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
  assert (i < _val.size());
  
  return _val[i]; 
}



template<typename T>
inline
T & DenseVector<T>::operator () (const unsigned int i)
{
  assert (i < _val.size());
  
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
void DenseVector<T>::add (const T factor,
			  const DenseVector<T>& vec)
{
  assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    (*this)(i) += factor*vec(i);
}



#ifdef USE_COMPLEX_NUMBERS

/*
 * For complex numbers, also offer a method
 * to add a real-valued vector to a complex-
 * valued vector (but not the other way around)!
 */
template<>
inline
void DenseVector<Complex>::add (const Complex factor,
				const DenseVector<Real>& vec)
{
  assert (this->size() == vec.size());

  for (unsigned int i=0; i<this->size(); i++)
    (*this)(i) += factor*vec(i);
}

#endif





#endif // #ifndef __dense_vector_h__

