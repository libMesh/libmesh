// $Id: dense_subvector.h,v 1.6 2003-09-25 21:46:55 benkirk Exp $

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



#ifndef __dense_subvector_h__
#define __dense_subvector_h__

// C++ includes

// Local Includes
#include "libmesh_common.h"
#include "dense_vector.h"



/**
 * Defines a dense subvector for use in Finite Element-type computations.
 * Useful for storing element load vectors  before summation
 * into a global vector, particularly when you have systems of equations.
 *
 * @author Benjamin S. Kirk, 2003
 */ 

// ------------------------------------------------------------
// DenseSubVector class definition
template<typename T>
class DenseSubVector : public DenseVectorBase<T>
{
public:

  /**
   * Constructor.  Creates a dense subvector of the vector
   * \p parent.  The subvector has dimensions \f$(m \times n)\f$,
   * and the \f$(0,0) entry of the subvector is located
   * at the \f$(ioff,joff)\f$ location in the parent vector.
   */
  DenseSubVector(DenseVector<T>& parent,
		 const unsigned int ioff=0,
		 const unsigned int n=0);

  /**
   * Destructor.  Does nothing.
   */     
  ~DenseSubVector() {}


  /**
   * Set every element in the subvector to 0.
   */
  virtual void zero();

  /**
   * @returns the \p (i,j) element of the subvector.
   */
  T operator() (const unsigned int i) const;

  /**
   * @returns the \p (i,j) element of the subvector as a writeable reference.
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
   * @returns the size of the subvector.
   */
  virtual unsigned int size() const { return _n; }

  /**
   * @returns the row offset into the parent vector.
   */
  unsigned int i_off() const { return _i_off; }

  /**
   * Changes the location of the subvector in the parent vector. 
   */
  void reposition(const unsigned int ioff,
		  const unsigned int n);
  
private:


  /**
   * The parent vector that contains this subvector.
   */
  DenseVector<T>& _parent_vector;
  
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
DenseSubVector<T>::DenseSubVector(DenseVector<T>& parent,
				  const unsigned int ioff,
				  const unsigned int n) :
  _parent_vector(parent)
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
  assert ((this->i_off() + this->size()) <= _parent_vector.size());
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
T DenseSubVector<T>::operator () (const unsigned int i) const
{
  assert (i < this->size());
  assert (i + this->i_off() < _parent_vector.size());
  
  return _parent_vector (i + this->i_off());
}


template<typename T>
inline
T & DenseSubVector<T>::operator () (const unsigned int i)
{
  assert (i < this->size());
  assert (i + this->i_off() < _parent_vector.size());
  
  return _parent_vector (i + this->i_off());
}



#endif // #ifndef __dense_vector_h__

