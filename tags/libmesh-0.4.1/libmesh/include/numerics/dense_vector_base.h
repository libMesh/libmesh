// $Id: dense_vector_base.h,v 1.1 2003-11-05 22:26:44 benkirk Exp $

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



#ifndef __dense_vector_base_h__
#define __dense_vector_base_h__

// C++ includes
#include <iomanip> // for std::setw

// Local Includes
#include "libmesh_common.h"


// Forward Declarations



/**
 * Defines an abstract dense vector base class for use in
 * Finite Element-type computations.
 * Specialized dense vectors, for example DenseSubVectors, can
 * be derived from this class.
 *
 * @author John W. Peterson, 2003
 */ 

// ------------------------------------------------------------
// DenseVectorBase class definition
template<typename T>
class DenseVectorBase
{
public:

  /**
   * Constructor.  Empty.
   */
  DenseVectorBase() {}
  
  /**
   * Destructor.  Does nothing.
   */     
  virtual ~DenseVectorBase() {}

  /**
   * Set every element in the vector to 0.  Needs to
   * be pure virtual since the storage method may
   * be different in derived classes.
   */
  virtual void zero() = 0;

  /**
   * @returns the \p (i) element of the vector.
   */
  virtual T el(const unsigned int i) const = 0;

  /**
   * @returns the \p (i) element of the vector as a writeable reference.
   */
  virtual T & el(const unsigned int i) = 0;

  /**
   * @returns the size of the vector.
   */
  virtual unsigned int size() const = 0; 

  /**
   * Pretty-print the vector to \p stdout.
   */
  void print() const;
};



// ------------------------------------------------------------
// DenseVectorBase member functions



template<typename T>
inline
void DenseVectorBase<T>::print () const
{  
  for (unsigned int i=0; i<this->size(); i++)
    std::cout << std::setw(8) << this->el(i) << std::endl;
}



#endif // #ifndef __dense_vector_base_h__

