// $Id: dense_matrix_base.h,v 1.9 2004-08-17 03:03:49 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __dense_matrix_base_h__
#define __dense_matrix_base_h__

// C++ includes
#include <iostream>
#include <iomanip> // for std::setw()

// Local Includes
#include "libmesh_common.h"






/**
 * Defines an abstract dense matrix base class for use in Finite Element-type
 * computations.  Specialized dense matrices, for example DenseSubMatrices,
 * can be derived from this class.
 *
 * @author John W. Peterson, 2003
 */ 
template<typename T>
class DenseMatrixBase
{
  
protected:
  /**
   * Constructor.  Creates a dense matrix of dimension \p m by \p n.
   * Protected so that there is no way the user can create one.
   */
  DenseMatrixBase(const unsigned int m=0,
		  const unsigned int n=0) : _m(m), _n(n) {};
  
public:
  /**
   * Destructor. Empty.
   */     
  virtual ~DenseMatrixBase() {};

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
  virtual void left_multiply (const DenseMatrixBase<T>& M2) = 0;

  /**
   * Performs the operation: (*this) <- (*this) * M3
   */
  virtual void right_multiply (const DenseMatrixBase<T>& M3) = 0;
  
  /**
   * @returns the row-dimension of the matrix.
   */
  unsigned int m() const { return _m; }
  
  /**
   * @returns the column-dimension of the matrix.
   */
  unsigned int n() const { return _n; }
  
  /**
   * Pretty-print the matrix to \p stdout.
   */
  void print() const;

  /**
   * Formatted print as above but allows you to do
   * DenseMatrix K;
   * std::cout << K << std::endl;
   */
  friend std::ostream& operator << (std::ostream& os, const DenseMatrixBase<T>& m)
  {
    m.print();
    return os;
  }
  
  /**
   * Prints the matrix entries with more decimal places in
   * scientific notation.
   */
  void print_scientific() const;
  
  /**
   * Adds \p factor to every element in the matrix.
   */
  void add (const T factor,
	    const DenseMatrixBase<T>& mat);
  
    
#ifdef USE_COMPLEX_NUMBERS
  
  /**
   * For a complex-valued matrix, let a real-valued
   * matrix being added to us (Note that the other
   * way around would be wrong!). 
   */
  void add (const Complex factor,
	    const DenseMatrixBase<Real>& mat);
  
#endif
   
protected:
  /**
   * Performs the computation M1 = M2 * M3 where:
   * M1 = (m x n)
   * M2 = (m x p)
   * M3 = (p x n)
   */
  void multiply (DenseMatrixBase<T>& M1,
		 const DenseMatrixBase<T>& M2,
		 const DenseMatrixBase<T>& M3);

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
inline
void DenseMatrixBase<T>::print_scientific () const
{
#ifndef BROKEN_IOSTREAM
  
  // save the initial format flags
  std::ios_base::fmtflags cout_flags = std::cout.flags();
  
  // Print the matrix entries.
  for (unsigned int i=0; i<this->m(); i++)
    {
      for (unsigned int j=0; j<this->n(); j++)
	std::cout << std::setw(15)
		  << std::scientific
		  << std::setprecision(8)
		  << this->el(i,j) << " ";

      std::cout << std::endl;
    }
  
  // reset the original format flags
  std::cout.flags(cout_flags);

#else
  
  // Print the matrix entries.
  for (unsigned int i=0; i<this->m(); i++)
    {
      for (unsigned int j=0; j<this->n(); j++)	
	std::cout << std::setprecision(8)
		  << this->el(i,j)
		  << " ";
      
      std::cout << std::endl;
    }
  
  
#endif
}



template<typename T>
inline
void DenseMatrixBase<T>::print () const
{  
  for (unsigned int i=0; i<this->m(); i++)
    {
      for (unsigned int j=0; j<this->n(); j++)
	std::cout << std::setw(8)
		  << this->el(i,j) << " ";

      std::cout << std::endl;
    }

  return;
}



template<typename T>
inline
void DenseMatrixBase<T>::add (const T factor,
			      const DenseMatrixBase<T>& mat)
{
  assert (this->m() == mat.m());
  assert (this->n() == mat.n());

  for (unsigned int j=0; j<this->n(); j++)
    for (unsigned int i=0; i<this->m(); i++)
      this->el(i,j) += factor*mat.el(i,j);
}



#ifdef USE_COMPLEX_NUMBERS

/*
 * For complex numbers, also offer a method
 * to add a real-valued matrix to a complex-de
 * valued matrix (but not the other way around)!
 */
template<>
inline
void DenseMatrixBase<Complex>::add (const Complex factor,
				    const DenseMatrixBase<Real>& mat)
{
  assert (this->m() == mat.m());
  assert (this->n() == mat.n());

  for (unsigned int j=0; j<this->n(); j++)
    for (unsigned int i=0; i<this->m(); i++)
      this->el(i,j) += factor*mat.el(i,j);
}


#endif







#endif // #ifndef __dense_matrix_base_h__

