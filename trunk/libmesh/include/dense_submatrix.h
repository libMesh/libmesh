// $Id: dense_submatrix.h,v 1.2 2003-03-07 04:44:38 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __dense_submatrix_h__
#define __dense_submatrix_h__

// C++ includes

// Local Includes
#include "mesh_common.h"
#include "dense_matrix_base.h"
#include "dense_matrix.h"




/**
 * Defines a dense submatrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation
 * into a global matrix, particularly when you have systems of equations.
 *
 * @author Benjamin S. Kirk, 2003
 */ 

// ------------------------------------------------------------
// DenseSubMatrix class definition
template<typename T>
class DenseSubMatrix : public DenseMatrixBase<T>
{
public:

  /**
   * Constructor.  Creates a dense submatrix of the matrix
   * \p parent.  The submatrix has dimensions \f$(m \times n)\f$,
   * and the \f$(0,0) entry of the submatrix is located
   * at the \f$(ioff,joff)\f$ location in the parent matrix.
   */
  DenseSubMatrix(DenseMatrix<T>& parent,
		 const unsigned int ioff=0,
		 const unsigned int joff=0,
		 const unsigned int m=0,
		 const unsigned int n=0);

  /**
   * Copy Constructor.
   */
  DenseSubMatrix (const DenseSubMatrix<T>& other_matrix);
  
  /**
   * Destructor.  Empty.
   */     
  virtual ~DenseSubMatrix() {};

  
  /**
   * Set every element in the submatrix to 0.
   */
  virtual void zero();

  /**
   * @returns the \p (i,j) element of the submatrix.
   */
  virtual T operator() (const unsigned int i,
			const unsigned int j) const;

  /**
   * @returns the \p (i,j) element of the submatrix as a writeable reference.
   */
  virtual T & operator() (const unsigned int i,
			  const unsigned int j);
  
  /**
   * Performs the operation: (*this) <- M2 * (*this) 
   */
  virtual void left_multiply (const DenseMatrixBase<T>& M2);

  /**
   * Performs the operation: (*this) <- (*this) * M3
   */
  virtual void right_multiply (const DenseMatrixBase<T>& M3);
  
  /**
   * Changes the location of the submatrix in the parent matrix. 
   */
  void reposition(const unsigned int ioff,
		  const unsigned int joff,
		  const unsigned int m,
		  const unsigned int n);

  /**
   * @returns the row offset into the parent matrix.
   */
  unsigned int i_off() const { return _i_off; }

  /**
   * @returns the column offset into the parent matrix.
   */
  unsigned int j_off() const { return _j_off; }

private:
  /**
   * The parent matrix that contains this submatrix.
   */
  DenseMatrix<T>& _parent_matrix;
  
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
DenseSubMatrix<T>::DenseSubMatrix(DenseMatrix<T>& parent,
				  const unsigned int ioff,
				  const unsigned int joff,
				  const unsigned int m,
				  const unsigned int n) 
  : DenseMatrixBase<T>(m,n),
    _parent_matrix(parent)
{
  this->reposition (ioff, joff, m, n);
}


// Copy Constructor
template<typename T>
inline
DenseSubMatrix<T>::DenseSubMatrix(const DenseSubMatrix<T>& other_matrix)
  : DenseMatrixBase<T>(other_matrix._m, other_matrix._n),
    _parent_matrix(other_matrix._parent_matrix)
{
  _i_off = other_matrix._i_off; 
  _j_off = other_matrix._j_off; 
}


template<typename T>
inline
void DenseSubMatrix<T>::reposition(const unsigned int ioff,
				   const unsigned int joff,
				   const unsigned int m,
				   const unsigned int n)
{				   
  _i_off = ioff;
  _j_off = joff;
  _m = m;
  _n = n;

  // Make sure we still fit in the parent matrix.
  assert ((this->i_off() + this->m()) <= _parent_matrix.m());
  assert ((this->j_off() + this->n()) <= _parent_matrix.n());
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
  assert (i < this->m());
  assert (j < this->n());
  assert (i + this->i_off() < _parent_matrix.m());
  assert (j + this->j_off() < _parent_matrix.n());
  
  return _parent_matrix (i + this->i_off(),
			 j + this->j_off());
}


template<typename T>
inline
T & DenseSubMatrix<T>::operator () (const unsigned int i,
				    const unsigned int j)
{
  assert (i < this->m());
  assert (j < this->n());
  assert (i + this->i_off() < _parent_matrix.m());
  assert (j + this->j_off() < _parent_matrix.n());
  
  return _parent_matrix (i + this->i_off(),
			 j + this->j_off());
}



#endif // #ifndef __dense_matrix_h__

