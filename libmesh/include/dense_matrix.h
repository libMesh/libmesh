// $Id: dense_matrix.h,v 1.14 2003-03-07 04:44:38 jwpeterson Exp $

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



#ifndef __dense_matrix_h__
#define __dense_matrix_h__

// C++ includes
#include <vector>

// Local Includes
#include "mesh_common.h"
#include "dense_matrix_base.h"


// Forward Declarations (for friends)




/**
 * Defines a dense matrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation
 * into a global matrix.
 *
 * @author Benjamin S. Kirk, 2002
 */ 

// ------------------------------------------------------------
// Dense Matrix class definition
template<typename T>
class DenseMatrix : public DenseMatrixBase<T>
{
public:
  
  /**
   * Constructor.  Creates a dense matrix of dimension \p m by \p n.
   */
  DenseMatrix(const unsigned int m=0,
	      const unsigned int n=0);
  
  /**
   * Copy-constructor.
   */
  DenseMatrix (const DenseMatrix<T>& other_matrix);
  
  /**
   * Destructor.  Empty.
   */     
  virtual ~DenseMatrix() {}
  

  /**
   * Set every element in the matrix to 0.
   */
  virtual void zero();

  /**
   * @returns the \p (i,j) element of the matrix.
   */
  virtual T operator() (const unsigned int i,
			const unsigned int j) const;

  /**
   * @returns the \p (i,j) element of the matrix as a writeable reference.
   */
  virtual T & operator() (const unsigned int i,
			  const unsigned int j);

  /**
   * Left multipliess by the matrix \p A.
   */
  virtual void left_multiply (const DenseMatrixBase<T>& M2);
  
  /**
   * Right multiplies by the matrix \p A.
   */
  virtual void right_multiply (const DenseMatrixBase<T>& M3);
  
  /**
   * Assignment operator.
   */
  DenseMatrix<T>& operator = (const DenseMatrix<T>& other_matrix);
  
  /**
   * Resize the matrix.  Will never free memory, but may
   * allocate more.  Sets all elements to 0.
   */
  void resize(const unsigned int m,
	      const unsigned int n);
  
  /**
   * Multiplies every element in the matrix by \p factor.
   */
  void scale (const T factor);

  /**
   * Left multiplies by the transpose of the matrix \p A.
   */
  void left_multiply_transpose (const DenseMatrix<T>& A);
  

  /**
   * Right multiplies by the transpose of the matrix \p A
   */
  void right_multiply_transpose (const DenseMatrix<T>& A);
  
  /**
   * @returns the \p (i,j) element of the transposed matrix.
   */
  T transpose (const unsigned int i,
	       const unsigned int j) const;
 
  /**
   * Access to the values array.  This should be used with
   * caution but can  be used to speed up code compilation
   * significantly.
   */
  std::vector<T>& get_values() { return _val; }

  /**
   * Return a constant reference to the matrix values.
   */
  const std::vector<T>& get_values() const { return _val; }

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,
		const unsigned int j,
		const T val,
		std::vector<T>& rhs);
  
  
private:
  /**
   * The actual data values, stored as a 1D array.
   */
  std::vector<T> _val;
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
   * Note that this typedef may be either
   * a real-only matrix, or a truly complex
   * matrix, depending on how \p Number
   * was defined in \p mesh_common.h.
   * Be also aware of the fact that \p DenseMatrix<T>
   * is likely to be more efficient for
   * real than for complex data.
   */
  typedef DenseMatrix<Complex> ComplexDenseMatrix;  

}



using namespace DenseMatrices;





// ------------------------------------------------------------
// Dense Matrix member functions
template<typename T>
inline
DenseMatrix<T>::DenseMatrix(const unsigned int m,
			    const unsigned int n)
  : DenseMatrixBase<T>(m,n)
{
  this->resize(m,n);
}



template<typename T>
inline
DenseMatrix<T>::DenseMatrix (const DenseMatrix<T>& other_matrix)
  : DenseMatrixBase<T>(other_matrix._m, other_matrix._n)
{
  _val = other_matrix._val;
}



template<typename T>
inline
void DenseMatrix<T>::resize(const unsigned int m,
			    const unsigned int n)
{
  if (m*n > _val.size())
    _val.resize(m*n);

  _m = m;
  _n = n;

  this->zero();
}



template<typename T>
inline
void DenseMatrix<T>::zero()
{
  for (unsigned int i=0; i<_val.size(); i++)
    _val[i] = 0.;
}



template<typename T>
inline
DenseMatrix<T>& DenseMatrix<T>::operator = (const DenseMatrix<T>& other_matrix)
{
  _m = other_matrix._m;
  _n = other_matrix._n;

  _val   = other_matrix._val;
  
  return *this;
}



template<typename T>
inline
T DenseMatrix<T>::operator () (const unsigned int i,
			       const unsigned int j) const
{
  assert (i*j<_val.size());
  
  //  return _val[(i) + (_m)*(j)]; // col-major
  return _val[(i)*(_n) + (j)]; // row-major
}



template<typename T>
inline
T & DenseMatrix<T>::operator () (const unsigned int i,
				 const unsigned int j)
{
  assert (i*j<_val.size());
  
  //return _val[(i) + (_m)*(j)]; // col-major
  return _val[(i)*(_n) + (j)]; // row-major
}


  
    


template<typename T>
inline
void DenseMatrix<T>::scale (const T factor)
{
  for (unsigned int i=0; i<_val.size(); i++)
    _val[i] *= factor;
}




template<typename T>
inline
T DenseMatrix<T>::transpose (const unsigned int i,
			     const unsigned int j) const
{
  // Be sure that we haven't exceeded
  // the bounds of the matrix.
  assert (i*j <= _m * _n);
  
  // Implement in terms of operator()
  return (*this)(j,i);
}





template<typename T>
inline
void DenseMatrix<T>::condense(const unsigned int iv,
			      const unsigned int jv,
			      const T val,
			      std::vector<T>& rhs)
{
  assert (_m == rhs.size());
  assert (iv == jv);


  // move the known value into the RHS
  // and zero the column
  for (unsigned int i=0; i<m(); i++)
    {
      rhs[i] -= ((*this)(i,jv))*val;
      (*this)(i,jv) = 0.;
    }

  // zero the row
  for (unsigned int j=0; j<n(); j++)
    (*this)(iv,j) = 0.;

  (*this)(iv,jv) = 1.;
  rhs[iv] = val;
  
}







#endif // #ifndef __dense_matrix_h__

