// $Id: dense_matrix.h,v 1.7 2003-02-07 18:07:45 ddreyer Exp $

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
#include <iomanip>
#include <vector>

// Local Includes
#include "mesh_common.h"


/**
 * Defines a dense matrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation
 * into a global matrix.
 *
 * @author Benjamin S. Kirk, 2002
 */ 

// ------------------------------------------------------------
// Dense Matrix class definition
template<typename Tp>
class DenseMatrix
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
  DenseMatrix (const DenseMatrix<Tp>& other_matrix);
  
  /**
   * Destructor.  Frees all associated memory.
   */     
  ~DenseMatrix();

  /**
   * Resize the matrix.  Will never free memory, but may
   * allocate more.  Sets all elements to 0.
   */
  void resize(const unsigned int m,
	      const unsigned int n);

  /**
   * Set every element in the matrix to 0.
   */
  void zero();

  /**
   * @returns the \p (i,j) element of the matrix.
   */
  Tp operator() (const unsigned int i,
		 const unsigned int j) const;

  /**
   * Assignment operator.
   */
  DenseMatrix<Tp>& operator = (const DenseMatrix<Tp>& other_matrix);
  
  /**
   * @returns the \p (i,j) element of the transposed matrix.
   */
  Tp transpose (const unsigned int i, const unsigned int j) const;
  
  /**
   * @returns the \p (i,j) element of the matrix as a writeable reference.
   */
  Tp & operator() (const unsigned int i,
		   const unsigned int j);

  /**
   * Multiplies every element in the matrix by \p factor.
   */
  void scale (const Tp factor);
  
  /**
   * Adds \p factor to every element in the matrix.
   */
  void add (const Tp factor,
	    const DenseMatrix<Tp>& mat);


#ifdef USE_COMPLEX_NUMBERS

  /**
   * For a complex-valued matrix, let a real-valued
   * matrix being added to us (Note that the other
   * way around would be wrong!). 
   */
  void add (const Complex factor,
	    const DenseMatrix<Real>& mat);

#endif


  /**
   * Left multipliess by the matrix \p M.
   * Optionally multiplies by the transpose.
   */
  void left_multiply (const DenseMatrix<Tp>& A, const bool transpose=false); 
  
  /**
   * Right multiplies by the matrix \p M.
   * Optionally multiplies by the transpose.
   */
  void right_multiply (const DenseMatrix<Tp>& A, const bool transpose=false); 
  
  /**
   * @returns the row-dimension of the matrix.
   */
  unsigned int m() const { return m_dim; };
  
  /**
   * @returns the column-dimension of the matrix.
   */
  unsigned int n() const { return n_dim; };

  /**
   * Access to the values array.  This should be used with
   * caution but can  be used to speed up code compilation
   * significantly.
   */
  std::vector<Tp>& get_values() { return val; };

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,
		const unsigned int j,
		const Tp val,
		std::vector<Tp>& rhs);
  
  /**
   * Pretty-print the matrix to \p stdout.
   */
  void print() const;
  
 private:

  /**
   * The actual data values, stored as a 1D array.
   */
  std::vector<Tp> val;

  /**
   * The row dimension.
   */
  unsigned int m_dim;

  /**
   * The column dimension.
   */
  unsigned int n_dim;

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
   * matrix, depending on how \p Complex
   * was defined in \p mesh_common.h.
   * Be also aware of the fact that \p DenseMatrix<Tp>
   * is likely to be more efficient for
   * real than for complex data.
   */
  typedef DenseMatrix<Complex> ComplexDenseMatrix;  

};



using namespace DenseMatrices;





// ------------------------------------------------------------
// Dense Matrix member functions
template<typename Tp>
inline
DenseMatrix<Tp>::DenseMatrix(const unsigned int m,
			     const unsigned int n)
{
  resize(m,n);

  zero();
};



template<typename Tp>
inline
DenseMatrix<Tp>::DenseMatrix (const DenseMatrix& other_matrix)
{
  m_dim = other_matrix.m_dim;
  n_dim = other_matrix.n_dim;
  
  val = other_matrix.val;
};



template<typename Tp>
inline
DenseMatrix<Tp>::~DenseMatrix()
{
};



template<typename Tp>
inline
void DenseMatrix<Tp>::resize(const unsigned int m,
			     const unsigned int n)
{
  if (m*n > val.size())
    val.resize(m*n);

  m_dim = m;
  n_dim = n;

  zero();
  
  return;
};



template<typename Tp>
inline
void DenseMatrix<Tp>::zero()
{
  for (unsigned int i=0; i<val.size(); i++)
    val[i] = 0.;
};



template<typename Tp>
inline
DenseMatrix<Tp>& DenseMatrix<Tp>::operator = (const DenseMatrix<Tp>& other_matrix)
{
  m_dim = other_matrix.m_dim;
  n_dim = other_matrix.n_dim;

  val   = other_matrix.val;
  
  return *this;
};



template<typename Tp>
inline
Tp DenseMatrix<Tp>::operator () (const unsigned int i,
				 const unsigned int j) const
{
  assert (i*j<val.size());
  
  //  return val[(i) + (m_dim)*(j)]; // row-major
  return val[(i)*(n_dim) + (j)]; // col-major
};



template<typename Tp>
inline
Tp & DenseMatrix<Tp>::operator () (const unsigned int i,
				   const unsigned int j)
{
  assert (i*j<val.size());
  
  //return val[(i) + (m_dim)*(j)]; // col-major
  return val[(i)*(n_dim) + (j)]; // row-major
};



template<typename Tp>
inline
Tp DenseMatrix<Tp>::transpose (const unsigned int i,
			       const unsigned int j) const
{
  assert (i*j<val.size());
  
  //return val[(j) + (m_dim)*(i)]; // col-major
  return val[(j)*(n_dim) + (i)]; // row-major
};
  
    


template<typename Tp>
inline
void DenseMatrix<Tp>::scale (const Tp factor)
{
  for (unsigned int i=0; i<val.size(); i++)
    val[i] *= factor;
};




// template<typename Tp>
// inline
// void DenseMatrix<Tp>::vmult(std::vector<Tp> &w,
// 			std::vector<Tp> &v,
// 			const bool adding)
// {
//   assert (n() == v.size());
//   assert (m() == w.size());

//   Tp beta = 0.;

//   if (adding)
//     beta = 1.;
//   /*
//     for (unsigned int i=0; i<m(); i++)
//     {
//     Tp sum = 0.;

//     for (unsigned int j=0; j<n(); j++)
//     sum += val[i + m_dim*j]*v[j];

//     w[i] = beta*w[i] + sum;
//     }
//   */
// #ifdef HAVE_BLAS
//   error();
//   //  cblas_dgemv(CblasColMajor, CblasNoTrans, m(), n(), 1.0, val, m(),
//   //	      &v[0], 1, beta, &w[0], 1);
// #else
//   error();
  
// #endif

// };
    
    

// template<typename Tp>
// inline
// void DenseMatrix<Tp>::Tvmult(std::vector<Tp> &w,
// 			 std::vector<Tp> &v,
// 			 const bool adding)
// {
//   assert (m() == v.size());
//   assert (n() == w.size());

//   Tp beta = 0.;

//   if (adding)
//     beta = 1.;
  
// #ifdef HAVE_BLAS
//   error(); 
//   //  cblas_dgemv(CblasColMajor, CblasTrans, m(), n(), 1.0, val, m(),
//   //	      &v[0], 1, beta, &w[0], 1);

// #else

//   error();

// #endif
// };




template<typename Tp>
inline
void DenseMatrix<Tp>::add (const Tp factor,
			   const DenseMatrix<Tp>& mat)
{
  assert (m() == mat.m());
  assert (n() == mat.n());

  for (unsigned int j=0; j<n(); j++)
    for (unsigned int i=0; i<m(); i++)
      (*this)(i,j) += factor*mat(i,j);

  return;
};



#ifdef USE_COMPLEX_NUMBERS

/*
 * For complex numbers, also offer a method
 * to add a real-valued matrix to a complex-
 * valued matrix (but not the other way around)!
 */

template<>
inline
void DenseMatrix<Complex>::add (const Complex factor,
				const DenseMatrix<Real>& mat)
{
  assert (m() == mat.m());
  assert (n() == mat.n());

  for (unsigned int j=0; j<n(); j++)
    for (unsigned int i=0; i<m(); i++)
      (*this)(i,j) += factor*mat(i,j);

  return;
};


#endif


template<typename Tp>
inline
void DenseMatrix<Tp>::condense(const unsigned int iv,
			       const unsigned int jv,
			       const Tp val,
			       std::vector<Tp>& rhs)
{
  assert (m() == rhs.size());
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
  
};




template<typename Tp>
inline
void DenseMatrix<Tp>::print () const
{  
  for (unsigned int i=0; i<m(); i++)
    {
      for (unsigned int j=0; j<n(); j++)
	std::cout << std::setw(8) << (*this)(i,j) << " ";

      std::cout << std::endl;
    }

  return;
};
  





#endif
