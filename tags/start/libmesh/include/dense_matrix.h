// $Id: dense_matrix.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
  DenseMatrix (const DenseMatrix& other_matrix);
  
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
  real operator() (const unsigned int i,
		   const unsigned int j) const;

  /**
   * Assignment operator.
   */
  DenseMatrix& operator = (const DenseMatrix& other_matrix);
  
  /**
   * @returns the \p (i,j) element of the transposed matrix.
   */
  real transpose (const unsigned int i, const unsigned int j) const;
  
  /**
   * @returns the \p (i,j) element of the matrix as a writeable reference.
   */
  real & operator() (const unsigned int i,
		     const unsigned int j);

  /**
   * Multiplies every element in the matrix by \p factor.
   */
  void scale (const real factor);
  
  /**
   * Adds \p factor to every element in the matrix.
   */
  void add (const real factor,
	    const DenseMatrix& mat);

  /**
   * Left multipliess by the matrix \p M.
   * Optionally multiplies by the transpose.
   */
  void left_multiply (const DenseMatrix& A, const bool transpose=false); 
  
  /**
   * Right multiplies by the matrix \p M.
   * Optionally multiplies by the transpose.
   */
  void right_multiply (const DenseMatrix& A, const bool transpose=false); 
  
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
  std::vector<real>& get_values() { return val; };

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,
		const unsigned int j,
		const real val,
		std::vector<real>& rhs);
  
  /**
   * Pretty-print the matrix to \p stdout.
   */
  void print() const;
  
 private:

  /**
   * The actual data values, stored as a 1D array.
   */
  std::vector<real> val;

  /**
   * The row dimension.
   */
  unsigned int m_dim;

  /**
   * The column dimension.
   */
  unsigned int n_dim;

};



/**
 * This class defines a coupling matrix.  A coupling
 * matrix is simply a matrix of ones and zeros describing
 * how different components in a system couple with each
 * other.  A coupling matrix is necessarily square but not
 * necessarily symmetric.
 */
class CouplingMatrix
{
public:
  
  /**
   * Constructor.
   */
  CouplingMatrix (const unsigned int n=0);

  /**
   * @returns the (i,j) entry of the matrix.
   */
  unsigned char operator() (const unsigned int i,
			    const unsigned int j) const;
  
  /**
   * @returns the (i,j) entry of the matrix as
   * a writeable reference.
   */
  unsigned char & operator() (const unsigned int i,
			      const unsigned int j);
  
  /**
   * @returns the size of the matrix, i.e. N for an
   * NxN matrix.
   */ 
  unsigned int size() const;

  /**
   * Resizes the matrix and initializes
   * all entries to be 0.
   */
  void resize(const unsigned int n);

  /**
   * Clears the matrix.
   */
  void clear();

  /**
   * @returns true if the matrix is empty.
   */
  bool empty() const;
  
private:
  
  /**
   * The actual matrix values.  These
   * are stored as unsigned chars because
   * a vector of bools is not what you
   * think.
   */
  std::vector<unsigned char> values;  

  /**
   * The size of the matrix.
   */
  unsigned int _size;
};



// ------------------------------------------------------------
// Dense Matrix member functions
inline
DenseMatrix::DenseMatrix(const unsigned int m,
			 const unsigned int n)
{
  resize(m,n);

  zero();
};



inline
DenseMatrix::DenseMatrix (const DenseMatrix& other_matrix)
{
  m_dim = other_matrix.m_dim;
  n_dim = other_matrix.n_dim;
  
  val = other_matrix.val;
};



inline
DenseMatrix::~DenseMatrix()
{
};



inline
void DenseMatrix::resize(const unsigned int m,
			 const unsigned int n)
{
  if (m*n > val.size())
    val.resize(m*n);

  m_dim = m;
  n_dim = n;

  zero();
  
  return;
};



inline
void DenseMatrix::zero()
{
  for (unsigned int i=0; i<val.size(); i++)
    val[i] = 0.;
};



inline
DenseMatrix& DenseMatrix::operator = (const DenseMatrix& other_matrix)
{
  m_dim = other_matrix.m_dim;
  n_dim = other_matrix.n_dim;

  val   = other_matrix.val;
  
  return *this;
};



inline
real DenseMatrix::operator () (const unsigned int i,
			       const unsigned int j) const
{
  assert (i*j<val.size());
  
  //  return val[(i) + (m_dim)*(j)]; // row-major
  return val[(i)*(n_dim) + (j)]; // col-major
};




inline
real & DenseMatrix::operator () (const unsigned int i,
				 const unsigned int j)
{
  assert (i*j<val.size());
  
  //return val[(i) + (m_dim)*(j)]; // col-major
  return val[(i)*(n_dim) + (j)]; // row-major
};



inline
real DenseMatrix::transpose (const unsigned int i,
			     const unsigned int j) const
{
  assert (i*j<val.size());
  
  //return val[(j) + (m_dim)*(i)]; // col-major
  return val[(j)*(n_dim) + (i)]; // row-major
};
  
    


inline
void DenseMatrix::scale (const real factor)
{
  for (unsigned int i=0; i<val.size(); i++)
    val[i] *= factor;
};




// inline
// void DenseMatrix::vmult(std::vector<real> &w,
// 			std::vector<real> &v,
// 			const bool adding)
// {
//   assert (n() == v.size());
//   assert (m() == w.size());

//   real beta = 0.;

//   if (adding)
//     beta = 1.;
//   /*
//     for (unsigned int i=0; i<m(); i++)
//     {
//     real sum = 0.;

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
    
    


// inline
// void DenseMatrix::Tvmult(std::vector<real> &w,
// 			 std::vector<real> &v,
// 			 const bool adding)
// {
//   assert (m() == v.size());
//   assert (n() == w.size());

//   real beta = 0.;

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




inline
void DenseMatrix::add (const real factor,
		       const DenseMatrix& mat)
{
  assert (m() == mat.m());
  assert (n() == mat.n());

  for (unsigned int j=0; j<n(); j++)
    for (unsigned int i=0; i<m(); i++)
      (*this)(i,j) += factor*mat(i,j);

  return;
};




inline
void DenseMatrix::condense(const unsigned int iv,
			   const unsigned int jv,
			   const real val,
			   std::vector<real>& rhs)
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




inline
void DenseMatrix::print () const
{  
  for (unsigned int i=0; i<m(); i++)
    {
      for (unsigned int j=0; j<n(); j++)
	std::cout << std::setw(8) << (*this)(i,j) << " ";

      std::cout << std::endl;
    }

  return;
};
  



//--------------------------------------------------
// CouplingMatrix inline methods
inline
CouplingMatrix::CouplingMatrix (const unsigned int n) :
  _size(n)
{
  resize(n);
};



inline
unsigned char CouplingMatrix::operator() (const unsigned int i,
					  const unsigned int j) const
{
  assert (i < _size);
  assert (j < _size);

  return values[i*_size + j];
};



inline
unsigned char & CouplingMatrix::operator() (const unsigned int i,
					    const unsigned int j)
{
  assert (i < _size);
  assert (j < _size);

  return values[i*_size + j];
};



inline
unsigned int CouplingMatrix::size() const
{
  return _size;
};



inline
void CouplingMatrix::resize(const unsigned int n)
{
  _size = n;

  values.resize(_size*_size);

  for (unsigned int i=0; i<values.size(); i++)
    values[i] = 0;
};



inline
void CouplingMatrix::clear()
{
  _size = 0;

  values.clear();
};



inline
bool CouplingMatrix::empty() const
{
  return (_size == 0);
};



#endif
