// $Id: laspack_matrix.h,v 1.1 2003-02-07 22:18:53 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __laspack_matrix_h__
#define __laspack_matrix_h__

#include "mesh_config.h"

#ifdef HAVE_LASPACK

// C++ includes

// Local includes
#include "sparse_matrix.h"


namespace Laspack {
extern "C" {
#include <qmatrix.h>
}
}






/**
 * Generic laspack matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to laspack matrices
 * from various different solver packages in
 * different formats.
 *
 * @author Benjamin S. Kirk, 2003
 */

class LaspackMatrix : public SparseMatrix
{

 public:
  /**
   * Constructor; initializes the matrix to
   * be empty, without any structure, i.e.
   * the matrix is not usable at all. This
   * constructor is therefore only useful
   * for matrices which are members of a
   * class. All other matrices should be
   * created at a point in the data flow
   * where all necessary information is
   * available.
   *
   * You have to initialize
   * the matrix before usage with
   * \p init(...).
   */
  LaspackMatrix ();

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~LaspackMatrix ();

  /**
   * Initialize a Laspack matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   */
  void init (const unsigned int m,
	     const unsigned int n,
	     const unsigned int m_l,
	     const unsigned int n_l,
	     const unsigned int nnz=30,
	     const unsigned int noz=10);

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */   
  void init (const DofMap& dof_map);
  
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  void clear ();

  /**
   * Set all entries to 0.
   */
  void zero ();
  
  /**
   * Call the Laspack assemble routines.
   * sends necessary messages to other
   * processors
   */
  void close () const;
  
  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  unsigned int m () const;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  unsigned int n () const;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  unsigned int row_start () const;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  unsigned int row_stop () const;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  void set (const unsigned int i,
		    const unsigned int j,
		    const Complex value);
    
  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  void add (const unsigned int i,
		    const unsigned int j,
		    const Complex value);

  /**
   * Add the full matrix to the
   * Laspack matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
    
  void add_matrix (const ComplexDenseMatrix &dm,
			   const std::vector<unsigned int> &rows,
			   const std::vector<unsigned int> &cols);
  
  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  void add_matrix (const ComplexDenseMatrix &dm,
			   const std::vector<unsigned int> &dof_indices);
      
  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation and you
   * should always take care where
   * to call this function.  In
   * order to avoid abuse, this
   * function throws an exception
   * if the required element does
   * not exist in the matrix.
   *
   * In case you want a function
   * that returns zero instead (for
   * entries that are not in the
   * sparsity pattern of the
   * matrix), use the \p el
   * function.
   */
  Complex operator () (const unsigned int i,
			       const unsigned int j) const;

  /**
   * Return the l1-norm of the matrix, that is
   * $|M|_1=max_{all columns j}\sum_{all 
   * rows i} |M_ij|$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * $|Mv|_1\leq |M|_1 |v|_1$.
   */
  Real l1_norm () const;

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * $|M|_infty=max_{all rows i}\sum_{all 
   * columns j} |M_ij|$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * $|Mv|_infty \leq |M|_infty |v|_infty$.
   */
  Real linfty_norm () const;

  /**
   * see if Laspack matrix has been closed
   * and fully assembled yet
   */
  bool closed() const;

  
private:
  
  /**
   *  The Laspack sparse matrix.
   */
  Laspack::QMatrix QMat;
};



//-----------------------------------------------------------------------
// LaspackMatrix class inline members
inline
LaspackMatrix::LaspackMatrix ()
{
  using namespace Laspack;
  
  Q_Constr(&QMat, ((char*) "Mat"), 0, False, Rowws, Normal, True);
};



inline
LaspackMatrix::~LaspackMatrix ()
{
  using namespace Laspack;
  
  Q_Destr(&QMat);
};



inline
unsigned int LaspackMatrix::m () const
{
  assert (initialized());

  return static_cast<unsigned int>(Q_GetDim(const_cast<Laspack::QMatrix*>(&QMat)));
};



inline
unsigned int LaspackMatrix::n () const
{
  assert (initialized());
  
  return static_cast<unsigned int>(Q_GetDim(const_cast<Laspack::QMatrix*>(&QMat)));
};



inline
unsigned int LaspackMatrix::row_start () const
{
  return 0;
};



inline
unsigned int LaspackMatrix::row_stop () const
{
  return n();
};



inline
void LaspackMatrix::set (const unsigned int i,
			 const unsigned int j,
			 const Complex value)
{
  assert (initialized());

  using namespace Laspack;

  error();
};



inline
void LaspackMatrix::add (const unsigned int i,
			 const unsigned int j,
			 const Complex value)
{
  assert (initialized());

  using namespace Laspack;

  error();
};



inline
void LaspackMatrix::add_matrix(const ComplexDenseMatrix& dm,
			       const std::vector<unsigned int>& dof_indices)
{
  add_matrix (dm, dof_indices, dof_indices);
};



inline
void LaspackMatrix::add_matrix(const ComplexDenseMatrix& dm,
			       const std::vector<unsigned int>& rows,
			       const std::vector<unsigned int>& cols)
		    
{
  assert (initialized());
  assert (dm.m() == rows.size());
  assert (dm.n() == cols.size());

  
  for (unsigned int i=0; i<rows.size(); i++)
    for (unsigned int j=0; j<cols.size(); j++)
      add(i,j,dm(i,j));
};



inline
Complex LaspackMatrix::operator () (const unsigned int i,
				    const unsigned int j) const
{
  assert (initialized());
  assert (i < m());
  assert (j < n());
  
  using namespace Laspack;

  return Q_GetEl (const_cast<Laspack::QMatrix*>(&QMat), i, j);
};



inline
bool LaspackMatrix::closed () const
{
  return true;
};


#endif // #ifdef HAVE_LASPACK


#endif // #ifdef __laspack_matrix_h__
