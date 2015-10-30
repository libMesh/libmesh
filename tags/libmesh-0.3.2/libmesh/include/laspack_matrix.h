// $Id: laspack_matrix.h,v 1.8 2003-02-20 23:18:06 benkirk Exp $

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
#include <algorithm>

// Local includes
#include "sparse_matrix.h"


namespace Laspack {
#include <qmatrix.h>
}


// Forward declarations
class LaspackInterface;



/**
 * Generic laspack matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to laspack matrices
 * from various different solver packages in
 * different formats.
 * Currently Laspack only supports real datatypes, so
 * this class is a full specialization of \p SparseMatrix<>
 * with \p T = \p Real

 *
 * @author Benjamin S. Kirk, 2003
 */

class LaspackMatrix : public SparseMatrix<Real>
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
   * Updates the matrix sparsity pattern. This will
   * tell the underlying matrix storage scheme how
   * to map the \f$ (i,j) \f$ elements.
   */
  void update_sparsity_pattern (std::vector<std::set<unsigned int> >&);
  
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
  void init ();
  
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
   * Close the matrix.  Dummy routine.  After calling
   * this method \p closed() is true and the matrix can
   * be used in computations.
   */
  void close () const { const_cast<LaspackMatrix*>(this)->_closed = true; }
  
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
	    const Real value);
    
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
	    const Real value);

  /**
   * Add the full matrix to the
   * Laspack matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
    
  void add_matrix (const DenseMatrix<Real> &dm,
		   const std::vector<unsigned int> &rows,
		   const std::vector<unsigned int> &cols);
  
  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  void add_matrix (const DenseMatrix<Real> &dm,
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
  Real operator () (const unsigned int i,
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
  Real l1_norm () const { error(); return 0.; }

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
  Real linfty_norm () const { error(); return 0.; }

  /**
   * see if Laspack matrix has been closed
   * and fully assembled yet
   */
  bool closed() const { return _closed; }

  
private:

  /**
   * @returns the position in the compressed row
   * storage scheme of the \f$ (i,j) \f$ element.
   */
  unsigned int pos (const unsigned int i,
		    const unsigned int j) const;
  
  /**
   *  The Laspack sparse matrix pointer.
   */
  Laspack::QMatrix _QMat;

  /**
   * The compressed row indices.
   */
  std::vector<unsigned int> _csr;

  /**
   * The start of each row in the compressed
   * row index data structure.
   */
  std::vector<unsigned int> _row_start;

  /**
   * Flag indicating if the matrix has been closed yet.
   */
  bool _closed;

  /**
   * Make other Laspack datatypes friends
   */
  friend class LaspackInterface;
};



//-----------------------------------------------------------------------
// LaspackMatrix class inline members
inline
LaspackMatrix::LaspackMatrix () :
  _closed (false)
{
}



inline
LaspackMatrix::~LaspackMatrix ()
{
  clear ();
}



inline
void LaspackMatrix::clear ()
{
  if (initialized())
    {
      Laspack::Q_Destr(&_QMat);
    }
  
  _csr.clear();
  _row_start.clear();
  _closed = false;
  _is_initialized = false;
}



inline
void LaspackMatrix::zero ()
{
  using namespace Laspack;

  const unsigned int n_rows = m();
  
  for (unsigned int row=0; row<n_rows; row++)
    {
      const unsigned int r_start = _row_start[row];
      const unsigned int len     = Q_GetLen(&_QMat, row+1);

      //       std::cout << "row=" << row << ", \t"
      // 		<< "len=" << len << ", \t"
      // 		<< "_row_start[row+1]-_row_start[row]="
      // 		<<  _row_start[row+1]-_row_start[row]
      // 		<< std::endl;
	
      // Make sure we agree on the row length
      assert (len == (_row_start[row+1]-_row_start[row]));
      
      for (unsigned int l=0; l<len; l++)
	{
	  const unsigned int j = _csr[r_start + l];

	  // Make sure the data structures are working
	  assert ((j+1) == Q_GetPos (&_QMat, row+1, l));
	  
	  Q_SetEntry (&_QMat, row+1, l, j+1, 0.);
	}
    }    
}



inline
unsigned int LaspackMatrix::m () const
{
  assert (initialized());

  return static_cast<unsigned int>(Laspack::Q_GetDim(const_cast<Laspack::QMatrix*>(&_QMat)));
}



inline
unsigned int LaspackMatrix::n () const
{
  assert (initialized());
  
  return static_cast<unsigned int>(Laspack::Q_GetDim(const_cast<Laspack::QMatrix*>(&_QMat)));
}



inline
unsigned int LaspackMatrix::row_start () const
{
  return 0;
}



inline
unsigned int LaspackMatrix::row_stop () const
{
  return m();
}



inline
void LaspackMatrix::set (const unsigned int i,
			 const unsigned int j,
			 const Real value)
{
  assert (initialized());
  assert (i < m());
  assert (j < n());
  
  const unsigned int position = pos(i,j);

  // Sanity check
  assert (_csr[_row_start[i]+position] == j);
  assert ((j+1) == Laspack::Q_GetPos (&_QMat, i+1, position));

  Laspack::Q_SetEntry (&_QMat, i+1, position, j+1, value);
}



inline
void LaspackMatrix::add (const unsigned int i,
			 const unsigned int j,
			 const Real value)
{
  assert (initialized());
  assert (i < m());
  assert (j < n());
  
  const unsigned int position = pos(i,j);

  // Sanity check
  assert (_csr[_row_start[i]+position] == j);

  Laspack::Q_AddVal (&_QMat, i+1, position, value);
}



inline
void LaspackMatrix::add_matrix(const DenseMatrix<Real>& dm,
			       const std::vector<unsigned int>& dof_indices)
{
  add_matrix (dm, dof_indices, dof_indices);
}



inline
void LaspackMatrix::add_matrix(const DenseMatrix<Real>& dm,
			       const std::vector<unsigned int>& rows,
			       const std::vector<unsigned int>& cols)
		    
{
  assert (initialized());
  assert (dm.m() == rows.size());
  assert (dm.n() == cols.size());

  
  for (unsigned int i=0; i<rows.size(); i++)
    for (unsigned int j=0; j<cols.size(); j++)
      add(rows[i],cols[j],dm(i,j));
}



inline
Real LaspackMatrix::operator () (const unsigned int i,
				 const unsigned int j) const
{
  assert (initialized());
  assert (i < m());
  assert (j < n());
  
  using namespace Laspack;

  return Q_GetEl (const_cast<Laspack::QMatrix*>(&_QMat), i+1, j+1);
}



inline
unsigned int LaspackMatrix::pos (const unsigned int i,
				 const unsigned int j) const
{
  //std::cout << "m()=" << m() << std::endl;
  assert (i < this->m());
  assert (j < this->n());
  assert (i+1 < _row_start.size());
  assert (_row_start.back() == _csr.size());

  std::pair<const unsigned int*,
    const unsigned int*> p =
    std::equal_range (&_csr[_row_start[i]],
		      &_csr[_row_start[i+1]],
		      j);

  // Make sure the row contains the element j
  assert (p.first != p.second);

  // Make sure the values match
  assert (*p.first == j);

  // Return the position in the compressed row
  return std::distance (&_csr[_row_start[i]], p.first);
}


#endif // #ifdef HAVE_LASPACK
#endif // #ifdef __laspack_matrix_h__
