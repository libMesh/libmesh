//    $Id: sparse_matrix.h,v 1.2 2003-02-10 03:55:51 benkirk Exp $

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



#ifndef __sparse_matrix_h__
#define __sparse_matrix_h__


// C++ includes
#include <vector>
#include <set>

// Local includes
#include "mesh_common.h"
#include "dof_map.h"
#include "reference_counted_object.h"


// forward declarations
class SparseMatrix;




/**
 * Generic sparse matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to sparse matrices
 * from various different solver packages in
 * different formats.
 *
 * @author Benjamin S. Kirk, 2003
 */

class SparseMatrix : public ReferenceCountedObject<SparseMatrix>
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
  SparseMatrix ();

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  virtual ~SparseMatrix ();

  /**
   * @returns true if the matrix has been initialized,
   * false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; };

  /**
   * Get a pointer to the \p DofMap to use.
   */
  void attach_dof_map (const DofMap& dof_map)
  { _dof_map = &dof_map; };
  

  /**
   * Updates the matrix sparsity pattern.  This method is
   * included because some sparse matrix storage schemes
   * (e.g. LASPACK) need it.  If your \p SparseMatrix
   * implementation does not need this data simply do
   * not overload this method.
   */
  virtual void update_sparsity_pattern (std::vector<std::set<unsigned int> >&) {};

  
  /**
   * Initialize a Sparse matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   */
  virtual void init (const unsigned int m,
		     const unsigned int n,
		     const unsigned int m_l,
		     const unsigned int n_l,
		     const unsigned int nnz=30,
		     const unsigned int noz=10) = 0;

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */   
  virtual void init () = 0;
  
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  virtual void clear () = 0;

  /**
   * Set all entries to 0.
   */
  virtual void zero () = 0;
  
  /**
   * Call the Sparse assemble routines.
   * sends necessary messages to other
   * processors
   */
  virtual void close () const = 0;
  
  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  virtual unsigned int m () const = 0;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  virtual unsigned int n () const = 0;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  virtual unsigned int row_start () const = 0;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  virtual unsigned int row_stop () const = 0;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  virtual void set (const unsigned int i,
		    const unsigned int j,
		    const Complex value) = 0;
    
  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  virtual void add (const unsigned int i,
		    const unsigned int j,
		    const Complex value) = 0;

  /**
   * Add the full matrix to the
   * Sparse matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
    
  virtual void add_matrix (const ComplexDenseMatrix &dm,
			   const std::vector<unsigned int> &rows,
			   const std::vector<unsigned int> &cols) = 0;
  
  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_matrix (const ComplexDenseMatrix &dm,
			   const std::vector<unsigned int> &dof_indices) = 0;
      
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
  virtual Complex operator () (const unsigned int i,
			       const unsigned int j) const = 0;

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
  virtual Real l1_norm () const = 0;

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
  virtual Real linfty_norm () const = 0;

  /**
   * see if Sparse matrix has been closed
   * and fully assembled yet
   */
  virtual bool closed() const = 0;

  /**
   * Print the contents of the matrix to the screen.
   */
  virtual void print() const;
  
  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string name="NULL") const
    {
      std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
      std::cerr << "ERROR writing MATLAB file " << name << std::endl;
      error();
    };
  
protected:
  
  /**
   * The \p DofMap object associated with this object.
   */
  DofMap const *_dof_map;
  
  /**
   * Flag indicating whether or not the matrix
   * has been initialized.
   */
  bool _is_initialized;
};



//-----------------------------------------------------------------------
// SparseMatrix inline members
inline
SparseMatrix::SparseMatrix () :
  _dof_map(NULL),
  _is_initialized(false)
{};



inline
SparseMatrix::~SparseMatrix ()
{};



inline
void SparseMatrix::print() const
{
  assert (initialized());

#ifndef USE_COMPLEX_NUMBERS

  for (unsigned int i=0; i<m(); i++)
    {
      for (unsigned int j=0; j<n(); j++)
	std::cout << std::setw(8) << (*this)(i,j) << " ";
      std::cout << std::endl;
    }

#else
  // std::complex<>::operator<<() is defined, but use this form

  std::cout << "Real part:" << std::endl;
  for (unsigned int i=0; i<m(); i++)
    {
      for (unsigned int j=0; j<n(); j++)
	std::cout << std::setw(14) << (*this)(i,j).real() << " ";
      std::cout << std::endl;
    }

  std::cout << std::endl << "Imaginary part:" << std::endl;
  for (unsigned int i=0; i<m(); i++)
    {
      for (unsigned int j=0; j<n(); j++)
	std::cout << std::setw(8) << (*this)(i,j).imag() << " ";
      std::cout << std::endl;
    }

#endif

};



#endif
