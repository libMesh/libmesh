//    $Id: petsc_matrix.h,v 1.9 2003-02-07 22:18:53 benkirk Exp $

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



#ifndef __petsc_matrix_h__
#define __petsc_matrix_h__

#include "mesh_common.h"

#ifdef HAVE_PETSC

// TODO:[BSK} This seems necessary to use petsc on IBM Power3 at NERSC, but only there?  
#include <cmath>


// C++ includes

// Local includes
#include "sparse_matrix.h"



/**
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
# ifndef USE_COMPLEX_NUMBERS

namespace Petsc {
extern "C" {
#include <petscmat.h>
}
} 
using namespace Petsc;

#else

#include <petscmat.h>

#endif



// Forward declarations
class PetscVector;
class PetscInterface;



/**
 * Petsc matrix. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 *
 * @author Benjamin S. Kirk, 2002
 */

class PetscMatrix : public SparseMatrix
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
  PetscMatrix ();

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~PetscMatrix ();

  /**
   * Initialize a Petsc matrix that is of global
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
   * Call the Petsc assemble routines.
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
   * Petsc matrix.  This is useful
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
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
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
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real linfty_norm () const;

  /**
   * see if Petsc matrix has been closed
   * and fully assembled yet
   */
  bool closed() const;

  /**
   * Print the contents of the matrix to the screen.
   */
  void print() const;
  
  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  void print_matlab(const std::string name="NULL") const;

  
 private:

  
  /**
   * Petsc matrix datatype to store values
   */				      
  Mat mat;

  /**
   * Make other Petsc datatypes friends
   */
  friend class PetscInterface;
};




//-----------------------------------------------------------------------
// PetscMatrix inline members
inline
PetscMatrix::PetscMatrix()
{};




inline
PetscMatrix::~PetscMatrix()
{
  clear();
};




inline
void PetscMatrix::close () const
{
  if (closed())
    return;
  
  int ierr=0;
 
  ierr = MatAssemblyBegin (mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd   (mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  return;
};




inline
unsigned int PetscMatrix::m () const
{
  assert (initialized());
  
  int petsc_m=0, petsc_n=0, ierr=0;

  ierr = MatGetSize (mat, &petsc_m, &petsc_n);

  return static_cast<unsigned int>(petsc_m);
};




inline
unsigned int PetscMatrix::n () const
{
  assert (initialized());
  
  int petsc_m=0, petsc_n=0, ierr=0;

  ierr = MatGetSize (mat, &petsc_m, &petsc_n);

  return static_cast<unsigned int>(petsc_n);
};




inline
unsigned int PetscMatrix::row_start () const
{
  assert (initialized());
  
  int start=0, stop=0, ierr=0;

  ierr = MatGetOwnershipRange(mat, &start, &stop); CHKERRQ(ierr);

  return start;
};




inline
unsigned int PetscMatrix::row_stop () const
{
  assert (initialized());
  
  int start=0, stop=0, ierr=0;

  ierr = MatGetOwnershipRange(mat, &start, &stop); CHKERRQ(ierr);

  return stop;
};




inline
void PetscMatrix::set (const unsigned int i,
		       const unsigned int j,
		       const Complex value)
{  
  assert (initialized());
  
  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(mat, 1, &i_val, 1, &j_val,
		      &petsc_value, INSERT_VALUES); CHKERRQ(ierr);

  return;
};




inline
void PetscMatrix::add (const unsigned int i,
		       const unsigned int j,
		       const Complex value)
{
  assert (initialized());
  
  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(mat, 1, &i_val, 1, &j_val,
		      &petsc_value, ADD_VALUES); CHKERRQ(ierr);
  
  return;
};



inline
void PetscMatrix::add_matrix(const ComplexDenseMatrix& dm,
			     const std::vector<unsigned int>& dof_indices)
{
  add_matrix (dm, dof_indices, dof_indices);
};


inline
void PetscMatrix::add_matrix(const ComplexDenseMatrix& dm,
			     const std::vector<unsigned int>& rows,
			     const std::vector<unsigned int>& cols)
		    
{
  assert (initialized());
  
  const unsigned int m = dm.m();
  const unsigned int n = dm.n();

  assert (rows.size() == m);
  assert (cols.size() == n);
  
  int ierr=0;

  // make this static to the function to aviod repeated allocations
  static std::vector<PetscScalar> values;

  values.resize (m*n);

  // notice values is row-major by default in Petsc
  for (unsigned int i=0; i<m; i++)
    for (unsigned int j=0; j<n; j++)
      values[(i)*(n) + (j)] = static_cast<PetscScalar>(dm(i,j)); 

  ierr = MatSetValues(mat,
		      m, (int*) &rows[0],
		      n, (int*) &cols[0],
		      &values[0],
		      ADD_VALUES);   CHKERRQ(ierr);

  return;
};



inline
Complex PetscMatrix::operator () (const unsigned int i,
				  const unsigned int j) const
{
  assert (initialized());
  
  PetscScalar *petsc_row;
  Complex value=0.;
  bool found=false;
  int ierr=0, ncols=0, *petsc_cols,
    i_val=static_cast<int>(i),
    j_val=static_cast<int>(j);
  

  // the matrix needs to be closed for this to work
  close();

  
  ierr = MatGetRow(mat, i_val, &ncols, &petsc_cols, &petsc_row); CHKERRQ(ierr);


  for (int entry=0; entry<ncols; entry++)
    if (petsc_cols[entry] == j_val)
      {
	found = true;
	  
	value = static_cast<Complex>(petsc_row[entry]);
	  
	ierr = MatRestoreRow(mat, i_val,
			     &ncols, &petsc_cols, &petsc_row); CHKERRQ(ierr);
	  
	return value;
      };
  
  // Otherwise the entry is not in the sparse matrix,
  // i.e. it is 0.
  
  return 0.;
};




inline
bool PetscMatrix::closed() const
{
  assert (initialized());
  
  int ierr=0;
  PetscTruth assembled;

  ierr = MatAssembled(mat, &assembled); CHKERRQ(ierr);

  if (assembled == PETSC_TRUE)
    return true;

  return false;
};




inline
void PetscMatrix::print() const
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
	std::cout << std::setw(8) << (*this)(i,j).real() << " ";
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
/*---------------------------   petsc_matrix.h     -------------------------*/
#endif
