//    $Id: petsc_matrix.h,v 1.26 2003-08-04 12:43:06 ddreyer Exp $

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
//#include <cmath>


// C++ includes

// Local includes
#include "sparse_matrix.h"
#include "dense_matrix.h"



/*
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
#ifndef USE_COMPLEX_NUMBERS

extern "C" {
#include <petscmat.h>
} 

#else

#include <petscmat.h>

#endif



// Forward declarations
template <typename T> class PetscVector;
template <typename T> class PetscInterface;



/**
 * Petsc matrix. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 *
 * @author Benjamin S. Kirk, 2002
 */

template <typename T>
class PetscMatrix : public SparseMatrix<T>
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
  void init ();
  
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  void clear ();

  /**
   * Set all entries to 0. This method retains 
   * sparsity structure.
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
	    const T value);
    
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
	    const T value);

  /**
   * Add the full matrix to the
   * Petsc matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
    
  void add_matrix (const DenseMatrix<T> &dm,
		   const std::vector<unsigned int> &rows,
		   const std::vector<unsigned int> &cols);	     

  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  void add_matrix (const DenseMatrix<T> &dm,
		   const std::vector<unsigned int> &dof_indices);	     
      
  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this: 
   * \f$\texttt{this} = a*X + \texttt{this} \f$.
   * Use this with caution, the sparse matrices need to have the
   * same nonzero pattern, otherwise \p PETSc will crash! 
   * It is advisable to not only allocate appropriate memory with 
   * \p init() , but also explicitly zero the terms of \p this
   * whenever you add a non-zero value to \p X.  Note: \p X will 
   * be closed, if not already done, before performing any work.
   */
  void add (const T a, SparseMatrix<T> &X);
    
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
  T operator () (const unsigned int i,
		 const unsigned int j) const;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all 
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real l1_norm () const;

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_infty=max_{all rows i}\sum_{all 
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_infty \leq |M|_infty |v|_infty\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real linfty_norm () const;

  /**
   * see if Petsc matrix has been closed
   * and fully assembled yet
   */
  bool closed() const;
  
  /**
   * Print the contents of the matrix to the screen
   * with the PETSc viewer.
   */
  void print_personal() const;

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
  friend class PetscInterface<T>;
  friend class PetscVector<T>;
};




//-----------------------------------------------------------------------
// PetscMatrix inline members
template <typename T>
inline
PetscMatrix<T>::PetscMatrix()
{}



template <typename T>
inline
PetscMatrix<T>::~PetscMatrix()
{
  clear();
}



template <typename T>
inline
void PetscMatrix<T>::close () const
{
  if (closed())
    return;
  
  int ierr=0;
 
  ierr = MatAssemblyBegin (mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd   (mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  return;
}



template <typename T>
inline
unsigned int PetscMatrix<T>::m () const
{
  assert (this->initialized());
  
  int petsc_m=0, petsc_n=0, ierr=0;

  ierr = MatGetSize (mat, &petsc_m, &petsc_n);

  return static_cast<unsigned int>(petsc_m);
}



template <typename T>
inline
unsigned int PetscMatrix<T>::n () const
{
  assert (this->initialized());
  
  int petsc_m=0, petsc_n=0, ierr=0;

  ierr = MatGetSize (mat, &petsc_m, &petsc_n);

  return static_cast<unsigned int>(petsc_n);
}



template <typename T>
inline
unsigned int PetscMatrix<T>::row_start () const
{
  assert (this->initialized());
  
  int start=0, stop=0, ierr=0;

  ierr = MatGetOwnershipRange(mat, &start, &stop); CHKERRQ(ierr);

  return static_cast<unsigned int>(start);
}



template <typename T>
inline
unsigned int PetscMatrix<T>::row_stop () const
{
  assert (this->initialized());
  
  int start=0, stop=0, ierr=0;

  ierr = MatGetOwnershipRange(mat, &start, &stop); CHKERRQ(ierr);

  return static_cast<unsigned int>(stop);
}



template <typename T>
inline
void PetscMatrix<T>::set (const unsigned int i,
			  const unsigned int j,
			  const T value)
{  
  assert (this->initialized());
  
  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(mat, 1, &i_val, 1, &j_val,
		      &petsc_value, INSERT_VALUES); CHKERRQ(ierr);

  return;
}



template <typename T>
inline
void PetscMatrix<T>::add (const unsigned int i,
			  const unsigned int j,
			  const T value)
{
  assert (this->initialized());
  
  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(mat, 1, &i_val, 1, &j_val,
		      &petsc_value, ADD_VALUES); CHKERRQ(ierr);
  
  return;
}



template <typename T>
inline
void PetscMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				const std::vector<unsigned int>& dof_indices)
{
  add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
inline
void PetscMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				const std::vector<unsigned int>& rows,
				const std::vector<unsigned int>& cols)
{
  assert (this->initialized());
  
  const unsigned int m = dm.m();
  const unsigned int n = dm.n();

  assert (rows.size() == m);
  assert (cols.size() == n);
  
  int ierr=0;

//   // make this static to the function to aviod repeated allocations
//   static std::vector<PetscScalar> values;

//   values.resize (m*n);

//   // notice values is row-major by default in Petsc
//   for (unsigned int i=0; i<m; i++)
//     for (unsigned int j=0; j<n; j++)
//       values[(i)*(n) + (j)] = static_cast<PetscScalar>(dm(i,j)); 

  ierr = MatSetValues(mat,
		      m, (int*) &rows[0],
		      n, (int*) &cols[0],
		      (PetscScalar*) &dm.get_values()[0],
		      ADD_VALUES);   CHKERRQ(ierr);

  return;
}




template <typename T>
inline
void PetscMatrix<T>::add (const T a_in, SparseMatrix<T> &X_in)
{
  assert (this->initialized());

  // sanity check. but this cannot avoid 
  // crash due to incompatible sparsity structure...
  assert (this->m() == X_in.m());
  assert (this->n() == X_in.n());

  PetscScalar     a = static_cast<PetscScalar>      (a_in);
  PetscMatrix<T>& X = dynamic_cast<PetscMatrix<T>&> (X_in);
  int ierr=0;

  // the matrix from which we copy the values has to be assembled/closed
  X.close ();

  ierr = MatAXPY(&a,  X.mat, mat, SAME_NONZERO_PATTERN);   CHKERRQ(ierr);
  return;
}




template <typename T>
inline
T PetscMatrix<T>::operator () (const unsigned int i,
			       const unsigned int j) const
{
  assert (this->initialized());
  
  PetscScalar *petsc_row;
  T value=0.;
  bool found=false;
  int ierr=0, ncols=0, *petsc_cols,
    i_val=static_cast<int>(i),
    j_val=static_cast<int>(j);
  

  // the matrix needs to be closed for this to work
  close();

  
  ierr = MatGetRow(mat, i_val, &ncols, &petsc_cols, &petsc_row); CHKERRQ(ierr);


  //TODO:[BSK] A binary search on petsc_cols would be faster!
  for (int entry=0; entry<ncols; entry++)
    if (petsc_cols[entry] == j_val)
      {
	found = true;
	  
	value = static_cast<T>(petsc_row[entry]);
	  
	ierr = MatRestoreRow(mat, i_val,
			     &ncols, &petsc_cols, &petsc_row); CHKERRQ(ierr);
	  
	return value;
      }
  
  // Otherwise the entry is not in the sparse matrix,
  // i.e. it is 0.
  
  return 0.;
}




template <typename T>
inline
bool PetscMatrix<T>::closed() const
{
  assert (this->initialized());
  
  int ierr=0;
  PetscTruth assembled;

  ierr = MatAssembled(mat, &assembled); CHKERRQ(ierr);

  return (assembled == PETSC_TRUE);
}




template <typename T>
inline
void PetscMatrix<T>::print_personal() const
{
  assert (this->initialized());
  
  int ierr=0;

  ierr = MatView(mat, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
}

#endif // #ifdef HAVE_PETSC
#endif // #ifdef __petsc_matrix_h__
