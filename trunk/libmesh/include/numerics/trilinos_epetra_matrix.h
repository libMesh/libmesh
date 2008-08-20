// $Id: petsc_matrix.h 2789 2008-04-13 02:24:40Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __trilinos_epetra_matrix_h__
#define __trilinos_epetra_matrix_h__

#include "libmesh_common.h"

#ifdef HAVE_TRILINOS

// Trilinos includes
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

// C++ includes
#include <algorithm>

// Local includes
#include "sparse_matrix.h"

// Forward Declarations
template <typename T> class DenseMatrix;



/**
 * Epetra matrix. Provides a nice interface to the
 * Epetra data structures for parallel,
 * sparse matrices.
 *
 * @author Benjamin S. Kirk, 2008
 */

template <typename T>
class EpetraMatrix : public SparseMatrix<T>
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
  EpetraMatrix ();

//   /**
//    * Constructor.  Creates a EpetraMatrix assuming you already
//    * have a valid Mat object.  In this case, m is NOT destroyed
//    * by the EpetraMatrix destructor when this object goes out of scope.
//    * This allows ownership of m to remain with the original creator,
//    * and to simply provide additional functionality with the EpetraMatrix.
//    */
//   EpetraMatrix (Mat m);
  
  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~EpetraMatrix ();
 
  /**
   * The \p EpetraMatrix needs the full sparsity pattern.
   */ 
  bool need_full_sparsity_pattern () const 
  { return true; }

  /**
   * Updates the matrix sparsity pattern.  This will
   * tell the underlying matrix storage scheme how
   * to map the \f$ (i,j) \f$ elements.
   */
  void update_sparsity_pattern (const SparsityPattern::Graph &);

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
   * with the PETSc viewer.  This function only allows
   * printing to standard out, this is because we have
   * limited ourselves to one PETSc implementation for
   * writing.
   */
  void print_personal(std::ostream& os=std::cout) const;

  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  void print_matlab(const std::string name="NULL") const;

  /**
   * Swaps the raw PETSc matrix context pointers.
   */
  void swap (EpetraMatrix<T> &);

//   /**
//    * Returns the raw PETSc matrix context pointer.  Note this is generally
//    * not required in user-level code. Just don't do anything crazy like
//    * calling MatDestroy()!
//    */
//   Mat mat () { libmesh_assert (_mat != NULL); return _mat; }

  
protected:

//   /**
//    * This function either creates or re-initializes
//    * a matrix called "submatrix" which is defined
//    * by the row and column indices given in the "rows" and "cols" entries.
//    * This function is implemented in terms of the MatGetSubMatrix()
//    * routine of PETSc.  The boolean reuse_submatrix parameter determines
//    * whether or not PETSc will treat "submatrix" as one which has already
//    * been used (had memory allocated) or as a new matrix.
//    */
//   virtual void _get_submatrix(SparseMatrix<T>& submatrix,
// 			      const std::vector<unsigned int>& rows,
// 			      const std::vector<unsigned int>& cols,
// 			      const bool reuse_submatrix) const;

private:
  
  /**
   * Actual Epetra datatype
   * to hold matrix entries
   */
  AutoPtr<Epetra_FECrsMatrix> _mat;

  /**
   * Holds the distributed Map
   */
  AutoPtr<Epetra_Map> _map;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Mat object. 
   */
  bool _destroy_mat_on_exit;
};




//-----------------------------------------------------------------------
// EpetraMatrix inline members
template <typename T>
inline
EpetraMatrix<T>::EpetraMatrix()
  : _destroy_mat_on_exit(true)
{}




// template <typename T>
// inline
// EpetraMatrix<T>::EpetraMatrix(Mat m)
//   : _destroy_mat_on_exit(false)
// {
//   this->_mat = m;
//   this->_is_initialized = true;
// }




template <typename T>
inline
EpetraMatrix<T>::~EpetraMatrix()
{
  this->clear();
}



template <typename T>
inline
void EpetraMatrix<T>::close () const
{
  libmesh_assert (_mat.get() != NULL);
  
  _mat->FillComplete();
}



template <typename T>
inline
unsigned int EpetraMatrix<T>::m () const
{
  libmesh_assert (this->initialized());
  
  return static_cast<unsigned int>(_mat->NumGlobalRows());
}



template <typename T>
inline
unsigned int EpetraMatrix<T>::n () const
{
  libmesh_assert (this->initialized());
  
  return static_cast<unsigned int>(_mat->NumGlobalCols());
}



template <typename T>
inline
unsigned int EpetraMatrix<T>::row_start () const
{
  libmesh_assert (this->initialized());
  libmesh_assert (_map.get() != NULL);

  return static_cast<unsigned int>(_map->MinMyGID());
}



template <typename T>
inline
unsigned int EpetraMatrix<T>::row_stop () const
{
  libmesh_assert (this->initialized());
  libmesh_assert (_map.get() != NULL);

  return static_cast<unsigned int>(_map->MaxMyGID())+1;
}



template <typename T>
inline
void EpetraMatrix<T>::set (const unsigned int i,
			   const unsigned int j,
			   const T value)
{  
  libmesh_assert (this->initialized());

  int 
    epetra_i = static_cast<int>(i),
    epetra_j = static_cast<int>(j);

  T epetra_value = value;

  if (_mat->Filled())    
    _mat->ReplaceGlobalValues (epetra_i, 1, &epetra_value, &epetra_j);
  else    
    _mat->InsertGlobalValues (epetra_i, 1, &epetra_value, &epetra_j);
}



template <typename T>
inline
void EpetraMatrix<T>::add (const unsigned int i,
			   const unsigned int j,
			   const T value)
{
  libmesh_assert (this->initialized());
  
  int 
    epetra_i = static_cast<int>(i),
    epetra_j = static_cast<int>(j);

  T epetra_value = value;
  
  _mat->SumIntoGlobalValues (epetra_i, 1, &epetra_value, &epetra_j);
}



template <typename T>
inline
void EpetraMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				 const std::vector<unsigned int>& dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
inline
void EpetraMatrix<T>::add (const T a_in, SparseMatrix<T> &X_in)
{
  libmesh_assert (this->initialized());

  // sanity check. but this cannot avoid 
  // crash due to incompatible sparsity structure...
  libmesh_assert (this->m() == X_in.m());
  libmesh_assert (this->n() == X_in.n());

  LIBMESH_THROW(libMesh::NotImplemented());

//   PetscScalar     a = static_cast<PetscScalar>      (a_in);
//   EpetraMatrix<T>* X = dynamic_cast<EpetraMatrix<T>*> (&X_in);

//   libmesh_assert (X != NULL);
  
//   int ierr=0;

//   // the matrix from which we copy the values has to be assembled/closed
//   X->close ();

// // 2.2.x & earlier style
// #if PETSC_VERSION_LESS_THAN(2,3,0)  
  
//   ierr = MatAXPY(&a,  X->_mat, _mat, SAME_NONZERO_PATTERN);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
// // 2.3.x & newer
// #else
  
//   ierr = MatAXPY(_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
// #endif
}




template <typename T>
inline
T EpetraMatrix<T>::operator () (const unsigned int i,
				const unsigned int j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->_mat.get() != NULL);
  libmesh_assert (this->_mat->MyGlobalRow(i));
  libmesh_assert (i >= this->row_start());
  libmesh_assert (i < this->row_stop());


  int row_length, *row_indices;
  double *values;

  _mat->ExtractMyRowView (i-this->row_start(),
			  row_length,
			  values,
			  row_indices);

  //std::cout << "row_length=" << row_length << std::endl;

  int *index = std::lower_bound (row_indices, row_indices+row_length, j);

  libmesh_assert (*index < row_length);
  libmesh_assert (row_indices[*index] == j);

  //std::cout << "val=" << values[*index] << std::endl;
  
  return values[*index];
}




template <typename T>
inline
bool EpetraMatrix<T>::closed() const
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->_mat.get() != NULL);
  
  return this->_mat->Filled();
}


template <typename T>
inline
void EpetraMatrix<T>::swap(EpetraMatrix<T> &m) 
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   std::swap(_mat, m._mat);
//   std::swap(_destroy_mat_on_exit, m._destroy_mat_on_exit);
}





template <typename T>
inline
void EpetraMatrix<T>::print_personal(std::ostream& os) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (_mat.get() != NULL);

  os << *_mat;
}

#endif // #ifdef HAVE_TRILINOS
#endif // #ifdef __trilinos_epetra_matrix_h__
