// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


// C++ includes
#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "petsc_matrix.h"
#include "dof_map.h"
#include "dense_matrix.h"
#include "petsc_vector.h"



//-----------------------------------------------------------------------
// PetscMatrix members
template <typename T>
void PetscMatrix<T>::init (const unsigned int m,
			   const unsigned int n,
			   const unsigned int m_l,
			   const unsigned int n_l,
			   const unsigned int nnz,
			   const unsigned int noz)
{
  // We allow 0x0 matrices now
  //if ((m==0) || (n==0))
  //  return;

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;

  
  int ierr     = 0;
  int m_global = static_cast<int>(m);
  int n_global = static_cast<int>(n);
  int m_local  = static_cast<int>(m_l);
  int n_local  = static_cast<int>(n_l);
  int n_nz     = static_cast<int>(nnz);
  int n_oz     = static_cast<int>(noz);
  
  // create a sequential matrix on one processor
  if (libMesh::n_processors() == 1)
    {
      libmesh_assert ((m_l == m) && (n_l == n));

      // Create matrix.  Revisit later to do preallocation and make more efficient
      ierr = MatCreateSeqAIJ (libMesh::COMM_WORLD, n_global, n_global,
			      n_nz, PETSC_NULL, &_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = MatSetFromOptions (_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  else
    {
      parallel_only();

      ierr = MatCreateMPIAIJ (libMesh::COMM_WORLD, m_local, n_local, m_global, n_global,
			      n_nz, PETSC_NULL, n_oz, PETSC_NULL, &_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = MatSetFromOptions (_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  this->zero ();
}




template <typename T>
void PetscMatrix<T>::init ()
{
  libmesh_assert (this->_dof_map != NULL);
  
  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;

  
  int proc_id = 0;

  MPI_Comm_rank (libMesh::COMM_WORLD, &proc_id);
  
  const unsigned int m   = this->_dof_map->n_dofs();
  const unsigned int n   = m;
  const unsigned int n_l = this->_dof_map->n_dofs_on_processor(proc_id); 
  const unsigned int m_l = n_l;


  const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
  const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();

  // Make sure the sparsity pattern isn't empty unless the matrix is 0x0
  libmesh_assert (n_nz.size() == n_l);
  libmesh_assert (n_oz.size() == n_l);
  
  // We allow 0x0 matrices now
  //if (m==0)
  //  return;
  
  int ierr     = 0;
  int m_global = static_cast<int>(m);
  int n_global = static_cast<int>(n);
  int m_local  = static_cast<int>(m_l);
  int n_local  = static_cast<int>(n_l);


  // create a sequential matrix on one processor
  if (libMesh::n_processors() == 1)
    {
      libmesh_assert ((m_l == m) && (n_l == n));
      if (n_nz.empty())
        ierr = MatCreateSeqAIJ (libMesh::COMM_WORLD, n_global, n_global,
			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
      else
        ierr = MatCreateSeqAIJ (libMesh::COMM_WORLD, n_global, n_global,
			        PETSC_NULL, (int*) &n_nz[0], &_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = MatSetFromOptions (_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  else
    {
      parallel_only();

      if (n_nz.empty())
        ierr = MatCreateMPIAIJ (libMesh::COMM_WORLD,
			        m_local, n_local,
			        m_global, n_global,
			        PETSC_NULL, (int*) PETSC_NULL,
			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
      else
        ierr = MatCreateMPIAIJ (libMesh::COMM_WORLD,
			        m_local, n_local,
			        m_global, n_global,
			        PETSC_NULL, (int*) &n_nz[0],
			        PETSC_NULL, (int*) &n_oz[0], &_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = MatSetFromOptions (_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  this->zero();
}



template <typename T>
void PetscMatrix<T>::zero ()
{
  libmesh_assert (this->initialized());

  semiparallel_only();
  
  int ierr=0;

  ierr = MatZeroEntries(_mat);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}



template <typename T>
void PetscMatrix<T>::clear ()
{
  int ierr=0;
  
  if ((this->initialized()) && (this->_destroy_mat_on_exit))
    {
      semiparallel_only();
  
      ierr = MatDestroy (_mat);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      this->_is_initialized = false;
    }
}



template <typename T>
Real PetscMatrix<T>::l1_norm () const
{
  libmesh_assert (this->initialized());
  
  semiparallel_only();

  int ierr=0;
  PetscReal petsc_value;
  Real value;
  
  libmesh_assert (this->closed());

  ierr = MatNorm(_mat, NORM_1, &petsc_value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
Real PetscMatrix<T>::linfty_norm () const
{
  libmesh_assert (this->initialized());
  
  semiparallel_only();

  int ierr=0;
  PetscReal petsc_value;
  Real value;
  
  libmesh_assert (this->closed());

  ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
void PetscMatrix<T>::print_matlab (const std::string name) const
{
  libmesh_assert (this->initialized());

  semiparallel_only();

  // libmesh_assert (this->closed());
  this->close();
  
  int ierr=0; 
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (libMesh::COMM_WORLD,
			    &petsc_viewer);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.  
   */
  if (name != "NULL")
    {
      ierr = PetscViewerASCIIOpen( libMesh::COMM_WORLD,
				   name.c_str(),
				   &petsc_viewer);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      ierr = PetscViewerSetFormat (petsc_viewer,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = MatView (_mat, petsc_viewer);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
    {
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (petsc_viewer);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}




template <typename T>
void PetscMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				const std::vector<unsigned int>& rows,
				const std::vector<unsigned int>& cols)
{
  libmesh_assert (this->initialized());
  
  const unsigned int m = dm.m();
  const unsigned int n = dm.n();

  libmesh_assert (rows.size() == m);
  libmesh_assert (cols.size() == n);
  
  int ierr=0;

  // These casts are required for PETSc <= 2.1.5
  ierr = MatSetValues(_mat,
		      m, (int*) &rows[0],
		      n, (int*) &cols[0],
		      (PetscScalar*) &dm.get_values()[0],
		      ADD_VALUES);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}





template <typename T>
void PetscMatrix<T>::_get_submatrix(SparseMatrix<T>& submatrix,
				    const std::vector<unsigned int> &rows,
				    const std::vector<unsigned int> &cols,
				    const bool reuse_submatrix) const
{
  // Can only extract submatrices from closed matrices
  this->close();
  
  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T>* petsc_submatrix = libmesh_cast_ptr<PetscMatrix<T>*>(&submatrix);

  // Construct row and column index sets.
  int ierr=0;
  IS isrow, iscol;

  ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			 rows.size(),
			 (int*) &rows[0],
			 &isrow); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			 cols.size(),
			 (int*) &cols[0],
			 &iscol); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Extract submatrix
  ierr = MatGetSubMatrix(_mat,
			 isrow,
			 iscol,
			 PETSC_DECIDE,
			 (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
			 &(petsc_submatrix->_mat));  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Specify that the new submatrix is initialized and close it.
  petsc_submatrix->_is_initialized = true;
  petsc_submatrix->close();

  // Clean up PETSc data structures
  ierr = ISDestroy(isrow); CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = ISDestroy(iscol); CHKERRABORT(libMesh::COMM_WORLD,ierr);
}



template <typename T>
void PetscMatrix<T>::get_diagonal (NumericVector<T>& dest) const
{
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T>& petsc_dest = libmesh_cast_ref<PetscVector<T>&>(dest);

  // Call PETSc function.

#if PETSC_VERSION_LESS_THAN(2,3,1)

  std::cout << "This method has been developed with PETSc 2.3.1.  "
	    << "No one has made it backwards compatible with older "
	    << "versions of PETSc so far; however, it might work "
	    << "without any change with some older version." << std::endl;
  libmesh_error();

#else

  // Needs a const_cast since PETSc does not work with const.
  int ierr =
    MatGetDiagonal(const_cast<PetscMatrix<T>*>(this)->mat(),petsc_dest.vec()); CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif

}



template <typename T>
void PetscMatrix<T>::get_transpose (SparseMatrix<T>& dest) const
{
  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T>& petsc_dest = libmesh_cast_ref<PetscMatrix<T>&>(dest);

  int ierr;
#if PETSC_VERSION_LESS_THAN(3,0,0)
  if (&petsc_dest == this)
    ierr = MatTranspose(_mat,PETSC_NULL);
  else
    ierr = MatTranspose(_mat,&petsc_dest._mat);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
  // FIXME - we can probably use MAT_REUSE_MATRIX in more situations
  if (&petsc_dest == this)
    ierr = MatTranspose(_mat,MAT_REUSE_MATRIX,&petsc_dest._mat);
  else
    ierr = MatTranspose(_mat,MAT_INITIAL_MATRIX,&petsc_dest._mat);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
}



//------------------------------------------------------------------
// Explicit instantiations
template class PetscMatrix<Number>;


#endif // #ifdef LIBMESH_HAVE_PETSC
