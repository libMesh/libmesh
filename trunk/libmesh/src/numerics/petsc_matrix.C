// $Id: petsc_matrix.C,v 1.18 2003-09-06 02:24:00 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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

// Local includes
#include "petsc_matrix.h"
#include "petsc_vector.h"

#ifdef HAVE_PETSC

#include "dof_map.h"




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
  if ((m==0) || (n==0))
    return;

  {
    // Clear initialized matrices
    if (this->initialized())
      this->clear();

    this->_is_initialized = true;
  }

  
  int ierr=0;
  int m_global=static_cast<int>(m);
  int n_global=static_cast<int>(n);
  int m_local=static_cast<int>(m_l);
  int n_local=static_cast<int>(n_l);
  int n_nz=static_cast<int>(nnz);
  int n_oz=static_cast<int>(noz);
  
  // create a sequential matrix on one processor
  if ((m_l == m) && (n_l == n))
    {
      // Create matrix.  Revisit later to do pReallocation and make more efficient
      ierr = MatCreateSeqAIJ (PETSC_COMM_WORLD, n_global, n_global,
			      n_nz, PETSC_NULL, &mat);                 CHKERRQ(ierr);
  
      ierr = MatSetFromOptions (mat);                                  CHKERRQ(ierr);
    }

  else
    {
      ierr = MatCreateMPIAIJ (PETSC_COMM_WORLD, m_local, n_local, m_global, n_global,
			      n_nz, PETSC_NULL, n_oz, PETSC_NULL, &mat); CHKERRQ(ierr);
  
      ierr = MatSetFromOptions (mat);                                  CHKERRQ(ierr);
    }

  zero ();
}




template <typename T>
void PetscMatrix<T>::init ()
{
  assert (this->_dof_map != NULL);
  
  {
    // Clear initialized matrices
    if (this->initialized())
      this->clear();

    this->_is_initialized = true;
  }

  
  int proc_id = 0;

  MPI_Comm_rank (PETSC_COMM_WORLD, &proc_id);
  
  const unsigned int m   = this->_dof_map->n_dofs();
  const unsigned int n   = m;
  const unsigned int n_l = this->_dof_map->n_dofs_on_processor(proc_id); 
  const unsigned int m_l = n_l;


  const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
  const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();

  // Make sure the sparsity pattern isn't empty
  assert (n_nz.size() == n_l);
  assert (n_oz.size() == n_l);
  
  if (m==0)
    return;
  
  int ierr=0;
  int m_global=static_cast<int>(m);
  int n_global=static_cast<int>(n);
  int m_local=static_cast<int>(m_l);
  int n_local=static_cast<int>(n_l);


  // create a sequential matrix on one processor
  if ((m_l == m) && (n_l == n))
    {
      ierr = MatCreateSeqAIJ (PETSC_COMM_WORLD, n_global, n_global,
			      PETSC_NULL, (int*) &n_nz[0], &mat);      CHKERRQ(ierr);
  
      ierr = MatSetFromOptions (mat);                                  CHKERRQ(ierr);
    }

  else
    {
      ierr = MatCreateMPIAIJ (PETSC_COMM_WORLD,
			      m_local, n_local,
			      m_global, n_global,
			      PETSC_NULL, (int*) &n_nz[0],
			      PETSC_NULL, (int*) &n_oz[0], &mat);      CHKERRQ(ierr);
  
      ierr = MatSetFromOptions (mat);                                  CHKERRQ(ierr);
    }

  zero();
}



template <typename T>
void PetscMatrix<T>::zero ()
{
  assert (this->initialized());
  
  int ierr=0;

  ierr = MatZeroEntries(mat); CHKERRQ(ierr);
}



template <typename T>
void PetscMatrix<T>::clear ()
{
  int ierr=0;
  
  if (this->initialized())
    {
      ierr = MatDestroy (mat); CHKERRQ(ierr);
      
      this->_is_initialized = false;
    }
}



template <typename T>
Real PetscMatrix<T>::l1_norm () const
{
  assert (this->initialized());
  
  int ierr=0;
  double petsc_value;
  Real value;
  
  assert (closed());

  ierr = MatNorm(mat, NORM_1, &petsc_value); CHKERRQ(ierr);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
Real PetscMatrix<T>::linfty_norm () const
{
  assert (this->initialized());
  
  int ierr=0;
  double petsc_value;
  Real value;
  
  assert (closed());

  ierr = MatNorm(mat, NORM_INFINITY, &petsc_value); CHKERRQ(ierr);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
void PetscMatrix<T>::print_matlab (const std::string name) const
{
  assert (this->initialized());
  assert (closed());
  
  int ierr=0; 
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (PETSC_COMM_WORLD,
			    &petsc_viewer);                    CHKERRQ(ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.  
   */
  if (name != "NULL")
    {
      ierr = PetscViewerASCIIOpen( PETSC_COMM_WORLD,
				   name.c_str(),
				   &petsc_viewer);             CHKERRQ(ierr);
      
      ierr = PetscViewerSetFormat (petsc_viewer,
				   PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
  
      ierr = MatView (mat, petsc_viewer);                      CHKERRQ(ierr);
    }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
    {
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
				   PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
  
      ierr = MatView (mat, PETSC_VIEWER_STDOUT_WORLD);         CHKERRQ(ierr);
    }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (petsc_viewer);                    CHKERRQ(ierr);
}



//------------------------------------------------------------------
// Explicit instantiations
template class PetscMatrix<Number>;


#endif // #ifdef HAVE_PETSC
