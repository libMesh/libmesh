//    $Id: petsc_matrix.C,v 1.3 2003-01-20 17:06:46 jwpeterson Exp $

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



#include "mesh_common.h"

#ifdef HAVE_PETSC

#include "petsc_matrix.h"
#include "petsc_vector.h"





PetscMatrix::PetscMatrix () :
  is_initialized(false)
{};




PetscMatrix::~PetscMatrix ()
{
  clear ();
};




void PetscMatrix::init (const unsigned int m,
			const unsigned int n,
			const unsigned int m_l,
			const unsigned int n_l,
			const unsigned int nnz,
			const unsigned int noz)
{
  if ((m==0) || (n==0))
    return;

  {
    if (initialized())
      {
	std::cerr << "ERROR: Matrix already initialized!"
		  << std::endl;
	
	error();
      }

    is_initialized = true;
  };

  
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
      // Create matrix.  Revisit later to do preallocation and make more efficient
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
};





void PetscMatrix::init (const DofMap& dof_map)
{
  {
    if (initialized())
      {
	std::cerr << "ERROR: Matrix already initialized!"
		  << std::endl;
	
	error();
      }

    is_initialized = true;
  };

  
  int proc_id = 0;

  MPI_Comm_rank (PETSC_COMM_WORLD, &proc_id);
  
  const unsigned int m   = dof_map.n_dofs();
  const unsigned int n   = m;
  const unsigned int n_l = dof_map.n_dofs_on_processor(proc_id); 
  const unsigned int m_l = n_l;


  const std::vector<unsigned int>& n_nz = dof_map.get_n_nz();
  const std::vector<unsigned int>& n_oz = dof_map.get_n_oz();

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
};





void PetscMatrix::init (const DofMap& dof_map,
			PetscMatrix& parent_matrix)
{
  {
    if (initialized())
      {
	std::cerr << "ERROR: Matrix already initialized!"
		  << std::endl;
	
	error();
      }

    is_initialized = true;
  };

  
  int proc_id = 0;

  MPI_Comm_rank (PETSC_COMM_WORLD, &proc_id);
  
  const unsigned int m   = dof_map.n_dofs();
  const unsigned int m_l = dof_map.n_dofs_on_processor(proc_id); 
  const unsigned int n_l = m_l;


  const std::vector<unsigned int>& n_nz = dof_map.get_n_nz();
  const std::vector<unsigned int>& n_oz = dof_map.get_n_oz();

  // Make sure the sparsity pattern isn't empty
  assert (n_nz.size() == n_l);
  assert (n_oz.size() == n_l);

  // Make sure our matrix will fit inside the parent matrix
  {
    assert (n_nz.size() <= (parent_matrix.row_stop() -
			    parent_matrix.row_start()) );
    
    for (unsigned int row=0; row<n_nz.size(); row++)
      assert (n_nz[row] <= (parent_matrix.row_stop() -
			    parent_matrix.row_start()) );
  };

  
  if (m==0)
    return;
  
  int ierr=0;
  int n_local=static_cast<int>(n_l);


  // create the index sets with which we will extract storage
  // from the parent matrix
  IS is_rows, is_cols;
  {
    std::vector<int> local_rows(m_l, 0);

    unsigned int next_row = parent_matrix.row_start();
  
    for (unsigned int row=0; row<m_l; row++)
      local_rows[row] = next_row++;

    ierr = ISCreateGeneral(PETSC_COMM_WORLD, m_l, &local_rows[0], &is_rows);
                                           CHKERRQ(ierr);

    ierr = ISAllGather(is_rows, &is_cols); CHKERRQ(ierr);

    ierr = MatGetSubMatrix(parent_matrix.mat, is_rows, is_cols,
			   n_local, MAT_INITIAL_MATRIX, &mat);
                                           CHKERRQ(ierr);

    ierr = ISDestroy(is_rows); CHKERRQ(ierr); 
    ierr = ISDestroy(is_cols); CHKERRQ(ierr); 					   
  };

  zero();
};



void PetscMatrix::zero ()
{
  assert (initialized());
  
  int ierr=0;

  ierr = MatZeroEntries(mat); CHKERRQ(ierr);
};



void PetscMatrix::clear ()
{
  int ierr=0;
  
  if (initialized())
    {
      ierr = MatDestroy (mat); CHKERRQ(ierr);
      
      is_initialized = false;
    };
};



real PetscMatrix::l1_norm () const
{
  assert (initialized());
  
  int ierr=0;
  double petsc_value;
  real value;
  
  assert (closed());

  ierr = MatNorm(mat, NORM_1, &petsc_value); CHKERRQ(ierr);

  value = static_cast<real>(petsc_value);

  return value;
};




real PetscMatrix::linfty_norm () const
{
  assert (initialized());
  
  int ierr=0;
  double petsc_value;
  real value;
  
  assert (closed());

  ierr = MatNorm(mat, NORM_INFINITY, &petsc_value); CHKERRQ(ierr);

  value = static_cast<real>(petsc_value);

  return value;
};


#endif
