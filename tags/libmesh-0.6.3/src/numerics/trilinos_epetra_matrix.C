// $Id: petsc_matrix.C 2789 2008-04-13 02:24:40Z roystgnr $

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


// C++ includes
#include "libmesh_config.h"

#ifdef HAVE_TRILINOS

// Local includes
#include "trilinos_epetra_matrix.h"
#include "dof_map.h"
#include "dense_matrix.h"
#include "parallel.h"



//-----------------------------------------------------------------------
// EpetraMatrix members
template <typename T>
void EpetraMatrix<T>::init (const unsigned int m,
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

    libmesh_assert (this->_mat.get() == NULL);
    libmesh_assert (this->_map.get() == NULL);

    this->_is_initialized = true;
  }

  // error checking
#ifndef NDEBUG
  {
    libmesh_assert (n == m);
    libmesh_assert (n_l == m_l);

    unsigned int 
      summed_m_l = m_l,
      summed_n_l = n_l;

    Parallel::sum (summed_m_l);
    Parallel::sum (summed_n_l);

    libmesh_assert (m == summed_m_l);
    libmesh_assert (n == summed_n_l);
  }
#endif

  // build a map defining the data distribution
  _map.reset (new Epetra_Map (m, 
			      m_l,
			      0,
			      Epetra_MpiComm (libMesh::COMM_WORLD))
	      );
  
  libmesh_assert (_map->NumGlobalPoints() == m);
  libmesh_assert (_map->MaxAllGID()+1 == m);

  _mat.reset (new Epetra_FECrsMatrix (Copy, *_map, nnz + noz));

//   this->zero ();
}




template <typename T>
void EpetraMatrix<T>::init ()
{
  libmesh_assert (this->_dof_map != NULL);
  
  {
    // Clear initialized matrices
    if (this->initialized())
      this->clear();

    this->_is_initialized = true;
  }

  
  LIBMESH_THROW(libMesh::NotImplemented());

//   int proc_id = 0;

//   MPI_Comm_rank (libMesh::COMM_WORLD, &proc_id);
  
//   const unsigned int m   = this->_dof_map->n_dofs();
//   const unsigned int n   = m;
//   const unsigned int n_l = this->_dof_map->n_dofs_on_processor(proc_id); 
//   const unsigned int m_l = n_l;


//   const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
//   const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();

//   // Make sure the sparsity pattern isn't empty
//   libmesh_assert (n_nz.size() == n_l);
//   libmesh_assert (n_oz.size() == n_l);
  
//   if (m==0)
//     return;
  
//   int ierr     = 0;
//   int m_global = static_cast<int>(m);
//   int n_global = static_cast<int>(n);
//   int m_local  = static_cast<int>(m_l);
//   int n_local  = static_cast<int>(n_l);


//   // create a sequential matrix on one processor
//   if ((m_l == m) && (n_l == n))
//     {
//       if (n_nz.empty())
//         ierr = MatCreateSeqAIJ (libMesh::COMM_WORLD, n_global, n_global,
// 			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
//       else
//         ierr = MatCreateSeqAIJ (libMesh::COMM_WORLD, n_global, n_global,
// 			        PETSC_NULL, (int*) &n_nz[0], &_mat);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//       ierr = MatSetFromOptions (_mat);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }

//   else
//     {
//       if (n_nz.empty())
//         ierr = MatCreateMPIAIJ (libMesh::COMM_WORLD,
// 			        m_local, n_local,
// 			        m_global, n_global,
// 			        PETSC_NULL, (int*) PETSC_NULL,
// 			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
//       else
//         ierr = MatCreateMPIAIJ (libMesh::COMM_WORLD,
// 			        m_local, n_local,
// 			        m_global, n_global,
// 			        PETSC_NULL, (int*) &n_nz[0],
// 			        PETSC_NULL, (int*) &n_oz[0], &_mat);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//       ierr = MatSetFromOptions (_mat);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }

//   this->zero();
}



template <typename T>
void EpetraMatrix<T>::zero ()
{
  libmesh_assert (this->initialized());
  
  LIBMESH_THROW(libMesh::NotImplemented());

//   int ierr=0;

//   ierr = MatZeroEntries(_mat);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
}



template <typename T>
void EpetraMatrix<T>::clear ()
{
  this->_mat.reset();
  this->_map.reset();

  this->_is_initialized = false;
  
  libmesh_assert (!this->initialized());
}



template <typename T>
Real EpetraMatrix<T>::l1_norm () const
{
  libmesh_assert (this->initialized());
  
  libmesh_assert (_mat.get() != NULL);

  return static_cast<Real>(_mat->NormOne());
}



template <typename T>
Real EpetraMatrix<T>::linfty_norm () const
{
  libmesh_assert (this->initialized());
  
  
  libmesh_assert (_mat.get() != NULL);

  return static_cast<Real>(_mat->NormInf());
}



template <typename T>
void EpetraMatrix<T>::print_matlab (const std::string name) const
{
  libmesh_assert (this->initialized());

  // libmesh_assert (this->closed());
  this->close();
  
  LIBMESH_THROW(libMesh::NotImplemented());

//   int ierr=0; 
//   PetscViewer petsc_viewer;


//   ierr = PetscViewerCreate (libMesh::COMM_WORLD,
// 			    &petsc_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   /**
//    * Create an ASCII file containing the matrix
//    * if a filename was provided.  
//    */
//   if (name != "NULL")
//     {
//       ierr = PetscViewerASCIIOpen( libMesh::COMM_WORLD,
// 				   name.c_str(),
// 				   &petsc_viewer);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
//       ierr = PetscViewerSetFormat (petsc_viewer,
// 				   PETSC_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//       ierr = MatView (_mat, petsc_viewer);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }

//   /**
//    * Otherwise the matrix will be dumped to the screen.
//    */
//   else
//     {
//       ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
// 				   PETSC_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//       ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }


//   /**
//    * Destroy the viewer.
//    */
//   ierr = PetscViewerDestroy (petsc_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
}




template <typename T>
void EpetraMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				 const std::vector<unsigned int>& rows,
				 const std::vector<unsigned int>& cols)
{
  libmesh_assert (this->initialized());
  
  LIBMESH_THROW(libMesh::NotImplemented());

//   const unsigned int m = dm.m();
//   const unsigned int n = dm.n();

//   libmesh_assert (rows.size() == m);
//   libmesh_assert (cols.size() == n);
  
//   int ierr=0;

//   // These casts are required for PETSc <= 2.1.5
//   ierr = MatSetValues(_mat,
// 		      m, (int*) &rows[0],
// 		      n, (int*) &cols[0],
// 		      (PetscScalar*) &dm.get_values()[0],
// 		      ADD_VALUES);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
}





// template <typename T>
// void EpetraMatrix<T>::_get_submatrix(SparseMatrix<T>& submatrix,
// 				     const std::vector<unsigned int> &rows,
// 				     const std::vector<unsigned int> &cols,
// 				     const bool reuse_submatrix) const
// {
//   // Can only extract submatrices from closed matrices
//   this->close();
  
//   LIBMESH_THROW(libMesh::NotImplemented());

// //   // Attempt to cast the input matrix to a EpetraMatrix*
// //   EpetraMatrix<T>* petsc_submatrix = dynamic_cast<EpetraMatrix<T>*>(&submatrix);
// //   libmesh_assert(petsc_submatrix != NULL);

// //   // Construct row and column index sets.
// //   int ierr=0;
// //   IS isrow, iscol;

// //   ierr = ISCreateGeneral(libMesh::COMM_WORLD,
// // 			 rows.size(),
// // 			 (int*) &rows[0],
// // 			 &isrow); CHKERRABORT(libMesh::COMM_WORLD,ierr);

// //   ierr = ISCreateGeneral(libMesh::COMM_WORLD,
// // 			 cols.size(),
// // 			 (int*) &cols[0],
// // 			 &iscol); CHKERRABORT(libMesh::COMM_WORLD,ierr);

// //   // Extract submatrix
// //   ierr = MatGetSubMatrix(_mat,
// // 			 isrow,
// // 			 iscol,
// // 			 PETSC_DECIDE,
// // 			 (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
// // 			 &(petsc_submatrix->_mat));  CHKERRABORT(libMesh::COMM_WORLD,ierr);

// //   // Specify that the new submatrix is initialized and close it.
// //   petsc_submatrix->_is_initialized = true;
// //   petsc_submatrix->close();

// //   // Clean up PETSc data structures
// //   ierr = ISDestroy(isrow); CHKERRABORT(libMesh::COMM_WORLD,ierr);
// //   ierr = ISDestroy(iscol); CHKERRABORT(libMesh::COMM_WORLD,ierr);
// }







//------------------------------------------------------------------
// Explicit instantiations
template class EpetraMatrix<Number>;


#endif // #ifdef HAVE_TRILINOS
