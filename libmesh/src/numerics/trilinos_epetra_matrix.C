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
//EpetraMatrix members

template <typename T> 
void EpetraMatrix<T>::update_sparsity_pattern (const SparsityPattern::Graph &sparsity_pattern)
{  
  // clear data, start over
  this->clear ();    

  // big trouble if this fails!
  libmesh_assert (this->_dof_map != NULL);
  
  const unsigned int n_rows = sparsity_pattern.size();

  const unsigned int m   = this->_dof_map->n_dofs();
  const unsigned int n   = m;
  const unsigned int n_l = this->_dof_map->n_dofs_on_processor(libMesh::processor_id()); 
  const unsigned int m_l = n_l;

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
  _map = new Epetra_Map (m, 
                         m_l,
                         0,
                         Epetra_MpiComm (libMesh::COMM_WORLD));

  libmesh_assert (_map->NumGlobalPoints() == m);
  libmesh_assert (_map->MaxAllGID()+1 == m);
  
  const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
  const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();
  
   // Make sure the sparsity pattern isn't empty
  libmesh_assert (n_nz.size() == n_l);
  libmesh_assert (n_oz.size() == n_l);

  // Epetra wants the total number of nonzeros, both local and remote.
  std::vector<int> n_nz_tot; /**/ n_nz_tot.reserve(n_nz.size());
  
  for (unsigned int i=0; i<n_nz.size(); i++)    
    n_nz_tot.push_back(std::min(n_nz[i] + n_oz[i], n));
  
  if (m==0)
    return;

  _graph = new Epetra_CrsGraph(Copy, *_map, &n_nz_tot[0]);

  // Tell the matrix about its structure.  Initialize it
  // to zero.
  for (unsigned int i=0; i<n_rows; i++)
    _graph->InsertGlobalIndices(_graph->GRID(i),
                                sparsity_pattern[i].size(),
                                const_cast<int *>((const int *)&sparsity_pattern[i][0]));

  _graph->FillComplete();

  //Initialize the matrix
  libmesh_assert (!this->initialized());
  this->init ();
  libmesh_assert (this->initialized());
}



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

    libmesh_assert (this->_mat == NULL);
    libmesh_assert (this->_map == NULL);

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
  _map = new Epetra_Map (m, 
                         m_l,
                         0,
                         Epetra_MpiComm (libMesh::COMM_WORLD));
  
  libmesh_assert (_map->NumGlobalPoints() == m);
  libmesh_assert (_map->MaxAllGID()+1 == m);

  _mat = new Epetra_FECrsMatrix (Copy, *_map, nnz + noz);
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
  

  _mat = new Epetra_FECrsMatrix (Copy, *_graph);
}



template <typename T>
void EpetraMatrix<T>::zero ()
{
  libmesh_assert (this->initialized());
  
  _mat->Scale(0.0);
}



template <typename T>
void EpetraMatrix<T>::clear ()
{
//  delete _mat;
//  delete _map;

  this->_is_initialized = false;
  
  libmesh_assert (!this->initialized());
}



template <typename T>
Real EpetraMatrix<T>::l1_norm () const
{
  libmesh_assert (this->initialized());
  
  libmesh_assert (_mat != NULL);

  return static_cast<Real>(_mat->NormOne());
}



template <typename T>
Real EpetraMatrix<T>::linfty_norm () const
{
  libmesh_assert (this->initialized());
  
  
  libmesh_assert (_mat != NULL);

  return static_cast<Real>(_mat->NormInf());
}



template <typename T>
void EpetraMatrix<T>::print_matlab (const std::string name) const
{
  libmesh_assert (this->initialized());

  // libmesh_assert (this->closed());
  this->close();
  
  libmesh_not_implemented();

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

  const unsigned int m = dm.m();
  const unsigned int n = dm.n();

  libmesh_assert (rows.size() == m);
  libmesh_assert (cols.size() == n);

  _mat->SumIntoGlobalValues(m, (int *)&rows[0], n, (int *)&cols[0], &dm.get_values()[0]);
}





// template <typename T>
// void EpetraMatrix<T>::_get_submatrix(SparseMatrix<T>& submatrix,
// 				     const std::vector<unsigned int> &rows,
// 				     const std::vector<unsigned int> &cols,
// 				     const bool reuse_submatrix) const
// {
//   // Can only extract submatrices from closed matrices
//   this->close();
  
//   libmesh_not_implemented();

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
