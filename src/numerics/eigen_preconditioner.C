// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// C++ includes

// Local Includes
#include "libmesh/eigen_preconditioner.h"
#include "libmesh/eigen_sparse_matrix.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/libmesh_common.h"

namespace libMesh
{

template <typename T>
void EigenPreconditioner<T>::apply(const NumericVector<T> & /* x */, NumericVector<T> & /* y */)
{
  libmesh_not_implemented();
  
  // EigenSparseVector<T> & x_pvec = libmesh_cast_ref<EigenSparseVector<T>&>(const_cast<NumericVector<T>&>(x));
  // EigenSparseVector<T> & y_pvec = libmesh_cast_ref<EigenSparseVector<T>&>(const_cast<NumericVector<T>&>(y));

  // Vec x_vec = x_pvec.vec();
  // Vec y_vec = y_pvec.vec();

  // int ierr = PCApply(_pc,x_vec,y_vec);
  // CHKERRABORT(libMesh::COMM_WORLD,ierr);
}




template <typename T>
void EigenPreconditioner<T>::init ()
{
  libmesh_not_implemented();
  
  // if(!this->_matrix)
  // {
  //   libMesh::err << "ERROR: No matrix set for EigenPreconditioner, but init() called" << std::endl;
  //   libmesh_error();
  // }

  // // Clear the preconditioner in case it has been created in the past
  // if (!this->_is_initialized)
  // {
  //   // Should probably use PCReset(), but it's not working at the moment so we'll destroy instead
  //   if (_pc)
  //   {
  //     int ierr = LibMeshPCDestroy(&_pc);
  //     CHKERRABORT(libMesh::COMM_WORLD,ierr);
  //   }

  //   int ierr = PCCreate(libMesh::COMM_WORLD,&_pc);
  //   CHKERRABORT(libMesh::COMM_WORLD,ierr);

  //   EigenSparseMatrix<T> * pmatrix = libmesh_cast_ptr<EigenSparseMatrix<T>*, SparseMatrix<T> >(this->_matrix);

  //   _mat = pmatrix->mat();
  // }

  // int ierr = PCSetOperators(_pc,_mat,_mat,SAME_NONZERO_PATTERN);
  // CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // // Set the PCType.  Note: this used to be done *before* the call to
  // // PCSetOperators(), and only when !_is_initialized, but
  // // 1.) Some preconditioners (those employing sub-preconditioners,
  // // for example) have to call PCSetUp(), and can only do this after
  // // the operators have been set.
  // // 2.) It should be safe to call set_petsc_preconditioner_type()
  // // multiple times.
  // set_petsc_preconditioner_type(this->_preconditioner_type, _pc);

  // this->_is_initialized = true;
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenPreconditioner<Number>;

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_EIGEN
