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

#include "libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

// C++ includes

// Local Includes
#include "auto_ptr.h"
#include "petsc_preconditioner.h"
#include "petsc_macro.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "petsc_macro.h"

#include "libmesh_common.h"

template <typename T>
void
PetscPreconditioner<T>::apply(const NumericVector<T> & x, NumericVector<T> & y)
{
  PetscVector<T> & x_pvec = libmesh_cast_ref<PetscVector<T>&>(const_cast<NumericVector<T>&>(x));
  PetscVector<T> & y_pvec = libmesh_cast_ref<PetscVector<T>&>(const_cast<NumericVector<T>&>(y));

  Vec x_vec = x_pvec.vec();
  Vec y_vec = y_pvec.vec();

  PCApply(_pc,x_vec,y_vec);
}

template <typename T>
void
PetscPreconditioner<T>::init ()
{
  if(!this->_matrix)
  {
    std::cerr << "ERROR: No matrix set for PetscPreconditioner, but init() called" << std::endl;
    libmesh_error();
  }

  //Clear the preconditioner in case it has been created in the past
  if(!this->_is_initialized)
  {
    //Create the preconditioning object
    PCCreate(libMesh::COMM_WORLD,&_pc);

    //Set the PCType
    set_petsc_preconditioner_type(this->_preconditioner_type, _pc);

#ifdef LIBMESH_HAVE_PETSC_HYPRE
    if(this->_preconditioner_type == AMG_PRECOND)
      PCHYPRESetType(this->_pc, "boomeramg");
#endif

    PetscMatrix<T> * pmatrix = libmesh_cast_ptr<PetscMatrix<T>*, SparseMatrix<T> >(this->_matrix);
    
    _mat = pmatrix->mat();
  }
  
  PCSetOperators(_pc,_mat,_mat,SAME_NONZERO_PATTERN);

  this->_is_initialized = true;
}

template <typename T>
void
PetscPreconditioner<T>::set_petsc_preconditioner_type (const PreconditionerType & preconditioner_type, PC & pc)
{
  int ierr = 0;
 
  switch (preconditioner_type)
  {
  case IDENTITY_PRECOND:
    ierr = PCSetType (pc, (char*) PCNONE);      CHKERRABORT(libMesh::COMM_WORLD,ierr); break;
	
  case CHOLESKY_PRECOND:
    ierr = PCSetType (pc, (char*) PCCHOLESKY);  CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case ICC_PRECOND:
    ierr = PCSetType (pc, (char*) PCICC);       CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case ILU_PRECOND:
    ierr = PCSetType (pc, (char*) PCILU);       CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case LU_PRECOND:
    ierr = PCSetType (pc, (char*) PCLU);        CHKERRABORT(libMesh::COMM_WORLD,ierr); break;
      
  case ASM_PRECOND:
    ierr = PCSetType (pc, (char*) PCASM);       CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case JACOBI_PRECOND:
    ierr = PCSetType (pc, (char*) PCJACOBI);    CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case BLOCK_JACOBI_PRECOND:
    ierr = PCSetType (pc, (char*) PCBJACOBI);   CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case SOR_PRECOND:
    ierr = PCSetType (pc, (char*) PCSOR);       CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case EISENSTAT_PRECOND:
    ierr = PCSetType (pc, (char*) PCEISENSTAT); CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  case AMG_PRECOND:
    ierr = PCSetType (pc, (char*) PCHYPRE);     CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

#if !(PETSC_VERSION_LESS_THAN(2,1,2))
    // Only available for PETSC >= 2.1.2      
  case USER_PRECOND:
    ierr = PCSetType (pc, (char*) PCMAT);       CHKERRABORT(libMesh::COMM_WORLD,ierr); break;
#endif

  case SHELL_PRECOND:
    ierr = PCSetType (pc, (char*) PCSHELL);     CHKERRABORT(libMesh::COMM_WORLD,ierr); break;

  default:
    std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
              << preconditioner_type       << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }

  //Let the commandline override stuff
  if( preconditioner_type != AMG_PRECOND )
    PCSetFromOptions(pc);
}

//------------------------------------------------------------------
// Explicit instantiations
template class PetscPreconditioner<Number>;

#endif // #ifdef LIBMESH_HAVE_PETSC
