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



#ifndef LIBMESH_PETSC_PRECONDITIONER_H
#define LIBMESH_PETSC_PRECONDITIONER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/preconditioner.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/petsc_macro.h"

// Petsc includes
EXTERN_C_FOR_PETSC_BEGIN
#include "petscpc.h"
EXTERN_C_FOR_PETSC_END

// C++ includes

namespace libMesh
{

// forward declarations
template <typename T> class AutoPtr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class ShellMatrix;

/**
 * This class provides an interface to the suite of preconditioners available
 * from Petsc.
 *
 * @author Derek Gaston, 2009
 */

template <typename T>
class PetscPreconditioner : public Preconditioner<T>
{
public:

  /**
   *  Constructor. Initializes PetscPreconditioner data structures
   */
  PetscPreconditioner (const libMesh::Parallel::Communicator &comm
		       LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~PetscPreconditioner ();

  /**
   * Computes the preconditioned vector "y" based on input "x".
   * Usually by solving Py=x to get the action of P^-1 x.
   */
  virtual void apply(const NumericVector<T> & x, NumericVector<T> & y);

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init ();

  /**
   * Returns the actual Petsc PC struct.  Useful for more advanced
   * purposes
   */
  PC pc() { return _pc; }

  /**
   * Tells PETSC to use the user-specified preconditioner
   */
  static void set_petsc_preconditioner_type (const PreconditionerType & preconditioner_type, PC & pc);

protected:

  /**
   * Preconditioner context
   */
  PC _pc;

  /**
   * Petsc Matrix that's been pulled out of the _matrix object.
   * This happens during init...
   */
  Mat _mat;

private:
  /**
   * Some PETSc preconditioners (ILU, LU) don't work in parallel.  This function
   * is called from set_petsc_preconditioner_type() to set additional options
   * for those so-called sub-preconditioners.  This method ends up being static
   * so that it can be called from set_petsc_preconditioner_type().  Not sure
   * why set_petsc_preconditioner_type() needs to be static though...
   */
#if PETSC_VERSION_LESS_THAN(3,0,0)
  // In Petsc 2.3.3, PCType was #define'd as const char*
  static void set_petsc_subpreconditioner_type(PCType type, PC& pc);
#else
  // In later versions, PCType is #define'd as char*, so we need the const
  static void set_petsc_subpreconditioner_type(const PCType type, PC& pc);
#endif
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
PetscPreconditioner<T>::PetscPreconditioner (const libMesh::Parallel::Communicator &comm) :
  Preconditioner<T>(comm),
  _pc(PETSC_NULL)
{
}



template <typename T>
inline
PetscPreconditioner<T>::~PetscPreconditioner ()
{
  this->clear ();
}

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_PRECONDITIONER_H
