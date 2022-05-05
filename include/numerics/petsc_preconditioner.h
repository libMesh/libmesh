// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libMesh includes
#include "libmesh/preconditioner.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/wrapped_petsc.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum PreconditionerType : int;
}
#else
#include "libmesh/enum_preconditioner_type.h"
#endif

// Petsc includes
#include "petscpc.h"

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class ShellMatrix;

/**
 * This class provides an interface to the suite of preconditioners
 * available from PETSc. All overridden virtual functions are
 * documented in preconditioner.h.
 *
 * \author Derek Gaston
 * \date 2009
 */
template <typename T>
class PetscPreconditioner : public Preconditioner<T>
{
public:

  /**
   *  Constructor. Initializes PetscPreconditioner data structures
   */
  PetscPreconditioner (const libMesh::Parallel::Communicator & comm_in);

  virtual ~PetscPreconditioner () = default;

  virtual void apply(const NumericVector<T> & x, NumericVector<T> & y) override;

  virtual void clear () override;

  virtual void init () override;

  /**
   * \returns The PETSc PC object.  Can be useful for implementing
   * more advanced algorithms.
   */
  PC pc();

  /**
   * Tells PETSc to use the user-specified preconditioner.
   */
  static void set_petsc_preconditioner_type (const PreconditionerType & preconditioner_type, PC & pc);

protected:

  /**
   * Preconditioner context
   */
  WrappedPetsc<PC> _pc;

  /**
   * PETSc Mat pulled out of the _matrix object during init(). We
   * aren't responsible for cleaning up this one.
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
  static void set_petsc_subpreconditioner_type(const PCType type, PC & pc);
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_PRECONDITIONER_H
