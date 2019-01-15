// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PRECONDITIONER_H
#define LIBMESH_PRECONDITIONER_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_object.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum SolverPackage : int;
enum PreconditionerType : int;
}
#else
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_preconditioner_type.h"
#endif


// C++ includes
#include <cstddef>

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class ShellMatrix;

/**
 * This class provides a uniform interface for preconditioners.  This base
 * class can be inherited from to wrap preconditioners from different packages
 * like PETSc or Trilinos.
 *
 * In the below comments P is the matrix to be preconditioned with Apply()
 * performing the equivalent of the matrix vector product P^-1 x.  This
 * can also be thought of as (usually approximately) solving for Py=x.
 *
 * \author Derek Gaston
 * \date 2009
 */
template <typename T>
class Preconditioner : public ReferenceCountedObject<Preconditioner<T>>,
                       public ParallelObject
{
public:

  /**
   * Constructor. Initializes Preconditioner data structures.
   */
  Preconditioner (const libMesh::Parallel::Communicator & comm);

  /**
   * Destructor.
   */
  virtual ~Preconditioner ();

  /**
   * Builds a \p Preconditioner using the linear solver package
   * specified by \p solver_package, returning the result wrapped in a
   * std::unique_ptr for safety.
   */
  static std::unique_ptr<Preconditioner<T>>
  build_preconditioner(const libMesh::Parallel::Communicator & comm,
                       const SolverPackage solver_package = libMesh::default_solver_package());

  /**
   * Builds a \p Preconditioner using the linear solver package specified by
   * \p solver_package
   *
   * \deprecated Use build_preconditioner() instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  static Preconditioner<T> *
  build(const libMesh::Parallel::Communicator & comm,
        const SolverPackage solver_package = libMesh::default_solver_package());
#endif

  /**
   * \returns \p true if the data structures are initialized, \p false
   * otherwise.
   */
  bool initialized () const { return _is_initialized; }

  /**
   * Computes the preconditioned vector \p y based on input vector \p
   * x. This is usually done by solving \f$ Py=x \f$ to get the
   * action of \f$ P^-1 x \f$.
   */
  virtual void apply(const NumericVector<T> & x, NumericVector<T> & y) = 0;

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () {}

  /**
   * Initialize data structures if not done so already.
   *
   * \note This MUST be called before the preconditioning object is used.
   */
  virtual void init () {}

  /**
   * This is called every time the "operator might have changed".
   *
   * This is where you need to fill in your preconditioning matrix.
   */
  virtual void setup () {}

  /**
   * Sets the matrix to be preconditioned.
   */
  void set_matrix(SparseMatrix<Number> & mat);

  /**
   * \returns The type of preconditioner to use.
   */
  PreconditionerType type () const { return _preconditioner_type; }

  /**
   * Sets the type of preconditioner to use.
   */
  void set_type (const PreconditionerType pct);

protected:

  /**
   * The matrix P... ie the matrix to be preconditioned.
   * This is often the actual system matrix of a linear system.
   */
  SparseMatrix<T> * _matrix;

  /**
   * Enum stating with type of preconditioner to use.
   */
  PreconditionerType _preconditioner_type;

  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
Preconditioner<T>::~Preconditioner ()
{
  this->clear ();
}

template <typename T>
void
Preconditioner<T>::set_matrix(SparseMatrix<Number> & mat)
{
  //If the matrix is changing then we (probably) need to reinitialize.
  _is_initialized = false;
  _matrix = &mat;
}

template <typename T>
void
Preconditioner<T>::set_type (const PreconditionerType pct)
{
  //If the preconditioner type changes we (probably) need to reinitialize.
  _is_initialized = false;
  _preconditioner_type = pct;
}

} // namespace libMesh


#endif // LIBMESH_PRECONDITIONER_H
