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



#ifndef __preconditioner_h__
#define __preconditioner_h__


// C++ includes

// Local includes
#include "libmesh_common.h"
#include "enum_solver_package.h"
#include "enum_solver_type.h"
#include "enum_preconditioner_type.h"
#include "reference_counted_object.h"
#include "libmesh.h"

namespace libMesh
{

// forward declarations
template <typename T> class AutoPtr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class ShellMatrix;





/**
 * This class provides a uniform interface for preconditioners.  This base
 * class is overloaded to provide linear solvers from different packages
 * like PETSC or Trilinos.
 *
 * In the below comments P is the matrix to be preconditioned with Apply()
 * performing the equivalent of the matrix vector product P^-1 x.  This
 * can also be thought of as (usually approximately) solving for Py=x.
 *
 * @author Derek Gaston, 2009
 */

template <typename T>
class Preconditioner : public ReferenceCountedObject<Preconditioner<T> >
{
public:

  /**
   *  Constructor. Initializes Preconditioner data structures
   */
  Preconditioner ();

  /**
   * Destructor.
   */
  virtual ~Preconditioner ();

  /**
   * Builds a \p Preconditioner using the linear solver package specified by
   * \p solver_package
   */
  static Preconditioner<T> * build(const SolverPackage solver_package =
                                   libMesh::default_solver_package());

  /**
   * @returns true if the data structures are
   * initialized, false otherwise.
   */
  bool initialized () const { return _is_initialized; }

  /**
   * Computes the preconditioned vector "y" based on input "x".
   * Usually by solving Py=x to get the action of P^-1 x.
   */
  virtual void apply(const NumericVector<T> & x, NumericVector<T> & y) = 0;

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () {}

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init () {};

  /**
   * Sets the matrix P to be preconditioned.
   */
  void set_matrix(SparseMatrix<Number> & mat);

  /**
   * Returns the type of preconditioner to use.
   */
  PreconditionerType type () const
  { return _preconditioner_type; }

  /**
   * Sets the type of preconditioner to use.
   */
  void set_type (const PreconditionerType pct);

protected:

  /**
   * The matrix P... ie the matrix to be preconditioned.
   * This is often the actual system matrix of a linear sytem.
   */
  SparseMatrix<T> * _matrix;

  /**
   * Enum statitng with type of preconditioner to use.
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
Preconditioner<T>::Preconditioner () :
  _matrix(NULL),
  _preconditioner_type (ILU_PRECOND),
  _is_initialized      (false)
{
}



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


#endif // #ifdef __preconditioner_h__
