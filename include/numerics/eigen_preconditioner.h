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



#ifndef LIBMESH_EIGEN_PRECONDITIONER_H
#define LIBMESH_EIGEN_PRECONDITIONER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_EIGEN

// Local includes
#include "libmesh/eigen_core_support.h"
#include "libmesh/preconditioner.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"

// Eigen includes

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
 * from Eigen.
 *
 * @author Benjamin Kirk, 2013.
 */

template <typename T>
class EigenPreconditioner : public Preconditioner<T>
{
public:

  /**
   *  Constructor. Initializes EigenPreconditioner data structures
   */
  EigenPreconditioner (const libMesh::Parallel::Communicator &comm
		       LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~EigenPreconditioner ();

  /**
   * Computes the preconditioned vector "y" based on input "x".
   * Usually by solving Py=x to get the action of P^-1 x.
   */
  virtual void apply(const NumericVector<T> & x, NumericVector<T> & y);

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () {}

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init ();

protected:
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
EigenPreconditioner<T>::EigenPreconditioner (const libMesh::Parallel::Communicator &comm) :
  Preconditioner<T>(comm)
{
}



template <typename T>
inline
EigenPreconditioner<T>::~EigenPreconditioner ()
{
  this->clear ();
}

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_EIGEN
#endif // LIBMESH_EIGEN_PRECONDITIONER_H
