// $Id: petsc_preconditioner.h 3158 2009-02-02 17:44:57Z drgasto $

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



#ifndef __petsc_preconditioner_h__
#define __petsc_preconditioner_h__

#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// C++ includes

// Local includes
#include "preconditioner.h"
#include "libmesh_common.h"
#include "enum_solver_package.h"
#include "enum_preconditioner_type.h"
#include "reference_counted_object.h"
#include "libmesh.h"

// Petsc includes
#include "petscpc.h"

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
  PetscPreconditioner ();
    
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
  virtual void clear () {}

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
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
PetscPreconditioner<T>::PetscPreconditioner () :
  Preconditioner<T>()
{
}



template <typename T>
inline
PetscPreconditioner<T>::~PetscPreconditioner ()
{
  this->clear ();
}

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
