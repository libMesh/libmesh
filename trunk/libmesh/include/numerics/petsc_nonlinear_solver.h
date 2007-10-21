// $Id: petsc_nonlinear_solver.h,v 1.4 2007-10-21 20:48:43 benkirk Exp $

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



#ifndef __petsc_nonlinear_solver_h__
#define __petsc_nonlinear_solver_h__

// C++ includes

// Local includes
#include "nonlinear_solver.h"


// Petsc include files.
#ifdef HAVE_PETSC

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscsnes.h>
}
#else
# include <petscsnes.h>
#endif



/**
 * This class provides an interface to PETSc
 * iterative solvers that is compatible with the \p libMesh
 * \p NonlinearSolver<>
 *
 * @author Benjamin Kirk, 2002-2007
 */

template <typename T>
class PetscNonlinearSolver : public NonlinearSolver<T>
{
public:
  /**
   *  Constructor. Initializes Petsc data structures
   */
  PetscNonlinearSolver ();
    
  /**
   * Destructor.
   */
  ~PetscNonlinearSolver ();
  
  /**
   * Release all memory and clear data structures.
   */
  virtual void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init ();
  
  /**
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T>&,    // System Jacobian Matrix
					       NumericVector<T>&,   // Solution vector
					       NumericVector<T>&,   // Residual vector
					       const double,        // Stopping tolerance
					       const unsigned int); // N. Iterations


  
private:

  /**
   * Nonlinear solver context
   */
  SNES _snes;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
PetscNonlinearSolver<T>::PetscNonlinearSolver ()
{
}



template <typename T>
inline
PetscNonlinearSolver<T>::~PetscNonlinearSolver ()
{
  this->clear ();
}



#endif // #ifdef HAVE_PETSC
#endif // #ifdef __petsc_nonlinear_solver_h__
