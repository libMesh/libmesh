// $Id: petsc_interface.h,v 1.9 2003-02-13 22:56:07 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __petsc_interface_h__
#define __petsc_interface_h__

#include "mesh_common.h"

#ifdef HAVE_PETSC


// C++ includes


// Local includes
#include "linear_solver_interface.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"


/**
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
# ifndef USE_COMPLEX_NUMBERS

namespace Petsc {
extern "C" {
#include <petscsles.h>
}
// for easy switching between Petsc 2.1.0/2.1.1
// typedef Scalar PetscScalar;
} 
using namespace Petsc;

#else

#include <petscsles.h>

#endif




/**
 * This class provides a deal.II interface to the Petsc
 * iterative solver library.
 *
 * @author Benjamin Kirk, 2002
 */

class PetscInterface : public LinearSolverInterface
{
 public:
  /**
   *  Constructor. Initializes Petsc data structures
   */
  PetscInterface ();
    
  /**
   * Destructor.
   */
  ~PetscInterface ();
  
  /**
   * Release all memory and clear data structures.
   */
  void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  void init ();
  
  /**
   * Call the Petsc solver
   */    
  std::pair<unsigned int, Real> 
    solve (SparseMatrix &matrix,
	   NumericVector &solution,
	   NumericVector &rhs,
	   const double tol,
	   const unsigned int m_its);
   
 private:

  /**
   * Tells PETSC to use the user-specified solver stored in
   * \p _solver_type
   */
  void set_petsc_solver_type ();

  /**
   * Tells PETSC to use the user-specified preconditioner stored in
   * \p _preconditioner_type
   */
  void set_petsc_preconditioner_type ();

  /**
   * Linear solver context
   */
  SLES _sles;
    
  /**
   * Preconditioner context
   */
  PC _pc; 

  /**
   * Krylov subspace context
   */
  KSP _ksp;
};


/*----------------------- functions ----------------------------------*/
inline
PetscInterface::PetscInterface ()
{
}



inline
PetscInterface::~PetscInterface ()
{
  clear ();
}



#endif
#endif
