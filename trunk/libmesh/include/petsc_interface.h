//    $Id: petsc_interface.h,v 1.4 2003-01-21 19:24:35 benkirk Exp $

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
#include "petsc_vector.h"
#include "petsc_matrix.h"


// Petsc include files
#ifdef USE_COMPLEX_NUMBERS
#  define PETSC_USE_COMPLEX 1
#endif

namespace Petsc {
extern "C" {
#include <petsc.h>
#include <petscsles.h>  
}
// for easy switching between Petsc 2.1.0/2.1.1
//typedef Scalar PetscScalar;
} 


using namespace Petsc;



/**
 * This class provides a deal.II interface to the Petsc
 * iterative solver library.
 *
 * @author Benjamin Kirk, 2002
 */

class PetscInterface
{
 public:
  /**
   *  Constructor. Initializes Petsc data structures
   **/
  PetscInterface ();
    
  /**
   * Destructor.
   **/
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
   *  Use no preconditioning
   **/
  void precond_identity ();

  /**
   * Use Jacobi preconditioning
   **/
  void precond_jacobi ();

  /**
   * Use block Jacobi preconditioning
   **/
  void precond_block_jacobi ();

  /**
   * Use SOR preconditioning
   **/
  void precond_sor ();

  /**
   * Use SOR with Eisenstat trick preconditioning
   **/
  void precond_eisenstat ();

  /**
   * Use Additive Schwarz preconditioning
   **/
  void precond_asm ();

  /**
   * Use Cholesky preconditioning
   **/
  void precond_cholesky ();

  /**
   * Use Incomplete Cholesky preconditioning
   **/
  void precond_icc ();

  /**
   * Use ILU (Incomplete LU factorization) preconditioning
   **/
  void precond_ilu ();

  /**
   * Use LU "preconditioning." Can be combined with the Richardson solver
   * (below) to provide a direct solver.
   **/
  void precond_lu ();

  // Methods to set the linear solver type --------------------------------
  
  /**
   * Use the Conjugate Gradient solver
   **/
  void solver_cg ();

  /**
   * Use the Conjugate Residual solver
   **/
  void solver_cr ();

  /**
   * Use the Conjugate Gradient Squared solver
   **/
  void solver_cgs ();

  /**
   * Use the Bi-Conjugate Gradient solver
   **/
  void solver_bicg ();

  /**
   * Use the Transpose-Free Quasi-Minimum residual (2) solver
   **/
  void solver_tcqmr ();

  /**
   * Use the Transpose-Free Quasi-Minimum residual (1) solver
   **/
  void solver_tfqmr ();

  /**
   * Use the Least-Squares method solver
   **/
  void solver_lsqr ();

  /**
   * Use the Bi-Conjugate Gradient Stabilized solver
   **/
  void solver_bcgstab ();

  /**
   * Use the Minimum-Residual solver
   **/
  void solver_minres ();

  /**
   * Use the GMRES solver
   **/
  void solver_gmres ();

  /**
   * Use the Richardson solver
   **/
  void solver_richardson ();

  /**
   * Use the Chebychev solver
   **/
  void solver_chebychev ();

  /**
   * Call the Petsc solver
   **/
    
  std::pair<unsigned int, real> 
    solve (PetscMatrix &matrix,
	   PetscVector &solution,
	   PetscVector &rhs,
	   const double tol,
	   const unsigned int m_its);
   
 private:

  /**
   * Linear solver context
   **/
  SLES sles;
    
  /**
   * Preconditioner context
   **/
  PC   pc; 

  /**
   * Krylov subspace context
   **/
  KSP  ksp;

  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool initialized;
};


/*----------------------- functions ----------------------------------*/
inline
PetscInterface::PetscInterface () :
  initialized (false)
{
};




inline
PetscInterface::~PetscInterface ()
{
  clear ();
};



inline
void PetscInterface::precond_identity()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCNONE); CHKERRQ(ierr);

  return;
};





inline
void PetscInterface::precond_cholesky()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCCHOLESKY); CHKERRQ(ierr);

  return;
};





inline
void PetscInterface::precond_icc()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCICC); CHKERRQ(ierr);

  return;
};





inline
void PetscInterface::precond_ilu()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCILU); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::precond_lu()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCLU); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::precond_asm()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCASM); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::precond_jacobi()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCJACOBI); CHKERRQ(ierr);

  return;
};



inline
void PetscInterface::precond_block_jacobi()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCBJACOBI); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::precond_sor()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCSOR); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::precond_eisenstat()
{
  init();
  
  int ierr=0;
  
  ierr = PCSetType (pc, (char*) PCEISENSTAT); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_cg()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPCG); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_cr()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPCR); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_cgs()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPCGS); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_bicg()
{
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPBICG); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_tcqmr()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPTCQMR); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_tfqmr()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPTFQMR); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_lsqr()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPLSQR); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_bcgstab()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPBCGS); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_minres()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPMINRES); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_gmres()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPGMRES); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_richardson()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPRICHARDSON); CHKERRQ(ierr);

  return;
};




inline
void PetscInterface::solver_chebychev()
{
  init();
  
  int ierr=0;
  
  ierr = KSPSetType (ksp, (char*) KSPCHEBYCHEV); CHKERRQ(ierr);

  return;
};



#endif
#endif
