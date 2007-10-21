// $Id: nonlinear_implicit_system.h,v 1.4 2007-10-21 20:48:44 benkirk Exp $

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



#ifndef __nonlinear_implicit_system_h__
#define __nonlinear_implicit_system_h__

// C++ includes

// Local Includes
#include "implicit_system.h"


// Forward declarations
template<typename T> class NonlinearSolver;


/**
 * This class provides a specific system class.  It aims
 * at implicit systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p ExplicitSystem.
 */

// ------------------------------------------------------------
// NonlinearImplicitSystem class definition

class NonlinearImplicitSystem : public ImplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  NonlinearImplicitSystem (EquationSystems& es,
			   const std::string& name,
			   const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~NonlinearImplicitSystem ();

  /**
   * The type of system.
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef ImplicitSystem Parent;
  
  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit ();
   
//   /**
//    * Prepares \p matrix and \p _dof_map for matrix assembly.
//    * Does not actually assemble anything.  For matrix assembly,
//    * use the \p assemble() in derived classes.
//    * @e Should be overloaded in derived classes.
//    */
//   virtual void assemble () { ImplicitSystem::assemble(); }
 
  /**
   * Assembles & solves the nonlinear system R(x) = 0.
   */
  virtual void solve ();
 
  /**
   * @returns \p "NonlinearImplicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "NonlinearImplicit"; }

  /**
   * The \p NonlinearSolver defines the interface used to
   * solve the nonlinear_implicit system.  This class handles all the
   * details of interfacing with various nonlinear algebra packages
   * like PETSc or LASPACK.
   */
  AutoPtr<NonlinearSolver<Number> > nonlinear_solver;
  
  /**
   * Returns  the number of iterations 
   * taken for the most recent nonlinear solve.
   */
  unsigned int n_nonlinear_iterations() const { return _n_nonlinear_iterations; }

  /**
   * Returns the final residual for the nonlinear system solve.
   */
  Real final_nonlinear_residual() const { return _final_nonlinear_residual; }
  
protected:
  
  /**
   * The number of nonlinear iterations required to solve the nonlinear
   * system R(x)=0.
   */
  unsigned int _n_nonlinear_iterations;

  /**
   * The final residual for the nonlinear system R(x)
   */
  Real _final_nonlinear_residual;
};



// ------------------------------------------------------------
// NonlinearImplicitSystem inline methods


#endif
