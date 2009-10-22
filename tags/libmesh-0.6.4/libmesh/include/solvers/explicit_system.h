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



#ifndef __explicit_system_h__
#define __explicit_system_h__

// C++ includes

// Local Includes
#include "system.h"


// Forward Declarations


/**
 * This class provides a specific system class.  It aims
 * at explicit systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p System.
 */

// ------------------------------------------------------------
// ExplicitSystem class definition

class ExplicitSystem : public System
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ExplicitSystem (EquationSystems& es,
		  const std::string& name,
		  const unsigned int number);

  /**
   * Destructor.
   */
  ~ExplicitSystem ();
 
  /**
   * The type of system.
   */
  typedef ExplicitSystem sys_type;

  /**
   * The type of the parent
   */
  typedef System Parent;
  
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
  
  /**
   * Prepares \p qoi for quantity of interest assembly, then calls
   * user qoi function.
   * @e Can be overloaded in derived classes.
   */
  virtual void assemble_qoi ();

  /**
   * Prepares \p rhs for quantity of interest derivative assembly,
   * then calls user qoi derivative function.
   * @e Can be overloaded in derived classes.
   */
  virtual void assemble_qoi_derivative ();
 
  /**
   * Assembles & solves the linear system Ax=b. 
   */
  virtual void solve ();
 
  /**
   * Yells at any user who tried to use an adjoint without
   * providing more than an explicit system assembly.
   */
  virtual void adjoint_solve () { libmesh_error(); }
 
  /**
   * Yells at any user who tried to calculate parameter sensitivities
   * without providing more than an explicit system assembly.
   */
  virtual void qoi_parameter_sensitivity (std::vector<Number *>&, std::vector<Number>&)
    { libmesh_error(); }

  /**
   * @returns \p "Explicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Explicit"; }

  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * right-hand-side vector b.
   */
  NumericVector<Number> * rhs;
  

protected:

  
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  
private:

  
  /**
   * Add the system right-hand-side vector to the \p _vectors data structure.
   * Useful in initialization.
   */
  void add_system_rhs ();
};



// ------------------------------------------------------------
// ExplicitSystem inline methods


#endif
