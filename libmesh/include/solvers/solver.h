// $Id: solver.h,v 1.1 2004-01-03 15:37:42 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __solver_h__
#define __solver_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "reference_counted_object.h"
#include "equation_systems.h"

// Forward Declarations
class Mesh;
class Solver;

/**
 * This is a generic class that defines a solver to be used in a
 * simulation.  A user can define a solver by deriving from this
 * class and implementing certain functions.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// Solver class definition
class Solver : public ReferenceCountedObject<Solver>
{
protected:
  
  /**
   * Constructor. Requires a reference to the system
   * to be solved. The constructor is protected since
   * it should not be instantiated by users.
   */
  Solver (EquationSystems& es);

  /**
   * Constructor.  Requires a reference to the \p EquationSystems
   * object, a name for the system, and the system number.
   */
  Solver (EquationSystems& es,
	  const std::string& name,
	  const unsigned int number);
  

public:
  
  /**
   * Destructor.
   */
  ~Solver ();

  /**
   * The type of system
   */
  typedef EquationSystems sys_type;
  
  /**
   * The initialization function.  This method is used to
   * initialize data structures befor a simulation begins.
   */
  virtual void init ();

  /**
   * This method may be called before each solve step in order
   * to perform any required pre-processing.
   */
  virtual void pre_process ();

  /**
   * This method performs a solve step.  What occurs in
   * this method will depend on the type of solver.  See
   * the example programs for more details.
   */
  virtual void solve ();

  /**
   * This method may be called after each solve step in order
   * to perform any required post-processing.
   */
  virtual void post_process ();

  /**
   * @returns a constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * @returns a reference to the \p Mesh.
   */
  const Mesh & mesh () const { return _mesh; }

  
protected:

  /**
   * @returns a writeable reference to the system we are solving.
   */
  sys_type & system () { return _system; }
  
  /**
   * @returns a reference to the \p Mesh.
   */
  Mesh & mesh () { return _mesh; }

  /**
   * A reference to the system we are solving.
   */
  sys_type& _system;

  /**
   * A reference to the \p Mesh for the system
   * we are solving.
   */
  Mesh& _mesh;
};



// ------------------------------------------------------------
// Solver inline members
Solver::Solver (EquationSystems& es) :
  _system (es),
  _mesh   (es.get_mesh())
{
}



Solver::~Solver ()
{
}



void Solver::init ()
{
  std::cout << "Initializing $Id: solver.h,v 1.1 2004-01-03 15:37:42 benkirk Exp $"
	    << std::endl;
  
  // Initialize the system.
  this->system().init ();
}



void Solver::pre_process ()
{
  std::cout << "Pre-processing"
	    << std::endl;
}



void Solver::solve ()
{
  // Perform any necessary pre-processing
  Solver::pre_process ();
  
  std::cout << "Solving $Id: solver.h,v 1.1 2004-01-03 15:37:42 benkirk Exp $"
	    << std::endl;

  // Solve the system
  this->system().solve ();
  
  // Perform any necessary post-processing
  Solver::post_process ();  
}



void Solver::post_process ()
{
  std::cout << "Post-processing"
	    << std::endl;
}


#endif // #define __solver_h__
