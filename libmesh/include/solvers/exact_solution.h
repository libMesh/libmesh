// $Id: exact_solution.h,v 1.1 2004-05-24 19:58:39 jwpeterson Exp $

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

#ifndef __exact_solution_h__
#define __exact_solution_h__


// C++ includes

// Local Includes
#include "equation_systems.h"

// Forward Declarations




/**
 * This class handles the computation of the L2 and/or H1
 * error for the Systems in the EquationSystems object
 * which is passed to it.  Note that for it to be useful,
 * the user must attach at least one, and possibly two functions
 * which can compute the exact solution and its derivative
 * for any component of any system.  These are the exact_value
 * and exact_deriv functions below.
 *
 * @author Benjamin S. Kirk w/ modifications for libmesh
 * by John W. Peterson
 */


class ExactSolution
{
  
public:
  /**
   * Constructor.  The ExactSolution object
   * must be initialized with an EquationSystems
   * object.
   */
  ExactSolution (EquationSystems& es) :
    _exact_value (NULL),
    _exact_deriv (NULL),
    _equation_systems(es)
    {}

  /**
   * Destructor.
   */
  ~ExactSolution() {}

  
  /**
   * Attach function similar to system.h which
   * allows the user to attach an arbitrary function
   * which computes the exact value of the solution
   * at any point, time.
   */
  void attach_exact_value ( Real fptr(const Point& p,
				      const Real time,
				      const std::string& sys_name,
				      const std::string& unknown_name));

  /**
   * Attach function similar to system.h which
   * allows the user to attach an arbitrary function
   * which computes the exact derivative of the solution
   * at any point, time.
   */
  void attach_exact_deriv ( Point fptr(const Point& p,
				       const Real time,
				       const std::string& sys_name,
				       const std::string& unknown_name));
				      
  
private:
  /**
   * This function computes the error for a single unknown in a single
   * system.  It is a private function since it is used by the implementation
   * when solving for several unknowns in several systems.
   */
  
  
  /**
   * Function pointer to user-provided function which
   * computes the exact value of the solution.
   */
  Real (* _exact_value) (const Point& p,
			 const Real   time,
			 const std::string& sys_name,
			 const std::string& unknown_name);

  /**
   * Function pointer to user-provided function which
   * computes the exact derivative of the solution.
   */
  Point (* _exact_deriv) (const Point& p,
			  const Real   time,
			  const std::string& sys_name,
			  const std::string& unknown_name);
  
  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems& _equation_systems;
};




#endif
