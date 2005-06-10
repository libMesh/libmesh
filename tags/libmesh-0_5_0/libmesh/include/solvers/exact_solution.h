// $Id: exact_solution.h,v 1.6 2005-06-06 14:53:18 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "vector_value.h" // for RealGradient

// Forward Declarations
class Point;



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
  ExactSolution (EquationSystems& es);

  /**
   * Destructor.
   */
  ~ExactSolution() {}

  
  /**
   * Attach function similar to system.h which
   * allows the user to attach an arbitrary function
   * which computes the exact value of the solution
   * at any point. 
   */
  void attach_exact_value ( Number fptr(const Point& p,
					const Parameters& Parameters,
					const std::string& sys_name,
					const std::string& unknown_name));

  /**
   * Attach function similar to system.h which
   * allows the user to attach an arbitrary function
   * which computes the exact derivative of the solution
   * at any point.
   */
  void attach_exact_deriv ( Gradient fptr(const Point& p,
					  const Parameters& parameters,
					  const std::string& sys_name,
					  const std::string& unknown_name));

  /**
   * Computes and stores the error in the solution value e = u-u_h
   * and the gradient grad(e) = grad(u) - grad(u_h).  Does not return
   * any value.  For that you need to call the l2_error() or h1_error()
   * functions respectively.
   */
  void compute_error(const std::string& sys_name,
		     const std::string& unknown_name);
  
  /**
   * This function returns the integrated L2 error for the system
   * sys_name for the unknown unknown_name.  Note that no error computations
   * are actually performed, you must call compute_error() for that.
   */
  Number l2_error(const std::string& sys_name,
		  const std::string& unknown_name);
  
  /**
   * This function computes and returns the H1 (energy) error for the system
   * sys_name for the unknown unknown_name.  Note that no error computations
   * are actually performed, you must call compute_error() for that.
   */
  Number h1_error(const std::string& sys_name,
		  const std::string& unknown_name);
  
private:
  
  /**
   * This function computes the error (in the solution and its first
   * derivative) for a single unknown in a single system.  It is a
   * private function since it is used by the implementation when
   * solving for several unknowns in several systems.
   */
  void _compute_error(const std::string& sys_name,
		      const std::string& unknown_name,
		      std::pair<Number, Number>& error_pair);

  /**
   * This function is responsible for checking the validity of
   * the sys_name and unknown_name inputs, and returning a
   * reference to the proper pair for storing the values.
   */
  std::pair<Number, Number>& _check_inputs(const std::string& sys_name,
					   const std::string& unknown_name);
  
  /**
   * Function pointer to user-provided function which
   * computes the exact value of the solution.
   */
  Number (* _exact_value) (const Point& p,
			   const Parameters& parameters,
			   const std::string& sys_name,
			   const std::string& unknown_name);

  /**
   * Function pointer to user-provided function which
   * computes the exact derivative of the solution.
   */
  Gradient (* _exact_deriv) (const Point& p,
			     const Parameters& parameters,
			     const std::string& sys_name,
			     const std::string& unknown_name);

  /**
   * Data structure which stores the errors:
   * ||e|| = ||u - u_h||
   * ||grad(e)|| = ||grad(u) - grad(u_h)||
   * for each unknown in a single system.
   * The name of the unknown is
   * the key for the map.
   */
  typedef std::map<std::string, std::pair<Number, Number> > SystemErrorMap;

  /**
   * A map of SystemErrorMaps, which contains entries
   * for each system in the EquationSystems object.
   * This is required, since it is possible for two
   * systems to have unknowns with the *same name*.
   */
  std::map<std::string, SystemErrorMap> _errors;
  
  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems& _equation_systems;

  /**
   * Constant reference to the mesh in the EquationSystems object.
   */
  const Mesh& _mesh;
};




#endif
