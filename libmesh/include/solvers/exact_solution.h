// $Id$

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

#ifndef __exact_solution_h__
#define __exact_solution_h__


// C++ includes
#include <map>
#include <vector>

// Local Includes
#include "libmesh_common.h" // for Number


// Forward Declarations
class Point;
class EquationSystems;
class Parameters;
class Mesh;

// Is there any way to simplify this?
// All we need are Tensor and Gradient. - RHS
template <typename T> class TensorValue;
template <typename T> class VectorValue;
typedef TensorValue<Number> NumberTensorValue;
typedef NumberTensorValue   Tensor;
typedef VectorValue<Number> NumberVectorValue;
typedef NumberVectorValue   Gradient;

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
   * allows the user to attach a second EquationSystems
   * object with a reference fine grid solution.
   */
  void attach_reference_solution (EquationSystems* es_fine);

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
   * Attach function similar to system.h which
   * allows the user to attach an arbitrary function
   * which computes the exact second derivatives of the solution
   * at any point.
   */
  void attach_exact_hessian ( Tensor fptr(const Point& p,
					  const Parameters& parameters,
					  const std::string& sys_name,
					  const std::string& unknown_name));

  /**
   * Increases or decreases the order of the quadrature rule used for numerical
   * integration.
   */
  void extra_quadrature_order (const int extraorder)
    { _extra_order = extraorder; }

  /**
   * Computes and stores the error in the solution value e = u-u_h,
   * the gradient grad(e) = grad(u) - grad(u_h), and possibly the hessian
   * grad(grad(e)) = grad(grad(u)) - grad(grad(u_h)).  Does not return
   * any value.  For that you need to call the l2_error(), h1_error()
   * or h2_error() functions respectively.
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
   * This function computes and returns the H1 error for the system
   * sys_name for the unknown unknown_name.  Note that no error computations
   * are actually performed, you must call compute_error() for that.
   */
  Number h1_error(const std::string& sys_name,
		const std::string& unknown_name);
  
  /**
   * This function computes and returns the H2 error for the system
   * sys_name for the unknown unknown_name.  Note that no error computations
   * are actually performed, you must call compute_error() for that.
   */
  Number h2_error(const std::string& sys_name,
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
		      std::vector<Number>& error_vals);

  /**
   * This function is responsible for checking the validity of
   * the sys_name and unknown_name inputs, and returning a
   * reference to the proper vector for storing the values.
   */
  std::vector<Number>& _check_inputs(const std::string& sys_name,
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
   * Function pointer to user-provided function which
   * computes the exact hessian of the solution.
   */
  Tensor (* _exact_hessian) (const Point& p,
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
  typedef std::map<std::string, std::vector<Number> > SystemErrorMap;

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
   * Constant pointer to the \p EquationSystems object
   * containing the fine grid solution.
   */
  EquationSystems* _equation_systems_fine;

  /**
   * Extra order to use for quadrature rule
   */
  int _extra_order;
};




#endif
