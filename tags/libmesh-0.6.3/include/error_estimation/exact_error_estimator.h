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



#ifndef __exact_error_estimator_h__
#define __exact_error_estimator_h__

// C++ includes
#include <vector>
#include <string>

// Local Includes
#include "error_estimator.h"

// Forward Declarations
class Elem;
class FEBase;
class MeshFunction;
class Point;
class Parameters;

template <typename T> class DenseVector;

// Is there any way to simplify this?
// All we need are Tensor and Gradient. - RHS
template <typename T> class TensorValue;
template <typename T> class VectorValue;
typedef TensorValue<Number> NumberTensorValue;
typedef NumberTensorValue   Tensor;
typedef VectorValue<Number> NumberVectorValue;
typedef NumberVectorValue   Gradient;



/**
 * This class implements an "error estimator"
 * based on the difference between the approximate
 * and exact solution.  In theory the quadrature error
 * in this estimate should be much lower than the
 * approximation error in other estimates, so this
 * estimator can be used to calculate effectivity.
 *
 * @author Roy Stogner, 2006.
 */
class ExactErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to NULL.
   */
  ExactErrorEstimator() : _exact_value(NULL), 
                          _exact_deriv(NULL),
                          _exact_hessian(NULL),
			  _extra_order(0)
  {}
  
  /**
   * Destructor.  
   */
  ~ExactErrorEstimator() {}


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
   * Attach function similar to system.h which
   * allows the user to attach a second EquationSystems
   * object with a reference fine grid solution.
   */
  void attach_reference_solution (EquationSystems* es_fine);
  

  /**
   * Increases or decreases the order of the quadrature rule used for numerical
   * integration.
   */
  void extra_quadrature_order (const int extraorder)
    { _extra_order = extraorder; }


  // Bring the base class functionality into the name lookup
  // procedure.  This allows for alternative calling formats
  // defined in the base class.  Thanks Wolfgang.
  // GCC 2.95.3 cannot compile such code.  Since it was not really
  // essential to the functioning of this class, it's been removed.
  // using ErrorEstimator::estimate_error;

  /**
   * This function uses the exact solution function
   * to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell,
			       bool estimate_parent_error = false);

private:

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
   * Constant pointer to the \p EquationSystems object
   * containing the fine grid solution.
   */
  EquationSystems* _equation_systems_fine;

  /**
   * Helper method for calculating on each element
   */
  Real find_squared_element_error (const System& system,
				   const std::string& var_name,
				   const Elem *elem,
                                   const DenseVector<Number> &Uelem,
                                   FEBase *fe,
				   MeshFunction *fine_values) const;

  /**
   * Extra order to use for quadrature rule
   */
  int _extra_order;
};


#endif

