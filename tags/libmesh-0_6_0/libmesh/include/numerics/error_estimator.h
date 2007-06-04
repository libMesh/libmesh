// $Id: error_estimator.h,v 1.17 2007-04-12 17:21:32 roystgnr Exp $

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



#ifndef __error_estimator_h__
#define __error_estimator_h__

// C++ includes
#include <map>
#include <string>
#include <vector>

// Local Includes
#include "libmesh_common.h"

// Forward Declarations
class ErrorVector;
class EquationSystems;
class System;

/**
 * This class holds functions that will estimate the error
 * in a finite element solution on a given mesh.  These error
 * estimates can be useful in their own right, or may be used
 * to guide adaptive mesh refinement.  Note that in general
 * the computed errors are stored as floats rather than doubles
 * since the required precision is low.
 *
 * @author Benjamin S. Kirk, 2003.
 */
class ErrorEstimator
{
public:

  /**
   * Constructor. Empty.
   */
  ErrorEstimator() : _sobolev_order(1) {}
  
  /**
   * Destructor.  
   */
  virtual ~ErrorEstimator() {}


  /**
   * This pure virtual function must be redefined
   * in derived classes to compute the error for each
   * cell and place it in the "error_per_cell" vector.
   */
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell,
			       bool estimate_parent_error = false) = 0;

  /**
   * This virtual function can be redefined
   * in derived classes, but by default computes the sum of
   * the error_per_cell for each system in the equation_systems.
   *
   * Currently this function ignores the component_scale member variable,
   * and uses the function argument component_scales instead.
   *
   * This function is named estimate_errors instead of estimate_error
   * because otherwise C++ can get confused.
   */
  virtual void estimate_errors (const EquationSystems& equation_systems,
			        ErrorVector& error_per_cell,
			        std::map<const System*, std::vector<float> >& component_scales,
			        bool estimate_parent_error = false);
  /**
   * When calculating many error vectors at once, we need a data structure to
   * hold them all
   */
  typedef std::map<std::pair<const System*, unsigned int>, ErrorVector*> ErrorMap;

  /**
   * This virtual function can be redefined
   * in derived classes, but by default it calls estimate_error
   * repeatedly to calculate the requested error vectors.
   *
   * Currently this function ignores the component_scale member variable,
   * because it calculates each error individually, unscaled.
   * 
   * The user selects which errors get computed by filling a map with error
   * vectors: If errors_per_cell[&system][v] exists, it will be filled with the
   * error values in variable \p v of \p system
   */
  virtual void estimate_errors (const EquationSystems& equation_systems,
			        ErrorMap& errors_per_cell,
			        bool estimate_parent_error = false);

  /**
   * When passed an EquationSystems reference, this function
   * first extracts the System named "name" and then calls the
   * function above.
   */
//   virtual void estimate_error (const EquationSystems& es,
// 			       const std::string& name,
// 			       ErrorVector& error_per_cell)
//   { this->estimate_error (es.get_system<SteadySystem>(name), error_per_cell); }

  /**
   * This vector can be used to "scale" certain
   * variables in a system.  If the mask is empty then the error
   * will be calculated for all variables in the system.  If
   * the mask is not empty the error for each component will be
   * scaled by component_scale[c].
   */
  std::vector<float> component_scale;

  /**
   * In prior versions of libMesh, component_scale was not available
   * and component_mask was used instead - a component_mask entry of false
   * corresponds to a component_scale of 0.0, and a component_mask entry of true
   * corresponds to a component_scale of 1.0.  This vector is now deprecated
   * and should not be used by new code - it will be removed in future libMesh
   * releases.
   */
  std::vector<bool> component_mask;

    /**
     * Returns or allows you to set the Sobolev order for error computations
     * e.g. 0 for H^0/L_2 error, 1 for H^1, 2 for H^2
     * This value defaults to 1, and its exact effects depend on the particular
     * ErrorEstimator subclass being used.
     */
    unsigned int & sobolev_order (void)
      { return _sobolev_order; }

protected:

  /**
   * This method takes the local error contributions in
   * \p error_per_cell from each processor and combines
   * them to get the global error vector.
   */
  void reduce_error (std::vector<float>& error_per_cell) const;

  /**
   * This method tests for a non-empty component_mask vector,
   * and if found checks it for validity and then uses it to
   * fill the newer component_scale vector.
   */
  void convert_component_mask_to_scale();

  /**
   * Sobolev order - e.g. 0 for H^0/L_2 error, 1 for H^1, 2 for H^2
   */
  unsigned int _sobolev_order;
};


#endif

