// $Id: error_estimator.h,v 1.10 2006-07-25 17:59:42 roystgnr Exp $

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
#include <vector>
#include <string>

// Local Includes
#include "libmesh_common.h"

// Forward Declarations
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
  ErrorEstimator() {}
  
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
			       std::vector<float>& error_per_cell) = 0;

  /**
   * When passed an EquationSystems reference, this function
   * first extracts the System named "name" and then calls the
   * function above.
   */
//   virtual void estimate_error (const EquationSystems& es,
// 			       const std::string& name,
// 			       std::vector<float>& error_per_cell)
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
};


#endif

