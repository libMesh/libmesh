// $Id: error_estimator.h,v 1.7 2003-09-25 21:46:55 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
   * Destructor.  Virtual, so feel free to derive your
   * own error estimator from this if you want.
   */
  virtual ~ErrorEstimator() {}

  
  /**
   * This function uses the Kelley Flux Jump error
   * estimate to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void flux_jump (const EquationSystems& es,
			  const std::string& name,
			  std::vector<float>& error_per_cell);


  /**
   * This vector can be used to "mask" certain
   * variables in a system.  If the mask is empty then the error
   * will be calculated for all variables in the system.  If
   * the mask is not empty the error will be calculated only
   * for those components for which \p component_mask[c] is
   * true.
   */
  std::vector<bool> component_mask;
};


#endif

