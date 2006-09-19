// $Id: discontinuity_measure.h,v 1.3 2006-09-19 17:50:52 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2006  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __discontinuity_measure_h__
#define __discontinuity_measure_h__

// C++ includes
#include <vector>
#include <string>

// Local Includes
#include "error_estimator.h"

// Forward Declarations
class Point;





/**
 * This class measures discontinuities between elements
 * for debugging purposes.  It derives from ErrorEstimator
 * just in case someone finds it useful in a DG framework.
 *
 * @author Roy H. Stogner, 2006
 */
class DiscontinuityMeasure : public ErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to NULL.
   */
  DiscontinuityMeasure() : scale_by_n_internal_sides(false),
			  _bc_function(NULL) {}
  
  /**
   * Destructor.  
   */
  ~DiscontinuityMeasure() {}

  /**
   * This function uses numerical quadrature to sum up the
   * discontinuities between each cell and its neighbors.
   * The integrated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell,
			       bool estimate_parent_error = false);

  /**
   * Register a user function to use in computing external boundary 
   * jumps.
   * The return value is std::pair<bool, Real>
   */
  void attach_flux_bc_function (std::pair<bool,Real> fptr(const System& system,
							  const Point& p,
							  const std::string& var_name));
  /**
   * This boolean flag allows you to scale the error indicator
   * result for each element by the number of internal sides the element
   * actually has.  This tends to weight more evenly cells which are
   * on the boundaries and thus have fewer contributions to jump terms
   * The value is initialized to false, simply set it to true if you
   * want to use the feature.
   */
  bool scale_by_n_internal_sides;
  
private:
  /**
   * Pointer to function that returns BC information.
   */
  std::pair<bool,Real> (* _bc_function) (const System& system,
					 const Point& p,
					 const std::string& var_name);
};


#endif

