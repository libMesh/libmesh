// $Id: kelly_error_estimator.h,v 1.1 2004-05-20 22:56:52 jwpeterson Exp $

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



#ifndef __kelly_error_estimator_h__
#define __kelly_error_estimator_h__

// C++ includes
#include <vector>
#include <string>

// Local Includes
#include "error_estimator.h"







/**
 * This class implements the Kelly error indicator
 * which is based on the flux jumps between elements.
 *
 * @author Benjamin S. Kirk, 2003.
 */
class KellyErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.
   */
  KellyErrorEstimator() {}
  
  /**
   * Destructor.  
   */
  ~KellyErrorEstimator() {}


  /**
   * This function uses the Kelly Flux Jump error
   * estimate to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const SteadySystem& system,
			       std::vector<float>& error_per_cell);


  // Bring the base class functionality into the name lookup
  // procedure.  This allows for alternative calling formats
  // defined in the base class.  Thanks Wolfgang.
  using ErrorEstimator::estimate_error;
};


#endif

