// $Id: patch_recovery_error_estimator.h,v 1.1 2004-06-02 20:32:04 benkirk Exp $

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



#ifndef __patch_recovery_error_estimator_h__
#define __patch_recovery_error_estimator_h__

// C++ includes
#include <vector>
#include <string>
#include <set>

// Local Includes
#include "error_estimator.h"







/**
 * This class implements the Patch Recovery error indicator.
 * 
 *
 * @author Varis Carey, Benjamin S. Kirk, 2004.
 */
class PatchRecoveryErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.
   */
  PatchRecoveryErrorEstimator() {}
  
  /**
   * Destructor.  
   */
  ~PatchRecoveryErrorEstimator() {}


  /**
   * This function uses the Patch Recovery error
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

private:

  /**
   * Builds a patch containing element \p elem by repeated addition
   * of neighbors.  This procedure is repeated until the \p target_patch_size
   * is met or exceeded.  
   */
  void build_patch_from_local_neighbors (const Elem* elem,
					 std::set<const Elem*> patch,
					 const unsigned int target_patch_size = 10);
};


#endif
