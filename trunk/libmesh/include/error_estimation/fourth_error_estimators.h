// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __fourth_error_estimators_h__
#define __fourth_error_estimators_h__

// C++ includes
#include <vector>
#include <string>

// Local Includes
#include "jump_error_estimator.h"

namespace libMesh
{







/**
 * This class is an error indicator based on laplacian jumps between
 * elements.
 * See the JumpErrorEstimator class for most user APIs
 *
 * @author Roy H. Stogner, 2005
 */
class LaplacianErrorEstimator : public JumpErrorEstimator
{
public:

  /**
   * Constructor.  Defaults to H2 seminorm; changes to error_norm are
   * ignored.
   */
  LaplacianErrorEstimator() { error_norm = H2_SEMINORM; }

  /**
   * Destructor.
   */
  ~LaplacianErrorEstimator() {}

protected:

  /**
   * An initialization function, for requesting specific data from the FE
   * objects
   */
  virtual void initialize(const System& system,
                          ErrorVector& error_per_cell,
                          bool estimate_parent_error);

  /**
   * The function which calculates a laplacian jump based error
   * term on an internal side
   */
  virtual void internal_side_integration();
};


} // namespace libMesh

#endif

