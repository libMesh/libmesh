// $Id: fourth_error_estimators.h,v 1.1 2005-05-27 14:45:37 roystgnr Exp $

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



#ifndef __fourth_error_estimators_h__
#define __fourth_error_estimators_h__

// C++ includes
#include <vector>
#include <string>

// Local Includes
#include "error_estimator.h"







/**
 * This class is an error indicator based on laplacian jumps between
 * elements.
 *
 * @author Roy H. Stogner, 2005
 */
class LaplacianErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.
   */
  LaplacianErrorEstimator() {}
  
  /**
   * Destructor.  
   */
  ~LaplacianErrorEstimator() {}


  /**
   * This function uses the Laplacian jump error
   * indicator to estimate the error on each cell.
   * The output is in the vector \p error_per_cell
   */
  virtual void estimate_error (const System& system,
			       std::vector<float>& error_per_cell);


};


#endif

