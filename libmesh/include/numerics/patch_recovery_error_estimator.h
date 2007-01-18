// $Id: patch_recovery_error_estimator.h,v 1.11 2007-01-18 22:39:55 roystgnr Exp $

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



#ifndef __patch_recovery_error_estimator_h__
#define __patch_recovery_error_estimator_h__

// C++ includes
#include <vector>

// Local Includes
#include "error_estimator.h"
#include "enum_order.h"
#include "point.h"

// Forward Declarations
class Elem;



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
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell,
			       bool estimate_parent_error = false);


  // Bring the base class functionality into the name lookup
  // procedure.  This allows for alternative calling formats
  // defined in the base class.  Thanks Wolfgang.
  // GCC 2.95.3 cannot compile such code.  Since it was not really
  // essential to the functioning of this class, it's been removed.
  // using ErrorEstimator::estimate_error;

private:

  /**
   * Computes the factorial of n-used for determining matrix size
   */  
  unsigned int factorial(const unsigned int n);

  /**
   * Returns the spectral polynomial basis function values at a point x,y,z
   */
  
  std::vector<Real> specpoly(const unsigned int dim,
			     const Order order,
			     const Point p,
			     const unsigned int matsize);
};



//-----------------------------------------------------------------
// PatchRecoveryErrorEstimator implementations
inline
unsigned int PatchRecoveryErrorEstimator::factorial (const unsigned int n)
{
  unsigned int fac=1;
  for (int i=n; i>0; i--) fac *= i;
  return fac;
}



#endif // #define __patch_recovery_error_estimator_h__
