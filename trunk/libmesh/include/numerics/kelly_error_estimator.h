// $Id: kelly_error_estimator.h,v 1.9 2006-09-19 15:46:34 roystgnr Exp $

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



#ifndef __kelly_error_estimator_h__
#define __kelly_error_estimator_h__

// C++ includes
#include <vector>
#include <string>

// Local Includes
#include "error_estimator.h"

// Forward Declarations
class Point;





/**
 * This class implements the Kelly error indicator
 * which is based on the flux jumps between elements.
 *
 * Full BibteX reference:
 * 
 * @Article{Kelly83error,
 * author = {D.~W.~Kelly and J.~P.~Gago and O.~C.~Zienkiewicz and I.~Babuska},
 * title  = {{A posteriori error analysis and adaptive
 *            processes in the finite element method: Part I Error analysis}},
 * journal = {Int. J. Num. Meth. Engng.},
 * volume  = {19},
 * pages   = {1593--1619},
 * year    = {1983}
 * }
 *
 * @author Benjamin S. Kirk, 2003.
 */
class KellyErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to NULL.
   */
  KellyErrorEstimator() : scale_by_n_flux_faces(false),
			  _bc_function(NULL) {}
  
  /**
   * Destructor.  
   */
  ~KellyErrorEstimator() {}


  // Bring the base class functionality into the name lookup
  // procedure.  This allows for alternative calling formats
  // defined in the base class.  Thanks Wolfgang.
  // GCC 2.95.3 cannot compile such code.  Since it was not really
  // essential to the functioning of this class, it's been removed.
  // using ErrorEstimator::estimate_error;

  /**
   * This function uses the Kelly Flux Jump error
   * estimate to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell);

  /**
   * Register a user function to use in computing the flux BCs.
   * The return value is std::pair<bool, Real>
   */
  void attach_flux_bc_function (std::pair<bool,Real> fptr(const System& system,
							  const Point& p,
							  const std::string& var_name));
  /**
   * This boolean flag allows you to scale the error indicator
   * result for each element by the number of "flux faces" the element
   * actually has.  This tends to weight more evenly cells which are
   * on the boundaries and thus have fewer contributions to their flux.
   * The value is initialized to false, simply set it to true if you
   * want to use the feature.
   */
  bool scale_by_n_flux_faces;
  
private:
  /**
   * Pointer to function that returns BC information.
   */
  std::pair<bool,Real> (* _bc_function) (const System& system,
					 const Point& p,
					 const std::string& var_name);
};


#endif

