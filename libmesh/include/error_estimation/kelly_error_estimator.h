// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "jump_error_estimator.h"

namespace libMesh
{

// Forward Declarations
class Point;





/**
 * This class implements the Kelly error indicator
 * which is based on the flux jumps between elements.
 * See the JumpErrorEstimator class for most user APIs
 *
 * Full BibteX reference:
\verbatim
@Article{Kelly83error,
  author = {D.~W.~Kelly and J.~P.~Gago and O.~C.~Zienkiewicz and I.~Babuska},
  title  = {{A posteriori error analysis and adaptive
             processes in the finite element method: Part I Error analysis}},
  journal = {Int. J. Num. Meth. Engng.},
  volume  = {19},
  pages   = {1593--1619},
  year    = {1983}
}
\endverbatim
 *
 * @author Benjamin S. Kirk, 2003.
 */
class KellyErrorEstimator : public JumpErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to NULL.  Defaults to H1 seminorm; changes to system norm
   * are ignored.
   */
  KellyErrorEstimator() : _bc_function(NULL)
  { error_norm = H1_SEMINORM; }

  /**
   * Destructor.
   */
  ~KellyErrorEstimator() {}

  /**
   * Register a user function to use in computing the flux BCs.
   * The return value is std::pair<bool, Real>
   */
  void attach_flux_bc_function (std::pair<bool,Real> fptr(const System& system,
							  const Point& p,
							  const std::string& var_name));

protected:

  /**
   * An initialization function, for requesting specific data from the FE
   * objects
   */
  virtual void initialize(const System& system,
                          ErrorVector& error_per_cell,
                          bool estimate_parent_error);

  /**
   * The function which calculates a normal derivative jump based error
   * term on an internal side
   */
  virtual void internal_side_integration();

  /**
   * The function which calculates a normal derivative jump based error
   * term on a boundary side.
   * Returns true if the flux bc function is in fact defined on the current side.
   */
  virtual bool boundary_side_integration();

  /**
   * A pointer to the current System
   */
  const System *my_system;

  /**
   * Pointer to function that returns BC information.
   */
  std::pair<bool,Real> (* _bc_function) (const System& system,
					 const Point& p,
					 const std::string& var_name);
};


} // namespace libMesh

#endif

