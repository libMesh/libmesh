// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_KELLY_ERROR_ESTIMATOR_H
#define LIBMESH_KELLY_ERROR_ESTIMATOR_H

// Local Includes
#include "libmesh/jump_error_estimator.h"

// C++ includes
#include <cstddef>
#include <string>
#include <vector>

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
 *
 * \verbatim
 * @Article{Kelly83error,
 * author = {D.~W.~Kelly and J.~P.~Gago and O.~C.~Zienkiewicz and I.~Babuska},
 * title  = {{A posteriori error analysis and adaptive
 *            processes in the finite element method: Part I Error analysis}},
 * journal = {Int. J. Num. Meth. Engng.},
 * volume  = {19},
 * pages   = {1593--1619},
 * year    = {1983}
 * }
 * \endverbatim
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
class KellyErrorEstimator : public JumpErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to libmesh_nullptr.  Defaults to H1 seminorm; changes to system norm
   * are ignored.
   */
  KellyErrorEstimator() :
    JumpErrorEstimator(),
    _bc_function(libmesh_nullptr)
  { error_norm = H1_SEMINORM; }

  /**
   * Destructor.
   */
  ~KellyErrorEstimator() {}

  /**
   * Register a user function to use in computing the flux BCs.
   * The return value is std::pair<bool, Real>
   */
  void attach_flux_bc_function (std::pair<bool,Real> fptr(const System & system,
                                                          const Point & p,
                                                          const std::string & var_name));

  virtual ErrorEstimatorType type() const libmesh_override
  { return KELLY;}

protected:

  /**
   * An initialization function, for requesting specific data from the FE
   * objects
   */
  virtual void init_context(FEMContext & c) libmesh_override;

  /**
   * The function which calculates a normal derivative jump based error
   * term on an internal side
   */
  virtual void internal_side_integration() libmesh_override;

  /**
   * The function which calculates a normal derivative jump based error
   * term on a boundary side.
   * Returns true if the flux bc function is in fact defined on the current side.
   */
  virtual bool boundary_side_integration() libmesh_override;

  /**
   * Pointer to function that returns BC information.
   */
  std::pair<bool,Real> (* _bc_function) (const System & system,
                                         const Point & p,
                                         const std::string & var_name);
};


} // namespace libMesh

#endif // LIBMESH_KELLY_ERROR_ESTIMATOR_H
