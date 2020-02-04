// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DISCONTINUITY_MEASURE_H
#define LIBMESH_DISCONTINUITY_MEASURE_H

// Local Includes
#include "libmesh/jump_error_estimator.h"

// C++ includes
#include <cstddef>
#include <string>
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename> class PointTempl;
typedef PointTempl<Real> Point;





/**
 * This class measures discontinuities between elements
 * for debugging purposes.  It derives from ErrorEstimator
 * just in case someone finds it useful in a DG framework.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class DiscontinuityMeasure : public JumpErrorEstimator
{
public:

  /**
   * Constructor.  Responsible for initializing the _bc_function function
   * pointer to nullptr.  Defaults to L2 norm; changes to system norm are
   * ignored.
   */
  DiscontinuityMeasure();

  /**
   * This class cannot be (default) copy constructed/assigned because
   * its base class has unique_ptr members.
   */
  DiscontinuityMeasure (const DiscontinuityMeasure &) = delete;
  DiscontinuityMeasure & operator= (const DiscontinuityMeasure &) = delete;

  /**
   * Defaulted move ctor, move assignment operator, and destructor.
   */
  DiscontinuityMeasure (DiscontinuityMeasure &&) = default;
  DiscontinuityMeasure & operator= (DiscontinuityMeasure &&) = default;
  virtual ~DiscontinuityMeasure() = default;

  /**
   * Register a user function to use in computing the essential BCs.
   */
  void attach_essential_bc_function (std::pair<bool,Real> fptr(const System & system,
                                                               const Point & p,
                                                               const std::string & var_name));

  virtual ErrorEstimatorType type() const override;

protected:

  /**
   * An initialization function, for requesting specific data from the FE
   * objects
   */
  virtual void init_context(FEMContext & c) override;

  /**
   * The function which calculates a normal derivative jump based error
   * term on an internal side
   */
  virtual void internal_side_integration() override;

  /**
   * The function which calculates a normal derivative jump based error
   * term on a boundary side.
   *
   * \returns \p true if the flux bc function is in fact defined on the current side.
   */
  virtual bool boundary_side_integration() override;

  /**
   * Pointer to function that provides BC information.
   */
  std::pair<bool,Real> (* _bc_function) (const System & system,
                                         const Point & p,
                                         const std::string & var_name);
};


} // namespace libMesh

#endif // LIBMESH_DISCONTINUITY_MEASURE_H
