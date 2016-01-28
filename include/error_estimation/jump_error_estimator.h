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



#ifndef LIBMESH_JUMP_ERROR_ESTIMATOR_H
#define LIBMESH_JUMP_ERROR_ESTIMATOR_H

// Local Includes
#include "libmesh/auto_ptr.h"
#include "libmesh/dense_vector.h"
#include "libmesh/error_estimator.h"
#include "libmesh/fem_context.h"

// C++ includes
#include <cstddef>
#include <string>
#include <vector>

namespace libMesh
{

// Forward Declarations
class Point;
class Elem;

/**
 * This abstract base class implements utility functions for error estimators
 * which are based on integrated jumps between elements.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class JumpErrorEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.
   */
  JumpErrorEstimator()
    : ErrorEstimator(),
      scale_by_n_flux_faces(false),
      integrate_boundary_sides(false),
      fine_context(),
      coarse_context(),
      fine_error(0),
      coarse_error(0) {}

  /**
   * Destructor.
   */
  virtual ~JumpErrorEstimator() {}


  /**
   * This function uses the derived class's jump error
   * estimate formula to estimate the error on each cell.
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & error_per_cell,
                               const NumericVector<Number> * solution_vector = libmesh_nullptr,
                               bool estimate_parent_error = false) libmesh_override;

  /**
   * This boolean flag allows you to scale the error indicator
   * result for each element by the number of "flux faces" the element
   * actually has.  This tends to weight more evenly cells which are
   * on the boundaries and thus have fewer contributions to their flux.
   * The value is initialized to false, simply set it to true if you
   * want to use the feature.
   */
  bool scale_by_n_flux_faces;

protected:
  /**
   * A utility function to reinit the finite element data on elements sharing a
   * side
   */
  void reinit_sides();

  /**
   * A utility function to correctly increase n_flux_faces for the coarse element
   */
  float coarse_n_flux_faces_increment();

  /**
   * An initialization function, to give derived classes a chance to
   * request specific data from the FE objects
   */
  virtual void init_context(FEMContext & c);

  /**
   * The function, to be implemented by derived classes, which calculates an error
   * term on an internal side
   */
  virtual void internal_side_integration() = 0;

  /**
   * The function, to be implemented by derived classes, which calculates an error
   * term on a boundary side
   * Returns true if the flux bc function is in fact defined on the current side
   */
  virtual bool boundary_side_integration() { return false; }

  /**
   * A boolean flag, by default false, to be set to true if integrations
   * with boundary_side_integration() should be performed
   */
  bool integrate_boundary_sides;

  /**
   * Context objects for integrating on the fine and coarse elements
   * sharing a face
   */
  UniquePtr<FEMContext> fine_context, coarse_context;

  /**
   * The fine and coarse error values to be set by each side_integration();
   */
  Real fine_error, coarse_error;

  /**
   * The variable number currently being evaluated
   */
  unsigned int var;
};


} // namespace libMesh

#endif // LIBMESH_JUMP_ERROR_ESTIMATOR_H
