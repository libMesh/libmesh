// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DIFF_QOI_H
#define LIBMESH_DIFF_QOI_H

// Local Includes
#include "libmesh/diff_context.h"

// C++ includes
#include <memory>

namespace libMesh
{

// Forward declarations
class DiffContext;
class QoISet;

namespace Parallel {
  class Communicator;
}

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class DifferentiableQoI
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DifferentiableQoI ();

  /**
   * Destructor.
   */
  virtual ~DifferentiableQoI () = default;

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * Initialize system qoi.  This version of the function required
   * direct vector access, and is now deprecated.
   */
  virtual void init_qoi( std::vector<Number> & /*sys_qoi*/){}
#else
  /**
   * Non-virtual, to try to help deprecated user code catch this
   * change at compile time (if they specified override)
   */
  void init_qoi( std::vector<Number> & /*sys_qoi*/){}
#endif

  /**
   * Initialize system qoi.  Often this will just call
   * sys.init_qois(some_desired_number_of_qois)
   */
  virtual void init_qoi_count( System & /*sys*/){}

  /**
   * Clear all the data structures associated with
   * the QoI.
   */
  virtual void clear_qoi () {}

  /**
   * If \p assemble_qoi_sides is true (it is false by default), the
   * assembly loop for a quantity of interest or its derivatives will
   * loop over domain boundary sides.  To add domain interior sides,
   * also set assemble_qoi_internal_sides to true.
   */
  bool assemble_qoi_sides;

  /**
   * If \p assemble_qoi_internal_sides is true (it is false by
   * default), the assembly loop for a quantity of interest or its
   * derivatives will loop over element sides which do not fall on
   * domain boundaries.
   */
  bool assemble_qoi_internal_sides;

  /**
   * If \p assemble_qoi_elements is false (it is true by default), the
   * assembly loop for a quantity of interest or its derivatives will
   * skip computing on mesh elements, and will only compute on mesh
   * sides.
   */
  bool assemble_qoi_elements;

  /**
   * Does any work that needs to be done on \p elem in a quantity of
   * interest assembly loop, outputting to elem_qoi.
   *
   * Only qois included in the supplied \p QoISet need to be
   * assembled.
   */
  virtual void element_qoi (DiffContext &,
                            const QoISet &)
  {}

  /**
   * Does any work that needs to be done on \p elem in a quantity of
   * interest derivative assembly loop, outputting to
   * elem_qoi_derivative
   *
   * Only qois included in the supplied \p QoISet need their
   * derivatives assembled.
   */
  virtual void element_qoi_derivative (DiffContext &,
                                       const QoISet &)
  {}

  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * quantity of interest assembly loop, outputting to elem_qoi.
   *
   * Only qois included in the supplied \p QoISet need to be
   * assembled.
   */
  virtual void side_qoi (DiffContext &,
                         const QoISet &)
  {}

  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * quantity of interest derivative assembly loop, outputting to
   * elem_qoi_derivative.
   *
   * Only qois included in the supplied \p QoISet need their
   * derivatives assembled.
   */
  virtual void side_qoi_derivative (DiffContext &,
                                    const QoISet &)
  {}

  /**
   * Prepares the result of a build_context() call for use.
   *
   * FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular QoI requires.  Trying to
   * evaluate a QoI without overriding init_context is both
   * inefficient and deprecated.
   */
  virtual void init_context(DiffContext &) { libmesh_deprecated(); }

  /**
   * Copy of this object. User should override to copy any needed state.
   */
  virtual std::unique_ptr<DifferentiableQoI> clone() =0;

  /**
   * Method to combine thread-local qois. By default, simply sums thread qois.
   */
  virtual void thread_join(std::vector<Number> & qoi,
                           const std::vector<Number> & other_qoi,
                           const QoISet & qoi_indices);

  /**
   * Method to populate system qoi data structure with process-local qoi. By default, simply
   * sums process qois into system qoi.
   */
  virtual void parallel_op(const Parallel::Communicator & communicator,
                           std::vector<Number> & sys_qoi,
                           std::vector<Number> & local_qoi,
                           const QoISet & qoi_indices);

  /**
   * Method to finalize qoi derivatives which require more than just a simple
   * sum of element contributions.
   */
  virtual void finalize_derivative(NumericVector<Number> & derivatives, std::size_t qoi_index);
};

} // namespace libMesh


#endif // LIBMESH_DIFF_QOI_H
