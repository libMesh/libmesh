
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



#ifndef LIBMESH_DIFF_QOI_H
#define LIBMESH_DIFF_QOI_H

// Local Includes
#include "libmesh/auto_ptr.h"
#include "libmesh/parallel.h"

// C++ includes

namespace libMesh
{

// Forward declarations
class DiffContext;
class QoISet;

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
  virtual ~DifferentiableQoI () {}

  /**
   * Initialize system qoi. By default, does nothing in order to maintain backward
   * compatibility for FEMSystem applications that control qoi.
   */
  virtual void init_qoi( std::vector<Number> & /*sys_qoi*/){}

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
   * Most FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular QoI requires.
   */
  virtual void init_context(DiffContext &) {}

  /**
   * Copy of this object. User should override to copy any needed state.
   */
  virtual UniquePtr<DifferentiableQoI> clone() =0;

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
};

} // namespace libMesh


#endif // LIBMESH_DIFF_QOI_H
