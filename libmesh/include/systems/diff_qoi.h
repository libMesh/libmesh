
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



#ifndef __diff_qoi_h__
#define __diff_qoi_h__

// C++ includes

// Local Includes
#include "diff_context.h"
#include "qoi_set.h"

namespace libMesh
{

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// DifferentiableQoI class definition

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
   * Clear all the data structures associated with
   * the QoI. 
   */
  virtual void clear_qoi () {}

  /**
   * Executes a postprocessing loop over all elements, and if
   * \p postprocess_sides is true over all sides.
   */
  virtual void postprocess () = 0;

  /**
   * If \p postprocess_sides is true (it is false by default), the
   * postprocessing loop will loop over all sides as well as all
   * elements.
   */
  bool postprocess_sides;

  /**
   * If \p assemble_qoi_sides is true (it is false by default), the
   * assembly loop for a quantity of interest or its derivatives will
   * loop over all sides as well as all elements.
   */
  bool assemble_qoi_sides;

  /**
   * Does any work that needs to be done on \p elem in a postprocessing loop.
   */
  virtual void element_postprocess (DiffContext &) {}
 
  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * postprocessing loop.
   */
  virtual void side_postprocess (DiffContext &) {}
 
  /**
   * Does any work that needs to be done on \p elem in a quantity of
   * interest assembly loop, outputting to elem_qoi.
   *
   * Only qois included in the supplied \p QoISet need to be
   * assembled.
   */
  virtual void element_qoi (DiffContext&, 
                            const QoISet&)
    {}
 
  /**
   * Does any work that needs to be done on \p elem in a quantity of
   * interest derivative assembly loop, outputting to
   * elem_qoi_derivative
   *
   * Only qois included in the supplied \p QoISet need their
   * derivatives assembled.
   */
  virtual void element_qoi_derivative (DiffContext&,
                                       const QoISet&)
    {}
 
  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * quantity of interest assembly loop, outputting to elem_qoi.
   *
   * Only qois included in the supplied \p QoISet need to be
   * assembled.
   */
  virtual void side_qoi (DiffContext&,
                         const QoISet&)
    {}
 
  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * quantity of interest derivative assembly loop, outputting to
   * elem_qoi_derivative.
   *
   * Only qois included in the supplied \p QoISet need their
   * derivatives assembled.
   */
  virtual void side_qoi_derivative (DiffContext&,
                                    const QoISet&)
    {}
 
  /*
   * Prepares the result of a build_context() call for use.
   * 
   * Most FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular QoI requires.
   */
  virtual void init_context(DiffContext &) {}
};

} // namespace libMesh


#endif
