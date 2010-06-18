// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef __qn_transient_rb_context_h__
#define __qn_transient_rb_context_h__

// Configuration data
#include "libmesh_config.h"

// Depends on the QNTransientRBSystem, which requires SLEPC
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

// Local Includes
#include "rb_context.h"

// Forward declaration
class QNTransientRBSystem;
class FEBase;
class QBase;

/**
 * This class is part of the rbOOmit framework.
 *
 * QNTransientRBContext provides extra context
 * information relevant to QNTransient problems.
 *
 * @author David J. Knezevic, 2009
 */

// ------------------------------------------------------------
// QNTransientRBContext class definition

class QNTransientRBContext : public RBContext
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  QNTransientRBContext (const QNTransientRBSystem &);

  /**
   * Destructor.
   */
  virtual ~QNTransientRBContext ();

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior
   */
  Number old_interior_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element side
   */
  Number old_side_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the solution variable \p var at the physical
   * point \p p on the current element
   */
  Number old_point_value(unsigned int var, Point &p);

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior
   */
  Gradient old_interior_gradient(unsigned int var, unsigned int qp);

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side
   */
  Gradient old_side_gradient(unsigned int var, unsigned int qp);

  /**
   * Reinitialize all the context data on a given
   * element for the given system. Overload to
   * also reinitialize the solution data from
   * the previous time step.
   */
  virtual void reinit(RBSystem&, Elem*);

  /**
   * Element by element components of nonlinear_solution
   * as adjusted by a time_solver
   */
  DenseVector<Number> elem_old_solution;
  std::vector<DenseSubVector<Number> *> elem_old_subsolutions;
};

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
