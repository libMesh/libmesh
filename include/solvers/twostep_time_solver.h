// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TWOSTEP_TIME_SOLVER_H
#define LIBMESH_TWOSTEP_TIME_SOLVER_H

// Local includes
#include "libmesh/adaptive_time_solver.h"

// C++ includes

namespace libMesh
{

// Forward declarations
class System;

// UPDATE THIS DESCRIPTION

/**
 * This class wraps another UnsteadySolver derived class, and compares
 * the results of timestepping with deltat and timestepping with
 * 2*deltat to adjust future timestep lengths.
 *
 * Currently this class only works on fully coupled Systems
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2007
 */
class TwostepTimeSolver : public AdaptiveTimeSolver
{
public:
  /**
   * The parent class
   */
  typedef AdaptiveTimeSolver Parent;

  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  TwostepTimeSolver (sys_type & s);

  /**
   * Destructor.
   */
  ~TwostepTimeSolver ();

  virtual void solve() override;

  virtual std::pair<unsigned int, Real> adjoint_solve (const QoISet & qoi_indices) override;

  /**
   * A method to integrate the adjoint sensitivity w.r.t a given parameter
   * vector. int_{tstep_start}^{tstep_end} dQ/dp dt = int_{tstep_start}^{tstep_end} (\partialQ / \partial p) - ( \partial R (u,z) / \partial p ) dt
   * The midpoint rule is used to integrate each substep
   */
  virtual void integrate_adjoint_sensitivity(const QoISet & qois, const ParameterVector & parameter_vector, SensitivityData & sensitivities) override;

};


} // namespace libMesh


#endif // LIBMESH_TWOSTEP_TIME_SOLVER_H
