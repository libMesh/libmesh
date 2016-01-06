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

#ifndef LIBMESH_RB_TEMPORAL_DISCRETIZATION_H
#define LIBMESH_RB_TEMPORAL_DISCRETIZATION_H

// libMesh includes
#include "libmesh/libmesh_common.h"

#include <vector>

namespace libMesh
{

/**
 * Define a class that encapsulates the details of a
 * "generalized Euler" temporal discretization to be
 * used in the rbOOmit framework.
 */
class RBTemporalDiscretization
{
public:

  /**
   * Constructor.
   */
  RBTemporalDiscretization();

  /**
   * Get/set delta_t, the time-step size.
   */
  Real get_delta_t() const;
  void set_delta_t(const Real delta_t_in);

  /**
   * Get/set euler_theta, parameter that determines
   * the temporal discretization.
   */
  Real get_euler_theta() const;
  void set_euler_theta(const Real euler_theta_in);

  /**
   * Get/set the current time-step.
   */
  unsigned int get_time_step() const;
  void set_time_step(const unsigned int k);

  /**
   * Get/set the total number of time-steps.
   */
  unsigned int get_n_time_steps() const;
  void set_n_time_steps(const unsigned int K);

  /**
   * Get/set the RHS control.
   */
  Real get_control(const unsigned int k) const;
  void set_control(const std::vector<Real> & control);

  /**
   * Read in and initialize parameters from \p parameters_filename.
   */
  void process_temporal_parameters_file (const std::string & parameters_filename);

  /**
   * Pull the temporal discretization data from \p other.
   */
  void pull_temporal_discretization_data(RBTemporalDiscretization & other);

private:

  /**
   * The time-step size.
   */
  Real _delta_t;

  /**
   * The parameter that determines the generalized Euler scheme
   * discretization that we employ.
   * euler_theta = 0   ---> Forward Euler
   * euler_theta = 0.5 ---> Crank-Nicolson
   * euler_theta = 1   ---> Backward Euler
   */
  Real _euler_theta;

  /**
   * The current time-step.
   */
  unsigned int _current_time_step;

  /**
   * The number of time-steps.
   */
  unsigned int _n_time_steps;

  /**
   * The RHS control (scalar function of time).
   * A function h(t) that is used in the RHS as h(t)*f(x,\f$ \mu \f$).
   * See Martin Grepl's thesis
   */
  std::vector<Real> _control;

};

}

#endif // LIBMESH_RB_TEMPORAL_DISCRETIZATION_H
