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

// rbOOmit includes
#include "libmesh/rb_temporal_discretization.h"

// libMesh includes
#include "libmesh/getpot.h"

namespace libMesh
{

RBTemporalDiscretization::RBTemporalDiscretization()
  : _delta_t(0.),
    _euler_theta(0.),
    _current_time_step(0),
    _n_time_steps(0)
{}

Real RBTemporalDiscretization::get_delta_t() const
{
  return _delta_t;
}

void RBTemporalDiscretization::set_delta_t(const Real delta_t_in)
{
  _delta_t = delta_t_in;
}

Real RBTemporalDiscretization::get_euler_theta() const
{
  return _euler_theta;
}

void RBTemporalDiscretization::set_euler_theta(const Real euler_theta_in)
{
  libmesh_assert((0. <= euler_theta_in ) && (euler_theta_in <= 1.));
  _euler_theta = euler_theta_in;
}

unsigned int RBTemporalDiscretization::get_time_step() const
{
  return _current_time_step;
}

void RBTemporalDiscretization::set_time_step(const unsigned int k)
{
  libmesh_assert_less_equal (k, get_n_time_steps());
  this->_current_time_step = k;
}

unsigned int RBTemporalDiscretization::get_n_time_steps() const
{
  return _n_time_steps;
}

void RBTemporalDiscretization::set_n_time_steps(const unsigned int K)
{
  _n_time_steps = K;
  _control.assign(_n_time_steps+1,1.0);
}

Real RBTemporalDiscretization::get_control(const unsigned int k) const
{
  libmesh_assert_less_equal (k, get_n_time_steps());
  return _control[k];
}

void RBTemporalDiscretization::set_control(const std::vector<Real> & control)
{
  libmesh_assert_less_equal(control.size(),_n_time_steps+1);
  _control = control;
  // If the input vector is smaller than the number of time steps (+1), we complete it with zeros
  _control.resize(_n_time_steps+1);
}

void RBTemporalDiscretization::process_temporal_parameters_file (const std::string & parameters_filename)
{
  // Read in data from parameters_filename
  GetPot infile(parameters_filename);

  // Read in parameters related to temporal discretization
  unsigned int n_time_steps_in = infile("n_time_steps", get_n_time_steps());
  const Real delta_t_in        = infile("delta_t",      get_delta_t());
  const Real euler_theta_in    = infile("euler_theta",  get_euler_theta());

  // and set the relevant member variables
  set_n_time_steps(n_time_steps_in);
  set_delta_t(delta_t_in);
  set_euler_theta(euler_theta_in);
  set_time_step(0);
}

void RBTemporalDiscretization::pull_temporal_discretization_data(RBTemporalDiscretization & other)
{
  this->set_delta_t( other.get_delta_t() );
  this->set_euler_theta( other.get_euler_theta() );
  this->set_n_time_steps( other.get_n_time_steps() );
  this->set_time_step( other.get_time_step() );
  this->set_control( other._control );
}

}
