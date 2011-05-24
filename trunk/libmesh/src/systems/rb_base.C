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

#include "libmesh_logging.h"
#include "equation_systems.h"
#include "parallel.h"
#include "rb_eim_system.h"
#include "rb_eim_evaluation.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBBase implementation

RBBase::RBBase()
{
  libmesh_experimental();
}

RBBase::~RBBase ()
{
}

unsigned int RBBase::get_n_params() const
{
  libmesh_assert( mu_min_vector.size() == mu_max_vector.size() );

  return mu_min_vector.size();
}

void RBBase::set_parameter_range(std::vector<Real> mu_min_in, std::vector<Real> mu_max_in)
{
  libmesh_assert( mu_min_in.size() == mu_max_in.size() );
  
  mu_min_vector = mu_min_in;
  mu_max_vector = mu_max_in;
}

Real RBBase::get_parameter_min(unsigned int i) const
{
  libmesh_assert(i < mu_min_vector.size());

  return mu_min_vector[i];
}

Real RBBase::get_parameter_max(unsigned int i) const
{
  libmesh_assert(i < mu_max_vector.size());

  return mu_max_vector[i];
}

void RBBase::set_current_parameters(const std::vector<Real>& params)
{
  libmesh_assert(params.size() == get_n_params());

  for(unsigned int i=0; i<params.size(); i++)
    libmesh_assert( (mu_min_vector[i] <= libmesh_real(params[i])) &&
                    (libmesh_real(params[i]) <= mu_max_vector[i]) );

  current_parameters = params;
}

void RBBase::print_current_parameters()
{
  for(unsigned int j=0; j<get_n_params(); j++)
  {
    libMesh::out << "mu[" << j << "] = " << current_parameters[j] << std::endl;
  }
  libMesh::out << std::endl;
}

bool RBBase::valid_params(const std::vector<Real>& params)
{
  bool valid = ( params.size() == get_n_params() );

  if(!valid)
  {
   return false;
  }
  else
  {
    for(unsigned int i=0; i<params.size(); i++)
    {
      valid = valid && ( (mu_min_vector[i] <= params[i]) &&
                         (params[i] <= mu_max_vector[i]) );
    }
  }

  return valid;
}

void RBBase::broadcast_current_parameters(unsigned int proc_id)
{
  libmesh_assert(proc_id < libMesh::n_processors());

  Parallel::broadcast(current_parameters, proc_id);
}

AutoPtr<RBThetaExpansion> RBBase::build_rb_theta_expansion(std::vector<Real>& parameters_ref)
{
  return AutoPtr<RBThetaExpansion>(new RBThetaExpansion( parameters_ref ));
}

void RBBase::init_extra_data_objects ()
{
  rb_theta_expansion = build_rb_theta_expansion(current_parameters);
}

} // namespace libMesh
