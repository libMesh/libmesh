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

#include "rb_parametrized_object.h"
#include "parallel.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBParametrizedObject implementation

RBParametrizedObject::RBParametrizedObject()
{
  libmesh_experimental();
}

RBParametrizedObject::~RBParametrizedObject ()
{
}

unsigned int RBParametrizedObject::get_n_params() const
{
  return current_parameters.size();
}

void RBParametrizedObject::set_n_params(unsigned int n_params_in)
{
  current_parameters.clear();
  current_parameters.resize(n_params_in);
}

void RBParametrizedObject::set_current_parameters(const std::vector<Real>& params)
{
  current_parameters = params;
}

std::vector<Real>& RBParametrizedObject::get_current_parameters()
{
  return current_parameters;
}

void RBParametrizedObject::print_current_parameters()
{
  for(unsigned int j=0; j<get_n_params(); j++)
  {
    libMesh::out << "mu[" << j << "] = " << current_parameters[j] << std::endl;
  }
  libMesh::out << std::endl;
}

void RBParametrizedObject::broadcast_current_parameters(unsigned int proc_id)
{
  libmesh_assert(proc_id < libMesh::n_processors());

  Parallel::broadcast(current_parameters, proc_id);
}

} // namespace libMesh
