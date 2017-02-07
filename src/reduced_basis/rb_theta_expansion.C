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
#include "libmesh/rb_theta_expansion.h"
#include "libmesh/rb_theta.h"
#include "libmesh/rb_parameters.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBThetaExpansion implementation

RBThetaExpansion::RBThetaExpansion()
{
}

unsigned int RBThetaExpansion::get_n_A_terms() const
{
  return cast_int<unsigned int>
    (_A_theta_vector.size());
}

unsigned int RBThetaExpansion::get_n_F_terms() const
{
  return cast_int<unsigned int>
    (_F_theta_vector.size());
}

unsigned int RBThetaExpansion::get_n_outputs() const
{
  return cast_int<unsigned int>
    (_output_theta_vector.size());
}

unsigned int RBThetaExpansion::get_n_output_terms(unsigned int index) const
{
  if(index >= get_n_outputs())
    libmesh_error_msg("Error: We must have index < n_outputs in get_Q_l.");

  return cast_int<unsigned int>
    (_output_theta_vector[index].size());
}

void RBThetaExpansion::attach_A_theta(RBTheta * theta_q_a)
{
  libmesh_assert(theta_q_a);

  _A_theta_vector.push_back(theta_q_a);
}

void RBThetaExpansion::attach_multiple_A_theta(std::vector<RBTheta *> theta_q_a)
{
  for (std::size_t i=0; i<theta_q_a.size(); i++)
    {
      libmesh_assert(theta_q_a[i]);
      _A_theta_vector.push_back(theta_q_a[i]);
    }
}

void RBThetaExpansion::attach_F_theta(RBTheta * theta_q_f)
{
  libmesh_assert(theta_q_f);

  _F_theta_vector.push_back(theta_q_f);
}

void RBThetaExpansion::attach_multiple_F_theta(std::vector<RBTheta *> theta_q_f)
{
  for (std::size_t i=0; i<theta_q_f.size(); i++)
    {
      libmesh_assert(theta_q_f[i]);
      _F_theta_vector.push_back(theta_q_f[i]);
    }
}

void RBThetaExpansion::attach_output_theta(std::vector<RBTheta *> theta_q_l)
{
  _output_theta_vector.push_back(theta_q_l);
}

void RBThetaExpansion::attach_output_theta(RBTheta * theta_q_l)
{
  libmesh_assert(theta_q_l);

  std::vector<RBTheta *> theta_l_vector(1);
  theta_l_vector[0] = theta_q_l;

  attach_output_theta(theta_l_vector);
}

Number RBThetaExpansion::eval_A_theta(unsigned int q,
                                      const RBParameters & mu)
{
  if(q >= get_n_A_terms())
    libmesh_error_msg("Error: We must have q < get_n_A_terms in eval_A_theta.");

  libmesh_assert(_A_theta_vector[q]);

  return _A_theta_vector[q]->evaluate( mu );
}

Number RBThetaExpansion::eval_F_theta(unsigned int q,
                                      const RBParameters & mu)
{
  if(q >= get_n_F_terms())
    libmesh_error_msg("Error: We must have q < get_n_F_terms in eval_F_theta.");

  libmesh_assert(_F_theta_vector[q]);

  return _F_theta_vector[q]->evaluate( mu );
}

Number RBThetaExpansion::eval_output_theta(unsigned int output_index,
                                           unsigned int q_l,
                                           const RBParameters & mu)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_n_output_terms(output_index)) )
    libmesh_error_msg("Error: We must have output_index < n_outputs and " \
                      << "q_l < get_n_output_terms(output_index) in eval_output_theta.");

  libmesh_assert(_output_theta_vector[output_index][q_l]);

  return _output_theta_vector[output_index][q_l]->evaluate( mu );
}


}
