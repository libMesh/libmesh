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

#include "rb_theta_expansion.h"
#include "rb_theta.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBThetaExpansion implementation

RBThetaExpansion::RBThetaExpansion()
{
  theta_q_a_vector.clear();
  theta_q_f_vector.clear();
  theta_q_l_vector.clear();
}

unsigned int RBThetaExpansion::get_Q_a()
{
  return theta_q_a_vector.size();
}

unsigned int RBThetaExpansion::get_Q_f() const
{
  return theta_q_f_vector.size();
}

unsigned int RBThetaExpansion::get_n_outputs() const
{
  return theta_q_l_vector.size();
}

unsigned int RBThetaExpansion::get_Q_l(unsigned int index) const
{
  if(index >= get_n_outputs())
  {
    libMesh::err << "Error: We must have index < n_outputs in get_Q_l."
                 << std::endl;
    libmesh_error();
  }
  return theta_q_l_vector[index].size();
}

void RBThetaExpansion::attach_theta_q_a(RBTheta* theta_q_a)
{
  libmesh_assert(theta_q_a != NULL);

  theta_q_a_vector.push_back(theta_q_a);
}

void RBThetaExpansion::attach_theta_q_f(RBTheta* theta_q_f)
{
  libmesh_assert(theta_q_f != NULL);

  theta_q_f_vector.push_back(theta_q_f);
}

void RBThetaExpansion::attach_output_theta(std::vector<RBTheta*> theta_q_l)
{
  theta_q_l_vector.push_back(theta_q_l);
}

void RBThetaExpansion::attach_output_theta(RBTheta* theta_q_l)
{
  libmesh_assert(theta_q_l != NULL);

  std::vector<RBTheta*> theta_l_vector(1);
  theta_l_vector[0] = theta_q_l;
  
  attach_output_theta(theta_l_vector);
}

Number RBThetaExpansion::eval_theta_q_a(unsigned int q, const std::vector<Real>& mu)
{
  if(q >= get_Q_a())
  {
    libMesh::err << "Error: We must have q < Q_a in eval_theta_q_a."
                 << std::endl;
    libmesh_error();
  }

  return theta_q_a_vector[q]->evaluate( mu );

}

Number RBThetaExpansion::eval_theta_q_f(unsigned int q, const std::vector<Real>& mu)
{
  if(q >= get_Q_f())
  {
    libMesh::err << "Error: We must have q < Q_f in eval_theta_q_f."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_f_vector[q] != NULL);

  return theta_q_f_vector[q]->evaluate( mu );
}

Number RBThetaExpansion::eval_theta_q_l(unsigned int output_index, unsigned int q_l, const std::vector<Real>& mu)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_Q_l(output_index)) )
  {
    libMesh::err << "Error: We must have output_index < n_outputs and "
                 << "q_l < get_Q_l(output_index) in eval_theta_q_l."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_l_vector[output_index][q_l] != NULL);

  return theta_q_l_vector[output_index][q_l]->evaluate( mu );
}


}