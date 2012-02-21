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

#include "transient_rb_theta_expansion.h"
#include "rb_theta.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBThetaExpansion implementation

TransientRBThetaExpansion::TransientRBThetaExpansion()
{
  theta_q_m_vector.clear();
}

Number TransientRBThetaExpansion::eval_theta_q_m(unsigned int q,
                                                 const std::vector<Real>& mu)
{
  if(q >= get_Q_m())
  {
    libMesh::err << "Error: We must have q < Q_m in eval_theta_q_m."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_m_vector[q] != NULL);

  return theta_q_m_vector[q]->evaluate( mu );
}

void TransientRBThetaExpansion::attach_theta_q_m(RBTheta* theta_q_m)
{
  libmesh_assert(theta_q_m != NULL);

  theta_q_m_vector.push_back(theta_q_m);
}


}