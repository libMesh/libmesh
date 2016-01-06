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

#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/rb_theta.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBThetaExpansion implementation

TransientRBThetaExpansion::TransientRBThetaExpansion()
{
}

Number TransientRBThetaExpansion::eval_M_theta(unsigned int q,
                                               const RBParameters & mu)
{
  if(q >= get_n_M_terms())
    libmesh_error_msg("Error: We must have q < get_n_M_terms in eval_M_theta.");

  libmesh_assert(_M_theta_vector[q]);

  return _M_theta_vector[q]->evaluate( mu );
}

void TransientRBThetaExpansion::attach_M_theta(RBTheta * theta_q_m)
{
  libmesh_assert(theta_q_m);

  _M_theta_vector.push_back(theta_q_m);
}


}
