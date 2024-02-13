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

// Local includes
#include "libmesh/rb_theta.h"
#include "libmesh/rb_parameters.h"
#include "libmesh/int_range.h"

namespace libMesh
{

Number RBTheta::evaluate(const RBParameters & mu)
{
  // The RBTheta::evaluate() API is not general enough to handle the
  // multi-sample RBParameters case, and you must therefore call
  // RBTheta::evaluate_vec() instead.
  libmesh_error_msg_if(mu.n_samples() > 1,
                       "You should only call the evaluate_vec() API when using multi-sample RBParameters objects.");

  return 1.;
}

std::vector<Number>
RBTheta::evaluate_vec(const std::vector<RBParameters> & mus)
{
  // Eventual return value
  std::vector<Number> result;

  for (const auto & mu : mus)
    {
      // Backwards-compatible behavior: for single-sample RBParameters objects, we fall back on
      // calling the scalar evaluate() function for this RBTheta object, which may have been
      // overridden by the user.
      if (mu.n_samples() == 1)
        result.push_back( this->evaluate(mu) );
      else
        {
          // For multi-sample RBParameters objects, all we can do is return
          // mu.n_samples() copies of 1 here at the base class level.
          result.insert(result.end(), /*count=*/mu.n_samples(), /*val=*/1.);
        }
    }

  return result;
}

}
