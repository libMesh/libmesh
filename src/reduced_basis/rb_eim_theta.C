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

#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parameters.h"
#include "libmesh/rb_eim_evaluation.h"

namespace libMesh
{

RBEIMTheta::RBEIMTheta(RBEIMEvaluation & rb_eim_eval_in, unsigned int index_in)
  :
  rb_eim_eval(rb_eim_eval_in),
  index(index_in)
{
}

Number RBEIMTheta::evaluate(const RBParameters & mu)
{
  if (mu.n_parameters() > rb_eim_eval.get_n_params())
    {
      // In this case the parameters related to the EIM are a subset of
      // the parameters from the associated RB problem, hence we need to "pull out"
      // the parameters related to the EIM
      RBParameters mu_eim;
      for (const auto & pr : rb_eim_eval.get_parameters())
        {
          const std::string & param_name = pr.first;
          mu_eim.set_value(param_name, mu.get_value(param_name));
        }
      rb_eim_eval.set_parameters(mu_eim);
    }
  else
    {
      rb_eim_eval.set_parameters(mu);
    }

  rb_eim_eval.rb_eim_solve(rb_eim_eval.get_n_basis_functions());

  const DenseVector<Number> & eim_solution = rb_eim_eval.get_rb_eim_solution();
  if(index >= eim_solution.size())
  {
    libmesh_error_msg("Error: Invalid EIM solution index, " + std::to_string(index));
  }

  return eim_solution(index);
}

}
