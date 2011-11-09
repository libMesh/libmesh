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

#ifndef __rb_eim_theta_h__
#define __rb_eim_theta_h__

#include "rb_theta.h"
#include "rb_eim_evaluation.h"

namespace libMesh
{

/**
 * This class provides functionality required to define an RBTheta
 * object that arises from an "Empirical Interpolation Method" (EIM)
 * approximation.
 *
 * @author David J. Knezevic, 2011
 */

class RBEIMTheta : public RBTheta
{
public:

  /**
   * Constructor.
   */
  RBEIMTheta(RBEIMEvaluation& rb_eim_eval_in, unsigned int index_in);

  /**
   * Evaluate this RBEIMTheta object at the parameter \p mu.
   * This entails solving the RB EIM approximation and picking
   * out the appropriate coefficient.
   */
  virtual Number evaluate(const std::vector<Real>& mu);
  
  /**
   * The RBEIMEvaluation object that this RBEIMTheta is based on.
   */
  RBEIMEvaluation& rb_eim_eval;

  /**
   * The index of the RB_solution vector that we pick out
   * from rb_eim_eval to provide the value of the evaluation.
   */
  unsigned int index;
};

RBEIMTheta::RBEIMTheta(RBEIMEvaluation& rb_eim_eval_in, unsigned int index_in)
  :
  rb_eim_eval(rb_eim_eval_in),
  index(index_in)
{}

Number RBEIMTheta::evaluate(const std::vector<Real>& mu)
{
  rb_eim_eval.set_current_parameters(mu);
  rb_eim_eval.rb_solve(rb_eim_eval.get_n_basis_functions());
    
  return rb_eim_eval.RB_solution(index);
}



}
 
#endif