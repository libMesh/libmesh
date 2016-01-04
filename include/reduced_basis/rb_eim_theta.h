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

#ifndef LIBMESH_RB_EIM_THETA_H
#define LIBMESH_RB_EIM_THETA_H

// rbOOmit includes
#include "libmesh/rb_theta.h"

// C++ includes

namespace libMesh
{

class RBParameters;
class RBEIMEvaluation;

/**
 * This class provides functionality required to define an RBTheta
 * object that arises from an "Empirical Interpolation Method" (EIM)
 * approximation.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class RBEIMTheta : public RBTheta
{
public:

  /**
   * Constructor.
   */
  RBEIMTheta(RBEIMEvaluation & rb_eim_eval_in, unsigned int index_in);

  /**
   * Evaluate this RBEIMTheta object at the parameter \p mu.
   * This entails solving the RB EIM approximation and picking
   * out the appropriate coefficient.
   */
  virtual Number evaluate(const RBParameters & mu) libmesh_override;

  /**
   * The RBEIMEvaluation object that this RBEIMTheta is based on.
   */
  RBEIMEvaluation & rb_eim_eval;

  /**
   * The index of the RB_solution vector that we pick out
   * from rb_eim_eval to provide the value of the evaluation.
   */
  unsigned int index;
};

}

#endif // LIBMESH_RB_EIM_THETA_H
