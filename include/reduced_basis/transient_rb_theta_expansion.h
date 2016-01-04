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

#ifndef LIBMESH_TRANSIENT_RB_THETA_EXPANSION_H
#define LIBMESH_TRANSIENT_RB_THETA_EXPANSION_H

// rbOOmit includes
#include "libmesh/rb_theta_expansion.h"

// libMesh includes

// C++ includes


namespace libMesh
{

/**
 * This class stores the set of RBTheta functor objects that define
 * the "parameter-dependent expansion" of a PDE.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class TransientRBThetaExpansion : public RBThetaExpansion
{
public:

  /**
   * Constructor.
   */
  TransientRBThetaExpansion();

  /**
   * The type of the parent.
   */
  typedef RBThetaExpansion Parent;

  /**
   * Evaluate theta at the current parameter. Overload
   * if the theta functions need to be treated differently
   * in subclasses.
   */
  virtual Number eval_M_theta(unsigned int q,
                              const RBParameters & mu);

  /**
   * Get Q_m, the number of terms in the affine
   * expansion for the mass operator.
   */
  virtual unsigned int get_n_M_terms()
  { return cast_int<unsigned int>(_M_theta_vector.size()); }

  /**
   * Attach a pointer to a functor object that defines one
   * of the theta_q_m terms.
   */
  virtual void attach_M_theta(RBTheta * theta_q_m);

private:

  /**
   * Vector storing the pointers to the RBTheta functors.
   */
  std::vector<RBTheta *> _M_theta_vector;
};

}

#endif // LIBMESH_TRANSIENT_RB_THETA_EXPANSION_H
