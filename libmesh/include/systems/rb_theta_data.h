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

#ifndef __rb_theta_data_h__
#define __rb_theta_data_h__

// C++ includes
#include <vector>

// Local includes
#include "reference_counted_object.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBThetaData provides a wrapper for any data required
 * in the evaluation of parameter dependent functions
 * (the theta_q functions). In the simplest case this
 * just provides access to a system's current_parameters.
 *
 * @author David J. Knezevic, 2009
 */
class RBThetaData : public ReferenceCountedObject<RBThetaData>
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  RBThetaData () {};
  
  /**
   * Destructor.
   */
  virtual ~RBThetaData () {};
  
  /**
   * Store a const pointer to an associated system's
   * parameter vector.
   */
  void point_to_parameters(const std::vector<Real>& mu_in);
  
  /**
   * Return a const reference to mu. Overload in order to
   * provide different behavior.
   */
  virtual const std::vector<Real>& get_mu();

  /**
   * A pointer to the parameters in an associated system.
   */
  const std::vector<Real>* mu;
};

inline void RBThetaData::point_to_parameters(const std::vector<Real>& mu_in)
{
  mu = &mu_in;
}

inline const std::vector<Real>& RBThetaData::get_mu()
{
  return *mu;
}
 
}

#endif
