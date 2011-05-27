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

#ifndef __rb_base_h__
#define __rb_base_h__

#include "auto_ptr.h"
#include "rb_theta.h"
#include "rb_theta_expansion.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * This is the base class for certified reduced basis (RB)
 * functionality. Here we store generic RB quantities such
 * as parameter ranges and functions in the `affine'
 * expansion of the bilinear form operators of the PDE.
 *
 * @author David J. Knezevic, 2011
 */


// ------------------------------------------------------------
// RBBase class definition
class RBBase
{
public:

  /**
   * Constructor.
   */
  RBBase ();

  /**
   * Destructor.
   */
  virtual ~RBBase ();

  /**
   * Get the number of parameters. Value is determined
   * by specifying the parameter ranges.
   */
  unsigned int get_n_params() const;

  /**
   * Set the parameter range for this RB object.
   */
  void set_parameter_range(std::vector<Real> mu_min, std::vector<Real> mu_max);

  /**
   * Get range of i^th parameter.
   */
  Real get_parameter_min(unsigned int i) const;
  Real get_parameter_max(unsigned int i) const;

  /**
   * Get/set current parameters.
   */
  std::vector<Real>& get_current_parameters() { return current_parameters; };
  void set_current_parameters(const std::vector<Real>& params);

  /**
   * Print the current parameters.
   */
  void print_current_parameters();

  /**
   * Broadcasts current_parameters on processor proc_id
   * to all processors.
   */
  void broadcast_current_parameters(unsigned int proc_id);
  
  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * A pointer to to the object that stores the theta expansion.
   * This is not an AutoPtr since we may want to share it between
   * multiple RBBases. (Note: a shared_ptr would be a good option here.)
   */
  RBThetaExpansion* rb_theta_expansion;
  
protected:

  /**
   * Helper function to indicate if the input parameters are valid.
   */
  bool valid_params(const std::vector<Real>& params);

  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * Vector of parameter ranges.
   */
  std::vector<Real> mu_min_vector;
  std::vector<Real> mu_max_vector;

  /**
   * Vector storing the current parameters.
   */
  std::vector<Real> current_parameters;

};

} // namespace libMesh


#endif
