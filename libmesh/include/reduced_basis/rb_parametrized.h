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

#ifndef __rb_parametrized_h__
#define __rb_parametrized_h__

// libMesh includes
#include "reference_counted_object.h"

// C++ includes
#include <vector>

namespace libMesh
{
    
/**
 * This class is part of the rbOOmit framework.
 *
 * This class defines basic functionality of
 * a parametrized object.
 *
 * @author David J. Knezevic, 2011
 */
// ------------------------------------------------------------
// RBParametrized class definition
class RBParametrized : public ReferenceCountedObject<RBParametrized>
{
public:

  /**
   * Constructor.
   */
  RBParametrized ();

  /**
   * Destructor.
   */
  virtual ~RBParametrized ();

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Initialize the parameter ranges and set current_parameters.
   */
  void initialize_parameters(std::vector<Real> mu_min_in,
                             std::vector<Real> mu_max_in,
                             std::vector<Real> mu_in);

  /**
   * Initialize the parameter ranges and set current_parameters.
   */
  void initialize_parameters(RBParametrized& rb_parametrized);
  
  /**
   * Initialize the parameter ranges and set current_parameters by
   * reading in data from the file \p input_filename
   */
  virtual void initialize_parameters(const std::string& parameters_filename);

  /**
   * Get the number of parameters.
   */
  unsigned int get_n_params() const;

  /**
   * Get the current parameters.
   */
  std::vector<Real> get_parameters() const;

  /**
   * Set the current parameters to \p params
   */
  void set_parameters(const std::vector<Real>& params);

  /**
   * Get a vector that specifies the minimum allowable value
   * for each parameter.
   */
  std::vector<Real> get_parameters_min_vector() const;

  /**
   * Get a vector that specifies the maximum allowable value
   * for each parameter.
   */
  std::vector<Real> get_parameters_max_vector() const;

  /**
   * Get minimum allowable value of parameter \p i.
   */
  Real get_parameter_min(unsigned int i) const;

  /**
   * Get maximum allowable value of parameter \p i.
   */
  Real get_parameter_max(unsigned int i) const;

  /**
   * Print the current parameters.
   */
  void print_parameters();

private:

  /**
   * Helper function to check that \p params is valid.
   */
  bool valid_params(const std::vector<Real>& params);

  //--------------- PRIVATE DATA MEMBERS ---------------//

  /**
   * Flag indicating whether the parameters have been initialized.
   */
  bool parameters_initialized;

  /**
   * Vector storing the current parameters.
   */
  std::vector<Real> parameters;

  /**
   * Vectors that define the ranges (min and max) for the parameters.
   */
  std::vector<Real> parameters_min_vector;
  std::vector<Real> parameters_max_vector;

};

} // namespace libMesh


#endif
