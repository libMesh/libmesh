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

#ifndef LIBMESH_RB_PARAMETRIZED_H
#define LIBMESH_RB_PARAMETRIZED_H

// rbOOmit includes
#include "libmesh/rb_parameters.h"

// libMesh includes
#include "libmesh/reference_counted_object.h"

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
  void initialize_parameters(const RBParameters& mu_min_in,
                             const RBParameters& mu_max_in,
                             const RBParameters& mu_in);

  /**
   * Initialize the parameter ranges and set current_parameters.
   */
  void initialize_parameters(const RBParametrized& rb_parametrized);

  /**
   * Get the number of parameters.
   */
  unsigned int get_n_params() const;

  /**
   * Get the current parameters.
   */
  const RBParameters& get_parameters() const;

  /**
   * Set the current parameters to \p params
   */
  void set_parameters(const RBParameters& params);

  /**
   * Get an RBParameters object that specifies the minimum allowable value
   * for each parameter.
   */
  const RBParameters& get_parameters_min() const;

  /**
   * Get an RBParameters object that specifies the maximum allowable value
   * for each parameter.
   */
  const RBParameters& get_parameters_max() const;

  /**
   * Get minimum allowable value of parameter \p param_name.
   */
  Real get_parameter_min(const std::string& param_name) const;

  /**
   * Get maximum allowable value of parameter \p param_name.
   */
  Real get_parameter_max(const std::string& param_name) const;

  /**
   * Print the current parameters.
   */
  void print_parameters() const;

  /**
   * Write out the parameter ranges to file.
   */
  void write_parameter_ranges_to_file(const std::string& file_name,
                                      const bool write_binary);

  /**
   * Read in the parameter ranges from file. Initialize parameters
   * to the "minimum" parameter values.
   */
  void read_parameter_ranges_from_file(const std::string& file_name,
                                       const bool read_binary);

  /**
   * Public boolean to toggle verbose mode.
   */
  bool verbose_mode;

private:

  /**
   * Helper function to check that \p params is valid.
   */
  bool valid_params(const RBParameters& params);

  //--------------- PRIVATE DATA MEMBERS ---------------//

  /**
   * Flag indicating whether the parameters have been initialized.
   */
  bool parameters_initialized;

  /**
   * Vector storing the current parameters.
   */
  RBParameters parameters;

  /**
   * Vectors that define the ranges (min and max) for the parameters.
   */
  RBParameters parameters_min;
  RBParameters parameters_max;

};

} // namespace libMesh


#endif // LIBMESH_RB_PARAMETRIZED_H
