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

#ifndef LIBMESH_RB_CONSTRUCTION_BASE_H
#define LIBMESH_RB_CONSTRUCTION_BASE_H

// rbOOmit includes
#include "libmesh/id_types.h"
#include "libmesh/rb_parametrized.h"
#include "libmesh/rb_theta_expansion.h"
#include "libmesh/rb_theta.h"

// libMesh includes
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_solver.h"
#include "libmesh/perf_log.h"

// C++ includes
#include <set>

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * This is the base class for the Construction stage
 * of the certified reduced basis (RB) method.
 * We template the Base class so that we can derive from
 * the appropriate libMesh System type (e.g. LinearImplicitSystem
 * for standard reduced basis, EigenSystem for SCM)
 * at compile time.
 *
 * \author David J. Knezevic
 * \date 2009
 */
template<class Base>
class RBConstructionBase : public Base, public RBParametrized
{
public:

  /**
   * Constructor.
   */
  RBConstructionBase (EquationSystems & es,
                      const std::string & name,
                      const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions as the union of its
   *   potential base classes (currently LinearImplicitSystem and
   *   EigenSystem).
   * - Destructor is defaulted out-of-line.
   */
  RBConstructionBase (RBConstructionBase &&) = default;
  RBConstructionBase & operator= (RBConstructionBase &&) = delete;
  RBConstructionBase (const RBConstructionBase &) = delete;
  RBConstructionBase & operator= (const RBConstructionBase &) = delete;
  virtual ~RBConstructionBase ();

  /**
   * The type of system.
   */
  typedef RBConstructionBase<Base> sys_type;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Set the quiet_mode flag. If quiet == false then
   * we print out a lot of extra information
   * during the Offline stage.
   */
  void set_quiet_mode(bool quiet_mode_in)
  { this->quiet_mode = quiet_mode_in; }

  /**
   * Is the system in quiet mode?
   */
  bool is_quiet() const
  { return this->quiet_mode; }

  /**
   * Set the boolean option that indicates if we normalization
   * solution snapshots or not.
   */
  void set_normalize_solution_snapshots(bool value);

  /**
   * Get the number of global training samples.
   */
  numeric_index_type get_n_training_samples() const;

  /**
   * Get the total number of training samples local to this processor.
   */
  numeric_index_type get_local_n_training_samples() const;

  /**
   * Get the first local index of the training parameters.
   */
  numeric_index_type get_first_local_training_index() const;

  /**
   * Get the last local index of the training parameters.
   */
  numeric_index_type get_last_local_training_index() const;

  /**
   * Initialize the parameter ranges and indicate whether deterministic
   * or random training parameters should be used and whether or
   * not we want the parameters to be scaled logarithmically.
   * \p n_global_training_samples is the total number of samples to
   * generate, which will be distributed across all the processors.
   */
  virtual void initialize_training_parameters(const RBParameters & mu_min,
                                              const RBParameters & mu_max,
                                              const unsigned int n_global_training_samples,
                                              const std::map<std::string, bool> & log_param_scale,
                                              const bool deterministic=true);

  /**
   * Overwrite the training parameters with new_training_set.
   * This training set is assumed to contain only the samples local to this processor.
   */
  virtual void load_training_set(const std::map<std::string, std::vector<RBParameter>> & new_training_set);

  /**
   * Overwrite the local training samples for \p param_name using \p values.
   * This assumes that values.size() matches get_local_n_training_samples().
   */
  void set_training_parameter_values(const std::string & param_name, const std::vector<RBParameter> & values);

  /**
   * Broadcasts parameters from processor proc_id to all processors.
   * This broadcasts the RBParameters object from .get_parameters(),
   * and then sets it on all processors with .set_parameters().
   */
  void broadcast_parameters(const unsigned int proc_id);

  /**
   * Set the seed that is used to randomly generate training parameters.
   */
  void set_training_random_seed(int seed);

  /**
   * In some cases we only want to allow discrete parameter values, instead
   * of parameters that may take any value in a specified interval.
   * Here we provide a method to set the d
   * Set the discrete values for parameter \p mu that are allowed in the
   * training set. This must be called before the training set is generated.
   */

  /**
   * Set the name of the parameter that we will generate deterministic training parameters for.
   * Defaults to "NONE".
   */
  void set_deterministic_training_parameter_name(const std::string & name);

  /**
   * Get the name of the parameter that we will generate deterministic training parameters for.
   */
  const std::string & get_deterministic_training_parameter_name() const;

  /**
   * Static helper function for generating a randomized set of parameters.
   * The parameter \p n_global_training_samples_in is the total number of parameters to
   * generate, and they will be split across all the processors (unless serial_training_set=true)
   * in the \p local_training_parameters_in map.
   * \return a pair of {first_local_index,last_local_index}
   */
  static std::pair<std::size_t, std::size_t>
  generate_training_parameters_random(const Parallel::Communicator & communicator,
                                      const std::map<std::string, bool> & log_param_scale,
                                      std::map<std::string, std::vector<RBParameter>> & local_training_parameters_in,
                                      const unsigned int n_global_training_samples_in,
                                      const RBParameters & min_parameters,
                                      const RBParameters & max_parameters,
                                      const int training_parameters_random_seed=-1,
                                      const bool serial_training_set=false);

  /**
   * Static helper function for generating a deterministic set of parameters. Only works with 1 or 2
   * parameters (as defined by the lengths of min/max parameters vectors), otherwise throws an error.
   * The parameter \p n_global_training_samples_in is the total number of parameters to
   * generate, and they will be split across all the processors (unless serial_training_set=true)
   * in the \p local_training_parameters_in map.
   * \return a pair of {first_local_index,last_local_index}
   */
  static std::pair<std::size_t, std::size_t>
  generate_training_parameters_deterministic(const Parallel::Communicator & communicator,
                                             const std::map<std::string, bool> & log_param_scale,
                                             std::map<std::string, std::vector<RBParameter>> & local_training_parameters_in,
                                             const unsigned int n_global_training_samples_in,
                                             const RBParameters & min_parameters,
                                             const RBParameters & max_parameters,
                                             const bool serial_training_set=false);

protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Return the RBParameters in index \p global_index of the global training set.
   * Why do we use an index here? RBParameters supports loading the full sample set.
   * This seems probably unnecessary now to load individually. Maybe it's a memory issue?
   */
  RBParameters get_params_from_training_set(unsigned int global_index);

  /**
   * Set parameters to the RBParameters stored in index \p global_index of the global training set.
   */
  void set_params_from_training_set(unsigned int global_index);

  /**
   * Load the specified training parameter and then broadcast to all processors.
   */
  virtual void set_params_from_training_set_and_broadcast(unsigned int global_index);

  /**
   * Static function to return the error pair (index,error)
   * that is corresponds to the largest error on all
   * processors.
   */
  static void get_global_max_error_pair(const Parallel::Communicator & communicator,
                                        std::pair<numeric_index_type, Real> & error_pair);


  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * Flag to indicate whether we print out extra information during
   * the Offline stage.
   */
  bool quiet_mode;

  /**
   * This boolean flag indicates whether or not the training set should
   * be the same on all processors. By default it is false, but in the
   * case of the Empirical Interpolation Method (RBEIMConstruction),
   * for example, we need the training set to be identical on all processors.
   */
  bool serial_training_set;

  /**
   * Set this boolean to true if we want to normalize solution snapshots
   * used in training to have norm of 1. This is relevant if snapshots
   * have differing magnitudes and we want to approximate them all with
   * equal accuracy.
   */
  bool _normalize_solution_snapshots;

  /**
   * We keep an extra temporary vector that is useful for
   * performing inner products (avoids unnecessary memory
   * allocation/deallocation).
   */
  std::unique_ptr<NumericVector<Number>> inner_product_storage_vector;


private:

  /**
   * Boolean flag to indicate whether or not the
   * parameter ranges have been initialized.
   */
  bool _training_parameters_initialized;

  /**
   * The training samples for each parameter.
   * When serial_training_set is true, the map contains all samples of all parameters.
   * Otherwise, the sample vectors will only contain the values for the local samples
   * as defined by _first_local_index and _n_local_training_samples.
   * Mapped from parameter_name -> sample_vector -> value_vector.
   */
  std::map<std::string, std::vector<RBParameter>> _training_parameters;

  /**
   * The first sample-vector index from the global vector which is stored in
   * the _training_parameters on this processor.
   * _n_local_training_samples is equivalent to the .size() of any vector in _training_parameters.
   */
  numeric_index_type _first_local_index;
  numeric_index_type _n_local_training_samples;
  numeric_index_type _n_global_training_samples;

  /**
   * If < 0, use std::time() * processor_id() to seed the random
   * number generator for the training parameters (default).  If
   * >= 0, use the provided value * processor_id() as the random
   * number generator seed.
   */
  int _training_parameters_random_seed;
};

} // namespace libMesh


#endif // LIBMESH_RB_CONSTRUCTION_BASE_H
