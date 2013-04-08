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
 * @author David J. Knezevic, 2009
 */


// ------------------------------------------------------------
// RBConstructionBase class definition
template<class Base>
class RBConstructionBase : public Base, public RBParametrized
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  RBConstructionBase (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBConstructionBase ();

  /**
   * The type of system.
   */
  typedef RBConstructionBase<Base> sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Get the total number of training samples.
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
   */
  virtual void initialize_training_parameters(const RBParameters& mu_min,
                                              const RBParameters& mu_max,
                                              unsigned int n_training_parameters,
                                              std::map<std::string, bool> log_param_scale,
                                              bool deterministic=true);

  /**
   * Overwrite the training parameters with new_training_set.
   */
  virtual void load_training_set(std::map< std::string, std::vector<Number> >& new_training_set);

  /**
   * Changes the current PC (and iterative solver, if desired) in the
   * passed-in LinearSolver object to an alternative solver specified
   * by the alternative_solver string stored in this class.  You might
   * use this to e.g. switch to a sparse direct solver for the multiple
   * RHS solves executed during the update_residual_terms function.
   * The return strings are names of the original PC and KSP objects,
   * you can reset these using the reset_alternative_solver() function below.
   */
  std::pair<std::string,std::string> set_alternative_solver(AutoPtr<LinearSolver<Number> >& ls);

  /**
   * Resets the PC (and iterative solver, if desired) in the passed-in
   * LinearSolver object to the values specified in the pair of strings
   * passed as the second argument.  If the "alternative_solver" string,
   * defined below, is "unchanged", this function does nothing.
   */
  void reset_alternative_solver(AutoPtr<LinearSolver<Number> >& ls,
				const std::pair<std::string,std::string>& orig);

  /**
   * Broadcasts parameters on processor proc_id
   * to all processors.
   */
  void broadcast_parameters(unsigned int proc_id);

  /**
   * Set the seed that is used to randomly generate training parameters.
   */
  void set_training_random_seed(unsigned int seed);

  /**
   * Set the name of the parameter that we will generate deterministic training parameters for.
   * Defaults to "NONE".
   */
  void set_deterministic_training_parameter_name(const std::string name);

  /**
   * Get the name of the parameter that we will generate deterministic training parameters for.
   */
  const std::string& get_deterministic_training_parameter_name() const;

  /**
   * Set the number of times each sample of the deterministic training parameter is repeated.
   */
  void set_deterministic_training_parameter_repeats(unsigned int repeats);

  /**
   * Get the number of times each sample of the deterministic training parameter is repeated.
   */
  unsigned int get_deterministic_training_parameter_repeats() const;

protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Return the RBParameters in index \p index of training set.
   */
  RBParameters get_params_from_training_set(unsigned int index);

  /**
   * Set parameters to the RBParameters stored in index \p index of the training set.
   */
  void set_params_from_training_set(unsigned int index);

  /**
   * Load the specified training parameter and then broadcast to all processors.
   */
  virtual void set_params_from_training_set_and_broadcast(unsigned int index);

  /**
   * Static function to return the error pair (index,error)
   * that is corresponds to the largest error on all
   * processors.
   */
  static void get_global_max_error_pair(const Parallel::Communicator &communicator,
					std::pair<unsigned int, Real>& error_pair);

  /**
   * Static helper function for generating a randomized set of parameters.
   */
  static void generate_training_parameters_random(const Parallel::Communicator &communicator,
						  std::map<std::string, bool> log_param_scale,
                                                  std::map< std::string, NumericVector<Number>* >& training_parameters_in,
                                                  unsigned int n_training_samples_in,
                                                  const RBParameters& min_parameters,
                                                  const RBParameters& max_parameters,
                                                  int training_parameters_random_seed=-1,
                                                  bool serial_training_set=false);

  /**
   * Static helper function for generating a "partially" random set of parameters, that is
   * the parameter indicated by this->get_deterministic_training_parameter() will be
   * deterministic.
   */
  static void generate_training_parameters_partially_random(const Parallel::Communicator &communicator,
							    const std::string& deterministic_parameter_name,
                                                            const unsigned int deterministic_parameter_repeats,
                                                            std::map<std::string, bool> log_param_scale,
                                                            std::map< std::string, NumericVector<Number>* >& training_parameters_in,
                                                            unsigned int n_deterministic_training_samples_in,
                                                            const RBParameters& min_parameters,
                                                            const RBParameters& max_parameters,
                                                            int training_parameters_random_seed=-1,
                                                            bool serial_training_set=false);

  /**
   * Static helper function for generating a deterministic set of parameters. Only works with 1 or 2
   * parameters (as defined by the lengths of min/max parameters vectors), otherwise throws an error.
   */
  static void generate_training_parameters_deterministic(const Parallel::Communicator &communicator,
							 std::map<std::string, bool> log_param_scale,
                                                         std::map< std::string, NumericVector<Number>* >& training_parameters_in,
                                                         unsigned int n_training_samples_in,
                                                         const RBParameters& min_parameters,
                                                         const RBParameters& max_parameters,
                                                         bool serial_training_set=false);


  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * This boolean flag indicates whether or not the training set should
   * be the same on all processors. By default it is false, but in the
   * case of the Empirical Interpolation Method (RBEIMConstruction),
   * for example, we need the training set to be identical on all processors.
   */
  bool serial_training_set;

  /**
   * We keep an extra temporary vector that is useful for
   * performing inner products (avoids unnecessary memory
   * allocation/deallocation).
   */
  AutoPtr< NumericVector<Number> > inner_product_storage_vector;

  /**
   * Set this string to specify an alternative solver used in the set_alternative_solver()
   * function above.  Currently-supported values are:
   * .) unchanged, to continue using the default truth solve solver
   * .) amg, to use the BoomerAMG from Hypre (NOT for indefinite problems!)
   * .) mumps, to use a sparse direct solver
   * Note1: mumps and amg will only be available if PETSc has been compiled with them.
   * Note2: RBConstruction::process_parameters_file() is responsible for reading in this value ("rb_alternative_solver")
   *        from file for RBConstruction-derived subclasses
   * Note3: RBSCMSystem::process_parameters_file() reads this value ("scm_alternative_solver")
   *        for RBSCMSystem-derived subclasses
   */
  std::string alternative_solver;


private:

  /**
   * Boolean flag to indicate whether or not the
   * parameter ranges have been initialized.
   */
  bool training_parameters_initialized;

  /**
   * The training samples.
   */
  std::map< std::string, NumericVector<Number>* > training_parameters;

  /**
   * If < 0, use std::time() * processor_id() to seed the random
   * number generator for the training parameters (default).  If
   * >= 0, use the provided value * processor_id() as the random
   * number generator seed.
   */
  int training_parameters_random_seed;

  /**
   * The name of the parameter that we will generate a deterministic
   * training parameters for in the case of a "partially random" training
   * set.
   */
  std::string _deterministic_training_parameter_name;

  /**
   * The number of times each sample of the deterministic training parameter
   * is repeated in generating the training set.
   */
  unsigned int _deterministic_training_parameter_repeats;

};

} // namespace libMesh


#endif // LIBMESH_RB_CONSTRUCTION_BASE_H
