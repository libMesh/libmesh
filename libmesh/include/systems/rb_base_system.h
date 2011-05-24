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

#ifndef __rb_base_system_h__
#define __rb_base_system_h__

#include "system.h"
#include "numeric_vector.h"
#include "linear_solver.h"
#include "perf_log.h"

#include "rb_base.h"
#include "rb_theta.h"

#include <set>

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * This is the base class for certified reduced basis (RB)
 * systems.
 * We template the Base class so that we can derive from
 * the appropriate System type (e.g. LinearImplicitSystem
 * for standard reduced basis, EigenSystem for SCM)
 * at compile time.
 * We also inherit from RBBase to define a parameter domain
 * and an affine expansion of a PDE.
 *
 * @author David J. Knezevic, 2009
 */


// ------------------------------------------------------------
// RBBaseSystem class definition
template<class Base>
class RBBaseSystem : public Base, public RBBase
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  RBBaseSystem (EquationSystems& es,
                const std::string& name,
                const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBBaseSystem ();

  /**
   * The type of system.
   */
  typedef RBBaseSystem<Base> sys_type;

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
  unsigned int get_n_training_samples() const
    { libmesh_assert(training_parameters_initialized);
       return training_parameters[0]->size(); }

  /**
   * Get the total number of training samples local to this processor.
   */
  unsigned int get_local_n_training_samples() const
    { libmesh_assert(training_parameters_initialized);
      return training_parameters[0]->local_size(); }

  /**
   * Get the first local index of the training parameters.
   */
  unsigned int get_first_local_training_index() const
    { libmesh_assert(training_parameters_initialized);
      return training_parameters[0]->first_local_index(); }

  /**
   * Get the last local index of the training parameters.
   */
  unsigned int get_last_local_training_index() const
    { libmesh_assert(training_parameters_initialized);
      return training_parameters[0]->last_local_index(); }

  /**
   * return the specified training parameter.
   */
  std::vector<Real> get_training_parameter(unsigned int index) const;

  /**
   * Load the specified training parameter into current_params.
   * Local to this processor.
   */
  virtual void load_training_parameter_locally(unsigned int index);

  /**
   * Load the specified training parameter into current_params.
   * Loads the parameter on all processors.
   */
  virtual void load_training_parameter_globally(unsigned int index);

  /**
   * Initialize the parameter ranges and indicate whether deterministic
   * or random training parameters should be used and whether or
   * not we want the parameters to be scaled logarithmically.
   */
  virtual void initialize_training_parameters(const std::vector<Real>& mu_min_vector,
                                              const std::vector<Real>& mu_max_vector,
                                              const unsigned int n_training_parameters,
                                              const std::vector<bool> log_param_scale,
                                              const bool deterministic=true);

  /**
   * Overwrite the training parameters with new_training_set.
   */
  virtual void load_training_set(std::vector< std::vector<Number> >& new_training_set);

  /**
   * Changes the current PC (and iterative solver, if desired) in the
   * passed-in LinearSolver object to an alternative solver specified
   * by the alternative_solver string stored in this class.  You might
   * use this to e.g. switch to a sparse direct solver for the multiple
   * RHS solves executed during the update_residual_terms function.
   * The return strings are names of the original PC and KSP objects,
   * you can reset these using the reset_alternative_solver() function below.
   */
  std::pair<std::string,std::string>
  set_alternative_solver(AutoPtr<LinearSolver<Number> >& ls);

  /**
   * Resets the PC (and iterative solver, if desired) in the passed-in
   * LinearSolver object to the values specified in the pair of strings
   * passed as the second argument.  If the "alternative_solver" string,
   * defined below, is "unchanged", this function does nothing.
   */
  void reset_alternative_solver(AutoPtr<LinearSolver<Number> >& ls,
				const std::pair<std::string,std::string>& orig);
  
  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * Boolean flag to indicate whether or not the
   * parameter ranges have been initialized.
   */
  bool training_parameters_initialized;

  /**
   * Boolean flag to indicate whether or not we initialize
   * mesh dependent matrices and vectors when init_data
   * is called. Default value is true.
   */
  bool initialize_mesh_dependent_data;

  /**
   * The training samples.
   */
  std::vector< NumericVector<Number>* > training_parameters;

  /**
   * If < 0, use std::time() * processor_id() to seed the random
   * number generator for the training parameters (default).  If
   * >= 0, use the provided value * processor_id() as the random
   * number generator seed.
   */
  int training_parameters_random_seed;
  
protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Static function to return the error pair (index,error)
   * that is corresponds to the largest error on all
   * processors.
   */
  static void get_global_max_error_pair(std::pair<unsigned int, Real>& error_pair);

  /**
   * Static helper function for generating a randomized set of parameters.
   */
  static void generate_training_parameters_random(const std::vector<bool> log_param_scale,
                                                  std::vector< NumericVector<Number>* >& training_parameters_in,
                                                  const unsigned int n_training_samples_in,
                                                  const std::vector<Real>& min_parameters,
                                                  const std::vector<Real>& max_parameters,
						  int training_parameters_random_seed=-1,
						  bool serial_training_set=false);

  /**
   * Static helper function for generating a deterministic set of parameters. Only works with 1 or 2
   * parameters (as defined by the lengths of min/max parameters vectors), otherwise throws an error.
   */
  static void generate_training_parameters_deterministic(const std::vector<bool> log_param_scale,
                                                         std::vector< NumericVector<Number>* >& training_parameters_in,
                                                         const unsigned int n_training_samples_in,
                                                         const std::vector<Real>& min_parameters,
                                                         const std::vector<Real>& max_parameters,
						         bool serial_training_set=false);


  //----------- PROTECTED DATA MEMBERS -----------//
  
  /**
   * This boolean flag indicates whether or not the training set should
   * be the same on all processors. By default it is false, but in the
   * case of the Empirical Interpolation Method (RBEIMSystem), for example,
   * we need the training set to be identical on all processors.
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
   * Note2: RBSystem::init_data() is responsible for reading in this value ("rb_alternative_solver")
   *        from file for RBSystem-derived subclasses
   * Note3: RBSCMSystem::init_data() reads this value ("scm_alternative_solver")
   *        for RBSCMSystem-derived subclasses
   */
  std::string alternative_solver;

};

} // namespace libMesh


#endif
