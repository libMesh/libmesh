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

#ifndef LIBMESH_RB_EIM_CONSTRUCTION_H
#define LIBMESH_RB_EIM_CONSTRUCTION_H

// rbOOmit includes
#include "libmesh/rb_construction.h"
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_eim_assembly.h"

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/coupling_matrix.h"

// C++ includes

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBEIMConstruction implements the Construction stage of the
 * Empirical Interpolation Method (EIM). This can be used to
 * generate an affine approximation to non-affine
 * operators.
 *
 * \author David J. Knezevic
 * \date 2010
 */
class RBEIMConstruction : public RBConstruction
{
public:

  enum BEST_FIT_TYPE { PROJECTION_BEST_FIT, EIM_BEST_FIT };

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBEIMConstruction (EquationSystems & es,
                     const std::string & name,
                     const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBEIMConstruction ();

  /**
   * The type of system.
   */
  typedef RBEIMConstruction sys_type;

  /**
   * The type of the parent.
   */
  typedef RBConstruction Parent;

  /**
   * Clear this object.
   */
  virtual void clear() libmesh_override;

  /**
   * Read parameters in from file and set up this system
   * accordingly.
   */
  virtual void process_parameters_file (const std::string & parameters_filename) libmesh_override;

  /**
   * Specify which type of "best fit" we use to guide the EIM
   * greedy algorithm.
   */
  void set_best_fit_type_flag (const std::string & best_fit_type_string);

  /**
   * Print out info that describes the current setup of this RBConstruction.
   */
  virtual void print_info() libmesh_override;

  /**
   * Initialize this system so that we can perform
   * the Construction stage of the RB method.
   */
  virtual void initialize_rb_construction(bool skip_matrix_assembly=false,
                                          bool skip_vector_assembly=false) libmesh_override;

  /**
   * Override train_reduced_basis to first initialize _parametrized_functions_in_training_set.
   */
  virtual Real train_reduced_basis(const bool resize_rb_eval_data=true) libmesh_override;

  /**
   * Load the truth representation of the parametrized function
   * at the current parameters into the solution vector.
   * The truth representation is the projection of
   * parametrized_function into the finite element space.
   * If \p plot_solution > 0 the solution will be plotted
   * to an output file.
   */
  virtual Real truth_solve(int plot_solution) libmesh_override;

  /**
   * We compute the best fit of parametrized_function
   * into the EIM space and then evaluate the error
   * in the norm defined by inner_product_matrix.
   *
   * @return the error in the best fit
   */
  virtual Real compute_best_fit_error();

  /**
   * Initialize \p c based on \p sys.
   */
  virtual void init_context_with_sys(FEMContext & c, System & sys);

  /**
   * Add variables to the ExplicitSystem that is used to store
   * the basis functions.
   */
  virtual void init_explicit_system() = 0;

  /**
   * Add one variable to the ImplicitSystem (i.e. this system) that is
   * used to perform L2 project solves.
   */
  virtual void init_implicit_system() = 0;

  /**
   * Evaluate the mesh function at the specified point and for the specified variable.
   */
  Number evaluate_mesh_function(unsigned int var_number,
                                Point p);

  /**
   * Build a vector of ElemAssembly objects that accesses the basis
   * functions stored in this RBEIMConstruction object. This is useful
   * for performing the Offline stage of the Reduced Basis method where
   * we want to use assembly functions based on this EIM approximation.
   */
  virtual void initialize_eim_assembly_objects();

  /**
   * @return the vector of assembly objects that point to this RBEIMConstruction.
   */
  std::vector<ElemAssembly *> get_eim_assembly_objects();

  /**
   * Build an element assembly object that will access basis function
   * \p bf_index.
   * This is pure virtual, override in subclasses to specify the appropriate
   * ElemAssembly object.
   */
  virtual UniquePtr<ElemAssembly> build_eim_assembly(unsigned int bf_index) = 0;

  /**
   * Get the ExplicitSystem associated with this system.
   */
  ExplicitSystem& get_explicit_system();

  /**
   * Load the i^th RB function into the RBConstruction
   * solution vector.
   * Override to load the basis function into the ExplicitSystem.
   */
  virtual void load_basis_function(unsigned int i) libmesh_override;

  /**
   * Load the RB solution from the most recent solve with rb_eval
   * into this system's solution vector.
   * Override to load the solution into the ExplicitSystem.
   */
  virtual void load_rb_solution() libmesh_override;

  /**
   * Load \p source into the subvector of \p dest corresponding
   * to var \p var.
   */
  void set_explicit_sys_subvector(
    NumericVector<Number>& dest, unsigned int var, NumericVector<Number>& source);

  /**
   * Load the subvector of \p localized_source corresponding to variable \p var into
   * \p dest. We require localized_source to be localized before we call this method.
   */
  void get_explicit_sys_subvector(
    NumericVector<Number>& dest, unsigned int var, NumericVector<Number>& localized_source);

  /**
   * Set up the index map between the implicit and explicit systems.
   */
  void init_dof_map_between_systems();

  /**
   * Plot all the parameterized functions that we are storing
   * in _parametrized_functions_in_training_set. \p pathname
   * provides the path to where the plot data will be saved.
   */
  void plot_parametrized_functions_in_training_set(const std::string& pathname);

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * Enum that indicates which type of "best fit" algorithm
   * we should use.
   * a) projection: Find the best fit in the inner product
   * b) eim: Use empirical interpolation to find a "best fit"
   *
   * @return the error associated with the "best fit" in the
   * norm induced by inner_product_matrix.
   */
  BEST_FIT_TYPE best_fit_type_flag;

protected:

  /**
   * Override to initialize the coupling matrix to decouple variables in this system.
   */
  virtual void init_data() libmesh_override;

  /**
   * Add a new basis function to the RB space. Overload
   * to enrich with the EIM basis functions.
   */
  virtual void enrich_RB_space() libmesh_override;

  /**
   * Update the system after enriching the RB space; this calls
   * a series of functions to update the system properly.
   */
  virtual void update_system() libmesh_override;

  /**
   * Compute the reduced basis matrices for the current basis.
   * Overload to update the inner product matrix that
   * is used to compute the best fit to parametrized_function.
   */
  virtual void update_RB_system_matrices() libmesh_override;

  /**
   * Overload to return the best fit error. This function is used in
   * the Greedy algorithm to select the next parameter.
   */
  virtual Real get_RB_error_bound() libmesh_override;

  /**
   * Function that indicates when to terminate the Greedy
   * basis training. Overload in subclasses to specialize.
   */
  virtual bool greedy_termination_test(
    Real abs_greedy_error, Real initial_greedy_error, int count) libmesh_override;

  /**
   * Loop over the training set and compute the parametrized function for each
   * training index.
   */
  void initialize_parametrized_functions_in_training_set();

  /**
   * Boolean flag to indicate whether or not we have called
   * compute_parametrized_functions_in_training_set() yet.
   */
  bool _parametrized_functions_in_training_set_initialized;

  /**
   * The libMesh vectors storing the finite element coefficients
   * of the RB basis functions.
   */
  std::vector< NumericVector<Number> * > _parametrized_functions_in_training_set;

private:

  /**
   * A mesh function to interpolate on the mesh.
   */
  MeshFunction * _mesh_function;

  /**
   * This flag indicates that we're in the process of
   * performing one extra Greedy step in order to compute
   * the data needed for the EIM a posteriori error bound
   * in the case that we use all of our basis functions.
   */
  bool _performing_extra_greedy_step;

  /**
   * We also need an extra vector in which we can store a ghosted
   * copy of the vector that we wish to use MeshFunction on.
   */
  UniquePtr< NumericVector<Number> > _ghosted_meshfunction_vector;

  /**
   * We initialize RBEIMConstruction so that it has an "empty" RBAssemblyExpansion,
   * because this isn't used at all in the EIM.
   */
  RBAssemblyExpansion _empty_rb_assembly_expansion;

  /**
   * The vector of assembly objects that are created to point to
   * this RBEIMConstruction.
   */
  std::vector<ElemAssembly *> _rb_eim_assembly_objects;

  /**
   * We use an ExplicitSystem to store the EIM basis functions.
   * This is because if we have an EIM system with many variables
   * we need to allocate a large matrix. Better to avoid this
   * and use an ExplicitSystem instead, and only use the ImplicitSystem
   * to deal with the per-variable L2 projections.
   */
  std::string _explicit_system_name;

  /**
   * The index map between the explicit system and the implicit system.
   */
  std::vector< std::vector<dof_id_type> > _dof_map_between_systems;

  /**
   * This vector is used to store inner_product_matrix * basis_function[i] for each i,
   * since we frequently use this data.
   */
  std::vector< NumericVector<Number>* > _matrix_times_bfs;

};

} // namespace libMesh

#endif // LIBMESH_RB_EIM_CONSTRUCTION_H
