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
 * @author David J. Knezevic, 2010
 */

// ------------------------------------------------------------
// RBEIMConstruction class definition

class RBEIMConstruction : public RBConstruction
{
public:

  enum BEST_FIT_TYPE { PROJECTION_BEST_FIT, EIM_BEST_FIT };

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBEIMConstruction (EquationSystems& es,
                     const std::string& name,
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
  virtual void clear();

  /**
   * Read parameters in from file and set up this system
   * accordingly.
   */
  virtual void process_parameters_file (const std::string& parameters_filename);

  /**
   * Specify which type of "best fit" we use to guide the EIM
   * greedy algorithm.
   */
  void set_best_fit_type_flag (const std::string& best_fit_type_string);

  /**
   * Print out info that describes the current setup of this RBConstruction.
   */
  virtual void print_info();

  /**
   * Initialize this system so that we can perform
   * the Construction stage of the RB method.
   */
  virtual void initialize_rb_construction();

  /**
   * Override train_reduced_basis to first initialize _parametrized_functions_in_training_set.
   */
  virtual Real train_reduced_basis(const std::string& directory_name = "offline_data",
                                   const bool resize_rb_eval_data=true);

  /**
   * Load the truth representation of the parametrized function
   * at the current parameters into the solution vector.
   * The truth representation is the projection of
   * parametrized_function into the finite element space.
   * If \p plot_solution > 0 the solution will be plotted
   * to an output file.
   */
  virtual Real truth_solve(int plot_solution);

  /**
   * We compute the best fit of parametrized_function
   * into the EIM space and then evaluate the error
   * in the norm defined by inner_product_matrix.
   *
   * @return the error in the best fit
   */
  virtual Real compute_best_fit_error();

  /**
   * Provide an implementation of init_context that is
   * relevant to the projection calculations in
   * load_calN_parametrized_function.
   */
  virtual void init_context(FEMContext &c);

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
  std::vector<ElemAssembly*> get_eim_assembly_objects();

  /**
   * Build an element assembly object that will access basis function
   * \p bf_index.
   * This is pure virtual, override in subclasses to specify the appropriate
   * ElemAssembly object.
   */
  virtual AutoPtr<ElemAssembly> build_eim_assembly(unsigned int bf_index) = 0;

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


  /**
   * The matrix associated with this system should be block diagonal
   * hence we use a CouplingMatrix to "decouple" all variables to
   * save memory. That is, _coupling_matrix is diagonal.
   */
  CouplingMatrix _coupling_matrix;

protected:

  /**
   * Override to initialize the coupling matrix to decouple variables in this system.
   */
  virtual void init_data();

  /**
   * Add a new basis function to the RB space. Overload
   * to enrich with the EIM basis functions.
   */
  virtual void enrich_RB_space();

  /**
   * Update the system after enriching the RB space; this calls
   * a series of functions to update the system properly.
   */
  virtual void update_system();

  /**
   * Compute the reduced basis matrices for the current basis.
   * Overload to update the inner product matrix that
   * is used to compute the best fit to parametrized_function.
   */
  virtual void update_RB_system_matrices();

  /**
   * Overload to return the best fit error. This function is used in
   * the Greedy algorithm to select the next parameter.
   */
  virtual Real get_RB_error_bound();

  /**
   * Function that indicates when to terminate the Greedy
   * basis training. Overload in subclasses to specialize.
   */
  virtual bool greedy_termination_test(Real training_greedy_error, int count);

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
  std::vector< NumericVector<Number>* > _parametrized_functions_in_training_set;

private:

  /**
   * A mesh function to interpolate on the mesh.
   */
  MeshFunction* _mesh_function;

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
  AutoPtr< NumericVector<Number> > _ghosted_meshfunction_vector;

  /**
   * We initialize RBEIMConstruction so that it has an "empty" RBAssemblyExpansion,
   * because this isn't used at all in the EIM.
   */
  RBAssemblyExpansion _empty_rb_assembly_expansion;

  /**
   * The vector of assembly objects that are created to point to
   * this RBEIMConstruction.
   */
  std::vector<ElemAssembly*> _rb_eim_assembly_objects;

};

} // namespace libMesh

#endif // LIBMESH_RB_EIM_CONSTRUCTION_H
