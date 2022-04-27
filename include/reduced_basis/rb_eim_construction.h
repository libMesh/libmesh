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
#include "libmesh/rb_eim_evaluation.h"

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/coupling_matrix.h"

// C++ includes
#include <unordered_map>
#include <map>
#include <string>
#include <memory>
#include <vector>

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBEIMConstruction implements the Construction stage of the
 * Empirical Interpolation Method (EIM). This can be used to
 * create an approximation to parametrized functions. In the context
 * of the reduced basis (RB) method, the EIM approximation is typically
 * used to create an affine approximation to non-affine operators,
 * so that the standard RB method can be applied in that case.
 */
class RBEIMConstruction : public RBConstructionBase<System>
{
public:

  enum BEST_FIT_TYPE { PROJECTION_BEST_FIT, EIM_BEST_FIT, POD_BEST_FIT };

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBEIMConstruction (EquationSystems & es,
                     const std::string & name,
                     const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - Destructor is defaulted out-of-line
   */
  RBEIMConstruction (RBEIMConstruction &&) = default;
  RBEIMConstruction (const RBEIMConstruction &) = delete;
  RBEIMConstruction & operator= (const RBEIMConstruction &) = delete;
  RBEIMConstruction & operator= (RBEIMConstruction &&) = delete;
  virtual ~RBEIMConstruction ();

  /**
   * Type of the data structure used to map from (elem id) -> [n_vars][n_qp] data.
   */
  typedef RBEIMEvaluation::QpDataMap QpDataMap;

  /**
   * Type of the data structure used to map from (elem id,side_index) -> [n_vars][n_qp] data.
   */
  typedef RBEIMEvaluation::SideQpDataMap SideQpDataMap;

  /**
   * Type of the data structure used to map from node id -> [n_vars] data.
   */
  typedef RBEIMEvaluation::NodeDataMap NodeDataMap;

  /**
   * Clear this object.
   */
  virtual void clear() override;

  /**
   * Set the RBEIMEvaluation object.
   */
  void set_rb_eim_evaluation(RBEIMEvaluation & rb_eim_eval_in);

  /**
   * Get a reference to the RBEvaluation object.
   */
  RBEIMEvaluation & get_rb_eim_evaluation();

  /**
   * Get a const reference to the RBEvaluation object.
   */
  const RBEIMEvaluation & get_rb_eim_evaluation() const;

  /**
   * Perform initialization of this object to prepare for running
   * train_eim_approximation().
   */
  void initialize_eim_construction();

  /**
   * Read parameters in from file and set up this system
   * accordingly.
   */
  virtual void process_parameters_file (const std::string & parameters_filename);

  /**
   * Set the state of this RBConstruction object based on the arguments
   * to this function.
   */
  void set_rb_construction_parameters(unsigned int n_training_samples_in,
                                      bool deterministic_training_in,
                                      unsigned int training_parameters_random_seed_in,
                                      bool quiet_mode_in,
                                      unsigned int Nmax_in,
                                      Real rel_training_tolerance_in,
                                      Real abs_training_tolerance_in,
                                      RBParameters mu_min_in,
                                      RBParameters mu_max_in,
                                      std::map<std::string, std::vector<Real>> discrete_parameter_values_in,
                                      std::map<std::string,bool> log_scaling,
                                      std::map<std::string, std::vector<Number>> * training_sample_list=nullptr);

  /**
   * Specify which type of "best fit" we use to guide the EIM
   * greedy algorithm.
   */
  virtual void set_best_fit_type_flag (const std::string & best_fit_type_string);

  /**
   * Print out info that describes the current setup of this RBConstruction.
   */
  virtual void print_info();

  /**
   * Generate the EIM approximation for the specified parametrized function
   * using either POD or the Greedy Algorithm. Return the final tolerance.
   */
  virtual Real train_eim_approximation();

  /**
   * Generate the EIM approximation for the specified parametrized function
   * using the Greedy Algorithm. Return the final tolerance.
   * Method is virtual so that behavior can be specialized further in
   * subclasses, if needed.
   */
  virtual Real train_eim_approximation_with_greedy();

  /**
   * Generate the EIM approximation for the specified parametrized function
   * using Proper Orthogonal Decomposition (POD). Return the final tolerance.
   * Method is virtual so that behavior can be specialized further in
   * subclasses, if needed.
   */
  virtual Real train_eim_approximation_with_POD();

  /**
   * Build a vector of ElemAssembly objects that accesses the basis
   * functions stored in this RBEIMConstruction object. This is useful
   * for performing the Offline stage of the Reduced Basis method where
   * we want to use assembly functions based on this EIM approximation.
   */
  virtual void initialize_eim_assembly_objects();

  /**
   * \returns The vector of assembly objects that point to this RBEIMConstruction.
   */
  std::vector<std::unique_ptr<ElemAssembly>> & get_eim_assembly_objects();

  /**
   * Build an element assembly object that will access basis function
   * \p bf_index.
   * This is pure virtual, override in subclasses to specify the appropriate
   * ElemAssembly object.
   */
  virtual std::unique_ptr<ElemAssembly> build_eim_assembly(unsigned int bf_index) = 0;

  /**
   * Pre-request FE data needed for calculations.
   */
  virtual void init_context(FEMContext &);

  /**
   * Get/set the relative tolerance for the basis training.
   */
  void set_rel_training_tolerance(Real new_training_tolerance);
  Real get_rel_training_tolerance();

  /**
   * Get/set the absolute tolerance for the basis training.
   */
  void set_abs_training_tolerance(Real new_training_tolerance);
  Real get_abs_training_tolerance();

  /**
   * Get/set Nmax, the maximum number of RB
   * functions we are willing to compute.
   */
  unsigned int get_Nmax() const;
  virtual void set_Nmax(unsigned int Nmax);

  /**
   * Get the maximum value (across all processors) from
   * the parametrized functions in the training set.
   */
  Real get_max_abs_value_in_training_set() const;

  /**
   * Get the EIM solution vector at all parametrized functions in the training
   * set. In some cases we want to store this data for future use. For example
   * this is useful in the case that the parametrized function is defined
   * based on a look-up table rather than an analytical function, since
   * if we store the EIM solution data, we can do Online solves without
   * initializing the look-up table data.
   */
  void store_eim_solutions_for_training_set();

  /**
   * Get a const reference to the specified parametrized function from
   * the training set.
   */
  const QpDataMap & get_parametrized_function_from_training_set(unsigned int training_index) const;
  const SideQpDataMap & get_side_parametrized_function_from_training_set(unsigned int training_index) const;
  const NodeDataMap & get_node_parametrized_function_from_training_set(unsigned int training_index) const;

  /**
   * Get the interior and side quadrature weights.
   */
  const std::unordered_map<dof_id_type, std::vector<Real> > & get_local_quad_point_JxW();
  const std::map<std::pair<dof_id_type,unsigned int>, std::vector<Real> > & get_local_side_quad_point_JxW();

  /**
   * Get the number of parametrized functions used for training.
   */
  unsigned int get_n_parametrized_functions_for_training() const;

  /**
   * Zero the _eim_projection_matrix and resize it to be get_Nmax() x get_Nmax().
   */
  void reinit_eim_projection_matrix();

  /**
   * Enum that indicates which type of "best fit" algorithm
   * we should use.
   * a) projection: Find the best fit in the inner product
   * b) eim: Use empirical interpolation to find a "best fit"
   */
  BEST_FIT_TYPE best_fit_type_flag;

protected:

  /**
   * Implementation of enrich_eim_approximation() for the case of element sides.
   */
  void enrich_eim_approximation_on_sides(const SideQpDataMap & side_pf);

  /**
   * Implementation of enrich_eim_approximation() for the case of element nodes.
   */
  void enrich_eim_approximation_on_nodes(const NodeDataMap & node_pf);

  /**
   * Implementation of enrich_eim_approximation() for the case of element interiors.
   */
  void enrich_eim_approximation_on_interiors(const QpDataMap & interior_pf,
                                             const std::vector<std::vector<Number>> & interior_pf_obs_values);

  /**
   * Update the matrices used in training the EIM approximation.
   */
  void update_eim_matrices();

private:

  /**
   * Find the training sample that has the largest EIM approximation error
   * based on the current EIM approximation. Return the maximum error, and
   * the training sample index at which it occured.
   */
  std::pair<Real, unsigned int> compute_max_eim_error();

  /**
   * Compute and store the parametrized function for each
   * parameter in the training set at all the stored qp locations.
   */
  void initialize_parametrized_functions_in_training_set();

  /**
   * Initialize the data associated with each quad point (location, JxW, etc.)
   * so that we can use this in evaluation of the parametrized functions.
   */
  void initialize_qp_data();

  /**
   * Initialize the \p elem_ids and \p sbd_ids associated with the observation
   * points so that we can subsequently evaluate parametrized functions at the
   * observations points.
   */
  void initialize_observation_points_data(
    std::vector<dof_id_type> & observation_points_elem_ids,
    std::vector<subdomain_id_type> & observation_points_sbd_ids);

  /**
   * Evaluate the inner product of vec1 and vec2 which specify values at
   * quadrature points. The inner product includes the JxW contributions
   * stored in _local_quad_point_JxW, so that this is equivalent to
   * computing w^t M v, where M is the mass matrix.
   */
  Number inner_product(const QpDataMap & v, const QpDataMap & w);

  /**
   * Same as inner_product() except for side data.
   */
  Number side_inner_product(const SideQpDataMap & v, const SideQpDataMap & w);

  /**
   * Same as inner_product() except for node data.
   */
  Number node_inner_product(const NodeDataMap & v, const NodeDataMap & w);

  /**
   * Get the maximum absolute value from a vector stored in the format that we use
   * for basis functions.
   */
  template <class DataMap>
  Real get_max_abs_value(const DataMap & v) const
  {
    Real max_value = 0.;

    for (const auto & pr : v)
      {
        const auto & v_comp_and_qp = pr.second;

        for (const auto & comp : index_range(v_comp_and_qp))
          {
            // If scale_components_in_enrichment() returns true then we
            // apply a scaling to give an approximately uniform scaling
            // for all components.
            Real comp_scaling = 1.;
            if (get_rb_eim_evaluation().scale_components_in_enrichment())
              {
                // Make sure that _component_scaling_in_training_set is initialized
                libmesh_error_msg_if(comp >= _component_scaling_in_training_set.size(),
                                    "Invalid vector index");
                comp_scaling = _component_scaling_in_training_set[comp];
              }

            const std::vector<Number> & v_qp = v_comp_and_qp[comp];
            for (Number value : v_qp)
              max_value = std::max(max_value, std::abs(value * comp_scaling));
          }
      }

    comm().max(max_value);
    return max_value;
  }

  /**
   * Add a new basis function to the EIM approximation.
   */
  void enrich_eim_approximation(unsigned int training_index);

  /**
   * We compute the best fit of parametrized_function
   * into the EIM space and then evaluate the error
   * in the norm defined by inner_product_matrix.
   *
   * \returns The error in the best fit
   */
  Real compute_best_fit_error();

  /**
   * Scale all values in \p pf by \p scaling_factor
   */
  template<class DataMap>
  static void scale_parametrized_function(DataMap & local_pf,
                                          Number scaling_factor)
  {
    for (auto & pr : local_pf)
      {
        auto & comp_and_qp = pr.second;

        for (unsigned int comp : index_range(comp_and_qp))
          {
            std::vector<Number> & qp_values = comp_and_qp[comp];

            for (unsigned int qp : index_range(qp_values))
              {
                qp_values[qp] *= scaling_factor;
              }
          }
      }
  }

  /**
   * Maximum number of EIM basis functions we are willing to use.
   */
  unsigned int _Nmax;

  /**
   * Relative and absolute tolerances for training the EIM approximation.
   */
  Real _rel_training_tolerance;
  Real _abs_training_tolerance;

  /**
   * The matrix we use in order to perform L2 projections of
   * parametrized functions as part of EIM training.
   */
  DenseMatrix<Number> _eim_projection_matrix;

  /**
   * The RBEIMEvaluation object that we use to perform the EIM training.
   */
  RBEIMEvaluation * _rb_eim_eval;

  /**
   * The vector of assembly objects that are created to point to
   * this RBEIMConstruction.
   */
  std::vector<std::unique_ptr<ElemAssembly>> _rb_eim_assembly_objects;

  /**
   * The parametrized functions that are used for training. We pre-compute and
   * store all of these functions, rather than recompute them at each iteration
   * of the training.
   *
   * We store values at quadrature points on elements that are local to this processor.
   * The indexing is as follows:
   *   basis function index --> element ID --> variable --> quadrature point --> value
   * We use a map to index the element ID, since the IDs on this processor in
   * generally will not start at zero.
   */
  std::vector<QpDataMap> _local_parametrized_functions_for_training;

  /**
   * Same as _local_parametrized_functions_for_training except for side data.
   * The indexing is as follows:
   *   basis function index --> (element ID,side index) --> variable --> quadrature point --> value
   */
  std::vector<SideQpDataMap> _local_side_parametrized_functions_for_training;

  /**
   * Same as _local_parametrized_functions_for_training except for node data.
   * The indexing is as follows:
   *   basis function index --> node ID --> variable --> value
   */
  std::vector<NodeDataMap> _local_node_parametrized_functions_for_training;

  /**
   * Maximum value in _local_parametrized_functions_for_training across all processors.
   * This can be used for normalization purposes, for example.
   */
  Real _max_abs_value_in_training_set;

  /**
   * The training sample index at which we found _max_abs_value_in_training_set.
   */
  unsigned int _max_abs_value_in_training_set_index;

  /**
   * Keep track of a scaling factor for each component of the parametrized functions in
   * the training set which "scales up" each component to have a similar magnitude as
   * the largest component encountered in the training set. This can give more uniform
   * scaling across all components and is helpful in cases where components have widely
   * varying magnitudes.
   */
  std::vector<Real> _component_scaling_in_training_set;

  /**
   * The quadrature point locations, quadrature point weights (JxW), and subdomain IDs
   * on every element local to this processor.
   *
   * The indexing is as follows:
   *   element ID --> quadrature point --> xyz
   *   element ID --> quadrature point --> JxW
   *   element ID --> subdomain_id
   * We use a map to index the element ID, since the IDs on this processor in
   * generally will not start at zero.
   */
  std::unordered_map<dof_id_type, std::vector<Point> > _local_quad_point_locations;
  std::unordered_map<dof_id_type, std::vector<Real> > _local_quad_point_JxW;
  std::unordered_map<dof_id_type, subdomain_id_type > _local_quad_point_subdomain_ids;

  /**
   * EIM approximations often arise when applying a geometric mapping to a Reduced Basis
   * formulation. In this context, we often need to approximate derivates of the mapping
   * function via EIM. In order to enable this, we also optionally store perturbations
   * about each point in _local_quad_point_locations to enable finite difference approximation
   * to the mapping function derivatives.
   */
  std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > _local_quad_point_locations_perturbations;

  /**
   * Same as above except for side data.
   */
  std::map<std::pair<dof_id_type,unsigned int>, std::vector<Point> > _local_side_quad_point_locations;
  std::map<std::pair<dof_id_type,unsigned int>, std::vector<Real> > _local_side_quad_point_JxW;
  std::map<std::pair<dof_id_type,unsigned int>, subdomain_id_type > _local_side_quad_point_subdomain_ids;
  std::map<std::pair<dof_id_type,unsigned int>, boundary_id_type > _local_side_quad_point_boundary_ids;
  std::map<std::pair<dof_id_type,unsigned int>, std::vector<std::vector<Point>> > _local_side_quad_point_locations_perturbations;

  /**
   * Same as above except for node data.
   */
  std::unordered_map<dof_id_type, Point > _local_node_locations;
  std::unordered_map<dof_id_type, boundary_id_type > _local_node_boundary_ids;

  /**
   * For side data, we also store "side type" info. This is used to distinguish between
   * data that is stored on a "shellface" vs. a "standard side". The convention we use
   * here is:
   *  0 --> standard side
   *  1 --> shellface
   */
  std::map<std::pair<dof_id_type,unsigned int>, unsigned int > _local_side_quad_point_side_types;

  /**
   * We also optionally store the values at the "observation points" for all parametrized functions
   * in the training set. These values are used to obtain the observation values that are stored in
   * RBEIMEvaluation.
   *
   * Indexing is: training_index --> observation point index --> component --> value.
   */
  std::vector<std::vector<std::vector<Number>>> _parametrized_functions_for_training_obs_values;

};

} // namespace libMesh

#endif // LIBMESH_RB_EIM_CONSTRUCTION_H
