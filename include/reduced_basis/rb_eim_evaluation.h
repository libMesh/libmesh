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

#ifndef LIBMESH_RB_EIM_EVALUATION_H
#define LIBMESH_RB_EIM_EVALUATION_H

// libMesh includes
#include "libmesh/point.h"
#include "libmesh/rb_theta_expansion.h"
#include "libmesh/rb_parametrized.h"
#include "libmesh/parallel_object.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/fe_type.h"

// C++ includes
#include <memory>
#include <map>
#include <vector>
#include <string>

namespace libMesh
{

class RBParameters;
class RBParametrizedFunction;
class RBTheta;
class System;
class EquationSystems;
class Elem;

/**
 * This struct encapsulates data that specifies how we will
 * perform plotting for EIM variable groups.
 */
struct EIMVarGroupPlottingInfo
{

  /**
   * Default constructor. Initialize values that would otherwise
   * be uninitialized. We do not initialize the FEType or string
   * members since they already have default constructors.
   */
  EIMVarGroupPlottingInfo()
  :
  first_eim_var_index(0),
  n_eim_vars(0),
  enforce_min_value(false),
  enforce_max_value(false),
  min_value(0.),
  max_value(0.)
  {
  }

  /**
   * The index for the first EIM variable in this variable group.
   */
  unsigned int first_eim_var_index;

  /**
   * The number of EIM variables in the group. The variables
   * are assumed to be numbered contiguously.
   */
  unsigned int n_eim_vars;

  /**
   * The name of the System we use for plotting this EIM variable group.
   */
  std::string eim_sys_name;

  /**
   * These booleans indicate if we should clamp the resulting output
   * to be above a min value or below a max value. This can be relevant
   * if we want to satisfy some physical constraints on the outputs, for
   * example, since these constraints may not be exactly satisfied by
   * the EIM output.
   */
  bool enforce_min_value;
  bool enforce_max_value;

  /**
   * The min (resp. max) value that we enforce if enforce_min_value
   * (resp. enforce_max_value) is true.
   */
  Real min_value;
  Real max_value;
};

/**
 * This class enables evaluation of an Empirical Interpolation Method (EIM)
 * approximation. RBEvaluation plays an analogous role in the context of
 * the regular reduced basis method.
 */
class RBEIMEvaluation : public RBParametrized,
                        public ParallelObject
{
public:

  /**
   * Constructor.
   */
  RBEIMEvaluation(const Parallel::Communicator & comm);

  /**
   * Special functions.
   * - This class contains unique_ptrs, so it can't be default copy
       constructed/assigned.
   * - The destructor is defaulted out of line.
   */
  RBEIMEvaluation (RBEIMEvaluation &&) = default;
  RBEIMEvaluation (const RBEIMEvaluation &) = delete;
  RBEIMEvaluation & operator= (const RBEIMEvaluation &) = delete;
  RBEIMEvaluation & operator= (RBEIMEvaluation &&) = default;
  virtual ~RBEIMEvaluation ();

  /**
   * Type of the data structure used to map from (elem id) -> [n_vars][n_qp] data.
   */
  typedef std::map<dof_id_type, std::vector<std::vector<Number>>> QpDataMap;

  /**
   * Type of the data structure used to map from (elem id, side index) -> [n_vars][n_qp] data.
   */
  typedef std::map<std::pair<dof_id_type,unsigned int>, std::vector<std::vector<Number>>> SideQpDataMap;

  /**
   * Type of the data structure used to map from (node id) -> [n_vars] data.
   */
  typedef std::map<dof_id_type, std::vector<Number>> NodeDataMap;

  /**
   * Clear this object.
   */
  virtual void clear() override;

  /**
   * Resize the data structures for storing data associated
   * with this object.
   */
  void resize_data_structures(const unsigned int Nmax);

  /**
   * Set the parametrized function that we will approximate
   * using the Empirical Interpolation Method. This object
   * will take ownership of the unique pointer.
   */
  void set_parametrized_function(std::unique_ptr<RBParametrizedFunction> pf);

  /**
   * Get a reference to the parametrized function.
   */
  RBParametrizedFunction & get_parametrized_function();

  /**
   * Get a const reference to the parametrized function.
   */
  const RBParametrizedFunction & get_parametrized_function() const;

  /**
   * Calculate the EIM approximation for the given
   * right-hand side vector \p EIM_rhs. Store the
   * solution coefficients in the member _eim_solution.
   */
  DenseVector<Number> rb_eim_solve(DenseVector<Number> & EIM_rhs);

  /**
   * Perform rb_eim_solves at each mu in \p mus and store the results
   * in _rb_eim_solutions.
   */
  void rb_eim_solves(const std::vector<RBParameters> & mus, unsigned int N);

  /**
   * Initialize _interpolation_points_spatial_indices. Once this data is
   * initialized, we can store it in the training data, and read it back
   * in during the Online stage to be used in solves.
   */
  void initialize_interpolation_points_spatial_indices();

  /**
   * The Online counterpart of initialize_interpolation_points_spatial_indices().
   * This is used to initialize the spatial indices data in _parametrized_function
   * so that the _parametrized_function can be used in the Online stage without
   * reconstructing the spatial indices on every element or node in the mesh.
   */
  void initialize_param_fn_spatial_indices();

  /**
   * Return the current number of EIM basis functions.
   */
  unsigned int get_n_basis_functions() const;

  /**
   * Return the number of interpolation points. If we're not using the EIM error
   * indicator, then this matches get_n_basis_functions(), but if we are using
   * the EIM error indicator then we should have one extra interpolation point.
   */
  unsigned int get_n_interpolation_points() const;

  /**
   * Return the number of unique elements containing interpolation points
   */
  unsigned int get_n_elems() const;

  /**
   * Set the number of basis functions. Useful when reading in
   * stored data.
   */
  void set_n_basis_functions(unsigned int n_bfs);

  /**
   * Subtract coeffs[i]*basis_function[i] from \p v.
   */
  void decrement_vector(QpDataMap & v,
                        const DenseVector<Number> & coeffs);

  /**
   * Same as decrement_vector() except for Side data.
   */
  void side_decrement_vector(SideQpDataMap & v,
                             const DenseVector<Number> & coeffs);

  /**
   * Same as decrement_vector() except for node data.
   */
  void node_decrement_vector(NodeDataMap & v,
                             const DenseVector<Number> & coeffs);

  /**
   * Build a vector of RBTheta objects that accesses the components
   * of the RB_solution member variable of this RBEvaluation.
   * Store these objects in the member vector rb_theta_objects.
   */
  void initialize_eim_theta_objects();

  /**
   * \returns The vector of theta objects that point to this RBEIMEvaluation.
   */
  std::vector<std::unique_ptr<RBTheta>> & get_eim_theta_objects();

  /**
   * Build a theta object corresponding to EIM index \p index.
   * The default implementation builds an RBEIMTheta object, possibly
   * override in subclasses if we need more specialized behavior.
   */
  virtual std::unique_ptr<RBTheta> build_eim_theta(unsigned int index);

  /**
   * Fill up values by evaluating the parametrized function \p pf for all quadrature
   * points on element \p elem_id and component \p comp.
   */
  static void get_parametrized_function_values_at_qps(
    const QpDataMap & pf,
    dof_id_type elem_id,
    unsigned int comp,
    std::vector<Number> & values);

  /**
   * Same as get_parametrized_function_values_at_qps() except for side data.
   */
  static void get_parametrized_function_side_values_at_qps(
    const SideQpDataMap & pf,
    dof_id_type elem_id,
    unsigned int side_index,
    unsigned int comp,
    std::vector<Number> & values);

  /**
   * Same as get_parametrized_function_values_at_qps() except for node data.
   * Note that this does not do any parallel communication, so it is only
   * applicable to looking up local values.
   */
  static Number get_parametrized_function_node_local_value(
    const NodeDataMap & pf,
    dof_id_type node_id,
    unsigned int comp);

  /**
   * Same as above, except that we just return the value at the qp^th
   * quadrature point.
   */
  static Number get_parametrized_function_value(
    const Parallel::Communicator & comm,
    const QpDataMap & pf,
    dof_id_type elem_id,
    unsigned int comp,
    unsigned int qp);

  /**
   * Same as get_parametrized_function_value() except for side data.
   */
  static Number get_parametrized_function_side_value(
    const Parallel::Communicator & comm,
    const SideQpDataMap & pf,
    dof_id_type elem_id,
    unsigned int side_index,
    unsigned int comp,
    unsigned int qp);

  /**
   * Same as get_parametrized_function_value() except for node data.
   * Unlike get_parametrized_function_node_local_value(), this does
   * parallel communication, and therefore if can be used to look up
   * values regardless of whether or not \p node_id is local.
   */
  static Number get_parametrized_function_node_value(
    const Parallel::Communicator & comm,
    const NodeDataMap & pf,
    dof_id_type node_id,
    unsigned int comp);

  /**
   * Fill up \p values with the basis function values for basis function
   * \p basis_function_index and variable \p var, at all quadrature points
   * on element \p elem_id. Each processor stores data for only the
   * elements local to that processor, so if elem_id is not on this processor
   * then \p values will be empty.
   */
  void get_eim_basis_function_values_at_qps(unsigned int basis_function_index,
                                            dof_id_type elem_id,
                                            unsigned int var,
                                            std::vector<Number> & values) const;

  /**
   * Same as get_eim_basis_function_values_at_qps() except for side data.
   */
  void get_eim_basis_function_side_values_at_qps(unsigned int basis_function_index,
                                                 dof_id_type elem_id,
                                                 unsigned int side_index,
                                                 unsigned int var,
                                                 std::vector<Number> & values) const;

  /**
   * Same as get_eim_basis_function_values_at_qps() except for node data.
   * Note that this does not do any parallel communication, it just looks
   * up the value from _local_node_eim_basis_functions.
   */
  Number get_eim_basis_function_node_local_value(unsigned int basis_function_index,
                                                 dof_id_type node_id,
                                                 unsigned int var) const;

  /**
   * Same as above, except that we just return the value at the qp^th
   * quadrature point.
   */
  Number get_eim_basis_function_value(unsigned int basis_function_index,
                                      dof_id_type elem_id,
                                      unsigned int comp,
                                      unsigned int qp) const;

  /**
   * Same as get_eim_basis_function_value() except for side data.
   */
  Number get_eim_basis_function_side_value(unsigned int basis_function_index,
                                          dof_id_type elem_id,
                                          unsigned int side_index,
                                          unsigned int comp,
                                          unsigned int qp) const;

  /**
   * Same as get_eim_basis_function_value() except for node data.
   * Note that unlike get_eim_basis_function_node_local_value(),
   * this does do parallel communication so that it can be called
   * on any processor regardless of whether \p node_id is local or not.
   */
  Number get_eim_basis_function_node_value(unsigned int basis_function_index,
                                           dof_id_type node_id,
                                           unsigned int var) const;

  /**
   * Get a reference to the i^th basis function.
   */
  const QpDataMap & get_basis_function(unsigned int i) const;

  /**
   * Get a reference to the i^th side basis function.
   */
  const SideQpDataMap & get_side_basis_function(unsigned int i) const;

  /**
   * Get a reference to the i^th node basis function.
   */
  const NodeDataMap & get_node_basis_function(unsigned int i) const;

  /**
   * Set _rb_eim_solutions. Normally we update _rb_eim_solutions by performing
   * and EIM solve, but in some cases we want to set the EIM solution coefficients
   * elsewhere, so this setter enables us to do that.
   */
  void set_rb_eim_solutions(const std::vector<DenseVector<Number>> & rb_eim_solutions);

  /**
   * Return the EIM solution coefficients from the most recent call to rb_eim_solves().
   */
  const std::vector<DenseVector<Number>> & get_rb_eim_solutions() const;

  /**
   * Return entry \p index for each solution in _rb_eim_solutions.
   */
  std::vector<Number> get_rb_eim_solutions_entries(unsigned int index) const;

  /**
   * Return a const reference to the EIM solutions for the parameters in the training set.
   */
  const std::vector<DenseVector<Number>> & get_eim_solutions_for_training_set() const;

  /**
   * Return a writeable reference to the EIM solutions for the parameters in the training set.
   */
  std::vector<DenseVector<Number>> & get_eim_solutions_for_training_set();

  /**
   * Return the EIM error indicator values from the most recent call to rb_eim_solves().
   * The first entry of each pair is the normalized error indicator, and the second
   * entry is the normalization factor.
   */
  const std::vector<std::pair<Real,Real>> & get_rb_eim_error_indicators() const;

  /**
   * Set the data associated with EIM interpolation points.
   */
  void add_interpolation_points_xyz(Point p);
  void add_interpolation_points_comp(unsigned int comp);
  void add_interpolation_points_subdomain_id(subdomain_id_type sbd_id);
  void add_interpolation_points_boundary_id(boundary_id_type b_id);
  void add_interpolation_points_xyz_perturbations(const std::vector<Point> & perturbs);
  void add_interpolation_points_elem_id(dof_id_type elem_id);
  void add_interpolation_points_side_index(unsigned int side_index);
  void add_interpolation_points_node_id(dof_id_type node_id);
  void add_interpolation_points_qp(unsigned int qp);
  void add_interpolation_points_elem_type(ElemType elem_type);
  void add_interpolation_points_phi_i_qp(const std::vector<Real> & phi_i_qp);
  void add_interpolation_points_JxW_all_qp(const std::vector<Real> & JxW_all_qp);
  void add_interpolation_points_phi_i_all_qp(const std::vector<std::vector<Real>> & phi_i_all_qp);
  void add_interpolation_points_qrule_order(Order qrule_order);
  void add_elem_center_dxyzdxi(const Point & dxyzdxi);
  void add_elem_center_dxyzdeta(const Point & dxyzdxi);
  void add_interpolation_points_spatial_indices(const std::vector<unsigned int> & spatial_indices);
  void add_elem_id_local_index_map_entry(const dof_id_type & elem_id, const unsigned int local_index);

  /**
   * Get the data associated with EIM interpolation points.
   */
  Point get_interpolation_points_xyz(unsigned int index) const;
  unsigned int get_interpolation_points_comp(unsigned int index) const;
  subdomain_id_type get_interpolation_points_subdomain_id(unsigned int index) const;
  boundary_id_type get_interpolation_points_boundary_id(unsigned int index) const;
  const std::vector<Point> & get_interpolation_points_xyz_perturbations(unsigned int index) const;
  dof_id_type get_interpolation_points_elem_id(unsigned int index) const;
  unsigned int get_interpolation_points_side_index(unsigned int index) const;
  dof_id_type get_interpolation_points_node_id(unsigned int index) const;
  unsigned int get_interpolation_points_qp(unsigned int index) const;
  ElemType get_interpolation_points_elem_type(unsigned int index) const;
  const std::vector<Real> & get_interpolation_points_phi_i_qp(unsigned int index) const;
  const std::vector<Real> & get_interpolation_points_JxW_all_qp(unsigned int index) const;
  const std::vector<std::vector<Real>> & get_interpolation_points_phi_i_all_qp(unsigned int index) const;
  Order get_interpolation_points_qrule_order(unsigned int index) const;
  const Point & get_elem_center_dxyzdxi(unsigned int index) const;
  const Point & get_elem_center_dxyzdeta(unsigned int index) const;
  const std::vector<unsigned int> & get_interpolation_points_spatial_indices(unsigned int index) const;
  const std::map<dof_id_type, unsigned int> & get_elem_id_to_local_index_map() const;

  /**
   * _interpolation_points_spatial_indices is optional data, so we need to be able to
   * check how many _interpolation_points_spatial_indices values have actually been set
   * since it may not match the number of interpolation points.
   */
  unsigned int get_n_interpolation_points_spatial_indices() const;

  /**
   * Set entry of the EIM interpolation matrix.
   */
  void set_interpolation_matrix_entry(unsigned int i, unsigned int j, Number value);

  /**
   * Get the EIM interpolation matrix.
   */
  const DenseMatrix<Number> & get_interpolation_matrix() const;

  /**
   * Add \p bf to our EIM basis.
   */
  void add_basis_function(
    const QpDataMap & bf);

  /**
   * Add interpolation data associated with a new basis function.
   */
  void add_interpolation_data(
    Point p,
    unsigned int comp,
    dof_id_type elem_id,
    subdomain_id_type subdomain_id,
    unsigned int qp,
    const std::vector<Point> & perturbs,
    const std::vector<Real> & phi_i_qp,
    ElemType elem_type,
    const std::vector<Real> & JxW_all_qp,
    const std::vector<std::vector<Real>> & phi_i_all_qp,
    Order qrule_order,
    const Point & dxyz_dxi_elem_center,
    const Point & dxyz_deta_elem_center);

  /**
   * Add \p side_bf to our EIM basis.
   */
  void add_side_basis_function(
    const SideQpDataMap & side_bf);

  /**
   * Add interpolation data associated with a new basis function.
   */
  void add_side_interpolation_data(
    Point p,
    unsigned int comp,
    dof_id_type elem_id,
    unsigned int side_index,
    subdomain_id_type subdomain_id,
    boundary_id_type boundary_id,
    unsigned int qp,
    const std::vector<Point> & perturbs,
    const std::vector<Real> & phi_i_qp);

  /**
   * Add \p node_bf to our EIM basis.
   */
  void add_node_basis_function(
    const NodeDataMap & node_bf);

  /**
   * Add interpolation data associated with a new basis function.
   */
  void add_node_interpolation_data(
    Point p,
    unsigned int comp,
    dof_id_type node_id,
    boundary_id_type boundary_id);

  /**
   * Set _preserve_rb_eim_solutions.
   */
  void set_preserve_rb_eim_solutions(bool preserve_rb_eim_solutions);

  /**
   * Get _preserve_rb_eim_solutions.
   */
  bool get_preserve_rb_eim_solutions() const;

  /**
   * Write out all the basis functions to file.
   * \p sys is used for file IO
   * \p directory_name specifies which directory to write files to
   * \p read_binary_basis_functions indicates whether to write
   * binary or ASCII data
   *
   * Note: this is not currently a virtual function and is not related
   * to the RBEvaluation function of the same name.
   */
  void write_out_basis_functions(const std::string & directory_name = "offline_data",
                                 bool write_binary_basis_functions = true);

  /**
   * Read in all the basis functions from file.
   *
   * \param sys The Mesh in this System determines the parallel distribution of the basis functions.
   * \param directory_name Specifies which directory to write files to.
   * \param read_binary_basis_functions Indicates whether to expect binary or ASCII data.
   *
   * Note: this is not a virtual function and is not related to the
   * RBEvaluation function of the same name.
   */
  void read_in_basis_functions(const System & sys,
                               const std::string & directory_name = "offline_data",
                               bool read_binary_basis_functions = true);

  /**
   * Project the EIM basis function data stored in \p bf_data onto
   * sys.solution. The intent of this function is to work with the
   * data format provided by get_interior_basis_functions_as_vecs().
   * That format can be easily serialized, if needed, and hence can
   * be used to provide \p bf_data after reading in data from disk,
   * for example.
   *
   * \p extra_options can be used to pass extra information to this
   * method, e.g. options related to how to perform the projection.
   *
   * This is a no-op by default, implement in sub-classes if needed.
   */
  virtual void project_qp_data_vector_onto_system(System & sys,
                                                  const std::vector<Number> & bf_data,
                                                  const EIMVarGroupPlottingInfo & eim_vargroup,
                                                  const std::map<std::string,std::string> & extra_options);

  /**
   * Get _eim_vars_to_project_and_write.
   */
  const std::vector<EIMVarGroupPlottingInfo> & get_eim_vars_to_project_and_write() const;

  /**
   * Get _scale_components_in_enrichment.
   */
  const std::set<unsigned int> & scale_components_in_enrichment() const;

  /**
   * Virtual function to indicate if we use the EIM error indicator in this case.
   * This indicates if we will generate the data during the Offline training for
   * the EIM error indicator. In EIM solves, the error indicator will only
   * be used if set_eim_error_indicator_active() is set to true, since we want
   * to be able to enable or disable the error indicator depending on the type
   * of solve we are doing (e.g. EIM solves during training do not need the
   * error indicator).
   */
  virtual bool use_eim_error_indicator() const;

  /**
   * Activate/decative the error indicator in EIM solves. We need this option since
   * in some cases (e.g. during EIM training) we do not want to activate the EIM
   * error indicator, whereas in "online solves" we do want to activate it.
   */
  void set_eim_error_indicator_active(bool is_active);

  /**
   * Get/set _extra_points_interpolation_matrix.
   */
  const DenseVector<Number> & get_error_indicator_interpolation_row() const;
  void set_error_indicator_interpolation_row(const DenseVector<Number> & error_indicator_row);

  /**
   * Evaluates the EIM error indicator based on \p error_indicator_rhs, \p eim_solution,
   * and _error_indicator_interpolation_row. We also pass in \p eim_rhs since this is
   * used to normalize the error indicator.
   *
   * We return a pair that specifies the relative error indicator, and the normalization
   * that was used to compute the relative error indicator. We can then recover the
   * absolute error indicator via rel. indicator x normalization.
   */
  std::pair<Real,Real> get_eim_error_indicator(
    Number error_indicator_rhs,
    const DenseVector<Number> & eim_solution,
    const DenseVector<Number> & eim_rhs);

  /**
   * Get the VectorizedEvalInput data.
   */
  const VectorizedEvalInput & get_vec_eval_input() const;

  /**
   * Get all interior basis functions in the form of std::vectors.
   * This can provide a convenient format for processing the basis
   * function data, or writing it to disk.
   *
   * Indexing is as follows:
   *  basis function index --> variable index --> data at qps per elems
   *
   * In parallel we gather basis functions to processor 0 and hence
   * we only return non-empty data on processor 0.
   */
  std::vector<std::vector<std::vector<Number>>> get_interior_basis_functions_as_vecs();

  /**
   * Get data that defines the sizes of interior EIM basis functions.
   */
  std::map<std::string,std::size_t>
    get_interior_basis_function_sizes(std::vector<unsigned int> & n_qp_per_elem);

  /**
   * Here we store an enum that defines the type of EIM error indicator
   * normalization that we use in get_eim_error_indicator(). The enum
   * is public so that it can be set in user code.
   *
   * RESIDUAL_SUM: Use the sum of the terms in the EIM residual to determine
   * the error indicator normalization. This ensures that the error indicator
   * value will be at most 1.0, which may be a desirable property of the
   * indicator.
   *
   * RHS: Use only the right-hand side value for the EIM residual to
   * determine the error indicator normalization.
   *
   * MAX_RHS: Use the maximum value in the EIM RHS vector to determine
   * the error indicator normalization (default). This is helpful when
   * the values at some EIM points are much larger than others, since in
   * this scenario we typically want to normalize the error indicator
   * based on the largest values in order to avoid overestimating the
   * error.
   */
  enum EimErrorIndicatorNormalization { RESIDUAL_SUM, RESIDUAL_RHS, MAX_RHS };

  EimErrorIndicatorNormalization eim_error_indicator_normalization;

  /**
   * If this boolean is true then we clamp EIM error indicator values to
   * be at most 1. This is often desirable since we typically think of
   * the error indicator as a percentage, and a value of 1 means 100%
   * error. We set this to true by default.
   */
  bool limit_eim_error_indicator_to_one;

protected:

  /**
   * This vector specifies which EIM variables we want to write to disk and/or
   * project to nodes for plotting purposes. By default this is an empty
   * set, but can be updated in subclasses to specify the EIM variables that
   * are relevant for visualization.
   *
   * We identify groups of variables with one or more variables in a group.
   * The purpose of using a group is often we plot multiple components of
   * a tensor-valued or vector-valued quantity, so it makes sense to refer
   * to the entire group of variables together in those cases.
   */
  std::vector<EIMVarGroupPlottingInfo> _eim_vars_to_project_and_write;

  /**
   * This set that specifies which EIM variables will be scaled during EIM
   * enrichment so that their maximum value matches the maximum value across
   * all variables. This is helpful in cases where some components are much
   * smaller in magnitude than others, since in those cases if we do not apply
   * component scaling to the small components then the accuracy of the EIM
   * approximation for those components will not be controlled well by the
   * EIM enrichment process.
   */
  std::set<unsigned int> _scale_components_in_enrichment;

private:

  /**
   * Method that writes out element interior EIM basis functions. This may be called by
   * write_out_basis_functions().
   */
  void write_out_interior_basis_functions(const std::string & directory_name,
                                          bool write_binary_basis_functions);

  /**
   * Method that writes out element side EIM basis functions. This may be called by
   * write_out_basis_functions().
   */
  void write_out_side_basis_functions(const std::string & directory_name,
                                      bool write_binary_basis_functions);

  /**
   * Method that writes out element node EIM basis functions. This may be called by
   * write_out_basis_functions().
   */
  void write_out_node_basis_functions(const std::string & directory_name,
                                      bool write_binary_basis_functions);

  /**
   * Method that reads in element interior EIM basis functions. This may be called by
   * read_in_basis_functions().
   */
  void read_in_interior_basis_functions(const System & sys,
                                        const std::string & directory_name,
                                        bool read_binary_basis_functions);

  /**
   * Method that reads in element side EIM basis functions. This may be called by
   * read_in_basis_functions().
   */
  void read_in_side_basis_functions(const System & sys,
                                    const std::string & directory_name,
                                    bool read_binary_basis_functions);

  /**
   * Method that reads in element node EIM basis functions. This may be called by
   * read_in_basis_functions().
   */
  void read_in_node_basis_functions(const System & sys,
                                    const std::string & directory_name,
                                    bool read_binary_basis_functions);

  /**
   * Helper function called by write_out_interior_basis_functions() to
   * get basis function \p bf_index stored as a std::vector per variable.
   */
  std::vector<std::vector<Number>> get_interior_basis_function_as_vec_helper(
    unsigned int n_vars,
    unsigned int n_qp_data,
    unsigned int bf_index);

  /**
   * The EIM solution coefficients from the most recent call to rb_eim_solves().
   */
  std::vector<DenseVector<Number>> _rb_eim_solutions;

  /**
   * If we're using the EIM error indicator, then we store the error indicator
   * values corresponding to _rb_eim_solutions here.
   */
  std::vector<std::pair<Real,Real>> _rb_eim_error_indicators;

  /**
   * Storage for EIM solutions from the training set. This is typically used in
   * the case that we have is_lookup_table==true in our RBParametrizedFunction,
   * since in that case we need to store all the EIM solutions on the training
   * set so that we do not always need to refer to the lookup table itself
   * (since in some cases, like in the Online stage, the lookup table is not
   * available).
   */
  std::vector<DenseVector<Number>> _eim_solutions_for_training_set;

  /**
   * The parameters and the number of basis functions that were used in the
   * most recent call to rb_eim_solves(). We store this so that we can
   * check if we can skip calling rb_eim_solves() again if the inputs
   * haven't changed.
   */
  std::vector<RBParameters> _rb_eim_solves_mus;
  unsigned int _rb_eim_solves_N;

  /**
   * Dense matrix that stores the lower triangular
   * interpolation matrix that can be used
   */
  DenseMatrix<Number> _interpolation_matrix;

  /**
   * We store the EIM interpolation point data in this object.
   */
  VectorizedEvalInput _vec_eval_input;

  /**
   * In the case of a "vector-valued" EIM, this vector determines which
   * component of the parameterized function we sample at each EIM point.
   */
  std::vector<unsigned int> _interpolation_points_comp;

  /**
   * Here we store the spatial indices that were initialized by
   * initialize_spatial_indices_at_interp_pts(). These are relevant
   * in the case that _parametrized_function is defined by indexing
   * into separate data based on the mesh-based data.
   */
  std::vector<std::vector<unsigned int>> _interpolation_points_spatial_indices;

  /**
   * Store the parametrized function that will be approximated
   * by this EIM system. Note that the parametrized function
   * may have more than one component, and each component is
   * approximated by a separate variable in the EIM system.
   */
  std::unique_ptr<RBParametrizedFunction> _parametrized_function;

  /**
   * The vector of RBTheta objects that are created to point to
   * this RBEIMEvaluation.
   */
  std::vector<std::unique_ptr<RBTheta>> _rb_eim_theta_objects;

  /**
   * The EIM basis functions. We store values at quadrature points
   * on elements that are local to this processor. The indexing
   * is as follows:
   *   basis function index --> element ID --> variable --> quadrature point --> value
   * We use a map to index the element ID, since the IDs on this processor in
   * general will not start at zero.
   */
  std::vector<QpDataMap> _local_eim_basis_functions;

  /**
   * The EIM basis functions on element sides. We store values at quadrature points
   * on elements that are local to this processor. The indexing
   * is as follows:
   *   basis function index --> (element ID,side index) --> variable --> quadrature point --> value
   * We use a map to index the element ID, since the IDs on this processor in
   * general will not start at zero.
   */
  std::vector<SideQpDataMap> _local_side_eim_basis_functions;

  /**
   * The EIM basis functions on element nodes (e.g. on a nodeset). We store values at nodes
   * that are local to this processor. The indexing
   * is as follows:
   *   basis function index --> node ID --> variable --> value
   * We use a map to index the node ID, since the IDs on this processor in
   * general will not start at zero.
   */
  std::vector<NodeDataMap> _local_node_eim_basis_functions;

  /**
   * Print the contents of _local_eim_basis_functions to libMesh::out.
   * Helper function mainly useful for debugging.
   */
  void print_local_eim_basis_functions() const;

  /**
   * Helper function that gathers the contents of
   * _local_eim_basis_functions to processor 0 in preparation for
   * printing to file.
   */
  void gather_bfs();

  /**
   * Same as gather_bfs() except for side data.
   */
  void side_gather_bfs();

  /**
   * Same as gather_bfs() except for node data.
   */
  void node_gather_bfs();

  /**
   * Helper function that distributes the entries of
   * _local_eim_basis_functions to their respective processors after
   * they are read in on processor 0.
   */
  void distribute_bfs(const System & sys);

  /**
   * Same as distribute_bfs() except for side data.
   */
  void side_distribute_bfs(const System & sys);

  /**
   * Same as distribute_bfs() except for node data.
   */
  void node_distribute_bfs(const System & sys);

  /**
   * Boolean to indicate if we skip updating _rb_eim_solutions in rb_eim_solves().
   * This is relevant for cases when we set up _rb_eim_solutions elsewhere and we
   * want to avoid changing it.
   */
  bool _preserve_rb_eim_solutions;

  /**
   * Indicate if the EIM error indicator is active in RB EIM solves. Note that
   * this is distinct from use_eim_error_indicator(), since use_eim_error_indicator()
   * indicates if this RBEIMEvaluation has an EIM error indicator defined,
   * whereas _is_eim_error_indicator_active is used to turn on or off the
   * error indicator. This primary purpose of _is_eim_error_indicator_active
   * is to turn the error indicator off during EIM training (when it is not relevant)
   * and to turn it on during "online solves".
   */
  bool _is_eim_error_indicator_active;

  /**
   * Here we store an extra row of the interpolation matrix which is used to
   * compute the EIM error indicator. This stores the EIM basis function
   * values at the extra point associated with the error indicator.
   */
  DenseVector<Number> _error_indicator_interpolation_row;

};

}

#endif // LIBMESH_RB_EIM_EVALUATION_H
