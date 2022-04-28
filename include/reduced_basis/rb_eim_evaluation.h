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
class Elem;

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
   * This is used to initalize the spatial indices data in _parametrized_function
   * so that the _parametrized_function can be used in the Online stage without
   * reconstructing the spatial indices on every element or node in the mesh.
   */
  void initialize_param_fn_spatial_indices();

  /**
   * Return the current number of EIM basis functions.
   */
  unsigned int get_n_basis_functions() const;

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
   */
  static Number get_parametrized_function_node_value(
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
   */
  Number get_eim_basis_function_node_value(unsigned int basis_function_index,
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
  void add_interpolation_points_phi_i_qp(const std::vector<Real> & phi_i_qp);
  void add_interpolation_points_spatial_indices(const std::vector<unsigned int> & spatial_indices);

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
  const std::vector<Real> & get_interpolation_points_phi_i_qp(unsigned int index) const;
  const std::vector<unsigned int> & get_interpolation_points_spatial_indices(unsigned int index) const;

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
  void add_basis_function_and_interpolation_data(
    const QpDataMap & bf,
    Point p,
    unsigned int comp,
    dof_id_type elem_id,
    subdomain_id_type subdomain_id,
    unsigned int qp,
    const std::vector<Point> & perturbs,
    const std::vector<Real> & phi_i_qp);

  /**
   * Add \p side_bf to our EIM basis.
   */
  void add_side_basis_function_and_interpolation_data(
    const SideQpDataMap & side_bf,
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
  void add_node_basis_function_and_interpolation_data(
    const NodeDataMap & node_bf,
    Point p,
    unsigned int comp,
    dof_id_type node_id,
    boundary_id_type boundary_id);

  /**
   * Set the observation points and components.
   */
  void set_observation_points(const std::vector<Point> & observation_points_xyz);

  /**
   * Get the number of observation points.
   */
  unsigned int get_n_observation_points() const;

  /**
   * Get the observation points.
   */
  const std::vector<Point> & get_observation_points() const;

  /**
   * Get the observation value for the specified basis function and observation point.
   */
  const std::vector<Number> & get_observation_values(unsigned int bf_index, unsigned int obs_pt_index) const;

  /**
   * Get a const reference to all the observation values, indexed as follows:
   *  basis_function index --> observation point index --> value.
   */
  const std::vector<std::vector<std::vector<Number>>> & get_observation_values() const;

  /**
   * Add values at the observation points for a new basis function.
   */
  void add_observation_values_for_basis_function(const std::vector<std::vector<Number>> & values);

  /**
   * Set all observation values.
   */
  void set_observation_values(const std::vector<std::vector<std::vector<Number>>> & values);

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
   * Project variable \p var of \p bf_data into the solution vector of System.
   */
  void project_qp_data_map_onto_system(System & sys,
                                       const QpDataMap & bf_data,
                                       unsigned int var);

  /**
   * Return a set that specifies which EIM variables will be projected
   * and written out in write_out_projected_basis_functions().
   * By default this returns an empty vector, but can be overridden in
   * subclasses to specify the EIM variables that are relevant for visualization.
   */
  virtual std::set<unsigned int> get_eim_vars_to_project_and_write() const;

  /**
   * Project all basis functions using project_qp_data_map_onto_system() and
   * then write out the resulting vectors.
   */
  void write_out_projected_basis_functions(System & sys,
                                           const std::string & directory_name = "offline_data");

  /**
   * Indicate whether we should apply scaling to the components of the parametrized
   * function during basis function enrichment in order give an approximately uniform
   * magnitude for all components. This is helpful in cases where the components vary
   * widely in magnitude.
   */
  virtual bool scale_components_in_enrichment() const;

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
   * The EIM solution coefficients from the most recent call to rb_eim_solves().
   */
  std::vector<DenseVector<Number>> _rb_eim_solutions;

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
   * We need to store interpolation point data in order to
   * evaluate parametrized functions at the interpolation points.
   * This requires the xyz locations, the components to evaluate,
   * and the subdomain IDs.
   */
  std::vector<Point> _interpolation_points_xyz;
  std::vector<unsigned int> _interpolation_points_comp;
  std::vector<subdomain_id_type> _interpolation_points_subdomain_id;

  /**
   * We also store perturbations of the xyz locations that may be
   * needed to evaluate finite difference approximations to derivatives.
   */
  std::vector<std::vector<Point>> _interpolation_points_xyz_perturbations;

  /**
   * We also store the element ID and qp index of each interpolation
   * point so that we can evaluate our basis functions at these
   * points by simply looking up the appropriate stored values.
   * This data is only needed during the EIM training.
   */
  std::vector<dof_id_type> _interpolation_points_elem_id;
  std::vector<unsigned int> _interpolation_points_qp;

  /**
   * If the EIM approximation applies to element sides, then we need to
   * store the side index and boundary ID for each quadrature point.
   */
  std::vector<unsigned int> _interpolation_points_side_index;
  std::vector<boundary_id_type> _interpolation_points_boundary_id;

  /**
   * If the EIM approximation applies to element nodes (e.g. from a nodeset),
   * then we need to store the ID for each node. We also store the associated
   * boundary ID in this case, using _interpolation_points_boundary_id.
   */
  std::vector<dof_id_type> _interpolation_points_node_id;

  /**
   * We store the shape function values at the qp as well. These values
   * allows us to evaluate parametrized functions that depend on nodal
   * data.
   */
  std::vector<std::vector<Real>> _interpolation_points_phi_i_qp;

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
   * Let {p_1,...,p_n} be a set of n "observation points", where we can
   * observe the values of our EIM basis functions. Also, let
   * {comp_k} be the components of the EIM basis function that
   * we will observe. Then the corresponding observation values, v_ijk,
   * are given by:
   *  v_ijk = eim_basis_function[i][p_j][comp_k].
   *
   * These observation values can be used to observe the EIM approximation
   * at specific points of interest, where the points of interest are defined
   * by the observation points.
   *
   * _observation_points_value is indexed as follows:
   *  basis_function index --> observation point index --> comp index --> value
   */
  std::vector<Point> _observation_points_xyz;
  std::vector<std::vector<std::vector<Number>>> _observation_points_values;

  /**
   * Boolean to indicate if we skip updating _rb_eim_solutions in rb_eim_solves().
   * This is relevant for cases when we set up _rb_eim_solutions elsewhere and we
   * want to avoid changing it.
   */
  bool _preserve_rb_eim_solutions;

};

}

#endif // LIBMESH_RB_EIM_EVALUATION_H
