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

#ifndef LIBMESH_RB_PARAMETRIZED_FUNCTION_H
#define LIBMESH_RB_PARAMETRIZED_FUNCTION_H

// libMesh includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <unordered_map>
#include <vector>
#include <map>
#include <set>

namespace libMesh
{

class RBParameters;
class Point;
class System;

/**
 * A simple functor class that provides a RBParameter-dependent function.
 *
 * \author David Knezevic
 * \date 2012
 * \brief Provides a reduced basis parameterized function.
 */
class RBParametrizedFunction
{
public:

  /**
   * Constructor.
   */
  RBParametrizedFunction();

  /**
   * Special functions.
   * - This class can be default copy/move assigned/constructed.
   * - The destructor is defaulted out-of-line.
   */
  RBParametrizedFunction (RBParametrizedFunction &&) = default;
  RBParametrizedFunction (const RBParametrizedFunction &) = default;
  RBParametrizedFunction & operator= (const RBParametrizedFunction &) = default;
  RBParametrizedFunction & operator= (RBParametrizedFunction &&) = default;
  virtual ~RBParametrizedFunction();

  /**
   * Specify the number of components in this parametrized function.
   * A scalar-valued function has one component, a vector-valued
   * function has more than one component.
   */
  virtual unsigned int get_n_components() const = 0;

  /**
   * Evaluate the parametrized function at the specified point for
   * parameter \p mu.  If requires_xyz_perturbations==false, then
   * xyz_perturb will not be used.
   *
   * In this case we return the value for component \p comp only, but
   * the base class implementation simply calls the vector-returning
   * evaluate() function below and returns the comp'th component, so
   * derived classes should provide a more efficient routine or just call
   * the vector-returning function instead.
   */
  virtual Number evaluate_comp(const RBParameters & mu,
                               unsigned int comp,
                               const Point & xyz,
                               dof_id_type elem_id,
                               unsigned int qp,
                               subdomain_id_type subdomain_id,
                               const std::vector<Point> & xyz_perturb,
                               const std::vector<Real> & phi_i_qp);

  /**
   * Same as evaluate_comp() but for element sides.
   */
  virtual Number side_evaluate_comp(const RBParameters & mu,
                                    unsigned int comp,
                                    const Point & xyz,
                                    dof_id_type elem_id,
                                    unsigned int side_index,
                                    unsigned int qp,
                                    subdomain_id_type subdomain_id,
                                    boundary_id_type boundary_id,
                                    const std::vector<Point> & xyz_perturb,
                                    const std::vector<Real> & phi_i_qp);

  /**
   * Same as evaluate_comp() but for element nodes.
   */
  virtual Number node_evaluate_comp(const RBParameters & mu,
                                    unsigned int comp,
                                    const Point & xyz,
                                    dof_id_type node_id,
                                    boundary_id_type boundary_id);

  /**
   * Evaluate the parametrized function at the specified point for
   * parameter \p mu.  If requires_xyz_perturbations==false, then
   * xyz_perturb will not be used.
   *
   * In this case we evaluate for all components.
   */
  virtual std::vector<Number> evaluate(const RBParameters & mu,
                                       const Point & xyz,
                                       dof_id_type elem_id,
                                       unsigned int qp,
                                       subdomain_id_type subdomain_id,
                                       const std::vector<Point> & xyz_perturb,
                                       const std::vector<Real> & phi_i_qp);

  /**
   * Same as evaluate() but for element sides.
   */
  virtual std::vector<Number> side_evaluate(const RBParameters & mu,
                                            const Point & xyz,
                                            dof_id_type elem_id,
                                            unsigned int side_index,
                                            unsigned int qp,
                                            subdomain_id_type subdomain_id,
                                            boundary_id_type boundary_id,
                                            const std::vector<Point> & xyz_perturb,
                                            const std::vector<Real> & phi_i_qp);

  /**
   * Same as evaluate() but for element nodes.
   */
  virtual std::vector<Number> node_evaluate(const RBParameters & mu,
                                            const Point & xyz,
                                            dof_id_type node_id,
                                            boundary_id_type boundary_id);

  /**
   * Vectorized version of evaluate. If requires_xyz_perturbations==false, then all_xyz_perturb will not be used.
   */
  virtual void vectorized_evaluate(const std::vector<RBParameters> & mus,
                                   const std::vector<Point> & all_xyz,
                                   const std::vector<dof_id_type> & elem_ids,
                                   const std::vector<unsigned int> & qps,
                                   const std::vector<subdomain_id_type> & sbd_ids,
                                   const std::vector<std::vector<Point>> & all_xyz_perturb,
                                   const std::vector<std::vector<Real>> & phi_i_qp,
                                   std::vector<std::vector<std::vector<Number>>> & output);

  /**
   * Same as vectorized_evaluate() but on element sides.
   */
  virtual void side_vectorized_evaluate(const std::vector<RBParameters> & mus,
                                        const std::vector<Point> & all_xyz,
                                        const std::vector<dof_id_type> & elem_ids,
                                        const std::vector<unsigned int> & side_indices,
                                        const std::vector<unsigned int> & qps,
                                        const std::vector<subdomain_id_type> & sbd_ids,
                                        const std::vector<boundary_id_type> & boundary_ids,
                                        const std::vector<std::vector<Point>> & all_xyz_perturb,
                                        const std::vector<std::vector<Real>> & phi_i_qp,
                                        std::vector<std::vector<std::vector<Number>>> & output);

  /**
   * Same as vectorized_evaluate() but on element nodes.
   */
  virtual void node_vectorized_evaluate(const std::vector<RBParameters> & mus,
                                        const std::vector<Point> & all_xyz,
                                        const std::vector<dof_id_type> & node_ids,
                                        const std::vector<boundary_id_type> & boundary_ids,
                                        std::vector<std::vector<std::vector<Number>>> & output);

  /**
   * Store the result of vectorized_evaluate. This is helpful during EIM training,
   * since we can pre-evaluate and store the parameterized function for each training
   * sample. If requires_xyz_perturbations==false, then all_xyz_perturb will not be used.
   */
  virtual void preevaluate_parametrized_function_on_mesh(const RBParameters & mu,
                                                         const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                         const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                         const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb,
                                                         const System & sys);

  /**
   * Same as preevaluate_parametrized_function_on_mesh() except for mesh sides.
   */
  virtual void preevaluate_parametrized_function_on_mesh_sides(const RBParameters & mu,
                                                               const std::map<std::pair<dof_id_type,unsigned int>, std::vector<Point>> & side_all_xyz,
                                                               const std::map<std::pair<dof_id_type,unsigned int>, subdomain_id_type> & sbd_ids,
                                                               const std::map<std::pair<dof_id_type,unsigned int>, boundary_id_type> & side_boundary_ids,
                                                               const std::map<std::pair<dof_id_type,unsigned int>, unsigned int> & side_types,
                                                               const std::map<std::pair<dof_id_type,unsigned int>, std::vector<std::vector<Point>> > & side_all_xyz_perturb,
                                                               const System & sys);

  /**
   * Same as preevaluate_parametrized_function_on_mesh() except for mesh nodes.
   */
  virtual void preevaluate_parametrized_function_on_mesh_nodes(const RBParameters & mu,
                                                               const std::unordered_map<dof_id_type, Point> & all_xyz,
                                                               const std::unordered_map<dof_id_type, boundary_id_type> & node_boundary_ids,
                                                               const System & sys);

  /**
   * Look up the preevaluate values of the parametrized function for
   * component \p comp, element \p elem_id, and quadrature point \p qp.
   */
  virtual Number lookup_preevaluated_value_on_mesh(unsigned int comp,
                                                   dof_id_type elem_id,
                                                   unsigned int qp) const;

  /**
   * Look up the preevaluated values of the parametrized function for
   * component \p comp, element \p elem_id, \p side_index, and quadrature point \p qp.
   */
  virtual Number lookup_preevaluated_side_value_on_mesh(unsigned int comp,
                                                        dof_id_type elem_id,
                                                        unsigned int side_index,
                                                        unsigned int qp) const;

  /**
   * Look up the preevaluate values of the parametrized function for
   * component \p comp, node \p node_id.
   */
  virtual Number lookup_preevaluated_node_value_on_mesh(unsigned int comp,
                                                        dof_id_type node_id) const;

  /**
   * If this parametrized function is defined based on a lookup table then
   * we can call this function to initialize the table. This is a no-op by
   * default, but it can be overridden in subclasses as needed.
   */
  virtual void initialize_lookup_table();

  /**
   * Get the value stored in _parameter_independent_data associated with
   * \p region_name and \p property_name.
   */
  Number get_parameter_independent_data(const std::string & property_name,
                                        subdomain_id_type sbd_id) const;

  /**
   * For RBParametrizedFunctions defined on element sides or nodes, we get/set the boundary
   * IDs that this parametrized function is defined on.
   */
  const std::set<boundary_id_type> & get_parametrized_function_boundary_ids() const;
  void set_parametrized_function_boundary_ids(const std::set<boundary_id_type> & boundary_ids, bool is_nodal_boundary);

  /**
   * @return true if this parametrized function is defined on mesh sides.
   */
  bool on_mesh_sides() const;

  /**
   * @return true if this parametrized function is defined on mesh nodes.
   */
  bool on_mesh_nodes() const;

  /**
   * Evaluate the parametrized function for the parameter \p mu at the set of
   * \p observation_points. We return a vector of values at each observation
   * point since we may want to evaluate more than one component of the
   * parametrized function. We also provide \p elem_ids and \p sbd_ids
   * since that info can be required for the evaluation in some cases.
   */
  virtual std::vector<std::vector<Number>> evaluate_at_observation_points(const RBParameters & mu,
                                                                          const std::vector<Point> & observation_points,
                                                                          const std::vector<dof_id_type> & elem_ids,
                                                                          const std::vector<subdomain_id_type> & sbd_ids,
                                                                          const System & sys);

  /**
   * In some cases a parametrized function is defined based on array data that
   * we index into based on the spatial data from the mesh (e.g. element, node,
   * or side indices). We refer to the indices that we use to index into this
   * array data as "spatial indices". This method sets \p spatial_indices based
   * on the provided mesh-based indices.
   *
   * Note that \p spatial_indices is defined as a doubly-nested vector so that we
   * can handle the case where the spatial function evaluation requires us to have
   * indices from all nodes of an element, since in that case we need a vector of
   * indices (one per node) for each point. Other cases, such as when we define
   * the paramtrized function based on the element index only, only require a
   * singly-nested vector which we handle as a special case of the doubly-nested
   * vector.
   *
   * This method is typically used in the Offline stage in order to generate
   * and store the relevant spatial indices.
   *
   * This method is a no-op by default, but it can be overridden in subclasses
   * to provide the relevant behavior.
   */
  virtual void get_spatial_indices(std::vector<std::vector<unsigned int>> & spatial_indices,
                                   const std::vector<dof_id_type> & elem_ids,
                                   const std::vector<unsigned int> & side_indices,
                                   const std::vector<dof_id_type> & node_ids,
                                   const std::vector<unsigned int> & qps,
                                   const std::vector<subdomain_id_type> & sbd_ids,
                                   const std::vector<boundary_id_type> & boundary_ids);

  /**
   * The Online stage counterpart of get_spatial_indices(). This method
   * is used to initialize the spatial index data in this object so that
   * we can evaluate it during an Online solve.
   */
  virtual void initialize_spatial_indices(const std::vector<std::vector<unsigned int>> & spatial_indices,
                                          const std::vector<dof_id_type> & elem_ids,
                                          const std::vector<unsigned int> & side_indices,
                                          const std::vector<dof_id_type> & node_ids,
                                          const std::vector<unsigned int> & qps,
                                          const std::vector<subdomain_id_type> & sbd_ids,
                                          const std::vector<boundary_id_type> & boundary_ids);

  /**
   * Storage for pre-evaluated values. The indexing is given by:
   *   parameter index --> point index --> component index --> value.
   */
  std::vector<std::vector<std::vector<Number>>> preevaluated_values;

  /**
   * Indexing into preevaluated_values for the case where the preevaluated values
   * were obtained from evaluations at elements/quadrature points on a mesh.
   * The indexing here is:
   *   elem_id --> qp --> point_index
   * Then preevaluated_values[0][point_index] provides the vector of component values at
   * that point.
   */
  std::unordered_map<dof_id_type, std::vector<unsigned int>> mesh_to_preevaluated_values_map;

  /**
   * Similar to the above except this map stores the data on element sides.
   * The indexing here is:
   *  (elem_id,side index) --> qp --> point_index
   */
  std::map<std::pair<dof_id_type,unsigned int>, std::vector<unsigned int>> mesh_to_preevaluated_side_values_map;

  /**
   * Indexing into preevaluated_values for the case where the preevaluated values
   * were obtained from evaluations at elements/quadrature points on a mesh.
   * The indexing here is:
   *   node_id --> point_index
   */
  std::unordered_map<dof_id_type, unsigned int> mesh_to_preevaluated_node_values_map;

  /**
   * Boolean to indicate whether this parametrized function requires xyz perturbations
   * in order to evaluate function values. An example of where perturbations are
   * required is when the parametrized function is based on finite difference
   * approximations to derivatives.
   */
  bool requires_xyz_perturbations;

  /**
   * Boolean to indicate if this parametrized function is defined based on a lookup
   * table or not. If it is defined based on a lookup table, then the evaluation
   * functions will access a discrete parameter to determine the index to lookup.
   */
  bool is_lookup_table;

  /**
   * If this is a lookup table, then lookup_table_param_name specifies the parameter
   * that is used to index into the lookup table.
   */
  std::string lookup_table_param_name;

  /**
   * The finite difference step size in the case that this function in the case
   * that this function uses finite differencing.
   */
  Real fd_delta;

protected:

  /**
   * In some cases we need to store parameter-independent data which is related
   * to this function but since it is parameter-indepedent should not be returned
   * as part of evaluate().
   *
   * We index this data by "property name" --> subdomain_id --> value.
   */
  std::map<std::string, std::map<subdomain_id_type, Number>> _parameter_independent_data;

  /**
   * In the case of an RBParametrizedFunction defined on element sides, this defines
   * the set of boundary IDs that the function is defined on.
   */
  std::set<boundary_id_type> _parametrized_function_boundary_ids;

  /**
   * In the case that _parametrized_function_boundary_ids is not empty, then this
   * parametrized function is defined on a mesh boundary. This boolean inidicates
   * if the mesh boundary under consideration is a set of sides, or a set of nodes.
   */
  bool _is_nodal_boundary;

};

}

#endif // LIBMESH_RB_PARAMETRIZED_FUNCTION_H
