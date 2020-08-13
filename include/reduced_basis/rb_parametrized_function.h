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

namespace libMesh
{

class RBParameters;
class Point;

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
                               subdomain_id_type subdomain_id,
                               const std::vector<Point> & xyz_perturb);

  /**
   * Evaluate the parametrized function at the specified point for
   * parameter \p mu.  If requires_xyz_perturbations==false, then
   * xyz_perturb will not be used.
   *
   * In this case we evaluate for all components.
   */
  virtual std::vector<Number> evaluate(const RBParameters & mu,
                                       const Point & xyz,
                                       subdomain_id_type subdomain_id,
                                       const std::vector<Point> & xyz_perturb) = 0;

  /**
   * Vectorized version of evaluate. If requires_xyz_perturbations==false, then all_xyz_perturb will not be used.
   */
  virtual void vectorized_evaluate(const std::vector<RBParameters> & mus,
                                   const std::vector<Point> & all_xyz,
                                   const std::vector<subdomain_id_type> & sbd_ids,
                                   const std::vector<std::vector<Point>> & all_xyz_perturb,
                                   std::vector<std::vector<std::vector<Number>>> & output);

  /**
   * Store the result of vectorized_evaluate. This is helpful during EIM training,
   * since we can pre-evaluate and store the parameterized function for each training
   * sample. If requires_xyz_perturbations==false, then all_xyz_perturb will not be used.
   */
  virtual void preevaluate_parametrized_function_on_mesh(const RBParameters & mu,
                                                         const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                         const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                         const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb);

  /**
   * Look up the preevaluate values of the parametrized function for
   * component \p comp, element \p elem_id, and quadrature point \p qp.
   */
  virtual Number lookup_preevaluated_value_on_mesh(unsigned int comp,
                                                   dof_id_type elem_id,
                                                   unsigned int qp) const;

  /**
   * If this parametrized function is defined based on a lookup table then
   * we can call this function to initialize the table. This is a no-op by
   * default, but it can be overridden in subclasses as needed.
   */
  virtual void initialize_lookup_table();

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
};

}

#endif // LIBMESH_RB_PARAMETRIZED_FUNCTION_H
