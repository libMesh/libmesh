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
#include <unordered_map>

namespace libMesh
{

class RBParameters;
class RBParametrizedFunction;
class Elem;
class RBTheta;

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
   * Destructor.
   */
  virtual ~RBEIMEvaluation ();

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
   * Get a const reference to the parametrized function.
   */
  RBParametrizedFunction & get_parametrized_function();

  /**
   * Calculate the EIM approximation to parametrized_function
   * using the first \p N EIM basis functions. Store the
   * solution coefficients in the member _eim_solution.
   * \returns The EIM a posteriori error bound.
   */
  virtual Real rb_eim_solve(unsigned int N);

  /**
   * Calculate the EIM approximation for the given
   * right-hand side vector \p EIM_rhs. Store the
   * solution coefficients in the member _eim_solution.
   */
  void rb_eim_solve(DenseVector<Number> & EIM_rhs);

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
  void decrement_vector(std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & v,
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
    const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & pf,
    dof_id_type elem_id,
    unsigned int comp,
    std::vector<Number> & values);

  /**
   * Same as above, except that we just return the value at the qp^th
   * quadrature point.
   */
  static Number get_parametrized_function_value(
    const Parallel::Communicator & comm,
    const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & pf,
    dof_id_type elem_id,
    unsigned int comp,
    unsigned int qp);

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
   * Same as above, except that we just return the value at the qp^th
   * quadrature point.
   */
  Number get_eim_basis_function_value(unsigned int basis_function_index,
                                      dof_id_type elem_id,
                                      unsigned int comp,
                                      unsigned int qp) const;

  /**
   * Get a reference to the i^th basis function.
   */
  const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> &
    get_basis_function(unsigned int i) const;

  /**
   * Return a const reference to the EIM solution coefficients from the most
   * recent solve.
   */
  const DenseVector<Number> & get_rb_eim_solution() const;

  /**
   * Set the data associated with EIM interpolation points.
   */
  void add_interpolation_points_xyz(Point p);
  void add_interpolation_points_comp(unsigned int comp);
  void add_interpolation_points_subdomain_id(subdomain_id_type sbd_id);
  void add_interpolation_points_xyz_perturbations(const std::vector<Point> & perturbs);
  void add_interpolation_points_elem_id(dof_id_type elem_id);
  void add_interpolation_points_qp(unsigned int qp);

  /**
   * Get the data associated with EIM interpolation points.
   */
  Point get_interpolation_points_xyz(unsigned int index) const;
  unsigned int get_interpolation_points_comp(unsigned int index) const;
  subdomain_id_type get_interpolation_points_subdomain_id(unsigned int index) const;
  const std::vector<Point> & get_interpolation_points_xyz_perturbations(unsigned int index) const;
  dof_id_type get_interpolation_points_elem_id(unsigned int index) const;
  unsigned int get_interpolation_points_qp(unsigned int index) const;

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
    const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & bf,
    Point p,
    unsigned int comp,
    dof_id_type elem_id,
    subdomain_id_type subdomain_id,
    unsigned int qp,
    const std::vector<Point> & perturbs);

  /**
   * Boolean to indicate whether we evaluate a posteriori error bounds
   * when eim_solve is called.
   */
  bool evaluate_eim_error_bound;

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

private:

  /**
   * The EIM solution coefficients from the most recent eim_solve().
   */
  DenseVector<Number> _rb_eim_solution;

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
   * generally will not start at zero.
   */
  std::vector<std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> > _local_eim_basis_functions;

};

}

#endif // LIBMESH_RB_EIM_EVALUATION_H
