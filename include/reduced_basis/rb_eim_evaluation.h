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
#include "libmesh/auto_ptr.h"
#include "libmesh/point.h"
#include "libmesh/rb_evaluation.h"
#include "libmesh/elem.h"
#include "libmesh/serial_mesh.h"

// C++ includes

namespace libMesh
{

class RBParameters;
class RBParametrizedFunction;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBEIMEvaluation extends RBEvaluation to
 * encapsulate the code and data required
 * to perform "online" evaluations for
 * EIM approximations.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBEIMEvaluation class definition

class RBEIMEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor.
   */
  RBEIMEvaluation (const libMesh::Parallel::Communicator &comm 
		   LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~RBEIMEvaluation ();

  /**
   * The type of the parent.
   */
  typedef RBEvaluation Parent;

  /**
   * Clear this object.
   */
  virtual void clear();

  /**
   * Resize the data structures for storing data associated
   * with this object.
   */
  virtual void resize_data_structures(const unsigned int Nmax,
                                      bool resize_error_bound_data=true);

  /**
   * Attach the parametrized function that we will approximate
   * using the Empirical Interpolation Method.
   */
  void attach_parametrized_function(RBParametrizedFunction* pf);


  /**
   * Get the number of parametrized functions that have
   * been attached to this system.
   */
  unsigned int get_n_parametrized_functions() const;

  /**
   * @return the value of the parametrized function that is being
   * approximated at the point \p p.
   * \p var_index specifies the
   * variable (i.e. the parametrized function index) to be evaluated.
   * \p elem specifies the element of the mesh that contains p.
   */
  Number evaluate_parametrized_function(unsigned int var_index,
                                        const Point& p,
                                        const Elem& elem);

  /**
   * Calculate the EIM approximation to parametrized_function
   * using the first \p N EIM basis functions. Store the
   * solution coefficients in the member RB_solution.
   * @return the EIM a posteriori error bound.
   */
  virtual Real rb_solve(unsigned int N);

  /**
   * Calculate the EIM approximation for the given
   * right-hand side vector \p EIM_rhs. Store the
   * solution coefficients in the member RB_solution.
   */
  void rb_solve(DenseVector<Number>& EIM_rhs);

  /**
   * Build a vector of RBTheta objects that accesses the components
   * of the RB_solution member variable of this RBEvaluation.
   * Store these objects in the member vector rb_theta_objects.
   */
  void initialize_eim_theta_objects();

  /**
   * @return the vector of theta objects that point to this RBEIMEvaluation.
   */
  std::vector<RBTheta*> get_eim_theta_objects();

  /**
   * Build a theta object corresponding to EIM index \p index.
   * The default implementation builds an RBEIMTheta object, possibly
   * override in subclasses if we need more specialized behavior.
   */
  virtual AutoPtr<RBTheta> build_eim_theta(unsigned int index);

  /**
   * Write out all the data to text files in order to segregate the
   * Offline stage from the Online stage.
   */
  virtual void write_offline_data_to_files(const std::string& directory_name = "offline_data",
                                           const bool write_binary_data=true);

  /**
   * Read in the saved Offline reduced basis data
   * to initialize the system for Online solves.
   */
  virtual void read_offline_data_from_files(const std::string& directory_name = "offline_data",
                                            bool read_error_bound_data=true,
                                            const bool read_binary_data=true);

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * Dense matrix that stores the lower triangular
   * interpolation matrix that can be used
   */
  DenseMatrix<Number> interpolation_matrix;

  /**
   * The list of interpolation points, i.e. locations at
   * which the basis functions are maximized.
   */
  std::vector<Point> interpolation_points;

  /**
   * The corresponding list of variables indices at which
   * the interpolation points were identified.
   */
  std::vector<unsigned int> interpolation_points_var;
  
  /**
   * The corresponding list of elements at which
   * the interpolation points were identified.
   */
  std::vector<Elem*> interpolation_points_elem;

  /**
   * We also need an extra interpolation point and associated
   * variable and Elem for the "extra" solve we do at the end of
   * the Greedy algorithm.
   */
  Point extra_interpolation_point;
  unsigned int extra_interpolation_point_var;
  Elem* extra_interpolation_point_elem;

  /**
   * We also need a DenseVector to represent the corresponding
   * "extra" row of the interpolation matrix.
   */
  DenseVector<Number> extra_interpolation_matrix_row;

private:

  /**
   * Write out interpolation_points_elem by putting the elements into
   * a mesh and writing out the mesh.
   */
  void write_out_interpolation_points_elem(const std::string directory_name);

  /**
   * Read int interpolation_points_elem from a mesh.
   */
  void read_in_interpolation_points_elem(const std::string directory_name);

  /**
   * This vector stores the parametrized functions
   * that will be approximated in this EIM system.
   */
  std::vector<RBParametrizedFunction*> _parametrized_functions;

  /**
   * The vector of RBTheta objects that are created to point to
   * this RBEIMEvaluation.
   */
  std::vector<RBTheta*> _rb_eim_theta_objects;

  /**
   * We initialize RBEIMEvaluation so that it has an "empty" RBThetaExpansion, because
   * this isn't used at all in the EIM.
   */
  RBThetaExpansion _empty_rb_theta_expansion;

  /**
   * Store the parameters at which the previous solve was performed (so we can avoid
   * an unnecessary repeat solve).
   */
  RBParameters _previous_parameters;

  /**
   * Store the number of basis functions used for the previous solve (so we can avoid
   * an unnecessary repeat solve).
   */
  unsigned int _previous_N;

  /**
   * Store the previous error bound returned by rb_solve (so we can return it if we
   * are avoiding an unnecessary repeat solve).
   */
  Real _previous_error_bound;

  /**
   * Mesh object that we use to store copies of the elements associated with
   * interpolation points.
   */
  SerialMesh _interpolation_points_mesh;

};

}

#endif // LIBMESH_RB_EIM_EVALUATION_H
