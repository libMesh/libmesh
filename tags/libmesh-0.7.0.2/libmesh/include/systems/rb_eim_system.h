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

#ifndef __rb_eim_system_h__
#define __rb_eim_system_h__

#include "rb_system.h"
#include "mesh_function.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBEIMSystem implements the Empirical Interpolation Method (EIM)
 * that can be used to generate an affine approximation to non-affine
 * operators.
 *
 * @author David J. Knezevic, 2010
 */

/**
 * The typedef for the function to be approximated.
 */
typedef Number (*parametrized_function_fptr)(const Point&, const RBSystem&);

// ------------------------------------------------------------
// RBEIMSystem class definition

class RBEIMSystem : public RBSystem
{
public:

  enum BEST_FIT_TYPE { PROJECTION_BEST_FIT, EIM_BEST_FIT };

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBEIMSystem (EquationSystems& es,
               const std::string& name,
               const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBEIMSystem ();

  /**
   * The type of system.
   */
  typedef RBEIMSystem sys_type;
  
  /**
   * The type of the parent.
   */
  typedef RBSystem Parent;
  
  /**
   * @returns a string indicating the type of the system.
   */
  virtual std::string system_type () const;
  
  /**
   * Initialize this system so that we can perform
   * RB calculations.
   */
  virtual void initialize_RB_system(bool online_mode);
 
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
   * Calculate the EIM approximation to parametrized_function
   * using the first \p N EIM basis functions. Store the
   * solution coefficients in the member RB_solution.
   * @return the EIM a posteriori error bound.
   */
  virtual Real RB_solve(unsigned int N);
  
  /**
   * Calculate the EIM approximation for the given
   * right-hand side vector \p EIM_rhs. Store the
   * solution coefficients in the member RB_solution.
   */
  void RB_solve(DenseVector<Number>& EIM_rhs);

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
   * Override attach_theta_q_a to just throw an error. Should
   * use attach_A_q in RBSystem and its subclasses.
   */
  void attach_paramerized_function(parametrized_function_fptr fptr)
    { parametrized_functions.push_back(fptr); }
  
  /**
   * Get the number of parametrized_functions that have
   * been attached to this system.
   */
  unsigned int get_n_parametrized_functions() const
    { return parametrized_functions.size(); }

  /**
   * @return the value of the parametrized function that is
   * being approximated.
   */
  Number evaluate_parametrized_function(unsigned int var, const Point& p);
  
  /**
   * @return the number of affine terms defined by the current EIM
   * approximation. Each function is typically used in an associated
   * reduced basis approximation.
   */
  virtual unsigned int get_n_affine_functions() const;
  
  /**
   * Evaluate the affine function at the specified points on element.
   */
  std::vector<Number> evaluate_current_affine_function(Elem& element,
                                                       const std::vector<Point>& qpoints);

  /**
   * Load a GHOSTED version of the basis function specified by function_index
   * into the vector current_ghosted_bf so that we can efficiently interpolate
   * this basis function.
   */
  void cache_ghosted_basis_function(unsigned int function_index);

  /**
   * Clear the basis functions and all basis-function-dependent data.
   * Overload in subclasses to clear any extra data.
   */
  virtual void clear_basis_function_dependent_data();
  
  /**
   * Write out all the data to text files in order to segregate the
   * Offline stage from the Online stage.
   */
  virtual void write_offline_data_to_files(const std::string& directory_name = "offline_data");
  
  /**
   * Read in the saved Offline reduced basis data
   * to initialize the system for Online solves.
   */
  virtual void read_offline_data_from_files(const std::string& directory_name = "offline_data");

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
   * We also need an extra interpolation point and associated
   * variable for the "extra" solve we do at the end of
   * the Greedy algorithm.
   */
  Point extra_interpolation_point;
  unsigned int extra_interpolation_point_var;
  
  /**
   * We also need a DenseVector to represent the corresponding
   * "extra" row of the interpolation matrix.
   */
  DenseVector<Number> extra_interpolation_matrix_row;
  
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
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

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
  virtual Real get_RB_error_bound() { return compute_best_fit_error(); }
  
  /**
   * Function that indicates when to terminate the Greedy
   * basis training. Overload in subclasses to specialize.
   */
  virtual bool greedy_termination_test(Real training_greedy_error, int count);
  
private:

  /**
   * A mesh function to interpolate on the mesh.
   */
  MeshFunction* mesh_function;
  
  /**
   * This flag allows us to perform one extra Greedy step
   * in order to compute the data needed for the EIM
   * a posteriori error bound in the case that we use
   * all of our basis functions.
   */
  bool performing_extra_greedy_step;
  
  /**
   * This vector stores the parametrized functions
   * that will be approximated in this EIM system.
   */
  std::vector<parametrized_function_fptr> parametrized_functions;
  
  /**
   * The current basis function that we sample to evaluate the
   * empirical interpolation approximation. This will be a GHOSTED
   * vector to facilitate interpolation in the case of multiple processors.
   */
  AutoPtr< NumericVector<Number> > current_ghosted_bf;
  
  /**
   * We also need to store a variable number associated wtih current_ghosted_bf
   * because we assume that each variable can correspond to a different
   * EIM approximation.
   */
  unsigned int current_variable_number;
  
  /**
   * We also need an extra vector in which we can store a serialized
   * copy of the solution vector so that we can use MeshFunction
   * in parallel.
   */
  AutoPtr< NumericVector<Number> > serialized_vector;
  
  /**
   * This flag indicates whether or not we evaluate the error
   * estimate in RB_solve. We need this to turn off error
   * estimation during the Greedy algorithm.
   */
  bool eval_error_estimate;

};

} // namespace libMesh

#endif
