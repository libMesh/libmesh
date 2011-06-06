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
 * Functor class in which we can define function to be approximated.
 */
class ParametrizedFunction
{
public:

  /**
   * Evaluate this parametrized function for the parameter value
   * \p mu at the point \p p.
   */
  virtual Number evaluate(std::vector<Real>& , const Point& ) { return 0.; }

};

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
   * the Construction stage of the RB method.
   */
  virtual void initialize_RB_system(RBEvaluation& rb_evaluation_in);

  /**
   * Read parameters in from file and set up this system
   * accordingly.
   */
  virtual void process_parameters_file (const std::string& parameters_filename);
 
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
   * Build a new RBEIMEvaluation object.
   */
  virtual AutoPtr<RBEvaluation> build_rb_evaluation();

  /**
   * Override attach_theta_q_a to just throw an error. Should
   * use attach_A_q in RBSystem and its subclasses.
   */
  void attach_paramerized_function(ParametrizedFunction* pf)
    { parametrized_functions.push_back(pf); }
  
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
  Number evaluate_parametrized_function(unsigned int index, const Point& p);
  
  /**
   * @return the number of affine terms defined by the current EIM
   * approximation. Each function is typically used in an associated
   * reduced basis approximation.
   */
  virtual unsigned int get_n_affine_functions() const;
  
  /**
   * Evaluate the basis function \p index at the points \p qpoints
   * on element \p element.
   * @return a vector of values corresponding to qpoints.
   */
  std::vector<Number> evaluate_basis_function(unsigned int index,
                                              Elem& element,
                                              const std::vector<Point>& qpoints);

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

  /**
   * Helper function to load a GHOSTED version of the basis function specified
   * by \p basis_function_index_in into the vector current_ghosted_bf so that
   * we can efficiently interpolate this basis function.
   * If basis_function_index_in == current_bf_index, then this function does nothing.
   */
  void set_current_basis_function(unsigned int basis_function_index_in);
  
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
  std::vector<ParametrizedFunction*> parametrized_functions;
  
  /**
   * The current basis function that we sample to evaluate the
   * empirical interpolation approximation. This will be a GHOSTED
   * vector to facilitate interpolation in the case of multiple processors.
   */
  AutoPtr< NumericVector<Number> > current_ghosted_bf;
  
  /**
   * We also need to store a basis function index to identify which basis function
   * is currently stored in current_ghosted_bf. This allows us to cache the basis
   * function and avoid unnecessarily reloading it all the time.
   */
  unsigned int current_bf_index;
  
  /**
   * We also need an extra vector in which we can store a serialized
   * copy of the solution vector so that we can use MeshFunction
   * in parallel.
   */
  AutoPtr< NumericVector<Number> > serialized_vector;

};

} // namespace libMesh

#endif
