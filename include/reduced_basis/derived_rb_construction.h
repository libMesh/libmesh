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

#ifndef LIBMESH_DERIVED_RB_CONSTRUCTION_H
#define LIBMESH_DERIVED_RB_CONSTRUCTION_H

// Local includes
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fem_context.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/rb_construction.h"

// C++ includes

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * DerivedRBConstruction implements the Construction
 * stage for the "two-level" Reduced Basis method.
 * In this context we obtain our "derived" basis functions
 * from a reference reduced basis approximation.
 *
 * @author David J. Knezevic, 2009
 */

template<class Base>
class DerivedRBConstruction : public Base
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DerivedRBConstruction (EquationSystems& es,
                         const std::string& name,
                         const unsigned int number);

  /**
   * The type of system.
   */
  typedef DerivedRBConstruction<Base> sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * @returns a string indicating the type of the system.
   */
  virtual std::string system_type () const;

  /**
   * Overload truth_solve so that it computes the
   * associated Unter RB solution.
   */
  virtual Real truth_solve(int plot_solution);

  /**
   * Overload train_reduced_basis to generate residual terms wrt truth
   * once the Greedy finishes.
   */
  virtual Real train_reduced_basis(const std::string& directory_name = "offline_data",
                                   const bool resize_rb_eval_data=true);

  /**
   * Set the uber_system's current_parameters to
   * match unter_system's current_parameters. We
   * require a virtual function here in case the
   * unter system has fewer parameters than the uber
   * system, for example.
   */
  virtual void set_uber_current_parameters();

  /**
   * Build a new DerivedRBEvaluation object.
   */
  virtual UniquePtr<RBEvaluation> build_rb_evaluation(const Parallel::Communicator &comm
                                                      LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Load the RB solution from the most recent solve
   * into the libMesh solution vector.
   */
  virtual void load_rb_solution();

  /**
   * Load the i^th derived basis function into vec.
   */
  virtual void load_basis_function(unsigned int i);

  /**
   * This function recomputes all the residual terms in order
   * to allow evaluation of the residual wrt the truth space.
   * Useful if we have performed the greedy algorithm with
   * respect to the uber space and need to replace with truth
   * residual terms.
   */
  void generate_residual_terms_wrt_truth();


  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The name of the uber RB system, i.e. the RB system
   * we use to develop the derived RB system.
   */
  std::string uber_system_name;

protected:

  /**
   * Get a copy of the specific derived basis function.
   */
  DenseVector<Number> get_derived_basis_function(unsigned int i);

  /**
   * Add a new basis function to the RB space. Overload
   * for uber-unter style enrichment, i.e. this function
   * is independent of the number of degrees of freedom
   * in the truth finite element discretization.
   */
  virtual void enrich_RB_space();

  /**
   * Compute the reduced basis matrices for the current basis.
   * This operation is based on the uber system and hence is
   * independent of the number of degrees of freedom
   * in the truth finite element discretization.
   */
  virtual void update_RB_system_matrices();

  /**
   * Compute the RHS terms that are combined `online'
   * to determine the dual norm of the residual. Overloaded
   * here for the two-stage RB method.
   */
  virtual void compute_Fq_representor_innerprods(bool compute_inner_products=true);

  /**
   * Compute the terms that are combined `online'
   * to determine the dual norm of the residual. Overload
   * here so that we perform an update based on the uber
   * system.
   */
  virtual void update_residual_terms(bool compute_inner_products=true);

private:

};

// And introduce convenient typedefs
typedef DerivedRBConstruction<RBConstruction> SteadyDerivedRBConstruction;

} // namespace libMesh


#endif // LIBMESH_DERIVED_RB_CONSTRUCTION_H
