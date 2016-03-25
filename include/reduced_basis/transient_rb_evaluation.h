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

#ifndef LIBMESH_TRANSIENT_RB_EVALUATION_H
#define LIBMESH_TRANSIENT_RB_EVALUATION_H

// rbOOmit includes
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_temporal_discretization.h"

// libMesh includes

// C++ includes

namespace libMesh
{

class TransientRBThetaExpansion;

/**
 * This class is part of the rbOOmit framework.
 *
 * TransientRBEvaluation extends RBEvaluation to
 * encapsulate the code and data required
 * to perform "online" RB evaluations for
 * Linear Time Invariant (LTI) transient problems.
 *
 * We can handle time controls on the RHS as h(t)*f(x,\f$ \mu \f$).
 * See Martin Grepl's thesis for more details.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class TransientRBEvaluation : public RBEvaluation, public RBTemporalDiscretization
{
public:

  /**
   * Constructor.
   */
  TransientRBEvaluation (const Parallel::Communicator & comm_in
                         LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  ~TransientRBEvaluation ();

  /**
   * The type of the parent.
   */
  typedef RBEvaluation Parent;

  /**
   * Clear this TransientRBEvaluation object.
   * Override to also clear the M_q representors
   */
  virtual void clear() libmesh_override;

  /**
   * Resize and clear the data vectors corresponding to the
   * value of \p Nmax. Optionally resize the data structures
   * required for the error bound.
   * Overridden to resize data relevant in the time-dependent
   * case.
   */
  virtual void resize_data_structures(const unsigned int Nmax,
                                      bool resize_error_bound_data=true) libmesh_override;

  /**
   * Perform online solve for current_params
   * with the N basis functions. Overridden
   * to perform a time-dependent solve.
   */
  virtual Real rb_solve(unsigned int N) libmesh_override;

  /**
   * If a solve has already been performed, then we cached some data
   * and we can perform a new solve much more rapidly
   * (with the same parameters but a possibly different initial condition/rhs control).
   */
  virtual Real rb_solve_again();

  /**
   * @return a scaling factor that we can use to provide a consistent
   * scaling of the RB error bound across different parameter values.
   */
  virtual Real get_error_bound_normalization() libmesh_override;

  /**
   * Specifies the residual scaling on the numerator to
   * be used in the a posteriori error bound. Override
   * in subclass in order to obtain the desired error bound.
   */
  virtual Real residual_scaling_numer(Real alpha_LB);

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution. This function uses the cached time-independent
   * data.
   */
  virtual Real compute_residual_dual_norm(const unsigned int N) libmesh_override;

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution. This function does not used the cached
   * data and therefore also works when the parameter changes as
   * a function of time.
   */
  virtual Real uncached_compute_residual_dual_norm(const unsigned int N);

  /**
   * Helper function for caching the terms in the
   * online residual assembly that do not change in time.
   * (This is only useful when the parameter is fixed in time.)
   */
  void cache_online_residual_terms(const unsigned int N);

  /**
   * Clear all the Riesz representors that are used to compute the RB residual
   * (and hence error bound). This is useful since once we complete the Greedy
   * we may not need the representors any more.
   * Override to clear the M_q representors.
   */
  virtual void clear_riesz_representors() libmesh_override;

  /**
   * Write out all the data to text files in order to segregate the
   * Offline stage from the Online stage.
   * Note: This is a legacy method, use RBDataSerialization instead.
   */
  virtual void legacy_write_offline_data_to_files(const std::string & directory_name = "offline_data",
                                                  const bool write_binary_data=true) libmesh_override;

  /**
   * Read in the saved Offline reduced basis data
   * to initialize the system for Online solves.
   * Note: This is a legacy method, use RBDataSerialization instead.
   */
  virtual void legacy_read_offline_data_from_files(const std::string & directory_name = "offline_data",
                                                   bool read_error_bound_data=true,
                                                   const bool read_binary_data=true) libmesh_override;

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * Dense RB L2 matrix.
   */
  DenseMatrix<Number> RB_L2_matrix;

  /**
   * Cached data for subsequent solves.
   */
  DenseMatrix<Number> RB_LHS_matrix;
  DenseMatrix<Number> RB_RHS_matrix;
  DenseVector<Number> RB_RHS_save;

  /**
   * Dense matrices for the RB mass matrices.
   */
  std::vector< DenseMatrix<Number> > RB_M_q_vector;

  /**
   * The RB outputs for all time-levels from the
   * most recent rb_solve.
   */
  std::vector< std::vector<Number> > RB_outputs_all_k;

  /**
   * The error bounds for each RB output for all
   * time-levels from the most recent rb_solve.
   */
  std::vector< std::vector<Real> > RB_output_error_bounds_all_k;

  /**
   * The RB solution at the previous time-level.
   */
  DenseVector<Number> old_RB_solution;

  /**
   * Array storing the solution data at each time level from the most recent solve.
   */
  std::vector< DenseVector<Number> > RB_temporal_solution_data;

  /**
   * The error bound data for all time-levels from the
   * most recent rb_solve.
   */
  std::vector< Real > error_bound_all_k;

  /**
   * Vector storing initial L2 error for all
   * 1 <= N <= RB_size.
   */
  std::vector<Real> initial_L2_error_all_N;

  /**
   * The RB initial conditions (i.e. L2 projection of the truth
   * initial condition) for each N.
   */
  std::vector< DenseVector<Number> > RB_initial_condition_all_N;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   */
  std::vector< std::vector< std::vector<Number> > > Fq_Mq_representor_innerprods;
  std::vector< std::vector< std::vector<Number> > > Mq_Mq_representor_innerprods;
  std::vector< std::vector< std::vector< std::vector<Number> > > > Aq_Mq_representor_innerprods;


  /**
   * Cached residual terms. These can be used to accelerate residual calculations
   * when we have an LTI system.
   */
  Number cached_Fq_term;
  DenseVector<Number> cached_Fq_Aq_vector;
  DenseMatrix<Number> cached_Aq_Aq_matrix;
  DenseVector<Number> cached_Fq_Mq_vector;
  DenseMatrix<Number> cached_Aq_Mq_matrix;
  DenseMatrix<Number> cached_Mq_Mq_matrix;

  /**
   * Vector storing the mass matrix representors.
   * These are basis dependent and hence stored here.
   */
  std::vector< std::vector< NumericVector<Number> * > > M_q_representor;

  /**
   * Check that the data has been cached in case of using rb_solve_again
   */
  bool _rb_solve_data_cached;
};

}

#endif // LIBMESH_TRANSIENT_RB_EVALUATION_H
