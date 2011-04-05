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

#ifndef __transient_rb_evaluation_h__
#define __transient_rb_evaluation_h__

#include "rb_evaluation.h"
#include "transient_rb_system.h"

namespace libMesh
{
        
/**
 * This class is part of the rbOOmit framework.
 *
 * TransientRBEvaluation extends RBEvaluation to
 * encapsulates the code and data required
 * to perform "online" RB evaluations for transient problems.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// QNTransientRBEvaluation class definition

class TransientRBEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor.
   */
  TransientRBEvaluation (TransientRBSystem& rb_sys_in);

  /**
   * The type of the parent.
   */
  typedef RBEvaluation Parent;

  /**
   * Clear this TransientRBEvaluation object.
   * Overload to also clear the M_q representors
   */
  virtual void clear();

  /**
   * Initialize this object by allocating the necessary data fields.
   */
  virtual void initialize();

  /**
   * Perform online solve for current_params
   * with the N basis functions. Overloaded
   * to perform a time-dependent solve.
   */
  virtual Real RB_solve(unsigned int N);

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution. This function uses the cached time-independent
   * data.
   */
  virtual Real compute_residual_dual_norm(const unsigned int N);
  
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
   * Dense RB L2 matrix.
   */
  DenseMatrix<Number> RB_L2_matrix;

  /**
   * Dense matrices for the RB mass matrices.
   */
  std::vector< DenseMatrix<Number> > RB_M_q_vector;

  /**
   * The RB outputs for all time-levels from the
   * most recent RB_solve.
   */
  std::vector< std::vector<Number> > RB_outputs_all_k;

  /**
   * The error bounds for each RB output for all
   * time-levels from the most recent RB_solve.
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
   * most recent RB_solve.
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
  std::vector< std::vector< std::vector<Number> > > Fq_Mq_representor_norms;
  std::vector< std::vector< std::vector<Number> > > Mq_Mq_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > Aq_Mq_representor_norms;
  

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
  std::vector< std::vector< NumericVector<Number>* > > M_q_representor;

};

}

#endif
