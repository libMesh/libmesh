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

#ifndef __rb_evaluation_h__
#define __rb_evaluation_h__

#include "reference_counted_object.h"
#include "dense_matrix.h"
#include "dense_vector.h"

namespace libMesh
{

class RBSystem;
template <typename T> class NumericVector;
        
/**
 * This class is part of the rbOOmit framework.
 *
 * RBEvaluation encapsulates the code and data required
 * to perform "online" RB evaluations.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBEvaluation class definition

class RBEvaluation : public ReferenceCountedObject<RBEvaluation>
{
public:

  /**
   * Constructor.
   */
  RBEvaluation (RBSystem& rb_sys_in);

  /**
   * Destructor.
   */
  virtual ~RBEvaluation ();

  /**
   * Clear this RBEvaluation object. Delete the basis functions
   * and clear and extra data in subclasses.
   */
  virtual void clear();

  /**
   * Initialize this object by allocating the necessary data fields.
   */
  virtual void initialize();

  /**
   * Perform online solve with the N RB basis functions, for the
   * set of parameters in current_params, where 1 <= N <= RB_size.
   */
  virtual Real RB_solve(unsigned int N);

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution_vector.
   */
  virtual Real compute_residual_dual_norm(const unsigned int N);

  /**
   * Get the current number of basis functions.
   */
  virtual unsigned int get_n_basis_functions() const { return basis_functions.size(); }

  /**
   * Set the number of basis functions. Useful when reading in
   * stored data.
   */
  virtual void set_n_basis_functions(unsigned int n_bfs) { basis_functions.resize(n_bfs); }

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

  /**
   * If store_basis_functions=true, write out all the basis functions to file.
   * Precision level specifies the number of significant digits to write out.
   */
  virtual void write_out_basis_functions(const std::string& directory_name, const unsigned int precision_level);
  
  /**
   * If store_basis_functions=true, read in all the basis functions from file.
   */
  virtual void read_in_basis_functions(const std::string& directory_name);
  
  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * A reference to the associated RBSystem.
   */
  RBSystem& rb_sys;

  /**
   * The libMesh vectors storing the finite element coefficients
   * of the RB basis functions.
   */
  std::vector< NumericVector<Number>* > basis_functions;

  /**
   * The inner product matrix. This should be close to the identity,
   * we need to calculate this rather than assume diagonality in order
   * to accurately perform projections since orthogonality degrades
   * with increasing N.
   */
  DenseMatrix<Number> RB_inner_product_matrix;

  /**
   * Dense matrices for the RB computations.
   */
  std::vector< DenseMatrix<Number> > RB_A_q_vector;

  /**
   * Dense vector for the RHS.
   */
  std::vector< DenseVector<Number> > RB_F_q_vector;

  /**
   * The RB solution vector.
   */
  DenseVector<Number> RB_solution;

  /**
   * The vectors storing the RB output vectors.
   */
  std::vector< std::vector< DenseVector<Number> > > RB_output_vectors;

  /**
   * The vectors storing the RB output values and
   * corresponding error bounds.
   */
  std::vector< Number > RB_outputs;
  std::vector< Real > RB_output_error_bounds;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   */
  std::vector<Number> Fq_representor_norms;
  std::vector< std::vector< std::vector<Number> > > Fq_Aq_representor_norms;
  std::vector< std::vector< std::vector<Number> > > Aq_Aq_representor_norms;

};

}

#endif