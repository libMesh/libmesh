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

#ifndef LIBMESH_RB_EVALUATION_H
#define LIBMESH_RB_EVALUATION_H

// rbOOmit includes
#include "libmesh/rb_parametrized.h"

// libMesh includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/parallel_object.h"

// C++ includes

namespace libMesh
{

class System;
template <typename T> class NumericVector;
class RBThetaExpansion;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBEvaluation encapsulates the functionality required
 * to _evaluate_ a given reduced basis model.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class RBEvaluation : public RBParametrized,
                     public ParallelObject
{
public:

  /**
   * Constructor.
   */
  RBEvaluation (const Parallel::Communicator & comm
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~RBEvaluation ();

  /**
   * Clear this RBEvaluation object. Delete the basis functions
   * and clear and extra data in subclasses.
   */
  virtual void clear() libmesh_override;

  /**
   * Set the RBThetaExpansion object.
   */
  void set_rb_theta_expansion(RBThetaExpansion & rb_theta_expansion_in);

  /**
   * Get a reference to the rb_theta_expansion.
   */
  RBThetaExpansion & get_rb_theta_expansion();

  /**
   * @return true if the theta expansion has been initialized.
   */
  bool is_rb_theta_expansion_initialized() const;

  /**
   * Resize and clear the data vectors corresponding to the
   * value of \p Nmax. Optionally resize the data structures
   * required for the error bound.
   * Overload to also clear and resize any extra
   * data in subclasses.
   */
  virtual void resize_data_structures(const unsigned int Nmax,
                                      bool resize_error_bound_data=true);

  /**
   * Get a reference to the i^th basis function.
   */
  NumericVector<Number> & get_basis_function(unsigned int i);

  /**
   * Perform online solve with the N RB basis functions, for the
   * set of parameters in current_params, where 0 <= N <= RB_size.
   * @return the (absolute) error bound associated with
   * the RB approximation.
   * With an empty RB space (N=0), our RB solution is zero, but we
   * still obtain a meaningful error bound associated with the
   * forcing terms.
   */
  virtual Real rb_solve(unsigned int N);

  /**
   * @return a scaling factor that we can use to provide a consistent
   * scaling of the RB error bound across different parameter values.
   */
  virtual Real get_error_bound_normalization();

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution_vector.
   */
  virtual Real compute_residual_dual_norm(const unsigned int N);

  /**
   * Specifies the residual scaling on the denominator to
   * be used in the a posteriori error bound. Overload
   * in subclass in order to obtain the desired error bound.
   */
  virtual Real residual_scaling_denom(Real alpha_LB);

  /**
   * Evaluate the dual norm of output \p n
   * for the current parameters.
   */
  Real eval_output_dual_norm(unsigned int n, const RBParameters & mu);

  /**
   * Get a lower bound for the stability constant (e.g. coercivity constant or
   * inf-sup constant) at the current parameter value.
   */
  virtual Real get_stability_lower_bound();

  /**
   * Get the current number of basis functions.
   */
  virtual unsigned int get_n_basis_functions() const
  { return cast_int<unsigned int>(basis_functions.size()); }

  /**
   * Set the number of basis functions. Useful when reading in
   * stored data.
   */
  virtual void set_n_basis_functions(unsigned int n_bfs) { basis_functions.resize(n_bfs); }

  /**
   * Clear all the Riesz representors that are used to compute the RB residual
   * (and hence error bound). This is useful since once we complete the Greedy
   * we may not need the representors any more.
   */
  virtual void clear_riesz_representors();

  /**
   * Write out all the data to text files in order to segregate the
   * Offline stage from the Online stage.
   * Note: This is a legacy method, use RBDataSerialization instead.
   */
  virtual void legacy_write_offline_data_to_files(const std::string & directory_name = "offline_data",
                                                  const bool write_binary_data=true);

  /**
   * Read in the saved Offline reduced basis data
   * to initialize the system for Online solves.
   * Note: This is a legacy method, use RBDataSerialization instead.
   */
  virtual void legacy_read_offline_data_from_files(const std::string & directory_name = "offline_data",
                                                   bool read_error_bound_data=true,
                                                   const bool read_binary_data=true);

  /**
   * Write out all the basis functions to file.
   * \p sys is used for file IO
   * \p directory_name specifies which directory to write files to
   * \p read_binary_basis_functions indicates whether to expect
   * binary or ASCII data
   */
  virtual void write_out_basis_functions(System & sys,
                                         const std::string & directory_name = "offline_data",
                                         const bool write_binary_basis_functions = true);

  /**
   * Same as write_out_basis_functions, except in this case we pass in the vectors to be
   * written.
   */
  virtual void write_out_vectors(System & sys,
                                 std::vector<NumericVector<Number> *> & vectors,
                                 const std::string & directory_name = "offline_data",
                                 const std::string & data_name = "bf",
                                 const bool write_binary_basis_functions = true);

  /**
   * Read in all the basis functions from file.
   * \p sys is used for file IO
   * \p directory_name specifies which directory to write files to
   * \p read_binary_basis_functions indicates whether to expect
   * binary or ASCII data
   */
  virtual void read_in_basis_functions(System & sys,
                                       const std::string & directory_name = "offline_data",
                                       const bool read_binary_basis_functions = true);

  /**
   * Same as read_in_basis_functions, except in this case we pass in the vectors to be
   * written. We assume that the size of vectors indicates the number of vectors
   * that need to be read in.
   */
  void read_in_vectors(System & sys,
                       std::vector<NumericVector<Number> *> & vectors,
                       const std::string & directory_name,
                       const std::string & data_name,
                       const bool read_binary_vectors);

  /**
   * Performs read_in_vectors for a list of directory names and data names.
   * Reading in vectors requires us to renumber the dofs in a partition-indepdent
   * way. This function only renumbers the dofs once at the start (and reverts
   * it at the end), which can save a lot of work compared to renumbering on every read.
   */
  void read_in_vectors_from_multiple_files(System & sys,
                                           std::vector<std::vector<NumericVector<Number> *> *> multiple_vectors,
                                           const std::vector<std::string> & multiple_directory_names,
                                           const std::vector<std::string> & multiple_data_names,
                                           const bool read_binary_vectors);

  /**
   * Version string that we need to use for writing/reading basis functions.
   */
  static std::string get_io_version_string();

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The libMesh vectors storing the finite element coefficients
   * of the RB basis functions.
   */
  std::vector<NumericVector<Number> *> basis_functions;

  /**
   * The list of parameters selected by the Greedy algorithm in generating
   * the Reduced Basis associated with this RBEvaluation object.
   */
  std::vector<RBParameters> greedy_param_list;

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
  std::vector<DenseMatrix<Number> > RB_Aq_vector;

  /**
   * Dense vector for the RHS.
   */
  std::vector<DenseVector<Number> > RB_Fq_vector;

  /**
   * The RB solution vector.
   */
  DenseVector<Number> RB_solution;

  /**
   * The vectors storing the RB output vectors.
   */
  std::vector<std::vector<DenseVector<Number> > > RB_output_vectors;

  /**
   * The vectors storing the RB output values and
   * corresponding error bounds.
   */
  std::vector< Number > RB_outputs;
  std::vector< Real > RB_output_error_bounds;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   * These values are independent of a basis, hence they can
   * be copied over directly from an RBSystem.
   */
  std::vector<Number> Fq_representor_innerprods;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   * We store the Aq-dependent representor inner products because they depend
   * on a reduced basis space. The basis independent representors
   * are stored in RBSystem.
   */
  std::vector< std::vector< std::vector<Number> > > Fq_Aq_representor_innerprods;
  std::vector< std::vector< std::vector<Number> > > Aq_Aq_representor_innerprods;

  /**
   * The vector storing the dual norm inner product terms
   * for each output.
   * These values are independent of a basis, hence they can
   * be copied over directly from an RBSystem.
   */
  std::vector< std::vector< Number > > output_dual_innerprods;

  /**
   * Vector storing the residual representors associated with the
   * left-hand side.
   * These are basis dependent and hence stored here, whereas
   * the Fq_representors are stored in RBSystem.
   */
  std::vector< std::vector< NumericVector<Number> *> > Aq_representor;

  /**
   * Boolean to indicate whether we evaluate a posteriori error bounds
   * when rb_solve is called.
   */
  bool evaluate_RB_error_bound;

  /**
   * Boolean flag to indicate whether we compute the RB_inner_product_matrix.
   */
  bool compute_RB_inner_product;

protected:

  /**
   * Helper function that checks if \p file_name exists.
   */
  void assert_file_exists(const std::string & file_name);

private:

  /**
   * A pointer to to the object that stores the theta expansion.
   * This is not a UniquePtr since we may want to share it.
   * (Note: a shared_ptr would be a good option here.)
   */
  RBThetaExpansion * rb_theta_expansion;

};

}

#endif // LIBMESH_RB_EVALUATION_H
