
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

#ifndef LIBMESH_RB_SCM_EVALUATION_H
#define LIBMESH_RB_SCM_EVALUATION_H

// RBSCMEvaluation should only be available
// if SLEPc and GLPK support is enabled.
#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

// rbOOmit includes
#include "libmesh/rb_parametrized.h"

// libMesh includes
#include "libmesh/parallel_object.h"

// C++ includes

namespace libMesh
{

// Forward declarations
class RBThetaExpansion;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBSCMEvaluation encapsulates the functionality required
 * to _evaluate_ the Successive Constraint Method for
 * associated with a reduced basis model.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class RBSCMEvaluation : public RBParametrized,
                        public ParallelObject
{
public:

  /**
   * Constructor.
   */
  RBSCMEvaluation (const Parallel::Communicator & comm
                   LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor.
   */
  virtual ~RBSCMEvaluation ();

  /**
   * Set the RBThetaExpansion object.
   */
  void set_rb_theta_expansion(RBThetaExpansion & rb_theta_expansion_in);

  /**
   * Get a reference to the rb_theta_expansion.
   */
  RBThetaExpansion & get_rb_theta_expansion();

  /**
   * Evaluate single SCM lower bound.
   */
  virtual Real get_SCM_LB();

  /**
   * Evaluate single SCM upper bound.
   */
  virtual Real get_SCM_UB();

  /**
   * Get stability constraints (i.e. the values of coercivity/
   * inf-sup/stability constants at the parameter values chosen
   * during the greedy); we store one constraint for each element
   * of C_J.
   */
  Real get_C_J_stability_constraint(unsigned int j) const;

  /**
   * Set stability constraints (i.e. the values of coercivity/
   * inf-sup/stability constants at the parameter values chosen
   * during the greedy); we store one constraint for each element
   * of C_J.
   */
  void set_C_J_stability_constraint(unsigned int j, Real stability_constraint_in);

  /**
   * Get entries of SCM_UB_vector, which stores the
   * vector y, corresponding to the minimizing eigenvectors
   * for the elements of C_J.
   */
  Real get_SCM_UB_vector(unsigned int j, unsigned int q);

  /**
   * Set entries of SCM_UB_vector, which stores the
   * vector y, corresponding to the minimizing eigenvectors
   * for the elements of C_J.
   */
  void set_SCM_UB_vector(unsigned int j, unsigned int q, Real y_q);

  /**
   * Get size of the set C_J.
   */
  unsigned int get_C_J_size()
  { return cast_int<unsigned int>(C_J.size()); }

  /**
   * Get entry of C_J.
   */
  const RBParameters & get_C_J_entry(unsigned int j);

  /**
   * Get entry of C_J_stability_vector.
   */
  Real get_C_J_stability_value(unsigned int j) { return C_J_stability_vector[j]; }

  /**
   * Get B_min and B_max.
   */
  Real get_B_min(unsigned int i) const;
  Real get_B_max(unsigned int i) const;

  /**
   * Set B_min and B_max.
   */
  void set_B_min(unsigned int i, Real B_min_val);
  void set_B_max(unsigned int i, Real B_max_val);

  /**
   * Helper function to save current_parameters in
   * saved_parameters.
   */
  virtual void save_current_parameters();

  /**
   * Helper functiont to (re)load current_parameters
   * from saved_parameters.
   */
  virtual void reload_current_parameters();

  /**
   * Set parameters based on values saved in "C_J"
   */
  virtual void set_current_parameters_from_C_J(unsigned int C_J_index);

  /**
   * Write out all the data to text files in order to segregate the
   * Offline stage from the Online stage.
   * Note: This is a legacy method, use RBDataSerialization instead.
   */
  virtual void legacy_write_offline_data_to_files(const std::string & directory_name = "offline_data",
                                                  const bool write_binary_data = true);

  /**
   * Read in the saved Offline reduced basis data
   * to initialize the system for Online solves.
   * Note: This is a legacy method, use RBDataSerialization instead.
   */
  virtual void legacy_read_offline_data_from_files(const std::string & directory_name = "offline_data",
                                                   const bool read_binary_data = true);

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * B_min, B_max define the bounding box.
   */
  std::vector<Real> B_min;
  std::vector<Real> B_max;

  /**
   * Vector storing the greedily selected parameters
   * during SCM training.
   */
  std::vector< RBParameters > C_J;

  /**
   * Vector storing the (truth) stability values
   * at the parameters in C_J.
   */
  std::vector<Real> C_J_stability_vector;

  /**
   * This matrix stores the infimizing vectors
   * y_1(\f$ \mu \f$),...,y_Q_a(\f$ \mu \f$), for each \f$ \mu \f$ in
   * C_J, which are used in computing the SCM
   * upper bounds.
   */
  std::vector< std::vector<Real> > SCM_UB_vectors;

private:

  /**
   * Vector in which to save a parameter set. Useful
   * in get_SCM_LB, for example.
   */
  RBParameters saved_parameters;

  /**
   * A pointer to to the object that stores the theta expansion.
   * This is not a UniquePtr since we may want to share it.
   * (Note: a shared_ptr would be a good option here.)
   */
  RBThetaExpansion * rb_theta_expansion;

};

}

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif // LIBMESH_RB_SCM_EVALUATION_H
