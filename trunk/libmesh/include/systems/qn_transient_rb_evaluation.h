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

#ifndef __qn_transient_rb_evaluation_h__
#define __qn_transient_rb_evaluation_h__

// Configuration data
#include "libmesh_config.h"

// This class requires QNTransientRBSystem, which is not
// defined if SLEPc is not present.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "transient_rb_evaluation.h"
#include "qn_transient_rb_system.h"

namespace libMesh
{
        
class RBSystem;

/**
 * This class is part of the rbOOmit framework.
 *
 * QNTransientRBEvaluation extends TransientRBEvaluation to
 * encapsulates the code and data required
 * to perform "online" RB evaluations for quadratically
 * nonlinear transient problems.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// QNTransientRBEvaluation class definition

class QNTransientRBEvaluation : public TransientRBEvaluation
{
public:

  /**
   * Constructor.
   */
  QNTransientRBEvaluation (QNTransientRBSystem& rb_sys_in);

  /**
   * The type of the parent.
   */
  typedef TransientRBEvaluation Parent;

  /**
   * Clear this RBEvaluation object.
   * Overload in subclasses to clear any extra data.
   */
  virtual void clear();

  /**
   * Initialize this object by allocating the necessary data fields.
   */
  virtual void initialize();

  /**
   * Perform online solve for current_params
   * with the N basis functions. Overload this
   * to solve the nonlinear RB system using
   * Newton's method.
   */
  virtual Real RB_solve(unsigned int N);

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution_vector.
   */
  virtual Real compute_residual_dual_norm(const unsigned int N);

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
   * Storage for the reduced basis trilinear form that arises
   * from the quadratic nonlinearity.
   */
  std::vector< std::vector< std::vector<Number> > > RB_trilinear_form;


  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   */
  std::vector< std::vector< std::vector<Number> > > Fq_C_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > Mq_C_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > Aq_C_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > C_C_representor_norms;

};

}

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
