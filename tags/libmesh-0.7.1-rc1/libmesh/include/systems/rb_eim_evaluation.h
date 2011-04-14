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

#ifndef __rb_eim_evaluation_h__
#define __rb_eim_evaluation_h__

#include "rb_evaluation.h"
#include "point.h"

namespace libMesh
{
        
/**
 * This class is part of the rbOOmit framework.
 *
 * RBEIMEvaluation extends RBEvaluation to
 * encapsulate the code and data required
 * to perform "online" RB evaluations for transient problems.
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
  RBEIMEvaluation (RBSystem& rb_sys_in);

  /**
   * The type of the parent.
   */
  typedef RBEvaluation Parent;

  /**
   * Clear this object.
   */
  virtual void clear();

  /**
   * Initialize this object by allocating the necessary data fields.
   */
  virtual void initialize();

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

};

}

#endif
