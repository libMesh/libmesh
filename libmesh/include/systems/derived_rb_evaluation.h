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

#ifndef __derived_rb_evaluation_h__
#define __derived_rb_evaluation_h__

#include "dense_vector.h"
#include "rb_evaluation.h"

namespace libMesh
{

class RBSystem;
class string;

/**
 * This class is part of the rbOOmit framework.
 *
 * DerivedRBEvaluation encapsulates the code and data required
 * to perform "online" evaluations with a DerivedRBSystem.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// DerivedRBEvaluation class definition

template<class Base>
class DerivedRBEvaluation : public Base
{
public:

  /**
   * Constructor.
   */
  DerivedRBEvaluation (RBSystem& rb_sys_in);

  /**
   * Get the current number of basis functions.
   */
  virtual unsigned int get_n_basis_functions() const;

  /**
   * Set the number of basis functions. Useful when reading in
   * stored data.
   */
  virtual void set_n_basis_functions(unsigned int n_bfs);

  /**
   * If store_basis_functions=true, write out all the basis functions to file.
   * Precision level specifies the number of significant digits to write out.
   */
  virtual void write_out_basis_functions(const std::string& directory_name,
                                        const unsigned int precision_level);
  
  /**
   * If store_basis_functions=true, read in all the basis functions from file.
   */
  virtual void read_in_basis_functions(const std::string& directory_name);
  
  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The dense vectors that define the derived basis
   * functions based on the uber system.
   */
   std::vector< DenseVector<Number> > derived_basis_functions;

};

// And introduce convenient typedefs
typedef DerivedRBEvaluation<RBEvaluation> SteadyDerivedRBEvaluation;

}

#endif