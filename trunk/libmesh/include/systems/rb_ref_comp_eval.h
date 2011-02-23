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

#ifndef __rb_ref_comp_eval_h__
#define __rb_ref_comp_eval_h__

#include "rb_evaluation.h"

namespace libMesh
{
        
/**
 * This class is part of the rbOOmit framework.
 *
 * RBReferenceComponentEvaluation extends RBEvaluation to
 * also encapsulate the "extra data" required to assemble
 * the global matrix in an RBComponentSystem.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBReferenceComponentEvaluation class definition

class RBReferenceComponentEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor.
   */
  RBReferenceComponentEvaluation (RBSystem& rb_sys_in);

  /**
   * The type of the parent.
   */
  typedef RBEvaluation Parent;

  /**
   * Initialize this object by allocating the necessary data fields.
   */
  virtual void initialize();

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
   * Vectors storing the data necessary to assemble the local
   * contribution to the global matrix
   * using the reduced basis approximation.
   */
  std::vector< std::vector< std::vector<Number> > > Aq_g_g;
  std::vector< std::vector< std::vector<Number> > > Aq_bubble_g;

};

}

#endif