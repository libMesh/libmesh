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

#ifndef __rb_assembly_expansion_h__
#define __rb_assembly_expansion_h__

// libMesh includes
#include "elem_assembly.h"
#include "reference_counted_object.h"

// misc includes
#include <vector>


namespace libMesh
{

/**
 * This class stores the set of ElemAssembly functor objects that define
 * the "parameter-independent expansion" of a PDE.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBAssemblyExpansion class definition
class RBAssemblyExpansion : public ReferenceCountedObject<RBAssemblyExpansion>
{
public:

  /**
   * Constructor.
   */
  RBAssemblyExpansion();

  /**
   * Destructor.
   */
  virtual ~RBAssemblyExpansion() {}

  /**
   * Attach ElemAssembly object for the left-hand side
   * (both interior and boundary assembly).
   */
  void attach_A_q_assembly(ElemAssembly* A_q_assembly);

  /**
   * Attach multiple ElemAssembly objects for the left-hand side
   * (both interior and boundary assembly).
   */
  void attach_multiple_A_q_assembly(std::vector<ElemAssembly*> A_q_assembly);

  /**
   * Attach ElemAssembly object for the right-hand side
   * (both interior and boundary assembly).
   */
  void attach_F_q_assembly(ElemAssembly* F_q_assembly);

  /**
   * Attach multiple ElemAssembly objects for the right-hand side
   * (both interior and boundary assembly).
   */
  void attach_multiple_F_q_assembly(std::vector<ElemAssembly*> F_q_assembly);

  /**
   * Attach ElemAssembly object for an output
   * (both interior and boundary assembly).
   * In this case we pass in vector arguments to allow for Q_l > 1.
   */
  virtual void attach_output_assembly(std::vector<ElemAssembly*> output_assembly);

  /**
   * Attach ElemAssembly object for an output
   * (both interior and boundary assembly).
   * This function provides simpler syntax in the case that Q_l = 1; we
   * do not need to use a vector in this case.
   */
  virtual void attach_output_assembly(ElemAssembly* output_assembly);

  // -------- Data members --------

  /**
   * Vectors storing the function pointers to the assembly
   * routines for the affine operators, both interior and boundary
   * assembly.
   */
  std::vector<ElemAssembly*> A_q_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the rhs affine vectors.
   */
  std::vector<ElemAssembly*> F_q_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the outputs. Element interior part.
   */
  std::vector< std::vector<ElemAssembly*> > output_assembly_vector;

};

}

#endif