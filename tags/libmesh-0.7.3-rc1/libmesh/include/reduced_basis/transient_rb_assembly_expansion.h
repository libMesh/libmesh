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

#ifndef __transient_rb_assembly_expansion_h__
#define __transient_rb_assembly_expansion_h__

// libMesh includes
#include "rb_assembly_expansion.h"


namespace libMesh
{

/**
 * This extends RBAssemblyExpansion to provide an assembly
 * expansion for the case of time-dependent PDEs.
 * This just requires an extra set of ElemAssembly functors
 * for the time-derivative term.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// TransientRBAssemblyExpansion class definition
class TransientRBAssemblyExpansion : public RBAssemblyExpansion
{
public:

  /**
   * Constructor.
   */
  TransientRBAssemblyExpansion();

  /**
   * Attach ElemAssembly object for the time-derivative
   * (both interior and boundary assembly).
   */
  void attach_M_q_assembly(ElemAssembly* A_q_assembly);

  // -------- Data members --------

  /**
   * Vectors storing the function pointers to the assembly
   * routines for the time-derivative operators, both interior and boundary
   * assembly.
   */
  std::vector<ElemAssembly*> M_q_assembly_vector;

};

}

#endif