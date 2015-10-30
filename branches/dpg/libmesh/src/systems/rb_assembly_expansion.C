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

// rbOOmit includes
#include "rb_assembly_expansion.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBAssemblyExpansion implementation

RBAssemblyExpansion::RBAssemblyExpansion()
{
  // Clear the assembly vectors so that we can push_back
  A_q_assembly_vector.clear();
  F_q_assembly_vector.clear();

  output_assembly_vector.clear();
}

void RBAssemblyExpansion::attach_A_q_assembly(ElemAssembly* A_q_assembly)
{
  A_q_assembly_vector.push_back(A_q_assembly);
}

void RBAssemblyExpansion::attach_F_q_assembly(ElemAssembly* F_q_assembly)
{
  F_q_assembly_vector.push_back(F_q_assembly);
}

void RBAssemblyExpansion::attach_output_assembly(std::vector<ElemAssembly*> output_assembly)
{
  output_assembly_vector.push_back(output_assembly);
}

void RBAssemblyExpansion::attach_output_assembly(ElemAssembly* output_assembly)
{
  std::vector<ElemAssembly*> L_vector(1); L_vector[0] = output_assembly;

  attach_output_assembly(L_vector);
}


}