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
#include "libmesh/transient_rb_assembly_expansion.h"
#include "libmesh/elem_assembly.h"

namespace libMesh
{

// ------------------------------------------------------------
// TransientRBAssemblyExpansion implementation

TransientRBAssemblyExpansion::TransientRBAssemblyExpansion()
{
}

void TransientRBAssemblyExpansion::perform_M_interior_assembly(unsigned int q,
                                                               FEMContext & context)
{
  if(q >= get_n_M_terms())
    libmesh_error_msg("Error: We must have q < get_n_M_terms in perform_M_interior_assembly.");

  libmesh_assert(_M_assembly_vector[q]);

  return _M_assembly_vector[q]->interior_assembly( context );
}

void TransientRBAssemblyExpansion::perform_M_boundary_assembly(unsigned int q,
                                                               FEMContext & context)
{
  if(q >= get_n_M_terms())
    libmesh_error_msg("Error: We must have q < get_n_M_terms in perform_M_boundary_assembly.");

  libmesh_assert(_M_assembly_vector[q]);

  return _M_assembly_vector[q]->boundary_assembly( context );
}

unsigned int TransientRBAssemblyExpansion::get_n_M_terms() const
{
  return cast_int<unsigned int>(_M_assembly_vector.size());
}

void TransientRBAssemblyExpansion::attach_M_assembly(ElemAssembly * M_q_assembly)
{
  _M_assembly_vector.push_back(M_q_assembly);
}

ElemAssembly & TransientRBAssemblyExpansion::get_M_assembly(unsigned int q)
{
  if(q >= get_n_M_terms())
    libmesh_error_msg("Error: We must have q < get_n_M_terms in get_M_assembly.");

  return *_M_assembly_vector[q];
}

}
