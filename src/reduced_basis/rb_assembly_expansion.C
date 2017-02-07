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
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/elem_assembly.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBAssemblyExpansion implementation

RBAssemblyExpansion::RBAssemblyExpansion()
{
}

void RBAssemblyExpansion::perform_A_interior_assembly(unsigned int q,
                                                      FEMContext & context)
{
  if(q >= get_n_A_terms())
    libmesh_error_msg("Error: We must have q < get_n_A_terms in perform_A_interior_assembly.");

  libmesh_assert(_A_assembly_vector[q]);

  return _A_assembly_vector[q]->interior_assembly( context );
}

void RBAssemblyExpansion::perform_A_boundary_assembly(unsigned int q,
                                                      FEMContext & context)
{
  if(q >= get_n_A_terms())
    libmesh_error_msg("Error: We must have q < get_n_A_terms in perform_A_boundary_assembly.");

  libmesh_assert(_A_assembly_vector[q]);

  return _A_assembly_vector[q]->boundary_assembly( context );
}

void RBAssemblyExpansion::perform_F_interior_assembly(unsigned int q,
                                                      FEMContext & context)
{
  if(q >= get_n_F_terms())
    libmesh_error_msg("Error: We must have q < get_n_F_terms in perform_F_interior_assembly.");

  libmesh_assert(_A_assembly_vector[q]);

  return _F_assembly_vector[q]->interior_assembly( context );
}

void RBAssemblyExpansion::perform_F_boundary_assembly(unsigned int q,
                                                      FEMContext & context)
{
  if(q >= get_n_F_terms())
    libmesh_error_msg("Error: We must have q < get_n_F_terms in perform_F_interior_assembly.");

  libmesh_assert(_A_assembly_vector[q]);

  return _F_assembly_vector[q]->boundary_assembly( context );
}

void RBAssemblyExpansion::perform_output_interior_assembly(unsigned int output_index,
                                                           unsigned int q_l,
                                                           FEMContext & context)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_n_output_terms(output_index)) )
    libmesh_error_msg("Error: We must have output_index < n_outputs and " \
                      << "q_l < get_n_output_terms(output_index) in perform_output_interior_assembly.");

  libmesh_assert(_output_assembly_vector[output_index][q_l]);

  return _output_assembly_vector[output_index][q_l]->interior_assembly(context);
}

void RBAssemblyExpansion::perform_output_boundary_assembly(unsigned int output_index,
                                                           unsigned int q_l,
                                                           FEMContext & context)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_n_output_terms(output_index)) )
    libmesh_error_msg("Error: We must have output_index < n_outputs and " \
                      << "q_l < get_n_output_terms(output_index) in perform_output_boundary_assembly.");

  libmesh_assert(_output_assembly_vector[output_index][q_l]);

  return _output_assembly_vector[output_index][q_l]->boundary_assembly(context);
}

unsigned int RBAssemblyExpansion::get_n_A_terms() const
{
  return cast_int<unsigned int>
    (_A_assembly_vector.size());
}

unsigned int RBAssemblyExpansion::get_n_F_terms() const
{
  return cast_int<unsigned int>
    (_F_assembly_vector.size());
}

unsigned int RBAssemblyExpansion::get_n_outputs() const
{
  return cast_int<unsigned int>
    (_output_assembly_vector.size());
}

unsigned int RBAssemblyExpansion::get_n_output_terms(unsigned int index) const
{
  if(index >= get_n_outputs())
    libmesh_error_msg("Error: We must have index < n_outputs in get_Q_l.");

  return cast_int<unsigned int>
    (_output_assembly_vector[index].size());
}

void RBAssemblyExpansion::attach_A_assembly(ElemAssembly * Aq_assembly)
{
  _A_assembly_vector.push_back(Aq_assembly);
}

void RBAssemblyExpansion::attach_multiple_A_assembly(std::vector<ElemAssembly *> Aq_assembly)
{
  for (std::size_t i=0; i<Aq_assembly.size(); i++)
    _A_assembly_vector.push_back(Aq_assembly[i]);
}

void RBAssemblyExpansion::attach_F_assembly(ElemAssembly * Fq_assembly)
{
  _F_assembly_vector.push_back(Fq_assembly);
}

void RBAssemblyExpansion::attach_multiple_F_assembly(std::vector<ElemAssembly *> Fq_assembly)
{
  for (std::size_t i=0; i<Fq_assembly.size(); i++)
    _F_assembly_vector.push_back(Fq_assembly[i]);
}

void RBAssemblyExpansion::attach_output_assembly(std::vector<ElemAssembly *> output_assembly)
{
  _output_assembly_vector.push_back(output_assembly);
}

void RBAssemblyExpansion::attach_output_assembly(ElemAssembly * output_assembly)
{
  std::vector<ElemAssembly *> L_vector(1); L_vector[0] = output_assembly;

  attach_output_assembly(L_vector);
}

ElemAssembly & RBAssemblyExpansion::get_A_assembly(unsigned int q)
{
  if(q >= get_n_A_terms())
    libmesh_error_msg("Error: We must have q < get_n_A_terms in get_A_assembly.");

  return *_A_assembly_vector[q];
}

ElemAssembly & RBAssemblyExpansion::get_F_assembly(unsigned int q)
{
  if(q >= get_n_F_terms())
    libmesh_error_msg("Error: We must have q < get_n_F_terms in get_F_assembly.");

  return *_F_assembly_vector[q];
}

ElemAssembly & RBAssemblyExpansion::get_output_assembly(unsigned int output_index,
                                                        unsigned int q_l)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_n_output_terms(output_index)) )
    libmesh_error_msg("Error: We must have output_index < n_outputs and " \
                      << "q_l < get_n_output_terms(output_index) in get_output_assembly.");

  return *_output_assembly_vector[output_index][q_l];
}

}
