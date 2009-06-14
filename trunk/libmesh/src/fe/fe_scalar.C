// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local includes
#include "dof_map.h"
#include "fe.h"
#include "fe_macro.h"
#include "elem.h"




// ------------------------------------------------------------
// SCALAR-specific implementations

template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<Number>& elem_soln,
			   std::vector<Number>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  nodal_soln.resize(n_nodes);

  // If the SCALAR order is CONSTANT, just set the nodal values
  // to zero, otherwise, set to the value of the first SCALAR dof
  for(unsigned int i=0; i<n_nodes; i++)
  {
    nodal_soln[i] = (order == CONSTANT) ? 0. : elem_soln[0];
  }
}

template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType, const Order o)
{
  // The Order indicates the number of SCALAR dofs
  return o;
}

template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType,
				       const Order,
				       const unsigned int)
{
  // SCALARs have no dofs at nodes

  return 0;
}

template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType,
					const Order)
{
  // SCALARs have no dofs per element

  return 0;
}

template <unsigned int Dim, FEFamily T>
FEContinuity FE<Dim,T>::get_continuity() const
{
  // This doesn't really make sense for a SCALAR...
  return C_ZERO;
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::is_hierarchic() const
{
  return false;
}

#ifdef LIBMESH_ENABLE_AMR
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_constraints (DofConstraints &constraints,
				     DofMap &dof_map,
				     const unsigned int variable_number,
				     const Elem* elem)
{
  return;
}
#endif // #ifdef LIBMESH_ENABLE_AMR

template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  return false;
}

//--------------------------------------------------------------
// Explicit instantiation of member functions
INSTANTIATE_MBRF(0,SCALAR);
INSTANTIATE_MBRF(1,SCALAR);
INSTANTIATE_MBRF(2,SCALAR);
INSTANTIATE_MBRF(3,SCALAR);