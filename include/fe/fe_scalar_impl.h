// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_FE_SCALAR_IMPL_H
#define LIBMESH_FE_SCALAR_IMPL_H

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{

// ------------------------------------------------------------
// SCALAR-specific implementations

// Anonymous namespace for local helper functions
namespace {

template <typename RealType>
void scalar_nodal_soln(const ElemTempl<RealType> * elem,
                       const Order order,
                       const std::vector<Number> & elem_soln,
                       std::vector<Number> &       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  nodal_soln.resize(n_nodes);

  // If the SCALAR order is CONSTANT, just set the nodal values
  // to zero, otherwise, set to the value of the first SCALAR dof
  for (unsigned int i=0; i<n_nodes; i++)
    nodal_soln[i] = (order == CONSTANT) ? 0. : elem_soln[0];
} // scalar_nodal_soln()

} // anonymous namespace




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,SCALAR,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <typename RealType>
void FEShim<1,SCALAR,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <typename RealType>
void FEShim<2,SCALAR,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <typename RealType>
void FEShim<3,SCALAR,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

// Full specialization of n_dofs() function for every dimension
// The Order indicates the number of SCALAR dofs
template <typename RealType> unsigned int FEShim<0,SCALAR,RealType>::n_dofs(const ElemType, const Order o) { return o; }
template <typename RealType> unsigned int FEShim<1,SCALAR,RealType>::n_dofs(const ElemType, const Order o) { return o; }
template <typename RealType> unsigned int FEShim<2,SCALAR,RealType>::n_dofs(const ElemType, const Order o) { return o; }
template <typename RealType> unsigned int FEShim<3,SCALAR,RealType>::n_dofs(const ElemType, const Order o) { return o; }

// Full specialization of n_dofs_at_node() function for every dimension.
// SCALARs have no dofs at nodes
template <typename RealType> unsigned int FEShim<0,SCALAR,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<1,SCALAR,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<2,SCALAR,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<3,SCALAR,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
// SCALARs have no dofs per element
template <typename RealType> unsigned int FEShim<0,SCALAR,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <typename RealType> unsigned int FEShim<1,SCALAR,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <typename RealType> unsigned int FEShim<2,SCALAR,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <typename RealType> unsigned int FEShim<3,SCALAR,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Scalar FEMs are discontinuous
template <typename RealType> FEContinuity FEShim<0,SCALAR,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<1,SCALAR,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<2,SCALAR,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<3,SCALAR,RealType>::get_continuity() { return DISCONTINUOUS; }

// Scalar FEMs are not hierarchic
template <typename RealType> bool FEShim<0,SCALAR,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<1,SCALAR,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<2,SCALAR,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<3,SCALAR,RealType>::is_hierarchic() { return false; }


#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() just returns for SCALAR FEMs
template <typename RealType>
void FEShim<2,SCALAR,RealType>::compute_constraints (DofConstraints &,
                                        DofMap &,
                                        const unsigned int,
                                                     const ElemTempl<Real> *)
{ }

template <typename RealType>
void FEShim<3,SCALAR,RealType>::compute_constraints (DofConstraints &,
                                        DofMap &,
                                        const unsigned int,
                                                     const ElemTempl<Real> *)
{ }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Scalar FEM shapes do not need reinit
template <typename RealType> bool FEShim<0,SCALAR,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<1,SCALAR,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<2,SCALAR,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<3,SCALAR,RealType>::shapes_need_reinit() { return false; }

} // namespace libMesh

#endif // LIBMESH_FE_SCALAR_IMPL_H
