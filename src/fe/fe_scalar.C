// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/dof_map.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{

// ------------------------------------------------------------
// SCALAR-specific implementations

// Anonymous namespace for local helper functions
namespace {

void scalar_nodal_soln(const Elem * elem,
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
template <>
void FE<0,SCALAR>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<1,SCALAR>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<2,SCALAR>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<3,SCALAR>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ scalar_nodal_soln(elem, order, elem_soln, nodal_soln); }

// Full specialization of n_dofs() function for every dimension
// The Order indicates the number of SCALAR dofs
template <> unsigned int FE<0,SCALAR>::n_dofs(const ElemType, const Order o) { return o; }
template <> unsigned int FE<1,SCALAR>::n_dofs(const ElemType, const Order o) { return o; }
template <> unsigned int FE<2,SCALAR>::n_dofs(const ElemType, const Order o) { return o; }
template <> unsigned int FE<3,SCALAR>::n_dofs(const ElemType, const Order o) { return o; }

// Full specialization of n_dofs_at_node() function for every dimension.
// SCALARs have no dofs at nodes
template <> unsigned int FE<0,SCALAR>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,SCALAR>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,SCALAR>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,SCALAR>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
// SCALARs have no dofs per element
template <> unsigned int FE<0,SCALAR>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<1,SCALAR>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<2,SCALAR>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<3,SCALAR>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Scalar FEMs are discontinuous
template <> FEContinuity FE<0,SCALAR>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,SCALAR>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,SCALAR>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,SCALAR>::get_continuity() const { return DISCONTINUOUS; }

// Scalar FEMs are not hierarchic
template <> bool FE<0,SCALAR>::is_hierarchic() const { return false; }
template <> bool FE<1,SCALAR>::is_hierarchic() const { return false; }
template <> bool FE<2,SCALAR>::is_hierarchic() const { return false; }
template <> bool FE<3,SCALAR>::is_hierarchic() const { return false; }


#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() just returns for SCALAR FEMs
template <>
void FE<2,SCALAR>::compute_constraints (DofConstraints &,
                                        DofMap &,
                                        const unsigned int,
                                        const Elem *)
{ }

template <>
void FE<3,SCALAR>::compute_constraints (DofConstraints &,
                                        DofMap &,
                                        const unsigned int,
                                        const Elem *)
{ }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Scalar FEM shapes do not need reinit
template <> bool FE<0,SCALAR>::shapes_need_reinit() const { return false; }
template <> bool FE<1,SCALAR>::shapes_need_reinit() const { return false; }
template <> bool FE<2,SCALAR>::shapes_need_reinit() const { return false; }
template <> bool FE<3,SCALAR>::shapes_need_reinit() const { return false; }

} // namespace libMesh
