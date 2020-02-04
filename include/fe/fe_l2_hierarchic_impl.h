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

#ifndef LIBMESH_L2_HIERARCHIC_IMPL_H
#define LIBMESH_L2_HIERARCHIC_IMPL_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// ------------------------------------------------------------
// Hierarchic-specific implementations

// Anonymous namespace for local helper functions
namespace {

void l2_hierarchic_nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> &       nodal_soln,
                              unsigned Dim)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(totalorder, L2_HIERARCHIC);

  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
        libmesh_assert_equal_to (elem_soln.size(), 1);

        const Number val = elem_soln[0];

        for (unsigned int n=0; n<n_nodes; n++)
          nodal_soln[n] = val;

        return;
      }


      // For other orders do interpolation at the nodes
      // explicitly.
    default:
      {

        const unsigned int n_sf =
          // FEShim<Dim,T,RealType>::n_shape_functions(elem_type, totalorder);
          FEInterface::n_shape_functions(Dim, fe_type, elem_type);

        std::vector<Point> refspace_nodes;
        FEBase::get_refspace_nodes(elem_type,refspace_nodes);
        libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);

        for (unsigned int n=0; n<n_nodes; n++)
          {
            libmesh_assert_equal_to (elem_soln.size(), n_sf);

            // Zero before summation
            nodal_soln[n] = 0;

            // u_i = Sum (alpha_i phi_i)
            for (unsigned int i=0; i<n_sf; i++)
              nodal_soln[n] += elem_soln[i] *
                // FEShim<Dim,T,RealType>::shape(elem, order, i, mapped_point);
                FEInterface::shape(Dim, fe_type, elem, i, refspace_nodes[n]);
          }

        return;
      }
    }
} // l2_hierarchic_nodal_soln()




unsigned int l2_hierarchic_n_dofs(const ElemType t, const Order o)
{
  libmesh_assert_greater (o, 0);
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
    case EDGE3:
      return (o+1);
    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      return ((o+1)*(o+1));
    case HEX8:
    case HEX20:
    case HEX27:
      return ((o+1)*(o+1)*(o+1));
    case TRI3:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case TRI6:
      return ((o+1)*(o+2)/2);
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for L2_HIERARCHIC FE family!");
    }
} // l2_hierarchic_n_dofs()


} // anonymous namespace


  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,L2_HIERARCHIC,RealType>::nodal_soln(const Elem * elem,
                                     const Order order,
                                     const std::vector<Number> & elem_soln,
                                     std::vector<Number> & nodal_soln)
{ l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <typename RealType>
void FEShim<1,L2_HIERARCHIC,RealType>::nodal_soln(const Elem * elem,
                                     const Order order,
                                     const std::vector<Number> & elem_soln,
                                     std::vector<Number> & nodal_soln)
{ l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <typename RealType>
void FEShim<2,L2_HIERARCHIC,RealType>::nodal_soln(const Elem * elem,
                                     const Order order,
                                     const std::vector<Number> & elem_soln,
                                     std::vector<Number> & nodal_soln)
{ l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,L2_HIERARCHIC,RealType>::nodal_soln(const Elem * elem,
                                     const Order order,
                                     const std::vector<Number> & elem_soln,
                                     std::vector<Number> & nodal_soln)
{ l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }

// Full specialization of n_dofs() function for every dimension
template <typename RealType> unsigned int FEShim<0,L2_HIERARCHIC,RealType>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,L2_HIERARCHIC,RealType>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,L2_HIERARCHIC,RealType>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,L2_HIERARCHIC,RealType>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
// Discontinuous L2 elements only have interior nodes
template <typename RealType> unsigned int FEShim<0,L2_HIERARCHIC,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<1,L2_HIERARCHIC,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<2,L2_HIERARCHIC,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<3,L2_HIERARCHIC,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <typename RealType> unsigned int FEShim<0,L2_HIERARCHIC,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,L2_HIERARCHIC,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,L2_HIERARCHIC,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,L2_HIERARCHIC,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }

// L2 Hierarchic FEMs are C^0 continuous
template <typename RealType> FEContinuity FEShim<0,L2_HIERARCHIC,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<1,L2_HIERARCHIC,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<2,L2_HIERARCHIC,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<3,L2_HIERARCHIC,RealType>::get_continuity() { return DISCONTINUOUS; }

// L2 Hierarchic FEMs are hierarchic (duh!)
template <typename RealType> bool FEShim<0,L2_HIERARCHIC,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<1,L2_HIERARCHIC,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<2,L2_HIERARCHIC,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<3,L2_HIERARCHIC,RealType>::is_hierarchic() { return true; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() is a NOOP for DISCONTINUOUS FE's
template <typename RealType>
void FEShim<2,L2_HIERARCHIC,RealType>::compute_constraints (DofConstraints &,
                                               DofMap &,
                                               const unsigned int,
                                                            const ElemTempl<Real> *)
{ }

template <typename RealType>
void FEShim<3,L2_HIERARCHIC,RealType>::compute_constraints (DofConstraints &,
                                               DofMap &,
                                               const unsigned int,
                                                            const ElemTempl<Real> *)
{ }
#endif // #ifdef LIBMESH_ENABLE_AMR

// L2-Hierarchic FEM shapes need reinit
template <typename RealType> bool FEShim<0,L2_HIERARCHIC,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<1,L2_HIERARCHIC,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<2,L2_HIERARCHIC,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<3,L2_HIERARCHIC,RealType>::shapes_need_reinit() { return true; }

} // namespace libMesh

#endif // LIBMESH_L2_HIERARCHIC_IMPL_H
