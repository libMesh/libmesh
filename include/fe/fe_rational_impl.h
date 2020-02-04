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

#ifndef LIBMESH_FE_RATIONAL_IMPL_H
#define LIBMESH_FE_RATIONAL_IMPL_H

// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

// Anonymous namespace for local helper functions
namespace {

using namespace libMesh;

#ifndef RATIONAL_UNDERLYING_FAMILY
#define RATIONAL_UNDERLYING_FAMILY

static const FEFamily _underlying_fe_family = BERNSTEIN;

#endif

template <typename RealType>
void rational_nodal_soln(const ElemTempl<RealType> * elem,
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
  FEType fe_type(totalorder, _underlying_fe_family);

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


      // For other bases do interpolation at the nodes
      // explicitly.
    default:
      {
        const unsigned int n_sf =
          FEInterface::n_shape_functions(Dim, fe_type, elem_type);

        std::vector<Point> refspace_nodes;
        FEBase::get_refspace_nodes(elem_type,refspace_nodes);
        libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);
        libmesh_assert_equal_to (n_sf, n_nodes);
        libmesh_assert_equal_to (elem_soln.size(), n_sf);

        std::vector<Real> node_weights(n_nodes);

        for (unsigned int n=0; n<n_nodes; n++)
          node_weights[n] = elem->node_ref(n).get_extra_integer(0);

        for (unsigned int n=0; n<n_nodes; n++)
          {
            std::vector<Real> weighted_shape(n_sf);
            Real weighted_sum = 0;

            for (unsigned int i=0; i<n_sf; i++)
              {
                weighted_shape[i] = node_weights[i] *
                  FEInterface::shape(Dim, fe_type, elem, i, refspace_nodes[n]);
                weighted_sum += weighted_shape[i];
              }

            // Zero before summation
            nodal_soln[n] = 0;

            // u_i = Sum (alpha_i w_i phi_i) / Sum (w_j phi_j)
            for (unsigned int i=0; i<n_sf; i++)
              nodal_soln[n] += elem_soln[i] * weighted_shape[i];
            nodal_soln[n] /= weighted_sum;
          }

        return;
      }
    }
} //  rational_nodal_soln()

} // anon namespace


namespace libMesh {

template <typename RealType>
void FEShim<0,RATIONAL_BERNSTEIN,RealType>::nodal_soln(const Elem * elem,
                                          const Order order,
                                          const std::vector<Number> & elem_soln,
                                          std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <typename RealType>
void FEShim<1,RATIONAL_BERNSTEIN,RealType>::nodal_soln(const Elem * elem,
                                            const Order order,
                                            const std::vector<Number> & elem_soln,
                                            std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <typename RealType>
void FEShim<2,RATIONAL_BERNSTEIN,RealType>::nodal_soln(const Elem * elem,
                                            const Order order,
                                            const std::vector<Number> & elem_soln,
                                            std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,RATIONAL_BERNSTEIN,RealType>::nodal_soln(const Elem * elem,
                                            const Order order,
                                            const std::vector<Number> & elem_soln,
                                            std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }


// Full specialization of n_dofs() function for every dimension
template <typename RealType> unsigned int FEShim<0,RATIONAL_BERNSTEIN,RealType>::n_dofs(const ElemType t, const Order o) { return FE<0,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,RATIONAL_BERNSTEIN,RealType>::n_dofs(const ElemType t, const Order o) { return FE<1,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,RATIONAL_BERNSTEIN,RealType>::n_dofs(const ElemType t, const Order o) { return FE<2,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,RATIONAL_BERNSTEIN,RealType>::n_dofs(const ElemType t, const Order o) { return FE<3,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
template <typename RealType> unsigned int FEShim<0,RATIONAL_BERNSTEIN,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<0,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<1,RATIONAL_BERNSTEIN,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<1,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<2,RATIONAL_BERNSTEIN,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<2,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<3,RATIONAL_BERNSTEIN,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<3,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <typename RealType> unsigned int FEShim<0,RATIONAL_BERNSTEIN,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<0,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<1,RATIONAL_BERNSTEIN,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<1,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<2,RATIONAL_BERNSTEIN,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<2,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<3,RATIONAL_BERNSTEIN,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<3,UnderlyingFEFamily<RATIONAL_BERNSTEIN>::value,RealType>::n_dofs_per_elem(t, o); }

// Our current rational FEMs are C^0 continuous
template <typename RealType> FEContinuity FEShim<0,RATIONAL_BERNSTEIN,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<1,RATIONAL_BERNSTEIN,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<2,RATIONAL_BERNSTEIN,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<3,RATIONAL_BERNSTEIN,RealType>::get_continuity() { return C_ZERO; }

// Rational FEMs are in general not hierarchic
template <typename RealType> bool FEShim<0,RATIONAL_BERNSTEIN,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<1,RATIONAL_BERNSTEIN,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<2,RATIONAL_BERNSTEIN,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<3,RATIONAL_BERNSTEIN,RealType>::is_hierarchic() { return false; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <typename RealType>
void FEShim<2,RATIONAL_BERNSTEIN,RealType>::compute_constraints (DofConstraints & constraints,
                                                      DofMap & dof_map,
                                                      const unsigned int variable_number,
                                                                 const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<RATIONAL_BERNSTEIN>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <typename RealType>
void FEShim<3,RATIONAL_BERNSTEIN,RealType>::compute_constraints (DofConstraints & constraints,
                                                      DofMap & dof_map,
                                                      const unsigned int variable_number,
                                                                 const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<RATIONAL_BERNSTEIN>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Bernstein shapes need reinit only for approximation orders >= 3,
// and at the moment we're only supporting equal-order shapes due to
// the need for nodal weights.
template <typename RealType> bool FEShim<0,RATIONAL_BERNSTEIN,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<1,RATIONAL_BERNSTEIN,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<2,RATIONAL_BERNSTEIN,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<3,RATIONAL_BERNSTEIN,RealType>::shapes_need_reinit() { return false; }

} // namespace libMesh

#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_RATIONAL_IMPL_H
