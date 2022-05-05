// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

// Anonymous namespace for local helper functions
namespace {

using namespace libMesh;

static const FEFamily _underlying_fe_family = BERNSTEIN;

void rational_nodal_soln(const Elem * elem,
                         const Order order,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> & nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(order, _underlying_fe_family);
  FEType p_refined_fe_type(totalorder, _underlying_fe_family);

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
          FEInterface::n_shape_functions(fe_type, elem);

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
                  FEInterface::shape(fe_type, elem, i, refspace_nodes[n]);
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

template <>
void FE<0,RATIONAL_BERNSTEIN>::nodal_soln(const Elem * elem,
                                          const Order order,
                                          const std::vector<Number> & elem_soln,
                                          std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<1,RATIONAL_BERNSTEIN>::nodal_soln(const Elem * elem,
                                            const Order order,
                                            const std::vector<Number> & elem_soln,
                                            std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<2,RATIONAL_BERNSTEIN>::nodal_soln(const Elem * elem,
                                            const Order order,
                                            const std::vector<Number> & elem_soln,
                                            std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<3,RATIONAL_BERNSTEIN>::nodal_soln(const Elem * elem,
                                            const Order order,
                                            const std::vector<Number> & elem_soln,
                                            std::vector<Number> & nodal_soln)
{ rational_nodal_soln(elem, order, elem_soln, nodal_soln); }


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,RATIONAL_BERNSTEIN>::n_dofs(const ElemType t, const Order o) { return FE<0,_underlying_fe_family>::n_dofs(t, o); }
template <> unsigned int FE<1,RATIONAL_BERNSTEIN>::n_dofs(const ElemType t, const Order o) { return FE<1,_underlying_fe_family>::n_dofs(t, o); }
template <> unsigned int FE<2,RATIONAL_BERNSTEIN>::n_dofs(const ElemType t, const Order o) { return FE<2,_underlying_fe_family>::n_dofs(t, o); }
template <> unsigned int FE<3,RATIONAL_BERNSTEIN>::n_dofs(const ElemType t, const Order o) { return FE<3,_underlying_fe_family>::n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
template <> unsigned int FE<0,RATIONAL_BERNSTEIN>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<0,_underlying_fe_family>::n_dofs_at_node(t, o, n); }
template <> unsigned int FE<1,RATIONAL_BERNSTEIN>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<1,_underlying_fe_family>::n_dofs_at_node(t, o, n); }
template <> unsigned int FE<2,RATIONAL_BERNSTEIN>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<2,_underlying_fe_family>::n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,RATIONAL_BERNSTEIN>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<3,_underlying_fe_family>::n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,RATIONAL_BERNSTEIN>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<0,_underlying_fe_family>::n_dofs_per_elem(t, o); }
template <> unsigned int FE<1,RATIONAL_BERNSTEIN>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<1,_underlying_fe_family>::n_dofs_per_elem(t, o); }
template <> unsigned int FE<2,RATIONAL_BERNSTEIN>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<2,_underlying_fe_family>::n_dofs_per_elem(t, o); }
template <> unsigned int FE<3,RATIONAL_BERNSTEIN>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<3,_underlying_fe_family>::n_dofs_per_elem(t, o); }

// Our current rational FEMs are C^0 continuous
template <> FEContinuity FE<0,RATIONAL_BERNSTEIN>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<1,RATIONAL_BERNSTEIN>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<2,RATIONAL_BERNSTEIN>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<3,RATIONAL_BERNSTEIN>::get_continuity() const { return C_ZERO; }

// Rational FEMs are in general not hierarchic
template <> bool FE<0,RATIONAL_BERNSTEIN>::is_hierarchic() const { return false; }
template <> bool FE<1,RATIONAL_BERNSTEIN>::is_hierarchic() const { return false; }
template <> bool FE<2,RATIONAL_BERNSTEIN>::is_hierarchic() const { return false; }
template <> bool FE<3,RATIONAL_BERNSTEIN>::is_hierarchic() const { return false; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <>
void FE<2,RATIONAL_BERNSTEIN>::compute_constraints (DofConstraints & constraints,
                                                      DofMap & dof_map,
                                                      const unsigned int variable_number,
                                                      const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <>
void FE<3,RATIONAL_BERNSTEIN>::compute_constraints (DofConstraints & constraints,
                                                      DofMap & dof_map,
                                                      const unsigned int variable_number,
                                                      const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Bernstein shapes need reinit only for approximation orders >= 3,
// but for *RATIONAL* shapes we could need reinit anywhere because our
// nodal weights can change from element to element.
template <> bool FE<0,RATIONAL_BERNSTEIN>::shapes_need_reinit() const { return true; }
template <> bool FE<1,RATIONAL_BERNSTEIN>::shapes_need_reinit() const { return true; }
template <> bool FE<2,RATIONAL_BERNSTEIN>::shapes_need_reinit() const { return true; }
template <> bool FE<3,RATIONAL_BERNSTEIN>::shapes_need_reinit() const { return true; }

} // namespace libMesh

#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
