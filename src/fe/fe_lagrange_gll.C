// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_lagrange_gll_shape_1D.h"
#include "libmesh/fe_macro.h"


namespace libMesh
{

/**
 * \returns The number of degrees of freedom that the LAGRANGE_GLL basis of order \p o carries
 * on an element of type \p t.
 *
 * The basis is a tensor product, so it applies only to tensor product element types. On
 * these types it has \f$(p+1)^d\f$ degrees of freedom.
 *
 * Each entity that owns degrees of freedom needs a node, and this limits the order on some
 * types. A vertex owns one degree of freedom and an edge owns p-1. Thus a type without mid-edge
 * nodes supports only order one. A face owns \f$(p-1)^2\f$, so in three dimensions a type
 * without face nodes supports only order one as well. The element owns the interior degrees of
 * freedom, so no type needs an interior node. HIERARCHIC has the same limits for the same
 * reason.
 */
unsigned int lagrange_gll_n_dofs(const ElemType t, const Order o)
{
  libmesh_assert_greater (o, 0);

  const unsigned int n = static_cast<unsigned int>(o) + 1;

  libmesh_error_msg_if(n > fe_lagrange_gll_max_points_1D,
                       "The LAGRANGE_GLL basis of order " << o << " needs " << n
                       << " Gauss-Lobatto points. The table has at most "
                       << fe_lagrange_gll_max_points_1D << ".");

  switch (t)
    {
    case NODEELEM:
      return 1;

      // An edge element has only vertices and an interior, so it needs no other nodes at any
      // order
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return n;

      // Above order one, the edge degrees of freedom need mid-edge nodes
    case QUAD4:
    case QUADSHELL4:
      libmesh_error_msg_if(o > FIRST,
                           "The LAGRANGE_GLL basis of order " << o << " puts " << o - 1
                           << " degrees of freedom on each edge. "
                           << Utility::enum_to_string(t) << " has no mid-edge nodes. "
                           "Use QUAD8 or QUAD9.");
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      return n*n;

      // In three dimensions, the face degrees of freedom also need face nodes
    case HEX8:
      libmesh_error_msg_if(o > FIRST,
                           "The LAGRANGE_GLL basis of order " << o
                           << " puts degrees of freedom on the edges and faces of a hexahedron. "
                           "HEX8 has no edge or face nodes. Use HEX27.");
      libmesh_fallthrough();
    case HEX20:
      libmesh_error_msg_if(o > FIRST,
                           "The LAGRANGE_GLL basis of order " << o << " puts " << (o - 1)*(o - 1)
                           << " degrees of freedom on each face. HEX20 has no face nodes. "
                           "Use HEX27.");
      libmesh_fallthrough();
    case HEX27:
      return n*n*n;

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("The LAGRANGE_GLL basis is a tensor product. It does not support "
                        << Utility::enum_to_string(t) << ".");
    }
}


unsigned int lagrange_gll_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return lagrange_gll_n_dofs(e->type(), o);
}

/**
 * \returns The number of degrees of freedom that node \p n of an element of type \p t owns.
 *
 * A vertex owns the degree of freedom at its own interpolation point. A mid-edge node owns
 * the p-1 points inside its edge. A mid-face node owns the \f$(p-1)^2\f$ points inside its face.
 * An interior node owns none. The element owns the interior points. This lets static
 * condensation eliminate them, and it lets types without an interior node use the basis.
 */
unsigned int lagrange_gll_n_dofs_at_node(const ElemType t,
                                         const Order o,
                                         const unsigned int n)
{
  libmesh_assert_greater (o, 0);

  const unsigned int per_edge = static_cast<unsigned int>(o) - 1;

  switch (t)
    {
    case NODEELEM:
      return 1;

    case EDGE2:
    case EDGE3:
    case EDGE4:
      return (n < 2) ? 1 : 0;

    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (n, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      if (n < 4)
        return 1;
      if (n < 8)
        return per_edge;
      return 0;

    case HEX8:
      libmesh_assert_less (n, 8);
      libmesh_fallthrough();
    case HEX20:
    case HEX27:
      if (n < 8)
        return 1;
      if (n < 20)
        return per_edge;
      if (n < 26)
        return per_edge * per_edge;
      return 0;

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("The LAGRANGE_GLL basis is a tensor product. It does not support "
                        << Utility::enum_to_string(t) << ".");
    }
}

/**
 * \returns The number of degrees of freedom that the element itself owns. These are the
 * \f$(p-1)^d\f$ interpolation points in its interior.
 */
unsigned int lagrange_gll_n_dofs_per_elem(const ElemType t, const Order o)
{
  libmesh_assert_greater (o, 0);

  const unsigned int per_edge = static_cast<unsigned int>(o) - 1;

  switch (t)
    {
    case NODEELEM:
    case INVALID_ELEM:
      return 0;

    case EDGE2:
    case EDGE3:
    case EDGE4:
      return per_edge;

    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      return per_edge * per_edge;

    case HEX8:
    case HEX20:
    case HEX27:
      return per_edge * per_edge * per_edge;

    default:
      libmesh_error_msg("The LAGRANGE_GLL basis is a tensor product. It does not support "
                        << Utility::enum_to_string(t) << ".");
    }
}


/**
 * Evaluates a solution at the nodes of \p elem. Output writers use these nodal values.
 *
 * Each vertex coefficient is the value at its node. Above order two, the other interpolation
 * points are not nodes. Thus this function evaluates the basis at each node.
 */
void lagrange_gll_nodal_soln(const Elem * elem,
                             const Order order,
                             const std::vector<Number> & elem_soln,
                             std::vector<Number> & nodal_soln,
                             const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  nodal_soln.resize(n_nodes);

  const FEType fe_type(order, LAGRANGE_GLL);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(fe_type, elem, add_p_level);

  std::vector<Point> refspace_nodes;
  FEBase::get_refspace_nodes(elem->type(), refspace_nodes);
  libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);
  libmesh_assert_equal_to (elem_soln.size(), n_sf);

  std::fill(nodal_soln.begin(), nodal_soln.end(), 0);

  for (const auto n : make_range(n_nodes))
    for (const auto i : make_range(n_sf))
      nodal_soln[n] += elem_soln[i] *
        FEInterface::shape(fe_type, elem, i, refspace_nodes[n], add_p_level);
}

LIBMESH_FE_NODAL_SOLN(LAGRANGE_GLL, lagrange_gll_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(LAGRANGE_GLL)


template <> unsigned int FE<0,LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<1,LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<2,LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<3,LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return lagrange_gll_n_dofs(t, o); }

template <> unsigned int FE<0,LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return lagrange_gll_n_dofs(e, o); }
template <> unsigned int FE<1,LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return lagrange_gll_n_dofs(e, o); }
template <> unsigned int FE<2,LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return lagrange_gll_n_dofs(e, o); }
template <> unsigned int FE<3,LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return lagrange_gll_n_dofs(e, o); }

template <> unsigned int FE<0,LAGRANGE_GLL>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<1,LAGRANGE_GLL>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<2,LAGRANGE_GLL>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,LAGRANGE_GLL>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(t, o, n); }

template <> unsigned int FE<0,LAGRANGE_GLL>::n_dofs_at_node(const Elem & e, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(e.type(), o, n); }
template <> unsigned int FE<1,LAGRANGE_GLL>::n_dofs_at_node(const Elem & e, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(e.type(), o, n); }
template <> unsigned int FE<2,LAGRANGE_GLL>::n_dofs_at_node(const Elem & e, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(e.type(), o, n); }
template <> unsigned int FE<3,LAGRANGE_GLL>::n_dofs_at_node(const Elem & e, const Order o, const unsigned int n) { return lagrange_gll_n_dofs_at_node(e.type(), o, n); }

template <> unsigned int FE<0,LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return lagrange_gll_n_dofs_per_elem(t, o); }
template <> unsigned int FE<1,LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return lagrange_gll_n_dofs_per_elem(t, o); }
template <> unsigned int FE<2,LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return lagrange_gll_n_dofs_per_elem(t, o); }
template <> unsigned int FE<3,LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return lagrange_gll_n_dofs_per_elem(t, o); }

template <> unsigned int FE<0,LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return lagrange_gll_n_dofs_per_elem(e.type(), o); }
template <> unsigned int FE<1,LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return lagrange_gll_n_dofs_per_elem(e.type(), o); }
template <> unsigned int FE<2,LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return lagrange_gll_n_dofs_per_elem(e.type(), o); }
template <> unsigned int FE<3,LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return lagrange_gll_n_dofs_per_elem(e.type(), o); }

// The basis is continuous across a shared entity, because the two elements put each
// shared degree of freedom at the same point
template <> FEContinuity FE<0,LAGRANGE_GLL>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<1,LAGRANGE_GLL>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<2,LAGRANGE_GLL>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<3,LAGRANGE_GLL>::get_continuity() const { return C_ZERO; }

// A nodal basis is not hierarchic
template <> bool FE<0,LAGRANGE_GLL>::is_hierarchic() const { return false; }
template <> bool FE<1,LAGRANGE_GLL>::is_hierarchic() const { return false; }
template <> bool FE<2,LAGRANGE_GLL>::is_hierarchic() const { return false; }
template <> bool FE<3,LAGRANGE_GLL>::is_hierarchic() const { return false; }

// Each shape function depends on the orientation of the entity that owns it. Elements
// with different orientations thus need new shape values.
template <> bool FE<0,LAGRANGE_GLL>::shapes_need_reinit() const { return true; }
template <> bool FE<1,LAGRANGE_GLL>::shapes_need_reinit() const { return true; }
template <> bool FE<2,LAGRANGE_GLL>::shapes_need_reinit() const { return true; }
template <> bool FE<3,LAGRANGE_GLL>::shapes_need_reinit() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,LAGRANGE_GLL>::compute_constraints (DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                          const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <>
void FE<3,LAGRANGE_GLL>::compute_constraints (DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                          const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // LIBMESH_ENABLE_AMR


} // namespace libMesh
