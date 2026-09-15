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
#include "libmesh/fe_macro.h"
#include "libmesh/quadrature_gauss_lobatto.h"


namespace libMesh
{

/**
 * \returns The number of degrees of freedom that the L2_LAGRANGE_GLL basis of order \p o
 * carries on an element of type \p t.
 *
 * The basis is a tensor product of one-dimensional bases, so it lives on the tensor
 * product element types alone, and it carries \f$(p+1)^d\f$ degrees of freedom on one of
 * those whatever its node count: every degree of freedom belongs to the element, so an
 * order of eight asks nothing of the element beyond its shape.
 */
unsigned int l2_lagrange_gll_n_dofs(const ElemType t, const Order o)
{
  libmesh_assert_greater (o, 0);

  const unsigned int n = static_cast<unsigned int>(o) + 1;

  libmesh_error_msg_if(n > QGaussLobatto::max_points_1D,
                       "The L2_LAGRANGE_GLL basis of order " << o << " would interpolate at "
                       << n << " Gauss-Lobatto points, and they are tabulated up to "
                       << QGaussLobatto::max_points_1D << ".");

  switch (t)
    {
    case NODEELEM:
      return 1;

    case EDGE2:
    case EDGE3:
    case EDGE4:
      return n;

    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      return n*n;

    case HEX8:
    case HEX20:
    case HEX27:
      return n*n*n;

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("The L2_LAGRANGE_GLL basis is a tensor product, so it has no form on "
                        << Utility::enum_to_string(t) << ".");
    }
}


unsigned int l2_lagrange_gll_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return l2_lagrange_gll_n_dofs(e->type(), o);
}

/**
 * Evaluates a solution at the nodes of \p elem, for the writers that report a variable
 * there.
 *
 * The interpolation points of this basis are the Gauss-Lobatto points, which are the nodes
 * of the element only at orders one and two, so beyond those the values are interpolated
 * at the nodes rather than read off the coefficients.
 */
void l2_lagrange_gll_nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln,
                                const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  nodal_soln.resize(n_nodes);

  const FEType fe_type(order, L2_LAGRANGE_GLL);

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

LIBMESH_FE_NODAL_SOLN(L2_LAGRANGE_GLL, l2_lagrange_gll_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(L2_LAGRANGE_GLL)


// Every degree of freedom of this basis belongs to the element
template <> unsigned int FE<0,L2_LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<1,L2_LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<2,L2_LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<3,L2_LAGRANGE_GLL>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }

template <> unsigned int FE<0,L2_LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return l2_lagrange_gll_n_dofs(e, o); }
template <> unsigned int FE<1,L2_LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return l2_lagrange_gll_n_dofs(e, o); }
template <> unsigned int FE<2,L2_LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return l2_lagrange_gll_n_dofs(e, o); }
template <> unsigned int FE<3,L2_LAGRANGE_GLL>::n_dofs(const Elem * e, const Order o) { return l2_lagrange_gll_n_dofs(e, o); }

template <> unsigned int FE<0,L2_LAGRANGE_GLL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_LAGRANGE_GLL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_LAGRANGE_GLL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_LAGRANGE_GLL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,L2_LAGRANGE_GLL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_LAGRANGE_GLL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_LAGRANGE_GLL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_LAGRANGE_GLL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,L2_LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<1,L2_LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<2,L2_LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }
template <> unsigned int FE<3,L2_LAGRANGE_GLL>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_gll_n_dofs(t, o); }

template <> unsigned int FE<0,L2_LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_lagrange_gll_n_dofs(&e, o); }
template <> unsigned int FE<1,L2_LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_lagrange_gll_n_dofs(&e, o); }
template <> unsigned int FE<2,L2_LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_lagrange_gll_n_dofs(&e, o); }
template <> unsigned int FE<3,L2_LAGRANGE_GLL>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_lagrange_gll_n_dofs(&e, o); }


template <> FEContinuity FE<0,L2_LAGRANGE_GLL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,L2_LAGRANGE_GLL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,L2_LAGRANGE_GLL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,L2_LAGRANGE_GLL>::get_continuity() const { return DISCONTINUOUS; }

// A nodal basis is not hierarchic
template <> bool FE<0,L2_LAGRANGE_GLL>::is_hierarchic() const { return false; }
template <> bool FE<1,L2_LAGRANGE_GLL>::is_hierarchic() const { return false; }
template <> bool FE<2,L2_LAGRANGE_GLL>::is_hierarchic() const { return false; }
template <> bool FE<3,L2_LAGRANGE_GLL>::is_hierarchic() const { return false; }

// The shape functions are the same on every element of a type
template <> bool FE<0,L2_LAGRANGE_GLL>::shapes_need_reinit() const { return false; }
template <> bool FE<1,L2_LAGRANGE_GLL>::shapes_need_reinit() const { return false; }
template <> bool FE<2,L2_LAGRANGE_GLL>::shapes_need_reinit() const { return false; }
template <> bool FE<3,L2_LAGRANGE_GLL>::shapes_need_reinit() const { return false; }

// A DISCONTINUOUS basis shares nothing across a hanging node, so it constrains nothing
#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,L2_LAGRANGE_GLL>::compute_constraints (DofConstraints &,
                                                 DofMap &,
                                                 const unsigned int,
                                                 const Elem *)
{ }

template <>
void FE<3,L2_LAGRANGE_GLL>::compute_constraints (DofConstraints &,
                                                 DofMap &,
                                                 const unsigned int,
                                                 const Elem *)
{ }
#endif // LIBMESH_ENABLE_AMR


} // namespace libMesh
