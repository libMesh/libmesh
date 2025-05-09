// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

// ------------------------------------------------------------
// Hierarchic-specific implementations

// Anonymous namespace for local helper functions
namespace {

void l2_hierarchic_nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln,
                              const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = order + add_p_level*elem->p_level();

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(order, L2_HIERARCHIC);

  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
        libmesh_assert_equal_to (elem_soln.size(), 1);

        std::fill(nodal_soln.begin(), nodal_soln.end(), elem_soln[0]);

        return;
      }


      // For other orders do interpolation at the nodes
      // explicitly.
    default:
      {
        const unsigned int n_sf =
          FEInterface::n_shape_functions(fe_type, elem);

        std::vector<Point> refspace_nodes;
        FEBase::get_refspace_nodes(elem_type,refspace_nodes);
        libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);
        libmesh_assert_equal_to (elem_soln.size(), n_sf);

        // Zero before summation
        std::fill(nodal_soln.begin(), nodal_soln.end(), 0);

        for (unsigned int n=0; n<n_nodes; n++)
          // u_i = Sum (alpha_i phi_i)
          for (unsigned int i=0; i<n_sf; i++)
            nodal_soln[n] += elem_soln[i] *
              FEInterface::shape(fe_type, elem, i, refspace_nodes[n]);

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
    case QUADSHELL9:
      return ((o+1)*(o+1));
    case HEX8:
    case HEX20:
    case HEX27:
      return ((o+1)*(o+1)*(o+1));
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      return ((o+1)*(o+1)*(o+2)/2);
    case TRI3:
    case TRI6:
    case TRI7:
      return ((o+1)*(o+2)/2);
    case TET4:
    case TET10:
    case TET14:
      return ((o+1)*(o+2)*(o+3)/6);
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for L2_HIERARCHIC FE family!");
    }
} // l2_hierarchic_n_dofs()



unsigned int l2_hierarchic_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return l2_hierarchic_n_dofs(e->type(), o);
}


} // anonymous namespace


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(L2_HIERARCHIC, l2_hierarchic_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(L2_HIERARCHIC)


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }

template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs(const Elem * e, const Order o) { return l2_hierarchic_n_dofs(e, o); }
template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs(const Elem * e, const Order o) { return l2_hierarchic_n_dofs(e, o); }
template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs(const Elem * e, const Order o) { return l2_hierarchic_n_dofs(e, o); }
template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs(const Elem * e, const Order o) { return l2_hierarchic_n_dofs(e, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
// Discontinuous L2 elements only have interior nodes
template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }

template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_hierarchic_n_dofs(&e, o); }
template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_hierarchic_n_dofs(&e, o); }
template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_hierarchic_n_dofs(&e, o); }
template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs_per_elem(const Elem & e, const Order o) { return l2_hierarchic_n_dofs(&e, o); }

// L2 Hierarchic FEMs are C^0 continuous
template <> FEContinuity FE<0,L2_HIERARCHIC>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,L2_HIERARCHIC>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,L2_HIERARCHIC>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,L2_HIERARCHIC>::get_continuity() const { return DISCONTINUOUS; }

// L2 Hierarchic FEMs are hierarchic (duh!)
template <> bool FE<0,L2_HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<1,L2_HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<2,L2_HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<3,L2_HIERARCHIC>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() is a NOOP for DISCONTINUOUS FE's
template <>
void FE<2,L2_HIERARCHIC>::compute_constraints (DofConstraints &,
                                               DofMap &,
                                               const unsigned int,
                                               const Elem *)
{ }

template <>
void FE<3,L2_HIERARCHIC>::compute_constraints (DofConstraints &,
                                               DofMap &,
                                               const unsigned int,
                                               const Elem *)
{ }
#endif // #ifdef LIBMESH_ENABLE_AMR

// L2-Hierarchic FEM shapes need reinit
template <> bool FE<0,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<1,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<2,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<3,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }

} // namespace libMesh
