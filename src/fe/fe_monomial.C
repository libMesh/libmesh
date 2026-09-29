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
#include "libmesh/fe.h"
#include "libmesh/fe_monomial_n_dofs.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"

namespace libMesh
{

unsigned int monomial_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return monomial_n_dofs(e->type(), o);
}


unsigned int monomial_n_dofs(const ElemType t, const Order o)
{
  // The counts live in a header so that device code, which cannot call in
  // here, reads the same ones; reporting the error stays here.
  const unsigned int n = FECounts::monomial_n_dofs(t, o);

  libmesh_error_msg_if(n == invalid_uint,
                       "ERROR: Bad ElemType = " << Utility::enum_to_string(t) <<
                       " for " << Utility::enum_to_string(o) << " order approximation!");

  return n;
} // monomial_n_dofs()


// Anonymous namespace for local helper functions
namespace {

void monomial_nodal_soln(const Elem * elem,
                         const Order order,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> & nodal_soln,
                         const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = order + add_p_level*elem->p_level();

  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
        libmesh_assert_equal_to (elem_soln.size(), 1);

        std::fill(nodal_soln.begin(), nodal_soln.end(), elem_soln[0]);

        return;
      }


      // For other orders, do interpolation at the nodes
      // explicitly.
    default:
      {
        // FEType object to be passed to various FEInterface functions below.
        FEType fe_type(order, MONOMIAL);

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
      } // default
    } // switch
} // monomial_nodal_soln()




} // anonymous namespace


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(MONOMIAL, monomial_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(MONOMIAL)


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

template <> unsigned int FE<0,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
// Monomials have no dofs at nodes, only element dofs.
template <> unsigned int FE<0,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

template <> unsigned int FE<0,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }


// Full specialization of get_continuity() function for every dimension.
template <> FEContinuity FE<0,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }

// Full specialization of is_hierarchic() function for every dimension.
// The monomials are hierarchic!
template <> bool FE<0,MONOMIAL>::is_hierarchic() const { return true; }
template <> bool FE<1,MONOMIAL>::is_hierarchic() const { return true; }
template <> bool FE<2,MONOMIAL>::is_hierarchic() const { return true; }
template <> bool FE<3,MONOMIAL>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR

// Full specialization of compute_constraints() function for 2D and
// 3D only.  There are no constraints for discontinuous elements, so
// we do nothing.
template <> void FE<2,MONOMIAL>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem *) {}
template <> void FE<3,MONOMIAL>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem *) {}

#endif // #ifdef LIBMESH_ENABLE_AMR

// Full specialization of shapes_need_reinit() function for every dimension.
template <> bool FE<0,MONOMIAL>::shapes_need_reinit() const { return false; }
template <> bool FE<1,MONOMIAL>::shapes_need_reinit() const { return false; }
template <> bool FE<2,MONOMIAL>::shapes_need_reinit() const { return false; }
template <> bool FE<3,MONOMIAL>::shapes_need_reinit() const { return false; }

} // namespace libMesh
