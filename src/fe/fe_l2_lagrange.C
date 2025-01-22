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
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"

namespace libMesh
{

// ------------------------------------------------------------
// Lagrange-specific implementations


// Anonymous namespace for local helper functions
namespace {
unsigned int l2_lagrange_n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {

      // linear Lagrange shape functions
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 2;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
            return 3;

          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 4;

          case TET4:
          case TET10:
          case TET14:
            return 4;

          case HEX8:
          case HEX20:
          case HEX27:
            return 8;

          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            return 6;

          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
            return 5;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 3;

          case TRI3:
          case TRI6:
          case TRI7:
            return 6;

          case QUAD8:
          case QUADSHELL8:
            return 8;

          case QUAD4:
          case QUAD9:
          case QUADSHELL9:
            return 9;

          case TET4:
          case TET10:
          case TET14:
            return 10;

          case HEX20:
            return 20;

          case HEX8:
          case HEX27:
            return 27;

          case PRISM15:
            return 15;

          case PRISM6:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            return 18;

          case PYRAMID13:
            return 13;

          case PYRAMID5:
          case PYRAMID14:
          case PYRAMID18:
            return 14;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 4;

          case TRI7:
            return 7;

          case TET14:
            return 14;

          case PRISM20:
            return 20;

          case PRISM21:
            return 21;

          case PYRAMID18:
            return 18;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for L2_LAGRANGE FE family!");
    }
}

} // anonymous namespace


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(L2_LAGRANGE, lagrange_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(L2_LAGRANGE)


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }


// L2 Lagrange elements have all dofs on elements (hence none on nodes)
template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }


// L2 Lagrange elements have all dofs on elements
template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }

// L2 Lagrange FEMs are DISCONTINUOUS
template <> FEContinuity FE<0,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }

// Lagrange FEMs are not hierarchic
template <> bool FE<0,L2_LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<1,L2_LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<2,L2_LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<3,L2_LAGRANGE>::is_hierarchic() const { return false; }

// Lagrange FEM shapes do not need reinit (is this always true?)
template <> bool FE<0,L2_LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<1,L2_LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<2,L2_LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<3,L2_LAGRANGE>::shapes_need_reinit() const { return false; }

// We don't need any constraints for this DISCONTINUOUS basis, hence these
// are NOOPs
#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,L2_LAGRANGE>::compute_constraints (DofConstraints &,
                                             DofMap &,
                                             const unsigned int,
                                             const Elem *)
{ }

template <>
void FE<3,L2_LAGRANGE>::compute_constraints (DofConstraints &,
                                             DofMap &,
                                             const unsigned int,
                                             const Elem *)
{ }
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
