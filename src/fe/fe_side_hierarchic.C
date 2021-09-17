// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// ------------------------------------------------------------
// Hierarchic-specific implementations

// Anonymous namespace for local helper functions
namespace {

void side_hierarchic_nodal_soln(const Elem * elem,
                                const Order /* order */,
                                const std::vector<Number> & /* elem_soln */,
                                std::vector<Number> & nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();

  nodal_soln.resize(n_nodes);

  libmesh_warning("Nodal solution requested for a side element; this makes no sense.");
} // side_hierarchic_nodal_soln()




unsigned int side_hierarchic_n_dofs_at_node(const ElemType t,
                                            const Order o,
                                            const unsigned int n)
{
  switch (t)
    {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      if (n < 2)
        return 1; // One per side
      else
        return 0;
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      if (n > 3 && n < 8)
        return o+1;
      else
        return 0;
    case HEX27:
      if (n > 19 && n < 26)
        return (o+1)*(o+1); // (o+1)^2 per side
      else
        return 0;
    case TRI6:
    case TRI7:
      if (n > 2 && n < 6)
        return o+1;
      else
        return 0;
    case INVALID_ELEM:
      return 0;
    // Without side nodes on all sides we can't support side elements
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SIDE_HIERARCHIC FE family!");
    }
} // side_hierarchic_n_dofs()



unsigned int side_hierarchic_n_dofs(const ElemType t, const Order o)
{
  switch (t)
    {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return 2; // One per side
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      return ((o+1)*4); // o+1 per side
    case HEX27:
      return ((o+1)*(o+1)*6); // (o+1)^2 per side
    case TRI6:
    case TRI7:
      return ((o+1)*3); // o+1 per side
    case INVALID_ELEM:
      return 0;
    // Without side nodes on all sides we can't support side elements
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for HIERARCHIC FE family!");
    }
} // side_hierarchic_n_dofs()


} // anonymous namespace




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <>
void FE<0,SIDE_HIERARCHIC>::nodal_soln(const Elem * elem,
                                       const Order order,
                                       const std::vector<Number> & elem_soln,
                                       std::vector<Number> & nodal_soln)
{ side_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<1,SIDE_HIERARCHIC>::nodal_soln(const Elem * elem,
                                       const Order order,
                                       const std::vector<Number> & elem_soln,
                                       std::vector<Number> & nodal_soln)
{ side_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<2,SIDE_HIERARCHIC>::nodal_soln(const Elem * elem,
                                       const Order order,
                                       const std::vector<Number> & elem_soln,
                                       std::vector<Number> & nodal_soln)
{ side_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<3,SIDE_HIERARCHIC>::nodal_soln(const Elem * elem,
                                       const Order order,
                                       const std::vector<Number> & elem_soln,
                                       std::vector<Number> & nodal_soln)
{ side_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln); }


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,SIDE_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return side_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<1,SIDE_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return side_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<2,SIDE_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return side_hierarchic_n_dofs(t, o); }
template <> unsigned int FE<3,SIDE_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return side_hierarchic_n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
template <> unsigned int FE<0,SIDE_HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return side_hierarchic_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<1,SIDE_HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return side_hierarchic_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<2,SIDE_HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return side_hierarchic_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,SIDE_HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return side_hierarchic_n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,SIDE_HIERARCHIC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<1,SIDE_HIERARCHIC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<2,SIDE_HIERARCHIC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<3,SIDE_HIERARCHIC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Side FEMs are discontinuous from side to side
template <> FEContinuity FE<0,SIDE_HIERARCHIC>::get_continuity() const { return SIDE_DISCONTINUOUS; }
template <> FEContinuity FE<1,SIDE_HIERARCHIC>::get_continuity() const { return SIDE_DISCONTINUOUS; }
template <> FEContinuity FE<2,SIDE_HIERARCHIC>::get_continuity() const { return SIDE_DISCONTINUOUS; }
template <> FEContinuity FE<3,SIDE_HIERARCHIC>::get_continuity() const { return SIDE_DISCONTINUOUS; }

// Side Hierarchic FEMs are hierarchic (duh!)
template <> bool FE<0,SIDE_HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<1,SIDE_HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<2,SIDE_HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<3,SIDE_HIERARCHIC>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <>
void FE<2,SIDE_HIERARCHIC>::compute_constraints (DofConstraints & constraints,
                                                 DofMap & dof_map,
                                                 const unsigned int variable_number,
                                                 const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <>
void FE<3,SIDE_HIERARCHIC>::compute_constraints (DofConstraints & constraints,
                                                 DofMap & dof_map,
                                                 const unsigned int variable_number,
                                                 const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Hierarchic FEM shapes need reinit
template <> bool FE<0,SIDE_HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<1,SIDE_HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<2,SIDE_HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<3,SIDE_HIERARCHIC>::shapes_need_reinit() const { return true; }

} // namespace libMesh
