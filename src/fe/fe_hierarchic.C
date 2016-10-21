// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/string_to_enum.h"

namespace libMesh
{

// ------------------------------------------------------------
// Hierarchic-specific implementations

// Anonymous namespace for local helper functions
namespace {

void hierarchic_nodal_soln(const Elem * elem,
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
  FEType fe_type(totalorder, HIERARCHIC);

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
          // FE<Dim,T>::n_shape_functions(elem_type, totalorder);
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
                // FE<Dim,T>::shape(elem, order, i, mapped_point);
                FEInterface::shape(Dim, fe_type, elem, i, refspace_nodes[n]);
          }

        return;
      }
    }
} // hierarchic_nodal_soln()




unsigned int hierarchic_n_dofs(const ElemType t, const Order o)
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
      libmesh_assert_less (o, 2);
    case QUAD8:
    case QUAD9:
      return ((o+1)*(o+1));
    case HEX8:
      libmesh_assert_less (o, 2);
    case HEX20:
      libmesh_assert_less (o, 2);
    case HEX27:
      return ((o+1)*(o+1)*(o+1));
    case TRI3:
      libmesh_assert_less (o, 2);
    case TRI6:
      return ((o+1)*(o+2)/2);
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for HIERARCHIC FE family!");
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
} // hierarchic_n_dofs()




unsigned int hierarchic_n_dofs_at_node(const ElemType t,
                                       const Order o,
                                       const unsigned int n)
{
  libmesh_assert_greater (o, 0);
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
    case EDGE3:
      switch (n)
        {
        case 0:
        case 1:
          return 1;
          // Internal DoFs are associated with the elem, not its nodes
        case 2:
          return 0;
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
        }
    case TRI3:
      libmesh_assert_less (n, 3);
      libmesh_assert_less (o, 2);
    case TRI6:
      switch (n)
        {
        case 0:
        case 1:
        case 2:
          return 1;

        case 3:
        case 4:
        case 5:
          return (o-1);

          // Internal DoFs are associated with the elem, not its nodes
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
        }
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (n, 4);
      libmesh_assert_less (o, 2);
    case QUAD8:
    case QUAD9:
      switch (n)
        {
        case 0:
        case 1:
        case 2:
        case 3:
          return 1;

        case 4:
        case 5:
        case 6:
        case 7:
          return (o-1);

          // Internal DoFs are associated with the elem, not its nodes
        case 8:
          return 0;

        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD4/8/9!");
        }
    case HEX8:
      libmesh_assert_less (n, 8);
      libmesh_assert_less (o, 2);
    case HEX20:
      libmesh_assert_less (n, 20);
      libmesh_assert_less (o, 2);
    case HEX27:
      switch (n)
        {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
          return 1;

        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
        case 18:
        case 19:
          return (o-1);

        case 20:
        case 21:
        case 22:
        case 23:
        case 24:
        case 25:
          return ((o-1)*(o-1));

          // Internal DoFs are associated with the elem, not its nodes
        case 26:
          return 0;
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for HEX8/20/27!");
        }

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
} // hierarchic_n_dofs_at_node()




unsigned int hierarchic_n_dofs_per_elem(const ElemType t,
                                        const Order o)
{
  libmesh_assert_greater (o, 0);
  switch (t)
    {
    case NODEELEM:
      return 0;
    case EDGE2:
    case EDGE3:
      return (o-1);
    case TRI3:
    case TRISHELL3:
    case QUAD4:
    case QUADSHELL4:
      return 0;
    case TRI6:
      return ((o-1)*(o-2)/2);
    case QUAD8:
    case QUAD9:
      return ((o-1)*(o-1));
    case HEX8:
    case HEX20:
      libmesh_assert_less (o, 2);
      return 0;
    case HEX27:
      return ((o-1)*(o-1)*(o-1));
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
} // hierarchic_n_dofs_per_elem()

} // anonymous namespace




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <>
void FE<0,HIERARCHIC>::nodal_soln(const Elem * elem,
                                  const Order order,
                                  const std::vector<Number> & elem_soln,
                                  std::vector<Number> & nodal_soln)
{ hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <>
void FE<1,HIERARCHIC>::nodal_soln(const Elem * elem,
                                  const Order order,
                                  const std::vector<Number> & elem_soln,
                                  std::vector<Number> & nodal_soln)
{ hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <>
void FE<2,HIERARCHIC>::nodal_soln(const Elem * elem,
                                  const Order order,
                                  const std::vector<Number> & elem_soln,
                                  std::vector<Number> & nodal_soln)
{ hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <>
void FE<3,HIERARCHIC>::nodal_soln(const Elem * elem,
                                  const Order order,
                                  const std::vector<Number> & elem_soln,
                                  std::vector<Number> & nodal_soln)
{ hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return hierarchic_n_dofs(t, o); }
template <> unsigned int FE<1,HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return hierarchic_n_dofs(t, o); }
template <> unsigned int FE<2,HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return hierarchic_n_dofs(t, o); }
template <> unsigned int FE<3,HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return hierarchic_n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
template <> unsigned int FE<0,HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hierarchic_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<1,HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hierarchic_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<2,HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hierarchic_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,HIERARCHIC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hierarchic_n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return hierarchic_n_dofs_per_elem(t, o); }
template <> unsigned int FE<1,HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return hierarchic_n_dofs_per_elem(t, o); }
template <> unsigned int FE<2,HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return hierarchic_n_dofs_per_elem(t, o); }
template <> unsigned int FE<3,HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return hierarchic_n_dofs_per_elem(t, o); }

// Hierarchic FEMs are C^0 continuous
template <> FEContinuity FE<0,HIERARCHIC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<1,HIERARCHIC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<2,HIERARCHIC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<3,HIERARCHIC>::get_continuity() const { return C_ZERO; }

// Hierarchic FEMs are hierarchic (duh!)
template <> bool FE<0,HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<1,HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<2,HIERARCHIC>::is_hierarchic() const { return true; }
template <> bool FE<3,HIERARCHIC>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <>
void FE<2,HIERARCHIC>::compute_constraints (DofConstraints & constraints,
                                            DofMap & dof_map,
                                            const unsigned int variable_number,
                                            const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <>
void FE<3,HIERARCHIC>::compute_constraints (DofConstraints & constraints,
                                            DofMap & dof_map,
                                            const unsigned int variable_number,
                                            const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Hierarchic FEM shapes need reinit
template <> bool FE<0,HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<1,HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<2,HIERARCHIC>::shapes_need_reinit() const { return true; }
template <> bool FE<3,HIERARCHIC>::shapes_need_reinit() const { return true; }

} // namespace libMesh
