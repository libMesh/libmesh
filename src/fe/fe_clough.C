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
// Clough-specific implementations

// Anonymous namespace for local helper functions
namespace {

void clough_nodal_soln(const Elem * elem,
                       const Order order,
                       const std::vector<Number> & elem_soln,
                       std::vector<Number> &       nodal_soln,
                       const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = order + add_p_level*elem->p_level();

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(order, CLOUGH);

  switch (totalorder)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      // Piecewise cubic shape functions
    case THIRD:
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

    default:
      libmesh_error_msg("ERROR: Invalid total order " << totalorder);
    }
} // clough_nodal_soln()



unsigned int CLOUGH_n_dofs(const ElemType t, const Order o)
{
  if (t == INVALID_ELEM)
    return 0;

  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      {
        switch (t)
          {
          case TRI6:
          case TRI7:
            return 9;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
      // Piecewise cubic Clough-Tocher element
    case THIRD:
      {
        switch (t)
          {
          case TRI6:
          case TRI7:
            return 12;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for CLOUGH FE family!");
    }
} // CLOUGH_n_dofs()



unsigned int CLOUGH_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return CLOUGH_n_dofs(e->type(), o);
}



unsigned int CLOUGH_n_dofs_at_node(const ElemType t,
                                   const Order o,
                                   const unsigned int n)
{
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      {
        switch (t)
          {
            // The 2D Clough-Tocher defined on a 6-noded triangle
          case TRI6:
          case TRI7:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 3;

                case 3:
                case 4:
                case 5:
                case 6:
                  return 0;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI!");
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
      // The third-order clough shape functions
    case THIRD:
      {
        switch (t)
          {
            // The 2D Clough-Tocher defined on a 6-noded triangle
          case TRI6:
          case TRI7:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 3;

                case 3:
                case 4:
                case 5:
                  return 1;

                case 6:
                  return 0;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for CLOUGH FE family!");
    }
} // CLOUGH_n_dofs_at_node()



unsigned int CLOUGH_n_dofs_at_node(const Elem & e,
                                   const Order o,
                                   const unsigned int n)
{
  return CLOUGH_n_dofs_at_node(e.type(), o, n);
}



unsigned int CLOUGH_n_dofs_per_elem(const ElemType t, const Order o)
{
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      // The third-order Clough-Tocher shape functions
    case THIRD:
      {
        switch (t)
          {
            // The 2D clough defined on a 6-noded triangle
          case TRI6:
          case TRI7:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for CLOUGH FE family!");
    }
} // CLOUGH_n_dofs_per_elem()



unsigned int CLOUGH_n_dofs_per_elem(const Elem & e, const Order o)
{
  return CLOUGH_n_dofs_per_elem(e.type(), o);
}

} // anonymous




// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(CLOUGH, clough_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(CLOUGH)


// Instantiate n_dofs*() functions for every dimension
LIBMESH_DEFAULT_NDOFS(CLOUGH)


// Clough FEMs are C^1 continuous
template <> FEContinuity FE<0,CLOUGH>::get_continuity() const { return C_ONE; }
template <> FEContinuity FE<1,CLOUGH>::get_continuity() const { return C_ONE; }
template <> FEContinuity FE<2,CLOUGH>::get_continuity() const { return C_ONE; }
template <> FEContinuity FE<3,CLOUGH>::get_continuity() const { return C_ONE; }

// Clough FEMs are not (currently) hierarchic
template <> bool FE<0,CLOUGH>::is_hierarchic() const { return false; } // FIXME - this will be changed
template <> bool FE<1,CLOUGH>::is_hierarchic() const { return false; } // FIXME - this will be changed
template <> bool FE<2,CLOUGH>::is_hierarchic() const { return false; } // FIXME - this will be changed
template <> bool FE<3,CLOUGH>::is_hierarchic() const { return false; } // FIXME - this will be changed

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <>
void FE<2,CLOUGH>::compute_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                        const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <>
void FE<3,CLOUGH>::compute_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                        const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Clough FEM shapes need reinit
template <> bool FE<0,CLOUGH>::shapes_need_reinit() const { return true; }
template <> bool FE<1,CLOUGH>::shapes_need_reinit() const { return true; }
template <> bool FE<2,CLOUGH>::shapes_need_reinit() const { return true; }
template <> bool FE<3,CLOUGH>::shapes_need_reinit() const { return true; }

} // namespace libMesh
