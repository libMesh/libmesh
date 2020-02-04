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

#ifndef LIBMESH_FE_HERMITE_IMPL_H
#define LIBMESH_FE_HERMITE_IMPL_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// ------------------------------------------------------------
// Hermite-specific implementations

// Anonymous namespace for local helper functions
namespace {

template <typename RealType>
void hermite_nodal_soln(const ElemTempl<RealType> * elem,
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
  FEType fe_type(totalorder, HERMITE);

  const unsigned int n_sf =
    // FEShim<Dim,T,RealType>::n_shape_functions(elem_type, totalorder);
    FEInterface::n_shape_functions(Dim, fe_type, elem_type);

  std::vector<PointTempl<RealType>> refspace_nodes;
  FEGenericBase<Real, RealType>::get_refspace_nodes(elem_type,refspace_nodes);
  libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);

  for (unsigned int n=0; n<n_nodes; n++)
    {
      libmesh_assert_equal_to (elem_soln.size(), n_sf);

      // Zero before summation
      nodal_soln[n] = 0;

      // u_i = Sum (alpha_i phi_i)
      for (unsigned int i=0; i<n_sf; i++)
        nodal_soln[n] += elem_soln[i] *
          // FEShim<Dim,T,RealType>::shape(elem, order, i, mapped_point);
          FEInterface::shape(Dim, fe_type, elem, i, refspace_nodes[n]);
    }
} // hermite_nodal_soln()



unsigned int hermite_n_dofs(const ElemType t, const Order o)
{
#ifdef DEBUG
  if (o < 3)
    libmesh_error_msg("Error: Hermite elements require order>=3, but you asked for order=" << o);
#endif

  // Piecewise (bi/tri)cubic C1 Hermite splines
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case EDGE3:
      return (o+1);

    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case QUAD9:
      return ((o+1)*(o+1));

    case HEX8:
    case HEX20:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case HEX27:
      return ((o+1)*(o+1)*(o+1));

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }
} // hermite_n_dofs()




unsigned int hermite_n_dofs_at_node(const ElemType t,
                                    const Order o,
                                    const unsigned int n)
{
  libmesh_assert_greater (o, 2);
  // Piecewise (bi/tri)cubic C1 Hermite splines
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
    case EDGE3:
      {
        switch (n)
          {
          case 0:
          case 1:
            return 2;
          case 2:
            //          Interior DoFs are carried on Elems
            //    return (o-3);
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
          }
      }

    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        switch (n)
          {
            // Vertices
          case 0:
          case 1:
          case 2:
          case 3:
            return 4;
            // Edges
          case 4:
          case 5:
          case 6:
          case 7:
            return (2*(o-3));
          case 8:
            //          Interior DoFs are carried on Elems
            //    return ((o-3)*(o-3));
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD4/8/9!");
          }
      }

    case HEX8:
    case HEX20:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case HEX27:
      {
        switch (n)
          {
            // Vertices
          case 0:
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
          case 6:
          case 7:
            return 8;
            // Edges
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
            return (4*(o-3));
            // Faces
          case 20:
          case 21:
          case 22:
          case 23:
          case 24:
          case 25:
            return (2*(o-3)*(o-3));
          case 26:
            // Interior DoFs are carried on Elems
            //    return ((o-3)*(o-3)*(o-3));
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for HEX8/20/27!");
          }
      }

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }
} // hermite_n_dofs_at_node()



unsigned int hermite_n_dofs_per_elem(const ElemType t,
                                     const Order o)
{
  libmesh_assert_greater (o, 2);

  switch (t)
    {
    case NODEELEM:
      return 0;
    case EDGE2:
    case EDGE3:
      return (o-3);
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      return ((o-3)*(o-3));
    case HEX8:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case HEX20:
    case HEX27:
      return ((o-3)*(o-3)*(o-3));

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }
} // hermite_n_dofs_per_elem()


} // anonymous namespace




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,HERMITE,RealType>::nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               std::vector<Number> & nodal_soln)
{ hermite_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <typename RealType>
void FEShim<1,HERMITE,RealType>::nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               std::vector<Number> & nodal_soln)
{ hermite_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <typename RealType>
void FEShim<2,HERMITE,RealType>::nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               std::vector<Number> & nodal_soln)
{ hermite_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,HERMITE,RealType>::nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               std::vector<Number> & nodal_soln)
{ hermite_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }




// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <typename RealType> unsigned int FEShim<0,HERMITE,RealType>::n_dofs(const ElemType t, const Order o) { return hermite_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,HERMITE,RealType>::n_dofs(const ElemType t, const Order o) { return hermite_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,HERMITE,RealType>::n_dofs(const ElemType t, const Order o) { return hermite_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,HERMITE,RealType>::n_dofs(const ElemType t, const Order o) { return hermite_n_dofs(t, o); }

// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <typename RealType> unsigned int FEShim<0,HERMITE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hermite_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<1,HERMITE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hermite_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<2,HERMITE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hermite_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<3,HERMITE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return hermite_n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <typename RealType> unsigned int FEShim<0,HERMITE,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return hermite_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<1,HERMITE,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return hermite_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<2,HERMITE,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return hermite_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<3,HERMITE,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return hermite_n_dofs_per_elem(t, o); }

// Hermite FEMs are C^1 continuous
template <typename RealType> FEContinuity FEShim<0,HERMITE,RealType>::get_continuity() { return C_ONE; }
template <typename RealType> FEContinuity FEShim<1,HERMITE,RealType>::get_continuity() { return C_ONE; }
template <typename RealType> FEContinuity FEShim<2,HERMITE,RealType>::get_continuity() { return C_ONE; }
template <typename RealType> FEContinuity FEShim<3,HERMITE,RealType>::get_continuity() { return C_ONE; }

// Hermite FEMs are hierarchic
template <typename RealType> bool FEShim<0,HERMITE,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<1,HERMITE,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<2,HERMITE,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<3,HERMITE,RealType>::is_hierarchic() { return true; }


#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <typename RealType>
void FEShim<2,HERMITE,RealType>::compute_constraints (DofConstraints & constraints,
                                         DofMap & dof_map,
                                         const unsigned int variable_number,
                                                      const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<HERMITE>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <typename RealType>
void FEShim<3,HERMITE,RealType>::compute_constraints (DofConstraints & constraints,
                                         DofMap & dof_map,
                                         const unsigned int variable_number,
                                                      const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<HERMITE>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Hermite FEM shapes need reinit
template <typename RealType> bool FEShim<0,HERMITE,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<1,HERMITE,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<2,HERMITE,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<3,HERMITE,RealType>::shapes_need_reinit() { return true; }

} // namespace libMesh

#endif // LIBMESH_FE_HERMITE_IMPL_H
