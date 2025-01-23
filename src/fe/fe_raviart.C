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
#include "libmesh/tensor_value.h"


namespace libMesh
{


LIBMESH_DEFAULT_VECTORIZED_FE(0,RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(1,RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(2,RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(3,RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(0,L2_RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(1,L2_RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(2,L2_RAVIART_THOMAS)
LIBMESH_DEFAULT_VECTORIZED_FE(3,L2_RAVIART_THOMAS)


// Anonymous namespace for local helper functions
namespace {
void raviart_thomas_nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               const int dim,
                               std::vector<Number> & nodal_soln,
                               const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();
  const ElemType elem_type   = elem->type();

  const Order totalorder = order + add_p_level*elem->p_level();

  nodal_soln.resize(n_nodes*dim);

  FEType p_refined_fe_type(totalorder, RAVIART_THOMAS);

  if (elem_type != TRI6  && elem_type != TRI7  &&
      elem_type != QUAD8 && elem_type != QUAD9 &&
      elem_type != TET14 && elem_type != HEX27)
    libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(elem_type) << " selected for RAVIART_THOMAS FE family!");

  libmesh_assert_equal_to (elem_soln.size(), FEInterface::n_dofs(p_refined_fe_type, elem, false));

  const unsigned int n_sf = FEInterface::n_shape_functions(p_refined_fe_type, elem, false);

  std::vector<Point> refspace_nodes;
  FEVectorBase::get_refspace_nodes(elem_type,refspace_nodes);
  libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);

  // Need to create new fe object so the shape function has the FETransformation
  // applied to it.
  std::unique_ptr<FEVectorBase> vis_fe = FEVectorBase::build(dim, p_refined_fe_type);

  const std::vector<std::vector<RealGradient>> & vis_phi = vis_fe->get_phi();

  vis_fe->reinit(elem,&refspace_nodes);

  for (unsigned int n = 0; n < n_nodes; n++)
    {
      libmesh_assert_equal_to (elem_soln.size(), n_sf);

      // Zero before summation
      for (int d = 0; d < dim; d++)
        nodal_soln[dim*n+d] = 0;

      // u = Sum (u_i phi_i)
      for (unsigned int i=0; i<n_sf; i++)
        for (int d = 0; d < dim; d++)
          nodal_soln[dim*n+d] += elem_soln[i]*(vis_phi[i][n](d));
    }

  return;
} // raviart_thomas_nodal_soln


unsigned int raviart_thomas_n_dofs(const ElemType t, const Order o)
{
  libmesh_assert_greater (o, 0);
  switch (t)
    {
    case TRI6:
    case TRI7:
      return o*(o+2);
    case QUAD8:
    case QUAD9:
      return 2*o*(o+1);
    case TET14:
      return o*(o+1)*(o+3)/2;
    case HEX27:
      return 3*o*o*(o+1);
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for RAVIART_THOMAS FE family!");
    }
}


unsigned int raviart_thomas_n_dofs_at_node(const ElemType t,
                                           const Order o,
                                           const unsigned int n)
{
  libmesh_assert_greater (o, 0);
  switch (t)
    {
    case TRI6:
    case TRI7:
      {
        switch (n)
          {
          case 0:
          case 1:
          case 2:
            return 0;
          case 3:
          case 4:
          case 5:
            return o;
          case 6:
            libmesh_assert_equal_to(t, TRI7);
            return 0;
          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n);
          }
      }
    case QUAD8:
    case QUAD9:
      {
        switch (n)
          {
          case 0:
          case 1:
          case 2:
          case 3:
            return 0;
          case 4:
          case 5:
          case 6:
          case 7:
            return o;
          case 8:
            libmesh_assert_equal_to(t, QUAD9);
            return 0;
          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n);
          }
      }
    case TET14:
      {
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
          case 8:
          case 9:
            return 0;
          case 10:
          case 11:
          case 12:
          case 13:
            return o*(o+1)/2;
          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n);
          }
      }
    case HEX27:
      {
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
            return 0;
          case 20:
          case 21:
          case 22:
          case 23:
          case 24:
          case 25:
            return o*o;
          case 26:
            return 0;
          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n);
          }
      }
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for RAVIART_THOMAS FE family!");
    }
}


unsigned int raviart_thomas_n_dofs_per_elem(const ElemType t,
                                            const Order o)
{
  libmesh_assert_greater (o, 0);
  switch (t)
    {
    case TRI6:
    case TRI7:
      return o*(o-1);
    case QUAD8:
    case QUAD9:
      return 2*o*(o-1);
    case TET14:
      return (o+1)*o*(o-1)/2;
    case HEX27:
      return 3*o*o*(o-1);
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for RAVIART_THOMAS FE family!");
    }
}


#ifdef LIBMESH_ENABLE_AMR
void raviart_thomas_compute_constraints (DofConstraints & /*constraints*/,
                                         DofMap & /*dof_map*/,
                                         const unsigned int /*variable_number*/,
                                         const Elem * libmesh_dbg_var(elem),
                                         const unsigned Dim)
{
  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  libmesh_assert(elem);

  libmesh_not_implemented();
} // raviart_thomas_compute_constraints()
#endif // #ifdef LIBMESH_ENABLE_AMR

} // anonymous namespace

#define RAVIART_LOW_D_ERROR_MESSAGE                                     \
  libmesh_error_msg("ERROR: This method makes no sense for low-D elements!");


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this file.
template <>
void FE<0,RAVIART_THOMAS>::nodal_soln(const Elem *,
                                      const Order,
                                      const std::vector<Number> &,
                                      std::vector<Number> &,
                                      bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<1,RAVIART_THOMAS>::nodal_soln(const Elem *,
                                      const Order,
                                      const std::vector<Number> &,
                                      std::vector<Number> &,
                                      bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<2,RAVIART_THOMAS>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln,
                                      const bool add_p_level)
{ raviart_thomas_nodal_soln(elem, order, elem_soln, 2 /*dim*/, nodal_soln, add_p_level); }

template <>
void FE<3,RAVIART_THOMAS>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln,
                                      const bool add_p_level)
{ raviart_thomas_nodal_soln(elem, order, elem_soln, 3 /*dim*/, nodal_soln, add_p_level); }

template <>
void FE<0,L2_RAVIART_THOMAS>::nodal_soln(const Elem *,
                                         const Order,
                                         const std::vector<Number> &,
                                         std::vector<Number> &,
                                         bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<1,L2_RAVIART_THOMAS>::nodal_soln(const Elem *,
                                         const Order,
                                         const std::vector<Number> &,
                                         std::vector<Number> &,
                                         bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<2,L2_RAVIART_THOMAS>::nodal_soln(const Elem * elem,
                                         const Order order,
                                         const std::vector<Number> & elem_soln,
                                         std::vector<Number> & nodal_soln,
                                         const bool add_p_level)
{ raviart_thomas_nodal_soln(elem, order, elem_soln, 2 /*dim*/, nodal_soln, add_p_level); }

template <>
void FE<3,L2_RAVIART_THOMAS>::nodal_soln(const Elem * elem,
                                         const Order order,
                                         const std::vector<Number> & elem_soln,
                                         std::vector<Number> & nodal_soln,
                                         const bool add_p_level)
{ raviart_thomas_nodal_soln(elem, order, elem_soln, 3 /*dim*/, nodal_soln, add_p_level); }

LIBMESH_FE_SIDE_NODAL_SOLN(RAVIART_THOMAS)
LIBMESH_FE_SIDE_NODAL_SOLN(L2_RAVIART_THOMAS)


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <> unsigned int FE<0,RAVIART_THOMAS>::n_dofs(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,RAVIART_THOMAS>::n_dofs(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,RAVIART_THOMAS>::n_dofs(const ElemType t, const Order o) { return raviart_thomas_n_dofs(t, o); }
template <> unsigned int FE<3,RAVIART_THOMAS>::n_dofs(const ElemType t, const Order o) { return raviart_thomas_n_dofs(t, o); }
template <> unsigned int FE<0,L2_RAVIART_THOMAS>::n_dofs(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,L2_RAVIART_THOMAS>::n_dofs(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,L2_RAVIART_THOMAS>::n_dofs(const ElemType t, const Order o) { return raviart_thomas_n_dofs(t, o); }
template <> unsigned int FE<3,L2_RAVIART_THOMAS>::n_dofs(const ElemType t, const Order o) { return raviart_thomas_n_dofs(t, o); }

template <> unsigned int FE<0,RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,RAVIART_THOMAS>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return raviart_thomas_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,RAVIART_THOMAS>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return raviart_thomas_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<0,L2_RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,L2_RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,L2_RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType t, const Order o) { return raviart_thomas_n_dofs_per_elem(t, o); }
template <> unsigned int FE<3,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType t, const Order o) { return raviart_thomas_n_dofs_per_elem(t, o); }

// L2 Raviart-Thomas elements have all their dofs per element
template <> unsigned int FE<0,L2_RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,L2_RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,L2_RAVIART_THOMAS>::n_dofs_per_elem(const ElemType t, const Order o) { return n_dofs(t, o); }
template <> unsigned int FE<3,L2_RAVIART_THOMAS>::n_dofs_per_elem(const ElemType t, const Order o) { return n_dofs(t, o); }

// Raviart-Thomas FEMs are always normally continuous
template <> FEContinuity FE<0,RAVIART_THOMAS>::get_continuity() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> FEContinuity FE<1,RAVIART_THOMAS>::get_continuity() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> FEContinuity FE<2,RAVIART_THOMAS>::get_continuity() const { return H_DIV; }
template <> FEContinuity FE<3,RAVIART_THOMAS>::get_continuity() const { return H_DIV; }

// L2 Raviart-Thomas FEMs are discontinuous
template <> FEContinuity FE<0,L2_RAVIART_THOMAS>::get_continuity() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> FEContinuity FE<1,L2_RAVIART_THOMAS>::get_continuity() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> FEContinuity FE<2,L2_RAVIART_THOMAS>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,L2_RAVIART_THOMAS>::get_continuity() const { return DISCONTINUOUS; }

// Raviart-Thomas FEMs are not hierarchic
template <> bool FE<0,RAVIART_THOMAS>::is_hierarchic() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<1,RAVIART_THOMAS>::is_hierarchic() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<2,RAVIART_THOMAS>::is_hierarchic() const { return false; }
template <> bool FE<3,RAVIART_THOMAS>::is_hierarchic() const { return false; }
template <> bool FE<0,L2_RAVIART_THOMAS>::is_hierarchic() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<1,L2_RAVIART_THOMAS>::is_hierarchic() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<2,L2_RAVIART_THOMAS>::is_hierarchic() const { return false; }
template <> bool FE<3,L2_RAVIART_THOMAS>::is_hierarchic() const { return false; }

// Raviart-Thomas FEM shapes always need to be reinit'ed (because of orientation dependence)
template <> bool FE<0,RAVIART_THOMAS>::shapes_need_reinit() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<1,RAVIART_THOMAS>::shapes_need_reinit() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<2,RAVIART_THOMAS>::shapes_need_reinit() const { return true; }
template <> bool FE<3,RAVIART_THOMAS>::shapes_need_reinit() const { return true; }

// Do L2 Raviart-Thomas FEM shapes always need to be reinit'ed because of orientation dependence?
template <> bool FE<0,L2_RAVIART_THOMAS>::shapes_need_reinit() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<1,L2_RAVIART_THOMAS>::shapes_need_reinit() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<2,L2_RAVIART_THOMAS>::shapes_need_reinit() const { return true; }
template <> bool FE<3,L2_RAVIART_THOMAS>::shapes_need_reinit() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<0,RAVIART_THOMAS>::compute_constraints (DofConstraints &,
                                                DofMap &,
                                                const unsigned int,
                                                const Elem *)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<1,RAVIART_THOMAS>::compute_constraints (DofConstraints &,
                                                DofMap &,
                                                const unsigned int,
                                                const Elem *)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<2,RAVIART_THOMAS>::compute_constraints (DofConstraints & constraints,
                                                DofMap & dof_map,
                                                const unsigned int variable_number,
                                                const Elem * elem)
{ raviart_thomas_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/2); }

template <>
void FE<3,RAVIART_THOMAS>::compute_constraints (DofConstraints & constraints,
                                                DofMap & dof_map,
                                                const unsigned int variable_number,
                                                const Elem * elem)
{ raviart_thomas_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/3); }

template <>
void FE<0,L2_RAVIART_THOMAS>::compute_constraints (DofConstraints &,
                                                   DofMap &,
                                                   const unsigned int,
                                                   const Elem *)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<1,L2_RAVIART_THOMAS>::compute_constraints (DofConstraints &,
                                                   DofMap &,
                                                   const unsigned int,
                                                   const Elem *)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<2,L2_RAVIART_THOMAS>::compute_constraints (DofConstraints & constraints,
                                                   DofMap & dof_map,
                                                   const unsigned int variable_number,
                                                   const Elem * elem)
{ raviart_thomas_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/2); }

template <>
void FE<3,L2_RAVIART_THOMAS>::compute_constraints (DofConstraints & constraints,
                                                   DofMap & dof_map,
                                                   const unsigned int variable_number,
                                                   const Elem * elem)
{ raviart_thomas_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/3); }
#endif // LIBMESH_ENABLE_AMR

// Specialize useless shape function methods
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape(const ElemType, const Order,const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape(const Elem *,const Order,const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,L2_RAVIART_THOMAS>::shape(const ElemType, const Order,const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,L2_RAVIART_THOMAS>::shape(const Elem *,const Order,const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape_deriv(const ElemType, const Order,const unsigned int,
                                               const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape_deriv(const Elem *,const Order,const unsigned int,
                                               const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,L2_RAVIART_THOMAS>::shape_deriv(const ElemType, const Order,const unsigned int,
                                                  const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,L2_RAVIART_THOMAS>::shape_deriv(const Elem *,const Order,const unsigned int,
                                                  const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape_second_deriv(const ElemType, const Order,const unsigned int,
                                                      const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape_second_deriv(const Elem *,const Order,const unsigned int,
                                                      const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,L2_RAVIART_THOMAS>::shape_second_deriv(const ElemType, const Order,const unsigned int,
                                                         const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,L2_RAVIART_THOMAS>::shape_second_deriv(const Elem *,const Order,const unsigned int,
                                                         const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

#endif

template <>
RealGradient FE<1,RAVIART_THOMAS>::shape(const ElemType, const Order,const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape(const Elem *,const Order,const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,L2_RAVIART_THOMAS>::shape(const ElemType, const Order,const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,L2_RAVIART_THOMAS>::shape(const Elem *,const Order,const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape_deriv(const ElemType, const Order,const unsigned int,
                                               const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape_deriv(const Elem *,const Order,const unsigned int,
                                               const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,L2_RAVIART_THOMAS>::shape_deriv(const ElemType, const Order,const unsigned int,
                                                  const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,L2_RAVIART_THOMAS>::shape_deriv(const Elem *,const Order,const unsigned int,
                                                  const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape_second_deriv(const ElemType, const Order,const unsigned int,
                                                      const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape_second_deriv(const Elem *,const Order,const unsigned int,
                                                      const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,L2_RAVIART_THOMAS>::shape_second_deriv(const ElemType, const Order,const unsigned int,
                                                         const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,L2_RAVIART_THOMAS>::shape_second_deriv(const Elem *,const Order,const unsigned int,
                                                         const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }

#endif

} // namespace libMesh
