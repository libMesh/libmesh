// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Anonymous namespace for local helper functions
namespace {
void raviart_thomas_nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               const int dim,
                               std::vector<Number> & nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  const ElemType elem_type   = elem->type();

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  nodal_soln.resize(n_nodes*dim);

  FEType fe_type(order, RAVIART_THOMAS);
  FEType p_refined_fe_type(totalorder, RAVIART_THOMAS);

  switch (totalorder)
    {
    case FIRST:
      {
        switch (elem_type)
          {
          case TRI6:
          case TRI7:
            {
              libmesh_assert_equal_to (elem_soln.size(), 3);

              if (elem_type == TRI6)
                libmesh_assert_equal_to (nodal_soln.size(), 6*2);
              else
                libmesh_assert_equal_to (nodal_soln.size(), 7*2);
              break;
            }
          case QUAD8:
          case QUAD9:
            {
              libmesh_assert_equal_to (elem_soln.size(), 4);

              if (elem_type == QUAD8)
                libmesh_assert_equal_to (nodal_soln.size(), 8*2);
              else
                libmesh_assert_equal_to (nodal_soln.size(), 9*2);
              break;
            }
          case TET14:
            {
              libmesh_assert_equal_to (elem_soln.size(), 4);
              libmesh_assert_equal_to (nodal_soln.size(), 14*3);
              break;
            }
          case HEX27:
            {
              libmesh_assert_equal_to (elem_soln.size(), 6);
              libmesh_assert_equal_to (nodal_soln.size(), 27*3);
              break;
            }

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(elem_type) << " selected for RAVIART_THOMAS FE family!");

          } // switch(elem_type)

        const unsigned int n_sf =
          FEInterface::n_shape_functions(fe_type, elem);

        std::vector<Point> refspace_nodes;
        FEVectorBase::get_refspace_nodes(elem_type,refspace_nodes);
        libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);


        // Need to create new fe object so the shape function as the FETransformation
        // applied to it.
        std::unique_ptr<FEVectorBase> vis_fe = FEVectorBase::build(dim, p_refined_fe_type);

        const std::vector<std::vector<RealGradient>> & vis_phi = vis_fe->get_phi();

        vis_fe->reinit(elem,&refspace_nodes);

        for (unsigned int n = 0; n < n_nodes; n++)
          {
            libmesh_assert_equal_to (elem_soln.size(), n_sf);

            // Zero before summation
            for (int d = 0; d < dim; d++)
              {
                nodal_soln[dim*n+d] = 0;
              }

            // u = Sum (u_i phi_i)
            for (unsigned int i=0; i<n_sf; i++)
              {
                for (int d = 0; d < dim; d++)
                  {
                    nodal_soln[dim*n+d]   += elem_soln[i]*(vis_phi[i][n](d));
                  }
              }
          }

        return;
      } // case FIRST

    default:
      libmesh_error_msg("ERROR: Invalid total order " << Utility::enum_to_string(totalorder) << " selected for RAVIART_THOMAS FE family!");

    }//switch (totalorder)

  return;
} // raviart_thomas_nodal_soln


unsigned int raviart_thomas_n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
    case FIRST:
      {
        switch (t)
          {
          case TRI6:
          case TRI7:
            return 3;

          case QUAD8:
          case QUAD9:
            return 4;

          case TET14:
            return 4;

          case HEX27:
            return 6;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for RAVIART_THOMAS FE family!");
    }
}




unsigned int raviart_thomas_n_dofs_at_node(const ElemType t,
                                           const Order o,
                                           const unsigned int n)
{
  switch (o)
    {
    case FIRST:
      {
        switch (t)
          {
          case TRI6:
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
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n);
                }
            }
          case TRI7:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 6:
                  return 0;
                case 3:
                case 4:
                case 5:
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n);
                }
            }
          case QUAD8:
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
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n);
                }
            }
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                case 8:
                  return 0;
                case 4:
                case 5:
                case 6:
                case 7:
                  return 1;

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
                  return 1;

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
                case 26:
                  return 0;
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n);
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for RAVIART_THOMAS FE family!");
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
                                      std::vector<Number> &)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<1,RAVIART_THOMAS>::nodal_soln(const Elem *,
                                      const Order,
                                      const std::vector<Number> &,
                                      std::vector<Number> &)
{ RAVIART_LOW_D_ERROR_MESSAGE }

template <>
void FE<2,RAVIART_THOMAS>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln)
{ raviart_thomas_nodal_soln(elem, order, elem_soln, 2 /*dim*/, nodal_soln); }

template <>
void FE<3,RAVIART_THOMAS>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln)
{ raviart_thomas_nodal_soln(elem, order, elem_soln, 3 /*dim*/, nodal_soln); }

LIBMESH_FE_SIDE_NODAL_SOLN(RAVIART_THOMAS)


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <> unsigned int FE<0,RAVIART_THOMAS>::n_dofs(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,RAVIART_THOMAS>::n_dofs(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,RAVIART_THOMAS>::n_dofs(const ElemType t, const Order o) { return raviart_thomas_n_dofs(t, o); }
template <> unsigned int FE<3,RAVIART_THOMAS>::n_dofs(const ElemType t, const Order o) { return raviart_thomas_n_dofs(t, o); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <> unsigned int FE<0,RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,RAVIART_THOMAS>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,RAVIART_THOMAS>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return raviart_thomas_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,RAVIART_THOMAS>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return raviart_thomas_n_dofs_at_node(t, o, n); }


// Raviart-Thomas elements have no dofs per element
template <> unsigned int FE<0,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<1,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { RAVIART_LOW_D_ERROR_MESSAGE }
template <> unsigned int FE<2,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<3,RAVIART_THOMAS>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Raviart-Thomas FEMs are always normally continuous
template <> FEContinuity FE<0,RAVIART_THOMAS>::get_continuity() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> FEContinuity FE<1,RAVIART_THOMAS>::get_continuity() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> FEContinuity FE<2,RAVIART_THOMAS>::get_continuity() const { return H_DIV; }
template <> FEContinuity FE<3,RAVIART_THOMAS>::get_continuity() const { return H_DIV; }

// Raviart-Thomas FEMs are not hierarchic
template <> bool FE<0,RAVIART_THOMAS>::is_hierarchic() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<1,RAVIART_THOMAS>::is_hierarchic() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<2,RAVIART_THOMAS>::is_hierarchic() const { return false; }
template <> bool FE<3,RAVIART_THOMAS>::is_hierarchic() const { return false; }

// Raviart-Thomas FEM shapes always need to be reinit'ed (because of orientation dependence)
template <> bool FE<0,RAVIART_THOMAS>::shapes_need_reinit() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<1,RAVIART_THOMAS>::shapes_need_reinit() const { RAVIART_LOW_D_ERROR_MESSAGE }
template <> bool FE<2,RAVIART_THOMAS>::shapes_need_reinit() const { return true; }
template <> bool FE<3,RAVIART_THOMAS>::shapes_need_reinit() const { return true; }

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
#endif // LIBMESH_ENABLE_AMR

// Specialize useless shape function methods
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape(const ElemType, const Order,const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape(const Elem *,const Order,const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape_deriv(const ElemType, const Order,const unsigned int,
                                               const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<0,RAVIART_THOMAS>::shape_deriv(const Elem *,const Order,const unsigned int,
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

#endif

template <>
RealGradient FE<1,RAVIART_THOMAS>::shape(const ElemType, const Order,const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape(const Elem *,const Order,const unsigned int,const Point &,const bool)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape_deriv(const ElemType, const Order,const unsigned int,
                                               const unsigned int,const Point &)
{ RAVIART_LOW_D_ERROR_MESSAGE }
template <>
RealGradient FE<1,RAVIART_THOMAS>::shape_deriv(const Elem *,const Order,const unsigned int,
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

#endif

} // namespace libMesh
