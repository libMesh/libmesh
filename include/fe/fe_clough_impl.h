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

#ifndef LIBMESH_FE_CLOUGH_IMPL_H
#define LIBMESH_FE_CLOUGH_IMPL_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// ------------------------------------------------------------
// Clough-specific implementations

// Anonymous namespace for local helper functions
namespace {

template <typename RealType>
void clough_nodal_soln(const ElemTempl<RealType> * elem,
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
  FEType fe_type(totalorder, CLOUGH);

  switch (totalorder)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      // Piecewise cubic shape functions
    case THIRD:
      {

        const unsigned int n_sf =
          // FEShim<Dim,T,RealType>::n_shape_functions(elem_type, totalorder);
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
                // FEShim<Dim,T,RealType>::shape(elem, order, i, mapped_point);
                FEInterface::shape(Dim, fe_type, elem, i, refspace_nodes[n]);
          }

        return;
      }

    default:
      libmesh_error_msg("ERROR: Invalid total order " << totalorder);
    }
} // clough_nodal_soln()


unsigned int clough_n_dofs(const ElemType t, const Order o);



unsigned int clough_n_dofs_at_node(const ElemType t,
                                   const Order o,
                                   const unsigned int n);


unsigned int clough_n_dofs_per_elem(const ElemType t, const Order o);

} // anonymous




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,CLOUGH,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ clough_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <typename RealType>
void FEShim<1,CLOUGH,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ clough_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <typename RealType>
void FEShim<2,CLOUGH,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ clough_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,CLOUGH,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ clough_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }


// Full specialization of n_dofs() function for every dimension
template <typename RealType> unsigned int FEShim<0,CLOUGH,RealType>::n_dofs(const ElemType t, const Order o) { return clough_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,CLOUGH,RealType>::n_dofs(const ElemType t, const Order o) { return clough_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,CLOUGH,RealType>::n_dofs(const ElemType t, const Order o) { return clough_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,CLOUGH,RealType>::n_dofs(const ElemType t, const Order o) { return clough_n_dofs(t, o); }


// Full specialization of n_dofs_at_node() function for every dimension.
template <typename RealType> unsigned int FEShim<0,CLOUGH,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return clough_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<1,CLOUGH,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return clough_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<2,CLOUGH,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return clough_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<3,CLOUGH,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return clough_n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <typename RealType> unsigned int FEShim<0,CLOUGH,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return clough_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<1,CLOUGH,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return clough_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<2,CLOUGH,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return clough_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<3,CLOUGH,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return clough_n_dofs_per_elem(t, o); }

// Clough FEMs are C^1 continuous
template <typename RealType> FEContinuity FEShim<0,CLOUGH,RealType>::get_continuity() { return C_ONE; }
template <typename RealType> FEContinuity FEShim<1,CLOUGH,RealType>::get_continuity() { return C_ONE; }
template <typename RealType> FEContinuity FEShim<2,CLOUGH,RealType>::get_continuity() { return C_ONE; }
template <typename RealType> FEContinuity FEShim<3,CLOUGH,RealType>::get_continuity() { return C_ONE; }

// Clough FEMs are not (currently) hierarchic
template <typename RealType> bool FEShim<0,CLOUGH,RealType>::is_hierarchic() { return false; } // FIXME - this will be changed
template <typename RealType> bool FEShim<1,CLOUGH,RealType>::is_hierarchic() { return false; } // FIXME - this will be changed
template <typename RealType> bool FEShim<2,CLOUGH,RealType>::is_hierarchic() { return false; } // FIXME - this will be changed
template <typename RealType> bool FEShim<3,CLOUGH,RealType>::is_hierarchic() { return false; } // FIXME - this will be changed

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <typename RealType>
void FEShim<2,CLOUGH,RealType>::compute_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                                     const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<CLOUGH>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <typename RealType>
void FEShim<3,CLOUGH,RealType>::compute_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                                     const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<CLOUGH>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Clough FEM shapes need reinit
template <typename RealType> bool FEShim<0,CLOUGH,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<1,CLOUGH,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<2,CLOUGH,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<3,CLOUGH,RealType>::shapes_need_reinit() { return true; }

} // namespace libMesh

#endif // LIBMESH_FE_CLOUGH_IMPL_H
