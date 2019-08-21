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

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"
#include "libmesh/tensor_value.h"

namespace libMesh
{

// ------------------------------------------------------------
// Vector monomial specific implementations

// Anonymous namespace for local helper functions
namespace
{
void
monomial_vec_nodal_soln(const Elem * elem,
                        const Order order,
                        const std::vector<Number> & elem_soln,
                        const int dim,
                        std::vector<Number> & nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(dim * n_nodes);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (totalorder)
  {
      // Constant shape functions
    case CONSTANT:
    {
      switch (dim)
      {
        case 2:
        {
          libmesh_assert_equal_to(elem_soln.size(), 2);
          for (unsigned int n = 0; n < n_nodes; n++)
          {
            nodal_soln[2 * n] = elem_soln[0];
            nodal_soln[1 + 2 * n] = elem_soln[1];
          }
          return;
        }
        case 3:
        {
          libmesh_assert_equal_to(elem_soln.size(), 3);
          for (unsigned int n = 0; n < n_nodes; n++)
          {
            nodal_soln[3 * n] = elem_soln[0];
            nodal_soln[1 + 3 * n] = elem_soln[1];
            nodal_soln[2 + 3 * n] = elem_soln[2];
          }
          return;
        }
        default:
          libmesh_error_msg(
              "The monomial_vec_nodal_soln helper should only be called for 2 and 3 dimensions");
      }
    }

    // For other orders, do interpolation at the nodes
    // explicitly.
    default:
    {
      // FEType object to be passed to various FEInterface functions below.
      FEType fe_type(totalorder, MONOMIAL);

      const unsigned int n_sf = FEInterface::n_shape_functions(dim, fe_type, elem_type);

      std::vector<Point> refspace_nodes;
      FEBase::get_refspace_nodes(elem_type, refspace_nodes);
      libmesh_assert_equal_to(refspace_nodes.size(), n_nodes);
      libmesh_assert_equal_to(elem_soln.size(), n_sf * dim);

      for (unsigned int d = 0; d < static_cast<unsigned int>(dim); d++)
        for (unsigned int n = 0; n < n_nodes; n++)
        {

          // Zero before summation
          nodal_soln[d + dim * n] = 0;

          // u_i = Sum (alpha_i phi_i)
          for (unsigned int i = 0; i < n_sf; i++)
            nodal_soln[d + dim * n] += elem_soln[d + dim * i] *
                                       FEInterface::shape(dim, fe_type, elem, i, refspace_nodes[n]);
        }

      return;
    } // default
  }   // switch
}
} // anonymous namespace

// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this file.
// This could be macro-ified so that it fits on one line...
template <>
void
FE<0, MONOMIAL_VEC>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{
  FE<0, MONOMIAL>::nodal_soln(elem, order, elem_soln, nodal_soln);
}

template <>
void
FE<1, MONOMIAL_VEC>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{
  FE<1, MONOMIAL>::nodal_soln(elem, order, elem_soln, nodal_soln);
}

template <>
void
FE<2, MONOMIAL_VEC>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{
  monomial_vec_nodal_soln(elem, order, elem_soln, 2 /*dimension*/, nodal_soln);
}

template <>
void
FE<3, MONOMIAL_VEC>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{
  monomial_vec_nodal_soln(elem, order, elem_soln, 3 /*dimension*/, nodal_soln);
}

// Specialize for shape function routines by leveraging scalar MONOMIAL elements

// 0-D
template <>
RealVectorValue
FE<0, MONOMIAL_VEC>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value = FE<0, MONOMIAL>::shape(type, order, i, p);
  return libMesh::RealVectorValue(value);
}
template <>
RealVectorValue
FE<0, MONOMIAL_VEC>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<0, MONOMIAL>::shape_deriv(type, order, i, j, p);
  return libMesh::RealVectorValue(value);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealVectorValue
FE<0, MONOMIAL_VEC>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<0, MONOMIAL>::shape_second_deriv(type, order, i, j, p);
  return libMesh::RealVectorValue(value);
}
#endif

// 1-D
template <>
RealVectorValue
FE<1, MONOMIAL_VEC>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value = FE<1, MONOMIAL>::shape(type, order, i, p);
  return libMesh::RealVectorValue(value);
}
template <>
RealVectorValue
FE<1, MONOMIAL_VEC>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<1, MONOMIAL>::shape_deriv(type, order, i, j, p);
  return libMesh::RealVectorValue(value);
}
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
RealVectorValue
FE<1, MONOMIAL_VEC>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<1, MONOMIAL>::shape_second_deriv(type, order, i, j, p);
  return libMesh::RealVectorValue(value);
}

#endif

// 2-D
template <>
RealVectorValue
FE<2, MONOMIAL_VEC>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value = FE<2, MONOMIAL>::shape(type, order, i / 2, p);

  switch (i % 2)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
  }

  // dummy
  return libMesh::RealVectorValue();
}
template <>
RealVectorValue
FE<2, MONOMIAL_VEC>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<2, MONOMIAL>::shape_deriv(type, order, i / 2, j, p);

  switch (i % 2)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
RealVectorValue
FE<2, MONOMIAL_VEC>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<2, MONOMIAL>::shape_second_deriv(type, order, i / 2, j, p);

  switch (i % 2)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#endif

// 3-D
template <>
RealVectorValue
FE<3, MONOMIAL_VEC>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value = FE<3, MONOMIAL>::shape(type, order, i / 3, p);

  switch (i % 3)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    case 2:
      return libMesh::RealVectorValue(Real(0), Real(0), value);

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
  }

  // dummy
  return libMesh::RealVectorValue();
}
template <>
RealVectorValue
FE<3, MONOMIAL_VEC>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<3, MONOMIAL>::shape_deriv(type, order, i / 3, j, p);

  switch (i % 3)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    case 2:
      return libMesh::RealVectorValue(Real(0), Real(0), value);

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealVectorValue
FE<3, MONOMIAL_VEC>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<3, MONOMIAL>::shape_second_deriv(type, order, i / 3, j, p);

  switch (i % 3)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    case 2:
      return libMesh::RealVectorValue(Real(0), Real(0), value);

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#endif

// 0-D
template <>
RealVectorValue
FE<0, MONOMIAL_VEC>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value =
      FE<0, MONOMIAL>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
  return libMesh::RealVectorValue(value);
}
template <>
RealVectorValue
FE<0, MONOMIAL_VEC>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<0, MONOMIAL>::shape_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealVectorValue(value);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealVectorValue
FE<0, MONOMIAL_VEC>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<0, MONOMIAL>::shape_second_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealVectorValue(value);
}

#endif

// 1-D
template <>
RealVectorValue
FE<1, MONOMIAL_VEC>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value =
      FE<1, MONOMIAL>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
  return libMesh::RealVectorValue(value);
}
template <>
RealVectorValue
FE<1, MONOMIAL_VEC>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<1, MONOMIAL>::shape_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealVectorValue(value);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
RealVectorValue
FE<1, MONOMIAL_VEC>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<1, MONOMIAL>::shape_second_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealVectorValue(value);
}

#endif

// 2-D
template <>
RealVectorValue
FE<2, MONOMIAL_VEC>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value =
      FE<2, MONOMIAL>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i / 2, p);

  switch (i % 2)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
  }

  // dummy
  return libMesh::RealVectorValue();
}
template <>
RealVectorValue
FE<2, MONOMIAL_VEC>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<2, MONOMIAL>::shape_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i / 2, j, p);

  switch (i % 2)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
RealVectorValue
FE<2, MONOMIAL_VEC>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<2, MONOMIAL>::shape_second_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i / 2, j, p);

  switch (i % 2)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#endif

// 3-D
template <>
RealVectorValue
FE<3, MONOMIAL_VEC>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  Real value =
      FE<3, MONOMIAL>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i / 3, p);

  switch (i % 3)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    case 2:
      return libMesh::RealVectorValue(Real(0), Real(0), value);

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
  }

  // dummy
  return libMesh::RealVectorValue();
}
template <>
RealVectorValue
FE<3, MONOMIAL_VEC>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  Real value = FE<3, MONOMIAL>::shape_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i / 3, j, p);

  switch (i % 3)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    case 2:
      return libMesh::RealVectorValue(Real(0), Real(0), value);

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealVectorValue
FE<3, MONOMIAL_VEC>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  Real value = FE<3, MONOMIAL>::shape_second_deriv(
      elem->type(), static_cast<Order>(order + elem->p_level()), i / 3, j, p);

  switch (i % 3)
  {
    case 0:
      return libMesh::RealVectorValue(value);

    case 1:
      return libMesh::RealVectorValue(Real(0), value);

    case 2:
      return libMesh::RealVectorValue(Real(0), Real(0), value);

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
  }

  // dummy
  return libMesh::RealVectorValue();
}

#endif

// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <>
unsigned int
FE<0, MONOMIAL_VEC>::n_dofs(const ElemType t, const Order o)
{
  return FE<0, MONOMIAL>::n_dofs(t, o);
}
template <>
unsigned int
FE<1, MONOMIAL_VEC>::n_dofs(const ElemType t, const Order o)
{
  return FE<1, MONOMIAL>::n_dofs(t, o);
}
template <>
unsigned int
FE<2, MONOMIAL_VEC>::n_dofs(const ElemType t, const Order o)
{
  return 2 * FE<2, MONOMIAL>::n_dofs(t, o);
}
template <>
unsigned int
FE<3, MONOMIAL_VEC>::n_dofs(const ElemType t, const Order o)
{
  return 3 * FE<3, MONOMIAL>::n_dofs(t, o);
}

template <>
unsigned int
FE<0, MONOMIAL_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int)
{
  return 0;
}
template <>
unsigned int
FE<1, MONOMIAL_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int)
{
  return 0;
}
template <>
unsigned int
FE<2, MONOMIAL_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int)
{
  return 0;
}
template <>
unsigned int
FE<3, MONOMIAL_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int)
{
  return 0;
}

template <>
unsigned int
FE<0, MONOMIAL_VEC>::n_dofs_per_elem(const ElemType t, const Order o)
{
  return FE<0, MONOMIAL>::n_dofs_per_elem(t, o);
}
template <>
unsigned int
FE<1, MONOMIAL_VEC>::n_dofs_per_elem(const ElemType t, const Order o)
{
  return FE<1, MONOMIAL>::n_dofs_per_elem(t, o);
}
template <>
unsigned int
FE<2, MONOMIAL_VEC>::n_dofs_per_elem(const ElemType t, const Order o)
{
  return 2 * FE<2, MONOMIAL>::n_dofs_per_elem(t, o);
}
template <>
unsigned int
FE<3, MONOMIAL_VEC>::n_dofs_per_elem(const ElemType t, const Order o)
{
  return 3 * FE<3, MONOMIAL>::n_dofs_per_elem(t, o);
}

// Monomial FEMs are always C^0 continuous
template <>
FEContinuity
FE<0, MONOMIAL_VEC>::get_continuity() const
{
  return DISCONTINUOUS;
}
template <>
FEContinuity
FE<1, MONOMIAL_VEC>::get_continuity() const
{
  return DISCONTINUOUS;
}
template <>
FEContinuity
FE<2, MONOMIAL_VEC>::get_continuity() const
{
  return DISCONTINUOUS;
}
template <>
FEContinuity
FE<3, MONOMIAL_VEC>::get_continuity() const
{
  return DISCONTINUOUS;
}

// Monomial FEMs are not hierarchic
template <>
bool
FE<0, MONOMIAL_VEC>::is_hierarchic() const
{
  return true;
}
template <>
bool
FE<1, MONOMIAL_VEC>::is_hierarchic() const
{
  return true;
}
template <>
bool
FE<2, MONOMIAL_VEC>::is_hierarchic() const
{
  return true;
}
template <>
bool
FE<3, MONOMIAL_VEC>::is_hierarchic() const
{
  return true;
}

// Monomial FEM shapes do not need reinit (is this always true?)
template <>
bool
FE<0, MONOMIAL_VEC>::shapes_need_reinit() const
{
  return false;
}
template <>
bool
FE<1, MONOMIAL_VEC>::shapes_need_reinit() const
{
  return false;
}
template <>
bool
FE<2, MONOMIAL_VEC>::shapes_need_reinit() const
{
  return false;
}
template <>
bool
FE<3, MONOMIAL_VEC>::shapes_need_reinit() const
{
  return false;
}

// we only need instantiations of this function for
// Dim==2 and 3. There are no constraints for discontinuous elements, so
// we do nothing.
#ifdef LIBMESH_ENABLE_AMR
template <>
void
FE<2, MONOMIAL_VEC>::compute_constraints(DofConstraints &,
                                         DofMap &,
                                         const unsigned int,
                                         const Elem *)
{
}

template <>
void
FE<3, MONOMIAL_VEC>::compute_constraints(DofConstraints &,
                                         DofMap &,
                                         const unsigned int,
                                         const Elem *)
{
}
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
