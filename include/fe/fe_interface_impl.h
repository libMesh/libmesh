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

#ifndef FE_INTERFACE_IMPL_H
#define FE_INTERFACE_IMPL_H

// Local includes
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
#define fe_family_switch(dim, func_and_args, prefix, suffix, RealType)           \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case CLOUGH:                                                      \
        prefix FE<dim,CLOUGH,RealType>::func_and_args suffix                     \
      case HERMITE:                                                     \
        prefix FE<dim,HERMITE,RealType>::func_and_args suffix                    \
      case HIERARCHIC:                                                  \
        prefix FE<dim,HIERARCHIC,RealType>::func_and_args suffix                 \
      case L2_HIERARCHIC:                                               \
        prefix FE<dim,L2_HIERARCHIC,RealType>::func_and_args suffix              \
      case LAGRANGE:                                                    \
        prefix FE<dim,LAGRANGE,RealType>::func_and_args suffix                   \
      case L2_LAGRANGE:                                                 \
        prefix FE<dim,L2_LAGRANGE,RealType>::func_and_args suffix                \
      case MONOMIAL:                                                    \
        prefix FE<dim,MONOMIAL,RealType>::func_and_args suffix                   \
      case SCALAR:                                                      \
        prefix FE<dim,SCALAR,RealType>::func_and_args suffix                     \
      case BERNSTEIN:                                                   \
        prefix FE<dim,BERNSTEIN,RealType>::func_and_args suffix                  \
      case SZABAB:                                                      \
        prefix FE<dim,SZABAB,RealType>::func_and_args suffix                     \
      case RATIONAL_BERNSTEIN:                                          \
        prefix FE<dim,RATIONAL_BERNSTEIN,RealType>::func_and_args suffix         \
      case XYZ:                                                         \
        prefix FEXYZ<dim,RealType>::func_and_args suffix                         \
      case SUBDIVISION:                                                 \
        libmesh_assert_equal_to (dim, 2);                               \
        prefix FE<2,SUBDIVISION,RealType>::func_and_args suffix                  \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)

#define fe_family_with_vec_switch(dim, func_and_args, prefix, suffix, RealType)  \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case CLOUGH:                                                      \
        prefix FE<dim,CLOUGH,RealType>::func_and_args suffix                     \
      case HERMITE:                                                     \
        prefix FE<dim,HERMITE,RealType>::func_and_args suffix                    \
      case HIERARCHIC:                                                  \
        prefix FE<dim,HIERARCHIC,RealType>::func_and_args suffix                 \
      case L2_HIERARCHIC:                                               \
        prefix FE<dim,L2_HIERARCHIC,RealType>::func_and_args suffix              \
      case LAGRANGE:                                                    \
        prefix FE<dim,LAGRANGE,RealType>::func_and_args suffix                   \
      case LAGRANGE_VEC:                                                \
        prefix FELagrangeVec<dim,RealType>::func_and_args suffix                 \
      case L2_LAGRANGE:                                                 \
        prefix FE<dim,L2_LAGRANGE,RealType>::func_and_args suffix                \
      case MONOMIAL:                                                    \
        prefix FE<dim,MONOMIAL,RealType>::func_and_args suffix                   \
      case MONOMIAL_VEC:                                                \
        prefix FEMonomialVec<dim,RealType>::func_and_args suffix                 \
      case SCALAR:                                                      \
        prefix FE<dim,SCALAR,RealType>::func_and_args suffix                     \
      case BERNSTEIN:                                                   \
        prefix FE<dim,BERNSTEIN,RealType>::func_and_args suffix                  \
      case SZABAB:                                                      \
        prefix FE<dim,SZABAB,RealType>::func_and_args suffix                     \
      case RATIONAL_BERNSTEIN:                                          \
        prefix FE<dim,RATIONAL_BERNSTEIN,RealType>::func_and_args suffix         \
      case XYZ:                                                         \
        prefix FEXYZ<dim,RealType>::func_and_args suffix                         \
      case SUBDIVISION:                                                 \
        libmesh_assert_equal_to (dim, 2);                               \
        prefix FE<2,SUBDIVISION,RealType>::func_and_args suffix                  \
      case NEDELEC_ONE:                                                 \
        prefix FENedelecOne<dim,RealType>::func_and_args suffix                  \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)

#define fe_scalar_vec_error_switch(dim, func_and_args, prefix, suffix, RealType) \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case CLOUGH:                                                      \
        prefix FE<dim,CLOUGH,RealType>::func_and_args suffix                     \
      case HERMITE:                                                     \
        prefix FE<dim,HERMITE,RealType>::func_and_args suffix                    \
      case HIERARCHIC:                                                  \
        prefix FE<dim,HIERARCHIC,RealType>::func_and_args suffix                 \
      case L2_HIERARCHIC:                                               \
        prefix FE<dim,L2_HIERARCHIC,RealType>::func_and_args suffix              \
      case LAGRANGE:                                                    \
        prefix FE<dim,LAGRANGE,RealType>::func_and_args suffix                   \
      case L2_LAGRANGE:                                                 \
        prefix FE<dim,L2_LAGRANGE,RealType>::func_and_args suffix                \
      case MONOMIAL:                                                    \
        prefix FE<dim,MONOMIAL,RealType>::func_and_args suffix                   \
      case SCALAR:                                                      \
        prefix FE<dim,SCALAR,RealType>::func_and_args suffix                     \
      case BERNSTEIN:                                                   \
        prefix FE<dim,BERNSTEIN,RealType>::func_and_args suffix                  \
      case RATIONAL_BERNSTEIN:                                          \
        prefix FE<dim,RATIONAL_BERNSTEIN,RealType>::func_and_args suffix         \
      case SZABAB:                                                      \
        prefix FE<dim,SZABAB,RealType>::func_and_args suffix                     \
      case XYZ:                                                         \
        prefix FEXYZ<dim,RealType>::func_and_args suffix                         \
      case SUBDIVISION:                                                 \
        libmesh_assert_equal_to (dim, 2);                               \
        prefix FE<2,SUBDIVISION,RealType>::func_and_args suffix                  \
      case LAGRANGE_VEC:                                                \
      case NEDELEC_ONE:                                                 \
      case MONOMIAL_VEC:                                                \
        libmesh_error_msg("Error: Can only request scalar valued elements for Real FEInterface::func_and_args"); \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)


#define fe_vector_scalar_error_switch(dim, func_and_args, prefix, suffix, RealType) \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case LAGRANGE_VEC:                                                \
        prefix FELagrangeVec<dim,RealType>::func_and_args suffix                 \
      case NEDELEC_ONE:                                                 \
        prefix FENedelecOne<dim,RealType>::func_and_args suffix                  \
      case MONOMIAL_VEC:                                                \
        prefix FEMonomialVec<dim,RealType>::func_and_args suffix                 \
      case HERMITE:                                                     \
      case HIERARCHIC:                                                  \
      case L2_HIERARCHIC:                                               \
      case LAGRANGE:                                                    \
      case L2_LAGRANGE:                                                 \
      case MONOMIAL:                                                    \
      case SCALAR:                                                      \
      case BERNSTEIN:                                                   \
      case SZABAB:                                                      \
      case RATIONAL_BERNSTEIN:                                          \
      case XYZ:                                                         \
      case SUBDIVISION:                                                 \
        libmesh_error_msg("Error: Can only request vector valued elements for RealGradient FEInterface::shape"); \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)

#else
#define fe_family_switch(dim, func_and_args, prefix, suffix, RealType)           \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case CLOUGH:                                                      \
        prefix FE<dim,CLOUGH,RealType>::func_and_args suffix                     \
      case HERMITE:                                                     \
        prefix FE<dim,HERMITE,RealType>::func_and_args suffix                    \
      case HIERARCHIC:                                                  \
        prefix FE<dim,HIERARCHIC,RealType>::func_and_args suffix                 \
      case L2_HIERARCHIC:                                               \
        prefix FE<dim,L2_HIERARCHIC,RealType>::func_and_args suffix              \
      case LAGRANGE:                                                    \
        prefix FE<dim,LAGRANGE,RealType>::func_and_args suffix                   \
      case L2_LAGRANGE:                                                 \
        prefix FE<dim,L2_LAGRANGE,RealType>::func_and_args suffix                \
      case MONOMIAL:                                                    \
        prefix FE<dim,MONOMIAL,RealType>::func_and_args suffix                   \
      case SCALAR:                                                      \
        prefix FE<dim,SCALAR,RealType>::func_and_args suffix                     \
      case XYZ:                                                         \
        prefix FEXYZ<dim,RealType>::func_and_args suffix                         \
      case SUBDIVISION:                                                 \
        libmesh_assert_equal_to (dim, 2);                               \
        prefix FE<2,SUBDIVISION,RealType>::func_and_args suffix                  \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)

#define fe_family_with_vec_switch(dim, func_and_args, prefix, suffix, RealType)  \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case CLOUGH:                                                      \
        prefix FE<dim,CLOUGH,RealType>::func_and_args suffix                     \
      case HERMITE:                                                     \
        prefix FE<dim,HERMITE,RealType>::func_and_args suffix                    \
      case HIERARCHIC:                                                  \
        prefix FE<dim,HIERARCHIC,RealType>::func_and_args suffix                 \
      case L2_HIERARCHIC:                                               \
        prefix FE<dim,L2_HIERARCHIC,RealType>::func_and_args suffix              \
      case LAGRANGE:                                                    \
        prefix FE<dim,LAGRANGE,RealType>::func_and_args suffix                   \
      case LAGRANGE_VEC:                                                \
        prefix FELagrangeVec<dim,RealType>::func_and_args suffix                 \
      case L2_LAGRANGE:                                                 \
        prefix FE<dim,L2_LAGRANGE,RealType>::func_and_args suffix                \
      case MONOMIAL:                                                    \
        prefix FE<dim,MONOMIAL,RealType>::func_and_args suffix                   \
      case MONOMIAL_VEC:                                                \
        prefix FEMonomialVec<dim,RealType>::func_and_args suffix                 \
      case SCALAR:                                                      \
        prefix FE<dim,SCALAR,RealType>::func_and_args suffix                     \
      case XYZ:                                                         \
        prefix FEXYZ<dim,RealType>::func_and_args suffix                         \
      case SUBDIVISION:                                                 \
        libmesh_assert_equal_to (dim, 2);                               \
        prefix FE<2,SUBDIVISION,RealType>::func_and_args suffix                  \
      case NEDELEC_ONE:                                                 \
        prefix FENedelecOne<dim,RealType>::func_and_args suffix                  \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)

#define fe_scalar_vec_error_switch(dim, func_and_args, prefix, suffix, RealType) \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case CLOUGH:                                                      \
        prefix  FE<dim,CLOUGH,RealType>::func_and_args suffix                    \
      case HERMITE:                                                     \
        prefix  FE<dim,HERMITE,RealType>::func_and_args suffix                   \
      case HIERARCHIC:                                                  \
        prefix  FE<dim,HIERARCHIC,RealType>::func_and_args suffix                \
      case L2_HIERARCHIC:                                               \
        prefix  FE<dim,L2_HIERARCHIC,RealType>::func_and_args suffix             \
      case LAGRANGE:                                                    \
        prefix  FE<dim,LAGRANGE,RealType>::func_and_args suffix                  \
      case L2_LAGRANGE:                                                 \
        prefix  FE<dim,L2_LAGRANGE,RealType>::func_and_args suffix               \
      case MONOMIAL:                                                    \
        prefix  FE<dim,MONOMIAL,RealType>::func_and_args suffix                  \
      case SCALAR:                                                      \
        prefix  FE<dim,SCALAR,RealType>::func_and_args suffix                    \
      case XYZ:                                                         \
        prefix  FEXYZ<dim,RealType>::func_and_args suffix                        \
      case SUBDIVISION:                                                 \
        libmesh_assert_equal_to (dim, 2);                               \
        prefix  FE<2,SUBDIVISION,RealType>::func_and_args suffix                 \
      case LAGRANGE_VEC:                                                \
      case NEDELEC_ONE:                                                 \
      case MONOMIAL_VEC:                                                \
        libmesh_error_msg("Error: Can only request scalar valued elements for Real FEInterface::func_and_args"); \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)


#define fe_vector_scalar_error_switch(dim, func_and_args, prefix, suffix, RealType) \
  do {                                                                  \
    switch (fe_t.family)                                                \
      {                                                                 \
      case LAGRANGE_VEC:                                                \
        prefix FELagrangeVec<dim,RealType>::func_and_args suffix                 \
      case NEDELEC_ONE:                                                 \
        prefix FENedelecOne<dim,RealType>::func_and_args suffix                  \
      case MONOMIAL_VEC:                                                \
        prefix FEMonomialVec<dim,RealType>::func_and_args suffix                 \
      case HERMITE:                                                     \
      case HIERARCHIC:                                                  \
      case L2_HIERARCHIC:                                               \
      case LAGRANGE:                                                    \
      case L2_LAGRANGE:                                                 \
      case MONOMIAL:                                                    \
      case SCALAR:                                                      \
      case XYZ:                                                         \
      case SUBDIVISION:                                                 \
        libmesh_error_msg("Error: Can only request vector valued elements for RealGradient FEInterface::func_and_args"); \
      default:                                                          \
        libmesh_error_msg("Unsupported family = " << fe_t.family);      \
      }                                                                 \
  } while (0)
#endif


#define fe_switch(func_and_args, RealType)                       \
  do {                                                  \
    switch (dim)                                        \
      {                                                 \
        /* 0D */                                        \
      case 0:                                           \
        fe_family_switch (0, func_and_args, return, ;, RealType);        \
        /* 1D */                                        \
      case 1:                                           \
        fe_family_switch (1, func_and_args, return, ;, RealType);        \
        /* 2D */                                        \
      case 2:                                           \
        fe_family_switch (2, func_and_args, return, ;, RealType);        \
        /* 3D */                                        \
      case 3:                                           \
        fe_family_switch (3, func_and_args, return, ;, RealType);        \
      default:                                          \
        libmesh_error_msg("Invalid dim = " << dim);     \
      }                                                 \
  } while (0)

#define fe_with_vec_switch(func_and_args, RealType)                              \
  do {                                                                  \
    switch (dim)                                                        \
      {                                                                 \
        /* 0D */                                                        \
      case 0:                                                           \
        fe_family_with_vec_switch (0, func_and_args, return, ;, RealType);       \
        /* 1D */                                                        \
      case 1:                                                           \
        fe_family_with_vec_switch (1, func_and_args, return, ;, RealType);       \
        /* 2D */                                                        \
      case 2:                                                           \
        fe_family_with_vec_switch (2, func_and_args, return, ;, RealType);       \
        /* 3D */                                                        \
      case 3:                                                           \
        fe_family_with_vec_switch (3, func_and_args, return, ;, RealType);       \
      default:                                                          \
        libmesh_error_msg("Invalid dim = " << dim);                     \
      }                                                                 \
  } while (0)


#define void_fe_switch(func_and_args, RealType)                          \
  do {                                                          \
    switch (dim)                                                \
      {                                                         \
        /* 0D */                                                \
      case 0:                                                   \
        fe_family_switch (0, func_and_args, ;, ; return;, RealType);      \
        /* 1D */                                                \
      case 1:                                                   \
        fe_family_switch (1, func_and_args, ;, ; return;, RealType);      \
        /* 2D */                                                \
      case 2:                                                   \
        fe_family_switch (2, func_and_args, ;, ; return;, RealType);      \
        /* 3D */                                                \
      case 3:                                                   \
        fe_family_switch (3, func_and_args, ;, ; return;, RealType);      \
      default:                                                  \
        libmesh_error_msg("Invalid dim = " << dim);             \
      }                                                         \
  } while (0)

#define void_fe_with_vec_switch(func_and_args, RealType)                          \
  do {                                                                  \
    switch (dim)                                                        \
      {                                                                 \
        /* 0D */                                                        \
      case 0:                                                           \
        fe_family_with_vec_switch (0, func_and_args, ;, ; return;, RealType);     \
        /* 1D */                                                        \
      case 1:                                                           \
        fe_family_with_vec_switch (1, func_and_args, ;, ; return;, RealType);     \
        /* 2D */                                                        \
      case 2:                                                           \
        fe_family_with_vec_switch (2, func_and_args, ;, ; return;, RealType);     \
        /* 3D */                                                        \
      case 3:                                                           \
        fe_family_with_vec_switch (3, func_and_args, ;, ; return;, RealType);     \
      default:                                                          \
        libmesh_error_msg("Invalid dim = " << dim);                     \
      }                                                                 \
  } while (0)


template <typename RealType>
unsigned int
FEInterface::n_dofs (const unsigned int dim,
                     const FEType & fe_t,
                     const ElemTempl<RealType> * elem)
{
  FEType p_refined_fe_t = fe_t;
  p_refined_fe_t.order = static_cast<Order>(p_refined_fe_t.order + elem->p_level());
  return FEInterface::n_dofs(dim, p_refined_fe_t, elem->type());
}


template <typename RealType>
void FEInterface::dofs_on_side(const ElemTempl<RealType> * const elem,
                               const unsigned int dim,
                               const FEType & fe_t,
                               unsigned int s,
                               std::vector<unsigned int> & di)
{
  const Order o = fe_t.order;

  void_fe_with_vec_switch(dofs_on_side(elem, o, s, di), RealType);
}



template <typename RealType>
void FEInterface::dofs_on_edge(const ElemTempl<RealType> * const elem,
                               const unsigned int dim,
                               const FEType & fe_t,
                               unsigned int e,
                               std::vector<unsigned int> & di)
{
  const Order o = fe_t.order;

  void_fe_with_vec_switch(dofs_on_edge(elem, o, e, di), RealType);
}




template <typename RealType>
void FEInterface::nodal_soln(const unsigned int dim,
                             const FEType & fe_t,
                             const ElemTempl<RealType> * elem,
                             const std::vector<Number> & elem_soln,
                             std::vector<Number> &       nodal_soln)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(elem->type()))
    {
      ifem_nodal_soln(dim, fe_t, elem, elem_soln, nodal_soln);
      return;
    }

#endif

  const Order order = fe_t.order;

  void_fe_with_vec_switch(nodal_soln(elem, order, elem_soln, nodal_soln), RealType);
}




template <typename RealType>
PointTempl<RealType> FEInterface::map(unsigned int dim,
                       const FEType & fe_t,
                       const ElemTempl<RealType> * elem,
                       const PointTempl<RealType> & p)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (is_InfFE_elem(elem->type()))
    return ifem_map(dim, fe_t, elem, p);
#endif
  fe_with_vec_switch(map(elem, p), RealType);
}





template <typename RealType>
PointTempl<RealType> FEInterface::inverse_map (const unsigned int dim,
                                const FEType & fe_t,
                                const ElemTempl<RealType> * elem,
                                const PointTempl<RealType> & p,
                                const Real tolerance,
                                const bool secure)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(elem->type()))
    return ifem_inverse_map(dim, fe_t, elem, p,tolerance, secure);

#endif

  fe_with_vec_switch(inverse_map(elem, p, tolerance, secure), RealType);
}




template <typename RealType>
void FEInterface::inverse_map (const unsigned int dim,
                               const FEType & fe_t,
                               const ElemTempl<RealType> * elem,
                               const std::vector<PointTempl<RealType>> & physical_points,
                               std::vector<PointTempl<RealType>> &       reference_points,
                               const Real tolerance,
                               const bool secure)
{
  const std::size_t n_pts = physical_points.size();

  // Resize the vector
  reference_points.resize(n_pts);

  if (n_pts == 0)
    {
      libMesh::err << "WARNING: empty vector physical_points!"
                   << std::endl;
      libmesh_here();
      return;
    }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(elem->type()))
    {
      ifem_inverse_map(dim, fe_t, elem, physical_points, reference_points, tolerance, secure);
      return;
      // libmesh_not_implemented();
    }

#endif

  void_fe_with_vec_switch(inverse_map(elem, physical_points, reference_points, tolerance, secure), RealType);
}



template <typename RealType>
bool FEInterface::on_reference_element(const PointTempl<RealType> & p,
                                       const ElemType t,
                                       const Real eps)
{
  return FEBase::on_reference_element(p,t,eps);
}




template <typename RealType>
RealType FEInterface::shape(const unsigned int dim,
                        const FEType & fe_t,
                        const ElemType t,
                        const unsigned int i,
                        const PointTempl<RealType> & p)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(t))
    return ifem_shape(dim, fe_t, t, i, p);

#endif

  const Order o = fe_t.order;

  fe_switch(shape(t,o,i,p), RealType);
}

template <typename RealType>
RealType FEInterface::shape(const unsigned int dim,
                        const FEType & fe_t,
                        const ElemTempl<RealType> * elem,
                        const unsigned int i,
                        const PointTempl<RealType> & p)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (elem && is_InfFE_elem(elem->type()))
    return ifem_shape(dim, fe_t, elem, i, p);

#endif

  const Order o = fe_t.order;

  fe_switch(shape(elem,o,i,p), RealType);
}

template <typename RealType>
void FEInterfaceShim<RealType, RealType>::shape(const unsigned int dim,
                                                const FEType & fe_t,
                                                const ElemType t,
                                                const unsigned int i,
                                                const PointTempl<RealType> & p,
                                                RealType & phi)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(t))
    {
      phi = ifem_shape(dim, fe_t, t, i, p);
      return;
    }

#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_scalar_vec_error_switch(0, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    case 1:
      fe_scalar_vec_error_switch(1, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    case 2:
      fe_scalar_vec_error_switch(2, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    case 3:
      fe_scalar_vec_error_switch(3, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }

  return;
}

template <typename RealType>
void FEInterfaceShim<RealType, RealType>::shape(const unsigned int dim,
                                                const FEType & fe_t,
                                                const ElemTempl<RealType> * elem,
                                                const unsigned int i,
                                                const PointTempl<RealType> & p,
                                                RealType & phi)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (elem && is_InfFE_elem(elem->type()))
    {
      phi = ifem_shape(dim, fe_t, elem, i, p);
      return;
    }
#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_scalar_vec_error_switch(0, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    case 1:
      fe_scalar_vec_error_switch(1, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    case 2:
      fe_scalar_vec_error_switch(2, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    case 3:
      fe_scalar_vec_error_switch(3, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }

  return;
}

template <typename RealType, template <typename> class Wrapper>
void FEInterfaceShim<RealType, Wrapper<RealType>>::shape(const unsigned int dim,
                                                         const FEType & fe_t,
                                                         const ElemType t,
                                                         const unsigned int i,
                                                         const PointTempl<RealType> & p,
                                                         Wrapper<RealType> & phi)
{
  // This even does not handle infinite elements at all!?
  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_vector_scalar_error_switch(0, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    case 1:
      fe_vector_scalar_error_switch(1, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    case 2:
      fe_vector_scalar_error_switch(2, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    case 3:
      fe_vector_scalar_error_switch(3, shape(t,o,i,p), phi = , ; break;, RealType);
      break;
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }

  return;
}

template <typename RealType, template <typename> class Wrapper>
void FEInterfaceShim<RealType, Wrapper<RealType>>::shape(const unsigned int dim,
                                                         const FEType & fe_t,
                                                         const ElemTempl<RealType> * elem,
                                                         const unsigned int i,
                                                         const PointTempl<RealType> & p,
                                                         Wrapper<RealType> & phi)

{
  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_vector_scalar_error_switch(0, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    case 1:
      fe_vector_scalar_error_switch(1, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    case 2:
      fe_vector_scalar_error_switch(2, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    case 3:
      fe_vector_scalar_error_switch(3, shape(elem,o,i,p), phi = , ; break;, RealType);
      break;
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }

  return;
}

template <typename RealType, typename OutputType>
void FEInterface::shape(const unsigned int dim,
                        const FEType & fe_t,
                        const ElemType t,
                        const unsigned int i,
                        const PointTempl<RealType> & p,
                        OutputType & phi)
{
  FEInterfaceShim<RealType, OutputType>::shape(dim, fe_t, t, i, p, phi);
}

template <typename RealType, typename OutputType>
void FEInterface::shape(const unsigned int dim,
                        const FEType & fe_t,
                        const ElemTempl<RealType> * elem,
                        const unsigned int i,
                        const PointTempl<RealType> & p,
                        OutputType & phi)
{
  FEInterfaceShim<RealType, OutputType>::shape(dim, fe_t, elem, i, p, phi);
}


template <typename RealType>
FEInterface::shape_ptr<RealType>
FEInterface::shape_function(const unsigned int dim,
                            const FEType & fe_t)
{
  fe_switch(shape, RealType);
}


template <typename RealType>
RealType FEInterface::shape_deriv(const unsigned int dim,
                              const FEType & fe_t,
                              const ElemType t,
                              const unsigned int i,
                              const unsigned int j,
                              const PointTempl<RealType> & p)
{
  libmesh_assert_greater (dim,j);
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(t)){
    return ifem_shape_deriv(dim, fe_t, t, i, j, p);
  }

#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_family_switch (0, shape_deriv(t, o, i, j, p), return , ;, RealType);
    case 1:
      fe_family_switch (1, shape_deriv(t, o, i, j, p), return , ;, RealType);
    case 2:
      fe_family_switch (2, shape_deriv(t, o, i, j, p), return  , ;, RealType);
    case 3:
      fe_family_switch (3, shape_deriv(t, o, i, j, p), return , ;, RealType);
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }
  return 0;
}


template <typename RealType>
RealType FEInterface::shape_deriv(const unsigned int dim,
                              const FEType & fe_t,
                              const ElemTempl<RealType> * elem,
                              const unsigned int i,
                              const unsigned int j,
                              const PointTempl<RealType> & p)
{
  libmesh_assert_greater (dim,j);
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (elem->infinite()){
    return ifem_shape_deriv(dim, fe_t, elem, i, j, p);
  }

#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_family_switch (0, shape_deriv(elem, o, i, j, p), return , ;, RealType);
    case 1:
      fe_family_switch (1, shape_deriv(elem, o, i, j, p), return , ;, RealType);
    case 2:
      fe_family_switch (2, shape_deriv(elem, o, i, j, p), return , ;, RealType);
    case 3:
      fe_family_switch (3, shape_deriv(elem, o, i, j, p), return , ;, RealType);
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }
  return 0;
}


template <typename RealType>
FEInterface::shape_deriv_ptr<RealType>
FEInterface::shape_deriv_function(const unsigned int dim,
                                  const FEType & fe_t)
{
  fe_switch(shape_deriv, RealType);
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <typename RealType>
RealType FEInterface::shape_second_deriv(const unsigned int dim,
                                     const FEType & fe_t,
                                     const ElemType t,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const PointTempl<RealType> & p)
{
  libmesh_assert_greater_equal (dim*(dim-1),j);
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (is_InfFE_elem(t))
    libmesh_not_implemented();
#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_family_switch (0, shape_second_deriv(t, o, i, j, p), return , ;, RealType);
    case 1:
      fe_family_switch (1, shape_second_deriv(t, o, i, j, p), return , ;, RealType);
    case 2:
      fe_family_switch (2, shape_second_deriv(t, o, i, j, p), return  , ;, RealType);
    case 3:
      fe_family_switch (3, shape_second_deriv(t, o, i, j, p), return , ;, RealType);
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }
  return 0;
}


template <typename RealType>
RealType FEInterface::shape_second_deriv(const unsigned int dim,
                                     const FEType & fe_t,
                                     const ElemTempl<RealType> * elem,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const PointTempl<RealType> & p)
{
  libmesh_assert_greater_equal (dim,j);
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    libmesh_not_implemented();
#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_family_switch (0, shape_second_deriv(elem, o, i, j, p), return , ;, RealType);
    case 1:
      fe_family_switch (1, shape_second_deriv(elem, o, i, j, p), return , ;, RealType);
    case 2:
      fe_family_switch (2, shape_second_deriv(elem, o, i, j, p), return , ;, RealType);
    case 3:
      fe_family_switch (3, shape_second_deriv(elem, o, i, j, p), return , ;, RealType);
    default:
      libmesh_error_msg("Invalid dimension = " << dim);
    }
  return 0;
}


template <typename RealType>
FEInterface::shape_second_deriv_ptr<RealType>
FEInterface::shape_second_deriv_function(const unsigned int dim,
                                         const FEType & fe_t)
{
  fe_switch(shape_second_deriv, RealType);
}
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES


template <typename RealType>
unsigned int FEInterface::n_vec_dim (const MeshBaseTempl<RealType> & mesh,
                                     const FEType & fe_type)
{
  switch (fe_type.family)
    {
      //FIXME: We currently assume that the number of vector components is tied
      //       to the mesh dimension. This will break for mixed-dimension meshes.
    case LAGRANGE_VEC:
    case NEDELEC_ONE:
    case MONOMIAL_VEC:
      return mesh.mesh_dimension();
    default:
      return 1;
    }
}


} // namespace libMesh

#define INSTANTIATE_FE_INTERFACE_METHODS0(RealType)                                                \
  template unsigned int FEInterface::n_dofs(                                                       \
      const unsigned int dim, const FEType & fe_t, const ElemTempl<RealType> * elem);              \
  template void FEInterface::dofs_on_side(const ElemTempl<RealType> * const elem,                  \
                                          const unsigned int dim,                                  \
                                          const FEType & fe_t,                                     \
                                          unsigned int s,                                          \
                                          std::vector<unsigned int> & di);                         \
  template void FEInterface::dofs_on_edge(const ElemTempl<RealType> * const elem,                  \
                                          const unsigned int dim,                                  \
                                          const FEType & fe_t,                                     \
                                          unsigned int e,                                          \
                                          std::vector<unsigned int> & di);                         \
  template void FEInterface::nodal_soln(const unsigned int dim,                                    \
                                        const FEType & fe_t,                                       \
                                        const ElemTempl<RealType> * elem,                          \
                                        const std::vector<Number> & elem_soln,                     \
                                        std::vector<Number> & nodal_soln);                         \
  template PointTempl<RealType> FEInterface::map(unsigned int dim,                                 \
                                                 const FEType & fe_t,                              \
                                                 const ElemTempl<RealType> * elem,                 \
                                                 const PointTempl<RealType> & p);                  \
  template PointTempl<RealType> FEInterface::inverse_map(const unsigned int dim,                   \
                                                         const FEType & fe_t,                      \
                                                         const ElemTempl<RealType> * elem,         \
                                                         const PointTempl<RealType> & p,           \
                                                         const Real tolerance = TOLERANCE,         \
                                                         const bool secure = true);                \
  template void FEInterface::inverse_map(                                                          \
      const unsigned int dim,                                                                      \
      const FEType & fe_t,                                                                         \
      const ElemTempl<RealType> * elem,                                                            \
      const std::vector<PointTempl<RealType>> & physical_points,                                   \
      std::vector<PointTempl<RealType>> & reference_points,                                        \
      const Real tolerance = TOLERANCE,                                                            \
      const bool secure = true);                                                                   \
  template bool FEInterface::on_reference_element(                                                 \
      const PointTempl<RealType> & p, const ElemType t, const Real eps = TOLERANCE);               \
  template RealType FEInterface::shape(const unsigned int dim,                                     \
                                       const FEType & fe_t,                                        \
                                       const ElemType t,                                           \
                                       const unsigned int i,                                       \
                                       const PointTempl<RealType> & p);                            \
  template RealType FEInterface::shape(const unsigned int dim,                                     \
                                       const FEType & fe_t,                                        \
                                       const ElemTempl<RealType> * elem,                           \
                                       const unsigned int i,                                       \
                                       const PointTempl<RealType> & p);                            \
  template void FEInterface::shape(const unsigned int dim,                                         \
                                   const FEType & fe_t,                                            \
                                   const ElemType t,                                               \
                                   const unsigned int i,                                           \
                                   const PointTempl<RealType> & p,                                 \
                                   RealType & phi);                                                \
  template void FEInterface::shape(const unsigned int dim,                                         \
                                   const FEType & fe_t,                                            \
                                   const ElemType t,                                               \
                                   const unsigned int i,                                           \
                                   const PointTempl<RealType> & p,                                 \
                                   VectorValue<RealType> & phi);                                   \
  template void FEInterface::shape(const unsigned int dim,                                         \
                                   const FEType & fe_t,                                            \
                                   const ElemTempl<RealType> * elem,                               \
                                   const unsigned int i,                                           \
                                   const PointTempl<RealType> & p,                                 \
                                   RealType & phi);                                                \
  template void FEInterface::shape(const unsigned int dim,                                         \
                                   const FEType & fe_t,                                            \
                                   const ElemTempl<RealType> * elem,                               \
                                   const unsigned int i,                                           \
                                   const PointTempl<RealType> & p,                                 \
                                   VectorValue<RealType> & phi);                                   \
  template FEInterface::shape_ptr<RealType> FEInterface::shape_function<RealType>(                 \
      const unsigned int dim, const FEType & fe_t);                                                \
  template RealType FEInterface::shape_deriv(const unsigned int dim,                               \
                                             const FEType & fe_t,                                  \
                                             const ElemType t,                                     \
                                             const unsigned int i,                                 \
                                             const unsigned int j,                                 \
                                             const PointTempl<RealType> & p);                      \
  template RealType FEInterface::shape_deriv(const unsigned int dim,                               \
                                             const FEType & fe_t,                                  \
                                             const ElemTempl<RealType> * elem,                     \
                                             const unsigned int i,                                 \
                                             const unsigned int j,                                 \
                                             const PointTempl<RealType> & p);                      \
  template FEInterface::shape_deriv_ptr<RealType> FEInterface::shape_deriv_function<RealType>(     \
      const unsigned int dim, const FEType & fe_t);                                                \
  template unsigned int FEInterface::n_vec_dim(const MeshBaseTempl<RealType> & mesh,               \
                                               const FEType & fe_type)

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#define INSTANTIATE_FE_INTERFACE_METHODS1(RealType)                                                \
  INSTANTIATE_FE_INTERFACE_METHODS0(RealType);                                                     \
  template RealType FEInterface::shape_second_deriv(const unsigned int dim,                        \
                                                    const FEType & fe_t,                           \
                                                    const ElemType t,                              \
                                                    const unsigned int i,                          \
                                                    const unsigned int j,                          \
                                                    const PointTempl<RealType> & p);               \
  template RealType FEInterface::shape_second_deriv(const unsigned int dim,                        \
                                                    const FEType & fe_t,                           \
                                                    const ElemTempl<RealType> * elem,              \
                                                    const unsigned int i,                          \
                                                    const unsigned int j,                          \
                                                    const PointTempl<RealType> & p);               \
  template FEInterface::shape_second_deriv_ptr<RealType>                                           \
  FEInterface::shape_second_deriv_function<RealType>(const unsigned int dim, const FEType & fe_t)
#else
#define INSTANTIATE_FE_INTERFACE_METHODS1(RealType) INSTANTIATE_FE_INTERFACE_METHODS0(RealType)
#endif

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#define INSTANTIATE_FE_INTERFACE_METHODS(RealType)                                                 \
  INSTANTIATE_FE_INTERFACE_METHODS1(RealType);                                                     \
  template void FEInterface::ifem_nodal_soln(const unsigned int dim,                               \
                                             const FEType & fe_t,                                  \
                                             const ElemTempl<RealType> * elem,                     \
                                             const std::vector<Number> & elem_soln,                \
                                             std::vector<Number> & nodal_soln);                    \
  template PointTempl<RealType> FEInterface::ifem_map(const unsigned int dim,                      \
                                                      const FEType & fe_t,                         \
                                                      const ElemTempl<RealType> * elem,            \
                                                      const PointTempl<RealType> & p);             \
  template PointTempl<RealType> FEInterface::ifem_inverse_map(const unsigned int dim,              \
                                                              const FEType & fe_t,                 \
                                                              const ElemTempl<RealType> * elem,    \
                                                              const PointTempl<RealType> & p,      \
                                                              const Real tolerance = TOLERANCE,    \
                                                              const bool secure = true);           \
  template void FEInterface::ifem_inverse_map(                                                     \
      const unsigned int dim,                                                                      \
      const FEType & fe_t,                                                                         \
      const ElemTempl<RealType> * elem,                                                            \
      const std::vector<PointTempl<RealType>> & physical_points,                                   \
      std::vector<PointTempl<RealType>> & reference_points,                                        \
      const Real tolerance = TOLERANCE,                                                            \
      const bool secure = true);                                                                   \
  template bool FEInterface::ifem_on_reference_element(                                            \
      const PointTempl<RealType> & p, const ElemType t, const Real eps);                           \
  template RealType FEInterface::ifem_shape(const unsigned int dim,                                \
                                            const FEType & fe_t,                                   \
                                            const ElemType t,                                      \
                                            const unsigned int i,                                  \
                                            const PointTempl<RealType> & p);                       \
  template RealType FEInterface::ifem_shape(const unsigned int dim,                                \
                                            const FEType & fe_t,                                   \
                                            const ElemTempl<RealType> * elem,                      \
                                            const unsigned int i,                                  \
                                            const PointTempl<RealType> & p);                       \
  template RealType FEInterface::ifem_shape_deriv(const unsigned int dim,                          \
                                                  const FEType & fe_t,                             \
                                                  const ElemType t,                                \
                                                  const unsigned int i,                            \
                                                  const unsigned int j,                            \
                                                  const PointTempl<RealType> & p);                 \
  template RealType FEInterface::ifem_shape_deriv(const unsigned int dim,                          \
                                                  const FEType & fe_t,                             \
                                                  const ElemTempl<RealType> * elem,                \
                                                  const unsigned int i,                            \
                                                  const unsigned int j,                            \
                                                  const PointTempl<RealType> & p);                 \
  template void FEInterface::ifem_compute_data(const unsigned int dim,                             \
                                               const FEType & fe_t,                                \
                                               const ElemTempl<RealType> * elem,                   \
                                               FEComputeData & data)

#else
#define INSTANTIATE_FE_INTERFACE_METHODS(RealType) INSTANTIATE_FE_INTERFACE_METHODS1(RealType)
#endif

#endif // FE_INTERFACE_IMPL_H
