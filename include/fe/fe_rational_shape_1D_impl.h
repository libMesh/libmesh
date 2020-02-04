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

#ifndef LIBMESH_FE_RATIONAL_SHAPE_1D_IMPL_H
#define LIBMESH_FE_RATIONAL_SHAPE_1D_IMPL_H

// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"


namespace {
extern const libMesh::FEFamily _underlying_fe_family;
}


namespace libMesh
{

template <typename RealType>
typename FEShim<1,RATIONAL_BERNSTEIN,RealType>::OutputShape FEShim<1,RATIONAL_BERNSTEIN,RealType>::shape(const ElemType,
                                     const Order,
                                     const unsigned int,
                                     const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,RATIONAL_BERNSTEIN,RealType>::OutputShape FEShim<1,RATIONAL_BERNSTEIN,RealType>::shape(const Elem * elem,
                                     const Order order,
                                     const unsigned int i,
                                     const Point & p,
                                     const bool add_p_level)
{
  libmesh_assert(elem);

  const ElemType elem_type = elem->type();

  const Order totalorder = static_cast<Order>(order + add_p_level * elem->p_level());

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(totalorder, _underlying_fe_family);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(1, fe_type, elem_type);

  const unsigned int n_nodes = elem->n_nodes();
  libmesh_assert_equal_to (n_sf, n_nodes);

  std::vector<Real> node_weights(n_nodes);

  const unsigned char datum_index = elem->mapping_data();
  for (unsigned int n=0; n<n_nodes; n++)
    node_weights[n] =
      elem->node_ref(n).template get_extra_datum<Real>(datum_index);

  Real weighted_shape_i = 0, weighted_sum = 0;

  for (unsigned int sf=0; sf<n_sf; sf++)
    {
      Real weighted_shape = node_weights[sf] *
        FEInterface::shape(1, fe_type, elem, sf, p);
      weighted_sum += weighted_shape;
      if (sf == i)
        weighted_shape_i = weighted_shape;
    }

  return weighted_shape_i / weighted_sum;
}



template <typename RealType>
typename FEShim<1,RATIONAL_BERNSTEIN,RealType>::OutputShape FEShim<1,RATIONAL_BERNSTEIN,RealType>::shape_deriv(const ElemType,
                                           const Order,
                                           const unsigned int,
                                           const unsigned int,
                                           const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,RATIONAL_BERNSTEIN,RealType>::OutputShape FEShim<1,RATIONAL_BERNSTEIN,RealType>::shape_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int libmesh_dbg_var(j),
                                           const Point & p,
                                           const bool add_p_level)
{
  // only d()/dxi in 1D!
  libmesh_assert_equal_to (j, 0);

  libmesh_assert(elem);

  const ElemType elem_type = elem->type();

  const Order totalorder = static_cast<Order>(order + add_p_level * elem->p_level());

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(totalorder, _underlying_fe_family);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(1, fe_type, elem_type);

  const unsigned int n_nodes = elem->n_nodes();
  libmesh_assert_equal_to (n_sf, n_nodes);

  std::vector<Real> node_weights(n_nodes);

  const unsigned char datum_index = elem->mapping_data();
  for (unsigned int n=0; n<n_nodes; n++)
    node_weights[n] =
      elem->node_ref(n).template get_extra_datum<Real>(datum_index);

  Real weighted_shape_i = 0, weighted_sum = 0,
       weighted_grad_i = 0, weighted_grad_sum = 0;

  for (unsigned int sf=0; sf<n_sf; sf++)
    {
      Real weighted_shape = node_weights[sf] *
        FEInterface::shape(1, fe_type, elem, sf, p);
      Real weighted_grad = node_weights[sf] *
        FEInterface::shape_deriv(1, fe_type, elem, sf, 0, p);
      weighted_sum += weighted_shape;
      weighted_grad_sum += weighted_grad;
      if (sf == i)
        {
          weighted_shape_i = weighted_shape;
          weighted_grad_i = weighted_grad;
        }
    }

  return (weighted_sum * weighted_grad_i - weighted_shape_i * weighted_grad_sum) /
         weighted_sum / weighted_sum;
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<1,RATIONAL_BERNSTEIN,RealType>::OutputShape FEShim<1,RATIONAL_BERNSTEIN,RealType>::shape_second_deriv(const ElemType,
                                                  const Order,
                                                  const unsigned int,
                                                  const unsigned int,
                                                  const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}




template <typename RealType>
typename FEShim<1,RATIONAL_BERNSTEIN,RealType>::OutputShape FEShim<1,RATIONAL_BERNSTEIN,RealType>::shape_second_deriv(const Elem * elem,
                                                  const Order order,
                                                  const unsigned int i,
                                                  const unsigned int libmesh_dbg_var(j),
                                                  const Point & p,
                                                  const bool add_p_level)
{
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!
  libmesh_assert_equal_to (j, 0);

  libmesh_assert(elem);

  const ElemType elem_type = elem->type();

  const Order totalorder = static_cast<Order>(order + add_p_level * elem->p_level());

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(totalorder, _underlying_fe_family);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(1, fe_type, elem_type);

  const unsigned int n_nodes = elem->n_nodes();
  libmesh_assert_equal_to (n_sf, n_nodes);

  std::vector<Real> node_weights(n_nodes);

  const unsigned char datum_index = elem->mapping_data();
  for (unsigned int n=0; n<n_nodes; n++)
    node_weights[n] =
      elem->node_ref(n).template get_extra_datum<Real>(datum_index);

  Real weighted_shape_i = 0, weighted_sum = 0,
       weighted_grad_i = 0, weighted_grad_sum = 0,
       weighted_hess_i = 0, weighted_hess_sum = 0;

  for (unsigned int sf=0; sf<n_sf; sf++)
    {
      Real weighted_shape = node_weights[sf] *
        FEInterface::shape(1, fe_type, elem, sf, p);
      Real weighted_grad = node_weights[sf] *
        FEInterface::shape_deriv(1, fe_type, elem, sf, 0, p);
      Real weighted_hess = node_weights[sf] *
        FEInterface::shape_second_deriv(1, fe_type, elem, sf, 0, p);
      weighted_sum += weighted_shape;
      weighted_grad_sum += weighted_grad;
      weighted_hess_sum += weighted_hess;
      if (sf == i)
        {
          weighted_shape_i = weighted_shape;
          weighted_grad_i = weighted_grad;
          weighted_hess_i = weighted_hess;
        }
    }

  return (weighted_sum * weighted_sum *
          (weighted_sum * weighted_hess_i - weighted_shape_i * weighted_hess_sum) -
          (weighted_sum * weighted_grad_i - weighted_shape_i * weighted_grad_sum) *
          2 * weighted_sum * weighted_grad_sum) /
         (weighted_sum * weighted_sum * weighted_sum * weighted_sum);
}

#endif

} // namespace libMesh


#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_RATIONAL_SHAPE_1D_IMPL_H
