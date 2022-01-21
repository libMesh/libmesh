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
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"


namespace {
using namespace libMesh;

static const FEFamily _underlying_fe_family = BERNSTEIN;


// shapes[i][q] is shape function phi_i at point p[q]
// derivs[j][i][q] is dphi_i/dxi_j at p[q]
void weighted_shapes_derivs(const Elem * elem,
                            FEType fe_type,
                            std::vector<std::vector<Real>> & shapes,
                            std::vector<std::vector<std::vector<Real>>> & derivs,
                            const std::vector<Point> & p,
                            const bool add_p_level)
{
  const int extra_order = add_p_level * elem->p_level();
  constexpr int dim = 3;

  const unsigned int n_sf =
    FEInterface::n_shape_functions(fe_type, extra_order, elem);

  libmesh_assert_equal_to (n_sf, elem->n_nodes());

  libmesh_assert_equal_to (dim, derivs.size());
  for (unsigned int d = 0; d != dim; ++d)
    derivs[d].resize(n_sf);

  std::vector<Real> node_weights(n_sf);

  const unsigned char datum_index = elem->mapping_data();
  for (unsigned int n=0; n<n_sf; n++)
    node_weights[n] =
      elem->node_ref(n).get_extra_datum<Real>(datum_index);

  const std::size_t n_p = p.size();

  shapes.resize(n_sf);
  for (unsigned int i=0; i != n_sf; ++i)
    {
      auto & shapes_i = shapes[i];

      shapes_i.resize(n_p, 0);

      FEInterface::shapes(dim, fe_type, elem, i, p, shapes_i, add_p_level);
      for (auto & s : shapes_i)
        s *= node_weights[i];

      for (unsigned int d = 0; d != dim; ++d)
        {
          auto & derivs_di = derivs[d][i];
          derivs_di.resize(n_p);
          FEInterface::shape_derivs(fe_type, elem, i, d, p,
                                    derivs_di, add_p_level);
          for (auto & dip : derivs_di)
            dip *= node_weights[i];
        }
    }
}


} // anonymous namespace



namespace libMesh
{


template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape(const Elem * elem,
                                     const Order order,
                                     const unsigned int i,
                                     const Point & p,
                                     const bool add_p_level)
{
  // FEType object for the non-rational basis underlying this one
  FEType fe_type(order, _underlying_fe_family);

  return rational_fe_shape(elem, fe_type, i, p, add_p_level);
}



template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape(const ElemType,
                                     const Order,
                                     const unsigned int,
                                     const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}

template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape(const FEType fet,
                                     const Elem * elem,
                                     const unsigned int i,
                                     const Point & p,
                                     const bool add_p_level)
{
  return FE<3,RATIONAL_BERNSTEIN>::shape
    (elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);

  int extra_order = add_p_level * elem->p_level();

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(order, _underlying_fe_family);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(fe_type, extra_order, elem);

  const unsigned int n_nodes = elem->n_nodes();
  libmesh_assert_equal_to (n_sf, n_nodes);

  std::vector<Real> node_weights(n_nodes);

  const unsigned char datum_index = elem->mapping_data();
  for (unsigned int n=0; n<n_nodes; n++)
    node_weights[n] =
      elem->node_ref(n).get_extra_datum<Real>(datum_index);

  Real weighted_shape_i = 0, weighted_sum = 0,
       weighted_grad_i = 0, weighted_grad_sum = 0;

  for (unsigned int sf=0; sf<n_sf; sf++)
    {
      Real weighted_shape = node_weights[sf] *
        FEInterface::shape(fe_type, extra_order, elem, sf, p);
      Real weighted_grad = node_weights[sf] *
        FEInterface::shape_deriv(fe_type, extra_order, elem, sf, j, p);
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


template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape_deriv(const ElemType,
                                           const Order,
                                           const unsigned int,
                                           const unsigned int,
                                           const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}


template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape_deriv(const FEType fet,
                                           const Elem * elem,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  return FE<3,RATIONAL_BERNSTEIN>::shape_deriv
    (elem, fet.order, i, j, p, add_p_level);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape_second_deriv(const Elem * elem,
                                                  const Order order,
                                                  const unsigned int i,
                                                  const unsigned int j,
                                                  const Point & p,
                                                  const bool add_p_level)
{
  unsigned int j1, j2;
  switch (j)
    {
    case 0:
      // j = 0 ==> d^2 phi / dxi^2
      j1 = j2 = 0;
      break;
    case 1:
      // j = 1 ==> d^2 phi / dxi deta
      j1 = 0;
      j2 = 1;
      break;
    case 2:
      // j = 2 ==> d^2 phi / deta^2
      j1 = j2 = 1;
      break;
    case 3:
      // j = 3 ==> d^2 phi / dxi dzeta
      j1 = 0;
      j2 = 2;
      break;
    case 4:
      // j = 4 ==> d^2 phi / deta dzeta
      j1 = 1;
      j2 = 2;
      break;
    case 5:
      // j = 5 ==> d^2 phi / dzeta^2
      j1 = j2 = 2;
      break;
    default:
      libmesh_error();
    }

  libmesh_assert(elem);

  int extra_order = add_p_level * elem->p_level();

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(order, _underlying_fe_family);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(fe_type, extra_order, elem);

  const unsigned int n_nodes = elem->n_nodes();
  libmesh_assert_equal_to (n_sf, n_nodes);

  std::vector<Real> node_weights(n_nodes);

  const unsigned char datum_index = elem->mapping_data();
  for (unsigned int n=0; n<n_nodes; n++)
    node_weights[n] =
      elem->node_ref(n).get_extra_datum<Real>(datum_index);

  Real weighted_shape_i = 0, weighted_sum = 0,
       weighted_grada_i = 0, weighted_grada_sum = 0,
       weighted_gradb_i = 0, weighted_gradb_sum = 0,
       weighted_hess_i = 0, weighted_hess_sum = 0;

  for (unsigned int sf=0; sf<n_sf; sf++)
    {
      Real weighted_shape = node_weights[sf] *
        FEInterface::shape(fe_type, extra_order, elem, sf, p);
      Real weighted_grada = node_weights[sf] *
        FEInterface::shape_deriv(fe_type, extra_order, elem, sf, j1, p);
      Real weighted_hess = node_weights[sf] *
        FEInterface::shape_second_deriv(fe_type, extra_order, elem, sf, j, p);
      weighted_sum += weighted_shape;
      weighted_grada_sum += weighted_grada;
      Real weighted_gradb = weighted_grada;
      if (j1 != j2)
        {
          weighted_gradb = (j1 == j2) ? weighted_grada :
            node_weights[sf] *
            FEInterface::shape_deriv(fe_type, extra_order, elem, sf, j2, p);
          weighted_grada_sum += weighted_grada;
        }
      weighted_hess_sum += weighted_hess;
      if (sf == i)
        {
          weighted_shape_i = weighted_shape;
          weighted_grada_i = weighted_grada;
          weighted_gradb_i = weighted_gradb;
          weighted_hess_i = weighted_hess;
        }
    }

  if (j1 == j2)
    weighted_gradb_sum = weighted_grada_sum;

  return (weighted_sum * weighted_hess_i - weighted_grada_i * weighted_gradb_sum -
          weighted_shape_i * weighted_hess_sum - weighted_gradb_i * weighted_grada_sum +
          2 * weighted_grada_sum * weighted_shape_i * weighted_gradb_sum / weighted_sum) /
         weighted_sum / weighted_sum;
}



template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape_second_deriv(const ElemType,
                                                  const Order,
                                                  const unsigned int,
                                                  const unsigned int,
                                                  const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}


template <>
Real FE<3,RATIONAL_BERNSTEIN>::shape_second_deriv(const FEType fet,
                                                  const Elem * elem,
                                                  const unsigned int i,
                                                  const unsigned int j,
                                                  const Point & p,
                                                  const bool add_p_level)
{
  return FE<3,RATIONAL_BERNSTEIN>::shape_second_deriv
    (elem, fet.order, i, j, p, add_p_level);
}


#endif


template<>
void FE<3,RATIONAL_BERNSTEIN>::shapes
  (const Elem * elem,
   const Order o,
   const unsigned int i,
   const std::vector<Point> & p,
   std::vector<OutputShape> & vi,
   const bool add_p_level)
{
  libmesh_assert_equal_to(p.size(), vi.size());
  for (auto j : index_range(vi))
    vi[j] = FE<3,RATIONAL_BERNSTEIN>::shape (elem, o, i, p[j], add_p_level);
}

template<>
void FE<3,RATIONAL_BERNSTEIN>::all_shapes
  (const Elem * elem,
   const Order o,
   const std::vector<Point> & p,
   std::vector<std::vector<OutputShape>> & v,
   const bool add_p_level)
{
  std::vector<std::vector<Real>> shapes;

  FEType underlying_fe_type(o, _underlying_fe_family);

  rational_fe_weighted_shapes(elem, underlying_fe_type, shapes, p,
                              add_p_level);

  std::vector<Real> shape_sums(p.size(), 0);

  for (auto i : index_range(v))
    {
      libmesh_assert_equal_to ( p.size(), shapes[i].size() );
      for (auto j : index_range(p))
        shape_sums[j] += shapes[i][j];
    }

  for (auto i : index_range(v))
    {
      libmesh_assert_equal_to ( p.size(), v[i].size() );
      for (auto j : index_range(v[i]))
        {
          v[i][j] = shapes[i][j] / shape_sums[j];

#ifdef DEBUG
          Real old_shape = FE<3,RATIONAL_BERNSTEIN>::shape (elem, o, i, p[j], add_p_level);
          libmesh_assert(std::abs(v[i][j] - old_shape) < TOLERANCE*TOLERANCE);
#endif
        }
    }
}


template<>
void FE<3,RATIONAL_BERNSTEIN>::shape_derivs
  (const Elem * elem,
   const Order o,
   const unsigned int i,
   const unsigned int j,
   const std::vector<Point> & p,
   std::vector<OutputShape> & v,
   const bool add_p_level)
{
  libmesh_assert_equal_to(p.size(), v.size());
  for (auto vi : index_range(v))
    v[vi] = FE<3,RATIONAL_BERNSTEIN>::shape_deriv (elem, o, i, j, p[vi], add_p_level);
}


template <>
void
FE<3,RATIONAL_BERNSTEIN>::all_shape_derivs (const Elem * elem,
                                            const Order o,
                                            const std::vector<Point> & p,
                                            const bool add_p_level)
{
  constexpr int my_dim = 3;

  std::vector<std::vector<Real>> shapes;
  std::vector<std::vector<std::vector<Real>>> derivs(my_dim);
  FEType underlying_fe_type(o, _underlying_fe_family);

  weighted_shapes_derivs(elem, underlying_fe_type, shapes, derivs, p,
                         add_p_level);

  std::vector<Real> shape_sums(p.size(), 0);
  std::vector<std::vector<Real>> shape_deriv_sums(my_dim);
  for (int d=0; d != my_dim; ++d)
    shape_deriv_sums[d].resize(p.size());

  for (auto i : index_range(shapes))
    {
      libmesh_assert_equal_to ( p.size(), shapes[i].size() );
      for (auto j : index_range(p))
        shape_sums[j] += shapes[i][j];

      for (int d=0; d != my_dim; ++d)
        for (auto j : index_range(p))
          shape_deriv_sums[d][j] += derivs[d][i][j];
    }


  std::vector<std::vector<OutputShape>> * comps[my_dim]
    { &this->dphidxi, &this->dphideta, &this->dphidzeta };
  for (unsigned int d=0; d != my_dim; ++d)
    {
      auto & comps_d = *comps[d];
      libmesh_assert_equal_to(comps_d.size(), elem->n_nodes());

      for (auto i : index_range(comps_d))
        {
          auto & comps_di = comps_d[i];
          auto & derivs_di = derivs[d][i];

          for (auto j : index_range(comps_di))
            {
              comps_di[j] = (shape_sums[j] * derivs_di[j] -
                shapes[i][j] * shape_deriv_sums[d][j]) /
                shape_sums[j] / shape_sums[j];

#ifdef DEBUG
              Real old_deriv = FE<3,RATIONAL_BERNSTEIN>::shape_deriv (elem, o, i, d, p[j], add_p_level);
              libmesh_assert(std::abs(comps_di[j] - old_deriv) < TOLERANCE*TOLERANCE);
#endif
            }
        }
    }
}


} // namespace libMesh


#endif// LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
