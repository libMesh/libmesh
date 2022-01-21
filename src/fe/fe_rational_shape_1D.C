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

#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"


namespace {
using namespace libMesh;

static const FEFamily _underlying_fe_family = BERNSTEIN;

} // anonymous namespace



namespace libMesh
{


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape(const Elem * elem,
                                     const Order order,
                                     const unsigned int i,
                                     const Point & p,
                                     const bool add_p_level)
{
  libmesh_assert(elem);

  // FEType object for the non-rational basis underlying this one
  FEType fe_type(order, _underlying_fe_family);

  return rational_fe_shape(*elem, fe_type, i, p, add_p_level);
}



template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape(const ElemType,
                                     const Order,
                                     const unsigned int,
                                     const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape(const FEType fet,
                                     const Elem * elem,
                                     const unsigned int i,
                                     const Point & p,
                                     const bool add_p_level)
{
  libmesh_assert(elem);
  return FE<1,RATIONAL_BERNSTEIN>::shape(elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);

  FEType underlying_fe_type(order, _underlying_fe_family);

  return rational_fe_shape_deriv(*elem, underlying_fe_type, i, j, p,
                                 add_p_level);
}


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape_deriv(const ElemType,
                                           const Order,
                                           const unsigned int,
                                           const unsigned int,
                                           const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape_deriv(const FEType fet,
                                           const Elem * elem,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);
  return FE<1,RATIONAL_BERNSTEIN>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}




#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape_second_deriv(const Elem * elem,
                                                  const Order order,
                                                  const unsigned int i,
                                                  const unsigned int j,
                                                  const Point & p,
                                                  const bool add_p_level)
{
  libmesh_assert(elem);

  // FEType object to be passed to various FEInterface functions below.
  FEType underlying_fe_type(order, _underlying_fe_family);

  return rational_fe_shape_second_deriv(*elem, underlying_fe_type, i,
                                        j, p, add_p_level);
}


template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape_second_deriv(const ElemType,
                                                  const Order,
                                                  const unsigned int,
                                                  const unsigned int,
                                                  const Point &)
{
  libmesh_error_msg("Rational bases require the real element \nto query nodal weighting.");
  return 0.;
}




template <>
Real FE<1,RATIONAL_BERNSTEIN>::shape_second_deriv(const FEType fet,
                                                  const Elem * elem,
                                                  const unsigned int i,
                                                  const unsigned int j,
                                                  const Point & p,
                                                  const bool add_p_level)
{
  libmesh_assert(elem);
  return FE<1,RATIONAL_BERNSTEIN>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}


#endif


template<>
void FE<1,RATIONAL_BERNSTEIN>::shapes
  (const Elem * elem,
   const Order o,
   const unsigned int i,
   const std::vector<Point> & p,
   std::vector<OutputShape> & vi,
   const bool add_p_level)
{
  libmesh_assert_equal_to(p.size(), vi.size());
  for (auto j : index_range(vi))
    vi[j] = FE<1,RATIONAL_BERNSTEIN>::shape (elem, o, i, p[j], add_p_level);
}

template<>
void FE<1,RATIONAL_BERNSTEIN>::all_shapes
  (const Elem * elem,
   const Order o,
   const std::vector<Point> & p,
   std::vector<std::vector<OutputShape>> & v,
   const bool add_p_level)
{
  std::vector<std::vector<Real>> shapes;

  FEType fe_type(o, _underlying_fe_family);

  rational_fe_weighted_shapes(elem, fe_type, shapes, p, add_p_level);

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
          Real old_shape = FE<1,RATIONAL_BERNSTEIN>::shape (elem, o, i, p[j], add_p_level);
          libmesh_assert(std::abs(v[i][j] - old_shape) < TOLERANCE*TOLERANCE);
#endif
        }
    }
}


template<>
void FE<1,RATIONAL_BERNSTEIN>::shape_derivs
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
    v[vi] = FE<1,RATIONAL_BERNSTEIN>::shape_deriv (elem, o, i, j, p[vi], add_p_level);
}


template <>
void
FE<1,RATIONAL_BERNSTEIN>::all_shape_derivs (const Elem * elem,
                                            const Order o,
                                            const std::vector<Point> & p,
                                            const bool add_p_level)
{
  std::vector<std::vector<OutputShape>> * comps[3]
    { &this->dphidxi, &this->dphideta, &this->dphidzeta };
  for (unsigned int d=0; d != 1; ++d)
    {
      libmesh_assert_equal_to(comps[d]->size(), elem->n_nodes());
      auto & comps_d = *comps[d];
      for (auto i : index_range(comps_d))
        FE<1,RATIONAL_BERNSTEIN>::shape_derivs
          (elem,o,i,d,p,comps_d[i],add_p_level);
    }
}


} // namespace libMesh


#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
