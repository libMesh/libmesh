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

#ifndef LIBMESH_FE_HERMITE_SHAPE_1D_IMPL_H
#define LIBMESH_FE_HERMITE_SHAPE_1D_IMPL_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/utility.h"

namespace
{
using namespace libMesh;

// Compute the static coefficients for an element
template <typename RealType>
void hermite_compute_coefs(const ElemTempl<RealType> * elem, RealType & d1xd1x, RealType & d2xd2x)
{
  typedef PointTempl<RealType> Point;

  const FEFamily mapping_family = FEMapTempl<RealType>::map_fe_type(*elem);
  const Order mapping_order        (elem->default_order());
  const ElemType mapping_elem_type (elem->type());

  const FEType map_fe_type(mapping_order, mapping_family);

  const int n_mapping_shape_functions =
    FEInterface::n_shape_functions (1, map_fe_type,
                                    mapping_elem_type);

  // Degrees of freedom are at vertices and edge midpoints
  std::vector<Point> dofpt;
  dofpt.push_back(Point(-1));
  dofpt.push_back(Point(1));

  // Mapping functions - first derivatives at each dofpt
  std::vector<Real> dxdxi(2);
  std::vector<Real> dxidx(2);

  FEInterface::shape_deriv_ptr<RealType> shape_deriv_ptr =
    FEInterface::shape_deriv_function<RealType>(1, map_fe_type);

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[p] = 0;
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const auto ddxi = shape_deriv_ptr
            (elem, mapping_order, i, 0, dofpt[p], false);
          dxdxi[p] += elem->point(i)(0) * ddxi;
        }
    }

  // Calculate derivative scaling factors

  d1xd1x = dxdxi[0];
  d2xd2x = dxdxi[1];
}

} // end anonymous namespace


namespace libMesh
{

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
RealType FEHermiteShim<1,RealType>::hermite_raw_shape_second_deriv (const unsigned int i, const RealType & xi)
{
  using Utility::pow;

  switch (i)
    {
    case 0:
      return 1.5 * xi;
    case 1:
      return -1.5 * xi;
    case 2:
      return 0.5 * (-1. + 3.*xi);
    case 3:
      return 0.5 * (1. + 3.*xi);
    case 4:
      return (8.*xi*xi + 4.*(xi*xi-1.))/24.;
    case 5:
      return (8.*xi*xi*xi + 12.*xi*(xi*xi-1.))/120.;
      //      case 6:
      //        return (8.*pow<4>(xi) + 20.*xi*xi*(xi*xi-1.) +
      //          2.*(xi*xi-1)*(xi*xi-1))/720.;
    default:
      Real denominator = 720., xipower = 1.;
      for (unsigned n=6; n != i; ++n)
        {
          xipower *= xi;
          denominator *= (n+1);
        }
      return (8.*pow<4>(xi)*xipower +
              (8.*(i-4)+4.)*xi*xi*xipower*(xi*xi-1.) +
              (i-4)*(i-5)*xipower*(xi*xi-1.)*(xi*xi-1.))/denominator;
    }
}

#endif


template <typename RealType>
RealType FEHermiteShim<1,RealType>::hermite_raw_shape_deriv(const unsigned int i, const RealType & xi)
{
  switch (i)
    {
    case 0:
      return 0.75 * (-1. + xi*xi);
    case 1:
      return 0.75 * (1. - xi*xi);
    case 2:
      return 0.25 * (-1. - 2.*xi + 3.*xi*xi);
    case 3:
      return 0.25 * (-1. + 2.*xi + 3.*xi*xi);
    case 4:
      return 4.*xi * (xi*xi-1.)/24.;
    case 5:
      return (4*xi*xi*(xi*xi-1.) + (xi*xi-1.)*(xi*xi-1.))/120.;
      //      case 6:
      //        return (4*xi*xi*xi*(xi*xi-1.) + 2*xi*(xi*xi-1.)*(xi*xi-1.))/720.;
    default:
      Real denominator = 720., xipower = 1.;
      for (unsigned n=6; n != i; ++n)
        {
          xipower *= xi;
          denominator *= (n+1);
        }
      return (4*xi*xi*xi*xipower*(xi*xi-1.) +
              (i-4)*xi*xipower*(xi*xi-1.)*(xi*xi-1.))/denominator;
    }
}

template <typename RealType>
RealType FEHermiteShim<1,RealType>::hermite_raw_shape(const unsigned int i, const RealType & xi)
{
  switch (i)
    {
    case 0:
      return 0.25 * (2. - 3.*xi + xi*xi*xi);
    case 1:
      return 0.25 * (2. + 3.*xi - xi*xi*xi);
    case 2:
      return 0.25 * (1. - xi - xi*xi + xi*xi*xi);
    case 3:
      return 0.25 * (-1. - xi + xi*xi + xi*xi*xi);
      // All high order terms have the form x^(p-4)(x^2-1)^2/p!
    case 4:
      return (xi*xi-1.) * (xi*xi-1.)/24.;
    case 5:
      return xi * (xi*xi-1.) * (xi*xi-1.)/120.;
      //      case 6:
      //        return xi*xi * (xi*xi-1.) * (xi*xi-1.)/720.;
    default:
      Real denominator = 720., xipower = 1.;
      for (unsigned n=6; n != i; ++n)
        {
          xipower *= xi;
          denominator *= (n+1);
        }
      return (xi*xi*xipower*(xi*xi-1.)*(xi*xi-1.))/denominator;

    }
}


template <typename RealType>
typename FEShim<1,HERMITE,RealType>::OutputShape FEShim<1,HERMITE,RealType>::shape(const ElemType,
                          const Order,
                          const unsigned int,
                          const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,HERMITE,RealType>::OutputShape FEShim<1,HERMITE,RealType>::shape(const Elem * elem,
                          const Order libmesh_dbg_var(order),
                          const unsigned int i,
                          const Point & p,
                          const bool libmesh_dbg_var(add_p_level))
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  RealType d1xd1x, d2xd2x;

  hermite_compute_coefs(elem, d1xd1x, d2xd2x);

  const ElemType type = elem->type();

#ifndef NDEBUG
  const unsigned int totalorder =
    order + add_p_level * elem->p_level();
#endif

  switch (type)
    {
      // C1 functions on the C1 cubic edge
    case EDGE2:
    case EDGE3:
      {
        libmesh_assert_less (i, totalorder+1);

        switch (i)
          {
          case 0:
            return FEHermite<1>::hermite_raw_shape(0, p(0));
          case 1:
            return d1xd1x * FEHermite<1>::hermite_raw_shape(2, p(0));
          case 2:
            return FEHermite<1>::hermite_raw_shape(1, p(0));
          case 3:
            return d2xd2x * FEHermite<1>::hermite_raw_shape(3, p(0));
          default:
            return FEHermite<1>::hermite_raw_shape(i, p(0));
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}



template <typename RealType>
typename FEShim<1,HERMITE,RealType>::OutputShape FEShim<1,HERMITE,RealType>::shape_deriv(const ElemType,
                                const Order,
                                const unsigned int,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,HERMITE,RealType>::OutputShape FEShim<1,HERMITE,RealType>::shape_deriv(const Elem * elem,
                                const Order libmesh_dbg_var(order),
                                const unsigned int i,
                                const unsigned int,
                                const Point & p,
                                const bool libmesh_dbg_var(add_p_level))
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  RealType d1xd1x, d2xd2x;

  hermite_compute_coefs(elem, d1xd1x, d2xd2x);

  const ElemType type = elem->type();

#ifndef NDEBUG
  const unsigned int totalorder =
    order + add_p_level * elem->p_level();
#endif

  switch (type)
    {
      // C1 functions on the C1 cubic edge
    case EDGE2:
    case EDGE3:
      {
        libmesh_assert_less (i, totalorder+1);

        switch (i)
          {
          case 0:
            return FEHermite<1>::hermite_raw_shape_deriv(0, p(0));
          case 1:
            return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2, p(0));
          case 2:
            return FEHermite<1>::hermite_raw_shape_deriv(1, p(0));
          case 3:
            return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3, p(0));
          default:
            return FEHermite<1>::hermite_raw_shape_deriv(i, p(0));
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<1,HERMITE,RealType>::OutputShape FEShim<1,HERMITE,RealType>::shape_second_deriv(const ElemType,
                                       const Order,
                                       const unsigned int,
                                       const unsigned int,
                                       const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <typename RealType>
typename FEShim<1,HERMITE,RealType>::OutputShape FEShim<1,HERMITE,RealType>::shape_second_deriv(const Elem * elem,
                                       const Order libmesh_dbg_var(order),
                                       const unsigned int i,
                                       const unsigned int,
                                       const Point & p,
                                       const bool libmesh_dbg_var(add_p_level))
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  RealType d1xd1x, d2xd2x;

  hermite_compute_coefs(elem, d1xd1x, d2xd2x);

  const ElemType type = elem->type();

#ifndef NDEBUG
  const unsigned int totalorder =
    order + add_p_level * elem->p_level();
#endif

  switch (type)
    {
      // C1 functions on the C1 cubic edge
    case EDGE2:
    case EDGE3:
      {
        libmesh_assert_less (i, totalorder+1);

        switch (i)
          {
          case 0:
            return FEHermite<1>::hermite_raw_shape_second_deriv(0, p(0));
          case 1:
            return d1xd1x * FEHermite<1>::hermite_raw_shape_second_deriv(2, p(0));
          case 2:
            return FEHermite<1>::hermite_raw_shape_second_deriv(1, p(0));
          case 3:
            return d2xd2x * FEHermite<1>::hermite_raw_shape_second_deriv(3, p(0));
          default:
            return FEHermite<1>::hermite_raw_shape_second_deriv(i, p(0));
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_HERMITE_SHAPE_1D_IMPL_H
