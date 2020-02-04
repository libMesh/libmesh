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

#ifndef LIBMESH_FE_CLOUGH_IMPL_SHAPE_1D_H
#define LIBMESH_FE_CLOUGH_IMPL_SHAPE_1D_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"


// Anonymous namespace for persistent variables.
// This allows us to cache the global-to-local mapping transformation
// This should also screw up multithreading royally
namespace
{
using namespace libMesh;

// Keep track of which element was most recently used to generate
// cached data
static dof_id_type old_elem_id = DofObject::invalid_id;
static const Elem * old_elem_ptr = nullptr;


// Compute the static coefficients for an element
template <typename RealType>
void clough_compute_coefs(const ElemTempl<RealType> * elem, RealType & d1xd1x, RealType & d2xd2x)
{
  typedef PointTempl<RealType> Point;

  // Using static globals for old_elem_id, etc. will fail
  // horribly with more than one thread.
  libmesh_assert_equal_to (libMesh::n_threads(), 1);

  // Coefficients are cached from old elements; we rely on that cache
  // except in dbg mode
#ifndef DEBUG
  if (elem->id() == old_elem_id &&
      elem == old_elem_ptr)
    return;
#endif

  const FEFamily mapping_family = FEMap::map_fe_type(*elem);

  const Order mapping_order        (elem->default_order());
  const ElemType mapping_elem_type (elem->type());

  const FEType map_fe_type(mapping_order, mapping_family);

  const int n_mapping_shape_functions =
    FEInterface::n_shape_functions(1, map_fe_type, mapping_elem_type);

  // Degrees of freedom are at vertices and edge midpoints
  std::vector<Point> dofpt;
  dofpt.push_back(Point(0));
  dofpt.push_back(Point(1));

  // Mapping functions - first derivatives at each dofpt
  std::vector<RealType> dxdxi(2);
  std::vector<RealType> dxidx(2);

  FEInterface::shape_deriv_ptr<RealType> shape_deriv_ptr =
    FEInterface::shape_deriv_function<RealType>(2, map_fe_type);

  for (int p = 0; p != 2; ++p)
    {
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const auto ddxi = shape_deriv_ptr
            (elem, mapping_order, i, 0, dofpt[p], false);
          dxdxi[p] += dofpt[p](0) * ddxi;
        }
    }

  // Calculate derivative scaling factors

#ifdef DEBUG
  // The cached factors should equal our calculations
  if (elem->id() == old_elem_id &&
      elem == old_elem_ptr)
    {
      libmesh_assert_equal_to(d1xd1x, dxdxi[0]);
      libmesh_assert_equal_to(d2xd2x, dxdxi[1]);
    }
#endif

  old_elem_id = elem->id();
  old_elem_ptr = elem;

  d1xd1x = dxdxi[0];
  d2xd2x = dxdxi[1];
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

// Return shape function second derivatives on the unit interval
template <typename RealType>
RealType clough_raw_shape_second_deriv(const unsigned int basis_num,
                                       const unsigned int deriv_type,
                                       const PointTempl<RealType> & p)
{
  auto xi = p(0);

  switch (deriv_type)
    {

      // second derivative in xi-xi direction
    case 0:
      switch (basis_num)
        {
        case 0:
          return -6 + 12*xi;
        case 1:
          return 6 - 12*xi;
        case 2:
          return -4 + 6*xi;
        case 3:
          return -2 + 6*xi;

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " <<
                        deriv_type);
    }
}

#endif


template <typename RealType>
RealType clough_raw_shape_deriv(const unsigned int basis_num,
                                const unsigned int deriv_type,
                                const PointTempl<RealType> & p)
{
  auto xi = p(0);

  switch (deriv_type)
    {
    case 0:
      switch (basis_num)
        {
        case 0:
          return -6*xi + 6*xi*xi;
        case 1:
          return 6*xi - 6*xi*xi;
        case 2:
          return 1 - 4*xi + 3*xi*xi;
        case 3:
          return -2*xi + 3*xi*xi;

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " <<
                        deriv_type);
    }
}

template <typename RealType>
RealType clough_raw_shape(const unsigned int basis_num,
                          const PointTempl<RealType> & p)
{
  auto xi = p(0);

  switch (basis_num)
    {
    case 0:
      return 1 - 3*xi*xi + 2*xi*xi*xi;
    case 1:
      return 3*xi*xi - 2*xi*xi*xi;
    case 2:
      return xi - 2*xi*xi + xi*xi*xi;
    case 3:
      return -xi*xi + xi*xi*xi;

    default:
      libmesh_error_msg("Invalid shape function index i = " <<
                        basis_num);
    }
}


} // end anonymous namespace


namespace libMesh
{


template <typename RealType>
typename FEShim<1,CLOUGH,RealType>::OutputShape FEShim<1,CLOUGH,RealType>::shape(const ElemType,
                         const Order,
                         const unsigned int,
                         const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,CLOUGH,RealType>::OutputShape FEShim<1,CLOUGH,RealType>::shape(const Elem * elem,
                         const Order order,
                         const unsigned int i,
                         const Point & p,
                         const bool add_p_level)
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static RealType d1xd1x, d2xd2x;

  clough_compute_coefs(elem, d1xd1x, d2xd2x);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 3rd-order C1 cubic element
    case THIRD:
      {
        switch (type)
          {
            // C1 functions on the C1 cubic edge
          case EDGE2:
          case EDGE3:
            {
              libmesh_assert_less (i, 4);

              switch (i)
                {
                case 0:
                  return clough_raw_shape(0, p);
                case 1:
                  return clough_raw_shape(1, p);
                case 2:
                  return d1xd1x * clough_raw_shape(2, p);
                case 3:
                  return d2xd2x * clough_raw_shape(3, p);
                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << totalorder);
    }
}



template <typename RealType>
typename FEShim<1,CLOUGH,RealType>::OutputShape FEShim<1,CLOUGH,RealType>::shape_deriv(const ElemType,
                               const Order,
                               const unsigned int,
                               const unsigned int,
                               const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,CLOUGH,RealType>::OutputShape FEShim<1,CLOUGH,RealType>::shape_deriv(const Elem * elem,
                               const Order order,
                               const unsigned int i,
                               const unsigned int j,
                               const Point & p,
                               const bool add_p_level)
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static RealType d1xd1x, d2xd2x;

  clough_compute_coefs(elem, d1xd1x, d2xd2x);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 3rd-order C1 cubic element
    case THIRD:
      {
        switch (type)
          {
            // C1 functions on the C1 cubic edge
          case EDGE2:
          case EDGE3:
            {
              switch (i)
                {
                case 0:
                  return clough_raw_shape_deriv(0, j, p);
                case 1:
                  return clough_raw_shape_deriv(1, j, p);
                case 2:
                  return d1xd1x * clough_raw_shape_deriv(2, j, p);
                case 3:
                  return d2xd2x * clough_raw_shape_deriv(3, j, p);
                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << totalorder);
    }
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<1,CLOUGH,RealType>::OutputShape FEShim<1,CLOUGH,RealType>::shape_second_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<1,CLOUGH,RealType>::OutputShape FEShim<1,CLOUGH,RealType>::shape_second_deriv(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static RealType d1xd1x, d2xd2x;

  clough_compute_coefs(elem, d1xd1x, d2xd2x);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 3rd-order C1 cubic element
    case THIRD:
      {
        switch (type)
          {
            // C1 functions on the C1 cubic edge
          case EDGE2:
          case EDGE3:
            {
              switch (i)
                {
                case 0:
                  return clough_raw_shape_second_deriv(0, j, p);
                case 1:
                  return clough_raw_shape_second_deriv(1, j, p);
                case 2:
                  return d1xd1x * clough_raw_shape_second_deriv(2, j, p);
                case 3:
                  return d2xd2x * clough_raw_shape_second_deriv(3, j, p);
                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << totalorder);
    }
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_CLOUGH_IMPL_SHAPE_1D_H
