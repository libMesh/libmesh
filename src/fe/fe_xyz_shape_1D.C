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


// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"


namespace libMesh
{


LIBMESH_DEFAULT_VECTORIZED_FE(1,XYZ)


template <>
Real FE<1,XYZ>::shape(const Elem * elem,
                      const Order libmesh_dbg_var(order),
                      const unsigned int i,
                      const Point & point_in,
                      const bool libmesh_dbg_var(add_p_level))
{
  libmesh_assert(elem);
  libmesh_assert_less_equal (i, order + add_p_level * elem->p_level());

  Point avg = elem->vertex_average();
  Real max_distance = 0.;
  for (const Point & p : elem->node_ref_range())
    {
      const Real distance = std::abs(avg(0) - p(0));
      max_distance = std::max(distance, max_distance);
    }

  const Real x  = point_in(0);
  const Real xc = avg(0);
  const Real dx = (x - xc)/max_distance;

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
    case 0:
      return 1.;

    case 1:
      return dx;

    case 2:
      return dx*dx;

    case 3:
      return dx*dx*dx;

    case 4:
      return dx*dx*dx*dx;

    default:
      Real val = 1.;
      for (unsigned int index = 0; index != i; ++index)
        val *= dx;
      return val;
    }
}



template <>
Real FE<1,XYZ>::shape(const ElemType,
                      const Order,
                      const unsigned int,
                      const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element.");
  return 0.;
}




template <>
Real FE<1,XYZ>::shape(const FEType fet,
                      const Elem * elem,
                      const unsigned int i,
                      const Point & p,
                      const bool add_p_level)
{
  return FE<1,XYZ>::shape(elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<1,XYZ>::shape_deriv(const Elem * elem,
                            const Order libmesh_dbg_var(order),
                            const unsigned int i,
                            const unsigned int libmesh_dbg_var(j),
                            const Point & point_in,
                            const bool libmesh_dbg_var(add_p_level))
{
  libmesh_assert(elem);
  libmesh_assert_less_equal (i, order + add_p_level * elem->p_level());

  // only d()/dxi in 1D!

  libmesh_assert_equal_to (j, 0);

  Point avg = elem->vertex_average();
  Real max_distance = 0.;
  for (const Point & p : elem->node_ref_range())
    {
      const Real distance = std::abs(avg(0) - p(0));
      max_distance = std::max(distance, max_distance);
    }

  const Real x  = point_in(0);
  const Real xc = avg(0);
  const Real dx = (x - xc)/max_distance;

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
    case 0:
      return 0.;

    case 1:
      return 1./max_distance;

    case 2:
      return 2.*dx/max_distance;

    case 3:
      return 3.*dx*dx/max_distance;

    case 4:
      return 4.*dx*dx*dx/max_distance;

    default:
      Real val = i;
      for (unsigned int index = 1; index != i; ++index)
        val *= dx;
      return val/max_distance;
    }
}



template <>
Real FE<1,XYZ>::shape_deriv(const ElemType,
                            const Order,
                            const unsigned int,
                            const unsigned int,
                            const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element.");
  return 0.;
}


template <>
Real FE<1,XYZ>::shape_deriv(const FEType fet,
                            const Elem * elem,
                            const unsigned int i,
                            const unsigned int j,
                            const Point & p,
                            const bool add_p_level)
{
  return FE<1,XYZ>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}





#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



template <>
Real FE<1,XYZ>::shape_second_deriv(const Elem * elem,
                                   const Order libmesh_dbg_var(order),
                                   const unsigned int i,
                                   const unsigned int libmesh_dbg_var(j),
                                   const Point & point_in,
                                   const bool libmesh_dbg_var(add_p_level))
{
  libmesh_assert(elem);
  libmesh_assert_less_equal (i, order + add_p_level * elem->p_level());

  // only d2()/dxi2 in 1D!

  libmesh_assert_equal_to (j, 0);

  Point avg = elem->vertex_average();
  Real max_distance = 0.;
  for (const Point & p : elem->node_ref_range())
    {
      const Real distance = std::abs(avg(0) - p(0));
      max_distance = std::max(distance, max_distance);
    }

  const Real x  = point_in(0);
  const Real xc = avg(0);
  const Real dx = (x - xc)/max_distance;
  const Real dist2 = pow(max_distance,2.);

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
    case 0:
    case 1:
      return 0.;

    case 2:
      return 2./dist2;

    case 3:
      return 6.*dx/dist2;

    case 4:
      return 12.*dx*dx/dist2;

    default:
      Real val = 2.;
      for (unsigned int index = 2; index != i; ++index)
        val *= (index+1) * dx;
      return val/dist2;
    }
}


template <>
Real FE<1,XYZ>::shape_second_deriv(const ElemType,
                                   const Order,
                                   const unsigned int,
                                   const unsigned int,
                                   const Point &)
{
  libmesh_error_msg("XYZ polynomials require the element.");
  return 0.;
}



template <>
Real FE<1,XYZ>::shape_second_deriv(const FEType fet,
                                   const Elem * elem,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p,
                                   const bool add_p_level)
{
  return FE<1,XYZ>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

#endif

} // namespace libMesh
