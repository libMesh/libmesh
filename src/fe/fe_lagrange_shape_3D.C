// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_lagrange_shape_1D.h"
#include "libmesh/fe_serendipity_lagrange.h"
#include "libmesh/fe_simplex_lagrange.h"
#include "libmesh/fe_tensor_product_lagrange.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/cell_c0polyhedron.h"
#include "libmesh/tensor_value.h"

// Anonymous namespace for functions shared by LAGRANGE and
// L2_LAGRANGE implementations. Implementations appear at the bottom
// of this file.
namespace
{
using namespace libMesh;

template <FEFamily T>
Real fe_lagrange_3D_shape(const ElemType,
                          const Order order,
                          const Elem * elem,
                          const unsigned int i,
                          const Point & p);

template <FEFamily T>
Real fe_lagrange_3D_shape_deriv(const ElemType type,
                                const Order order,
                                const Elem * elem,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <FEFamily T>
Real fe_lagrange_3D_shape_second_deriv(const ElemType type,
                                       const Order order,
                                       const Elem * elem,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p);

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // anonymous namespace

namespace libMesh
{


// TODO: If optimizations for LAGRANGE work well we should do
// L2_LAGRANGE too...
LIBMESH_DEFAULT_VECTORIZED_FE(3,L2_LAGRANGE)


template<>
void FE<3,LAGRANGE>::all_shapes
  (const Elem * elem,
   const Order o,
   const std::vector<Point> & p,
   std::vector<std::vector<OutputShape>> & v,
   const bool add_p_level)
{
  const ElemType type = elem->type();

  // Just loop on the harder-to-optimize cases
  if (type != HEX8 && type != HEX27)
    {
      FE<3,LAGRANGE>::default_all_shapes
        (elem,o,p,v,add_p_level);
      return;
    }

#if LIBMESH_DIM == 3

  const unsigned int n_sf = v.size();

  switch (o)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // trilinear hexahedral shape functions
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less_equal (n_sf, 8);

              for (auto qp : index_range(p))
                {
                  const Point & q_point = p[qp];

                  for (unsigned int i : make_range(n_sf))
                    v[i][qp] = libMesh::detail::fe_lagrange_hex8_shape(i, q_point(0), q_point(1), q_point(2));
                }
              return;
            }

          default:
            libmesh_error(); // How did we get here?
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {
            // triquadratic hexahedral shape functions
          case HEX8:
// TODO: refactor to optimize this
//            libmesh_assert_msg(T == L2_LAGRANGE,
//                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case HEX27:
            {
              libmesh_assert_less_equal (n_sf, 27);

              for (auto qp : index_range(p))
                {
                  const Point & q_point = p[qp];

                  for (unsigned int i : make_range(n_sf))
                    v[i][qp] = libMesh::detail::fe_lagrange_hex27_shape(i, q_point(0), q_point(1), q_point(2));
                }
              return;
            }

          default:
            libmesh_error(); // How did we get here?
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order on HEX!: " << o);
    }
#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, o, p, v, add_p_level);
  libmesh_not_implemented();
#endif // LIBMESH_DIM == 3
}

template<>
void FE<3,LAGRANGE>::shapes
  (const Elem * elem,
   const Order o,
   const unsigned int i,
   const std::vector<Point> & p,
   std::vector<OutputShape> & v,
   const bool add_p_level)
{
  FE<3,LAGRANGE>::default_shapes
    (elem,o,i,p,v,add_p_level);
}

template<>
void FE<3,LAGRANGE>::shape_derivs
  (const Elem * elem,
   const Order o,
   const unsigned int i,
   const unsigned int j,
   const std::vector<Point> & p,
   std::vector<OutputShape> & v,
   const bool add_p_level)
{
  FE<3,LAGRANGE>::default_shape_derivs
    (elem,o,i,j,p,v,add_p_level);
}

template<>
void FE<3,LAGRANGE>::all_shape_derivs
  (const Elem * elem,
   const Order o,
   const std::vector<Point> & p,
   std::vector<std::vector<OutputShape>> * comps[3],
   const bool add_p_level)
{
  const ElemType type = elem->type();

  // Just loop on the harder-to-optimize cases
  if (type != HEX8 && type != HEX27)
    {
      FE<3,LAGRANGE>::default_all_shape_derivs
        (elem,o,p,comps,add_p_level);
      return;
    }

#if LIBMESH_DIM == 3

  libmesh_assert(comps[0]);
  libmesh_assert(comps[1]);
  libmesh_assert(comps[2]);
  const unsigned int n_sf = comps[0]->size();

  switch (o)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // trilinear hexahedral shape functions
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_equal_to (n_sf, 8);

              for (auto qp : index_range(p))
                {
                  const Point & q_point = p[qp];
                  for (unsigned int i : make_range(n_sf))
                    {
                      (*comps[0])[i][qp] = libMesh::detail::fe_lagrange_hex8_shape_deriv(i, 0, q_point(0), q_point(1), q_point(2));
                      (*comps[1])[i][qp] = libMesh::detail::fe_lagrange_hex8_shape_deriv(i, 1, q_point(0), q_point(1), q_point(2));
                      (*comps[2])[i][qp] = libMesh::detail::fe_lagrange_hex8_shape_deriv(i, 2, q_point(0), q_point(1), q_point(2));
                    }
                }
              return;
            }

          default:
            libmesh_error(); // How did we get here?
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {
            // triquadratic hexahedral shape functions
          case HEX8:
// TODO: refactor to optimize this
//            libmesh_assert_msg(T == L2_LAGRANGE,
//                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case HEX27:
            {
              libmesh_assert_less_equal (n_sf, 27);

              for (auto qp : index_range(p))
                {
                  const Point & q_point = p[qp];
                  for (unsigned int i : make_range(n_sf))
                    {
                      (*comps[0])[i][qp] = libMesh::detail::fe_lagrange_hex27_shape_deriv(i, 0, q_point(0), q_point(1), q_point(2));
                      (*comps[1])[i][qp] = libMesh::detail::fe_lagrange_hex27_shape_deriv(i, 1, q_point(0), q_point(1), q_point(2));
                      (*comps[2])[i][qp] = libMesh::detail::fe_lagrange_hex27_shape_deriv(i, 2, q_point(0), q_point(1), q_point(2));
                    }
                }
              return;
            }

          default:
            libmesh_error(); // How did we get here?
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order on HEX!: " << o);
    }
#else // LIBMESH_DIM != 3
  libmesh_ignore(elem, o, p, v, add_p_level);
  libmesh_not_implemented();
#endif // LIBMESH_DIM == 3
}




template <>
Real FE<3,LAGRANGE>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  return fe_lagrange_3D_shape<LAGRANGE>(type, order, nullptr, i, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape(const ElemType type,
                              const Order order,
                              const unsigned int i,
                              const Point & p)
{
  return fe_lagrange_3D_shape<L2_LAGRANGE>(type, order, nullptr, i, p);
}



template <>
Real FE<3,LAGRANGE>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_3D_shape<LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), elem, i, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape(const Elem * elem,
                              const Order order,
                              const unsigned int i,
                              const Point & p,
                              const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape functions
  return fe_lagrange_3D_shape<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), elem, i, p);
}



template <>
Real FE<3,LAGRANGE>::shape(const FEType fet,
                           const Elem * elem,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level)
{
  libmesh_assert(elem);
  return fe_lagrange_3D_shape<LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), elem, i, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape(const FEType fet,
                              const Elem * elem,
                              const unsigned int i,
                              const Point & p,
                              const bool add_p_level)
{
  libmesh_assert(elem);
  return fe_lagrange_3D_shape<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), elem, i, p);
}

template <>
Real FE<3,LAGRANGE>::shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  return fe_lagrange_3D_shape_deriv<LAGRANGE>(type, order, nullptr, i, j, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape_deriv(const ElemType type,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  return fe_lagrange_3D_shape_deriv<L2_LAGRANGE>(type, order, nullptr, i, j, p);
}



template <>
Real FE<3,LAGRANGE>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_deriv<LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), elem, i, j, p);
}


template <>
Real FE<3,L2_LAGRANGE>::shape_deriv(const Elem * elem,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p,
                                    const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_deriv<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), elem, i, j, p);
}


template <>
Real FE<3,LAGRANGE>::shape_deriv(const FEType fet,
                                 const Elem * elem,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level)
{
  libmesh_assert(elem);
  return fe_lagrange_3D_shape_deriv<LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), elem, i, j, p);
}


template <>
Real FE<3,L2_LAGRANGE>::shape_deriv(const FEType fet,
                                    const Elem * elem,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p,
                                    const bool add_p_level)
{
  libmesh_assert(elem);
  return fe_lagrange_3D_shape_deriv<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), elem, i, j, p);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
Real FE<3,LAGRANGE>::shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  return fe_lagrange_3D_shape_second_deriv<LAGRANGE>(type, order, nullptr, i, j, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape_second_deriv(const ElemType type,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p)
{
  return fe_lagrange_3D_shape_second_deriv<L2_LAGRANGE>(type, order, nullptr, i, j, p);
}



template <>
Real FE<3,LAGRANGE>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_second_deriv<LAGRANGE>
    (elem->type(), order + add_p_level*elem->p_level(), elem, i, j, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape_second_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);

  // call the orientation-independent shape function derivatives
  return fe_lagrange_3D_shape_second_deriv<L2_LAGRANGE>
    (elem->type(), order + add_p_level*elem->p_level(), elem, i, j, p);
}


template <>
Real FE<3,LAGRANGE>::shape_second_deriv(const FEType fet,
                                        const Elem * elem,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level)
{
  libmesh_assert(elem);
  return fe_lagrange_3D_shape_second_deriv<LAGRANGE>
    (elem->type(), fet.order + add_p_level*elem->p_level(), elem, i, j, p);
}



template <>
Real FE<3,L2_LAGRANGE>::shape_second_deriv(const FEType fet,
                                           const Elem * elem,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level)
{
  libmesh_assert(elem);
  return fe_lagrange_3D_shape_second_deriv<L2_LAGRANGE>
    (elem->type(), fet.order + add_p_level*elem->p_level(), elem, i, j, p);
}


#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh



namespace
{
using namespace libMesh;

template <FEFamily T>
Real fe_lagrange_3D_shape(const ElemType type,
                          const Order order,
                          const Elem * elem,
                          const unsigned int i,
                          const Point & p)
{
#if LIBMESH_DIM == 3

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // trilinear hexahedral shape functions
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              return libMesh::detail::fe_lagrange_hex8_shape(i, p(0), p(1), p(2));
            }

            // linear tetrahedral shape functions
          case TET4:
          case TET10:
          case TET14:
            {
              libmesh_assert_less (i, 4);
              return libMesh::detail::fe_lagrange_tet4_shape(i, p(0), p(1), p(2));
            }

            // linear prism shape functions
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              libmesh_assert_less (i, 6);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              return (FE<2,LAGRANGE>::shape(TRI3,  FIRST, i1[i], p2d)*
                      fe_lagrange_1D_linear_shape(i0[i], p1d));
            }

            // linear pyramid shape functions
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
            {
              libmesh_assert_less (i, 5);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              switch(i)
                {
                case 0:
                  return .25*(zeta + xi - 1.)*(zeta + eta - 1.)/((1. - zeta) + eps);

                case 1:
                  return .25*(zeta - xi - 1.)*(zeta + eta - 1.)/((1. - zeta) + eps);

                case 2:
                  return .25*(zeta - xi - 1.)*(zeta - eta - 1.)/((1. - zeta) + eps);

                case 3:
                  return .25*(zeta + xi - 1.)*(zeta - eta - 1.)/((1. - zeta) + eps);

                case 4:
                  return zeta;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case C0POLYHEDRON:
            {
              // Polyhedra require using newer FE APIs
              if (!elem)
                libmesh_error_msg("Code (see stack trace) used an outdated FE function overload.\n"
                                  "Shape functions on a polyhedron are not defined by ElemType alone.");

              libmesh_assert(elem->type() == C0POLYHEDRON);

              const C0Polyhedron & poly = *cast_ptr<const C0Polyhedron *>(elem);

              // We can't use a small tolerance here, because in
              // inverse_map() Newton might hand us intermediate
              // iterates outside the polyhedron.
              const auto [s, a, b, c] = poly.subelement_coordinates(p, 100);
              if (s == invalid_uint)
                return 0;
              libmesh_assert_less(s, poly.n_subelements());

              const auto subtet = poly.subelement(s);

              // Avoid signed/unsigned comparison warnings
              const int nodei = i;
              if (nodei == subtet[0])
                return 1-a-b-c;
              if (nodei == subtet[1])
                return a;
              if (nodei == subtet[2])
                return b;
              if (nodei == subtet[3])
                return c;

              // Basis function i is not supported on p's subtet
              return 0;
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {

            // serendipity hexahedral quadratic shape functions
          case HEX20:
            {
              libmesh_assert_less (i, 20);
              return libMesh::detail::fe_lagrange_hex20_shape(i, p(0), p(1), p(2));
            }

            // triquadratic hexahedral shape functions
          case HEX8:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case HEX27:
            {
              libmesh_assert_less (i, 27);

              return libMesh::detail::fe_lagrange_hex27_shape(i, p(0), p(1), p(2));
            }

            // quadratic tetrahedral shape functions
          case TET4:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case TET10:
            libmesh_assert_less (i, 10);
            libmesh_fallthrough();
          case TET14:
            {
              libmesh_assert_less (i, 14);
              return libMesh::detail::fe_lagrange_tet10_shape(i, p(0), p(1), p(2));
            }

            // "serendipity" prism
          case PRISM15:
            {
              libmesh_assert_less (i, 15);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              switch(i)
                {
                case 0:
                  return (1. - zeta)*(xi + eta - 1.)*(xi + eta + 0.5*zeta);

                case 1:
                  return (1. - zeta)*xi*(xi - 1. - 0.5*zeta);

                case 2: // phi1 with xi <- eta
                  return (1. - zeta)*eta*(eta - 1. - 0.5*zeta);

                case 3: // phi0 with zeta <- (-zeta)
                  return (1. + zeta)*(xi + eta - 1.)*(xi + eta - 0.5*zeta);

                case 4: // phi1 with zeta <- (-zeta)
                  return (1. + zeta)*xi*(xi - 1. + 0.5*zeta);

                case 5: // phi4 with xi <- eta
                  return (1. + zeta)*eta*(eta - 1. + 0.5*zeta);

                case 6:
                  return 2.*(1. - zeta)*xi*(1. - xi - eta);

                case 7:
                  return 2.*(1. - zeta)*xi*eta;

                case 8:
                  return 2.*(1. - zeta)*eta*(1. - xi - eta);

                case 9:
                  return (1. - zeta)*(1. + zeta)*(1. - xi - eta);

                case 10:
                  return (1. - zeta)*(1. + zeta)*xi;

                case 11: // phi10 with xi <-> eta
                  return (1. - zeta)*(1. + zeta)*eta;

                case 12: // phi6 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*xi*(1. - xi - eta);

                case 13: // phi7 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*xi*eta;

                case 14: // phi8 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*eta*(1. - xi - eta);

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // quadratic prism shape functions
          case PRISM6:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              libmesh_assert_less (i, 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              return (FE<2,LAGRANGE>::shape(TRI6,  SECOND, i1[i], p2d)*
                      fe_lagrange_1D_quadratic_shape(i0[i], p1d));
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              libmesh_assert_less (i, 13);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              Real den = (1. - zeta + eps);

              switch(i)
                {
                case 0:
                  return 0.25*(-xi - eta - 1.)*((1. - xi)*(1. - eta) - zeta + xi*eta*zeta/den);

                case 1:
                  return 0.25*(-eta + xi - 1.)*((1. + xi)*(1. - eta) - zeta - xi*eta*zeta/den);

                case 2:
                  return 0.25*(xi + eta - 1.)*((1. + xi)*(1. + eta) - zeta + xi*eta*zeta/den);

                case 3:
                  return 0.25*(eta - xi - 1.)*((1. - xi)*(1. + eta) - zeta - xi*eta*zeta/den);

                case 4:
                  return zeta*(2.*zeta - 1.);

                case 5:
                  return 0.5*(1. + xi - zeta)*(1. - xi - zeta)*(1. - eta - zeta)/den;

                case 6:
                  return 0.5*(1. + eta - zeta)*(1. - eta - zeta)*(1. + xi - zeta)/den;

                case 7:
                  return 0.5*(1. + xi - zeta)*(1. - xi - zeta)*(1. + eta - zeta)/den;

                case 8:
                  return 0.5*(1. + eta - zeta)*(1. - eta - zeta)*(1. - xi - zeta)/den;

                case 9:
                  return zeta*(1. - xi - zeta)*(1. - eta - zeta)/den;

                case 10:
                  return zeta*(1. + xi - zeta)*(1. - eta - zeta)/den;

                case 11:
                  return zeta*(1. + eta - zeta)*(1. + xi - zeta)/den;

                case 12:
                  return zeta*(1. - xi - zeta)*(1. + eta - zeta)/den;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID5:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case PYRAMID14:
          case PYRAMID18:
            {
              libmesh_assert_less (i, 14);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              Real
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              Real
                den = (-1. + zeta + eps),
                den2 = den*den;

              switch(i)
                {
                case 0:
                  return p4*p1*(xi*eta - zeta + zeta*zeta)/den2;

                case 1:
                  return -p1*p2*(xi*eta + zeta - zeta*zeta)/den2;

                case 2:
                  return p2*p3*(xi*eta - zeta + zeta*zeta)/den2;

                case 3:
                  return -p3*p4*(xi*eta + zeta - zeta*zeta)/den2;

                case 4:
                  return zeta*(2.*zeta - 1.);

                case 5:
                  return -4.*p2*p1*p4*eta/den2;

                case 6:
                  return 4.*p1*p2*p3*xi/den2;

                case 7:
                  return 4.*p2*p3*p4*eta/den2;

                case 8:
                  return -4.*p3*p4*p1*xi/den2;

                case 9:
                  return -4.*p1*p4*zeta/den;

                case 10:
                  return -4.*p2*p1*zeta/den;

                case 11:
                  return -4.*p3*p2*zeta/den;

                case 12:
                  return -4.*p4*p3*zeta/den;

                case 13:
                  return 16.*p1*p2*p3*p4/den2;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }

    case THIRD:
      {
        switch (type)
          {
            // quadratic Lagrange shape functions with cubic bubbles
          case PRISM20:
            {
              libmesh_assert_less (i, 20);

              // Compute Prism21 shape functions as a tensor-product
              // of a triangle and an edge, then redistribute the
              // central bubble function over the other Tri6 nodes
              // around it (in a way consistent with the Tri7 shape
              // function definitions).
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1};

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              const Real mainval = FE<3,LAGRANGE>::shape(PRISM21, THIRD, i, p);

              if (i0[i] != 2)
                return mainval;

              const Real bubbleval =
                FE<2,LAGRANGE>::shape(TRI7, THIRD, 6, p2d) *
                fe_lagrange_1D_quadratic_shape(2, p1d);

              if (i < 12) // vertices
                return mainval - bubbleval / 9;

              return mainval + bubbleval * (Real(4) / 9);
            }

            // quadratic Lagrange shape functions with cubic bubbles
          case PRISM21:
            {
              libmesh_assert_less (i, 21);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 6, 6};

              return (FE<2,LAGRANGE>::shape(TRI7, THIRD, i1[i], p2d)*
                      fe_lagrange_1D_quadratic_shape(i0[i], p1d));
            }

            // Weird rational shape functions with weirder bubbles...
          case PYRAMID18:
            {
              libmesh_assert_less (i, 18);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              const Real
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              const Real
                den = (-1. + zeta + eps),
                den2 = den*den;

              // Bubble functions on triangular sides.  We actually
              // have a degree of freedom to play with here, and I'm
              // not certain how best to use it, so let's leave it
              // as a variable in case we figure that out later.
              constexpr Real alpha = 0.5;
              const Real
                bub_f1 = ((1-alpha)*(1-zeta) + alpha*(-eta)),
                bub_f2 = ((1-alpha)*(1-zeta) + alpha*(xi)),
                bub_f3 = ((1-alpha)*(1-zeta) + alpha*(eta)),
                bub_f4 = ((1-alpha)*(1-zeta) + alpha*(-xi));

              const Real
                bub1 = bub_f1*p1*p2*p4*zeta/den2,
                bub2 = bub_f2*p1*p2*p3*zeta/den2,
                bub3 = bub_f3*p2*p3*p4*zeta/den2,
                bub4 = bub_f4*p1*p3*p4*zeta/den2;

              switch(i)
                {
                case 0:
                  return p4*p1*(xi*eta - zeta + zeta*zeta)/den2 + 3*(bub1+bub4);

                case 1:
                  return -p1*p2*(xi*eta + zeta - zeta*zeta)/den2 + 3*(bub1+bub2);

                case 2:
                  return p2*p3*(xi*eta - zeta + zeta*zeta)/den2 + 3*(bub2+bub3);

                case 3:
                  return -p3*p4*(xi*eta + zeta - zeta*zeta)/den2 + 3*(bub3+bub4);

                case 4:
                  return zeta*(2.*zeta - 1.) + 3*(bub1+bub2+bub3+bub4);

                case 5:
                  return -4.*p2*p1*p4*eta/den2 - 12*bub1;

                case 6:
                  return 4.*p1*p2*p3*xi/den2 - 12*bub2;

                case 7:
                  return 4.*p2*p3*p4*eta/den2 - 12*bub3;

                case 8:
                  return -4.*p3*p4*p1*xi/den2 - 12*bub4;

                case 9:
                  return -4.*p1*p4*zeta/den - 12*(bub1+bub4);

                case 10:
                  return -4.*p2*p1*zeta/den - 12*(bub1+bub2);

                case 11:
                  return -4.*p3*p2*zeta/den - 12*(bub2+bub3);

                case 12:
                  return -4.*p4*p3*zeta/den - 12*(bub3+bub4);

                case 13:
                  return 16.*p1*p2*p3*p4/den2;

                case 14:
                  return 27*bub1;

                case 15:
                  return 27*bub2;

                case 16:
                  return 27*bub3;

                case 17:
                  return 27*bub4;

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

            // quadratic Lagrange shape functions with cubic bubbles
          case TET14:
            {
              libmesh_assert_less (i, 14);
              return libMesh::detail::fe_lagrange_tet14_shape(i, p(0), p(1), p(2));
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << order);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(type, order, elem, i, p);
  libmesh_not_implemented();
#endif
}



template <FEFamily T>
Real fe_lagrange_3D_shape_deriv(const ElemType type,
                                const Order order,
                                const Elem * elem,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p)
{
#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 3);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // trilinear hexahedral shape functions
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);

              return libMesh::detail::fe_lagrange_hex8_shape_deriv(i, j, p(0), p(1), p(2));
            }

            // linear tetrahedral shape functions
          case TET4:
          case TET10:
          case TET14:
            {
              libmesh_assert_less (i, 4);
              return libMesh::detail::fe_lagrange_tet4_shape_deriv(i, j);
            }

            // linear prism shape functions
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              libmesh_assert_less (i, 6);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI3,  FIRST, i1[i], 0, p2d)*
                          fe_lagrange_1D_linear_shape(i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI3,  FIRST, i1[i], 1, p2d)*
                          fe_lagrange_1D_linear_shape(i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FE<2,LAGRANGE>::shape(TRI3,  FIRST, i1[i], p2d)*
                          fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // linear pyramid shape functions
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
            {
              libmesh_assert_less (i, 5);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return  .25*(zeta + eta - 1.)/((1. - zeta) + eps);

                    case 1:
                      return -.25*(zeta + eta - 1.)/((1. - zeta) + eps);

                    case 2:
                      return -.25*(zeta - eta - 1.)/((1. - zeta) + eps);

                    case 3:
                      return  .25*(zeta - eta - 1.)/((1. - zeta) + eps);

                    case 4:
                      return 0;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return  .25*(zeta + xi - 1.)/((1. - zeta) + eps);

                    case 1:
                      return  .25*(zeta - xi - 1.)/((1. - zeta) + eps);

                    case 2:
                      return -.25*(zeta - xi - 1.)/((1. - zeta) + eps);

                    case 3:
                      return -.25*(zeta + xi - 1.)/((1. - zeta) + eps);

                    case 4:
                      return 0;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/dzeta
                case 2:
                  {
                    // We computed the derivatives with general eps and
                    // then let eps tend to zero in the numerators...
                    Real
                      num = zeta*(2. - zeta) - 1.,
                      den = (1. - zeta + eps)*(1. - zeta + eps);

                    switch(i)
                      {
                      case 0:
                      case 2:
                        return .25*(num + xi*eta)/den;

                      case 1:
                      case 3:
                        return .25*(num - xi*eta)/den;

                      case 4:
                        return 1.;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case C0POLYHEDRON:
            {
              // Polyhedra require using newer FE APIs
              if (!elem)
                libmesh_error_msg("Code (see stack trace) used an outdated FE function overload.\n"
                                  "Shape functions on a polyhedron are not defined by ElemType alone.");

              libmesh_assert(elem->type() == C0POLYHEDRON);

              const C0Polyhedron & poly = *cast_ptr<const C0Polyhedron *>(elem);

              // We can't use a small tolerance here, because in
              // inverse_map() Newton might hand us intermediate
              // iterates outside the polyhedron.
              const auto [s, a, b, c] = poly.subelement_coordinates(p, 100);
              if (s == invalid_uint)
                return 0;
              libmesh_assert_less(s, poly.n_subelements());

              const auto subtet = poly.subelement(s);

              // Find derivatives w.r.t. subtriangle barycentric
              // coordinates
              Real du_da = 0, du_db = 0, du_dc = 0;

              // Avoid signed/unsigned comparison warnings
              const int nodei = i;
              if (nodei == subtet[0])
                du_da = du_db = du_dc = -1;
              else if (nodei == subtet[1])
                du_da = 1;
              else if (nodei == subtet[2])
                du_db = 1;
              else if (nodei == subtet[3])
                du_dc = 1;
              else
                // Basis function i is not supported on p's subtet
                return 0;

              // We want to return derivatives with respect to
              // xi/eta/zeta for the polyhedron, but what we
              // calculated above are with respect to xi and eta
              // coordinates for a master *triangle*.  We need to
              // convert from one to the other.

              const auto master_points = poly.master_subelement(s);

              const RealTensor dXi_dA(
                master_points[1](0) - master_points[0](0), master_points[2](0) - master_points[0](0), master_points[3](0) - master_points[0](0),
                master_points[1](1) - master_points[0](1), master_points[2](1) - master_points[0](1), master_points[3](1) - master_points[0](1),
                master_points[1](2) - master_points[0](2), master_points[2](2) - master_points[0](2), master_points[3](2) - master_points[0](2));

              // When we vectorize this we'll want a full inverse, but
              // when we're querying one component at a time it's
              // cheaper to manually compute a single column.
              // const RealTensor dabc_dxietazeta_dabc = dxietazeta_dabc.inverse();
              const Real jac = dXi_dA.det();

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    const Real da_dxi =  (dXi_dA(2,2)*dXi_dA(1,1) - dXi_dA(2,1)*dXi_dA(1,2)) / jac;
                    const Real db_dxi = -(dXi_dA(2,2)*dXi_dA(1,0) - dXi_dA(2,0)*dXi_dA(1,2)) / jac;
                    const Real dc_dxi =  (dXi_dA(2,1)*dXi_dA(1,0) - dXi_dA(2,0)*dXi_dA(1,1)) / jac;
                    return du_da*da_dxi + du_db*db_dxi + du_dc*dc_dxi;
                  }
                  // d()/deta
                case 1:
                  {
                    const Real da_deta = -(dXi_dA(2,2)*dXi_dA(0,1) - dXi_dA(2,1)*dXi_dA(0,2) ) / jac;
                    const Real db_deta =  (dXi_dA(2,2)*dXi_dA(0,0) - dXi_dA(2,0)*dXi_dA(0,2)) / jac;
                    const Real dc_deta = -(dXi_dA(2,1)*dXi_dA(0,0) - dXi_dA(2,0)*dXi_dA(0,1)) / jac;
                    return du_da*da_deta + du_db*db_deta + du_dc*dc_deta;
                  }
                  // d()/dzeta
                case 2:
                  {
                    const Real da_dzeta =  (dXi_dA(1,2)*dXi_dA(0,1) - dXi_dA(1,1)*dXi_dA(0,2) ) / jac;
                    const Real db_dzeta = -(dXi_dA(1,2)*dXi_dA(0,0) - dXi_dA(1,0)*dXi_dA(0,2)) / jac;
                    const Real dc_dzeta =  (dXi_dA(1,1)*dXi_dA(0,0) - dXi_dA(1,0)*dXi_dA(0,1)) / jac;
                    return du_da*da_dzeta + du_db*db_dzeta + du_dc*dc_dzeta;
                  }
                default:
                  libmesh_error_msg("ERROR: Invalid derivative index j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {

            // serendipity hexahedral quadratic shape functions
          case HEX20:
            {
              libmesh_assert_less (i, 20);
              return libMesh::detail::fe_lagrange_hex20_shape_deriv(i, j, p(0), p(1), p(2));
            }

            // triquadratic hexahedral shape functions
          case HEX8:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case HEX27:
            {
              libmesh_assert_less (i, 27);

              return libMesh::detail::fe_lagrange_hex27_shape_deriv(i, j, p(0), p(1), p(2));
            }

            // quadratic tetrahedral shape functions
          case TET4:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case TET10:
          case TET14:
            {
              libmesh_assert_less (i, 10);
              return libMesh::detail::fe_lagrange_tet10_shape_deriv(i, j, p(0), p(1), p(2));
            }


            // "serendipity" prism
          case PRISM15:
            {
              libmesh_assert_less (i, 15);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return (2.*xi + 2.*eta + 0.5*zeta - 1.)*(1. - zeta);
                      case 1:
                        return (2.*xi - 1. - 0.5*zeta)*(1. - zeta);
                      case 2:
                        return 0.;
                      case 3:
                        return (2.*xi + 2.*eta - 0.5*zeta - 1.)*(1. + zeta);
                      case 4:
                        return (2.*xi - 1. + 0.5*zeta)*(1. + zeta);
                      case 5:
                        return 0.;
                      case 6:
                        return (4.*xi + 2.*eta - 2.)*(zeta - 1.);
                      case 7:
                        return -2.*(zeta - 1.)*eta;
                      case 8:
                        return 2.*(zeta - 1.)*eta;
                      case 9:
                        return (zeta - 1.)*(1. + zeta);
                      case 10:
                        return (1. - zeta)*(1. + zeta);
                      case 11:
                        return 0.;
                      case 12:
                        return (-4.*xi - 2.*eta + 2.)*(1. + zeta);
                      case 13:
                        return 2.*(1. + zeta)*eta;
                      case 14:
                        return -2.*(1. + zeta)*eta;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return (2.*xi + 2.*eta + 0.5*zeta - 1.)*(1. - zeta);
                      case 1:
                        return 0.;
                      case 2:
                        return (2.*eta - 1. - 0.5*zeta)*(1. - zeta);
                      case 3:
                        return (2.*xi + 2.*eta - 0.5*zeta - 1.)*(1. + zeta);
                      case 4:
                        return 0.;
                      case 5:
                        return (2.*eta - 1. + 0.5*zeta)*(1. + zeta);
                      case 6:
                        return 2.*(zeta - 1.)*xi;
                      case 7:
                        return 2.*(1. - zeta)*xi;
                      case 8:
                        return (2.*xi + 4.*eta - 2.)*(zeta - 1.);
                      case 9:
                        return (zeta - 1.)*(1. + zeta);
                      case 10:
                        return 0.;
                      case 11:
                        return (1. - zeta)*(1. + zeta);
                      case 12:
                        return -2.*(1. + zeta)*xi;
                      case 13:
                        return 2.*(1. + zeta)*xi;
                      case 14:
                        return (-2.*xi - 4.*eta + 2.)*(1. + zeta);

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return (-xi - eta - zeta + 0.5)*(xi + eta - 1.);
                      case 1:
                        return -0.5*xi*(2.*xi - 1. - 2.*zeta);
                      case 2:
                        return -0.5*eta*(2.*eta - 1. - 2.*zeta);
                      case 3:
                        return (xi + eta - zeta - 0.5)*(xi + eta - 1.);
                      case 4:
                        return 0.5*xi*(2.*xi - 1. + 2.*zeta);
                      case 5:
                        return 0.5*eta*(2.*eta - 1. + 2.*zeta);
                      case 6:
                        return 2.*xi*(xi + eta - 1.);
                      case 7:
                        return -2.*xi*eta;
                      case 8:
                        return 2.*eta*(xi + eta - 1.);
                      case 9:
                        return 2.*zeta*(xi + eta - 1.);
                      case 10:
                        return -2.*xi*zeta;
                      case 11:
                        return -2.*eta*zeta;
                      case 12:
                        return 2.*xi*(1. - xi - eta);
                      case 13:
                        return 2.*xi*eta;
                      case 14:
                        return 2.*eta*(1. - xi - eta);

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



            // quadratic prism shape functions
          case PRISM6:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              libmesh_assert_less (i, 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FE<2,LAGRANGE>::shape(TRI6,  SECOND, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              libmesh_assert_less (i, 13);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              Real
                den = (-1. + zeta + eps),
                den2 = den*den,
                xi2 = xi*xi,
                eta2 = eta*eta,
                zeta2 = zeta*zeta,
                zeta3 = zeta2*zeta;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return 0.25*(-zeta - eta + 2.*eta*zeta - 2.*xi + 2.*zeta*xi + 2.*eta*xi + zeta2 + eta2)/den;

                    case 1:
                      return -0.25*(-zeta - eta + 2.*eta*zeta + 2.*xi - 2.*zeta*xi - 2.*eta*xi + zeta2 + eta2)/den;

                    case 2:
                      return -0.25*(-zeta + eta - 2.*eta*zeta + 2.*xi - 2.*zeta*xi + 2.*eta*xi + zeta2 + eta2)/den;

                    case 3:
                      return 0.25*(-zeta + eta - 2.*eta*zeta - 2.*xi + 2.*zeta*xi - 2.*eta*xi + zeta2 + eta2)/den;

                    case 4:
                      return 0.;

                    case 5:
                      return -(-1. + eta + zeta)*xi/den;

                    case 6:
                      return 0.5*(-1. + eta + zeta)*(1. + eta - zeta)/den;

                    case 7:
                      return (1. + eta - zeta)*xi/den;

                    case 8:
                      return -0.5*(-1. + eta + zeta)*(1. + eta - zeta)/den;

                    case 9:
                      return -(-1. + eta + zeta)*zeta/den;

                    case 10:
                      return (-1. + eta + zeta)*zeta/den;

                    case 11:
                      return -(1. + eta - zeta)*zeta/den;

                    case 12:
                      return (1. + eta - zeta)*zeta/den;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return 0.25*(-zeta - 2.*eta + 2.*eta*zeta - xi + 2.*zeta*xi + 2.*eta*xi + zeta2 + xi2)/den;

                    case 1:
                      return -0.25*(zeta + 2.*eta - 2.*eta*zeta - xi + 2.*zeta*xi + 2.*eta*xi - zeta2 - xi2)/den;

                    case 2:
                      return -0.25*(-zeta + 2.*eta - 2.*eta*zeta + xi - 2.*zeta*xi + 2.*eta*xi + zeta2 + xi2)/den;

                    case 3:
                      return 0.25*(zeta - 2.*eta + 2.*eta*zeta + xi - 2.*zeta*xi + 2.*eta*xi - zeta2 - xi2)/den;

                    case 4:
                      return 0.;

                    case 5:
                      return -0.5*(-1. + xi + zeta)*(1. + xi - zeta)/den;

                    case 6:
                      return (1. + xi - zeta)*eta/den;

                    case 7:
                      return 0.5*(-1. + xi + zeta)*(1. + xi - zeta)/den;

                    case 8:
                      return -(-1. + xi + zeta)*eta/den;

                    case 9:
                      return -(-1. + xi + zeta)*zeta/den;

                    case 10:
                      return (1. + xi - zeta)*zeta/den;

                    case 11:
                      return -(1. + xi - zeta)*zeta/den;

                    case 12:
                      return (-1. + xi + zeta)*zeta/den;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                  // d/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.25*(xi + eta + 1.)*(-1. + 2.*zeta - zeta2 + eta*xi)/den2;

                      case 1:
                        return 0.25*(eta - xi + 1.)*(1. - 2.*zeta + zeta2 + eta*xi)/den2;

                      case 2:
                        return 0.25*(xi + eta - 1.)*(-1. + 2.*zeta - zeta2 + eta*xi)/den2;

                      case 3:
                        return -0.25*(eta - xi - 1.)*(1. - 2.*zeta + zeta2 + eta*xi)/den2;

                      case 4:
                        return 4.*zeta - 1.;

                      case 5:
                        return 0.5*(-2 + eta + 6.*zeta + eta*xi2 + eta*zeta2 - 6.*zeta2 + 2.*zeta3 - 2.*eta*zeta)/den2;

                      case 6:
                        return -0.5*(2 - 6.*zeta + xi + xi*zeta2 + eta2*xi + 6.*zeta2 - 2.*zeta3 - 2.*zeta*xi)/den2;

                      case 7:
                        return -0.5*(2 + eta - 6.*zeta + eta*xi2 + eta*zeta2 + 6.*zeta2 - 2.*zeta3 - 2.*eta*zeta)/den2;

                      case 8:
                        return 0.5*(-2 + 6.*zeta + xi + xi*zeta2 + eta2*xi - 6.*zeta2 + 2.*zeta3 - 2.*zeta*xi)/den2;

                      case 9:
                        return (1. - eta - 4.*zeta - xi - xi*zeta2 - eta*zeta2 + eta*xi + 5.*zeta2 - 2.*zeta3 + 2.*eta*zeta + 2.*zeta*xi)/den2;

                      case 10:
                        return -(-1. + eta + 4.*zeta - xi - xi*zeta2 + eta*zeta2 + eta*xi - 5.*zeta2 + 2.*zeta3 - 2.*eta*zeta + 2.*zeta*xi)/den2;

                      case 11:
                        return (1. + eta - 4.*zeta + xi + xi*zeta2 + eta*zeta2 + eta*xi + 5.*zeta2 - 2.*zeta3 - 2.*eta*zeta - 2.*zeta*xi)/den2;

                      case 12:
                        return -(-1. - eta + 4.*zeta + xi + xi*zeta2 - eta*zeta2 + eta*xi - 5.*zeta2 + 2.*zeta3 + 2.*eta*zeta - 2.*zeta*xi)/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID5:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case PYRAMID14:
          case PYRAMID18:
            {
              libmesh_assert_less (i, 14);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              Real
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              Real
                den = (-1. + zeta + eps),
                den2 = den*den,
                den3 = den2*den;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return 0.5*p1*(-xi*eta + zeta - zeta*zeta + 2.*p4*eta)/den2;

                    case 1:
                      return -0.5*p1*(xi*eta + zeta - zeta*zeta + 2.*p2*eta)/den2;

                    case 2:
                      return 0.5*p3*(xi*eta - zeta + zeta*zeta + 2.*p2*eta)/den2;

                    case 3:
                      return -0.5*p3*(-xi*eta - zeta + zeta*zeta + 2.*p4*eta)/den2;

                    case 4:
                      return 0.;

                    case 5:
                      return 2.*p1*eta*xi/den2;

                    case 6:
                      return 2.*p1*p3*(xi + 2.*p2)/den2;

                    case 7:
                      return -2.*p3*eta*xi/den2;

                    case 8:
                      return -2.*p1*p3*(-xi + 2.*p4)/den2;

                    case 9:
                      return 2.*p1*zeta/den;

                    case 10:
                      return -2.*p1*zeta/den;

                    case 11:
                      return -2.*p3*zeta/den;

                    case 12:
                      return 2.*p3*zeta/den;

                    case 13:
                      return -8.*p1*p3*xi/den2;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }

                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return -0.5*p4*(xi*eta - zeta + zeta*zeta - 2.*p1*xi)/den2;

                    case 1:
                      return 0.5*p2*(xi*eta + zeta - zeta*zeta - 2.*p1*xi)/den2;

                    case 2:
                      return 0.5*p2*(xi*eta - zeta + zeta*zeta + 2.*p3*xi)/den2;

                    case 3:
                      return -0.5*p4*(xi*eta + zeta - zeta*zeta + 2.*p3*xi)/den2;

                    case 4:
                      return 0.;

                    case 5:
                      return 2.*p2*p4*(eta - 2.*p1)/den2;

                    case 6:
                      return -2.*p2*xi*eta/den2;

                    case 7:
                      return 2.*p2*p4*(eta + 2.*p3)/den2;

                    case 8:
                      return 2.*p4*xi*eta/den2;

                    case 9:
                      return 2.*p4*zeta/den;

                    case 10:
                      return 2.*p2*zeta/den;

                    case 11:
                      return -2.*p2*zeta/den;

                    case 12:
                      return -2.*p4*zeta/den;

                    case 13:
                      return -8.*p2*p4*eta/den2;

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    }


                  // d/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.5*p1*(xi*eta - zeta + zeta*zeta)/den2
                          - 0.5*p4*(xi*eta - zeta + zeta*zeta)/den2
                          + p4*p1*(2.*zeta - 1)/den2
                          - 2.*p4*p1*(xi*eta - zeta + zeta*zeta)/den3;

                      case 1:
                        return 0.5*p2*(xi*eta + zeta - zeta*zeta)/den2
                          + 0.5*p1*(xi*eta + zeta - zeta*zeta)/den2
                          - p1*p2*(1 - 2.*zeta)/den2
                          + 2.*p1*p2*(xi*eta + zeta - zeta*zeta)/den3;

                      case 2:
                        return -0.5*p3*(xi*eta - zeta + zeta*zeta)/den2
                          - 0.5*p2*(xi*eta - zeta + zeta*zeta)/den2
                          + p2*p3*(2.*zeta - 1)/den2
                          - 2.*p2*p3*(xi*eta - zeta + zeta*zeta)/den3;

                      case 3:
                        return 0.5*p4*(xi*eta + zeta - zeta*zeta)/den2
                          + 0.5*p3*(xi*eta + zeta - zeta*zeta)/den2
                          - p3*p4*(1 - 2.*zeta)/den2
                          + 2.*p3*p4*(xi*eta + zeta - zeta*zeta)/den3;

                      case 4:
                        return 4.*zeta - 1.;

                      case 5:
                        return 2.*p4*p1*eta/den2
                          + 2.*p2*p4*eta/den2
                          + 2.*p1*p2*eta/den2
                          + 8.*p2*p1*p4*eta/den3;

                      case 6:
                        return -2.*p2*p3*xi/den2
                          - 2.*p1*p3*xi/den2
                          - 2.*p1*p2*xi/den2
                          - 8.*p1*p2*p3*xi/den3;

                      case 7:
                        return -2.*p3*p4*eta/den2
                          - 2.*p2*p4*eta/den2
                          - 2.*p2*p3*eta/den2
                          - 8.*p2*p3*p4*eta/den3;

                      case 8:
                        return 2.*p4*p1*xi/den2
                          + 2.*p1*p3*xi/den2
                          + 2.*p3*p4*xi/den2
                          + 8.*p3*p4*p1*xi/den3;

                      case 9:
                        return 2.*p4*zeta/den
                          + 2.*p1*zeta/den
                          - 4.*p1*p4/den
                          + 4.*p1*p4*zeta/den2;

                      case 10:
                        return 2.*p1*zeta/den
                          + 2.*p2*zeta/den
                          - 4.*p2*p1/den
                          + 4.*p2*p1*zeta/den2;

                      case 11:
                        return 2.*p2*zeta/den
                          + 2.*p3*zeta/den
                          - 4.*p3*p2/den
                          + 4.*p3*p2*zeta/den2;

                      case 12:
                        return 2.*p3*zeta/den
                          + 2.*p4*zeta/den
                          - 4.*p4*p3/den
                          + 4.*p4*p3*zeta/den2;

                      case 13:
                        return -8.*p2*p3*p4/den2
                          - 8.*p3*p4*p1/den2
                          - 8.*p2*p1*p4/den2
                          - 8.*p1*p2*p3/den2
                          - 32.*p1*p2*p3*p4/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }

    case THIRD:
      {
        switch (type)
          {
            // quadratic Lagrange shape functions with a cubic bubble
          case TET14:
            {
              libmesh_assert_less (i, 14);
              return libMesh::detail::fe_lagrange_tet14_shape_deriv(i, j, p(0), p(1), p(2));
            }

          case PRISM20:
            {
              libmesh_assert_less (i, 20);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              const Real mainval = FE<3,LAGRANGE>::shape_deriv(PRISM21, THIRD, i, j, p);

              if (i0[i] != 2)
                return mainval;

              Real bubbleval = 0;

              switch (j)
                {
                  // d()/dxi
                case 0:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, 6, 0, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d()/deta
                case 1:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, 6, 1, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d()/dzeta
                case 2:
                  bubbleval =
                    FE<2,LAGRANGE>::shape(TRI7, THIRD, 6, p2d)*
                    fe_lagrange_1D_quadratic_shape_deriv(2, 0, p1d);
                  break;

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }

              if (i < 12) // vertices
                return mainval - bubbleval / 9;

              return mainval + bubbleval * (Real(4) / 9);
            }

          case PRISM21:
            {
              libmesh_assert_less (i, 21);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 6, 6};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/deta
                case 1:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (FE<2,LAGRANGE>::shape(TRI7, THIRD, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

          case PYRAMID18:
            {
              return fe_fdm_deriv(type, order, elem, i, j, p, fe_lagrange_3D_shape<T>);
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << order);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(type, order, elem, i, j, p);
  libmesh_not_implemented();
#endif
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <FEFamily T>
Real fe_lagrange_3D_shape_second_deriv(const ElemType type,
                                       const Order order,
                                       const Elem * elem,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p)
{
#if LIBMESH_DIM == 3

  libmesh_assert_less (j, 6);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
            // Linear tets have all second derivatives = 0
          case TET4:
          case TET10:
          case TET14:
            {
              return 0.;
            }

            // The following elements use either tensor product or
            // rational basis functions, and therefore probably have
            // second derivatives, but we have not implemented them
            // yet...
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              libmesh_assert_less (i, 6);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

              switch (j)
                {
                  // All repeated second derivatives and the xi-eta derivative are zero on PRISMs
                case 0: // d^2()/dxi^2
                case 1: // d^2()/dxideta
                case 2: // d^2()/deta^2
                case 5: // d^2()/dzeta^2
                  {
                    return 0.;
                  }

                case 3: // d^2()/dxidzeta
                  return (FE<2,LAGRANGE>::shape_deriv(TRI3,  FIRST, i1[i], 0, p2d)*
                          fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                case 4: // d^2()/detadzeta
                  return (FE<2,LAGRANGE>::shape_deriv(TRI3,  FIRST, i1[i], 1, p2d)*
                          fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
            {
              libmesh_assert_less (i, 5);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              switch (j)
                {
                  // xi-xi and eta-eta derivatives are all zero for PYRAMID5.
                case 0: // d^2()/dxi^2
                case 2: // d^2()/deta^2
                  return 0.;

                case 1: // d^2()/dxideta
                  {
                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.25/(1. - zeta + eps);
                      case 1:
                      case 3:
                        return -0.25/(1. - zeta + eps);
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 3: // d^2()/dxidzeta
                  {
                    Real den = (1. - zeta + eps)*(1. - zeta + eps);

                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.25*eta/den;
                      case 1:
                      case 3:
                        return -0.25*eta/den;
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    Real den = (1. - zeta + eps)*(1. - zeta + eps);

                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.25*xi/den;
                      case 1:
                      case 3:
                        return -0.25*xi/den;
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 5: // d^2()/dzeta^2
                  {
                    Real den = (1. - zeta + eps)*(1. - zeta + eps)*(1. - zeta + eps);

                    switch (i)
                      {
                      case 0:
                      case 2:
                        return 0.5*xi*eta/den;
                      case 1:
                      case 3:
                        return -0.5*xi*eta/den;
                      case 4:
                        return 0.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // Trilinear shape functions on HEX8s have nonzero mixed second derivatives
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 8);
              return libMesh::detail::fe_lagrange_hex8_shape_second_deriv(i, j, p(0), p(1), p(2));
            }

          // All second derivatives for piecewise-linear polyhedra are
          // zero or dirac-type distributions, but we can't put the
          // latter in a Real, so beware when integrating...
          case C0POLYHEDRON:
            {
              return 0.;
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }

      }

      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (type)
          {

            // serendipity hexahedral quadratic shape functions
          case HEX20:
            {
              libmesh_assert_less (i, 20);
              return libMesh::detail::fe_lagrange_hex20_shape_second_deriv(i, j, p(0), p(1), p(2));
            }

            // triquadratic hexahedral shape functions
          case HEX8:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case HEX27:
            {
              libmesh_assert_less (i, 27);
              return libMesh::detail::fe_lagrange_hex27_shape_second_deriv(i, j, p(0), p(1), p(2));
            }

            // quadratic tetrahedral shape functions
          case TET4:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case TET10:
          case TET14:
            {
              libmesh_assert_less (i, 10);
              return libMesh::detail::fe_lagrange_tet10_shape_second_deriv(i, j);
            }



            // "serendipity" prism
          case PRISM15:
            {
              libmesh_assert_less (i, 15);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return 2.*(1. - zeta);
                      case 2:
                      case 5:
                      case 7:
                      case 8:
                      case 9:
                      case 10:
                      case 11:
                      case 13:
                      case 14:
                        return 0.;
                      case 3:
                      case 4:
                        return 2.*(1. + zeta);
                      case 6:
                        return 4.*(zeta - 1);
                      case 12:
                        return -4.*(1. + zeta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/dxideta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                      case 7:
                        return 2.*(1. - zeta);
                      case 1:
                      case 2:
                      case 4:
                      case 5:
                      case 9:
                      case 10:
                      case 11:
                        return 0.;
                      case 3:
                      case 13:
                        return 2.*(1. + zeta);
                      case 6:
                      case 8:
                        return 2.*(zeta - 1.);
                      case 12:
                      case 14:
                        return -2.*(1. + zeta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/deta^2
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                      case 2:
                        return 2.*(1. - zeta);
                      case 1:
                      case 4:
                      case 6:
                      case 7:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                      case 13:
                        return 0.;
                      case 3:
                      case 5:
                        return 2.*(1. + zeta);
                      case 8:
                        return 4.*(zeta - 1.);
                      case 14:
                        return -4.*(1. + zeta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/dxidzeta
                case 3:
                  {
                    switch(i)
                      {
                      case 0:
                        return 1.5 - zeta - 2.*xi - 2.*eta;
                      case 1:
                        return 0.5 + zeta - 2.*xi;
                      case 2:
                      case 5:
                      case 11:
                        return 0.;
                      case 3:
                        return -1.5 - zeta + 2.*xi + 2.*eta;
                      case 4:
                        return -0.5 + zeta + 2.*xi;
                      case 6:
                        return 4.*xi + 2.*eta - 2.;
                      case 7:
                        return -2.*eta;
                      case 8:
                        return 2.*eta;
                      case 9:
                        return 2.*zeta;
                      case 10:
                        return -2.*zeta;
                      case 12:
                        return -4.*xi - 2.*eta + 2.;
                      case 13:
                        return 2.*eta;
                      case 14:
                        return -2.*eta;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/detadzeta
                case 4:
                  {
                    switch(i)
                      {
                      case 0:
                        return 1.5 - zeta - 2.*xi - 2.*eta;
                      case 1:
                      case 4:
                      case 10:
                        return 0.;
                      case 2:
                        return .5 + zeta - 2.*eta;
                      case 3:
                        return -1.5 - zeta + 2.*xi + 2.*eta;
                      case 5:
                        return -.5 + zeta + 2.*eta;
                      case 6:
                        return 2.*xi;
                      case 7:
                        return -2.*xi;
                      case 8:
                        return 2.*xi + 4.*eta - 2.;
                      case 9:
                        return 2.*zeta;
                      case 11:
                        return -2.*zeta;
                      case 12:
                        return -2.*xi;
                      case 13:
                        return 2.*xi;
                      case 14:
                        return -2.*xi - 4.*eta + 2.;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                  // d^2()/dzeta^2
                case 5:
                  {
                    switch(i)
                      {
                      case 0:
                      case 3:
                        return 1. - xi - eta;
                      case 1:
                      case 4:
                        return xi;
                      case 2:
                      case 5:
                        return eta;
                      case 6:
                      case 7:
                      case 8:
                      case 12:
                      case 13:
                      case 14:
                        return 0.;
                      case 9:
                        return 2.*xi + 2.*eta - 2.;
                      case 10:
                        return -2.*xi;
                      case 11:
                        return -2.*eta;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }



            // quadratic prism shape functions
          case PRISM6:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              libmesh_assert_less (i, 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  return (FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/dxideta
                case 1:
                  return (FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/deta^2
                case 2:
                  return (FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, i1[i], 2, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/dxidzeta
                case 3:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                  // d^2()/detadzeta
                case 4:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                  // d^2()/dzeta^2
                case 5:
                  return (FE<2,LAGRANGE>::shape(TRI6,  SECOND, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }


            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID5:
            libmesh_assert_msg(T == L2_LAGRANGE,
                               "High order on first order elements only supported for L2 families");
            libmesh_fallthrough();
          case PYRAMID14:
          case PYRAMID18:
            {
              libmesh_assert_less (i, 14);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              Real
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              Real
                den = (-1. + zeta + eps),
                den2 = den*den,
                den3 = den2*den,
                den4 = den2*den2;

              // These terms are used in several of the derivatives
              Real
                numer_mp = xi*eta - zeta + zeta*zeta,
                numer_pm = xi*eta + zeta - zeta*zeta;

              switch (j)
                {
                case 0: // d^2()/dxi^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return -p1*eta/den2;
                      case 2:
                      case 3:
                        return p3*eta/den2;
                      case 4:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;
                      case 5:
                        return 2.*p1*eta/den2;
                      case 6:
                      case 8:
                        return 4.*p1*p3/den2;
                      case 7:
                        return -2.*p3*eta/den2;
                      case 13:
                        return -8.*p1*p3/den2;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 1: // d^2()/dxideta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p1*xi/den2
                          - 0.5*p4*eta/den2
                          + p4*p1/den2;

                      case 1:
                        return 0.25*numer_pm/den2
                          - 0.5*p1*xi/den2
                          + 0.5*p2*eta/den2
                          - p1*p2/den2;

                      case 2:
                        return 0.25*numer_mp/den2
                          + 0.5*p3*xi/den2
                          + 0.5*p2*eta/den2
                          + p2*p3/den2;

                      case 3:
                        return 0.25*numer_pm/den2
                          + 0.5*p3*xi/den2
                          - 0.5*p4*eta/den2
                          - p3*p4/den2;

                      case 4:
                        return 0.;

                      case 5:
                        return p4*eta/den2
                          - 2.*p4*p1/den2
                          - p2*eta/den2
                          + 2.*p1*p2/den2;

                      case 6:
                        return -p3*xi/den2
                          + p1*xi/den2
                          - 2.*p2*p3/den2
                          + 2.*p1*p2/den2;

                      case 7:
                        return p4*eta/den2
                          + 2.*p3*p4/den2
                          - p2*eta/den2
                          - 2.*p2*p3/den2;

                      case 8:
                        return -p3*xi/den2
                          + p1*xi/den2
                          - 2.*p4*p1/den2
                          + 2.*p3*p4/den2;

                      case 9:
                      case 11:
                        return -zeta/den;

                      case 10:
                      case 12:
                        return zeta/den;

                      case 13:
                        return 4.*p4*p1/den2
                          - 4.*p3*p4/den2
                          + 4.*p2*p3/den2
                          - 4.*p1*p2/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 2: // d^2()/deta^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 3:
                        return -p4*xi/den2;
                      case 1:
                      case 2:
                        return p2*xi/den2;
                      case 4:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;
                      case 5:
                      case 7:
                        return 4.*p2*p4/den2;
                      case 6:
                        return -2.*p2*xi/den2;
                      case 8:
                        return 2.*p4*xi/den2;
                      case 13:
                        return -8.*p2*p4/den2;
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 3: // d^2()/dxidzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p1*(2.*zeta - 1.)/den2
                          + p1*numer_mp/den3
                          - 0.5*p1*eta/den2
                          - 0.5*p4*eta/den2
                          - 2.*p4*p1*eta/den3;

                      case 1:
                        return 0.25*numer_pm/den2
                          - 0.5*p1*(1 - 2.*zeta)/den2
                          + p1*numer_pm/den3
                          + 0.5*p2*eta/den2
                          + 0.5*p1*eta/den2
                          + 2.*p1*p2*eta/den3;

                      case 2:
                        return -0.25*numer_mp/den2
                          + 0.5*p3*(2.*zeta - 1.)/den2
                          - p3*numer_mp/den3
                          - 0.5*p3*eta/den2
                          - 0.5*p2*eta/den2
                          - 2.*p2*p3*eta/den3;

                      case 3:
                        return -0.25*numer_pm/den2
                          + 0.5*p3*(1 - 2.*zeta)/den2
                          - p3*numer_pm/den3
                          + 0.5*p4*eta/den2
                          + 0.5*p3*eta/den2
                          + 2.*p3*p4*eta/den3;

                      case 4:
                        return 0.;

                      case 5:
                        return p4*eta/den2
                          + 4.*p4*p1*eta/den3
                          - p2*eta/den2
                          - 4.*p1*p2*eta/den3;

                      case 6:
                        return -p3*xi/den2
                          - p1*xi/den2
                          - 4.*p1*p3*xi/den3
                          - 2.*p2*p3/den2
                          - 2.*p1*p3/den2
                          - 2.*p1*p2/den2
                          - 8.*p1*p2*p3/den3;

                      case 7:
                        return -p4*eta/den2
                          - 4.*p3*p4*eta/den3
                          + p2*eta/den2
                          + 4.*p2*p3*eta/den3;

                      case 8:
                        return -p3*xi/den2
                          - p1*xi/den2
                          - 4.*p1*p3*xi/den3
                          + 2.*p4*p1/den2
                          + 2.*p1*p3/den2
                          + 2.*p3*p4/den2
                          + 8.*p3*p4*p1/den3;

                      case 9:
                        return -zeta/den
                          + 2.*p1/den
                          - 2.*p1*zeta/den2;

                      case 10:
                        return zeta/den
                          - 2.*p1/den
                          + 2.*p1*zeta/den2;

                      case 11:
                        return zeta/den
                          - 2.*p3/den
                          + 2.*p3*zeta/den2;

                      case 12:
                        return -zeta/den
                          + 2.*p3/den
                          - 2.*p3*zeta/den2;

                      case 13:
                        return -4.*p4*p1/den2
                          - 4.*p3*p4/den2
                          - 16.*p3*p4*p1/den3
                          + 4.*p2*p3/den2
                          + 4.*p1*p2/den2
                          + 16.*p1*p2*p3/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p4*(2.*zeta - 1.)/den2
                          + p4*numer_mp/den3
                          - 0.5*p1*xi/den2
                          - 0.5*p4*xi/den2
                          - 2.*p4*p1*xi/den3;

                      case 1:
                        return -0.25*numer_pm/den2
                          + 0.5*p2*(1. - 2.*zeta)/den2
                          - p2*numer_pm/den3
                          + 0.5*p2*xi/den2
                          + 0.5*p1*xi/den2
                          + 2.*p1*p2*xi/den3;

                      case 2:
                        return -0.25*numer_mp/den2
                          + 0.5*p2*(2.*zeta - 1.)/den2
                          - p2*numer_mp/den3
                          - 0.5*p3*xi/den2
                          - 0.5*p2*xi/den2
                          - 2.*p2*p3*xi/den3;

                      case 3:
                        return 0.25*numer_pm/den2
                          - 0.5*p4*(1. - 2.*zeta)/den2
                          + p4*numer_pm/den3
                          + 0.5*p4*xi/den2
                          + 0.5*p3*xi/den2
                          + 2.*p3*p4*xi/den3;

                      case 4:
                        return 0.;

                      case 5:
                        return -p4*eta/den2
                          - p2*eta/den2
                          - 4.*p2*p4*eta/den3
                          + 2.*p4*p1/den2
                          + 2.*p2*p4/den2
                          + 2.*p1*p2/den2
                          + 8.*p2*p1*p4/den3;

                      case 6:
                        return p3*xi/den2
                          + 4.*p2*p3*xi/den3
                          - p1*xi/den2
                          - 4.*p1*p2*xi/den3;

                      case 7:
                        return -p4*eta/den2
                          - p2*eta/den2
                          - 4.*p2*p4*eta/den3
                          - 2.*p3*p4/den2
                          - 2.*p2*p4/den2
                          - 2.*p2*p3/den2
                          - 8.*p2*p3*p4/den3;

                      case 8:
                        return p1*xi/den2
                          + 4.*p4*p1*xi/den3
                          - p3*xi/den2
                          - 4.*p3*p4*xi/den3;

                      case 9:
                        return -zeta/den
                          + 2.*p4/den
                          - 2.*p4*zeta/den2;

                      case 10:
                        return -zeta/den
                          + 2.*p2/den
                          - 2.*p2*zeta/den2;

                      case 11:
                        return zeta/den
                          - 2.*p2/den
                          + 2.*p2*zeta/den2;

                      case 12:
                        return zeta/den
                          - 2.*p4/den
                          + 2.*p4*zeta/den2;

                      case 13:
                        return 4.*p3*p4/den2
                          + 4.*p2*p3/den2
                          + 16.*p2*p3*p4/den3
                          - 4.*p4*p1/den2
                          - 4.*p1*p2/den2
                          - 16.*p2*p1*p4/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 5: // d^2()/dzeta^2
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.5*numer_mp/den2
                          - p1*(2.*zeta - 1.)/den2
                          + 2.*p1*numer_mp/den3
                          - p4*(2.*zeta - 1.)/den2
                          + 2.*p4*numer_mp/den3
                          + 2.*p4*p1/den2
                          - 4.*p4*p1*(2.*zeta - 1.)/den3
                          + 6.*p4*p1*numer_mp/den4;

                      case 1:
                        return -0.5*numer_pm/den2
                          + p2*(1 - 2.*zeta)/den2
                          - 2.*p2*numer_pm/den3
                          + p1*(1 - 2.*zeta)/den2
                          - 2.*p1*numer_pm/den3
                          + 2.*p1*p2/den2
                          + 4.*p1*p2*(1 - 2.*zeta)/den3
                          - 6.*p1*p2*numer_pm/den4;

                      case 2:
                        return 0.5*numer_mp/den2
                          - p3*(2.*zeta - 1.)/den2
                          + 2.*p3*numer_mp/den3
                          - p2*(2.*zeta - 1.)/den2
                          + 2.*p2*numer_mp/den3
                          + 2.*p2*p3/den2
                          - 4.*p2*p3*(2.*zeta - 1.)/den3
                          + 6.*p2*p3*numer_mp/den4;

                      case 3:
                        return -0.5*numer_pm/den2
                          + p4*(1 - 2.*zeta)/den2
                          - 2.*p4*numer_pm/den3
                          + p3*(1 - 2.*zeta)/den2
                          - 2.*p3*numer_pm/den3
                          + 2.*p3*p4/den2
                          + 4.*p3*p4*(1 - 2.*zeta)/den3
                          - 6.*p3*p4*numer_pm/den4;

                      case 4:
                        return 4.;

                      case 5:
                        return -2.*p1*eta/den2
                          - 2.*p4*eta/den2
                          - 8.*p4*p1*eta/den3
                          - 2.*p2*eta/den2
                          - 8.*p2*p4*eta/den3
                          - 8.*p1*p2*eta/den3
                          - 24.*p2*p1*p4*eta/den4;

                      case 6:
                        return 2.*p3*xi/den2
                          + 2.*p2*xi/den2
                          + 8.*p2*p3*xi/den3
                          + 2.*p1*xi/den2
                          + 8.*p1*p3*xi/den3
                          + 8.*p1*p2*xi/den3
                          + 24.*p1*p2*p3*xi/den4;

                      case 7:
                        return 2.*p4*eta/den2
                          + 2.*p3*eta/den2
                          + 8.*p3*p4*eta/den3
                          + 2.*p2*eta/den2
                          + 8.*p2*p4*eta/den3
                          + 8.*p2*p3*eta/den3
                          + 24.*p2*p3*p4*eta/den4;

                      case 8:
                        return -2.*p1*xi/den2
                          - 2.*p4*xi/den2
                          - 8.*p4*p1*xi/den3
                          - 2.*p3*xi/den2
                          - 8.*p1*p3*xi/den3
                          - 8.*p3*p4*xi/den3
                          - 24.*p3*p4*p1*xi/den4;

                      case 9:
                        return -2.*zeta/den
                          + 4.*p4/den
                          - 4.*p4*zeta/den2
                          + 4.*p1/den
                          - 4.*p1*zeta/den2
                          + 8.*p4*p1/den2
                          - 8.*p1*p4*zeta/den3;

                      case 10:
                        return -2.*zeta/den
                          + 4.*p1/den
                          - 4.*p1*zeta/den2
                          + 4.*p2/den
                          - 4.*p2*zeta/den2
                          + 8.*p1*p2/den2
                          - 8.*p2*p1*zeta/den3;

                      case 11:
                        return -2.*zeta/den
                          + 4.*p2/den
                          - 4.*p2*zeta/den2
                          + 4.*p3/den
                          - 4.*p3*zeta/den2
                          + 8.*p2*p3/den2
                          - 8.*p3*p2*zeta/den3;

                      case 12:
                        return -2.*zeta/den
                          + 4.*p3/den
                          - 4.*p3*zeta/den2
                          + 4.*p4/den
                          - 4.*p4*zeta/den2
                          + 8.*p3*p4/den2
                          - 8.*p4*p3*zeta/den3;

                      case 13:
                        return 8.*p3*p4/den2
                          + 8.*p2*p4/den2
                          + 8.*p2*p3/den2
                          + 32.*p2*p3*p4/den3
                          + 8.*p4*p1/den2
                          + 8.*p1*p3/den2
                          + 32.*p3*p4*p1/den3
                          + 8.*p1*p2/den2
                          + 32.*p2*p1*p4/den3
                          + 32.*p1*p2*p3/den3
                          + 96.*p1*p2*p3*p4/den4;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              libmesh_assert_less (i, 13);

              const Real xi   = p(0);
              const Real eta  = p(1);
              const Real zeta = p(2);
              const Real eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              Real
                den = (-1. + zeta + eps),
                den2 = den*den,
                den3 = den2*den,
                xi2 = xi*xi,
                eta2 = eta*eta,
                zeta2 = zeta*zeta,
                zeta3 = zeta2*zeta;

              switch (j)
                {
                case 0: // d^2()/dxi^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 1:
                        return 0.5*(-1. + zeta + eta)/den;

                      case 2:
                      case 3:
                        return 0.5*(-1. + zeta - eta)/den;

                      case 4:
                      case 6:
                      case 8:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;

                      case 5:
                        return (1. - eta - zeta)/den;

                      case 7:
                        return (1. + eta - zeta)/den;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 1: // d^2()/dxideta
                  {
                    switch(i)
                      {
                      case 0:
                        return  0.25*(-1. + 2.*zeta + 2.*xi + 2.*eta)/den;

                      case 1:
                        return -0.25*(-1. + 2.*zeta - 2.*xi + 2.*eta)/den;

                      case 2:
                        return -0.25*(1. - 2.*zeta + 2.*xi + 2.*eta)/den;

                      case 3:
                        return  0.25*(1. - 2.*zeta - 2.*xi + 2.*eta)/den;

                      case 4:
                        return 0.;

                      case 5:
                        return -xi/den;

                      case 6:
                        return eta/den;

                      case 7:
                        return xi/den;

                      case 8:
                        return -eta/den;

                      case 9:
                        return -zeta/den;

                      case 10:
                        return zeta/den;

                      case 11:
                        return -zeta/den;

                      case 12:
                        return zeta/den;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 2: // d^2()/deta^2
                  {
                    switch(i)
                      {
                      case 0:
                      case 3:
                        return 0.5*(-1. + zeta + xi)/den;

                      case 1:
                      case 2:
                        return 0.5*(-1. + zeta - xi)/den;

                      case 4:
                      case 5:
                      case 7:
                      case 9:
                      case 10:
                      case 11:
                      case 12:
                        return 0.;

                      case 6:
                        return (1. + xi - zeta)/den;

                      case 8:
                        return (1. - xi - zeta)/den;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }


                case 3: // d^2()/dxidzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.25*(-1. + 2.*zeta - zeta2 + eta + 2.*eta*xi + eta2)/den2;

                      case 1:
                        return 0.25*(-1. + 2.*zeta - zeta2 + eta - 2.*eta*xi + eta2)/den2;

                      case 2:
                        return 0.25*(-1. + 2.*zeta - zeta2 - eta + 2.*eta*xi + eta2)/den2;

                      case 3:
                        return -0.25*(-1. + 2.*zeta - zeta2 - eta - 2.*eta*xi + eta2)/den2;

                      case 4:
                        return 0.;

                      case 5:
                        return eta*xi/den2;

                      case 6:
                        return -0.5*(1. + zeta2 + eta2 - 2.*zeta)/den2;

                      case 7:
                        return -eta*xi/den2;

                      case 8:
                        return 0.5*(1. + zeta2 + eta2 - 2.*zeta)/den2;

                      case 9:
                        return (-1. - zeta2 + eta + 2.*zeta)/den2;

                      case 10:
                        return -(-1. - zeta2 + eta + 2.*zeta)/den2;

                      case 11:
                        return (1. + zeta2 + eta - 2.*zeta)/den2;

                      case 12:
                        return -(1. + zeta2 + eta - 2.*zeta)/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.25*(-1. + 2.*zeta - zeta2 + xi + 2.*eta*xi + xi2)/den2;

                      case 1:
                        return 0.25*(1. - 2.*zeta + zeta2 + xi + 2.*eta*xi - xi2)/den2;

                      case 2:
                        return 0.25*(-1. + 2.*zeta - zeta2 - xi + 2.*eta*xi + xi2)/den2;

                      case 3:
                        return -0.25*(1. - 2.*zeta + zeta2 - xi + 2.*eta*xi - xi2)/den2;

                      case 4:
                        return 0.;

                      case 5:
                        return 0.5*(1. + xi2 + zeta2 - 2.*zeta)/den2;

                      case 6:
                        return -eta*xi/den2;

                      case 7:
                        return -0.5*(1. + xi2 + zeta2 - 2.*zeta)/den2;

                      case 8:
                        return eta*xi/den2;

                      case 9:
                        return (-1. - zeta2 + xi + 2.*zeta)/den2;

                      case 10:
                        return -(1. + zeta2 + xi - 2.*zeta)/den2;

                      case 11:
                        return (1. + zeta2 + xi - 2.*zeta)/den2;

                      case 12:
                        return -(-1. - zeta2 + xi + 2.*zeta)/den2;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                case 5: // d^2()/dzeta^2
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.5*(xi + eta + 1.)*eta*xi/den3;

                      case 1:
                        return -0.5*(eta - xi + 1.)*eta*xi/den3;

                      case 2:
                        return -0.5*(xi + eta - 1.)*eta*xi/den3;

                      case 3:
                        return 0.5*(eta - xi - 1.)*eta*xi/den3;

                      case 4:
                        return 4.;

                      case 5:
                        return -(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi2)/den3;

                      case 6:
                        return (-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta2*xi)/den3;

                      case 7:
                        return (-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi2)/den3;

                      case 8:
                        return -(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta2*xi)/den3;

                      case 9:
                        return -2.*(-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi)/den3;

                      case 10:
                        return 2.*(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi)/den3;

                      case 11:
                        return -2.*(-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi)/den3;

                      case 12:
                        return 2.*(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi)/den3;

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  }

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }

    case THIRD:
      {
        switch (type)
          {
            // quadratic Lagrange shape functions with a cubic bubble
          case TET14:
            {
              libmesh_assert_less (i, 14);
              return libMesh::detail::fe_lagrange_tet14_shape_second_deriv(i, j, p(0), p(1), p(2));
            }

          case PRISM20:
            {
              libmesh_assert_less (i, 20);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1};

              const Real mainval = FE<3,LAGRANGE>::shape_second_deriv(PRISM21, THIRD, i, j, p);

              if (i0[i] != 2)
                return mainval;

              Real bubbleval = 0;

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_second_deriv(TRI7, THIRD, 6, 0, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d^2()/dxideta
                case 1:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_second_deriv(TRI7, THIRD, 6, 1, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d^2()/deta^2
                case 2:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_second_deriv(TRI7, THIRD, 6, 2, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d^2()/dxidzeta
                case 3:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, 6, 0, p2d)*
                    fe_lagrange_1D_quadratic_shape_deriv(2, 0, p1d);
                  break;

                  // d^2()/detadzeta
                case 4:
                  bubbleval =
                    FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, 6, 1, p2d)*
                    fe_lagrange_1D_quadratic_shape_deriv(2, 0, p1d);
                  break;

                  // d^2()/dzeta^2
                case 5:
                  bubbleval =
                    FE<2,LAGRANGE>::shape(TRI7, THIRD, 6, p2d)*
                    fe_lagrange_1D_quadratic_shape_second_deriv(2, 0, p1d);
                  break;

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }

              if (i < 12) // vertices
                return mainval - bubbleval / 9;

              return mainval + bubbleval * (Real(4) / 9);
            }

          case PRISM21:
            {
              libmesh_assert_less (i, 21);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge

              Point p2d(p(0),p(1));
              Real p1d = p(2);

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 6, 6};

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  return (FE<2,LAGRANGE>::shape_second_deriv(TRI7, THIRD, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/dxideta
                case 1:
                  return (FE<2,LAGRANGE>::shape_second_deriv(TRI7, THIRD, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/deta^2
                case 2:
                  return (FE<2,LAGRANGE>::shape_second_deriv(TRI7, THIRD, i1[i], 2, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d^2()/dxidzeta
                case 3:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                  // d^2()/detadzeta
                case 4:
                  return (FE<2,LAGRANGE>::shape_deriv(TRI7, THIRD, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                  // d^2()/dzeta^2
                case 5:
                  return (FE<2,LAGRANGE>::shape(TRI7, THIRD, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, p1d));

                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }
            }

          case PYRAMID18:
            {
              return fe_fdm_second_deriv(type, order, elem, i, j, p,
                                         fe_lagrange_3D_shape_deriv<T>);
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(type));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << order);
    }

#else // LIBMESH_DIM != 3
  libmesh_ignore(type, order, elem, i, j, p);
  libmesh_not_implemented();
#endif
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES


} // anonymous namespace
