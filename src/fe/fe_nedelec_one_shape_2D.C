// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// An excellent discussion of Nedelec shape functions is given in
// https://www.dealii.org/reports/nedelec/nedelec.pdf
// An excellent summary of Nedelec shape functions is also given in
// https://defelement.com/elements/nedelec1.html
template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  libmesh_assert_less(i, n_dofs(elem->type(), total_order));

  const char sign = i >= total_order * elem->n_edges() || elem->point(i / total_order) > elem->point((i / total_order + 1) % elem->n_vertices()) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / total_order * 2 + 1) * total_order - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape functions
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              switch(ii)
                {
                case 0:
                  return sign * RealGradient( -0.25*(1.0-eta), 0.0 );
                case 1:
                  return sign * RealGradient( 0.0, -0.25*(1.0+xi) );
                case 2:
                  return sign * RealGradient( 0.25*(1.0+eta), 0.0 );
                case 3:
                  return sign * RealGradient( 0.0, -0.25*(xi-1.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( -1.0+eta, -xi );
                case 1:
                  return sign * RealGradient( eta, -xi );
                case 2:
                  return sign * RealGradient( eta, -xi+1.0 );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // quadratic Nedelec (first kind) shape functions
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              const Real x = 0.5 * (xi + 1.0);
              const Real y = 0.5 * (eta + 1.0);

              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 0.5*(-18.0*x*y*y+24.0*x*y-6.0*x+12.0*y*y-16.0*y+4.0), 0.0 );
                case 1:
                  return sign * RealGradient( 0.5*( 18.0*x*y*y-24.0*x*y+6.0*x-6.0*y*y+8.0*y-2.0), 0.0 );
                case 2:
                  return sign * RealGradient( 0.0, x*(-9.0*x*y+6.0*x+6.0*y-4.0) );
                case 3:
                  return sign * RealGradient( 0.0, x*( 9.0*x*y-3.0*x-6.0*y+2.0) );
                case 4:
                  return sign * RealGradient( y*(-9.0*x*y+6.0*x+3.0*y-2.0), 0.0 );
                case 5:
                  return sign * RealGradient( y*( 9.0*x*y-6.0*x-6.0*y+4.0), 0.0 );
                case 6:
                  return sign * RealGradient( 0.0, 0.5*(-18.0*x*x*y+6.0*x*x+24.0*x*y-8.0*x-6.0*y+2.0) );
                case 7:
                  return sign * RealGradient( 0.0, 0.5*( 18.0*x*x*y-12.0*x*x-24.0*x*y+16.0*x+6.0*y-4.0) );
                case 8:
                  return RealGradient( 0.0, 3.0*x*(3*x*y-2.0*x-3.0*y+2.0) );
                case 9:
                  return RealGradient( 3.0*y*(-3.0*x*y+3.0*x+2.0*y-2.0), 0.0 );
                case 10:
                  return RealGradient( 3.0*y*(3.0*x*y-3.0*x-y+1.0), 0.0 );
                case 11:
                  return RealGradient( 0.0, 3.0*x*(-3.0*x*y+x+3.0*y-1.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 8.0*xi*eta-6.0*xi+8.0*eta*eta-12.0*eta+4.0, 2.0*xi*(-4.0*xi-4.0*eta+3.0) );
                case 1:
                  return sign * RealGradient( -8.0*xi*eta+6.0*xi+2.0*eta-2.0, 4.0*xi*(2.0*xi-1.0) );
                case 2:
                  return sign * RealGradient( 2.0*eta*(1.0-4.0*xi), 4.0*xi*(2.0*xi-1.0) );
                case 3:
                  return sign * RealGradient( 4.0*eta*(1.0-2.0*eta), 2.0*xi*(4.0*eta-1.0) );
                case 4:
                  return sign * RealGradient( 4.0*eta*(1.0-2.0*eta), 8.0*xi*eta-2.0*xi-6.0*eta+2.0 );
                case 5:
                  return sign * RealGradient( 2.0*eta*(4.0*xi+4*eta-3.0), -8.0*xi*xi-8.0*xi*eta+12.0*xi+6.0*eta-4.0 );
                case 6:
                  return RealGradient( 8.0*eta*(-xi-2.0*eta+2.0), 8.0*xi*(xi+2.0*eta-1.0) );
                case 7:
                  return RealGradient( 8.0*eta*(2.0*xi+eta-1.0), 8.0*xi*(-2.0*xi-eta+2.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // cubic Nedelec (first kind) shape functions
    case THIRD:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient(-20.25*eta - 9.0*xi + 162.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 33.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 81.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 22.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 24.75 + 20.25*((eta + 1)*(eta + 1)) + 16.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.75*((xi + 1)*(xi + 1)) - 5.625*(eta + 1)*(eta + 1)*(eta + 1) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                case 1:
                  return sign * RealGradient(3.375*eta + 3.75*xi - 67.5*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 16.875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 33.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 9.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 6.375 - 3.375*(eta + 1)*(eta + 1) - 8.4375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 1.875*(xi + 1)*(xi + 1) + 0.9375*((eta + 1)*(eta + 1)*(eta + 1)) + 2.34375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)), 0);
                case 2:
                  return sign * RealGradient(-6.75*eta - 6.0*xi + 108.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 33.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 54.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 15.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 11.25 + 6.75*((eta + 1)*(eta + 1)) + 16.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.75*((xi + 1)*(xi + 1)) - 1.875*(eta + 1)*(eta + 1)*(eta + 1) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                case 3:
                  return sign * RealGradient(0, 6.75*xi - 54.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 54.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 22.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 11.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 6.75 - 11.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 13.5*(xi + 1)*(xi + 1) + 5.625*((xi + 1)*(xi + 1)*(xi + 1)) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)));
                case 4:
                  return sign * RealGradient(0, -1.125*xi + 22.5*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 22.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 9.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 5.625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 1.125 + 5.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 2.25*((xi + 1)*(xi + 1)) - 0.9375*(xi + 1)*(xi + 1)*(xi + 1) - 2.34375*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1));
                case 5:
                  return sign * RealGradient(0, 2.25*xi - 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 36.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 15.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 11.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 2.25 - 11.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 4.5*(xi + 1)*(xi + 1) + 1.875*((xi + 1)*(xi + 1)*(xi + 1)) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)));
                case 6:
                  return sign * RealGradient(-2.25*eta + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 11.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 36.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 15.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 2.25 + 4.5*((eta + 1)*(eta + 1)) + 11.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1.875*(eta + 1)*(eta + 1)*(eta + 1) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                case 7:
                  return sign * RealGradient(1.125*eta - 22.5*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 5.625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 22.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 9.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1.125 - 2.25*(eta + 1)*(eta + 1) - 5.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 0.9375*((eta + 1)*(eta + 1)*(eta + 1)) + 2.34375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)), 0);
                case 8:
                  return sign * RealGradient(-6.75*eta + 54.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 11.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 54.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 22.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.75 + 13.5*((eta + 1)*(eta + 1)) + 11.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 5.625*(eta + 1)*(eta + 1)*(eta + 1) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                case 9:
                  return sign * RealGradient(0, 6.0*eta + 6.75*xi - 108.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 54.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 15.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 33.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 11.25 - 3.75*(eta + 1)*(eta + 1) - 16.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 6.75*(xi + 1)*(xi + 1) + 1.875*((xi + 1)*(xi + 1)*(xi + 1)) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)));
                case 10:
                  return sign * RealGradient(0, -3.75*eta - 3.375*xi + 67.5*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 33.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 9.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 16.875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 6.375 + 1.875*((eta + 1)*(eta + 1)) + 8.4375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.375*((xi + 1)*(xi + 1)) - 0.9375*(xi + 1)*(xi + 1)*(xi + 1) - 2.34375*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1));
                case 11:
                  return sign * RealGradient(0, 9.0*eta + 20.25*xi - 162.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 81.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 22.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 33.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 24.75 - 3.75*(eta + 1)*(eta + 1) - 16.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 20.25*(xi + 1)*(xi + 1) + 5.625*((xi + 1)*(xi + 1)*(xi + 1)) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)));
                case 12:
                  return RealGradient(0, 18.0*xi - 144.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 81.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 22.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 18.0 - 16.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 20.25*(xi + 1)*(xi + 1) + 5.625*((xi + 1)*(xi + 1)*(xi + 1)) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)));
                case 13:
                  return RealGradient(0, -4.5*xi + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 54.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 22.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 7.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 4.5 + 11.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 13.5*((xi + 1)*(xi + 1)) - 5.625*(xi + 1)*(xi + 1)*(xi + 1) - 4.6875*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1));
                case 14:
                  return RealGradient(-18.0*eta + 144.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 81.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 22.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 18.0 + 20.25*((eta + 1)*(eta + 1)) + 16.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 5.625*(eta + 1)*(eta + 1)*(eta + 1) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                case 15:
                  return RealGradient(4.5*eta - 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 7.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 54.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 22.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 4.5 - 13.5*(eta + 1)*(eta + 1) - 11.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5.625*((eta + 1)*(eta + 1)*(eta + 1)) + 4.6875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)), 0);
                case 16:
                  return RealGradient(-6.0*eta + 96.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 54.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 15.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.0 + 6.75*((eta + 1)*(eta + 1)) + 16.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1.875*(eta + 1)*(eta + 1)*(eta + 1) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                case 17:
                  return RealGradient(1.5*eta - 24.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 7.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 36.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 15.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1.5 - 4.5*(eta + 1)*(eta + 1) - 11.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1.875*((eta + 1)*(eta + 1)*(eta + 1)) + 4.6875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)), 0);
                case 18:
                  return RealGradient(0, 6.0*xi - 96.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 54.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 15.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 6.0 - 16.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 6.75*(xi + 1)*(xi + 1) + 1.875*((xi + 1)*(xi + 1)*(xi + 1)) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)));
                case 19:
                  return RealGradient(0, -1.5*xi + 24.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 36.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 15.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 7.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 1.5 + 11.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 4.5*((xi + 1)*(xi + 1)) - 1.875*(xi + 1)*(xi + 1)*(xi + 1) - 4.6875*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1));
                case 20:
                  return RealGradient(2.0*eta + 2.0 - 2.25*(eta + 1)*(eta + 1) + 0.625*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 21:
                  return RealGradient(0, 2.0*xi + 2.0 - 2.25*(xi + 1)*(xi + 1) + 0.625*((xi + 1)*(xi + 1)*(xi + 1)));
                case 22:
                  return RealGradient(0, -0.5*xi - 0.5 + 1.5*((xi + 1)*(xi + 1)) - 0.625*(xi + 1)*(xi + 1)*(xi + 1));
                case 23:
                  return RealGradient(-0.5*eta - 0.5 + 1.5*((eta + 1)*(eta + 1)) - 0.625*(eta + 1)*(eta + 1)*(eta + 1), 0);
                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient(120.0*eta*xi - 54.0*eta - 45.0*eta*xi*xi - 36.0*xi - 90.0*xi*eta*eta + 9.0 + 90.0*(eta*eta) + 30.0*(xi*xi) - 45.0*eta*eta*eta, -60.0*eta*xi + 90.0*eta*(xi*xi) + 18.0*xi + 45.0*xi*(eta*eta) - 60.0*xi*xi + 45.0*(xi*xi*xi));
                case 1:
                  return sign * RealGradient(-37.5*eta*xi - 2.25*eta + 22.5*eta*(xi*xi) + 15.0*xi + 22.5*xi*(eta*eta) - 1.5 + 15.0*(eta*eta) - 15.0*xi*xi - 11.25*eta*eta*eta, -22.5*eta*xi*xi - 5.25*xi + 11.25*xi*(eta*eta) + 26.25*(xi*xi) - 22.5*xi*xi*xi);
                case 2:
                  return sign * RealGradient(30.0*eta*xi - 3.0*eta - 45.0*eta*xi*xi - 24.0*xi + 3.0 + 30.0*(xi*xi), 9.0*xi - 45.0*xi*xi + 45.0*(xi*xi*xi));
                case 3:
                  return sign * RealGradient(30.0*eta*xi - 3.0*eta - 45.0*eta*xi*xi, 9.0*xi - 45.0*xi*xi + 45.0*(xi*xi*xi));
                case 4:
                  return sign * RealGradient(22.5*eta*xi - 6.0*eta - 11.25*eta*xi*xi - 45.0*xi*eta*eta + 18.75*(eta*eta) - 11.25*eta*eta*eta, -22.5*eta*xi + 45.0*eta*(xi*xi) + 6.0*xi + 11.25*xi*(eta*eta) - 18.75*xi*xi + 11.25*(xi*xi*xi));
                case 5:
                  return sign * RealGradient(-9.0*eta + 45.0*(eta*eta) - 45.0*eta*eta*eta, -30.0*eta*xi + 3.0*xi + 45.0*xi*(eta*eta));
                case 6:
                  return sign * RealGradient(-9.0*eta + 45.0*(eta*eta) - 45.0*eta*eta*eta, -30.0*eta*xi + 24.0*eta + 3.0*xi + 45.0*xi*(eta*eta) - 3.0 - 30.0*eta*eta);
                case 7:
                  return sign * RealGradient(5.25*eta - 11.25*eta*xi*xi + 22.5*xi*(eta*eta) - 26.25*eta*eta + 22.5*(eta*eta*eta), 37.5*eta*xi - 15.0*eta - 22.5*eta*xi*xi + 2.25*xi - 22.5*xi*eta*eta + 1.5 + 15.0*(eta*eta) - 15.0*xi*xi + 11.25*(xi*xi*xi));
                case 8:
                  return sign * RealGradient(60.0*eta*xi - 18.0*eta - 45.0*eta*xi*xi - 90.0*xi*eta*eta + 60.0*(eta*eta) - 45.0*eta*eta*eta, -120.0*eta*xi + 36.0*eta + 90.0*eta*(xi*xi) + 54.0*xi + 45.0*xi*(eta*eta) - 9.0 - 30.0*eta*eta - 90.0*xi*xi + 45.0*(xi*xi*xi));
                case 9:
                  return RealGradient(-300.0*eta*xi + 180.0*eta + 90.0*eta*(xi*xi) + 360.0*xi*(eta*eta) - 450.0*eta*eta + 270.0*(eta*eta*eta), 300.0*eta*xi - 360.0*eta*xi*xi - 60.0*xi - 270.0*xi*eta*eta + 150.0*(xi*xi) - 90.0*xi*xi*xi);
                case 10:
                  return RealGradient(300.0*eta*xi - 60.0*eta - 270.0*eta*xi*xi - 360.0*xi*eta*eta + 150.0*(eta*eta) - 90.0*eta*eta*eta, -300.0*eta*xi + 360.0*eta*(xi*xi) + 180.0*xi + 90.0*xi*(eta*eta) - 450.0*xi*xi + 270.0*(xi*xi*xi));
                case 11:
                  return RealGradient(360.0*eta*xi - 60.0*eta - 180.0*eta*xi*xi - 360.0*xi*eta*eta + 60.0*(eta*eta), -120.0*eta*xi + 360.0*eta*(xi*xi) + 60.0*xi - 240.0*xi*xi + 180.0*(xi*xi*xi));
                case 12:
                  return RealGradient(-240.0*eta*xi + 30.0*eta + 270.0*eta*(xi*xi) + 180.0*xi*(eta*eta) - 30.0*eta*eta, 60.0*eta*xi - 180.0*eta*xi*xi - 90.0*xi + 360.0*(xi*xi) - 270.0*xi*xi*xi);
                case 13:
                  return RealGradient(60.0*eta*xi - 90.0*eta - 180.0*xi*eta*eta + 360.0*(eta*eta) - 270.0*eta*eta*eta, -240.0*eta*xi + 180.0*eta*(xi*xi) + 30.0*xi + 270.0*xi*(eta*eta) - 30.0*xi*xi);
                case 14:
                  return RealGradient(-120.0*eta*xi + 60.0*eta + 360.0*xi*(eta*eta) - 240.0*eta*eta + 180.0*(eta*eta*eta), 360.0*eta*xi - 360.0*eta*xi*xi - 60.0*xi - 180.0*xi*eta*eta + 60.0*(xi*xi));
                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          } // end switch (type)
      } // end case THIRD

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}




template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const FEType fet,
                                      const Elem * elem,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape(elem, fet.order, i, p, add_p_level);
}

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const Elem * elem,
                                            const Order order,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);
  libmesh_assert_less (j, 2);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  libmesh_assert_less(i, n_dofs(elem->type(), total_order));

  const char sign = i >= total_order * elem->n_edges() || elem->point(i / total_order) > elem->point((i / total_order + 1) % elem->n_vertices()) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / total_order * 2 + 1) * total_order - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape function first derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                      case 2:
                        return RealGradient();
                      case 1:
                      case 3:
                        return sign * RealGradient( 0.0, -0.25 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 1:
                      case 3:
                        return RealGradient();
                      case 0:
                      case 2:
                        return sign * RealGradient( 0.25 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  return sign * RealGradient( 0.0, -1.0 );

                  // d()/deta
                case 1:
                  return sign * RealGradient( 1.0 );

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // quadratic Nedelec (first kind) shape function first derivatives
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              const Real x = 0.5 * (xi + 1.0);
              const Real y = 0.5 * (eta + 1.0);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.25*(-18.0*y*y+24.0*y-6.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.25*( 18.0*y*y-24.0*y+6.0), 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 0.25*(-36.0*x*y+24.0*x+12.0*y-8.0) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.25*( 36.0*x*y-12.0*x-12.0*y+4.0) );
                      case 4:
                        return sign * RealGradient( 0.25*(-18.0*y*y+12.0*y), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.25*( 18.0*y*y-12.0*y), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.25*(-36.0*x*y+12.0*x+24.0*y-8.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.25*( 36.0*x*y-24.0*x-24.0*y+16.0) );
                      case 8:
                        return RealGradient( 0.0, 1.5*(6.0*x*y-4.0*x-3.0*y+2.0) );
                      case 9:
                        return RealGradient( 1.5*y*(-3.0*y+3.0), 0.0 );
                      case 10:
                        return RealGradient( 1.5*y*(3.0*y-3.0), 0.0 );
                      case 11:
                        return RealGradient( 0.0, 1.5*(-6.0*x*y+2.0*x+3.0*y-1.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.25*(-36.0*x*y+24.0*x+24.0*y-16.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.25*( 36.0*x*y-24.0*x-12.0*y+8.0), 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 0.25*x*(-18.0*x+12.0) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.25*x*( 18.0*x-12.0) );
                      case 4:
                        return sign * RealGradient( 0.25*(-36.0*x*y+12.0*x+12.0*y-4.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.25*( 36.0*x*y-12.0*x-24.0*y+8.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.25*(-18.0*x*x+24.0*x-6.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.25*( 18.0*x*x-24.0*x+6.0) );
                      case 8:
                        return RealGradient( 0.0, 1.5*x*(3.0*x-3.0) );
                      case 9:
                        return RealGradient( 1.5*(-6.0*x*y+3.0*x+4.0*y-2.0), 0.0 );
                      case 10:
                        return RealGradient( 1.5*(6.0*x*y-3.0*x-2.0*y+1.0), 0.0 );
                      case 11:
                        return RealGradient( 0.0, 1.5*x*(-3.0*x+3.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 8.0*eta-6.0, -16.0*xi-8.0*eta+6.0 );
                      case 1:
                        return sign * RealGradient( -8.0*eta+6.0, 16.0*xi-4.0 );
                      case 2:
                        return sign * RealGradient( -8.0*eta, 16.0*xi-4.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 8.0*eta-2.0 );
                      case 4:
                        return sign * RealGradient( 0.0, 8.0*eta-2.0 );
                      case 5:
                        return sign * RealGradient( 8.0*eta, -16.0*xi-8.0*eta+12.0 );
                      case 6:
                        return RealGradient( -8.0*eta, 16.0*xi+16.0*eta-8.0 );
                      case 7:
                        return RealGradient( 16.0*eta, -32.0*xi-8.0*eta+16.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 8.0*xi+16.0*eta-12.0, -8.0*xi );
                      case 1:
                        return sign * RealGradient( -8.0*xi+2.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 2.0-8.0*xi, 0.0 );
                      case 3:
                        return sign * RealGradient( 4.0-16.0*eta, 8.0*xi );
                      case 4:
                        return sign * RealGradient( 4.0-16.0*eta,8.0*xi-6.0 );
                      case 5:
                        return sign * RealGradient( 8.0*xi+16*eta-6.0,-8.0*xi+6.0 );
                      case 6:
                        return RealGradient( -8.0*xi-32.0*eta+16.0, 16.0*xi );
                      case 7:
                        return RealGradient( 16.0*xi+16.0*eta-8.0,-8.0*xi );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // cubic Nedelec (first kind) shape function first derivatives
    case THIRD:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(40.5*eta + 7.5*xi - 33.75*(0.5*eta + 0.5)*(2*xi + 2) + 16.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 39.0 - 40.5*(eta + 1)*(eta + 1) + 11.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 1:
                        return sign * RealGradient(-16.875*eta - 3.75*xi + 16.875*(0.5*eta + 0.5)*(2*xi + 2) - 8.4375*(2*xi + 2)*(eta + 1)*(eta + 1) + 2.34375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 16.875 + 16.875*((eta + 1)*(eta + 1)) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 2:
                        return sign * RealGradient(27.0*eta + 7.5*xi - 33.75*(0.5*eta + 0.5)*(2*xi + 2) + 16.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 28.5 - 27.0*(eta + 1)*(eta + 1) + 7.5*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 3:
                        return sign * RealGradient(0, -13.5*eta - 27.0*xi + 54.0*(0.5*eta + 0.5)*(2*xi + 2) - 67.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 11.25*(2*xi + 2)*(eta + 1)*(eta + 1) - 33.75 + 5.625*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 16.875*((xi + 1)*(xi + 1)));
                      case 4:
                        return sign * RealGradient(0, 5.625*eta + 4.5*xi - 22.5*(0.5*eta + 0.5)*(2*xi + 2) + 28.125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 5.625*(2*xi + 2)*((eta + 1)*(eta + 1)) + 9.0 - 2.8125*(eta + 1)*(eta + 1) - 7.03125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 2.8125*(xi + 1)*(xi + 1));
                      case 5:
                        return sign * RealGradient(0, -9.0*eta - 9.0*xi + 36.0*(0.5*eta + 0.5)*(2*xi + 2) - 45.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 11.25*(2*xi + 2)*(eta + 1)*(eta + 1) - 15.75 + 5.625*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 5.625*((xi + 1)*(xi + 1)));
                      case 6:
                        return sign * RealGradient(9.0*eta - 11.25*(0.5*eta + 0.5)*(2*xi + 2) + 11.25*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 9.0 - 18.0*(eta + 1)*(eta + 1) + 7.5*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 7:
                        return sign * RealGradient(-5.625*eta + 5.625*(0.5*eta + 0.5)*(2*xi + 2) - 5.625*(2*xi + 2)*(eta + 1)*(eta + 1) + 2.34375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 5.625 + 11.25*((eta + 1)*(eta + 1)) - 4.6875*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 8:
                        return sign * RealGradient(13.5*eta - 11.25*(0.5*eta + 0.5)*(2*xi + 2) + 11.25*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 13.5 - 27.0*(eta + 1)*(eta + 1) + 11.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 9:
                        return sign * RealGradient(0, -27.0*eta - 13.5*xi + 54.0*(0.5*eta + 0.5)*(2*xi + 2) - 45.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 16.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 33.75 + 16.875*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 5.625*((xi + 1)*(xi + 1)));
                      case 10:
                        return sign * RealGradient(0, 16.875*eta + 6.75*xi - 33.75*(0.5*eta + 0.5)*(2*xi + 2) + 28.125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 8.4375*(2*xi + 2)*((eta + 1)*(eta + 1)) + 20.25 - 8.4375*(eta + 1)*(eta + 1) - 7.03125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 2.8125*(xi + 1)*(xi + 1));
                      case 11:
                        return sign * RealGradient(0, -40.5*eta - 40.5*xi + 81.0*(0.5*eta + 0.5)*(2*xi + 2) - 67.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 16.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 60.75 + 16.875*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 16.875*((xi + 1)*(xi + 1)));
                      case 12:
                        return RealGradient(0, -36.0*eta - 40.5*xi + 81.0*(0.5*eta + 0.5)*(2*xi + 2) - 67.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 16.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 58.5 + 15.0*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 16.875*((xi + 1)*(xi + 1)));
                      case 13:
                        return RealGradient(0, 9.0*eta + 27.0*xi - 54.0*(0.5*eta + 0.5)*(2*xi + 2) + 67.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 11.25*(2*xi + 2)*((eta + 1)*(eta + 1)) + 31.5 - 3.75*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 16.875*(xi + 1)*(xi + 1));
                      case 14:
                        return RealGradient(36.0*eta - 30.0*(0.5*eta + 0.5)*(2*xi + 2) + 16.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 36.0 - 40.5*(eta + 1)*(eta + 1) + 11.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 15:
                        return RealGradient(-9.0*eta + 7.5*(0.5*eta + 0.5)*(2*xi + 2) - 11.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 4.6875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 9.0 + 27.0*((eta + 1)*(eta + 1)) - 11.25*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 16:
                        return RealGradient(24.0*eta - 30.0*(0.5*eta + 0.5)*(2*xi + 2) + 16.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 24.0 - 27.0*(eta + 1)*(eta + 1) + 7.5*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 17:
                        return RealGradient(-6.0*eta + 7.5*(0.5*eta + 0.5)*(2*xi + 2) - 11.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 4.6875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.0 + 18.0*((eta + 1)*(eta + 1)) - 7.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 18:
                        return RealGradient(0, -24.0*eta - 13.5*xi + 54.0*(0.5*eta + 0.5)*(2*xi + 2) - 45.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 16.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 31.5 + 15.0*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 5.625*((xi + 1)*(xi + 1)));
                      case 19:
                        return RealGradient(0, 6.0*eta + 9.0*xi - 36.0*(0.5*eta + 0.5)*(2*xi + 2) + 45.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 11.25*(2*xi + 2)*((eta + 1)*(eta + 1)) + 13.5 - 3.75*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 5.625*(xi + 1)*(xi + 1));
                      case 20:
                        return RealGradient(0, 0);
                      case 21:
                        return RealGradient(0, -4.5*xi - 2.5 + 1.875*((xi + 1)*(xi + 1)));
                      case 22:
                        return RealGradient(0, 3.0*xi + 2.5 - 1.875*(xi + 1)*(xi + 1));
                      case 23:
                        return RealGradient(0, 0);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(40.5*eta + 40.5*xi - 81.0*(2*eta + 2)*(0.5*xi + 0.5) + 16.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 67.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 60.75 - 16.875*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 16.875*(xi + 1)*(xi + 1), 0);
                      case 1:
                        return sign * RealGradient(-6.75*eta - 16.875*xi + 33.75*(2*eta + 2)*(0.5*xi + 0.5) - 8.4375*(2*eta + 2)*(xi + 1)*(xi + 1) - 28.125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 20.25 + 2.8125*((eta + 1)*(eta + 1)) + 7.03125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 8.4375*((xi + 1)*(xi + 1)), 0);
                      case 2:
                        return sign * RealGradient(13.5*eta + 27.0*xi - 54.0*(2*eta + 2)*(0.5*xi + 0.5) + 16.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 45.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 33.75 - 5.625*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 16.875*(xi + 1)*(xi + 1), 0);
                      case 3:
                        return sign * RealGradient(0, -13.5*xi + 11.25*(2*eta + 2)*(0.5*xi + 0.5) - 11.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 13.5 + 27.0*((xi + 1)*(xi + 1)) - 11.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 4:
                        return sign * RealGradient(0, 5.625*xi - 5.625*(2*eta + 2)*(0.5*xi + 0.5) + 5.625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 2.34375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 5.625 - 11.25*(xi + 1)*(xi + 1) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 5:
                        return sign * RealGradient(0, -9.0*xi + 11.25*(2*eta + 2)*(0.5*xi + 0.5) - 11.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 9.0 + 18.0*((xi + 1)*(xi + 1)) - 7.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 6:
                        return sign * RealGradient(9.0*eta + 9.0*xi - 36.0*(2*eta + 2)*(0.5*xi + 0.5) + 11.25*(2*eta + 2)*((xi + 1)*(xi + 1)) + 45.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 15.75 - 5.625*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 5.625*(xi + 1)*(xi + 1), 0);
                      case 7:
                        return sign * RealGradient(-4.5*eta - 5.625*xi + 22.5*(2*eta + 2)*(0.5*xi + 0.5) - 5.625*(2*eta + 2)*(xi + 1)*(xi + 1) - 28.125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 9.0 + 2.8125*((eta + 1)*(eta + 1)) + 7.03125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 2.8125*((xi + 1)*(xi + 1)), 0);
                      case 8:
                        return sign * RealGradient(27.0*eta + 13.5*xi - 54.0*(2*eta + 2)*(0.5*xi + 0.5) + 11.25*(2*eta + 2)*((xi + 1)*(xi + 1)) + 67.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 33.75 - 16.875*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 5.625*(xi + 1)*(xi + 1), 0);
                      case 9:
                        return sign * RealGradient(0, -7.5*eta - 27.0*xi + 33.75*(2*eta + 2)*(0.5*xi + 0.5) - 16.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 28.5 + 27.0*((xi + 1)*(xi + 1)) - 7.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 10:
                        return sign * RealGradient(0, 3.75*eta + 16.875*xi - 16.875*(2*eta + 2)*(0.5*xi + 0.5) + 8.4375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 2.34375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 16.875 - 16.875*(xi + 1)*(xi + 1) + 4.6875*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 11:
                        return sign * RealGradient(0, -7.5*eta - 40.5*xi + 33.75*(2*eta + 2)*(0.5*xi + 0.5) - 16.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 39.0 + 40.5*((xi + 1)*(xi + 1)) - 11.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 12:
                        return RealGradient(0, -36.0*xi + 30.0*(2*eta + 2)*(0.5*xi + 0.5) - 16.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 36.0 + 40.5*((xi + 1)*(xi + 1)) - 11.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 13:
                        return RealGradient(0, 9.0*xi - 7.5*(2*eta + 2)*(0.5*xi + 0.5) + 11.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 4.6875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 9.0 - 27.0*(xi + 1)*(xi + 1) + 11.25*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 14:
                        return RealGradient(40.5*eta + 36.0*xi - 81.0*(2*eta + 2)*(0.5*xi + 0.5) + 16.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 67.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 58.5 - 16.875*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1), 0);
                      case 15:
                        return RealGradient(-27.0*eta - 9.0*xi + 54.0*(2*eta + 2)*(0.5*xi + 0.5) - 11.25*(2*eta + 2)*(xi + 1)*(xi + 1) - 67.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 31.5 + 16.875*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.75*((xi + 1)*(xi + 1)), 0);
                      case 16:
                        return RealGradient(13.5*eta + 24.0*xi - 54.0*(2*eta + 2)*(0.5*xi + 0.5) + 16.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 45.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 31.5 - 5.625*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1), 0);
                      case 17:
                        return RealGradient(-9.0*eta - 6.0*xi + 36.0*(2*eta + 2)*(0.5*xi + 0.5) - 11.25*(2*eta + 2)*(xi + 1)*(xi + 1) - 45.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 13.5 + 5.625*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.75*((xi + 1)*(xi + 1)), 0);
                      case 18:
                        return RealGradient(0, -24.0*xi + 30.0*(2*eta + 2)*(0.5*xi + 0.5) - 16.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 24.0 + 27.0*((xi + 1)*(xi + 1)) - 7.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 19:
                        return RealGradient(0, 6.0*xi - 7.5*(2*eta + 2)*(0.5*xi + 0.5) + 11.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 4.6875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 6.0 - 18.0*(xi + 1)*(xi + 1) + 7.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 20:
                        return RealGradient(-4.5*eta - 2.5 + 1.875*((eta + 1)*(eta + 1)), 0);
                      case 21:
                        return RealGradient(0, 0);
                      case 22:
                        return RealGradient(0, 0);
                      case 23:
                        return RealGradient(3.0*eta + 2.5 - 1.875*(eta + 1)*(eta + 1), 0);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-90.0*eta*xi + 120.0*eta + 60.0*xi - 36.0 - 90.0*eta*eta, 180.0*eta*xi - 60.0*eta - 120.0*xi + 18.0 + 45.0*(eta*eta) + 135.0*(xi*xi));
                      case 1:
                        return sign * RealGradient(45.0*eta*xi - 37.5*eta - 30.0*xi + 15.0 + 22.5*(eta*eta), -45.0*eta*xi + 52.5*xi - 5.25 + 11.25*(eta*eta) - 67.5*xi*xi);
                      case 2:
                        return sign * RealGradient(-90.0*eta*xi + 30.0*eta + 60.0*xi - 24.0, -90.0*xi + 9.0 + 135.0*(xi*xi));
                      case 3:
                        return sign * RealGradient(-90.0*eta*xi + 30.0*eta, -90.0*xi + 9.0 + 135.0*(xi*xi));
                      case 4:
                        return sign * RealGradient(-22.5*eta*xi + 22.5*eta - 45.0*eta*eta, 90.0*eta*xi - 22.5*eta - 37.5*xi + 6.0 + 11.25*(eta*eta) + 33.75*(xi*xi));
                      case 5:
                        return sign * RealGradient(0, -30.0*eta + 3.0 + 45.0*(eta*eta));
                      case 6:
                        return sign * RealGradient(0, -30.0*eta + 3.0 + 45.0*(eta*eta));
                      case 7:
                        return sign * RealGradient(-22.5*eta*xi + 22.5*(eta*eta), -45.0*eta*xi + 37.5*eta - 30.0*xi + 2.25 - 22.5*eta*eta + 33.75*(xi*xi));
                      case 8:
                        return sign * RealGradient(-90.0*eta*xi + 60.0*eta - 90.0*eta*eta, 180.0*eta*xi - 120.0*eta - 180.0*xi + 54.0 + 45.0*(eta*eta) + 135.0*(xi*xi));
                      case 9:
                        return RealGradient(180.0*eta*xi - 300.0*eta + 360.0*(eta*eta), -720.0*eta*xi + 300.0*eta + 300.0*xi - 60.0 - 270.0*eta*eta - 270.0*xi*xi);
                      case 10:
                        return RealGradient(-540.0*eta*xi + 300.0*eta - 360.0*eta*eta, 720.0*eta*xi - 300.0*eta - 900.0*xi + 180.0 + 90.0*(eta*eta) + 810.0*(xi*xi));
                      case 11:
                        return RealGradient(-360.0*eta*xi + 360.0*eta - 360.0*eta*eta, 720.0*eta*xi - 120.0*eta - 480.0*xi + 60.0 + 540.0*(xi*xi));
                      case 12:
                        return RealGradient(540.0*eta*xi - 240.0*eta + 180.0*(eta*eta), -360.0*eta*xi + 60.0*eta + 720.0*xi - 90.0 - 810.0*xi*xi);
                      case 13:
                        return RealGradient(60.0*eta - 180.0*eta*eta, 360.0*eta*xi - 240.0*eta - 60.0*xi + 30.0 + 270.0*(eta*eta));
                      case 14:
                        return RealGradient(-120.0*eta + 360.0*(eta*eta), -720.0*eta*xi + 360.0*eta + 120.0*xi - 60.0 - 180.0*eta*eta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-180.0*eta*xi + 180.0*eta + 120.0*xi - 54.0 - 135.0*eta*eta - 45.0*xi*xi, 90.0*eta*xi - 60.0*xi + 90.0*(xi*xi));
                      case 1:
                        return sign * RealGradient(45.0*eta*xi + 30.0*eta - 37.5*xi - 2.25 - 33.75*eta*eta + 22.5*(xi*xi), 22.5*eta*xi - 22.5*xi*xi);
                      case 2:
                        return sign * RealGradient(30.0*xi - 3.0 - 45.0*xi*xi, 0);
                      case 3:
                        return sign * RealGradient(30.0*xi - 3.0 - 45.0*xi*xi, 0);
                      case 4:
                        return sign * RealGradient(-90.0*eta*xi + 37.5*eta + 22.5*xi - 6.0 - 33.75*eta*eta - 11.25*xi*xi, 22.5*eta*xi - 22.5*xi + 45.0*(xi*xi));
                      case 5:
                        return sign * RealGradient(90.0*eta - 9.0 - 135.0*eta*eta, 90.0*eta*xi - 30.0*xi);
                      case 6:
                        return sign * RealGradient(90.0*eta - 9.0 - 135.0*eta*eta, 90.0*eta*xi - 60.0*eta - 30.0*xi + 24.0);
                      case 7:
                        return sign * RealGradient(45.0*eta*xi - 52.5*eta + 5.25 + 67.5*(eta*eta) - 11.25*xi*xi, -45.0*eta*xi + 30.0*eta + 37.5*xi - 15.0 - 22.5*xi*xi);
                      case 8:
                        return sign * RealGradient(-180.0*eta*xi + 120.0*eta + 60.0*xi - 18.0 - 135.0*eta*eta - 45.0*xi*xi, 90.0*eta*xi - 60.0*eta - 120.0*xi + 36.0 + 90.0*(xi*xi));
                      case 9:
                        return RealGradient(720.0*eta*xi - 900.0*eta - 300.0*xi + 180.0 + 810.0*(eta*eta) + 90.0*(xi*xi), -540.0*eta*xi + 300.0*xi - 360.0*xi*xi);
                      case 10:
                        return RealGradient(-720.0*eta*xi + 300.0*eta + 300.0*xi - 60.0 - 270.0*eta*eta - 270.0*xi*xi, 180.0*eta*xi - 300.0*xi + 360.0*(xi*xi));
                      case 11:
                        return RealGradient(-720.0*eta*xi + 120.0*eta + 360.0*xi - 60.0 - 180.0*xi*xi, -120.0*xi + 360.0*(xi*xi));
                      case 12:
                        return RealGradient(360.0*eta*xi - 60.0*eta - 240.0*xi + 30.0 + 270.0*(xi*xi), 60.0*xi - 180.0*xi*xi);
                      case 13:
                        return RealGradient(-360.0*eta*xi + 720.0*eta + 60.0*xi - 90.0 - 810.0*eta*eta, 540.0*eta*xi - 240.0*xi + 180.0*(xi*xi));
                      case 14:
                        return RealGradient(720.0*eta*xi - 480.0*eta - 120.0*xi + 60.0 + 540.0*(eta*eta), -360.0*eta*xi + 360.0*xi - 360.0*xi*xi);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          } // end switch (type)
      } // end case THIRD

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(elem, order, i, j, add_p_level);
  libmesh_not_implemented();
#endif
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const ElemType,
                                            const Order,
                                            const unsigned int,
                                            const unsigned int,
                                            const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const FEType fet,
                                            const Elem * elem,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}





#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const Elem * elem,
                                                   const Order order,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point & p,
                                                   const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  libmesh_assert_less (j, 3);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  libmesh_assert_less(i, n_dofs(elem->type(), total_order));

  const char sign = i >= total_order * elem->n_edges() || elem->point(i / total_order) > elem->point((i / total_order + 1) % elem->n_vertices()) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / total_order * 2 + 1) * total_order - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape function second derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
          case TRI6:
          case TRI7:
            // All second derivatives for linear quads and triangles are zero.
            return RealGradient();

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));

          } // end switch (type)
      } // end case FIRST

      // quadratic Nedelec (first kind) shape function second derivatives
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              const Real x = 0.5 * (xi + 1.0);
              const Real y = 0.5 * (eta + 1.0);

              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                      case 1:
                      case 4:
                      case 5:
                      case 9:
                      case 10:
                        return RealGradient();
                      case 2:
                        return sign * RealGradient( 0.0, 0.125*(-36.0*y+24.0) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.125*( 36.0*y-12.0) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.125*(-36.0*y+12.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.125*( 36.0*y-24.0) );
                      case 8:
                        return RealGradient( 0.0,  0.75*(6.0*y-4.0) );
                      case 11:
                        return RealGradient( 0.0,  0.75*(-6.0*y+2.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d^2()/dxideta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.125*(-36.0*y+24.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.125*( 36.0*y-24.0), 0.0 );
                      case 2:
                       return sign * RealGradient( 0.0, 0.125*(-36.0*x+12.0) );
                      case 3:
                       return sign * RealGradient( 0.0, 0.125*( 36.0*x-12.0) );
                      case 4:
                        return sign * RealGradient( 0.125*(-36.0*y+12.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.125*( 36.0*y-12.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.125*(-36.0*x+24.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.125*( 36.0*x-24.0) );
                      case 8:
                        return RealGradient( 0.0,  0.75*(6.0*x-3.0) );
                      case 9:
                        return RealGradient( 0.75*(-6.0*y), 0.0 );
                      case 10:
                        return RealGradient( 0.75*(6.0*y), 0.0 );
                      case 11:
                        return RealGradient( 0.0, 0.75*(-6.0*x+3.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                  // d^2()/deta^2
                case 2:
                  {
                    switch(ii)
                      {
                      case 2:
                      case 3:
                      case 6:
                      case 7:
                      case 8:
                      case 11:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( 0.125*(-36.0*x+24.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.125*( 36.0*x-12.0), 0.0 );
                      case 4:
                        return sign * RealGradient( 0.125*(-36.0*x+12.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.125*( 36.0*x-24.0), 0.0 );
                      case 9:
                        return RealGradient( 0.75*(-6.0*x+4.0), 0.0 );
                      case 10:
                        return RealGradient( 0.75*( 6.0*x-2.0), 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    switch(ii)
                      {
                      case 3:
                      case 4:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( 0.0, -16.0 );
                      case 1:
                        return sign * RealGradient( 0.0, 16.0 );
                      case 2:
                        return sign * RealGradient( 0.0, 16.0 );
                      case 5:
                        return sign * RealGradient( 0.0, -16.0 );
                      case 6:
                        return RealGradient( 0.0, 16.0 );
                      case 7:
                        return RealGradient( 0.0,-32.0 );
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d^2()/dxideta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 8.0, -8.0 );
                      case 1:
                        return sign * RealGradient( -8.0, 0.0 );
                      case 2:
                        return sign * RealGradient( -8.0, 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 8.0 );
                      case 4:
                        return sign * RealGradient( 0.0, 8.0 );
                      case 5:
                        return sign * RealGradient( 8.0, -8.0 );
                      case 6:
                        return RealGradient( -8.0, 16.0 );
                      case 7:
                        return RealGradient( 16.0, -8.0 );
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                  // d^2()/deta^2
                case 2:
                  {
                    switch(ii)
                      {
                      case 1:
                      case 2:
                        return RealGradient();
                      case 0:
                        return sign * RealGradient( 16.0, 0.0 );
                      case 3:
                        return sign * RealGradient( -16.0, 0.0 );
                      case 4:
                        return sign * RealGradient( -16.0, 0.0 );
                      case 5:
                        return sign * RealGradient( 16.0, 0.0 );
                      case 6:
                        return RealGradient( -32.0, 0.0 );
                      case 7:
                        return RealGradient( 16.0, 0.0 );
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));

          } // end switch (type)
      } // end case SECOND

      // cubic Nedelec (first kind) shape function second derivatives
    case THIRD:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-33.75*eta - 26.25 + 33.75*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 1:
                        return sign * RealGradient(16.875*eta + 13.125 - 16.875*(eta + 1)*(eta + 1) + 4.6875*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 2:
                        return sign * RealGradient(-33.75*eta - 26.25 + 33.75*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 3:
                        return sign * RealGradient(0, 54.0*eta + 33.75*xi - 67.5*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 60.75 - 22.5*(eta + 1)*(eta + 1));
                      case 4:
                        return sign * RealGradient(0, -22.5*eta - 5.625*xi + 28.125*(eta + 1)*(xi + 1) - 14.0625*(xi + 1)*(eta + 1)*(eta + 1) - 23.625 + 11.25*((eta + 1)*(eta + 1)));
                      case 5:
                        return sign * RealGradient(0, 36.0*eta + 11.25*xi - 45.0*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 38.25 - 22.5*(eta + 1)*(eta + 1));
                      case 6:
                        return sign * RealGradient(-11.25*eta - 11.25 + 22.5*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 7:
                        return sign * RealGradient(5.625*eta + 5.625 - 11.25*(eta + 1)*(eta + 1) + 4.6875*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 8:
                        return sign * RealGradient(-11.25*eta - 11.25 + 22.5*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 9:
                        return sign * RealGradient(0, 54.0*eta + 11.25*xi - 45.0*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 51.75 - 33.75*(eta + 1)*(eta + 1));
                      case 10:
                        return sign * RealGradient(0, -33.75*eta - 5.625*xi + 28.125*(eta + 1)*(xi + 1) - 14.0625*(xi + 1)*(eta + 1)*(eta + 1) - 32.625 + 16.875*((eta + 1)*(eta + 1)));
                      case 11:
                        return sign * RealGradient(0, 81.0*eta + 33.75*xi - 67.5*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 74.25 - 33.75*(eta + 1)*(eta + 1));
                      case 12:
                        return RealGradient(0, 81.0*eta + 33.75*xi - 67.5*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 74.25 - 33.75*(eta + 1)*(eta + 1));
                      case 13:
                        return RealGradient(0, -54.0*eta - 33.75*xi + 67.5*(eta + 1)*(xi + 1) - 28.125*(xi + 1)*(eta + 1)*(eta + 1) - 60.75 + 22.5*((eta + 1)*(eta + 1)));
                      case 14:
                        return RealGradient(-30.0*eta - 30.0 + 33.75*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 15:
                        return RealGradient(7.5*eta + 7.5 - 22.5*(eta + 1)*(eta + 1) + 9.375*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 16:
                        return RealGradient(-30.0*eta - 30.0 + 33.75*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 17:
                        return RealGradient(7.5*eta + 7.5 - 22.5*(eta + 1)*(eta + 1) + 9.375*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 18:
                        return RealGradient(0, 54.0*eta + 11.25*xi - 45.0*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 51.75 - 33.75*(eta + 1)*(eta + 1));
                      case 19:
                        return RealGradient(0, -36.0*eta - 11.25*xi + 45.0*(eta + 1)*(xi + 1) - 28.125*(xi + 1)*(eta + 1)*(eta + 1) - 38.25 + 22.5*((eta + 1)*(eta + 1)));
                      case 20:
                        return RealGradient(0, 0);
                      case 21:
                        return RealGradient(0, 3.75*xi - 0.75);
                      case 22:
                        return RealGradient(0, -3.75*xi - 0.75);
                      case 23:
                        return RealGradient(0, 0);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d^2()/dxideta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-81.0*eta - 33.75*xi + 16.875*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 74.25 + 33.75*((eta + 1)*(eta + 1)), 0);
                      case 1:
                        return sign * RealGradient(33.75*eta + 16.875*xi - 8.4375*(2*eta + 2)*(2*xi + 2) + 7.03125*(2*xi + 2)*((eta + 1)*(eta + 1)) + 33.75 - 14.0625*(eta + 1)*(eta + 1), 0);
                      case 2:
                        return sign * RealGradient(-54.0*eta - 33.75*xi + 16.875*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 60.75 + 22.5*((eta + 1)*(eta + 1)), 0);
                      case 3:
                        return sign * RealGradient(0, 11.25*eta + 54.0*xi - 11.25*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 51.75 - 33.75*(xi + 1)*(xi + 1));
                      case 4:
                        return sign * RealGradient(0, -5.625*eta - 22.5*xi + 5.625*(2*eta + 2)*(2*xi + 2) - 7.03125*(2*eta + 2)*(xi + 1)*(xi + 1) - 22.5 + 14.0625*((xi + 1)*(xi + 1)));
                      case 5:
                        return sign * RealGradient(0, 11.25*eta + 36.0*xi - 11.25*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 38.25 - 22.5*(xi + 1)*(xi + 1));
                      case 6:
                        return sign * RealGradient(-36.0*eta - 11.25*xi + 11.25*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 38.25 + 22.5*((eta + 1)*(eta + 1)), 0);
                      case 7:
                        return sign * RealGradient(22.5*eta + 5.625*xi - 5.625*(2*eta + 2)*(2*xi + 2) + 7.03125*(2*xi + 2)*((eta + 1)*(eta + 1)) + 22.5 - 14.0625*(eta + 1)*(eta + 1), 0);
                      case 8:
                        return sign * RealGradient(-54.0*eta - 11.25*xi + 11.25*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 51.75 + 33.75*((eta + 1)*(eta + 1)), 0);
                      case 9:
                        return sign * RealGradient(0, 33.75*eta + 54.0*xi - 16.875*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 60.75 - 22.5*(xi + 1)*(xi + 1));
                      case 10:
                        return sign * RealGradient(0, -16.875*eta - 33.75*xi + 8.4375*(2*eta + 2)*(2*xi + 2) - 7.03125*(2*eta + 2)*(xi + 1)*(xi + 1) - 33.75 + 14.0625*((xi + 1)*(xi + 1)));
                      case 11:
                        return sign * RealGradient(0, 33.75*eta + 81.0*xi - 16.875*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 74.25 - 33.75*(xi + 1)*(xi + 1));
                      case 12:
                        return RealGradient(0, 30.0*eta + 81.0*xi - 16.875*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 75.0 - 33.75*(xi + 1)*(xi + 1));
                      case 13:
                        return RealGradient(0, -7.5*eta - 54.0*xi + 11.25*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*eta + 2)*(xi + 1)*(xi + 1) - 52.5 + 33.75*((xi + 1)*(xi + 1)));
                      case 14:
                        return RealGradient(-81.0*eta - 30.0*xi + 16.875*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 75.0 + 33.75*((eta + 1)*(eta + 1)), 0);
                      case 15:
                        return RealGradient(54.0*eta + 7.5*xi - 11.25*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*xi + 2)*((eta + 1)*(eta + 1)) + 52.5 - 33.75*(eta + 1)*(eta + 1), 0);
                      case 16:
                        return RealGradient(-54.0*eta - 30.0*xi + 16.875*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 60.0 + 22.5*((eta + 1)*(eta + 1)), 0);
                      case 17:
                        return RealGradient(36.0*eta + 7.5*xi - 11.25*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*xi + 2)*((eta + 1)*(eta + 1)) + 37.5 - 22.5*(eta + 1)*(eta + 1), 0);
                      case 18:
                        return RealGradient(0, 30.0*eta + 54.0*xi - 16.875*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 60.0 - 22.5*(xi + 1)*(xi + 1));
                      case 19:
                        return RealGradient(0, -7.5*eta - 36.0*xi + 11.25*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*eta + 2)*(xi + 1)*(xi + 1) - 37.5 + 22.5*((xi + 1)*(xi + 1)));
                      case 20:
                        return RealGradient(0, 0);
                      case 21:
                        return RealGradient(0, 0);
                      case 22:
                        return RealGradient(0, 0);
                      case 23:
                        return RealGradient(0, 0);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                  // d^2()/deta^2
                case 2:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-33.75*eta - 81.0*xi + 67.5*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 74.25 + 33.75*((xi + 1)*(xi + 1)), 0);
                      case 1:
                        return sign * RealGradient(5.625*eta + 33.75*xi - 28.125*(eta + 1)*(xi + 1) + 14.0625*(eta + 1)*((xi + 1)*(xi + 1)) + 32.625 - 16.875*(xi + 1)*(xi + 1), 0);
                      case 2:
                        return sign * RealGradient(-11.25*eta - 54.0*xi + 45.0*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 51.75 + 33.75*((xi + 1)*(xi + 1)), 0);
                      case 3:
                        return sign * RealGradient(0, 11.25*xi + 11.25 - 22.5*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 4:
                        return sign * RealGradient(0, -5.625*xi - 5.625 + 11.25*((xi + 1)*(xi + 1)) - 4.6875*(xi + 1)*(xi + 1)*(xi + 1));
                      case 5:
                        return sign * RealGradient(0, 11.25*xi + 11.25 - 22.5*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 6:
                        return sign * RealGradient(-11.25*eta - 36.0*xi + 45.0*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 38.25 + 22.5*((xi + 1)*(xi + 1)), 0);
                      case 7:
                        return sign * RealGradient(5.625*eta + 22.5*xi - 28.125*(eta + 1)*(xi + 1) + 14.0625*(eta + 1)*((xi + 1)*(xi + 1)) + 23.625 - 11.25*(xi + 1)*(xi + 1), 0);
                      case 8:
                        return sign * RealGradient(-33.75*eta - 54.0*xi + 67.5*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 60.75 + 22.5*((xi + 1)*(xi + 1)), 0);
                      case 9:
                        return sign * RealGradient(0, 33.75*xi + 26.25 - 33.75*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 10:
                        return sign * RealGradient(0, -16.875*xi - 13.125 + 16.875*((xi + 1)*(xi + 1)) - 4.6875*(xi + 1)*(xi + 1)*(xi + 1));
                      case 11:
                        return sign * RealGradient(0, 33.75*xi + 26.25 - 33.75*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 12:
                        return RealGradient(0, 30.0*xi + 30.0 - 33.75*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 13:
                        return RealGradient(0, -7.5*xi - 7.5 + 22.5*((xi + 1)*(xi + 1)) - 9.375*(xi + 1)*(xi + 1)*(xi + 1));
                      case 14:
                        return RealGradient(-33.75*eta - 81.0*xi + 67.5*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 74.25 + 33.75*((xi + 1)*(xi + 1)), 0);
                      case 15:
                        return RealGradient(33.75*eta + 54.0*xi - 67.5*(eta + 1)*(xi + 1) + 28.125*(eta + 1)*((xi + 1)*(xi + 1)) + 60.75 - 22.5*(xi + 1)*(xi + 1), 0);
                      case 16:
                        return RealGradient(-11.25*eta - 54.0*xi + 45.0*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 51.75 + 33.75*((xi + 1)*(xi + 1)), 0);
                      case 17:
                        return RealGradient(11.25*eta + 36.0*xi - 45.0*(eta + 1)*(xi + 1) + 28.125*(eta + 1)*((xi + 1)*(xi + 1)) + 38.25 - 22.5*(xi + 1)*(xi + 1), 0);
                      case 18:
                        return RealGradient(0, 30.0*xi + 30.0 - 33.75*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 19:
                        return RealGradient(0, -7.5*xi - 7.5 + 22.5*((xi + 1)*(xi + 1)) - 9.375*(xi + 1)*(xi + 1)*(xi + 1));
                      case 20:
                        return RealGradient(3.75*eta - 0.75, 0);
                      case 21:
                        return RealGradient(0, 0);
                      case 22:
                        return RealGradient(0, 0);
                      case 23:
                        return RealGradient(-3.75*eta - 0.75, 0);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(60.0 - 90.0*eta, 180.0*eta + 270.0*xi - 120.0);
                      case 1:
                        return sign * RealGradient(45.0*eta - 30.0, -45.0*eta - 135.0*xi + 52.5);
                      case 2:
                        return sign * RealGradient(60.0 - 90.0*eta, 270.0*xi - 90.0);
                      case 3:
                        return sign * RealGradient(-90.0*eta, 270.0*xi - 90.0);
                      case 4:
                        return sign * RealGradient(-22.5*eta, 90.0*eta + 67.5*xi - 37.5);
                      case 5:
                        return sign * RealGradient(0, 0);
                      case 6:
                        return sign * RealGradient(0, 0);
                      case 7:
                        return sign * RealGradient(-22.5*eta, -45.0*eta + 67.5*xi - 30.0);
                      case 8:
                        return sign * RealGradient(-90.0*eta, 180.0*eta + 270.0*xi - 180.0);
                      case 9:
                        return RealGradient(180.0*eta, -720.0*eta - 540.0*xi + 300.0);
                      case 10:
                        return RealGradient(-540.0*eta, 720.0*eta + 1620.0*xi - 900.0);
                      case 11:
                        return RealGradient(-360.0*eta, 720.0*eta + 1080.0*xi - 480.0);
                      case 12:
                        return RealGradient(540.0*eta, -360.0*eta - 1620.0*xi + 720.0);
                      case 13:
                        return RealGradient(0, 360.0*eta - 60.0);
                      case 14:
                        return RealGradient(0, 120.0 - 720.0*eta);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d^2()/dxideta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-180.0*eta - 90.0*xi + 120.0, 90.0*eta + 180.0*xi - 60.0);
                      case 1:
                        return sign * RealGradient(45.0*eta + 45.0*xi - 37.5, 22.5*eta - 45.0*xi);
                      case 2:
                        return sign * RealGradient(30.0 - 90.0*xi, 0);
                      case 3:
                        return sign * RealGradient(30.0 - 90.0*xi, 0);
                      case 4:
                        return sign * RealGradient(-90.0*eta - 22.5*xi + 22.5, 22.5*eta + 90.0*xi - 22.5);
                      case 5:
                        return sign * RealGradient(0, 90.0*eta - 30.0);
                      case 6:
                        return sign * RealGradient(0, 90.0*eta - 30.0);
                      case 7:
                        return sign * RealGradient(45.0*eta - 22.5*xi, -45.0*eta - 45.0*xi + 37.5);
                      case 8:
                        return sign * RealGradient(-180.0*eta - 90.0*xi + 60.0, 90.0*eta + 180.0*xi - 120.0);
                      case 9:
                        return RealGradient(720.0*eta + 180.0*xi - 300.0, -540.0*eta - 720.0*xi + 300.0);
                      case 10:
                        return RealGradient(-720.0*eta - 540.0*xi + 300.0, 180.0*eta + 720.0*xi - 300.0);
                      case 11:
                        return RealGradient(-720.0*eta - 360.0*xi + 360.0, 720.0*xi - 120.0);
                      case 12:
                        return RealGradient(360.0*eta + 540.0*xi - 240.0, 60.0 - 360.0*xi);
                      case 13:
                        return RealGradient(60.0 - 360.0*eta, 540.0*eta + 360.0*xi - 240.0);
                      case 14:
                        return RealGradient(720.0*eta - 120.0, -360.0*eta - 720.0*xi + 360.0);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                  // d^2()/deta^2
                case 2:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient(-270.0*eta - 180.0*xi + 180.0, 90.0*xi);
                      case 1:
                        return sign * RealGradient(-67.5*eta + 45.0*xi + 30.0, 22.5*xi);
                      case 2:
                        return sign * RealGradient(0, 0);
                      case 3:
                        return sign * RealGradient(0, 0);
                      case 4:
                        return sign * RealGradient(-67.5*eta - 90.0*xi + 37.5, 22.5*xi);
                      case 5:
                        return sign * RealGradient(90.0 - 270.0*eta, 90.0*xi);
                      case 6:
                        return sign * RealGradient(90.0 - 270.0*eta, 90.0*xi - 60.0);
                      case 7:
                        return sign * RealGradient(135.0*eta + 45.0*xi - 52.5, 30.0 - 45.0*xi);
                      case 8:
                        return sign * RealGradient(-270.0*eta - 180.0*xi + 120.0, 90.0*xi - 60.0);
                      case 9:
                        return RealGradient(1620.0*eta + 720.0*xi - 900.0, -540.0*xi);
                      case 10:
                        return RealGradient(-540.0*eta - 720.0*xi + 300.0, 180.0*xi);
                      case 11:
                        return RealGradient(120.0 - 720.0*xi, 0);
                      case 12:
                        return RealGradient(360.0*xi - 60.0, 0);
                      case 13:
                        return RealGradient(-1620.0*eta - 360.0*xi + 720.0, 540.0*xi);
                      case 14:
                        return RealGradient(1080.0*eta + 720.0*xi - 480.0, -360.0*xi);
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          } // end switch (type)
      } // end case THIRD

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);

    } // end switch (order)

#else // LIBMESH_DIM > 1
  libmesh_assert(true || i || j);
  libmesh_ignore(elem, order, add_p_level);
  libmesh_not_implemented();
#endif
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const ElemType,
                                                   const Order,
                                                   const unsigned int,
                                                   const unsigned int,
                                                   const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const FEType fet,
                                                   const Elem * elem,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point & p,
                                                   const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}



#endif

} // namespace libMesh
