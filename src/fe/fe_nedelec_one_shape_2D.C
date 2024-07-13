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

      // quartic Nedelec (first kind) shape functions
    case FOURTH:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient(-64.0*eta - 30.0*xi + 960.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 480.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 140.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 900.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 600.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 131.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 86.0 + 120.0*((eta + 1)*(eta + 1)) + 450.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 131.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0*((xi + 1)*(xi + 1)) - 300.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 80.0*(eta + 1)*(eta + 1)*(eta + 1) + 87.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 8.75*(xi + 1)*(xi + 1)*(xi + 1) - 19.140625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 1:
                  return sign * RealGradient(10.0740740740741*eta + 10.5555555555556*xi - 337.777777777778*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 191.111111111111*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 57.037037037037*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 316.666666666667*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 211.111111111111*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 46.1805555555556*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 19.3703703703704 - 18.8888888888889*(eta + 1)*(eta + 1) - 179.166666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 53.4722222222222*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 11.9444444444444*(xi + 1)*(xi + 1) + 119.444444444444*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 26.1284722222222*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 12.5925925925926*((eta + 1)*(eta + 1)*(eta + 1)) - 35.6481481481481*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3.56481481481481*((xi + 1)*(xi + 1)*(xi + 1)) + 7.79803240740741*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2.75462962962963*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 2:
                  return sign * RealGradient(-4.74074074074074*eta - 5.55555555555556*xi + 177.777777777778*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 151.111111111111*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 57.037037037037*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 166.666666666667*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 111.111111111111*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 24.3055555555556*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 9.7037037037037 + 8.88888888888889*((eta + 1)*(eta + 1)) + 141.666666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 53.4722222222222*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9.44444444444444*((xi + 1)*(xi + 1)) - 94.4444444444444*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 20.6597222222222*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 5.92592592592593*(eta + 1)*(eta + 1)*(eta + 1) + 35.6481481481481*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 3.56481481481481*(xi + 1)*(xi + 1)*(xi + 1) - 7.79803240740741*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1.2962962962963*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 3:
                  return sign * RealGradient(16.0*eta + 15.0*xi - 480.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 360.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 140.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 450.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 300.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 29.0 - 30.0*(eta + 1)*(eta + 1) - 337.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 131.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 22.5*(xi + 1)*(xi + 1) + 225.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 20.0*((eta + 1)*(eta + 1)*(eta + 1)) - 87.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 8.75*((xi + 1)*(xi + 1)*(xi + 1)) + 19.140625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 4.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 4:
                  return sign * RealGradient(0, -16.0*xi + 240.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 450.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 450.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 131.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 120.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 35.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 16.0 + 225.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 225.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 60.0*((xi + 1)*(xi + 1)) - 65.625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 60.0*(xi + 1)*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 5:
                  return sign * RealGradient(0, 2.51851851851852*xi - 84.4444444444444*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 158.333333333333*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 158.333333333333*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 46.1805555555556*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 47.7777777777778*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 14.2592592592593*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 2.51851851851852 - 89.5833333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 89.5833333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 26.1284722222222*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 9.44444444444444*(xi + 1)*(xi + 1) + 26.7361111111111*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 26.7361111111111*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.79803240740741*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 9.44444444444444*((xi + 1)*(xi + 1)*(xi + 1)) - 2.75462962962963*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 6:
                  return sign * RealGradient(0, -1.18518518518519*xi + 44.4444444444444*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 83.3333333333333*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 83.3333333333333*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 24.3055555555556*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 37.7777777777778*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 14.2592592592593*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 1.18518518518519 + 70.8333333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 70.8333333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 20.6597222222222*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 4.44444444444444*((xi + 1)*(xi + 1)) - 26.7361111111111*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 26.7361111111111*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 7.79803240740741*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 4.44444444444444*(xi + 1)*(xi + 1)*(xi + 1) + 1.2962962962963*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 7:
                  return sign * RealGradient(0, 4.0*xi - 120.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 225.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 225.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 90.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 35.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 4.0 - 168.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 168.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 49.21875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 15.0*((xi + 1)*(xi + 1)*(xi + 1)) - 4.375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 8:
                  return sign * RealGradient(-4.0*eta + 120.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 90.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 35.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 225.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 225.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 4.0 + 15.0*((eta + 1)*(eta + 1)) + 168.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 65.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 168.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 49.21875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 15.0*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4.375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 9:
                  return sign * RealGradient(1.18518518518519*eta - 44.4444444444444*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 37.7777777777778*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 14.2592592592593*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 83.3333333333333*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 83.3333333333333*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 24.3055555555556*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1.18518518518519 - 4.44444444444444*(eta + 1)*(eta + 1) - 70.8333333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 26.7361111111111*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 70.8333333333333*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 20.6597222222222*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4.44444444444444*((eta + 1)*(eta + 1)*(eta + 1)) - 26.7361111111111*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.79803240740741*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1.2962962962963*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 10:
                  return sign * RealGradient(-2.51851851851852*eta + 84.4444444444444*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 47.7777777777778*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 14.2592592592593*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 158.333333333333*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 158.333333333333*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 46.1805555555556*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 2.51851851851852 + 9.44444444444444*((eta + 1)*(eta + 1)) + 89.5833333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 26.7361111111111*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 89.5833333333333*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 26.1284722222222*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 9.44444444444444*(eta + 1)*(eta + 1)*(eta + 1) + 26.7361111111111*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 7.79803240740741*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2.75462962962963*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 11:
                  return sign * RealGradient(16.0*eta - 240.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 120.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 35.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 450.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 450.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 16.0 - 60.0*(eta + 1)*(eta + 1) - 225.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 65.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 225.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 60.0*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 17.5*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 12:
                  return sign * RealGradient(0, -15.0*eta - 16.0*xi + 480.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 450.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 300.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 360.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 140.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 29.0 + 22.5*((eta + 1)*(eta + 1)) + 337.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 225.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 30.0*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 8.75*(eta + 1)*(eta + 1)*(eta + 1) + 87.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 20.0*(xi + 1)*(xi + 1)*(xi + 1) + 4.375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 13:
                  return sign * RealGradient(0, 5.55555555555556*eta + 4.74074074074074*xi - 177.777777777778*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 166.666666666667*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 111.111111111111*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 24.3055555555556*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 151.111111111111*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 57.037037037037*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 9.7037037037037 - 9.44444444444444*(eta + 1)*(eta + 1) - 141.666666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 94.4444444444444*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 20.6597222222222*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 8.88888888888889*(xi + 1)*(xi + 1) + 53.4722222222222*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 3.56481481481481*((eta + 1)*(eta + 1)*(eta + 1)) - 35.6481481481481*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.79803240740741*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 5.92592592592593*((xi + 1)*(xi + 1)*(xi + 1)) - 1.2962962962963*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 14:
                  return sign * RealGradient(0, -10.5555555555556*eta - 10.0740740740741*xi + 337.777777777778*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 316.666666666667*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 211.111111111111*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 46.1805555555556*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 191.111111111111*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 57.037037037037*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 19.3703703703704 + 11.9444444444444*((eta + 1)*(eta + 1)) + 179.166666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 119.444444444444*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 26.1284722222222*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 18.8888888888889*((xi + 1)*(xi + 1)) - 53.4722222222222*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3.56481481481481*(eta + 1)*(eta + 1)*(eta + 1) + 35.6481481481481*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 7.79803240740741*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 12.5925925925926*(xi + 1)*(xi + 1)*(xi + 1) + 2.75462962962963*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 15:
                  return sign * RealGradient(0, 30.0*eta + 64.0*xi - 960.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 900.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 600.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 480.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 140.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 86.0 - 30.0*(eta + 1)*(eta + 1) - 450.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 300.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 120.0*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 8.75*((eta + 1)*(eta + 1)*(eta + 1)) - 87.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 80.0*((xi + 1)*(xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 16:
                  return RealGradient(0, 52.0*xi - 780.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 870.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 600.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 390.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 113.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 52.0 - 435.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 300.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 116.0*(xi + 1)*(xi + 1) + 126.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 87.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 80.0*((xi + 1)*(xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 17:
                  return RealGradient(0, 12.0*xi - 180.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 420.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 450.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 90.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 26.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 12.0 - 210.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 225.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 56.0*(xi + 1)*(xi + 1) + 61.25*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 60.0*((xi + 1)*(xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 18:
                  return RealGradient(0, 16.0*xi - 240.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 60.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 120.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 35.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 16.0 - 30.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 8.0*(xi + 1)*(xi + 1) + 8.75*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)));
                case 19:
                  return RealGradient(-52.0*eta + 780.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 390.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 113.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 870.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 600.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 131.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 52.0 + 116.0*((eta + 1)*(eta + 1)) + 435.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 126.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 300.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 80.0*(eta + 1)*(eta + 1)*(eta + 1) + 87.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 20:
                  return RealGradient(-12.0*eta + 180.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 90.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 26.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 420.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 450.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 131.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 12.0 + 56.0*((eta + 1)*(eta + 1)) + 210.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 61.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 225.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 60.0*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 21:
                  return RealGradient(-16.0*eta + 240.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 120.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 35.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 60.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 16.0 + 8.0*((eta + 1)*(eta + 1)) + 30.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 8.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                case 22:
                  return RealGradient(13.0*eta - 390.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 292.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 113.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 435.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 300.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 13.0 - 29.0*(eta + 1)*(eta + 1) - 326.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 126.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 225.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 20.0*((eta + 1)*(eta + 1)*(eta + 1)) - 87.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 4.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 23:
                  return RealGradient(3.0*eta - 90.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 67.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 26.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 210.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 225.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 3.0 - 14.0*(eta + 1)*(eta + 1) - 157.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 61.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 168.75*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 15.0*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 19.140625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 4.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 24:
                  return RealGradient(4.0*eta - 120.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 90.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 35.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 4.0 - 2.0*(eta + 1)*(eta + 1) - 22.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 8.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                case 25:
                  return RealGradient(0, -13.0*xi + 390.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 435.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 300.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 292.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 113.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 13.0 + 326.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 225.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 29.0*((xi + 1)*(xi + 1)) - 126.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 87.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 20.0*(xi + 1)*(xi + 1)*(xi + 1) + 4.375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 26:
                  return RealGradient(0, -3.0*xi + 90.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 210.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 225.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 67.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 26.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 3.0 + 157.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 168.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 14.0*((xi + 1)*(xi + 1)) - 61.25*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 19.140625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1)*(xi + 1) + 4.375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 27:
                  return RealGradient(0, -4.0*xi + 120.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 90.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 35.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 4.0 + 22.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 2.0*((xi + 1)*(xi + 1)) - 8.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                case 28:
                  return RealGradient(12.0*eta - 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 42.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 30.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 6.5625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 12.0 - 28.5*(eta + 1)*(eta + 1) + 20.0*((eta + 1)*(eta + 1)*(eta + 1)) - 4.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 29:
                  return RealGradient(-6.0*eta + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 42.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 30.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.5625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 6.0 + 14.25*((eta + 1)*(eta + 1)) - 10.0*(eta + 1)*(eta + 1)*(eta + 1) + 2.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 30:
                  return RealGradient(0, 12.0*xi - 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 42.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 6.5625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 12.0 - 28.5*(xi + 1)*(xi + 1) + 20.0*((xi + 1)*(xi + 1)*(xi + 1)) - 4.375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 31:
                  return RealGradient(0, -6.0*xi + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 42.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 30.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 6.5625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6.0 + 14.25*((xi + 1)*(xi + 1)) - 10.0*(xi + 1)*(xi + 1)*(xi + 1) + 2.1875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 32:
                  return RealGradient(0, 2.0*xi - 6.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 20.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 22.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 6.5625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2.0 - 13.5*(xi + 1)*(xi + 1) + 15.0*((xi + 1)*(xi + 1)*(xi + 1)) - 4.375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 33:
                  return RealGradient(0, -1.0*xi + 6.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 20.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 22.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 6.5625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1.0 + 6.75*((xi + 1)*(xi + 1)) - 7.5*(xi + 1)*(xi + 1)*(xi + 1) + 2.1875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 34:
                  return RealGradient(2.0*eta - 6.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 20.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 22.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 6.5625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2.0 - 13.5*(eta + 1)*(eta + 1) + 15.0*((eta + 1)*(eta + 1)*(eta + 1)) - 4.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 35:
                  return RealGradient(-1.0*eta + 6.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 20.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 22.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.5625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1.0 + 6.75*((eta + 1)*(eta + 1)) - 7.5*(eta + 1)*(eta + 1)*(eta + 1) + 2.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 36:
                  return RealGradient(0, 6.0*xi - 18.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 4.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 6.0 - 3.0*(xi + 1)*(xi + 1));
                case 37:
                  return RealGradient(-6.0*eta + 18.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 4.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 6.0 + 3.0*((eta + 1)*(eta + 1)), 0);
                case 38:
                  return RealGradient(3.0*eta - 18.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 4.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 3.0 - 1.5*(eta + 1)*(eta + 1), 0);
                case 39:
                  return RealGradient(0, -3.0*xi + 18.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 4.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 3.0 + 1.5*((xi + 1)*(xi + 1)));
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
                  return sign * RealGradient(720.0*eta*xi - 160.0*eta - 840.0*eta*xi*xi + 224.0*eta*(xi*xi*xi) - 120.0*xi - 1260.0*xi*eta*eta + 672.0*xi*(eta*eta*eta) + 16.0 + 480.0*(eta*eta) + 672.0*(eta*eta)*(xi*xi) + 240.0*(xi*xi) - 560.0*eta*eta*eta - 140.0*xi*xi*xi + 224.0*(eta*eta*eta*eta), -240.0*eta*xi + 840.0*eta*(xi*xi) - 672.0*eta*xi*xi*xi + 40.0*xi + 420.0*xi*(eta*eta) - 224.0*xi*eta*eta*eta - 672.0*eta*eta*xi*xi - 240.0*xi*xi + 420.0*(xi*xi*xi) - 224.0*xi*xi*xi*xi);
                case 1:
                  return sign * RealGradient(-200.0*eta*xi - 2.81481481481481*eta + 320.444444444444*eta*(xi*xi) - 91.2592592592593*eta*xi*xi*xi + 42.2222222222222*xi + 261.333333333333*xi*(eta*eta) - 99.5555555555556*xi*eta*eta*eta - 2.51851851851852 + 67.5555555555556*(eta*eta) - 248.888888888889*eta*eta*xi*xi - 95.5555555555556*xi*xi - 128.592592592593*eta*eta*eta + 57.037037037037*(xi*xi*xi) + 66.3703703703704*(eta*eta*eta*eta), 19.5555555555556*eta*xi - 236.444444444444*eta*xi*xi + 248.888888888889*eta*(xi*xi*xi) - 10.962962962963*xi + 49.7777777777778*xi*(eta*eta) - 66.3703703703704*xi*eta*eta*eta + 99.5555555555556*(eta*eta)*(xi*xi) + 88.8888888888889*(xi*xi) - 168.0*xi*xi*xi + 91.2592592592593*(xi*xi*xi*xi));
                case 2:
                  return sign * RealGradient(-16.0*eta*xi + 4.14814814814815*eta - 124.444444444444*eta*xi*xi + 91.2592592592593*eta*(xi*xi*xi) - 22.2222222222222*xi + 158.666666666667*xi*(eta*eta) - 124.444444444444*xi*eta*eta*eta + 1.18518518518519 - 11.5555555555556*eta*eta + 24.8888888888889*(eta*eta)*(xi*xi) + 75.5555555555556*(xi*xi) - 2.07407407407407*eta*eta*eta - 57.037037037037*xi*xi*xi + 8.2962962962963*(eta*eta*eta*eta), 12.4444444444444*eta*xi - 43.5555555555556*eta*xi*xi - 24.8888888888889*eta*xi*xi*xi + 5.62962962962963*xi - 21.7777777777778*xi*eta*eta - 8.2962962962963*xi*eta*eta*eta + 124.444444444444*(eta*eta)*(xi*xi) - 56.8888888888889*xi*xi + 140.0*(xi*xi*xi) - 91.2592592592593*xi*xi*xi*xi);
                case 3:
                  return sign * RealGradient(-72.0*eta*xi + 4.0*eta + 252.0*eta*(xi*xi) - 224.0*eta*xi*xi*xi + 60.0*xi - 4.0 - 180.0*xi*xi + 140.0*(xi*xi*xi), -16.0*xi + 144.0*(xi*xi) - 336.0*xi*xi*xi + 224.0*(xi*xi*xi*xi));
                case 4:
                  return sign * RealGradient(-72.0*eta*xi + 4.0*eta + 252.0*eta*(xi*xi) - 224.0*eta*xi*xi*xi, -16.0*xi + 144.0*(xi*xi) - 336.0*xi*xi*xi + 224.0*(xi*xi*xi*xi));
                case 5:
                  return sign * RealGradient(-80.0*eta*xi + 8.0*eta + 149.333333333333*eta*(xi*xi) - 66.3703703703704*eta*xi*xi*xi + 261.333333333333*xi*(eta*eta) - 149.333333333333*xi*eta*eta*eta - 32.0*eta*eta - 298.666666666667*eta*eta*xi*xi + 31.1111111111111*(eta*eta*eta) - 8.2962962962963*eta*eta*eta*eta, 56.0*eta*xi - 298.666666666667*eta*xi*xi + 298.666666666667*eta*(xi*xi*xi) - 12.0*xi - 46.6666666666667*xi*eta*eta + 8.2962962962963*xi*(eta*eta*eta) + 149.333333333333*(eta*eta)*(xi*xi) + 80.0*(xi*xi) - 136.888888888889*xi*xi*xi + 66.3703703703704*(xi*xi*xi*xi));
                case 6:
                  return sign * RealGradient(-56.0*eta*xi + 12.0*eta + 46.6666666666667*eta*(xi*xi) - 8.2962962962963*eta*xi*xi*xi + 298.666666666667*xi*(eta*eta) - 298.666666666667*xi*eta*eta*eta - 80.0*eta*eta - 149.333333333333*eta*eta*xi*xi + 136.888888888889*(eta*eta*eta) - 66.3703703703704*eta*eta*eta*eta, 80.0*eta*xi - 261.333333333333*eta*xi*xi + 149.333333333333*eta*(xi*xi*xi) - 8.0*xi - 149.333333333333*xi*eta*eta + 66.3703703703704*xi*(eta*eta*eta) + 298.666666666667*(eta*eta)*(xi*xi) + 32.0*(xi*xi) - 31.1111111111111*xi*xi*xi + 8.2962962962963*(xi*xi*xi*xi));
                case 7:
                  return sign * RealGradient(16.0*eta - 144.0*eta*eta + 336.0*(eta*eta*eta) - 224.0*eta*eta*eta*eta, 72.0*eta*xi - 4.0*xi - 252.0*xi*eta*eta + 224.0*xi*(eta*eta*eta));
                case 8:
                  return sign * RealGradient(16.0*eta - 144.0*eta*eta + 336.0*(eta*eta*eta) - 224.0*eta*eta*eta*eta, 72.0*eta*xi - 60.0*eta - 4.0*xi - 252.0*xi*eta*eta + 224.0*xi*(eta*eta*eta) + 4.0 + 180.0*(eta*eta) - 140.0*eta*eta*eta);
                case 9:
                  return sign * RealGradient(-12.4444444444444*eta*xi - 5.62962962962963*eta + 21.7777777777778*eta*(xi*xi) + 8.2962962962963*eta*(xi*xi*xi) + 43.5555555555556*xi*(eta*eta) + 24.8888888888889*xi*(eta*eta*eta) + 56.8888888888889*(eta*eta) - 124.444444444444*eta*eta*xi*xi - 140.0*eta*eta*eta + 91.2592592592593*(eta*eta*eta*eta), 16.0*eta*xi + 22.2222222222222*eta - 158.666666666667*eta*xi*xi + 124.444444444444*eta*(xi*xi*xi) - 4.14814814814815*xi + 124.444444444444*xi*(eta*eta) - 91.2592592592593*xi*eta*eta*eta - 1.18518518518519 - 75.5555555555556*eta*eta - 24.8888888888889*eta*eta*xi*xi + 11.5555555555556*(xi*xi) + 57.037037037037*(eta*eta*eta) + 2.07407407407407*(xi*xi*xi) - 8.2962962962963*xi*xi*xi*xi);
                case 10:
                  return sign * RealGradient(-19.5555555555556*eta*xi + 10.962962962963*eta - 49.7777777777778*eta*xi*xi + 66.3703703703704*eta*(xi*xi*xi) + 236.444444444444*xi*(eta*eta) - 248.888888888889*xi*eta*eta*eta - 88.8888888888889*eta*eta - 99.5555555555556*eta*eta*xi*xi + 168.0*(eta*eta*eta) - 91.2592592592593*eta*eta*eta*eta, 200.0*eta*xi - 42.2222222222222*eta - 261.333333333333*eta*xi*xi + 99.5555555555556*eta*(xi*xi*xi) + 2.81481481481481*xi - 320.444444444444*xi*eta*eta + 91.2592592592593*xi*(eta*eta*eta) + 2.51851851851852 + 95.5555555555556*(eta*eta) + 248.888888888889*(eta*eta)*(xi*xi) - 67.5555555555556*xi*xi - 57.037037037037*eta*eta*eta + 128.592592592593*(xi*xi*xi) - 66.3703703703704*xi*xi*xi*xi);
                case 11:
                  return sign * RealGradient(240.0*eta*xi - 40.0*eta - 420.0*eta*xi*xi + 224.0*eta*(xi*xi*xi) - 840.0*xi*eta*eta + 672.0*xi*(eta*eta*eta) + 240.0*(eta*eta) + 672.0*(eta*eta)*(xi*xi) - 420.0*eta*eta*eta + 224.0*(eta*eta*eta*eta), -720.0*eta*xi + 120.0*eta + 1260.0*eta*(xi*xi) - 672.0*eta*xi*xi*xi + 160.0*xi + 840.0*xi*(eta*eta) - 224.0*xi*eta*eta*eta - 16.0 - 240.0*eta*eta - 672.0*eta*eta*xi*xi - 480.0*xi*xi + 140.0*(eta*eta*eta) + 560.0*(xi*xi*xi) - 224.0*xi*xi*xi*xi);
                case 12:
                  return RealGradient(-3240.0*eta*xi + 960.0*eta + 3024.0*eta*(xi*xi) - 672.0*eta*xi*xi*xi + 9072.0*xi*(eta*eta) - 6048.0*xi*eta*eta*eta - 4320.0*eta*eta - 4032.0*eta*eta*xi*xi + 6048.0*(eta*eta*eta) - 2688.0*eta*eta*eta*eta, 2160.0*eta*xi - 6048.0*eta*xi*xi + 4032.0*eta*(xi*xi*xi) - 240.0*xi - 4536.0*xi*eta*eta + 2688.0*xi*(eta*eta*eta) + 6048.0*(eta*eta)*(xi*xi) + 1080.0*(xi*xi) - 1512.0*xi*xi*xi + 672.0*(xi*xi*xi*xi));
                case 13:
                  return RealGradient(2160.0*eta*xi - 240.0*eta - 4536.0*eta*xi*xi + 2688.0*eta*(xi*xi*xi) - 6048.0*xi*eta*eta + 4032.0*xi*(eta*eta*eta) + 1080.0*(eta*eta) + 6048.0*(eta*eta)*(xi*xi) - 1512.0*eta*eta*eta + 672.0*(eta*eta*eta*eta), -3240.0*eta*xi + 9072.0*eta*(xi*xi) - 6048.0*eta*xi*xi*xi + 960.0*xi + 3024.0*xi*(eta*eta) - 672.0*xi*eta*eta*eta - 4032.0*eta*eta*xi*xi - 4320.0*xi*xi + 6048.0*(xi*xi*xi) - 2688.0*xi*xi*xi*xi);
                case 14:
                  return RealGradient(-1944.0*eta*xi + 144.0*eta + 4536.0*eta*(xi*xi) - 2016.0*eta*xi*xi*xi + 2016.0*xi*(eta*eta) - 144.0*eta*eta - 4032.0*eta*eta*xi*xi, 432.0*eta*xi - 3024.0*eta*xi*xi + 4032.0*eta*(xi*xi*xi) - 216.0*xi + 1728.0*(xi*xi) - 3528.0*xi*xi*xi + 2016.0*(xi*xi*xi*xi));
                case 15:
                  return RealGradient(1152.0*eta*xi - 72.0*eta - 3528.0*eta*xi*xi + 2688.0*eta*(xi*xi*xi) - 1008.0*xi*eta*eta + 72.0*(eta*eta) + 2016.0*(eta*eta)*(xi*xi), -216.0*eta*xi + 1512.0*eta*(xi*xi) - 2016.0*eta*xi*xi*xi + 288.0*xi - 2304.0*xi*xi + 4704.0*(xi*xi*xi) - 2688.0*xi*xi*xi*xi);
                case 16:
                  return RealGradient(-216.0*eta*xi + 288.0*eta + 1512.0*xi*(eta*eta) - 2016.0*xi*eta*eta*eta - 2304.0*eta*eta + 4704.0*(eta*eta*eta) - 2688.0*eta*eta*eta*eta, 1152.0*eta*xi - 1008.0*eta*xi*xi - 72.0*xi - 3528.0*xi*eta*eta + 2688.0*xi*(eta*eta*eta) + 2016.0*(eta*eta)*(xi*xi) + 72.0*(xi*xi));
                case 17:
                  return RealGradient(432.0*eta*xi - 216.0*eta - 3024.0*xi*eta*eta + 4032.0*xi*(eta*eta*eta) + 1728.0*(eta*eta) - 3528.0*eta*eta*eta + 2016.0*(eta*eta*eta*eta), -1944.0*eta*xi + 2016.0*eta*(xi*xi) + 144.0*xi + 4536.0*xi*(eta*eta) - 2016.0*xi*eta*eta*eta - 4032.0*eta*eta*xi*xi - 144.0*xi*xi);
                case 18:
                  return RealGradient(-1332.0*eta*xi + 216.0*eta + 1638.0*eta*(xi*xi) - 504.0*eta*xi*xi*xi + 4788.0*xi*(eta*eta) - 3528.0*xi*eta*eta*eta - 1098.0*eta*eta - 3024.0*eta*eta*xi*xi + 1554.0*(eta*eta*eta) - 672.0*eta*eta*eta*eta, 1044.0*eta*xi - 4032.0*eta*xi*xi + 3024.0*eta*(xi*xi*xi) - 144.0*xi - 1638.0*xi*eta*eta + 672.0*xi*(eta*eta*eta) + 3528.0*(eta*eta)*(xi*xi) + 774.0*(xi*xi) - 1134.0*xi*xi*xi + 504.0*(xi*xi*xi*xi));
                case 19:
                  return RealGradient(1044.0*eta*xi - 144.0*eta - 1638.0*eta*xi*xi + 672.0*eta*(xi*xi*xi) - 4032.0*xi*eta*eta + 3024.0*xi*(eta*eta*eta) + 774.0*(eta*eta) + 3528.0*(eta*eta)*(xi*xi) - 1134.0*eta*eta*eta + 504.0*(eta*eta*eta*eta), -1332.0*eta*xi + 4788.0*eta*(xi*xi) - 3528.0*eta*xi*xi*xi + 216.0*xi + 1638.0*xi*(eta*eta) - 504.0*xi*eta*eta*eta - 3024.0*eta*eta*xi*xi - 1098.0*xi*xi + 1554.0*(xi*xi*xi) - 672.0*xi*xi*xi*xi);
                case 20:
                  return RealGradient(-216.0*eta*xi - 48.0*eta + 504.0*eta*(xi*xi) - 168.0*eta*xi*xi*xi - 756.0*xi*eta*eta + 1008.0*xi*(eta*eta*eta) + 720.0*(eta*eta) - 1344.0*eta*eta*eta + 672.0*(eta*eta*eta*eta), -360.0*eta*xi + 504.0*eta*(xi*xi) + 12.0*xi + 1008.0*xi*(eta*eta) - 672.0*xi*eta*eta*eta - 1008.0*eta*eta*xi*xi + 72.0*(xi*xi) - 252.0*xi*xi*xi + 168.0*(xi*xi*xi*xi));
                case 21:
                  return RealGradient(-216.0*eta*xi + 66.0*eta - 378.0*eta*xi*xi + 672.0*eta*(xi*xi*xi) + 2268.0*xi*(eta*eta) - 2016.0*xi*eta*eta*eta - 486.0*eta*eta - 1512.0*eta*eta*xi*xi + 756.0*(eta*eta*eta) - 336.0*eta*eta*eta*eta, 1188.0*eta*xi - 2772.0*eta*xi*xi + 1512.0*eta*(xi*xi*xi) + 6.0*xi - 1512.0*xi*eta*eta + 336.0*xi*(eta*eta*eta) + 2016.0*(eta*eta)*(xi*xi) - 468.0*xi*xi + 1134.0*(xi*xi*xi) - 672.0*xi*xi*xi*xi);
                case 22:
                  return RealGradient(1188.0*eta*xi + 6.0*eta - 1512.0*eta*xi*xi + 336.0*eta*(xi*xi*xi) - 2772.0*xi*eta*eta + 1512.0*xi*(eta*eta*eta) - 468.0*eta*eta + 2016.0*(eta*eta)*(xi*xi) + 1134.0*(eta*eta*eta) - 672.0*eta*eta*eta*eta, -216.0*eta*xi + 2268.0*eta*(xi*xi) - 2016.0*eta*xi*xi*xi + 66.0*xi - 378.0*xi*eta*eta + 672.0*xi*(eta*eta*eta) - 1512.0*eta*eta*xi*xi - 486.0*xi*xi + 756.0*(xi*xi*xi) - 336.0*xi*xi*xi*xi);
                case 23:
                  return RealGradient(-360.0*eta*xi + 12.0*eta + 1008.0*eta*(xi*xi) - 672.0*eta*xi*xi*xi + 504.0*xi*(eta*eta) + 72.0*(eta*eta) - 1008.0*eta*eta*xi*xi - 252.0*eta*eta*eta + 168.0*(eta*eta*eta*eta), -216.0*eta*xi - 756.0*eta*xi*xi + 1008.0*eta*(xi*xi*xi) - 48.0*xi + 504.0*xi*(eta*eta) - 168.0*xi*eta*eta*eta + 720.0*(xi*xi) - 1344.0*xi*xi*xi + 672.0*(xi*xi*xi*xi));
                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          } // end switch (type)
      } // end case FOURTH

      // quintic Nedelec (first kind) shape functions
    case FIFTH:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient(-156.25*eta - 75.0*xi + 3750.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 3281.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 2187.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 492.1875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5625.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 6562.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 3281.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 218.75 + 468.75*((eta + 1)*(eta + 1)) + 4921.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3281.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 738.28125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 131.25*((xi + 1)*(xi + 1)) - 5742.1875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2871.09375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 516.796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 546.875*(eta + 1)*(eta + 1)*(eta + 1) + 3828.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 861.328125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 87.5*(xi + 1)*(xi + 1)*(xi + 1) - 1914.0625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 344.53125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 273.4375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 430.6640625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 19.6875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 49.21875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 1:
                  return sign * RealGradient(23.681640625*eta + 22.3828125*xi - 1119.140625*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 1030.517578125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 664.794921875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 142.27294921875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1678.7109375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1958.49609375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 979.248046875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 176.2646484375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 44.169921875 - 71.044921875*(eta + 1)*(eta + 1) - 1545.7763671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 997.1923828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 213.409423828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 41.220703125*(xi + 1)*(xi + 1) + 1803.40576171875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 901.702880859375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 162.306518554688*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 82.8857421875*((eta + 1)*(eta + 1)*(eta + 1)) - 1163.39111328125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 248.977661132813*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 26.591796875*((xi + 1)*(xi + 1)*(xi + 1)) + 581.695556640625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 104.705200195313*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 41.44287109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 124.488830566406*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5.69091796875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 22.4079895019531*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 7.459716796875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 2:
                  return sign * RealGradient(-11.71875*eta - 13.125*xi + 656.25*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 902.34375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 738.28125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 184.5703125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 984.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1148.4375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 574.21875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 103.359375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 23.90625 + 35.15625*((eta + 1)*(eta + 1)) + 1353.515625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1107.421875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 276.85546875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 36.09375*((xi + 1)*(xi + 1)) - 1579.1015625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 789.55078125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 142.119140625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 41.015625*(eta + 1)*(eta + 1)*(eta + 1) + 1291.9921875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 322.998046875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 29.53125*(xi + 1)*(xi + 1)*(xi + 1) - 645.99609375*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 116.279296875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 20.5078125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 161.4990234375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 7.3828125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 29.06982421875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3.69140625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 3:
                  return sign * RealGradient(4.150390625*eta + 5.5078125*xi - 275.390625*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 456.298828125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 473.388671875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 142.27294921875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 413.0859375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 481.93359375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 240.966796875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 43.3740234375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 9.326171875 - 12.451171875*(eta + 1)*(eta + 1) - 684.4482421875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 710.0830078125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 213.409423828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 18.251953125*(xi + 1)*(xi + 1) + 798.52294921875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 399.261474609375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 71.8670654296875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 14.5263671875*((eta + 1)*(eta + 1)*(eta + 1)) - 828.43017578125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 248.977661132813*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 18.935546875*((xi + 1)*(xi + 1)*(xi + 1)) + 414.215087890625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 74.5587158203125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 7.26318359375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 124.488830566406*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5.69091796875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 22.4079895019531*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1.307373046875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 4:
                  return sign * RealGradient(-31.25*eta - 30.0*xi + 1500.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 1968.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1750.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 492.1875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2250.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 2625.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 1312.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 236.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 58.75 + 93.75*((eta + 1)*(eta + 1)) + 2953.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2625.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 738.28125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 78.75*((xi + 1)*(xi + 1)) - 3445.3125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1722.65625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 109.375*(eta + 1)*(eta + 1)*(eta + 1) + 3062.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 861.328125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 70.0*(xi + 1)*(xi + 1)*(xi + 1) - 1531.25*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 275.625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 54.6875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 430.6640625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 19.6875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 9.84375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 5:
                  return sign * RealGradient(0, 31.25*xi - 750.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 2250.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 3937.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 2625.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 656.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 437.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 98.4375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 31.25 - 1968.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3445.3125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2296.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 187.5*(xi + 1)*(xi + 1) + 1312.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 295.3125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 2296.875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1531.25*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 328.125*((xi + 1)*(xi + 1)*(xi + 1)) + 516.796875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 218.75*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 6:
                  return sign * RealGradient(0, -4.736328125*xi + 223.828125*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 671.484375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1175.09765625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 783.3984375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 176.2646484375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 206.103515625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 132.958984375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 28.45458984375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 4.736328125 + 618.310546875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1082.04345703125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 721.3623046875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 162.306518554688*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 28.41796875*((xi + 1)*(xi + 1)) - 398.876953125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 85.36376953125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 698.03466796875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 465.3564453125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 104.705200195313*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 49.7314453125*(xi + 1)*(xi + 1)*(xi + 1) - 149.386596679688*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 99.591064453125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 22.4079895019531*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 33.154296875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 7.459716796875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 7:
                  return sign * RealGradient(0, 2.34375*xi - 131.25*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 393.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 689.0625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 459.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 103.359375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 180.46875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 147.65625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 36.9140625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2.34375 - 541.40625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 947.4609375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 631.640625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 142.119140625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.0625*(xi + 1)*(xi + 1) + 442.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 110.7421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 775.1953125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 116.279296875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 24.609375*((xi + 1)*(xi + 1)*(xi + 1)) + 193.798828125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 129.19921875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 29.06982421875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 16.40625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3.69140625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 8:
                  return sign * RealGradient(0, -0.830078125*xi + 55.078125*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 165.234375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 289.16015625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 192.7734375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 43.3740234375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 91.259765625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 94.677734375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 28.45458984375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 0.830078125 + 273.779296875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 479.11376953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 319.4091796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 71.8670654296875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4.98046875*((xi + 1)*(xi + 1)) - 284.033203125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 85.36376953125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 497.05810546875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 331.3720703125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 74.5587158203125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 8.7158203125*(xi + 1)*(xi + 1)*(xi + 1) - 149.386596679688*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 99.591064453125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 22.4079895019531*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5.810546875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1.307373046875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 9:
                  return sign * RealGradient(0, 6.25*xi - 300.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 900.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1575.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1050.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 236.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 393.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 350.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 98.4375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 6.25 - 1181.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2067.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 37.5*(xi + 1)*(xi + 1) + 1050.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 295.3125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1837.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1225.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 275.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1)*(xi + 1)) + 516.796875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 43.75*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9.84375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 10:
                  return sign * RealGradient(-6.25*eta + 300.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 393.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 350.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 98.4375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 900.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1575.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 1050.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 236.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 6.25 + 37.5*((eta + 1)*(eta + 1)) + 1181.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1050.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 295.3125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2067.1875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1378.125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 65.625*(eta + 1)*(eta + 1)*(eta + 1) + 1837.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1225.0*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 275.625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 43.75*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 344.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 9.84375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 11:
                  return sign * RealGradient(0.830078125*eta - 55.078125*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 91.259765625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 94.677734375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 28.45458984375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 165.234375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 289.16015625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 192.7734375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 43.3740234375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 0.830078125 - 4.98046875*(eta + 1)*(eta + 1) - 273.779296875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 284.033203125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 85.36376953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 479.11376953125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 319.4091796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 71.8670654296875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 8.7158203125*((eta + 1)*(eta + 1)*(eta + 1)) - 497.05810546875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 149.386596679688*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 331.3720703125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 74.5587158203125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 5.810546875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 99.591064453125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 22.4079895019531*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1.307373046875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 12:
                  return sign * RealGradient(-2.34375*eta + 131.25*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 180.46875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 147.65625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 36.9140625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 393.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 689.0625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 459.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 103.359375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2.34375 + 14.0625*((eta + 1)*(eta + 1)) + 541.40625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 442.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 110.7421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 947.4609375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 631.640625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 142.119140625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 24.609375*(eta + 1)*(eta + 1)*(eta + 1) + 775.1953125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 193.798828125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 516.796875*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 116.279296875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 16.40625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 129.19921875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 29.06982421875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3.69140625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 13:
                  return sign * RealGradient(4.736328125*eta - 223.828125*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 206.103515625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 132.958984375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 28.45458984375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 671.484375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1175.09765625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 783.3984375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 176.2646484375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4.736328125 - 28.41796875*(eta + 1)*(eta + 1) - 618.310546875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.876953125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 85.36376953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1082.04345703125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 721.3623046875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 162.306518554688*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 49.7314453125*((eta + 1)*(eta + 1)*(eta + 1)) - 698.03466796875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 149.386596679688*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 465.3564453125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 104.705200195313*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 33.154296875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 99.591064453125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 22.4079895019531*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 7.459716796875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 14:
                  return sign * RealGradient(-31.25*eta + 750.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 656.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 437.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 98.4375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2250.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 3937.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 2625.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 31.25 + 187.5*((eta + 1)*(eta + 1)) + 1968.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1312.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 295.3125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 3445.3125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2296.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 516.796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 328.125*(eta + 1)*(eta + 1)*(eta + 1) + 2296.875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1531.25*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 344.53125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 218.75*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 344.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 49.21875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 15:
                  return sign * RealGradient(0, 30.0*eta + 31.25*xi - 1500.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 2250.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 2625.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1312.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 236.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1968.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1750.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 492.1875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 58.75 - 78.75*(eta + 1)*(eta + 1) - 2953.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3445.3125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 93.75*(xi + 1)*(xi + 1) + 2625.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 738.28125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 70.0*((eta + 1)*(eta + 1)*(eta + 1)) - 3062.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1531.25*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 275.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 109.375*((xi + 1)*(xi + 1)*(xi + 1)) + 861.328125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 19.6875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 430.6640625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 54.6875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9.84375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 16:
                  return sign * RealGradient(0, -5.5078125*eta - 4.150390625*xi + 275.390625*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 413.0859375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 481.93359375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 240.966796875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 43.3740234375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 456.298828125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 473.388671875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 142.27294921875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 9.326171875 + 18.251953125*((eta + 1)*(eta + 1)) + 684.4482421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 798.52294921875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 399.261474609375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 71.8670654296875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 12.451171875*((xi + 1)*(xi + 1)) - 710.0830078125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 213.409423828125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 18.935546875*(eta + 1)*(eta + 1)*(eta + 1) + 828.43017578125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 414.215087890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 74.5587158203125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.5263671875*(xi + 1)*(xi + 1)*(xi + 1) - 248.977661132813*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5.69091796875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 124.488830566406*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 22.4079895019531*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.26318359375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1.307373046875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 17:
                  return sign * RealGradient(0, 13.125*eta + 11.71875*xi - 656.25*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 984.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1148.4375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 574.21875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 103.359375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 902.34375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 738.28125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 184.5703125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 23.90625 - 36.09375*(eta + 1)*(eta + 1) - 1353.515625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1579.1015625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 789.55078125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 142.119140625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 35.15625*(xi + 1)*(xi + 1) + 1107.421875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 276.85546875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 29.53125*((eta + 1)*(eta + 1)*(eta + 1)) - 1291.9921875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 645.99609375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 116.279296875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 41.015625*((xi + 1)*(xi + 1)*(xi + 1)) + 322.998046875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 7.3828125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 161.4990234375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 29.06982421875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 20.5078125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3.69140625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 18:
                  return sign * RealGradient(0, -22.3828125*eta - 23.681640625*xi + 1119.140625*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 1678.7109375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1958.49609375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 979.248046875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 176.2646484375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1030.517578125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 664.794921875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 142.27294921875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 44.169921875 + 41.220703125*((eta + 1)*(eta + 1)) + 1545.7763671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1803.40576171875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 901.702880859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 162.306518554688*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 71.044921875*((xi + 1)*(xi + 1)) - 997.1923828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 213.409423828125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 26.591796875*(eta + 1)*(eta + 1)*(eta + 1) + 1163.39111328125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 581.695556640625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 104.705200195313*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 82.8857421875*(xi + 1)*(xi + 1)*(xi + 1) - 248.977661132813*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5.69091796875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 124.488830566406*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 22.4079895019531*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 41.44287109375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 7.459716796875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 19:
                  return sign * RealGradient(0, 75.0*eta + 156.25*xi - 3750.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 5625.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 6562.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 3281.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3281.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 2187.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 492.1875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 218.75 - 131.25*(eta + 1)*(eta + 1) - 4921.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5742.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2871.09375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 468.75*(xi + 1)*(xi + 1) + 3281.25*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 738.28125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 87.5*((eta + 1)*(eta + 1)*(eta + 1)) - 3828.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1914.0625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 546.875*((xi + 1)*(xi + 1)*(xi + 1)) + 861.328125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 19.6875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 430.6640625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 273.4375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 20:
                  return RealGradient(0, 121.875*xi - 2925.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 5287.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 6478.125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 3281.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2559.375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1706.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 383.90625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 121.875 - 4626.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5668.359375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2871.09375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 440.625*(xi + 1)*(xi + 1) + 3084.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 693.984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3778.90625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1914.0625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 539.84375*((xi + 1)*(xi + 1)*(xi + 1)) + 850.25390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 430.6640625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 273.4375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 21:
                  return RealGradient(0, -25.0*xi + 600.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 2081.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 3853.125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 2625.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 590.625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 525.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 350.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 78.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 25.0 + 1821.09375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3371.484375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2296.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 173.4375*((xi + 1)*(xi + 1)) - 1214.0625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 273.1640625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2247.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1531.25*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 344.53125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 321.09375*(xi + 1)*(xi + 1)*(xi + 1) - 505.72265625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 344.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 218.75*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 49.21875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 22:
                  return RealGradient(0, 56.25*xi - 1350.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 843.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 253.125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1181.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 787.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 177.1875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 56.25 - 738.28125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 221.484375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 70.3125*(xi + 1)*(xi + 1) + 492.1875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 110.7421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 147.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 21.09375*((xi + 1)*(xi + 1)*(xi + 1)) + 33.22265625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                case 23:
                  return RealGradient(0, -28.125*xi + 675.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 675.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 253.125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 393.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 88.59375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 28.125 + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 221.484375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 56.25*((xi + 1)*(xi + 1)) - 393.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 147.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 21.09375*(xi + 1)*(xi + 1)*(xi + 1) - 33.22265625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                case 24:
                  return RealGradient(-121.875*eta + 2925.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 2559.375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1706.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 383.90625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5287.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 6478.125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 3281.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 121.875 + 440.625*((eta + 1)*(eta + 1)) + 4626.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3084.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 693.984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 5668.359375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2871.09375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 516.796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 539.84375*(eta + 1)*(eta + 1)*(eta + 1) + 3778.90625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 850.25390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1914.0625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 344.53125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 273.4375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 430.6640625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 49.21875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 25:
                  return RealGradient(25.0*eta - 600.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 525.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 350.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 78.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2081.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 3853.125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 2625.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 590.625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 25.0 - 173.4375*(eta + 1)*(eta + 1) - 1821.09375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1214.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 273.1640625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3371.484375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 2296.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 516.796875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 321.09375*((eta + 1)*(eta + 1)*(eta + 1)) - 2247.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 505.72265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1531.25*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 344.53125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 218.75*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 49.21875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 26:
                  return RealGradient(-56.25*eta + 1350.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 1181.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 787.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 177.1875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 843.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 253.125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 56.25 + 70.3125*((eta + 1)*(eta + 1)) + 738.28125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 492.1875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 110.7421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 221.484375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 21.09375*(eta + 1)*(eta + 1)*(eta + 1) + 147.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 33.22265625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                case 27:
                  return RealGradient(28.125*eta - 675.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 590.625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 393.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 88.59375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 675.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 253.125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 28.125 - 56.25*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 393.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 221.484375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 21.09375*((eta + 1)*(eta + 1)*(eta + 1)) - 147.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 33.22265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                case 28:
                  return RealGradient(-24.375*eta + 1170.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 1535.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1365.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 383.90625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2115.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 2591.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 1312.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 236.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 24.375 + 88.125*((eta + 1)*(eta + 1)) + 2775.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2467.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 693.984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 3401.015625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1722.65625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 107.96875*(eta + 1)*(eta + 1)*(eta + 1) + 3023.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 850.25390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1531.25*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 275.625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 54.6875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 430.6640625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 9.84375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 29:
                  return RealGradient(5.0*eta - 240.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 315.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 280.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 78.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 832.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1541.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1050.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 236.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5.0 - 34.6875*(eta + 1)*(eta + 1) - 1092.65625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 971.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 273.1640625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2022.890625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1378.125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 310.078125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 64.21875*((eta + 1)*(eta + 1)*(eta + 1)) - 1798.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 505.72265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1225.0*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 275.625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 43.75*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 344.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 9.84375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 30:
                  return RealGradient(-11.25*eta + 540.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 708.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 630.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 177.1875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 337.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 101.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 11.25 + 14.0625*((eta + 1)*(eta + 1)) + 442.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 393.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 110.7421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 132.890625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 4.21875*(eta + 1)*(eta + 1)*(eta + 1) + 118.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 33.22265625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                case 31:
                  return RealGradient(5.625*eta - 270.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 354.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 315.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 88.59375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 270.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 101.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 5.625 - 11.25*(eta + 1)*(eta + 1) - 354.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 315.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 132.890625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 4.21875*((eta + 1)*(eta + 1)*(eta + 1)) - 118.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 33.22265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                case 32:
                  return RealGradient(0, 24.375*xi - 1170.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 2115.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 2591.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1312.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 236.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1535.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1365.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 383.90625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 24.375 - 2775.9375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3401.015625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 88.125*(xi + 1)*(xi + 1) + 2467.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 693.984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3023.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1531.25*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 275.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 107.96875*((xi + 1)*(xi + 1)*(xi + 1)) + 850.25390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 430.6640625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 77.51953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 54.6875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9.84375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 33:
                  return RealGradient(0, -5.0*xi + 240.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 832.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1541.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1050.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 236.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 315.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 280.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 78.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 5.0 + 1092.65625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2022.890625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 310.078125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 34.6875*((xi + 1)*(xi + 1)) - 971.25*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 273.1640625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1798.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1225.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 275.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 64.21875*(xi + 1)*(xi + 1)*(xi + 1) - 505.72265625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 344.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 77.51953125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 43.75*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 9.84375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 34:
                  return RealGradient(0, 11.25*xi - 540.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 337.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 101.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 708.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 630.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 177.1875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 11.25 - 442.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 132.890625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 14.0625*(xi + 1)*(xi + 1) + 393.75*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 110.7421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 118.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4.21875*((xi + 1)*(xi + 1)*(xi + 1)) + 33.22265625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                case 35:
                  return RealGradient(0, -5.625*xi + 270.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 270.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 101.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 354.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 315.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 88.59375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 5.625 + 354.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 132.890625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 11.25*((xi + 1)*(xi + 1)) - 315.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 118.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4.21875*(xi + 1)*(xi + 1)*(xi + 1) - 33.22265625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                case 36:
                  return RealGradient(36.0*eta - 288.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 60.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 594.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 765.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 393.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 70.875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 36.0 - 148.5*(eta + 1)*(eta + 1) - 123.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 159.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 82.03125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 191.25*((eta + 1)*(eta + 1)*(eta + 1)) - 98.4375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 17.71875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 37:
                  return RealGradient(12.0*eta - 192.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 60.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 396.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 510.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 262.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 47.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 12.0 - 49.5*(eta + 1)*(eta + 1) - 123.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 159.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 82.03125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 63.75*((eta + 1)*(eta + 1)*(eta + 1)) - 32.8125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5.90625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 38:
                  return RealGradient(-6.0*eta + 120.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 247.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 318.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 164.0625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 29.53125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 6.0 + 24.75*((eta + 1)*(eta + 1)) + 61.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 79.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 41.015625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 7.3828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 31.875*(eta + 1)*(eta + 1)*(eta + 1) + 16.40625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2.953125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 39:
                  return RealGradient(0, 36.0*xi - 288.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 594.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 765.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 393.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 70.875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 60.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 36.0 - 123.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 159.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 82.03125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 148.5*(xi + 1)*(xi + 1) + 191.25*((xi + 1)*(xi + 1)*(xi + 1)) - 98.4375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 17.71875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 40:
                  return RealGradient(0, 12.0*xi - 192.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 396.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 510.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 262.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 47.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 60.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 12.0 - 123.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 159.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 82.03125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 49.5*(xi + 1)*(xi + 1) + 63.75*((xi + 1)*(xi + 1)*(xi + 1)) - 32.8125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5.90625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 41:
                  return RealGradient(0, -6.0*xi + 120.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 247.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 318.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 164.0625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 29.53125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 30.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 6.0 + 61.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 79.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 41.015625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 7.3828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 24.75*((xi + 1)*(xi + 1)) - 31.875*(xi + 1)*(xi + 1)*(xi + 1) + 16.40625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2.953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 42:
                  return RealGradient(0, -4.5*xi + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 216.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 450.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 315.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 70.875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 7.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 4.5 + 45.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 54.0*((xi + 1)*(xi + 1)) - 112.5*(xi + 1)*(xi + 1)*(xi + 1) + 78.75*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 17.71875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 43:
                  return RealGradient(0, -1.5*xi + 24.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 144.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 300.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 210.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 47.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 7.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 1.5 + 45.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 18.0*((xi + 1)*(xi + 1)) - 37.5*(xi + 1)*(xi + 1)*(xi + 1) + 26.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 5.90625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                case 44:
                  return RealGradient(0, 0.75*xi - 15.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 90.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 187.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 29.53125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 0.75 - 22.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 46.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 32.8125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.3828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 9.0*(xi + 1)*(xi + 1) + 18.75*((xi + 1)*(xi + 1)*(xi + 1)) - 13.125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2.953125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                case 45:
                  return RealGradient(-4.5*eta + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 7.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 216.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 450.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 315.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 70.875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 4.5 + 54.0*((eta + 1)*(eta + 1)) + 45.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 112.5*(eta + 1)*(eta + 1)*(eta + 1) + 78.75*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 17.71875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 46:
                  return RealGradient(-1.5*eta + 24.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 7.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 144.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 300.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 210.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 47.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1.5 + 18.0*((eta + 1)*(eta + 1)) + 45.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 37.5*(eta + 1)*(eta + 1)*(eta + 1) + 26.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 5.90625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 47:
                  return RealGradient(0.75*eta - 15.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 3.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 90.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 187.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 29.53125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 0.75 - 9.0*(eta + 1)*(eta + 1) - 22.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 46.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 32.8125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 7.3828125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 18.75*((eta + 1)*(eta + 1)*(eta + 1)) - 13.125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2.953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 48:
                  return RealGradient(0, 18.0*xi - 144.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 81.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 22.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 18.0 - 16.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 4.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 20.25*(xi + 1)*(xi + 1) + 5.625*((xi + 1)*(xi + 1)*(xi + 1)));
                case 49:
                  return RealGradient(0, -4.5*xi + 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 54.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 22.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 7.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 4.5 + 11.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 4.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 13.5*((xi + 1)*(xi + 1)) - 5.625*(xi + 1)*(xi + 1)*(xi + 1));
                case 50:
                  return RealGradient(-18.0*eta + 144.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 81.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 22.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 18.0 + 20.25*((eta + 1)*(eta + 1)) + 16.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 4.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 5.625*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 51:
                  return RealGradient(4.5*eta - 36.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 7.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 54.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 22.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 4.5 - 13.5*(eta + 1)*(eta + 1) - 11.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 4.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 5.625*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 52:
                  return RealGradient(-6.0*eta + 96.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 30.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 54.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 15.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.0 + 6.75*((eta + 1)*(eta + 1)) + 16.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 4.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1.875*(eta + 1)*(eta + 1)*(eta + 1), 0);
                case 53:
                  return RealGradient(1.5*eta - 24.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 7.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 36.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 15.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1.5 - 4.5*(eta + 1)*(eta + 1) - 11.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 4.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 1.875*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 54:
                  return RealGradient(0, 6.0*xi - 96.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) + 54.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 15.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 6.0 - 16.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 4.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 6.75*(xi + 1)*(xi + 1) + 1.875*((xi + 1)*(xi + 1)*(xi + 1)));
                case 55:
                  return RealGradient(0, -1.5*xi + 24.0*(0.5*eta + 0.5)*(0.5*xi + 0.5) - 36.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 15.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 7.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 1.5 + 11.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 4.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4.5*((xi + 1)*(xi + 1)) - 1.875*(xi + 1)*(xi + 1)*(xi + 1));
                case 56:
                  return RealGradient(2.0*eta + 2.0 - 2.25*(eta + 1)*(eta + 1) + 0.625*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                case 57:
                  return RealGradient(0, 2.0*xi + 2.0 - 2.25*(xi + 1)*(xi + 1) + 0.625*((xi + 1)*(xi + 1)*(xi + 1)));
                case 58:
                  return RealGradient(0, -0.5*xi - 0.5 + 1.5*((xi + 1)*(xi + 1)) - 0.625*(xi + 1)*(xi + 1)*(xi + 1));
                case 59:
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
                  return sign * RealGradient(2800.0*eta*xi - 375.0*eta - 6300.0*eta*xi*xi + 5040.0*eta*(xi*xi*xi) - 1050.0*eta*xi*xi*xi*xi - 300.0*xi - 8400.0*xi*eta*eta + 10080.0*xi*(eta*eta*eta) - 4200.0*xi*eta*eta*eta*eta + 25.0 + 1750.0*(eta*eta) + 11340.0*(eta*eta)*(xi*xi) - 4200.0*eta*eta*xi*xi*xi + 1050.0*(xi*xi) - 6300.0*xi*xi*eta*eta*eta - 3500.0*eta*eta*eta - 1400.0*xi*xi*xi + 3150.0*(eta*eta*eta*eta) + 630.0*(xi*xi*xi*xi) - 1050.0*eta*eta*eta*eta*eta, -700.0*eta*xi + 4200.0*eta*(xi*xi) - 7560.0*eta*xi*xi*xi + 4200.0*eta*(xi*xi*xi*xi) + 75.0*xi + 2100.0*xi*(eta*eta) - 2520.0*xi*eta*eta*eta + 1050.0*xi*(eta*eta*eta*eta) - 7560.0*eta*eta*xi*xi + 6300.0*(eta*eta)*(xi*xi*xi) - 700.0*xi*xi + 4200.0*(xi*xi)*(eta*eta*eta) + 2100.0*(xi*xi*xi) - 2520.0*xi*xi*xi*xi + 1050.0*(xi*xi*xi*xi*xi));
                case 1:
                  return sign * RealGradient(-759.0625*eta*xi + 5.4296875*eta + 2127.890625*eta*(xi*xi) - 1706.25*eta*xi*xi*xi + 303.515625*eta*(xi*xi*xi*xi) + 89.53125*xi + 1811.25*xi*(eta*eta) - 1594.6875*xi*eta*eta*eta + 442.96875*xi*(eta*eta*eta*eta) - 3.7890625 + 175.546875*(eta*eta) - 3765.234375*eta*eta*xi*xi + 1525.78125*(eta*eta)*(xi*xi*xi) - 329.765625*xi*xi + 1993.359375*(xi*xi)*(eta*eta*eta) - 664.453125*eta*eta*eta + 425.46875*(xi*xi*xi) + 819.4921875*(eta*eta*eta*eta) - 182.109375*xi*xi*xi*xi - 332.2265625*eta*eta*eta*eta*eta, 95.15625*eta*xi - 1230.46875*eta*xi*xi + 2628.28125*eta*(xi*xi*xi) - 1525.78125*eta*xi*xi*xi*xi - 18.7109375*xi + 44.296875*xi*(eta*eta) - 442.96875*xi*eta*eta*eta + 332.2265625*xi*(eta*eta*eta*eta) + 1594.6875*(eta*eta)*(xi*xi) - 1993.359375*eta*eta*xi*xi*xi + 214.375*(xi*xi) - 442.96875*xi*xi*eta*eta*eta - 652.421875*xi*xi*xi + 759.609375*(xi*xi*xi*xi) - 303.515625*xi*xi*xi*xi*xi);
                case 2:
                  return sign * RealGradient(6.875*eta - 813.75*eta*xi*xi + 1260.0*eta*(xi*xi*xi) - 393.75*eta*xi*xi*xi*xi - 52.5*xi + 840.0*xi*(eta*eta) - 1575.0*xi*eta*eta*eta + 787.5*xi*(eta*eta*eta*eta) + 1.875 - 43.75*eta*eta + 236.25*(eta*eta)*(xi*xi) - 787.5*eta*eta*xi*xi*xi + 288.75*(xi*xi) + 393.75*(xi*xi)*(eta*eta*eta) + 8.75*(eta*eta*eta) - 472.5*xi*xi*xi + 91.875*(eta*eta*eta*eta) + 236.25*(xi*xi*xi*xi) - 65.625*eta*eta*eta*eta*eta, 17.5*eta*xi + 52.5*eta*(xi*xi) - 787.5*eta*xi*xi*xi + 787.5*eta*(xi*xi*xi*xi) + 10.625*xi - 131.25*xi*eta*eta + 52.5*xi*(eta*eta*eta) + 65.625*xi*(eta*eta*eta*eta) + 945.0*(eta*eta)*(xi*xi) - 393.75*eta*eta*xi*xi*xi - 157.5*xi*xi - 787.5*xi*xi*eta*eta*eta + 621.25*(xi*xi*xi) - 866.25*xi*xi*xi*xi + 393.75*(xi*xi*xi*xi*xi));
                case 3:
                  return sign * RealGradient(115.9375*eta*xi - 6.4453125*eta - 260.859375*eta*xi*xi - 236.25*eta*xi*xi*xi + 303.515625*eta*(xi*xi*xi*xi) + 22.03125*xi - 288.75*xi*eta*eta - 19.6875*xi*eta*eta*eta + 180.46875*xi*(eta*eta*eta*eta) - 0.6640625 + 9.296875*(eta*eta) + 1196.015625*(eta*eta)*(xi*xi) - 311.71875*eta*eta*xi*xi*xi - 146.015625*xi*xi - 762.890625*xi*xi*eta*eta*eta + 9.296875*(eta*eta*eta) + 302.96875*(xi*xi*xi) - 7.3828125*eta*eta*eta*eta - 182.109375*xi*xi*xi*xi - 4.1015625*eta*eta*eta*eta*eta, -27.34375*eta*xi + 292.03125*eta*(xi*xi) - 679.21875*eta*xi*xi*xi + 311.71875*eta*(xi*xi*xi*xi) - 4.3359375*xi + 18.046875*xi*(eta*eta) + 29.53125*xi*(eta*eta*eta) + 4.1015625*xi*(eta*eta*eta*eta) - 295.3125*eta*eta*xi*xi + 762.890625*(eta*eta)*(xi*xi*xi) + 74.375*(xi*xi) - 180.46875*xi*xi*eta*eta*eta - 346.171875*xi*xi*xi + 575.859375*(xi*xi*xi*xi) - 303.515625*xi*xi*xi*xi*xi);
                case 4:
                  return sign * RealGradient(140.0*eta*xi - 5.0*eta - 840.0*eta*xi*xi + 1680.0*eta*(xi*xi*xi) - 1050.0*eta*xi*xi*xi*xi - 120.0*xi + 5.0 + 630.0*(xi*xi) - 1120.0*xi*xi*xi + 630.0*(xi*xi*xi*xi), 25.0*xi - 350.0*xi*xi + 1400.0*(xi*xi*xi) - 2100.0*xi*xi*xi*xi + 1050.0*(xi*xi*xi*xi*xi));
                case 5:
                  return sign * RealGradient(140.0*eta*xi - 5.0*eta - 840.0*eta*xi*xi + 1680.0*eta*(xi*xi*xi) - 1050.0*eta*xi*xi*xi*xi, 25.0*xi - 350.0*xi*xi + 1400.0*(xi*xi*xi) - 2100.0*xi*xi*xi*xi + 1050.0*(xi*xi*xi*xi*xi));
                case 6:
                  return sign * RealGradient(183.75*eta*xi - 10.0*eta - 708.75*eta*xi*xi + 885.9375*eta*(xi*xi*xi) - 332.2265625*eta*xi*xi*xi*xi - 787.5*xi*eta*eta + 767.8125*xi*(eta*eta*eta) - 196.875*xi*eta*eta*eta*eta + 48.125*(eta*eta) + 2392.03125*(eta*eta)*(xi*xi) - 1771.875*eta*eta*xi*xi*xi - 1328.90625*xi*xi*eta*eta*eta - 61.25*eta*eta*eta + 27.890625*(eta*eta*eta*eta) - 4.1015625*eta*eta*eta*eta*eta, -113.75*eta*xi + 1102.5*eta*(xi*xi) - 2657.8125*eta*xi*xi*xi + 1771.875*eta*(xi*xi*xi*xi) + 20.0*xi + 131.25*xi*(eta*eta) - 45.9375*xi*eta*eta*eta + 4.1015625*xi*(eta*eta*eta*eta) - 974.53125*eta*eta*xi*xi + 1328.90625*(eta*eta)*(xi*xi*xi) - 223.125*xi*xi + 196.875*(xi*xi)*(eta*eta*eta) + 708.75*(xi*xi*xi) - 841.640625*xi*xi*xi*xi + 332.2265625*(xi*xi*xi*xi*xi));
                case 7:
                  return sign * RealGradient(175.0*eta*xi - 15.0*eta - 420.0*eta*xi*xi + 315.0*eta*(xi*xi*xi) - 65.625*eta*xi*xi*xi*xi - 1260.0*xi*eta*eta + 2205.0*xi*(eta*eta*eta) - 1050.0*xi*eta*eta*eta*eta + 122.5*(eta*eta) + 2362.5*(eta*eta)*(xi*xi) - 1050.0*eta*eta*xi*xi*xi - 2362.5*xi*xi*eta*eta*eta - 280.0*eta*eta*eta + 236.25*(eta*eta*eta*eta) - 65.625*eta*eta*eta*eta*eta, -175.0*eta*xi + 1260.0*eta*(xi*xi) - 2205.0*eta*xi*xi*xi + 1050.0*eta*(xi*xi*xi*xi) + 15.0*xi + 420.0*xi*(eta*eta) - 315.0*xi*eta*eta*eta + 65.625*xi*(eta*eta*eta*eta) - 2362.5*eta*eta*xi*xi + 2362.5*(eta*eta)*(xi*xi*xi) - 122.5*xi*xi + 1050.0*(xi*xi)*(eta*eta*eta) + 280.0*(xi*xi*xi) - 236.25*xi*xi*xi*xi + 65.625*(xi*xi*xi*xi*xi));
                case 8:
                  return sign * RealGradient(113.75*eta*xi - 20.0*eta - 131.25*eta*xi*xi + 45.9375*eta*(xi*xi*xi) - 4.1015625*eta*xi*xi*xi*xi - 1102.5*xi*eta*eta + 2657.8125*xi*(eta*eta*eta) - 1771.875*xi*eta*eta*eta*eta + 223.125*(eta*eta) + 974.53125*(eta*eta)*(xi*xi) - 196.875*eta*eta*xi*xi*xi - 1328.90625*xi*xi*eta*eta*eta - 708.75*eta*eta*eta + 841.640625*(eta*eta*eta*eta) - 332.2265625*eta*eta*eta*eta*eta, -183.75*eta*xi + 787.5*eta*(xi*xi) - 767.8125*eta*xi*xi*xi + 196.875*eta*(xi*xi*xi*xi) + 10.0*xi + 708.75*xi*(eta*eta) - 885.9375*xi*eta*eta*eta + 332.2265625*xi*(eta*eta*eta*eta) - 2392.03125*eta*eta*xi*xi + 1328.90625*(eta*eta)*(xi*xi*xi) - 48.125*xi*xi + 1771.875*(xi*xi)*(eta*eta*eta) + 61.25*(xi*xi*xi) - 27.890625*xi*xi*xi*xi + 4.1015625*(xi*xi*xi*xi*xi));
                case 9:
                  return sign * RealGradient(-25.0*eta + 350.0*(eta*eta) - 1400.0*eta*eta*eta + 2100.0*(eta*eta*eta*eta) - 1050.0*eta*eta*eta*eta*eta, -140.0*eta*xi + 5.0*xi + 840.0*xi*(eta*eta) - 1680.0*xi*eta*eta*eta + 1050.0*xi*(eta*eta*eta*eta));
                case 10:
                  return sign * RealGradient(-25.0*eta + 350.0*(eta*eta) - 1400.0*eta*eta*eta + 2100.0*(eta*eta*eta*eta) - 1050.0*eta*eta*eta*eta*eta, -140.0*eta*xi + 120.0*eta + 5.0*xi + 840.0*xi*(eta*eta) - 1680.0*xi*eta*eta*eta + 1050.0*xi*(eta*eta*eta*eta) - 5.0 - 630.0*eta*eta + 1120.0*(eta*eta*eta) - 630.0*eta*eta*eta*eta);
                case 11:
                  return sign * RealGradient(27.34375*eta*xi + 4.3359375*eta - 18.046875*eta*xi*xi - 29.53125*eta*xi*xi*xi - 4.1015625*eta*xi*xi*xi*xi - 292.03125*xi*eta*eta + 679.21875*xi*(eta*eta*eta) - 311.71875*xi*eta*eta*eta*eta - 74.375*eta*eta + 295.3125*(eta*eta)*(xi*xi) + 180.46875*(eta*eta)*(xi*xi*xi) - 762.890625*xi*xi*eta*eta*eta + 346.171875*(eta*eta*eta) - 575.859375*eta*eta*eta*eta + 303.515625*(eta*eta*eta*eta*eta), -115.9375*eta*xi - 22.03125*eta + 288.75*eta*(xi*xi) + 19.6875*eta*(xi*xi*xi) - 180.46875*eta*xi*xi*xi*xi + 6.4453125*xi + 260.859375*xi*(eta*eta) + 236.25*xi*(eta*eta*eta) - 303.515625*xi*eta*eta*eta*eta + 0.6640625 + 146.015625*(eta*eta) - 1196.015625*eta*eta*xi*xi + 762.890625*(eta*eta)*(xi*xi*xi) - 9.296875*xi*xi + 311.71875*(xi*xi)*(eta*eta*eta) - 302.96875*eta*eta*eta - 9.296875*xi*xi*xi + 182.109375*(eta*eta*eta*eta) + 7.3828125*(xi*xi*xi*xi) + 4.1015625*(xi*xi*xi*xi*xi));
                case 12:
                  return sign * RealGradient(-17.5*eta*xi - 10.625*eta + 131.25*eta*(xi*xi) - 52.5*eta*xi*xi*xi - 65.625*eta*xi*xi*xi*xi - 52.5*xi*eta*eta + 787.5*xi*(eta*eta*eta) - 787.5*xi*eta*eta*eta*eta + 157.5*(eta*eta) - 945.0*eta*eta*xi*xi + 787.5*(eta*eta)*(xi*xi*xi) + 393.75*(xi*xi)*(eta*eta*eta) - 621.25*eta*eta*eta + 866.25*(eta*eta*eta*eta) - 393.75*eta*eta*eta*eta*eta, 52.5*eta - 840.0*eta*xi*xi + 1575.0*eta*(xi*xi*xi) - 787.5*eta*xi*xi*xi*xi - 6.875*xi + 813.75*xi*(eta*eta) - 1260.0*xi*eta*eta*eta + 393.75*xi*(eta*eta*eta*eta) - 1.875 - 288.75*eta*eta - 236.25*eta*eta*xi*xi - 393.75*eta*eta*xi*xi*xi + 43.75*(xi*xi) + 787.5*(xi*xi)*(eta*eta*eta) + 472.5*(eta*eta*eta) - 8.75*xi*xi*xi - 236.25*eta*eta*eta*eta - 91.875*xi*xi*xi*xi + 65.625*(xi*xi*xi*xi*xi));
                case 13:
                  return sign * RealGradient(-95.15625*eta*xi + 18.7109375*eta - 44.296875*eta*xi*xi + 442.96875*eta*(xi*xi*xi) - 332.2265625*eta*xi*xi*xi*xi + 1230.46875*xi*(eta*eta) - 2628.28125*xi*eta*eta*eta + 1525.78125*xi*(eta*eta*eta*eta) - 214.375*eta*eta - 1594.6875*eta*eta*xi*xi + 442.96875*(eta*eta)*(xi*xi*xi) + 1993.359375*(xi*xi)*(eta*eta*eta) + 652.421875*(eta*eta*eta) - 759.609375*eta*eta*eta*eta + 303.515625*(eta*eta*eta*eta*eta), 759.0625*eta*xi - 89.53125*eta - 1811.25*eta*xi*xi + 1594.6875*eta*(xi*xi*xi) - 442.96875*eta*xi*xi*xi*xi - 5.4296875*xi - 2127.890625*xi*eta*eta + 1706.25*xi*(eta*eta*eta) - 303.515625*xi*eta*eta*eta*eta + 3.7890625 + 329.765625*(eta*eta) + 3765.234375*(eta*eta)*(xi*xi) - 1993.359375*eta*eta*xi*xi*xi - 175.546875*xi*xi - 1525.78125*xi*xi*eta*eta*eta - 425.46875*eta*eta*eta + 664.453125*(xi*xi*xi) + 182.109375*(eta*eta*eta*eta) - 819.4921875*xi*xi*xi*xi + 332.2265625*(xi*xi*xi*xi*xi));
                case 14:
                  return sign * RealGradient(700.0*eta*xi - 75.0*eta - 2100.0*eta*xi*xi + 2520.0*eta*(xi*xi*xi) - 1050.0*eta*xi*xi*xi*xi - 4200.0*xi*eta*eta + 7560.0*xi*(eta*eta*eta) - 4200.0*xi*eta*eta*eta*eta + 700.0*(eta*eta) + 7560.0*(eta*eta)*(xi*xi) - 4200.0*eta*eta*xi*xi*xi - 6300.0*xi*xi*eta*eta*eta - 2100.0*eta*eta*eta + 2520.0*(eta*eta*eta*eta) - 1050.0*eta*eta*eta*eta*eta, -2800.0*eta*xi + 300.0*eta + 8400.0*eta*(xi*xi) - 10080.0*eta*xi*xi*xi + 4200.0*eta*(xi*xi*xi*xi) + 375.0*xi + 6300.0*xi*(eta*eta) - 5040.0*xi*eta*eta*eta + 1050.0*xi*(eta*eta*eta*eta) - 25.0 - 1050.0*eta*eta - 11340.0*eta*eta*xi*xi + 6300.0*(eta*eta)*(xi*xi*xi) - 1750.0*xi*xi + 4200.0*(xi*xi)*(eta*eta*eta) + 1400.0*(eta*eta*eta) + 3500.0*(xi*xi*xi) - 630.0*eta*eta*eta*eta - 3150.0*xi*xi*xi*xi + 1050.0*(xi*xi*xi*xi*xi));
                case 15:
                  return RealGradient(-19600.0*eta*xi + 3500.0*eta + 35280.0*eta*(xi*xi) - 23520.0*eta*xi*xi*xi + 4200.0*eta*(xi*xi*xi*xi) + 94080.0*xi*(eta*eta) - 141120.0*xi*eta*eta*eta + 67200.0*xi*(eta*eta*eta*eta) - 24500.0*eta*eta - 105840.0*eta*eta*xi*xi + 33600.0*(eta*eta)*(xi*xi*xi) + 75600.0*(xi*xi)*(eta*eta*eta) + 58800.0*(eta*eta*eta) - 58800.0*eta*eta*eta*eta + 21000.0*(eta*eta*eta*eta*eta), 9800.0*eta*xi - 47040.0*eta*xi*xi + 70560.0*eta*(xi*xi*xi) - 33600.0*eta*xi*xi*xi*xi - 700.0*xi - 35280.0*xi*eta*eta + 47040.0*xi*(eta*eta*eta) - 21000.0*xi*eta*eta*eta*eta + 105840.0*(eta*eta)*(xi*xi) - 75600.0*eta*eta*xi*xi*xi + 4900.0*(xi*xi) - 67200.0*xi*xi*eta*eta*eta - 11760.0*xi*xi*xi + 11760.0*(xi*xi*xi*xi) - 4200.0*xi*xi*xi*xi*xi);
                case 16:
                  return RealGradient(9800.0*eta*xi - 700.0*eta - 35280.0*eta*xi*xi + 47040.0*eta*(xi*xi*xi) - 21000.0*eta*xi*xi*xi*xi - 47040.0*xi*eta*eta + 70560.0*xi*(eta*eta*eta) - 33600.0*xi*eta*eta*eta*eta + 4900.0*(eta*eta) + 105840.0*(eta*eta)*(xi*xi) - 67200.0*eta*eta*xi*xi*xi - 75600.0*xi*xi*eta*eta*eta - 11760.0*eta*eta*eta + 11760.0*(eta*eta*eta*eta) - 4200.0*eta*eta*eta*eta*eta, -19600.0*eta*xi + 94080.0*eta*(xi*xi) - 141120.0*eta*xi*xi*xi + 67200.0*eta*(xi*xi*xi*xi) + 3500.0*xi + 35280.0*xi*(eta*eta) - 23520.0*xi*eta*eta*eta + 4200.0*xi*(eta*eta*eta*eta) - 105840.0*eta*eta*xi*xi + 75600.0*(eta*eta)*(xi*xi*xi) - 24500.0*xi*xi + 33600.0*(xi*xi)*(eta*eta*eta) + 58800.0*(xi*xi*xi) - 58800.0*xi*xi*xi*xi + 21000.0*(xi*xi*xi*xi*xi));
                case 17:
                  return RealGradient(6440.0*eta*xi - 280.0*eta - 30240.0*eta*xi*xi + 43680.0*eta*(xi*xi*xi) - 16800.0*eta*xi*xi*xi*xi - 6720.0*xi*eta*eta + 280.0*(eta*eta) + 30240.0*(eta*eta)*(xi*xi) - 33600.0*eta*eta*xi*xi*xi, -1120.0*eta*xi + 13440.0*eta*(xi*xi) - 40320.0*eta*xi*xi*xi + 33600.0*eta*(xi*xi*xi*xi) + 560.0*xi - 7280.0*xi*xi + 26880.0*(xi*xi*xi) - 36960.0*xi*xi*xi*xi + 16800.0*(xi*xi*xi*xi*xi));
                case 18:
                  return RealGradient(-3640.0*eta*xi + 140.0*eta + 20160.0*eta*(xi*xi) - 36960.0*eta*xi*xi*xi + 21000.0*eta*(xi*xi*xi*xi) + 3360.0*xi*(eta*eta) - 140.0*eta*eta - 15120.0*eta*eta*xi*xi + 16800.0*(eta*eta)*(xi*xi*xi), 560.0*eta*xi - 6720.0*eta*xi*xi + 20160.0*eta*(xi*xi*xi) - 16800.0*eta*xi*xi*xi*xi - 700.0*xi + 9100.0*(xi*xi) - 33600.0*xi*xi*xi + 46200.0*(xi*xi*xi*xi) - 21000.0*xi*xi*xi*xi*xi);
                case 19:
                  return RealGradient(560.0*eta*xi - 700.0*eta - 6720.0*xi*eta*eta + 20160.0*xi*(eta*eta*eta) - 16800.0*xi*eta*eta*eta*eta + 9100.0*(eta*eta) - 33600.0*eta*eta*eta + 46200.0*(eta*eta*eta*eta) - 21000.0*eta*eta*eta*eta*eta, -3640.0*eta*xi + 3360.0*eta*(xi*xi) + 140.0*xi + 20160.0*xi*(eta*eta) - 36960.0*xi*eta*eta*eta + 21000.0*xi*(eta*eta*eta*eta) - 15120.0*eta*eta*xi*xi - 140.0*xi*xi + 16800.0*(xi*xi)*(eta*eta*eta));
                case 20:
                  return RealGradient(-1120.0*eta*xi + 560.0*eta + 13440.0*xi*(eta*eta) - 40320.0*xi*eta*eta*eta + 33600.0*xi*(eta*eta*eta*eta) - 7280.0*eta*eta + 26880.0*(eta*eta*eta) - 36960.0*eta*eta*eta*eta + 16800.0*(eta*eta*eta*eta*eta), 6440.0*eta*xi - 6720.0*eta*xi*xi - 280.0*xi - 30240.0*xi*eta*eta + 43680.0*xi*(eta*eta*eta) - 16800.0*xi*eta*eta*eta*eta + 30240.0*(eta*eta)*(xi*xi) + 280.0*(xi*xi) - 33600.0*xi*xi*eta*eta*eta);
                case 21:
                  return RealGradient(5973.33333333333*eta*xi - 420.0*eta - 17920.0*eta*xi*xi + 17422.2222222222*eta*(xi*xi*xi) - 4977.77777777778*eta*xi*xi*xi*xi - 29120.0*xi*eta*eta + 38826.6666666667*xi*(eta*eta*eta) - 15555.5555555556*xi*eta*eta*eta*eta + 2473.33333333333*(eta*eta) + 62720.0*(eta*eta)*(xi*xi) - 32355.5555555556*eta*eta*xi*xi*xi - 44800.0*xi*xi*eta*eta*eta - 4480.0*eta*eta*eta + 3204.44444444444*(eta*eta*eta*eta) - 777.777777777778*eta*eta*eta*eta*eta, -3453.33333333333*eta*xi + 28000.0*eta*(xi*xi) - 56746.6666666667*eta*xi*xi*xi + 32355.5555555556*eta*(xi*xi*xi*xi) + 420.0*xi + 6720.0*xi*(eta*eta) - 4355.55555555556*xi*eta*eta*eta + 777.777777777778*xi*(eta*eta*eta*eta) - 40880.0*eta*eta*xi*xi + 44800.0*(eta*eta)*(xi*xi*xi) - 4153.33333333333*xi*xi + 15555.5555555556*(xi*xi)*(eta*eta*eta) + 11946.6666666667*(xi*xi*xi) - 13191.1111111111*xi*xi*xi*xi + 4977.77777777778*(xi*xi*xi*xi*xi));
                case 22:
                  return RealGradient(-4293.33333333333*eta*xi + 280.0*eta + 14560.0*eta*(xi*xi) - 16924.4444444444*eta*xi*xi*xi + 6222.22222222222*eta*(xi*xi*xi*xi) + 22400.0*xi*(eta*eta) - 30613.3333333333*xi*eta*eta*eta + 12444.4444444444*xi*(eta*eta*eta*eta) - 1773.33333333333*eta*eta - 54880.0*eta*eta*xi*xi + 34844.4444444444*(eta*eta)*(xi*xi*xi) + 39200.0*(xi*xi)*(eta*eta*eta) + 3360.0*(eta*eta*eta) - 2488.88888888889*eta*eta*eta*eta + 622.222222222222*(eta*eta*eta*eta*eta), 4013.33333333333*eta*xi - 31360.0*eta*xi*xi + 61973.3333333333*eta*(xi*xi*xi) - 34844.4444444444*eta*xi*xi*xi*xi - 560.0*xi - 6720.0*xi*eta*eta + 3857.77777777778*xi*(eta*eta*eta) - 622.222222222222*xi*eta*eta*eta*eta + 38080.0*(eta*eta)*(xi*xi) - 39200.0*eta*eta*xi*xi*xi + 5413.33333333333*(xi*xi) - 12444.4444444444*xi*xi*eta*eta*eta - 15306.6666666667*xi*xi*xi + 16675.5555555556*(xi*xi*xi*xi) - 6222.22222222222*xi*xi*xi*xi*xi);
                case 23:
                  return RealGradient(4013.33333333333*eta*xi - 560.0*eta - 6720.0*eta*xi*xi + 3857.77777777778*eta*(xi*xi*xi) - 622.222222222222*eta*xi*xi*xi*xi - 31360.0*xi*eta*eta + 61973.3333333333*xi*(eta*eta*eta) - 34844.4444444444*xi*eta*eta*eta*eta + 5413.33333333333*(eta*eta) + 38080.0*(eta*eta)*(xi*xi) - 12444.4444444444*eta*eta*xi*xi*xi - 39200.0*xi*xi*eta*eta*eta - 15306.6666666667*eta*eta*eta + 16675.5555555556*(eta*eta*eta*eta) - 6222.22222222222*eta*eta*eta*eta*eta, -4293.33333333333*eta*xi + 22400.0*eta*(xi*xi) - 30613.3333333333*eta*xi*xi*xi + 12444.4444444444*eta*(xi*xi*xi*xi) + 280.0*xi + 14560.0*xi*(eta*eta) - 16924.4444444444*xi*eta*eta*eta + 6222.22222222222*xi*(eta*eta*eta*eta) - 54880.0*eta*eta*xi*xi + 39200.0*(eta*eta)*(xi*xi*xi) - 1773.33333333333*xi*xi + 34844.4444444444*(xi*xi)*(eta*eta*eta) + 3360.0*(xi*xi*xi) - 2488.88888888889*xi*xi*xi*xi + 622.222222222222*(xi*xi*xi*xi*xi));
                case 24:
                  return RealGradient(-3453.33333333333*eta*xi + 420.0*eta + 6720.0*eta*(xi*xi) - 4355.55555555556*eta*xi*xi*xi + 777.777777777778*eta*(xi*xi*xi*xi) + 28000.0*xi*(eta*eta) - 56746.6666666667*xi*eta*eta*eta + 32355.5555555556*xi*(eta*eta*eta*eta) - 4153.33333333333*eta*eta - 40880.0*eta*eta*xi*xi + 15555.5555555556*(eta*eta)*(xi*xi*xi) + 44800.0*(xi*xi)*(eta*eta*eta) + 11946.6666666667*(eta*eta*eta) - 13191.1111111111*eta*eta*eta*eta + 4977.77777777778*(eta*eta*eta*eta*eta), 5973.33333333333*eta*xi - 29120.0*eta*xi*xi + 38826.6666666667*eta*(xi*xi*xi) - 15555.5555555556*eta*xi*xi*xi*xi - 420.0*xi - 17920.0*xi*eta*eta + 17422.2222222222*xi*(eta*eta*eta) - 4977.77777777778*xi*eta*eta*eta*eta + 62720.0*(eta*eta)*(xi*xi) - 44800.0*eta*eta*xi*xi*xi + 2473.33333333333*(xi*xi) - 32355.5555555556*xi*xi*eta*eta*eta - 4480.0*xi*xi*xi + 3204.44444444444*(xi*xi*xi*xi) - 777.777777777778*xi*xi*xi*xi*xi);
                case 25:
                  return RealGradient(-1431.11111111111*eta*xi - 77.7777777777778*eta + 5600.0*eta*(xi*xi) - 5475.55555555556*eta*xi*xi*xi + 1244.44444444444*eta*(xi*xi*xi*xi) - 5226.66666666667*xi*eta*eta + 17173.3333333333*xi*(eta*eta*eta) - 10577.7777777778*xi*eta*eta*eta*eta + 2877.77777777778*(eta*eta) - 1120.0*eta*eta*xi*xi + 2488.88888888889*(eta*eta)*(xi*xi*xi) - 5600.0*xi*xi*eta*eta*eta - 9333.33333333333*eta*eta*eta + 10422.2222222222*(eta*eta*eta*eta) - 3888.88888888889*eta*eta*eta*eta*eta, -1151.11111111111*eta*xi + 2613.33333333333*eta*(xi*xi) + 746.666666666667*eta*(xi*xi*xi) - 2488.88888888889*eta*xi*xi*xi*xi + 15.5555555555556*xi + 5600.0*xi*(eta*eta) - 8337.77777777778*xi*eta*eta*eta + 3888.88888888889*xi*(eta*eta*eta*eta) - 12880.0*eta*eta*xi*xi + 5600.0*(eta*eta)*(xi*xi*xi) + 357.777777777778*(xi*xi) + 10577.7777777778*(xi*xi)*(eta*eta*eta) - 1866.66666666667*xi*xi*xi + 2737.77777777778*(xi*xi*xi*xi) - 1244.44444444444*xi*xi*xi*xi*xi);
                case 26:
                  return RealGradient(-1057.77777777778*eta*xi + 155.555555555556*eta - 1120.0*eta*xi*xi + 7964.44444444444*eta*(xi*xi*xi) - 6222.22222222222*eta*xi*xi*xi*xi + 14933.3333333333*xi*(eta*eta) - 27626.6666666667*xi*eta*eta*eta + 13688.8888888889*xi*(eta*eta*eta*eta) - 1648.88888888889*eta*eta - 25760.0*eta*eta*xi*xi + 9955.55555555555*(eta*eta)*(xi*xi*xi) + 28000.0*(xi*xi)*(eta*eta*eta) + 4106.66666666667*(eta*eta*eta) - 3857.77777777778*eta*eta*eta*eta + 1244.44444444444*(eta*eta*eta*eta*eta), 5755.55555555556*eta*xi - 23146.6666666667*eta*xi*xi + 27626.6666666667*eta*(xi*xi*xi) - 9955.55555555555*eta*xi*xi*xi*xi + 62.2222222222222*xi - 13440.0*xi*eta*eta + 8835.55555555555*xi*(eta*eta*eta) - 1244.44444444444*xi*eta*eta*eta*eta + 40320.0*(eta*eta)*(xi*xi) - 28000.0*eta*eta*xi*xi*xi - 2675.55555555556*xi*xi - 13688.8888888889*xi*xi*eta*eta*eta + 10826.6666666667*(xi*xi*xi) - 14435.5555555556*xi*xi*xi*xi + 6222.22222222222*(xi*xi*xi*xi*xi));
                case 27:
                  return RealGradient(311.111111111111*eta*xi + 77.7777777777778*eta - 560.0*eta*xi*xi - 124.444444444444*eta*xi*xi*xi + 155.555555555556*eta*(xi*xi*xi*xi) - 1493.33333333333*xi*eta*eta - 3733.33333333333*xi*eta*eta*eta + 4977.77777777778*xi*(eta*eta*eta*eta) - 1477.77777777778*eta*eta + 6160.0*(eta*eta)*(xi*xi) - 2488.88888888889*eta*eta*xi*xi*xi - 2800.0*xi*xi*eta*eta*eta + 6533.33333333333*(eta*eta*eta) - 9022.22222222222*eta*eta*eta*eta + 3888.88888888889*(eta*eta*eta*eta*eta), 591.111111111111*eta*xi + 746.666666666667*eta*(xi*xi) - 4106.66666666667*eta*xi*xi*xi + 2488.88888888889*eta*(xi*xi*xi*xi) - 15.5555555555556*xi - 3920.0*xi*eta*eta + 7217.77777777778*xi*(eta*eta*eta) - 3888.88888888889*xi*eta*eta*eta*eta + 2800.0*(eta*eta)*(xi*xi) + 2800.0*(eta*eta)*(xi*xi*xi) - 77.7777777777778*xi*xi - 4977.77777777778*xi*xi*eta*eta*eta + 186.666666666667*(xi*xi*xi) + 62.2222222222222*(xi*xi*xi*xi) - 155.555555555556*xi*xi*xi*xi*xi);
                case 28:
                  return RealGradient(31.1111111111111*eta*xi - 108.888888888889*eta + 1680.0*eta*(xi*xi) - 1244.44444444444*eta*xi*xi*xi - 777.777777777778*eta*xi*xi*xi*xi - 3733.33333333333*xi*eta*eta + 17546.6666666667*xi*(eta*eta*eta) - 13688.8888888889*xi*eta*eta*eta*eta + 1508.88888888889*(eta*eta) - 9520.0*eta*eta*xi*xi + 12444.4444444444*(eta*eta)*(xi*xi*xi) - 2800.0*xi*xi*eta*eta*eta - 5413.33333333333*eta*eta*eta + 6657.77777777778*(eta*eta*eta*eta) - 2644.44444444444*eta*eta*eta*eta*eta, -995.555555555556*eta*xi - 5973.33333333333*eta*xi*xi + 19413.3333333333*eta*(xi*xi*xi) - 12444.4444444444*eta*xi*xi*xi*xi - 15.5555555555556*xi + 9520.0*xi*(eta*eta) - 11075.5555555556*xi*eta*eta*eta + 2644.44444444444*xi*(eta*eta*eta*eta) - 15120.0*eta*eta*xi*xi + 2800.0*(eta*eta)*(xi*xi*xi) + 482.222222222222*(xi*xi) + 13688.8888888889*(xi*xi)*(eta*eta*eta) - 560.0*xi*xi*xi - 684.444444444444*xi*xi*xi*xi + 777.777777777778*(xi*xi*xi*xi*xi));
                case 29:
                  return RealGradient(5755.55555555556*eta*xi + 62.2222222222222*eta - 13440.0*eta*xi*xi + 8835.55555555555*eta*(xi*xi*xi) - 1244.44444444444*eta*xi*xi*xi*xi - 23146.6666666667*xi*eta*eta + 27626.6666666667*xi*(eta*eta*eta) - 9955.55555555555*xi*eta*eta*eta*eta - 2675.55555555556*eta*eta + 40320.0*(eta*eta)*(xi*xi) - 13688.8888888889*eta*eta*xi*xi*xi - 28000.0*xi*xi*eta*eta*eta + 10826.6666666667*(eta*eta*eta) - 14435.5555555556*eta*eta*eta*eta + 6222.22222222222*(eta*eta*eta*eta*eta), -1057.77777777778*eta*xi + 14933.3333333333*eta*(xi*xi) - 27626.6666666667*eta*xi*xi*xi + 13688.8888888889*eta*(xi*xi*xi*xi) + 155.555555555556*xi - 1120.0*xi*eta*eta + 7964.44444444444*xi*(eta*eta*eta) - 6222.22222222222*xi*eta*eta*eta*eta - 25760.0*eta*eta*xi*xi + 28000.0*(eta*eta)*(xi*xi*xi) - 1648.88888888889*xi*xi + 9955.55555555555*(xi*xi)*(eta*eta*eta) + 4106.66666666667*(xi*xi*xi) - 3857.77777777778*xi*xi*xi*xi + 1244.44444444444*(xi*xi*xi*xi*xi));
                case 30:
                  return RealGradient(-1151.11111111111*eta*xi + 15.5555555555556*eta + 5600.0*eta*(xi*xi) - 8337.77777777778*eta*xi*xi*xi + 3888.88888888889*eta*(xi*xi*xi*xi) + 2613.33333333333*xi*(eta*eta) + 746.666666666667*xi*(eta*eta*eta) - 2488.88888888889*xi*eta*eta*eta*eta + 357.777777777778*(eta*eta) - 12880.0*eta*eta*xi*xi + 10577.7777777778*(eta*eta)*(xi*xi*xi) + 5600.0*(xi*xi)*(eta*eta*eta) - 1866.66666666667*eta*eta*eta + 2737.77777777778*(eta*eta*eta*eta) - 1244.44444444444*eta*eta*eta*eta*eta, -1431.11111111111*eta*xi - 5226.66666666667*eta*xi*xi + 17173.3333333333*eta*(xi*xi*xi) - 10577.7777777778*eta*xi*xi*xi*xi - 77.7777777777778*xi + 5600.0*xi*(eta*eta) - 5475.55555555556*xi*eta*eta*eta + 1244.44444444444*xi*(eta*eta*eta*eta) - 1120.0*eta*eta*xi*xi - 5600.0*eta*eta*xi*xi*xi + 2877.77777777778*(xi*xi) + 2488.88888888889*(xi*xi)*(eta*eta*eta) - 9333.33333333333*xi*xi*xi + 10422.2222222222*(xi*xi*xi*xi) - 3888.88888888889*xi*xi*xi*xi*xi);
                case 31:
                  return RealGradient(-995.555555555556*eta*xi - 15.5555555555556*eta + 9520.0*eta*(xi*xi) - 11075.5555555556*eta*xi*xi*xi + 2644.44444444444*eta*(xi*xi*xi*xi) - 5973.33333333333*xi*eta*eta + 19413.3333333333*xi*(eta*eta*eta) - 12444.4444444444*xi*eta*eta*eta*eta + 482.222222222222*(eta*eta) - 15120.0*eta*eta*xi*xi + 13688.8888888889*(eta*eta)*(xi*xi*xi) + 2800.0*(xi*xi)*(eta*eta*eta) - 560.0*eta*eta*eta - 684.444444444444*eta*eta*eta*eta + 777.777777777778*(eta*eta*eta*eta*eta), 31.1111111111111*eta*xi - 3733.33333333333*eta*xi*xi + 17546.6666666667*eta*(xi*xi*xi) - 13688.8888888889*eta*xi*xi*xi*xi - 108.888888888889*xi + 1680.0*xi*(eta*eta) - 1244.44444444444*xi*eta*eta*eta - 777.777777777778*xi*eta*eta*eta*eta - 9520.0*eta*eta*xi*xi - 2800.0*eta*eta*xi*xi*xi + 1508.88888888889*(xi*xi) + 12444.4444444444*(xi*xi)*(eta*eta*eta) - 5413.33333333333*xi*xi*xi + 6657.77777777778*(xi*xi*xi*xi) - 2644.44444444444*xi*xi*xi*xi*xi);
                case 32:
                  return RealGradient(591.111111111111*eta*xi - 15.5555555555556*eta - 3920.0*eta*xi*xi + 7217.77777777778*eta*(xi*xi*xi) - 3888.88888888889*eta*xi*xi*xi*xi + 746.666666666667*xi*(eta*eta) - 4106.66666666667*xi*eta*eta*eta + 2488.88888888889*xi*(eta*eta*eta*eta) - 77.7777777777778*eta*eta + 2800.0*(eta*eta)*(xi*xi) - 4977.77777777778*eta*eta*xi*xi*xi + 2800.0*(xi*xi)*(eta*eta*eta) + 186.666666666667*(eta*eta*eta) + 62.2222222222222*(eta*eta*eta*eta) - 155.555555555556*eta*eta*eta*eta*eta, 311.111111111111*eta*xi - 1493.33333333333*eta*xi*xi - 3733.33333333333*eta*xi*xi*xi + 4977.77777777778*eta*(xi*xi*xi*xi) + 77.7777777777778*xi - 560.0*xi*eta*eta - 124.444444444444*xi*eta*eta*eta + 155.555555555556*xi*(eta*eta*eta*eta) + 6160.0*(eta*eta)*(xi*xi) - 2800.0*eta*eta*xi*xi*xi - 1477.77777777778*xi*xi - 2488.88888888889*xi*xi*eta*eta*eta + 6533.33333333333*(xi*xi*xi) - 9022.22222222222*xi*xi*xi*xi + 3888.88888888889*(xi*xi*xi*xi*xi));
                case 33:
                  return RealGradient(-715.555555555556*eta*xi + 31.1111111111111*eta + 1680.0*eta*(xi*xi) - 1493.33333333333*eta*xi*xi*xi + 466.666666666667*eta*(xi*xi*xi*xi) + 6346.66666666667*xi*(eta*eta) - 11200.0*xi*eta*eta*eta + 5600.0*xi*(eta*eta*eta*eta) - 311.111111111111*eta*eta - 9520.0*eta*eta*xi*xi + 3733.33333333333*(eta*eta)*(xi*xi*xi) + 8400.0*(xi*xi)*(eta*eta*eta) + 560.0*(eta*eta*eta) - 280.0*eta*eta*eta*eta, 684.444444444444*eta*xi - 4853.33333333333*eta*xi*xi + 7840.0*eta*(xi*xi*xi) - 3733.33333333333*eta*xi*xi*xi*xi - 62.2222222222222*xi - 1680.0*xi*eta*eta + 1120.0*xi*(eta*eta*eta) + 10080.0*(eta*eta)*(xi*xi) - 8400.0*eta*eta*xi*xi*xi + 528.888888888889*(xi*xi) - 5600.0*xi*xi*eta*eta*eta - 1306.66666666667*xi*xi*xi + 1306.66666666667*(xi*xi*xi*xi) - 466.666666666667*xi*xi*xi*xi*xi);
                case 34:
                  return RealGradient(684.444444444444*eta*xi - 62.2222222222222*eta - 1680.0*eta*xi*xi + 1120.0*eta*(xi*xi*xi) - 4853.33333333333*xi*eta*eta + 7840.0*xi*(eta*eta*eta) - 3733.33333333333*xi*eta*eta*eta*eta + 528.888888888889*(eta*eta) + 10080.0*(eta*eta)*(xi*xi) - 5600.0*eta*eta*xi*xi*xi - 8400.0*xi*xi*eta*eta*eta - 1306.66666666667*eta*eta*eta + 1306.66666666667*(eta*eta*eta*eta) - 466.666666666667*eta*eta*eta*eta*eta, -715.555555555556*eta*xi + 6346.66666666667*eta*(xi*xi) - 11200.0*eta*xi*xi*xi + 5600.0*eta*(xi*xi*xi*xi) + 31.1111111111111*xi + 1680.0*xi*(eta*eta) - 1493.33333333333*xi*eta*eta*eta + 466.666666666667*xi*(eta*eta*eta*eta) - 9520.0*eta*eta*xi*xi + 8400.0*(eta*eta)*(xi*xi*xi) - 311.111111111111*xi*xi + 3733.33333333333*(xi*xi)*(eta*eta*eta) + 560.0*(xi*xi*xi) - 280.0*xi*xi*xi*xi);
                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          } // end switch (type)
      } // end case FIFTH

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

      // quartic Nedelec (first kind) shape function first derivatives
    case FOURTH:
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
                        return sign * RealGradient(240.0*eta + 60.0*xi - 480.0*(0.5*eta + 0.5)*(2*xi + 2) + 420.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 450.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 300.0*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 270.0 - 450.0*(eta + 1)*(eta + 1) - 393.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 26.25*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 57.421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 300.0*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 1:
                        return sign * RealGradient(-84.4444444444444*eta - 23.8888888888889*xi + 191.111111111111*(0.5*eta + 0.5)*(2*xi + 2) - 171.111111111111*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 179.166666666667*(2*xi + 2)*(eta + 1)*(eta + 1) + 119.444444444444*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 26.1284722222222*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 97.7777777777778 + 158.333333333333*((eta + 1)*(eta + 1)) + 160.416666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 10.6944444444444*((xi + 1)*(xi + 1)) - 106.944444444444*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 23.3940972222222*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 105.555555555556*(eta + 1)*(eta + 1)*(eta + 1) + 23.0902777777778*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 2:
                        return sign * RealGradient(44.4444444444444*eta + 18.8888888888889*xi - 151.111111111111*(0.5*eta + 0.5)*(2*xi + 2) + 171.111111111111*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 141.666666666667*(2*xi + 2)*((eta + 1)*(eta + 1)) - 94.4444444444444*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 20.6597222222222*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 57.7777777777778 - 83.3333333333333*(eta + 1)*(eta + 1) - 160.416666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 10.6944444444444*(xi + 1)*(xi + 1) + 106.944444444444*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 23.3940972222222*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 55.5555555555556*((eta + 1)*(eta + 1)*(eta + 1)) - 12.1527777777778*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 3:
                        return sign * RealGradient(-120.0*eta - 45.0*xi + 360.0*(0.5*eta + 0.5)*(2*xi + 2) - 420.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 337.5*(2*xi + 2)*(eta + 1)*(eta + 1) + 225.0*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 150.0 + 225.0*((eta + 1)*(eta + 1)) + 393.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 26.25*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 57.421875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 150.0*(eta + 1)*(eta + 1)*(eta + 1) + 32.8125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 4:
                        return sign * RealGradient(0, 60.0*eta + 120.0*xi - 450.0*(0.5*eta + 0.5)*(2*xi + 2) + 1350.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 525.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 225.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 65.625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 164.0 - 60.0*(eta + 1)*(eta + 1) - 675.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 180.0*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 70.0*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 5:
                        return sign * RealGradient(0, -21.1111111111111*eta - 18.8888888888889*xi + 158.333333333333*(0.5*eta + 0.5)*(2*xi + 2) - 475.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 184.722222222222*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 89.5833333333333*(2*xi + 2)*(eta + 1)*(eta + 1) + 26.7361111111111*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 37.4814814814815 + 23.8888888888889*((eta + 1)*(eta + 1)) + 268.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 104.513888888889*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 28.3333333333333*((xi + 1)*(xi + 1)) - 80.2083333333333*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 7.12962962962963*(eta + 1)*(eta + 1)*(eta + 1) + 31.1921296296296*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 11.0185185185185*(xi + 1)*(xi + 1)*(xi + 1));
                      case 6:
                        return sign * RealGradient(0, 11.1111111111111*eta + 8.88888888888889*xi - 83.3333333333333*(0.5*eta + 0.5)*(2*xi + 2) + 250.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 97.2222222222222*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 70.8333333333333*(2*xi + 2)*((eta + 1)*(eta + 1)) - 26.7361111111111*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 18.8148148148148 - 18.8888888888889*(eta + 1)*(eta + 1) - 212.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 82.6388888888889*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 13.3333333333333*(xi + 1)*(xi + 1) + 80.2083333333333*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 7.12962962962963*((eta + 1)*(eta + 1)*(eta + 1)) - 31.1921296296296*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5.18518518518519*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 7:
                        return sign * RealGradient(0, -30.0*eta - 30.0*xi + 225.0*(0.5*eta + 0.5)*(2*xi + 2) - 675.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 262.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 168.75*(2*xi + 2)*(eta + 1)*(eta + 1) + 65.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 56.0 + 45.0*((eta + 1)*(eta + 1)) + 506.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 45.0*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 17.5*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 8:
                        return sign * RealGradient(30.0*eta - 90.0*(0.5*eta + 0.5)*(2*xi + 2) + 105.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 168.75*(2*xi + 2)*((eta + 1)*(eta + 1)) - 168.75*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 49.21875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 30.0 - 112.5*(eta + 1)*(eta + 1) - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 57.421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 112.5*((eta + 1)*(eta + 1)*(eta + 1)) - 32.8125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 9:
                        return sign * RealGradient(-11.1111111111111*eta + 37.7777777777778*(0.5*eta + 0.5)*(2*xi + 2) - 42.7777777777778*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 70.8333333333333*(2*xi + 2)*(eta + 1)*(eta + 1) + 70.8333333333333*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 20.6597222222222*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 11.1111111111111 + 41.6666666666667*((eta + 1)*(eta + 1)) + 80.2083333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 80.2083333333333*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 23.3940972222222*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 41.6666666666667*(eta + 1)*(eta + 1)*(eta + 1) + 12.1527777777778*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 10:
                        return sign * RealGradient(21.1111111111111*eta - 47.7777777777778*(0.5*eta + 0.5)*(2*xi + 2) + 42.7777777777778*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 89.5833333333333*(2*xi + 2)*((eta + 1)*(eta + 1)) - 89.5833333333333*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 26.1284722222222*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 21.1111111111111 - 79.1666666666667*(eta + 1)*(eta + 1) - 80.2083333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 80.2083333333333*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 23.3940972222222*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 79.1666666666667*((eta + 1)*(eta + 1)*(eta + 1)) - 23.0902777777778*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 11:
                        return sign * RealGradient(-60.0*eta + 120.0*(0.5*eta + 0.5)*(2*xi + 2) - 105.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 225.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 225.0*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 60.0 + 225.0*((eta + 1)*(eta + 1)) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 57.421875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 225.0*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 12:
                        return sign * RealGradient(0, 120.0*eta + 60.0*xi - 450.0*(0.5*eta + 0.5)*(2*xi + 2) + 900.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 262.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 337.5*(2*xi + 2)*((eta + 1)*(eta + 1)) - 131.25*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 164.0 - 180.0*(eta + 1)*(eta + 1) - 675.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 60.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 70.0*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 13:
                        return sign * RealGradient(0, -44.4444444444444*eta - 17.7777777777778*xi + 166.666666666667*(0.5*eta + 0.5)*(2*xi + 2) - 333.333333333333*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 97.2222222222222*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 141.666666666667*(2*xi + 2)*(eta + 1)*(eta + 1) + 53.4722222222222*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 57.4814814814815 + 75.5555555555556*((eta + 1)*(eta + 1)) + 283.333333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 82.6388888888889*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 17.7777777777778*((xi + 1)*(xi + 1)) - 106.944444444444*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 28.5185185185185*(eta + 1)*(eta + 1)*(eta + 1) + 31.1921296296296*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5.18518518518519*(xi + 1)*(xi + 1)*(xi + 1));
                      case 14:
                        return sign * RealGradient(0, 84.4444444444444*eta + 37.7777777777778*xi - 316.666666666667*(0.5*eta + 0.5)*(2*xi + 2) + 633.333333333333*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 184.722222222222*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 179.166666666667*(2*xi + 2)*((eta + 1)*(eta + 1)) - 53.4722222222222*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 112.148148148148 - 95.5555555555556*(eta + 1)*(eta + 1) - 358.333333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 104.513888888889*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 37.7777777777778*(xi + 1)*(xi + 1) + 106.944444444444*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 28.5185185185185*((eta + 1)*(eta + 1)*(eta + 1)) - 31.1921296296296*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 11.0185185185185*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 15:
                        return sign * RealGradient(0, -240.0*eta - 240.0*xi + 900.0*(0.5*eta + 0.5)*(2*xi + 2) - 1800.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 525.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 450.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 131.25*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 416.0 + 240.0*((eta + 1)*(eta + 1)) + 900.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 240.0*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 70.0*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 70.0*(xi + 1)*(xi + 1)*(xi + 1));
                      case 16:
                        return RealGradient(0, -195.0*eta - 232.0*xi + 870.0*(0.5*eta + 0.5)*(2*xi + 2) - 1800.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 525.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 435.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 126.875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 375.0 + 195.0*((eta + 1)*(eta + 1)) + 900.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 240.0*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 56.875*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 70.0*(xi + 1)*(xi + 1)*(xi + 1));
                      case 17:
                        return RealGradient(0, -45.0*eta - 112.0*xi + 420.0*(0.5*eta + 0.5)*(2*xi + 2) - 1350.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 525.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 210.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 61.25*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 145.0 + 45.0*((eta + 1)*(eta + 1)) + 675.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 180.0*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 13.125*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 70.0*(xi + 1)*(xi + 1)*(xi + 1));
                      case 18:
                        return RealGradient(0, -60.0*eta - 16.0*xi + 60.0*(0.5*eta + 0.5)*(2*xi + 2) - 30.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 8.75*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 60.0 + 60.0*((eta + 1)*(eta + 1)) - 17.5*(eta + 1)*(eta + 1)*(eta + 1));
                      case 19:
                        return RealGradient(195.0*eta - 390.0*(0.5*eta + 0.5)*(2*xi + 2) + 341.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 435.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 300.0*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 195.0 - 435.0*(eta + 1)*(eta + 1) - 380.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 57.421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 300.0*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 20:
                        return RealGradient(45.0*eta - 90.0*(0.5*eta + 0.5)*(2*xi + 2) + 78.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 210.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 225.0*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 45.0 - 210.0*(eta + 1)*(eta + 1) - 183.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 57.421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 225.0*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 21:
                        return RealGradient(60.0*eta - 120.0*(0.5*eta + 0.5)*(2*xi + 2) + 105.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 30.0*(2*xi + 2)*((eta + 1)*(eta + 1)) + 60.0 - 30.0*(eta + 1)*(eta + 1) - 26.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1), 0);
                      case 22:
                        return RealGradient(-97.5*eta + 292.5*(0.5*eta + 0.5)*(2*xi + 2) - 341.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 326.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 225.0*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 97.5 + 217.5*((eta + 1)*(eta + 1)) + 380.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 57.421875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 150.0*(eta + 1)*(eta + 1)*(eta + 1) + 32.8125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 23:
                        return RealGradient(-22.5*eta + 67.5*(0.5*eta + 0.5)*(2*xi + 2) - 78.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 157.5*(2*xi + 2)*(eta + 1)*(eta + 1) + 168.75*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 22.5 + 105.0*((eta + 1)*(eta + 1)) + 183.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 57.421875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 112.5*(eta + 1)*(eta + 1)*(eta + 1) + 32.8125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 24:
                        return RealGradient(-30.0*eta + 90.0*(0.5*eta + 0.5)*(2*xi + 2) - 105.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 22.5*(2*xi + 2)*(eta + 1)*(eta + 1) - 30.0 + 15.0*((eta + 1)*(eta + 1)) + 26.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)), 0);
                      case 25:
                        return RealGradient(0, 97.5*eta + 58.0*xi - 435.0*(0.5*eta + 0.5)*(2*xi + 2) + 900.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 262.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 326.25*(2*xi + 2)*((eta + 1)*(eta + 1)) - 126.875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 142.5 - 146.25*(eta + 1)*(eta + 1) - 675.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 60.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 56.875*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 26:
                        return RealGradient(0, 22.5*eta + 28.0*xi - 210.0*(0.5*eta + 0.5)*(2*xi + 2) + 675.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 262.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 157.5*(2*xi + 2)*((eta + 1)*(eta + 1)) - 61.25*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 47.5 - 33.75*(eta + 1)*(eta + 1) - 506.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 45.0*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 13.125*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 27:
                        return RealGradient(0, 30.0*eta + 4.0*xi - 30.0*(0.5*eta + 0.5)*(2*xi + 2) + 22.5*(2*xi + 2)*((eta + 1)*(eta + 1)) - 8.75*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 30.0 - 45.0*(eta + 1)*(eta + 1) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 28:
                        return RealGradient(-9.0*eta - 9.0 + 21.375*((eta + 1)*(eta + 1)) - 15.0*(eta + 1)*(eta + 1)*(eta + 1) + 3.28125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 29:
                        return RealGradient(9.0*eta + 9.0 - 21.375*(eta + 1)*(eta + 1) + 15.0*((eta + 1)*(eta + 1)*(eta + 1)) - 3.28125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 30:
                        return RealGradient(0, -9.0*eta - 57.0*xi + 42.75*(0.5*eta + 0.5)*(2*xi + 2) - 90.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 26.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 54.0 + 60.0*((xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 31:
                        return RealGradient(0, 9.0*eta + 28.5*xi - 42.75*(0.5*eta + 0.5)*(2*xi + 2) + 90.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 26.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 31.5 - 30.0*(xi + 1)*(xi + 1) + 8.75*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 32:
                        return RealGradient(0, -1.5*eta - 27.0*xi + 20.25*(0.5*eta + 0.5)*(2*xi + 2) - 67.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 26.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 26.5 + 45.0*((xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 33:
                        return RealGradient(0, 1.5*eta + 13.5*xi - 20.25*(0.5*eta + 0.5)*(2*xi + 2) + 67.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 26.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 14.0 - 22.5*(xi + 1)*(xi + 1) + 8.75*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 34:
                        return RealGradient(-1.5*eta - 1.5 + 10.125*((eta + 1)*(eta + 1)) - 11.25*(eta + 1)*(eta + 1)*(eta + 1) + 3.28125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 35:
                        return RealGradient(1.5*eta + 1.5 - 10.125*(eta + 1)*(eta + 1) + 11.25*((eta + 1)*(eta + 1)*(eta + 1)) - 3.28125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 36:
                        return RealGradient(0, -4.5*eta - 6.0*xi + 4.5*(0.5*eta + 0.5)*(2*xi + 2) - 4.5);
                      case 37:
                        return RealGradient(4.5*eta + 4.5 - 2.25*(eta + 1)*(eta + 1), 0);
                      case 38:
                        return RealGradient(-4.5*eta - 4.5 + 2.25*((eta + 1)*(eta + 1)), 0);
                      case 39:
                        return RealGradient(0, 4.5*eta + 3.0*xi - 4.5*(0.5*eta + 0.5)*(2*xi + 2) + 4.5);
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
                        return sign * RealGradient(240.0*eta + 240.0*xi - 900.0*(2*eta + 2)*(0.5*xi + 0.5) + 450.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 131.25*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1800.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 525.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 416.0 - 240.0*(eta + 1)*(eta + 1) - 900.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 240.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 70.0*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 70.0*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 1:
                        return sign * RealGradient(-37.7777777777778*eta - 84.4444444444444*xi + 316.666666666667*(2*eta + 2)*(0.5*xi + 0.5) - 179.166666666667*(2*eta + 2)*(xi + 1)*(xi + 1) + 53.4722222222222*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 633.333333333333*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 184.722222222222*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 112.148148148148 + 37.7777777777778*((eta + 1)*(eta + 1)) + 358.333333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 106.944444444444*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 95.5555555555556*((xi + 1)*(xi + 1)) - 104.513888888889*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 11.0185185185185*(eta + 1)*(eta + 1)*(eta + 1) + 31.1921296296296*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 28.5185185185185*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 2:
                        return sign * RealGradient(17.7777777777778*eta + 44.4444444444444*xi - 166.666666666667*(2*eta + 2)*(0.5*xi + 0.5) + 141.666666666667*(2*eta + 2)*((xi + 1)*(xi + 1)) - 53.4722222222222*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 333.333333333333*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 97.2222222222222*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 57.4814814814815 - 17.7777777777778*(eta + 1)*(eta + 1) - 283.333333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 106.944444444444*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 75.5555555555556*(xi + 1)*(xi + 1) + 82.6388888888889*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 5.18518518518519*((eta + 1)*(eta + 1)*(eta + 1)) - 31.1921296296296*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 28.5185185185185*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 3:
                        return sign * RealGradient(-60.0*eta - 120.0*xi + 450.0*(2*eta + 2)*(0.5*xi + 0.5) - 337.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 131.25*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 900.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 262.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 164.0 + 60.0*((eta + 1)*(eta + 1)) + 675.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 180.0*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 17.5*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 70.0*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 4:
                        return sign * RealGradient(0, 60.0*xi - 120.0*(2*eta + 2)*(0.5*xi + 0.5) + 225.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 225.0*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 105.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 60.0 - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 57.421875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 225.0*(xi + 1)*(xi + 1) + 225.0*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 5:
                        return sign * RealGradient(0, -21.1111111111111*xi + 47.7777777777778*(2*eta + 2)*(0.5*xi + 0.5) - 89.5833333333333*(2*eta + 2)*(xi + 1)*(xi + 1) + 89.5833333333333*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 26.1284722222222*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 42.7777777777778*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 21.1111111111111 + 80.2083333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 80.2083333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 23.3940972222222*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 79.1666666666667*((xi + 1)*(xi + 1)) - 79.1666666666667*(xi + 1)*(xi + 1)*(xi + 1) + 23.0902777777778*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 6:
                        return sign * RealGradient(0, 11.1111111111111*xi - 37.7777777777778*(2*eta + 2)*(0.5*xi + 0.5) + 70.8333333333333*(2*eta + 2)*((xi + 1)*(xi + 1)) - 70.8333333333333*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 20.6597222222222*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 42.7777777777778*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 11.1111111111111 - 80.2083333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 80.2083333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 23.3940972222222*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 41.6666666666667*(xi + 1)*(xi + 1) + 41.6666666666667*((xi + 1)*(xi + 1)*(xi + 1)) - 12.1527777777778*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 7:
                        return sign * RealGradient(0, -30.0*xi + 90.0*(2*eta + 2)*(0.5*xi + 0.5) - 168.75*(2*eta + 2)*(xi + 1)*(xi + 1) + 168.75*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 49.21875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 105.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 30.0 + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 57.421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 112.5*((xi + 1)*(xi + 1)) - 112.5*(xi + 1)*(xi + 1)*(xi + 1) + 32.8125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 8:
                        return sign * RealGradient(30.0*eta + 30.0*xi - 225.0*(2*eta + 2)*(0.5*xi + 0.5) + 168.75*(2*eta + 2)*((xi + 1)*(xi + 1)) - 65.625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 675.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 262.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 56.0 - 45.0*(eta + 1)*(eta + 1) - 506.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 45.0*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 9:
                        return sign * RealGradient(-8.88888888888889*eta - 11.1111111111111*xi + 83.3333333333333*(2*eta + 2)*(0.5*xi + 0.5) - 70.8333333333333*(2*eta + 2)*(xi + 1)*(xi + 1) + 26.7361111111111*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 250.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 97.2222222222222*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 18.8148148148148 + 13.3333333333333*((eta + 1)*(eta + 1)) + 212.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 80.2083333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 18.8888888888889*((xi + 1)*(xi + 1)) - 82.6388888888889*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 5.18518518518519*(eta + 1)*(eta + 1)*(eta + 1) + 31.1921296296296*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 7.12962962962963*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 10:
                        return sign * RealGradient(18.8888888888889*eta + 21.1111111111111*xi - 158.333333333333*(2*eta + 2)*(0.5*xi + 0.5) + 89.5833333333333*(2*eta + 2)*((xi + 1)*(xi + 1)) - 26.7361111111111*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 475.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 184.722222222222*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 37.4814814814815 - 28.3333333333333*(eta + 1)*(eta + 1) - 268.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 80.2083333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 23.8888888888889*(xi + 1)*(xi + 1) + 104.513888888889*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 11.0185185185185*((eta + 1)*(eta + 1)*(eta + 1)) - 31.1921296296296*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.12962962962963*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 11:
                        return sign * RealGradient(-120.0*eta - 60.0*xi + 450.0*(2*eta + 2)*(0.5*xi + 0.5) - 225.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 65.625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 1350.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 525.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 164.0 + 180.0*((eta + 1)*(eta + 1)) + 675.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 60.0*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 70.0*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 12:
                        return sign * RealGradient(0, 45.0*eta + 120.0*xi - 360.0*(2*eta + 2)*(0.5*xi + 0.5) + 337.5*(2*eta + 2)*((xi + 1)*(xi + 1)) - 225.0*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 420.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 150.0 - 26.25*(eta + 1)*(eta + 1) - 393.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 57.421875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 225.0*(xi + 1)*(xi + 1) + 150.0*((xi + 1)*(xi + 1)*(xi + 1)) - 32.8125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 13:
                        return sign * RealGradient(0, -18.8888888888889*eta - 44.4444444444444*xi + 151.111111111111*(2*eta + 2)*(0.5*xi + 0.5) - 141.666666666667*(2*eta + 2)*(xi + 1)*(xi + 1) + 94.4444444444444*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 20.6597222222222*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 171.111111111111*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 57.7777777777778 + 10.6944444444444*((eta + 1)*(eta + 1)) + 160.416666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 106.944444444444*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 23.3940972222222*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 83.3333333333333*((xi + 1)*(xi + 1)) - 55.5555555555556*(xi + 1)*(xi + 1)*(xi + 1) + 12.1527777777778*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 14:
                        return sign * RealGradient(0, 23.8888888888889*eta + 84.4444444444444*xi - 191.111111111111*(2*eta + 2)*(0.5*xi + 0.5) + 179.166666666667*(2*eta + 2)*((xi + 1)*(xi + 1)) - 119.444444444444*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 26.1284722222222*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 171.111111111111*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 97.7777777777778 - 10.6944444444444*(eta + 1)*(eta + 1) - 160.416666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 106.944444444444*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 23.3940972222222*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 158.333333333333*(xi + 1)*(xi + 1) + 105.555555555556*((xi + 1)*(xi + 1)*(xi + 1)) - 23.0902777777778*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 15:
                        return sign * RealGradient(0, -60.0*eta - 240.0*xi + 480.0*(2*eta + 2)*(0.5*xi + 0.5) - 450.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 300.0*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 420.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 270.0 + 26.25*((eta + 1)*(eta + 1)) + 393.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 57.421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 450.0*((xi + 1)*(xi + 1)) - 300.0*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 16:
                        return RealGradient(0, -195.0*xi + 390.0*(2*eta + 2)*(0.5*xi + 0.5) - 435.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 300.0*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 341.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 195.0 + 380.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 57.421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 435.0*((xi + 1)*(xi + 1)) - 300.0*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 17:
                        return RealGradient(0, -45.0*xi + 90.0*(2*eta + 2)*(0.5*xi + 0.5) - 210.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 225.0*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 78.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 45.0 + 183.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 57.421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 210.0*((xi + 1)*(xi + 1)) - 225.0*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 18:
                        return RealGradient(0, -60.0*xi + 120.0*(2*eta + 2)*(0.5*xi + 0.5) - 30.0*(2*eta + 2)*(xi + 1)*(xi + 1) - 105.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 60.0 + 26.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 30.0*((xi + 1)*(xi + 1)));
                      case 19:
                        return RealGradient(232.0*eta + 195.0*xi - 870.0*(2*eta + 2)*(0.5*xi + 0.5) + 435.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 126.875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1800.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 525.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 375.0 - 240.0*(eta + 1)*(eta + 1) - 900.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 195.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 70.0*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 56.875*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 20:
                        return RealGradient(112.0*eta + 45.0*xi - 420.0*(2*eta + 2)*(0.5*xi + 0.5) + 210.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 61.25*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1350.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 525.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 145.0 - 180.0*(eta + 1)*(eta + 1) - 675.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 45.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 70.0*((eta + 1)*(eta + 1)*(eta + 1)) - 76.5625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 13.125*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 21:
                        return RealGradient(16.0*eta + 60.0*xi - 60.0*(2*eta + 2)*(0.5*xi + 0.5) + 30.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 8.75*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 60.0 - 60.0*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 22:
                        return RealGradient(-58.0*eta - 97.5*xi + 435.0*(2*eta + 2)*(0.5*xi + 0.5) - 326.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 126.875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 900.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 262.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 142.5 + 60.0*((eta + 1)*(eta + 1)) + 675.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 262.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 146.25*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 17.5*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 56.875*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 23:
                        return RealGradient(-28.0*eta - 22.5*xi + 210.0*(2*eta + 2)*(0.5*xi + 0.5) - 157.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 61.25*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 675.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 262.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 47.5 + 45.0*((eta + 1)*(eta + 1)) + 506.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 196.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 33.75*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 17.5*(eta + 1)*(eta + 1)*(eta + 1) + 76.5625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 13.125*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 24:
                        return RealGradient(-4.0*eta - 30.0*xi + 30.0*(2*eta + 2)*(0.5*xi + 0.5) - 22.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 8.75*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 30.0 + 45.0*((xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 25:
                        return RealGradient(0, 97.5*xi - 292.5*(2*eta + 2)*(0.5*xi + 0.5) + 326.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 225.0*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 341.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 97.5 - 380.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 57.421875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 217.5*(xi + 1)*(xi + 1) + 150.0*((xi + 1)*(xi + 1)*(xi + 1)) - 32.8125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 26:
                        return RealGradient(0, 22.5*xi - 67.5*(2*eta + 2)*(0.5*xi + 0.5) + 157.5*(2*eta + 2)*((xi + 1)*(xi + 1)) - 168.75*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 49.21875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 78.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 22.5 - 183.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 196.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 57.421875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 105.0*(xi + 1)*(xi + 1) + 112.5*((xi + 1)*(xi + 1)*(xi + 1)) - 32.8125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 27:
                        return RealGradient(0, 30.0*xi - 90.0*(2*eta + 2)*(0.5*xi + 0.5) + 22.5*(2*eta + 2)*((xi + 1)*(xi + 1)) + 105.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 30.0 - 26.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1));
                      case 28:
                        return RealGradient(-57.0*eta - 9.0*xi + 42.75*(2*eta + 2)*(0.5*xi + 0.5) - 90.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 26.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 54.0 + 60.0*((eta + 1)*(eta + 1)) - 17.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 29:
                        return RealGradient(28.5*eta + 9.0*xi - 42.75*(2*eta + 2)*(0.5*xi + 0.5) + 90.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 26.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 31.5 - 30.0*(eta + 1)*(eta + 1) + 8.75*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 30:
                        return RealGradient(0, -9.0*xi - 9.0 + 21.375*((xi + 1)*(xi + 1)) - 15.0*(xi + 1)*(xi + 1)*(xi + 1) + 3.28125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 31:
                        return RealGradient(0, 9.0*xi + 9.0 - 21.375*(xi + 1)*(xi + 1) + 15.0*((xi + 1)*(xi + 1)*(xi + 1)) - 3.28125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 32:
                        return RealGradient(0, -1.5*xi - 1.5 + 10.125*((xi + 1)*(xi + 1)) - 11.25*(xi + 1)*(xi + 1)*(xi + 1) + 3.28125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 33:
                        return RealGradient(0, 1.5*xi + 1.5 - 10.125*(xi + 1)*(xi + 1) + 11.25*((xi + 1)*(xi + 1)*(xi + 1)) - 3.28125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 34:
                        return RealGradient(-27.0*eta - 1.5*xi + 20.25*(2*eta + 2)*(0.5*xi + 0.5) - 67.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 26.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 26.5 + 45.0*((eta + 1)*(eta + 1)) - 17.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 35:
                        return RealGradient(13.5*eta + 1.5*xi - 20.25*(2*eta + 2)*(0.5*xi + 0.5) + 67.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 26.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 14.0 - 22.5*(eta + 1)*(eta + 1) + 8.75*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 36:
                        return RealGradient(0, -4.5*xi - 4.5 + 2.25*((xi + 1)*(xi + 1)));
                      case 37:
                        return RealGradient(6.0*eta + 4.5*xi - 4.5*(2*eta + 2)*(0.5*xi + 0.5) + 4.5, 0);
                      case 38:
                        return RealGradient(-3.0*eta - 4.5*xi + 4.5*(2*eta + 2)*(0.5*xi + 0.5) - 4.5, 0);
                      case 39:
                        return RealGradient(0, 4.5*xi + 4.5 - 2.25*(xi + 1)*(xi + 1));
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
                        return sign * RealGradient(-1680.0*eta*xi + 720.0*eta + 672.0*eta*(xi*xi) + 480.0*xi + 1344.0*xi*(eta*eta) - 120.0 - 1260.0*eta*eta - 420.0*xi*xi + 672.0*(eta*eta*eta), 1680.0*eta*xi - 240.0*eta - 2016.0*eta*xi*xi - 480.0*xi - 1344.0*xi*eta*eta + 40.0 + 420.0*(eta*eta) + 1260.0*(xi*xi) - 224.0*eta*eta*eta - 896.0*xi*xi*xi);
                      case 1:
                        return sign * RealGradient(640.888888888889*eta*xi - 200.0*eta - 273.777777777778*eta*xi*xi - 191.111111111111*xi - 497.777777777778*xi*eta*eta + 42.2222222222222 + 261.333333333333*(eta*eta) + 171.111111111111*(xi*xi) - 99.5555555555556*eta*eta*eta, -472.888888888889*eta*xi + 19.5555555555556*eta + 746.666666666667*eta*(xi*xi) + 177.777777777778*xi + 199.111111111111*xi*(eta*eta) - 10.962962962963 + 49.7777777777778*(eta*eta) - 504.0*xi*xi - 66.3703703703704*eta*eta*eta + 365.037037037037*(xi*xi*xi));
                      case 2:
                        return sign * RealGradient(-248.888888888889*eta*xi - 16.0*eta + 273.777777777778*eta*(xi*xi) + 151.111111111111*xi + 49.7777777777778*xi*(eta*eta) - 22.2222222222222 + 158.666666666667*(eta*eta) - 171.111111111111*xi*xi - 124.444444444444*eta*eta*eta, -87.1111111111111*eta*xi + 12.4444444444444*eta - 74.6666666666667*eta*xi*xi - 113.777777777778*xi + 248.888888888889*xi*(eta*eta) + 5.62962962962963 - 21.7777777777778*eta*eta + 420.0*(xi*xi) - 8.2962962962963*eta*eta*eta - 365.037037037037*xi*xi*xi);
                      case 3:
                        return sign * RealGradient(504.0*eta*xi - 72.0*eta - 672.0*eta*xi*xi - 360.0*xi + 60.0 + 420.0*(xi*xi), 288.0*xi - 16.0 - 1008.0*xi*xi + 896.0*(xi*xi*xi));
                      case 4:
                        return sign * RealGradient(504.0*eta*xi - 72.0*eta - 672.0*eta*xi*xi, 288.0*xi - 16.0 - 1008.0*xi*xi + 896.0*(xi*xi*xi));
                      case 5:
                        return sign * RealGradient(298.666666666667*eta*xi - 80.0*eta - 199.111111111111*eta*xi*xi - 597.333333333333*xi*eta*eta + 261.333333333333*(eta*eta) - 149.333333333333*eta*eta*eta, -597.333333333333*eta*xi + 56.0*eta + 896.0*eta*(xi*xi) + 160.0*xi + 298.666666666667*xi*(eta*eta) - 12.0 - 46.6666666666667*eta*eta - 410.666666666667*xi*xi + 8.2962962962963*(eta*eta*eta) + 265.481481481481*(xi*xi*xi));
                      case 6:
                        return sign * RealGradient(93.3333333333333*eta*xi - 56.0*eta - 24.8888888888889*eta*xi*xi - 298.666666666667*xi*eta*eta + 298.666666666667*(eta*eta) - 298.666666666667*eta*eta*eta, -522.666666666667*eta*xi + 80.0*eta + 448.0*eta*(xi*xi) + 64.0*xi + 597.333333333333*xi*(eta*eta) - 8.0 - 149.333333333333*eta*eta - 93.3333333333333*xi*xi + 66.3703703703704*(eta*eta*eta) + 33.1851851851852*(xi*xi*xi));
                      case 7:
                        return sign * RealGradient(0, 72.0*eta - 4.0 - 252.0*eta*eta + 224.0*(eta*eta*eta));
                      case 8:
                        return sign * RealGradient(0, 72.0*eta - 4.0 - 252.0*eta*eta + 224.0*(eta*eta*eta));
                      case 9:
                        return sign * RealGradient(43.5555555555556*eta*xi - 12.4444444444444*eta + 24.8888888888889*eta*(xi*xi) - 248.888888888889*xi*eta*eta + 43.5555555555556*(eta*eta) + 24.8888888888889*(eta*eta*eta), -317.333333333333*eta*xi + 16.0*eta + 373.333333333333*eta*(xi*xi) + 23.1111111111111*xi - 49.7777777777778*xi*eta*eta - 4.14814814814815 + 124.444444444444*(eta*eta) + 6.22222222222222*(xi*xi) - 91.2592592592593*eta*eta*eta - 33.1851851851852*xi*xi*xi);
                      case 10:
                        return sign * RealGradient(-99.5555555555556*eta*xi - 19.5555555555556*eta + 199.111111111111*eta*(xi*xi) - 199.111111111111*xi*eta*eta + 236.444444444444*(eta*eta) - 248.888888888889*eta*eta*eta, -522.666666666667*eta*xi + 200.0*eta + 298.666666666667*eta*(xi*xi) - 135.111111111111*xi + 497.777777777778*xi*(eta*eta) + 2.81481481481481 - 320.444444444444*eta*eta + 385.777777777778*(xi*xi) + 91.2592592592593*(eta*eta*eta) - 265.481481481481*xi*xi*xi);
                      case 11:
                        return sign * RealGradient(-840.0*eta*xi + 240.0*eta + 672.0*eta*(xi*xi) + 1344.0*xi*(eta*eta) - 840.0*eta*eta + 672.0*(eta*eta*eta), 2520.0*eta*xi - 720.0*eta - 2016.0*eta*xi*xi - 960.0*xi - 1344.0*xi*eta*eta + 160.0 + 840.0*(eta*eta) + 1680.0*(xi*xi) - 224.0*eta*eta*eta - 896.0*xi*xi*xi);
                      case 12:
                        return RealGradient(6048.0*eta*xi - 3240.0*eta - 2016.0*eta*xi*xi - 8064.0*xi*eta*eta + 9072.0*(eta*eta) - 6048.0*eta*eta*eta, -12096.0*eta*xi + 2160.0*eta + 12096.0*eta*(xi*xi) + 2160.0*xi + 12096.0*xi*(eta*eta) - 240.0 - 4536.0*eta*eta - 4536.0*xi*xi + 2688.0*(eta*eta*eta) + 2688.0*(xi*xi*xi));
                      case 13:
                        return RealGradient(-9072.0*eta*xi + 2160.0*eta + 8064.0*eta*(xi*xi) + 12096.0*xi*(eta*eta) - 6048.0*eta*eta + 4032.0*(eta*eta*eta), 18144.0*eta*xi - 3240.0*eta - 18144.0*eta*xi*xi - 8640.0*xi - 8064.0*xi*eta*eta + 960.0 + 3024.0*(eta*eta) + 18144.0*(xi*xi) - 672.0*eta*eta*eta - 10752.0*xi*xi*xi);
                      case 14:
                        return RealGradient(9072.0*eta*xi - 1944.0*eta - 6048.0*eta*xi*xi - 8064.0*xi*eta*eta + 2016.0*(eta*eta), -6048.0*eta*xi + 432.0*eta + 12096.0*eta*(xi*xi) + 3456.0*xi - 216.0 - 10584.0*xi*xi + 8064.0*(xi*xi*xi));
                      case 15:
                        return RealGradient(-7056.0*eta*xi + 1152.0*eta + 8064.0*eta*(xi*xi) + 4032.0*xi*(eta*eta) - 1008.0*eta*eta, 3024.0*eta*xi - 216.0*eta - 6048.0*eta*xi*xi - 4608.0*xi + 288.0 + 14112.0*(xi*xi) - 10752.0*xi*xi*xi);
                      case 16:
                        return RealGradient(-216.0*eta + 1512.0*(eta*eta) - 2016.0*eta*eta*eta, -2016.0*eta*xi + 1152.0*eta + 144.0*xi + 4032.0*xi*(eta*eta) - 72.0 - 3528.0*eta*eta + 2688.0*(eta*eta*eta));
                      case 17:
                        return RealGradient(432.0*eta - 3024.0*eta*eta + 4032.0*(eta*eta*eta), 4032.0*eta*xi - 1944.0*eta - 288.0*xi - 8064.0*xi*eta*eta + 144.0 + 4536.0*(eta*eta) - 2016.0*eta*eta*eta);
                      case 18:
                        return RealGradient(3276.0*eta*xi - 1332.0*eta - 1512.0*eta*xi*xi - 6048.0*xi*eta*eta + 4788.0*(eta*eta) - 3528.0*eta*eta*eta, -8064.0*eta*xi + 1044.0*eta + 9072.0*eta*(xi*xi) + 1548.0*xi + 7056.0*xi*(eta*eta) - 144.0 - 1638.0*eta*eta - 3402.0*xi*xi + 672.0*(eta*eta*eta) + 2016.0*(xi*xi*xi));
                      case 19:
                        return RealGradient(-3276.0*eta*xi + 1044.0*eta + 2016.0*eta*(xi*xi) + 7056.0*xi*(eta*eta) - 4032.0*eta*eta + 3024.0*(eta*eta*eta), 9576.0*eta*xi - 1332.0*eta - 10584.0*eta*xi*xi - 2196.0*xi - 6048.0*xi*eta*eta + 216.0 + 1638.0*(eta*eta) + 4662.0*(xi*xi) - 504.0*eta*eta*eta - 2688.0*xi*xi*xi);
                      case 20:
                        return RealGradient(1008.0*eta*xi - 216.0*eta - 504.0*eta*xi*xi - 756.0*eta*eta + 1008.0*(eta*eta*eta), 1008.0*eta*xi - 360.0*eta + 144.0*xi - 2016.0*xi*eta*eta + 12.0 + 1008.0*(eta*eta) - 756.0*xi*xi - 672.0*eta*eta*eta + 672.0*(xi*xi*xi));
                      case 21:
                        return RealGradient(-756.0*eta*xi - 216.0*eta + 2016.0*eta*(xi*xi) - 3024.0*xi*eta*eta + 2268.0*(eta*eta) - 2016.0*eta*eta*eta, -5544.0*eta*xi + 1188.0*eta + 4536.0*eta*(xi*xi) - 936.0*xi + 4032.0*xi*(eta*eta) + 6.0 - 1512.0*eta*eta + 3402.0*(xi*xi) + 336.0*(eta*eta*eta) - 2688.0*xi*xi*xi);
                      case 22:
                        return RealGradient(-3024.0*eta*xi + 1188.0*eta + 1008.0*eta*(xi*xi) + 4032.0*xi*(eta*eta) - 2772.0*eta*eta + 1512.0*(eta*eta*eta), 4536.0*eta*xi - 216.0*eta - 6048.0*eta*xi*xi - 972.0*xi - 3024.0*xi*eta*eta + 66.0 - 378.0*eta*eta + 2268.0*(xi*xi) + 672.0*(eta*eta*eta) - 1344.0*xi*xi*xi);
                      case 23:
                        return RealGradient(2016.0*eta*xi - 360.0*eta - 2016.0*eta*xi*xi - 2016.0*xi*eta*eta + 504.0*(eta*eta), -1512.0*eta*xi - 216.0*eta + 3024.0*eta*(xi*xi) + 1440.0*xi - 48.0 + 504.0*(eta*eta) - 4032.0*xi*xi - 168.0*eta*eta*eta + 2688.0*(xi*xi*xi));
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
                        return sign * RealGradient(-2520.0*eta*xi + 960.0*eta + 1344.0*eta*(xi*xi) + 720.0*xi + 2016.0*xi*(eta*eta) - 160.0 - 1680.0*eta*eta - 840.0*xi*xi + 896.0*(eta*eta*eta) + 224.0*(xi*xi*xi), 840.0*eta*xi - 1344.0*eta*xi*xi - 240.0*xi - 672.0*xi*eta*eta + 840.0*(xi*xi) - 672.0*xi*xi*xi);
                      case 1:
                        return sign * RealGradient(522.666666666667*eta*xi + 135.111111111111*eta - 497.777777777778*eta*xi*xi - 200.0*xi - 298.666666666667*xi*eta*eta - 2.81481481481481 - 385.777777777778*eta*eta + 320.444444444444*(xi*xi) + 265.481481481481*(eta*eta*eta) - 91.2592592592593*xi*xi*xi, 99.5555555555556*eta*xi + 199.111111111111*eta*(xi*xi) + 19.5555555555556*xi - 199.111111111111*xi*eta*eta - 236.444444444444*xi*xi + 248.888888888889*(xi*xi*xi));
                      case 2:
                        return sign * RealGradient(317.333333333333*eta*xi - 23.1111111111111*eta + 49.7777777777778*eta*(xi*xi) - 16.0*xi - 373.333333333333*xi*eta*eta + 4.14814814814815 - 6.22222222222222*eta*eta - 124.444444444444*xi*xi + 33.1851851851852*(eta*eta*eta) + 91.2592592592593*(xi*xi*xi), -43.5555555555556*eta*xi + 248.888888888889*eta*(xi*xi) + 12.4444444444444*xi - 24.8888888888889*xi*eta*eta - 43.5555555555556*xi*xi - 24.8888888888889*xi*xi*xi);
                      case 3:
                        return sign * RealGradient(-72.0*xi + 4.0 + 252.0*(xi*xi) - 224.0*xi*xi*xi, 0);
                      case 4:
                        return sign * RealGradient(-72.0*xi + 4.0 + 252.0*(xi*xi) - 224.0*xi*xi*xi, 0);
                      case 5:
                        return sign * RealGradient(522.666666666667*eta*xi - 64.0*eta - 597.333333333333*eta*xi*xi - 80.0*xi - 448.0*xi*eta*eta + 8.0 + 93.3333333333333*(eta*eta) + 149.333333333333*(xi*xi) - 33.1851851851852*eta*eta*eta - 66.3703703703704*xi*xi*xi, -93.3333333333333*eta*xi + 298.666666666667*eta*(xi*xi) + 56.0*xi + 24.8888888888889*xi*(eta*eta) - 298.666666666667*xi*xi + 298.666666666667*(xi*xi*xi));
                      case 6:
                        return sign * RealGradient(597.333333333333*eta*xi - 160.0*eta - 298.666666666667*eta*xi*xi - 56.0*xi - 896.0*xi*eta*eta + 12.0 + 410.666666666667*(eta*eta) + 46.6666666666667*(xi*xi) - 265.481481481481*eta*eta*eta - 8.2962962962963*xi*xi*xi, -298.666666666667*eta*xi + 597.333333333333*eta*(xi*xi) + 80.0*xi + 199.111111111111*xi*(eta*eta) - 261.333333333333*xi*xi + 149.333333333333*(xi*xi*xi));
                      case 7:
                        return sign * RealGradient(-288.0*eta + 16.0 + 1008.0*(eta*eta) - 896.0*eta*eta*eta, -504.0*eta*xi + 72.0*xi + 672.0*xi*(eta*eta));
                      case 8:
                        return sign * RealGradient(-288.0*eta + 16.0 + 1008.0*(eta*eta) - 896.0*eta*eta*eta, -504.0*eta*xi + 360.0*eta + 72.0*xi + 672.0*xi*(eta*eta) - 60.0 - 420.0*eta*eta);
                      case 9:
                        return sign * RealGradient(87.1111111111111*eta*xi + 113.777777777778*eta - 248.888888888889*eta*xi*xi - 12.4444444444444*xi + 74.6666666666667*xi*(eta*eta) - 5.62962962962963 - 420.0*eta*eta + 21.7777777777778*(xi*xi) + 365.037037037037*(eta*eta*eta) + 8.2962962962963*(xi*xi*xi), 248.888888888889*eta*xi - 151.111111111111*eta - 49.7777777777778*eta*xi*xi + 16.0*xi - 273.777777777778*xi*eta*eta + 22.2222222222222 + 171.111111111111*(eta*eta) - 158.666666666667*xi*xi + 124.444444444444*(xi*xi*xi));
                      case 10:
                        return sign * RealGradient(472.888888888889*eta*xi - 177.777777777778*eta - 199.111111111111*eta*xi*xi - 19.5555555555556*xi - 746.666666666667*xi*eta*eta + 10.962962962963 + 504.0*(eta*eta) - 49.7777777777778*xi*xi - 365.037037037037*eta*eta*eta + 66.3703703703704*(xi*xi*xi), -640.888888888889*eta*xi + 191.111111111111*eta + 497.777777777778*eta*(xi*xi) + 200.0*xi + 273.777777777778*xi*(eta*eta) - 42.2222222222222 - 171.111111111111*eta*eta - 261.333333333333*xi*xi + 99.5555555555556*(xi*xi*xi));
                      case 11:
                        return sign * RealGradient(-1680.0*eta*xi + 480.0*eta + 1344.0*eta*(xi*xi) + 240.0*xi + 2016.0*xi*(eta*eta) - 40.0 - 1260.0*eta*eta - 420.0*xi*xi + 896.0*(eta*eta*eta) + 224.0*(xi*xi*xi), 1680.0*eta*xi - 480.0*eta - 1344.0*eta*xi*xi - 720.0*xi - 672.0*xi*eta*eta + 120.0 + 420.0*(eta*eta) + 1260.0*(xi*xi) - 672.0*xi*xi*xi);
                      case 12:
                        return RealGradient(18144.0*eta*xi - 8640.0*eta - 8064.0*eta*xi*xi - 3240.0*xi - 18144.0*xi*eta*eta + 960.0 + 18144.0*(eta*eta) + 3024.0*(xi*xi) - 10752.0*eta*eta*eta - 672.0*xi*xi*xi, -9072.0*eta*xi + 12096.0*eta*(xi*xi) + 2160.0*xi + 8064.0*xi*(eta*eta) - 6048.0*xi*xi + 4032.0*(xi*xi*xi));
                      case 13:
                        return RealGradient(-12096.0*eta*xi + 2160.0*eta + 12096.0*eta*(xi*xi) + 2160.0*xi + 12096.0*xi*(eta*eta) - 240.0 - 4536.0*eta*eta - 4536.0*xi*xi + 2688.0*(eta*eta*eta) + 2688.0*(xi*xi*xi), 6048.0*eta*xi - 8064.0*eta*xi*xi - 3240.0*xi - 2016.0*xi*eta*eta + 9072.0*(xi*xi) - 6048.0*xi*xi*xi);
                      case 14:
                        return RealGradient(4032.0*eta*xi - 288.0*eta - 8064.0*eta*xi*xi - 1944.0*xi + 144.0 + 4536.0*(xi*xi) - 2016.0*xi*xi*xi, 432.0*xi - 3024.0*xi*xi + 4032.0*(xi*xi*xi));
                      case 15:
                        return RealGradient(-2016.0*eta*xi + 144.0*eta + 4032.0*eta*(xi*xi) + 1152.0*xi - 72.0 - 3528.0*xi*xi + 2688.0*(xi*xi*xi), -216.0*xi + 1512.0*(xi*xi) - 2016.0*xi*xi*xi);
                      case 16:
                        return RealGradient(3024.0*eta*xi - 4608.0*eta - 216.0*xi - 6048.0*xi*eta*eta + 288.0 + 14112.0*(eta*eta) - 10752.0*eta*eta*eta, -7056.0*eta*xi + 4032.0*eta*(xi*xi) + 1152.0*xi + 8064.0*xi*(eta*eta) - 1008.0*xi*xi);
                      case 17:
                        return RealGradient(-6048.0*eta*xi + 3456.0*eta + 432.0*xi + 12096.0*xi*(eta*eta) - 216.0 - 10584.0*eta*eta + 8064.0*(eta*eta*eta), 9072.0*eta*xi - 8064.0*eta*xi*xi - 1944.0*xi - 6048.0*xi*eta*eta + 2016.0*(xi*xi));
                      case 18:
                        return RealGradient(9576.0*eta*xi - 2196.0*eta - 6048.0*eta*xi*xi - 1332.0*xi - 10584.0*xi*eta*eta + 216.0 + 4662.0*(eta*eta) + 1638.0*(xi*xi) - 2688.0*eta*eta*eta - 504.0*xi*xi*xi, -3276.0*eta*xi + 7056.0*eta*(xi*xi) + 1044.0*xi + 2016.0*xi*(eta*eta) - 4032.0*xi*xi + 3024.0*(xi*xi*xi));
                      case 19:
                        return RealGradient(-8064.0*eta*xi + 1548.0*eta + 7056.0*eta*(xi*xi) + 1044.0*xi + 9072.0*xi*(eta*eta) - 144.0 - 3402.0*eta*eta - 1638.0*xi*xi + 2016.0*(eta*eta*eta) + 672.0*(xi*xi*xi), 3276.0*eta*xi - 6048.0*eta*xi*xi - 1332.0*xi - 1512.0*xi*eta*eta + 4788.0*(xi*xi) - 3528.0*xi*xi*xi);
                      case 20:
                        return RealGradient(-1512.0*eta*xi + 1440.0*eta - 216.0*xi + 3024.0*xi*(eta*eta) - 48.0 - 4032.0*eta*eta + 504.0*(xi*xi) + 2688.0*(eta*eta*eta) - 168.0*xi*xi*xi, 2016.0*eta*xi - 2016.0*eta*xi*xi - 360.0*xi - 2016.0*xi*eta*eta + 504.0*(xi*xi));
                      case 21:
                        return RealGradient(4536.0*eta*xi - 972.0*eta - 3024.0*eta*xi*xi - 216.0*xi - 6048.0*xi*eta*eta + 66.0 + 2268.0*(eta*eta) - 378.0*xi*xi - 1344.0*eta*eta*eta + 672.0*(xi*xi*xi), -3024.0*eta*xi + 4032.0*eta*(xi*xi) + 1188.0*xi + 1008.0*xi*(eta*eta) - 2772.0*xi*xi + 1512.0*(xi*xi*xi));
                      case 22:
                        return RealGradient(-5544.0*eta*xi - 936.0*eta + 4032.0*eta*(xi*xi) + 1188.0*xi + 4536.0*xi*(eta*eta) + 6.0 + 3402.0*(eta*eta) - 1512.0*xi*xi - 2688.0*eta*eta*eta + 336.0*(xi*xi*xi), -756.0*eta*xi - 3024.0*eta*xi*xi - 216.0*xi + 2016.0*xi*(eta*eta) + 2268.0*(xi*xi) - 2016.0*xi*xi*xi);
                      case 23:
                        return RealGradient(1008.0*eta*xi + 144.0*eta - 2016.0*eta*xi*xi - 360.0*xi + 12.0 - 756.0*eta*eta + 1008.0*(xi*xi) + 672.0*(eta*eta*eta) - 672.0*xi*xi*xi, 1008.0*eta*xi - 216.0*xi - 504.0*xi*eta*eta - 756.0*xi*xi + 1008.0*(xi*xi*xi));
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
      } // end case FOURTH

      // quintic Nedelec (first kind) shape function first derivatives
    case FIFTH:
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
                        return sign * RealGradient(937.5*eta + 262.5*xi - 3281.25*(0.5*eta + 0.5)*(2*xi + 2) + 6562.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1968.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 4921.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 5742.1875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 2871.09375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 516.796875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1125.0 - 2812.5*(eta + 1)*(eta + 1) - 9843.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2953.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1) + 11484.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5742.1875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1033.59375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 3281.25*((eta + 1)*(eta + 1)*(eta + 1)) - 3445.3125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 78.75*((xi + 1)*(xi + 1)*(xi + 1)) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1640.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 295.3125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 1:
                        return sign * RealGradient(-279.78515625*eta - 82.44140625*xi + 1030.517578125*(0.5*eta + 0.5)*(2*xi + 2) - 1994.384765625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 569.091796875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1545.7763671875*(2*xi + 2)*(eta + 1)*(eta + 1) + 1803.40576171875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 901.702880859375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 162.306518554688*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 339.84375 + 839.35546875*((eta + 1)*(eta + 1)) + 2991.5771484375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 853.6376953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 79.775390625*((xi + 1)*(xi + 1)) - 3490.17333984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1745.08666992188*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 314.115600585938*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 979.248046875*(eta + 1)*(eta + 1)*(eta + 1) + 995.91064453125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 22.763671875*(xi + 1)*(xi + 1)*(xi + 1) - 497.955322265625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 89.6319580078125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 489.6240234375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 88.13232421875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 2:
                        return sign * RealGradient(164.0625*eta + 72.1875*xi - 902.34375*(0.5*eta + 0.5)*(2*xi + 2) + 2214.84375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 738.28125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1353.515625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1579.1015625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 789.55078125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 142.119140625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 223.125 - 492.1875*(eta + 1)*(eta + 1) - 3322.265625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1107.421875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(xi + 1)*(xi + 1) + 3875.9765625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1937.98828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 348.837890625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 574.21875*((eta + 1)*(eta + 1)*(eta + 1)) - 1291.9921875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 29.53125*((xi + 1)*(xi + 1)*(xi + 1)) + 645.99609375*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 116.279296875*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 287.109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 51.6796875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 3:
                        return sign * RealGradient(-68.84765625*eta - 36.50390625*xi + 456.298828125*(0.5*eta + 0.5)*(2*xi + 2) - 1420.166015625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 569.091796875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 684.4482421875*(2*xi + 2)*(eta + 1)*(eta + 1) + 798.52294921875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 399.261474609375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 71.8670654296875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 99.84375 + 206.54296875*((eta + 1)*(eta + 1)) + 2130.2490234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 853.6376953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 56.806640625*((xi + 1)*(xi + 1)) - 2485.29052734375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1242.64526367188*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 223.676147460938*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 240.966796875*(eta + 1)*(eta + 1)*(eta + 1) + 995.91064453125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 22.763671875*(xi + 1)*(xi + 1)*(xi + 1) - 497.955322265625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 89.6319580078125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 120.4833984375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 21.68701171875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 4:
                        return sign * RealGradient(375.0*eta + 157.5*xi - 1968.75*(0.5*eta + 0.5)*(2*xi + 2) + 5250.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1968.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 2953.125*(2*xi + 2)*((eta + 1)*(eta + 1)) - 3445.3125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1722.65625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 502.5 - 1125.0*(eta + 1)*(eta + 1) - 7875.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2953.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 210.0*(xi + 1)*(xi + 1) + 9187.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4593.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 826.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1312.5*((eta + 1)*(eta + 1)*(eta + 1)) - 3445.3125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 78.75*((xi + 1)*(xi + 1)*(xi + 1)) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 656.25*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 118.125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 5:
                        return sign * RealGradient(0, -187.5*eta - 375.0*xi + 2250.0*(0.5*eta + 0.5)*(2*xi + 2) - 11812.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 10500.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 2953.125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1968.75*(2*xi + 2)*(eta + 1)*(eta + 1) + 1312.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 295.3125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 531.25 + 328.125*((eta + 1)*(eta + 1)) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 9187.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2583.984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 984.375*((xi + 1)*(xi + 1)) - 6890.625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1550.390625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 218.75*(eta + 1)*(eta + 1)*(eta + 1) + 6125.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 875.0*(xi + 1)*(xi + 1)*(xi + 1) - 1378.125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 49.21875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 246.09375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 6:
                        return sign * RealGradient(0, 55.95703125*eta + 56.8359375*xi - 671.484375*(0.5*eta + 0.5)*(2*xi + 2) + 3525.29296875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 3133.59375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 881.3232421875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 618.310546875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 398.876953125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 85.36376953125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 108.056640625 - 103.0517578125*(eta + 1)*(eta + 1) - 3246.13037109375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2885.44921875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 811.532592773438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 149.1943359375*(xi + 1)*(xi + 1) + 2094.10400390625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 448.159790039063*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 66.4794921875*((eta + 1)*(eta + 1)*(eta + 1)) - 1861.42578125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 523.526000976563*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 132.6171875*((xi + 1)*(xi + 1)*(xi + 1)) + 398.3642578125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.227294921875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 112.039947509766*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 37.298583984375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 7:
                        return sign * RealGradient(0, -32.8125*eta - 28.125*xi + 393.75*(0.5*eta + 0.5)*(2*xi + 2) - 2067.1875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1837.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 541.40625*(2*xi + 2)*(eta + 1)*(eta + 1) + 442.96875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 110.7421875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 58.59375 + 90.234375*((eta + 1)*(eta + 1)) + 2842.3828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2526.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 710.595703125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 73.828125*((xi + 1)*(xi + 1)) - 2325.5859375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 581.396484375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 73.828125*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 581.396484375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 65.625*(xi + 1)*(xi + 1)*(xi + 1) - 516.796875*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 18.45703125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 145.34912109375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 18.45703125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 8:
                        return sign * RealGradient(0, 13.76953125*eta + 9.9609375*xi - 165.234375*(0.5*eta + 0.5)*(2*xi + 2) + 867.48046875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 771.09375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 216.8701171875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 273.779296875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 284.033203125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 85.36376953125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 22.900390625 - 45.6298828125*(eta + 1)*(eta + 1) - 1437.34130859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1277.63671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 359.335327148438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 26.1474609375*(xi + 1)*(xi + 1) + 1491.17431640625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 448.159790039063*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 47.3388671875*((eta + 1)*(eta + 1)*(eta + 1)) - 1325.48828125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 372.793579101563*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 23.2421875*((xi + 1)*(xi + 1)*(xi + 1)) + 398.3642578125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.227294921875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 112.039947509766*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6.536865234375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 9:
                        return sign * RealGradient(0, -75.0*eta - 75.0*xi + 900.0*(0.5*eta + 0.5)*(2*xi + 2) - 4725.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 4200.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1181.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1181.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 1050.0*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 295.3125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 143.75 + 196.875*((eta + 1)*(eta + 1)) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 5512.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 196.875*((xi + 1)*(xi + 1)) - 5512.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1550.390625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 175.0*(eta + 1)*(eta + 1)*(eta + 1) + 4900.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 175.0*(xi + 1)*(xi + 1)*(xi + 1) - 1378.125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 49.21875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 49.21875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 10:
                        return sign * RealGradient(75.0*eta - 393.75*(0.5*eta + 0.5)*(2*xi + 2) + 1050.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 393.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1181.25*(2*xi + 2)*((eta + 1)*(eta + 1)) - 2067.1875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1378.125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 75.0 - 450.0*(eta + 1)*(eta + 1) - 3150.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1181.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 5512.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 3675.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 826.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1)*(eta + 1)) - 2067.1875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 525.0*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 118.125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 11:
                        return sign * RealGradient(-13.76953125*eta + 91.259765625*(0.5*eta + 0.5)*(2*xi + 2) - 284.033203125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 113.818359375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 273.779296875*(2*xi + 2)*(eta + 1)*(eta + 1) + 479.11376953125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 319.4091796875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 71.8670654296875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 13.76953125 + 82.6171875*((eta + 1)*(eta + 1)) + 852.099609375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 341.455078125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1491.17431640625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 994.1162109375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 223.676147460938*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 144.580078125*(eta + 1)*(eta + 1)*(eta + 1) + 597.54638671875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 398.3642578125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 89.6319580078125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 96.38671875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 21.68701171875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 12:
                        return sign * RealGradient(32.8125*eta - 180.46875*(0.5*eta + 0.5)*(2*xi + 2) + 442.96875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 147.65625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 541.40625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 947.4609375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 631.640625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 142.119140625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 32.8125 - 196.875*(eta + 1)*(eta + 1) - 1328.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 442.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 2325.5859375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 348.837890625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 344.53125*((eta + 1)*(eta + 1)*(eta + 1)) - 775.1953125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 116.279296875*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 229.6875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 51.6796875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 13:
                        return sign * RealGradient(-55.95703125*eta + 206.103515625*(0.5*eta + 0.5)*(2*xi + 2) - 398.876953125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 113.818359375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 618.310546875*(2*xi + 2)*(eta + 1)*(eta + 1) + 1082.04345703125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 721.3623046875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 162.306518554688*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 55.95703125 + 335.7421875*((eta + 1)*(eta + 1)) + 1196.630859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 341.455078125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2094.10400390625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1396.0693359375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 314.115600585938*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 587.548828125*(eta + 1)*(eta + 1)*(eta + 1) + 597.54638671875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 398.3642578125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 89.6319580078125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 391.69921875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 88.13232421875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 14:
                        return sign * RealGradient(187.5*eta - 656.25*(0.5*eta + 0.5)*(2*xi + 2) + 1312.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 393.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1968.75*(2*xi + 2)*((eta + 1)*(eta + 1)) - 3445.3125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 2296.875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 516.796875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 187.5 - 1125.0*(eta + 1)*(eta + 1) - 3937.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1181.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 6890.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4593.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1033.59375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1968.75*((eta + 1)*(eta + 1)*(eta + 1)) - 2067.1875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1312.5*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 295.3125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 15:
                        return sign * RealGradient(0, -375.0*eta - 187.5*xi + 2250.0*(0.5*eta + 0.5)*(2*xi + 2) - 7875.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 5250.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1181.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2953.125*(2*xi + 2)*(eta + 1)*(eta + 1) + 2625.0*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 738.28125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 531.25 + 984.375*((eta + 1)*(eta + 1)) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6890.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 328.125*((xi + 1)*(xi + 1)) - 9187.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2583.984375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 875.0*(eta + 1)*(eta + 1)*(eta + 1) + 6125.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 218.75*(xi + 1)*(xi + 1)*(xi + 1) - 1722.65625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 246.09375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 49.21875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 16:
                        return sign * RealGradient(0, 68.84765625*eta + 24.90234375*xi - 413.0859375*(0.5*eta + 0.5)*(2*xi + 2) + 1445.80078125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 963.8671875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 216.8701171875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 684.4482421875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 710.0830078125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 213.409423828125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 89.599609375 - 228.1494140625*(eta + 1)*(eta + 1) - 2395.56884765625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1597.0458984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 359.335327148438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 43.5791015625*(xi + 1)*(xi + 1) + 2485.29052734375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 746.932983398438*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 236.6943359375*((eta + 1)*(eta + 1)*(eta + 1)) - 1656.8603515625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 372.793579101563*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 29.052734375*((xi + 1)*(xi + 1)*(xi + 1)) + 497.955322265625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 71.136474609375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 112.039947509766*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6.536865234375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 17:
                        return sign * RealGradient(0, -164.0625*eta - 70.3125*xi + 984.375*(0.5*eta + 0.5)*(2*xi + 2) - 3445.3125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 2296.875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1353.515625*(2*xi + 2)*(eta + 1)*(eta + 1) + 1107.421875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 276.85546875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 222.65625 + 451.171875*((eta + 1)*(eta + 1)) + 4737.3046875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3158.203125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 710.595703125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 123.046875*((xi + 1)*(xi + 1)) - 3875.9765625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 968.994140625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 369.140625*(eta + 1)*(eta + 1)*(eta + 1) + 2583.984375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 581.396484375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 82.03125*(xi + 1)*(xi + 1)*(xi + 1) - 645.99609375*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 92.28515625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 145.34912109375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 18.45703125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 18:
                        return sign * RealGradient(0, 279.78515625*eta + 142.08984375*xi - 1678.7109375*(0.5*eta + 0.5)*(2*xi + 2) + 5875.48828125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 3916.9921875*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 881.3232421875*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1545.7763671875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 997.1923828125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 213.409423828125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 398.193359375 - 515.2587890625*(eta + 1)*(eta + 1) - 5410.21728515625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3606.8115234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 811.532592773438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 248.6572265625*(xi + 1)*(xi + 1) + 3490.17333984375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 746.932983398438*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 332.3974609375*((eta + 1)*(eta + 1)*(eta + 1)) - 2326.7822265625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 523.526000976563*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 165.771484375*((xi + 1)*(xi + 1)*(xi + 1)) + 497.955322265625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 71.136474609375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 112.039947509766*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 37.298583984375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 19:
                        return sign * RealGradient(0, -937.5*eta - 937.5*xi + 5625.0*(0.5*eta + 0.5)*(2*xi + 2) - 19687.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 13125.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 2953.125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 4921.875*(2*xi + 2)*(eta + 1)*(eta + 1) + 3281.25*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 738.28125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1718.75 + 1640.625*((eta + 1)*(eta + 1)) + 17226.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 11484.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2583.984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1640.625*((xi + 1)*(xi + 1)) - 11484.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2583.984375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1093.75*(eta + 1)*(eta + 1)*(eta + 1) + 7656.25*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1093.75*(xi + 1)*(xi + 1)*(xi + 1) - 1722.65625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 246.09375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 246.09375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 20:
                        return RealGradient(0, -731.25*eta - 881.25*xi + 5287.5*(0.5*eta + 0.5)*(2*xi + 2) - 19434.375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 13125.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 2953.125*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 4626.5625*(2*xi + 2)*(eta + 1)*(eta + 1) + 3084.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 693.984375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1490.625 + 1279.6875*((eta + 1)*(eta + 1)) + 17005.078125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 11484.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2583.984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1619.53125*((xi + 1)*(xi + 1)) - 11336.71875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2550.76171875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 853.125*(eta + 1)*(eta + 1)*(eta + 1) + 7656.25*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1093.75*(xi + 1)*(xi + 1)*(xi + 1) - 1722.65625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 191.953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 246.09375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 21:
                        return RealGradient(0, 150.0*eta + 346.875*xi - 2081.25*(0.5*eta + 0.5)*(2*xi + 2) + 11559.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 10500.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 2953.125*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1821.09375*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1214.0625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 273.1640625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 471.875 - 262.5*(eta + 1)*(eta + 1) - 10114.453125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 9187.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2583.984375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 963.28125*(xi + 1)*(xi + 1) + 6742.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1517.16796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 175.0*((eta + 1)*(eta + 1)*(eta + 1)) - 6125.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 875.0*((xi + 1)*(xi + 1)*(xi + 1)) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 39.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 246.09375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 22:
                        return RealGradient(0, -337.5*eta - 140.625*xi + 843.75*(0.5*eta + 0.5)*(2*xi + 2) - 759.375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 738.28125*(2*xi + 2)*(eta + 1)*(eta + 1) + 492.1875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 110.7421875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 421.875 + 590.625*((eta + 1)*(eta + 1)) + 664.453125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 63.28125*((xi + 1)*(xi + 1)) - 442.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 99.66796875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 393.75*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 23:
                        return RealGradient(0, 168.75*eta + 112.5*xi - 675.0*(0.5*eta + 0.5)*(2*xi + 2) + 759.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 590.625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 393.75*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 253.125 - 295.3125*(eta + 1)*(eta + 1) - 664.453125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 63.28125*(xi + 1)*(xi + 1) + 442.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 99.66796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 196.875*((eta + 1)*(eta + 1)*(eta + 1)) - 44.296875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 24:
                        return RealGradient(731.25*eta - 2559.375*(0.5*eta + 0.5)*(2*xi + 2) + 5118.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1535.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 4626.5625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 5668.359375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 2871.09375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 516.796875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 731.25 - 2643.75*(eta + 1)*(eta + 1) - 9253.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2775.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 11336.71875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5742.1875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1033.59375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 3239.0625*((eta + 1)*(eta + 1)*(eta + 1)) - 3401.015625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1640.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 295.3125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 25:
                        return RealGradient(-150.0*eta + 525.0*(0.5*eta + 0.5)*(2*xi + 2) - 1050.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 315.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1821.09375*(2*xi + 2)*(eta + 1)*(eta + 1) + 3371.484375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 2296.875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 516.796875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 150.0 + 1040.625*((eta + 1)*(eta + 1)) + 3642.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1092.65625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6742.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4593.75*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1033.59375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1926.5625*(eta + 1)*(eta + 1)*(eta + 1) + 2022.890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 310.078125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1312.5*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 295.3125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 26:
                        return RealGradient(337.5*eta - 1181.25*(0.5*eta + 0.5)*(2*xi + 2) + 2362.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 708.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 738.28125*(2*xi + 2)*((eta + 1)*(eta + 1)) - 221.484375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 337.5 - 421.875*(eta + 1)*(eta + 1) - 1476.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 442.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 442.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 126.5625*((eta + 1)*(eta + 1)*(eta + 1)) - 132.890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 27:
                        return RealGradient(-168.75*eta + 590.625*(0.5*eta + 0.5)*(2*xi + 2) - 1181.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 354.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(2*xi + 2)*(eta + 1)*(eta + 1) + 221.484375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 168.75 + 337.5*((eta + 1)*(eta + 1)) + 1181.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 354.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 442.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 126.5625*(eta + 1)*(eta + 1)*(eta + 1) + 132.890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 28:
                        return RealGradient(292.5*eta - 1535.625*(0.5*eta + 0.5)*(2*xi + 2) + 4095.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1535.625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 2775.9375*(2*xi + 2)*((eta + 1)*(eta + 1)) - 3401.015625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1722.65625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 292.5 - 1057.5*(eta + 1)*(eta + 1) - 7402.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2775.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 9069.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4593.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 826.875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1295.625*((eta + 1)*(eta + 1)*(eta + 1)) - 3401.015625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 310.078125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 656.25*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 118.125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 29:
                        return RealGradient(-60.0*eta + 315.0*(0.5*eta + 0.5)*(2*xi + 2) - 840.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 315.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1092.65625*(2*xi + 2)*(eta + 1)*(eta + 1) + 2022.890625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1378.125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 310.078125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 60.0 + 416.25*((eta + 1)*(eta + 1)) + 2913.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1092.65625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5394.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3675.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 826.875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 770.625*(eta + 1)*(eta + 1)*(eta + 1) + 2022.890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 310.078125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 525.0*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 118.125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 30:
                        return RealGradient(135.0*eta - 708.75*(0.5*eta + 0.5)*(2*xi + 2) + 1890.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 708.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 442.96875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 132.890625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 135.0 - 168.75*(eta + 1)*(eta + 1) - 1181.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 442.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) + 354.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 50.625*((eta + 1)*(eta + 1)*(eta + 1)) - 132.890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 31:
                        return RealGradient(-67.5*eta + 354.375*(0.5*eta + 0.5)*(2*xi + 2) - 945.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 354.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 354.375*(2*xi + 2)*(eta + 1)*(eta + 1) + 132.890625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 67.5 + 135.0*((eta + 1)*(eta + 1)) + 945.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 354.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 354.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 50.625*(eta + 1)*(eta + 1)*(eta + 1) + 132.890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 32:
                        return RealGradient(0, -292.5*eta - 176.25*xi + 2115.0*(0.5*eta + 0.5)*(2*xi + 2) - 7773.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 5250.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 1181.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2775.9375*(2*xi + 2)*(eta + 1)*(eta + 1) + 2467.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 693.984375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 444.375 + 767.8125*((eta + 1)*(eta + 1)) + 10203.046875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6890.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 323.90625*((xi + 1)*(xi + 1)) - 9069.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2550.76171875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 682.5*(eta + 1)*(eta + 1)*(eta + 1) + 6125.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 218.75*(xi + 1)*(xi + 1)*(xi + 1) - 1722.65625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 191.953125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 49.21875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 33:
                        return RealGradient(0, 60.0*eta + 69.375*xi - 832.5*(0.5*eta + 0.5)*(2*xi + 2) + 4623.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 4200.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 1181.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1092.65625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 971.25*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 273.1640625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 124.375 - 157.5*(eta + 1)*(eta + 1) - 6068.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5512.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 192.65625*(xi + 1)*(xi + 1) + 5394.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1517.16796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 140.0*((eta + 1)*(eta + 1)*(eta + 1)) - 4900.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 175.0*((xi + 1)*(xi + 1)*(xi + 1)) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 39.375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 49.21875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 34:
                        return RealGradient(0, -135.0*eta - 28.125*xi + 337.5*(0.5*eta + 0.5)*(2*xi + 2) - 303.75*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 442.96875*(2*xi + 2)*(eta + 1)*(eta + 1) + 393.75*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 110.7421875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 151.875 + 354.375*((eta + 1)*(eta + 1)) + 398.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 12.65625*((xi + 1)*(xi + 1)) - 354.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 99.66796875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 315.0*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 35:
                        return RealGradient(0, 67.5*eta + 22.5*xi - 270.0*(0.5*eta + 0.5)*(2*xi + 2) + 303.75*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 354.375*(2*xi + 2)*((eta + 1)*(eta + 1)) - 315.0*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 84.375 - 177.1875*(eta + 1)*(eta + 1) - 398.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 12.65625*(xi + 1)*(xi + 1) + 354.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 99.66796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 157.5*((eta + 1)*(eta + 1)*(eta + 1)) - 44.296875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 36:
                        return RealGradient(-72.0*eta + 60.0*(0.5*eta + 0.5)*(2*xi + 2) - 123.75*(2*xi + 2)*(eta + 1)*(eta + 1) + 159.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 82.03125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 72.0 + 297.0*((eta + 1)*(eta + 1)) - 382.5*(eta + 1)*(eta + 1)*(eta + 1) + 196.875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 35.4375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 37:
                        return RealGradient(-48.0*eta + 60.0*(0.5*eta + 0.5)*(2*xi + 2) - 123.75*(2*xi + 2)*(eta + 1)*(eta + 1) + 159.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 82.03125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 48.0 + 198.0*((eta + 1)*(eta + 1)) - 255.0*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 23.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 38:
                        return RealGradient(30.0*eta - 30.0*(0.5*eta + 0.5)*(2*xi + 2) + 61.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 79.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 41.015625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 7.3828125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 30.0 - 123.75*(eta + 1)*(eta + 1) + 159.375*((eta + 1)*(eta + 1)*(eta + 1)) - 82.03125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 39:
                        return RealGradient(0, -72.0*eta - 297.0*xi + 594.0*(0.5*eta + 0.5)*(2*xi + 2) - 2295.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1575.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 354.375*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 123.75*(2*xi + 2)*(eta + 1)*(eta + 1) - 333.0 + 30.0*((eta + 1)*(eta + 1)) + 478.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 328.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 573.75*((xi + 1)*(xi + 1)) - 393.75*(xi + 1)*(xi + 1)*(xi + 1) + 88.59375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 40:
                        return RealGradient(0, -48.0*eta - 99.0*xi + 396.0*(0.5*eta + 0.5)*(2*xi + 2) - 1530.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 1050.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 236.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 123.75*(2*xi + 2)*(eta + 1)*(eta + 1) - 135.0 + 30.0*((eta + 1)*(eta + 1)) + 478.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 328.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 191.25*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1) + 29.53125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 41:
                        return RealGradient(0, 30.0*eta + 49.5*xi - 247.5*(0.5*eta + 0.5)*(2*xi + 2) + 956.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 656.25*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 147.65625*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 61.875*(2*xi + 2)*((eta + 1)*(eta + 1)) + 73.5 - 15.0*(eta + 1)*(eta + 1) - 239.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 164.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 36.9140625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 95.625*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 42:
                        return RealGradient(0, 9.0*eta + 108.0*xi - 216.0*(0.5*eta + 0.5)*(2*xi + 2) + 1350.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 1260.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 354.375*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 45.0*(2*xi + 2)*((eta + 1)*(eta + 1)) + 112.5 - 3.75*(eta + 1)*(eta + 1) - 281.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 73.828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 337.5*(xi + 1)*(xi + 1) + 315.0*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 43:
                        return RealGradient(0, 6.0*eta + 36.0*xi - 144.0*(0.5*eta + 0.5)*(2*xi + 2) + 900.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) - 840.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1) + 236.25*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 45.0*(2*xi + 2)*((eta + 1)*(eta + 1)) + 40.5 - 3.75*(eta + 1)*(eta + 1) - 281.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 262.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 73.828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 112.5*(xi + 1)*(xi + 1) + 105.0*((xi + 1)*(xi + 1)*(xi + 1)) - 29.53125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 44:
                        return RealGradient(0, -3.75*eta - 18.0*xi + 90.0*(0.5*eta + 0.5)*(2*xi + 2) - 562.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) + 525.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)*(xi + 1)) - 147.65625*(0.5*eta + 0.5)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 22.5*(2*xi + 2)*(eta + 1)*(eta + 1) - 21.0 + 1.875*((eta + 1)*(eta + 1)) + 140.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 131.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 36.9140625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 56.25*((xi + 1)*(xi + 1)) - 52.5*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 45:
                        return RealGradient(9.0*eta - 7.5*(0.5*eta + 0.5)*(2*xi + 2) + 45.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 93.75*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 9.0 - 108.0*(eta + 1)*(eta + 1) + 225.0*((eta + 1)*(eta + 1)*(eta + 1)) - 157.5*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 35.4375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 46:
                        return RealGradient(6.0*eta - 7.5*(0.5*eta + 0.5)*(2*xi + 2) + 45.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 93.75*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 6.0 - 72.0*(eta + 1)*(eta + 1) + 150.0*((eta + 1)*(eta + 1)*(eta + 1)) - 105.0*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 23.625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 47:
                        return RealGradient(-3.75*eta + 3.75*(0.5*eta + 0.5)*(2*xi + 2) - 22.5*(2*xi + 2)*(eta + 1)*(eta + 1) + 46.875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 32.8125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 7.3828125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 3.75 + 45.0*((eta + 1)*(eta + 1)) - 93.75*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 48:
                        return RealGradient(0, -36.0*eta - 40.5*xi + 81.0*(0.5*eta + 0.5)*(2*xi + 2) - 67.5*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 16.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 58.5 + 15.0*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 16.875*((xi + 1)*(xi + 1)));
                      case 49:
                        return RealGradient(0, 9.0*eta + 27.0*xi - 54.0*(0.5*eta + 0.5)*(2*xi + 2) + 67.5*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 11.25*(2*xi + 2)*((eta + 1)*(eta + 1)) + 31.5 - 3.75*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 16.875*(xi + 1)*(xi + 1));
                      case 50:
                        return RealGradient(36.0*eta - 30.0*(0.5*eta + 0.5)*(2*xi + 2) + 16.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 36.0 - 40.5*(eta + 1)*(eta + 1) + 11.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 51:
                        return RealGradient(-9.0*eta + 7.5*(0.5*eta + 0.5)*(2*xi + 2) - 11.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 4.6875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 9.0 + 27.0*((eta + 1)*(eta + 1)) - 11.25*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 52:
                        return RealGradient(24.0*eta - 30.0*(0.5*eta + 0.5)*(2*xi + 2) + 16.875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 4.6875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 24.0 - 27.0*(eta + 1)*(eta + 1) + 7.5*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 53:
                        return RealGradient(-6.0*eta + 7.5*(0.5*eta + 0.5)*(2*xi + 2) - 11.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 4.6875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 6.0 + 18.0*((eta + 1)*(eta + 1)) - 7.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 54:
                        return RealGradient(0, -24.0*eta - 13.5*xi + 54.0*(0.5*eta + 0.5)*(2*xi + 2) - 45.0*(0.5*eta + 0.5)*(xi + 1)*(xi + 1) - 16.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 31.5 + 15.0*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 5.625*((xi + 1)*(xi + 1)));
                      case 55:
                        return RealGradient(0, 6.0*eta + 9.0*xi - 36.0*(0.5*eta + 0.5)*(2*xi + 2) + 45.0*(0.5*eta + 0.5)*((xi + 1)*(xi + 1)) + 11.25*(2*xi + 2)*((eta + 1)*(eta + 1)) + 13.5 - 3.75*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 5.625*(xi + 1)*(xi + 1));
                      case 56:
                        return RealGradient(0, 0);
                      case 57:
                        return RealGradient(0, -4.5*xi - 2.5 + 1.875*((xi + 1)*(xi + 1)));
                      case 58:
                        return RealGradient(0, 3.0*xi + 2.5 - 1.875*(xi + 1)*(xi + 1));
                      case 59:
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
                        return sign * RealGradient(937.5*eta + 937.5*xi - 5625.0*(2*eta + 2)*(0.5*xi + 0.5) + 4921.875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 3281.25*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 738.28125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 19687.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 13125.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 2953.125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1718.75 - 1640.625*(eta + 1)*(eta + 1) - 17226.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 11484.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2583.984375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1640.625*(xi + 1)*(xi + 1) + 11484.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1093.75*((eta + 1)*(eta + 1)*(eta + 1)) - 7656.25*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1093.75*((xi + 1)*(xi + 1)*(xi + 1)) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 246.09375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 246.09375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 1:
                        return sign * RealGradient(-142.08984375*eta - 279.78515625*xi + 1678.7109375*(2*eta + 2)*(0.5*xi + 0.5) - 1545.7763671875*(2*eta + 2)*(xi + 1)*(xi + 1) + 997.1923828125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 213.409423828125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5875.48828125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 3916.9921875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 881.3232421875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 398.193359375 + 248.6572265625*((eta + 1)*(eta + 1)) + 5410.21728515625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3490.17333984375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 746.932983398438*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 515.2587890625*((xi + 1)*(xi + 1)) - 3606.8115234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 811.532592773438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 165.771484375*(eta + 1)*(eta + 1)*(eta + 1) + 2326.7822265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 497.955322265625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 332.3974609375*(xi + 1)*(xi + 1)*(xi + 1) - 523.526000976563*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 37.298583984375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 112.039947509766*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 71.136474609375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 2:
                        return sign * RealGradient(70.3125*eta + 164.0625*xi - 984.375*(2*eta + 2)*(0.5*xi + 0.5) + 1353.515625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1107.421875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 276.85546875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 3445.3125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 2296.875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 516.796875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 222.65625 - 123.046875*(eta + 1)*(eta + 1) - 4737.3046875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3875.9765625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 968.994140625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 451.171875*(xi + 1)*(xi + 1) + 3158.203125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 710.595703125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 82.03125*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 645.99609375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 369.140625*((xi + 1)*(xi + 1)*(xi + 1)) + 581.396484375*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 18.45703125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 145.34912109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 92.28515625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 3:
                        return sign * RealGradient(-24.90234375*eta - 68.84765625*xi + 413.0859375*(2*eta + 2)*(0.5*xi + 0.5) - 684.4482421875*(2*eta + 2)*(xi + 1)*(xi + 1) + 710.0830078125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 213.409423828125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1445.80078125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 963.8671875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 216.8701171875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 89.599609375 + 43.5791015625*((eta + 1)*(eta + 1)) + 2395.56884765625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2485.29052734375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 746.932983398438*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 228.1494140625*((xi + 1)*(xi + 1)) - 1597.0458984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 359.335327148438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 29.052734375*(eta + 1)*(eta + 1)*(eta + 1) + 1656.8603515625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 497.955322265625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 236.6943359375*(xi + 1)*(xi + 1)*(xi + 1) - 372.793579101563*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 6.536865234375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 112.039947509766*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 71.136474609375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 4:
                        return sign * RealGradient(187.5*eta + 375.0*xi - 2250.0*(2*eta + 2)*(0.5*xi + 0.5) + 2953.125*(2*eta + 2)*((xi + 1)*(xi + 1)) - 2625.0*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 738.28125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 7875.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 5250.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1181.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 531.25 - 328.125*(eta + 1)*(eta + 1) - 10335.9375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 9187.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2583.984375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 984.375*(xi + 1)*(xi + 1) + 6890.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 218.75*((eta + 1)*(eta + 1)*(eta + 1)) - 6125.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 875.0*((xi + 1)*(xi + 1)*(xi + 1)) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 246.09375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 5:
                        return sign * RealGradient(0, -187.5*xi + 656.25*(2*eta + 2)*(0.5*xi + 0.5) - 1968.75*(2*eta + 2)*(xi + 1)*(xi + 1) + 3445.3125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 2296.875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1312.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 393.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 187.5 + 3937.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6890.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4593.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1033.59375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1125.0*((xi + 1)*(xi + 1)) - 1181.25*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1968.75*(xi + 1)*(xi + 1)*(xi + 1) + 1312.5*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 295.3125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 6:
                        return sign * RealGradient(0, 55.95703125*xi - 206.103515625*(2*eta + 2)*(0.5*xi + 0.5) + 618.310546875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1082.04345703125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 721.3623046875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 162.306518554688*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 398.876953125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 113.818359375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 55.95703125 - 1196.630859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2094.10400390625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1396.0693359375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 314.115600585938*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 335.7421875*(xi + 1)*(xi + 1) + 341.455078125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 597.54638671875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 398.3642578125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 89.6319580078125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 587.548828125*((xi + 1)*(xi + 1)*(xi + 1)) - 391.69921875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 88.13232421875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 7:
                        return sign * RealGradient(0, -32.8125*xi + 180.46875*(2*eta + 2)*(0.5*xi + 0.5) - 541.40625*(2*eta + 2)*(xi + 1)*(xi + 1) + 947.4609375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 631.640625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 142.119140625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 442.96875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 147.65625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 32.8125 + 1328.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2325.5859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 348.837890625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1)) - 442.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 775.1953125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 116.279296875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 344.53125*(xi + 1)*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 51.6796875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 8:
                        return sign * RealGradient(0, 13.76953125*xi - 91.259765625*(2*eta + 2)*(0.5*xi + 0.5) + 273.779296875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 479.11376953125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 319.4091796875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 71.8670654296875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 284.033203125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 113.818359375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 13.76953125 - 852.099609375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1491.17431640625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 994.1162109375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 223.676147460938*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 82.6171875*(xi + 1)*(xi + 1) + 341.455078125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 597.54638671875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 398.3642578125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 89.6319580078125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 144.580078125*((xi + 1)*(xi + 1)*(xi + 1)) - 96.38671875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 21.68701171875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 9:
                        return sign * RealGradient(0, -75.0*xi + 393.75*(2*eta + 2)*(0.5*xi + 0.5) - 1181.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 2067.1875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1050.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 393.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 75.0 + 3150.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 5512.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3675.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 826.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 450.0*((xi + 1)*(xi + 1)) - 1181.25*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 787.5*(xi + 1)*(xi + 1)*(xi + 1) + 525.0*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 118.125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 10:
                        return sign * RealGradient(75.0*eta + 75.0*xi - 900.0*(2*eta + 2)*(0.5*xi + 0.5) + 1181.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1050.0*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 295.3125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 4725.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 4200.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1181.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 143.75 - 196.875*(eta + 1)*(eta + 1) - 6201.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5512.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 196.875*(xi + 1)*(xi + 1) + 5512.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 175.0*((eta + 1)*(eta + 1)*(eta + 1)) - 4900.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 175.0*((xi + 1)*(xi + 1)*(xi + 1)) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 49.21875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 11:
                        return sign * RealGradient(-9.9609375*eta - 13.76953125*xi + 165.234375*(2*eta + 2)*(0.5*xi + 0.5) - 273.779296875*(2*eta + 2)*(xi + 1)*(xi + 1) + 284.033203125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 85.36376953125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 867.48046875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 771.09375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 216.8701171875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 22.900390625 + 26.1474609375*((eta + 1)*(eta + 1)) + 1437.34130859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1491.17431640625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 448.159790039063*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 45.6298828125*((xi + 1)*(xi + 1)) - 1277.63671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 359.335327148438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 23.2421875*(eta + 1)*(eta + 1)*(eta + 1) + 1325.48828125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 398.3642578125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 47.3388671875*(xi + 1)*(xi + 1)*(xi + 1) - 372.793579101563*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 6.536865234375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 112.039947509766*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 14.227294921875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 12:
                        return sign * RealGradient(28.125*eta + 32.8125*xi - 393.75*(2*eta + 2)*(0.5*xi + 0.5) + 541.40625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 442.96875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 110.7421875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2067.1875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1837.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 516.796875*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 58.59375 - 73.828125*(eta + 1)*(eta + 1) - 2842.3828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2325.5859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 581.396484375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 90.234375*(xi + 1)*(xi + 1) + 2526.5625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 710.595703125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1)) - 2067.1875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 73.828125*((xi + 1)*(xi + 1)*(xi + 1)) + 581.396484375*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 18.45703125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 145.34912109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 18.45703125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 13:
                        return sign * RealGradient(-56.8359375*eta - 55.95703125*xi + 671.484375*(2*eta + 2)*(0.5*xi + 0.5) - 618.310546875*(2*eta + 2)*(xi + 1)*(xi + 1) + 398.876953125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 85.36376953125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3525.29296875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 3133.59375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 881.3232421875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 108.056640625 + 149.1943359375*((eta + 1)*(eta + 1)) + 3246.13037109375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2094.10400390625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 448.159790039063*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 103.0517578125*((xi + 1)*(xi + 1)) - 2885.44921875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 811.532592773438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 132.6171875*(eta + 1)*(eta + 1)*(eta + 1) + 1861.42578125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 398.3642578125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 66.4794921875*(xi + 1)*(xi + 1)*(xi + 1) - 523.526000976563*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 37.298583984375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 112.039947509766*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 14.227294921875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 14:
                        return sign * RealGradient(375.0*eta + 187.5*xi - 2250.0*(2*eta + 2)*(0.5*xi + 0.5) + 1968.75*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1312.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 295.3125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 11812.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 10500.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 2953.125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 531.25 - 984.375*(eta + 1)*(eta + 1) - 10335.9375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6890.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 328.125*(xi + 1)*(xi + 1) + 9187.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 875.0*((eta + 1)*(eta + 1)*(eta + 1)) - 6125.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 218.75*((xi + 1)*(xi + 1)*(xi + 1)) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 246.09375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 49.21875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 15:
                        return sign * RealGradient(0, -157.5*eta - 375.0*xi + 1968.75*(2*eta + 2)*(0.5*xi + 0.5) - 2953.125*(2*eta + 2)*(xi + 1)*(xi + 1) + 3445.3125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 5250.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1968.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 502.5 + 210.0*((eta + 1)*(eta + 1)) + 7875.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 9187.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4593.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 826.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1125.0*((xi + 1)*(xi + 1)) - 2953.125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 78.75*(eta + 1)*(eta + 1)*(eta + 1) + 3445.3125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1312.5*(xi + 1)*(xi + 1)*(xi + 1) + 656.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 118.125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 16:
                        return sign * RealGradient(0, 36.50390625*eta + 68.84765625*xi - 456.298828125*(2*eta + 2)*(0.5*xi + 0.5) + 684.4482421875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 798.52294921875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 399.261474609375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 71.8670654296875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1420.166015625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 569.091796875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 99.84375 - 56.806640625*(eta + 1)*(eta + 1) - 2130.2490234375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2485.29052734375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1242.64526367188*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 223.676147460938*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 206.54296875*(xi + 1)*(xi + 1) + 853.6376953125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 22.763671875*((eta + 1)*(eta + 1)*(eta + 1)) - 995.91064453125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 497.955322265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 89.6319580078125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 240.966796875*((xi + 1)*(xi + 1)*(xi + 1)) - 120.4833984375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 21.68701171875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 17:
                        return sign * RealGradient(0, -72.1875*eta - 164.0625*xi + 902.34375*(2*eta + 2)*(0.5*xi + 0.5) - 1353.515625*(2*eta + 2)*(xi + 1)*(xi + 1) + 1579.1015625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 789.55078125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 142.119140625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2214.84375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 738.28125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 223.125 + 88.59375*((eta + 1)*(eta + 1)) + 3322.265625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3875.9765625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1937.98828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 348.837890625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 492.1875*((xi + 1)*(xi + 1)) - 1107.421875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 29.53125*(eta + 1)*(eta + 1)*(eta + 1) + 1291.9921875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 645.99609375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 116.279296875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 574.21875*(xi + 1)*(xi + 1)*(xi + 1) + 287.109375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 51.6796875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 18:
                        return sign * RealGradient(0, 82.44140625*eta + 279.78515625*xi - 1030.517578125*(2*eta + 2)*(0.5*xi + 0.5) + 1545.7763671875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1803.40576171875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 901.702880859375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 162.306518554688*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1994.384765625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 569.091796875*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 339.84375 - 79.775390625*(eta + 1)*(eta + 1) - 2991.5771484375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3490.17333984375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1745.08666992188*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 314.115600585938*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 839.35546875*(xi + 1)*(xi + 1) + 853.6376953125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 22.763671875*((eta + 1)*(eta + 1)*(eta + 1)) - 995.91064453125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 497.955322265625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 89.6319580078125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 979.248046875*((xi + 1)*(xi + 1)*(xi + 1)) - 489.6240234375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 88.13232421875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 19:
                        return sign * RealGradient(0, -262.5*eta - 937.5*xi + 3281.25*(2*eta + 2)*(0.5*xi + 0.5) - 4921.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 5742.1875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 2871.09375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 6562.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1968.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 1125.0 + 262.5*((eta + 1)*(eta + 1)) + 9843.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 11484.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5742.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1033.59375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2812.5*((xi + 1)*(xi + 1)) - 2953.125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 78.75*(eta + 1)*(eta + 1)*(eta + 1) + 3445.3125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 3281.25*(xi + 1)*(xi + 1)*(xi + 1) + 1640.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 295.3125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 20:
                        return RealGradient(0, -731.25*xi + 2559.375*(2*eta + 2)*(0.5*xi + 0.5) - 4626.5625*(2*eta + 2)*(xi + 1)*(xi + 1) + 5668.359375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 2871.09375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 516.796875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 5118.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1535.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 731.25 + 9253.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 11336.71875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5742.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1033.59375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2643.75*((xi + 1)*(xi + 1)) - 2775.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3401.015625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 3239.0625*(xi + 1)*(xi + 1)*(xi + 1) + 1640.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 295.3125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 21:
                        return RealGradient(0, 150.0*xi - 525.0*(2*eta + 2)*(0.5*xi + 0.5) + 1821.09375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 3371.484375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 2296.875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 516.796875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1050.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 315.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 150.0 - 3642.1875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6742.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4593.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1033.59375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1040.625*(xi + 1)*(xi + 1) + 1092.65625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 2022.890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 310.078125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1926.5625*((xi + 1)*(xi + 1)*(xi + 1)) - 1312.5*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 295.3125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 22:
                        return RealGradient(0, -337.5*xi + 1181.25*(2*eta + 2)*(0.5*xi + 0.5) - 738.28125*(2*eta + 2)*(xi + 1)*(xi + 1) + 221.484375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 2362.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 708.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 337.5 + 1476.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 442.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 421.875*((xi + 1)*(xi + 1)) - 442.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 132.890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 126.5625*(xi + 1)*(xi + 1)*(xi + 1));
                      case 23:
                        return RealGradient(0, 168.75*xi - 590.625*(2*eta + 2)*(0.5*xi + 0.5) + 590.625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 221.484375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1181.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 354.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 168.75 - 1181.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 442.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 337.5*(xi + 1)*(xi + 1) + 354.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 132.890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 126.5625*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 24:
                        return RealGradient(881.25*eta + 731.25*xi - 5287.5*(2*eta + 2)*(0.5*xi + 0.5) + 4626.5625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 3084.375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 693.984375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 19434.375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 13125.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 2953.125*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1490.625 - 1619.53125*(eta + 1)*(eta + 1) - 17005.078125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 11336.71875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2550.76171875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1279.6875*(xi + 1)*(xi + 1) + 11484.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1093.75*((eta + 1)*(eta + 1)*(eta + 1)) - 7656.25*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 853.125*((xi + 1)*(xi + 1)*(xi + 1)) + 1722.65625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 246.09375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 191.953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 25:
                        return RealGradient(-346.875*eta - 150.0*xi + 2081.25*(2*eta + 2)*(0.5*xi + 0.5) - 1821.09375*(2*eta + 2)*(xi + 1)*(xi + 1) + 1214.0625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 273.1640625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 11559.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 10500.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 2953.125*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 471.875 + 963.28125*((eta + 1)*(eta + 1)) + 10114.453125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6742.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1517.16796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 262.5*((xi + 1)*(xi + 1)) - 9187.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2583.984375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 875.0*(eta + 1)*(eta + 1)*(eta + 1) + 6125.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 175.0*(xi + 1)*(xi + 1)*(xi + 1) - 1722.65625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 246.09375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 39.375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 26:
                        return RealGradient(140.625*eta + 337.5*xi - 843.75*(2*eta + 2)*(0.5*xi + 0.5) + 738.28125*(2*eta + 2)*((xi + 1)*(xi + 1)) - 492.1875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 110.7421875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 759.375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 421.875 - 63.28125*(eta + 1)*(eta + 1) - 664.453125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 442.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 99.66796875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 590.625*(xi + 1)*(xi + 1) + 393.75*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 27:
                        return RealGradient(-112.5*eta - 168.75*xi + 675.0*(2*eta + 2)*(0.5*xi + 0.5) - 590.625*(2*eta + 2)*(xi + 1)*(xi + 1) + 393.75*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 759.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 253.125 + 63.28125*((eta + 1)*(eta + 1)) + 664.453125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 442.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 99.66796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 295.3125*((xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(xi + 1) + 44.296875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 28:
                        return RealGradient(176.25*eta + 292.5*xi - 2115.0*(2*eta + 2)*(0.5*xi + 0.5) + 2775.9375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 2467.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 693.984375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 7773.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 5250.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 1181.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 444.375 - 323.90625*(eta + 1)*(eta + 1) - 10203.046875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 9069.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2550.76171875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 767.8125*(xi + 1)*(xi + 1) + 6890.625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 218.75*((eta + 1)*(eta + 1)*(eta + 1)) - 6125.0*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1722.65625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 682.5*((xi + 1)*(xi + 1)*(xi + 1)) + 1378.125*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 49.21875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 387.59765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 191.953125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 29:
                        return RealGradient(-69.375*eta - 60.0*xi + 832.5*(2*eta + 2)*(0.5*xi + 0.5) - 1092.65625*(2*eta + 2)*(xi + 1)*(xi + 1) + 971.25*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 273.1640625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 4623.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 4200.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 1181.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 124.375 + 192.65625*((eta + 1)*(eta + 1)) + 6068.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 5394.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1517.16796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 157.5*((xi + 1)*(xi + 1)) - 5512.5*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1550.390625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 175.0*(eta + 1)*(eta + 1)*(eta + 1) + 4900.0*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1378.125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 140.0*(xi + 1)*(xi + 1)*(xi + 1) - 1378.125*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 49.21875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 387.59765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 39.375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 30:
                        return RealGradient(28.125*eta + 135.0*xi - 337.5*(2*eta + 2)*(0.5*xi + 0.5) + 442.96875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 393.75*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 110.7421875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 303.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 151.875 - 12.65625*(eta + 1)*(eta + 1) - 398.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 354.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 99.66796875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 354.375*(xi + 1)*(xi + 1) + 315.0*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 31:
                        return RealGradient(-22.5*eta - 67.5*xi + 270.0*(2*eta + 2)*(0.5*xi + 0.5) - 354.375*(2*eta + 2)*(xi + 1)*(xi + 1) + 315.0*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 88.59375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 303.75*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 84.375 + 12.65625*((eta + 1)*(eta + 1)) + 398.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 354.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 99.66796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 177.1875*((xi + 1)*(xi + 1)) - 157.5*(xi + 1)*(xi + 1)*(xi + 1) + 44.296875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 32:
                        return RealGradient(0, -292.5*xi + 1535.625*(2*eta + 2)*(0.5*xi + 0.5) - 2775.9375*(2*eta + 2)*(xi + 1)*(xi + 1) + 3401.015625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 4095.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1535.625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 292.5 + 7402.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 9069.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4593.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 826.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1057.5*((xi + 1)*(xi + 1)) - 2775.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3401.015625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1722.65625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 310.078125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1295.625*(xi + 1)*(xi + 1)*(xi + 1) + 656.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 118.125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 33:
                        return RealGradient(0, 60.0*xi - 315.0*(2*eta + 2)*(0.5*xi + 0.5) + 1092.65625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 2022.890625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 310.078125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 840.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 315.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 60.0 - 2913.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5394.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 3675.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 826.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 416.25*(xi + 1)*(xi + 1) + 1092.65625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 2022.890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1378.125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 310.078125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 770.625*((xi + 1)*(xi + 1)*(xi + 1)) - 525.0*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 118.125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 34:
                        return RealGradient(0, -135.0*xi + 708.75*(2*eta + 2)*(0.5*xi + 0.5) - 442.96875*(2*eta + 2)*(xi + 1)*(xi + 1) + 132.890625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 1890.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 708.75*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 135.0 + 1181.25*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 354.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 168.75*((xi + 1)*(xi + 1)) - 442.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 132.890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 50.625*(xi + 1)*(xi + 1)*(xi + 1));
                      case 35:
                        return RealGradient(0, 67.5*xi - 354.375*(2*eta + 2)*(0.5*xi + 0.5) + 354.375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 132.890625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 945.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 354.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 67.5 - 945.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 354.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 135.0*(xi + 1)*(xi + 1) + 354.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 132.890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 50.625*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 36:
                        return RealGradient(-297.0*eta - 72.0*xi + 594.0*(2*eta + 2)*(0.5*xi + 0.5) - 123.75*(2*eta + 2)*(xi + 1)*(xi + 1) - 2295.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1575.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 354.375*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 333.0 + 573.75*((eta + 1)*(eta + 1)) + 478.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 30.0*((xi + 1)*(xi + 1)) - 328.125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 73.828125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 393.75*(eta + 1)*(eta + 1)*(eta + 1) + 88.59375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 37:
                        return RealGradient(-99.0*eta - 48.0*xi + 396.0*(2*eta + 2)*(0.5*xi + 0.5) - 123.75*(2*eta + 2)*(xi + 1)*(xi + 1) - 1530.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 1050.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 236.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 135.0 + 191.25*((eta + 1)*(eta + 1)) + 478.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 30.0*((xi + 1)*(xi + 1)) - 328.125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 73.828125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 131.25*(eta + 1)*(eta + 1)*(eta + 1) + 29.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 38:
                        return RealGradient(49.5*eta + 30.0*xi - 247.5*(2*eta + 2)*(0.5*xi + 0.5) + 61.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 956.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 656.25*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 147.65625*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 73.5 - 95.625*(eta + 1)*(eta + 1) - 239.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1) + 164.0625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 36.9140625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 65.625*((eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 39:
                        return RealGradient(0, -72.0*xi + 60.0*(2*eta + 2)*(0.5*xi + 0.5) - 123.75*(2*eta + 2)*(xi + 1)*(xi + 1) + 159.375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 82.03125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 72.0 + 297.0*((xi + 1)*(xi + 1)) - 382.5*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 35.4375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 40:
                        return RealGradient(0, -48.0*xi + 60.0*(2*eta + 2)*(0.5*xi + 0.5) - 123.75*(2*eta + 2)*(xi + 1)*(xi + 1) + 159.375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 82.03125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 48.0 + 198.0*((xi + 1)*(xi + 1)) - 255.0*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 23.625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 41:
                        return RealGradient(0, 30.0*xi - 30.0*(2*eta + 2)*(0.5*xi + 0.5) + 61.875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 79.6875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 41.015625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 7.3828125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 30.0 - 123.75*(xi + 1)*(xi + 1) + 159.375*((xi + 1)*(xi + 1)*(xi + 1)) - 82.03125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 42:
                        return RealGradient(0, 9.0*xi - 7.5*(2*eta + 2)*(0.5*xi + 0.5) + 45.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 93.75*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9.0 - 108.0*(xi + 1)*(xi + 1) + 225.0*((xi + 1)*(xi + 1)*(xi + 1)) - 157.5*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 35.4375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 43:
                        return RealGradient(0, 6.0*xi - 7.5*(2*eta + 2)*(0.5*xi + 0.5) + 45.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 93.75*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 6.0 - 72.0*(xi + 1)*(xi + 1) + 150.0*((xi + 1)*(xi + 1)*(xi + 1)) - 105.0*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 23.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 44:
                        return RealGradient(0, -3.75*xi + 3.75*(2*eta + 2)*(0.5*xi + 0.5) - 22.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 46.875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 32.8125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7.3828125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 3.75 + 45.0*((xi + 1)*(xi + 1)) - 93.75*(xi + 1)*(xi + 1)*(xi + 1) + 65.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 45:
                        return RealGradient(108.0*eta + 9.0*xi - 216.0*(2*eta + 2)*(0.5*xi + 0.5) + 45.0*(2*eta + 2)*((xi + 1)*(xi + 1)) + 1350.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 1260.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 354.375*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 112.5 - 337.5*(eta + 1)*(eta + 1) - 281.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 3.75*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 73.828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 315.0*((eta + 1)*(eta + 1)*(eta + 1)) - 88.59375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 46:
                        return RealGradient(36.0*eta + 6.0*xi - 144.0*(2*eta + 2)*(0.5*xi + 0.5) + 45.0*(2*eta + 2)*((xi + 1)*(xi + 1)) + 900.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) - 840.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1) + 236.25*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 40.5 - 112.5*(eta + 1)*(eta + 1) - 281.25*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 3.75*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 73.828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 105.0*((eta + 1)*(eta + 1)*(eta + 1)) - 29.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 47:
                        return RealGradient(-18.0*eta - 3.75*xi + 90.0*(2*eta + 2)*(0.5*xi + 0.5) - 22.5*(2*eta + 2)*(xi + 1)*(xi + 1) - 562.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) + 525.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)*(eta + 1)) - 147.65625*(0.5*xi + 0.5)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 21.0 + 56.25*((eta + 1)*(eta + 1)) + 140.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 1.875*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 36.9140625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 52.5*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 48:
                        return RealGradient(0, -36.0*xi + 30.0*(2*eta + 2)*(0.5*xi + 0.5) - 16.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 36.0 + 40.5*((xi + 1)*(xi + 1)) - 11.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 49:
                        return RealGradient(0, 9.0*xi - 7.5*(2*eta + 2)*(0.5*xi + 0.5) + 11.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 4.6875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 9.0 - 27.0*(xi + 1)*(xi + 1) + 11.25*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 50:
                        return RealGradient(40.5*eta + 36.0*xi - 81.0*(2*eta + 2)*(0.5*xi + 0.5) + 16.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 67.5*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 58.5 - 16.875*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1), 0);
                      case 51:
                        return RealGradient(-27.0*eta - 9.0*xi + 54.0*(2*eta + 2)*(0.5*xi + 0.5) - 11.25*(2*eta + 2)*(xi + 1)*(xi + 1) - 67.5*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 31.5 + 16.875*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.75*((xi + 1)*(xi + 1)), 0);
                      case 52:
                        return RealGradient(13.5*eta + 24.0*xi - 54.0*(2*eta + 2)*(0.5*xi + 0.5) + 16.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 45.0*(0.5*xi + 0.5)*((eta + 1)*(eta + 1)) + 31.5 - 5.625*(eta + 1)*(eta + 1) - 14.0625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.0*(xi + 1)*(xi + 1), 0);
                      case 53:
                        return RealGradient(-9.0*eta - 6.0*xi + 36.0*(2*eta + 2)*(0.5*xi + 0.5) - 11.25*(2*eta + 2)*(xi + 1)*(xi + 1) - 45.0*(0.5*xi + 0.5)*(eta + 1)*(eta + 1) - 13.5 + 5.625*((eta + 1)*(eta + 1)) + 14.0625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 3.75*((xi + 1)*(xi + 1)), 0);
                      case 54:
                        return RealGradient(0, -24.0*xi + 30.0*(2*eta + 2)*(0.5*xi + 0.5) - 16.875*(2*eta + 2)*(xi + 1)*(xi + 1) + 4.6875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 24.0 + 27.0*((xi + 1)*(xi + 1)) - 7.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 55:
                        return RealGradient(0, 6.0*xi - 7.5*(2*eta + 2)*(0.5*xi + 0.5) + 11.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 4.6875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 6.0 - 18.0*(xi + 1)*(xi + 1) + 7.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 56:
                        return RealGradient(-4.5*eta - 2.5 + 1.875*((eta + 1)*(eta + 1)), 0);
                      case 57:
                        return RealGradient(0, 0);
                      case 58:
                        return RealGradient(0, 0);
                      case 59:
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
                        return sign * RealGradient(-12600.0*eta*xi + 2800.0*eta + 15120.0*eta*(xi*xi) - 4200.0*eta*xi*xi*xi + 2100.0*xi + 22680.0*xi*(eta*eta) - 12600.0*xi*eta*eta*eta - 300.0 - 8400.0*eta*eta - 12600.0*eta*eta*xi*xi - 4200.0*xi*xi + 10080.0*(eta*eta*eta) + 2520.0*(xi*xi*xi) - 4200.0*eta*eta*eta*eta, 8400.0*eta*xi - 700.0*eta - 22680.0*eta*xi*xi + 16800.0*eta*(xi*xi*xi) - 1400.0*xi - 15120.0*xi*eta*eta + 8400.0*xi*(eta*eta*eta) + 75.0 + 2100.0*(eta*eta) + 18900.0*(eta*eta)*(xi*xi) + 6300.0*(xi*xi) - 2520.0*eta*eta*eta - 10080.0*xi*xi*xi + 1050.0*(eta*eta*eta*eta) + 5250.0*(xi*xi*xi*xi));
                      case 1:
                        return sign * RealGradient(4255.78125*eta*xi - 759.0625*eta - 5118.75*eta*xi*xi + 1214.0625*eta*(xi*xi*xi) - 659.53125*xi - 7530.46875*xi*eta*eta + 3986.71875*xi*(eta*eta*eta) + 89.53125 + 1811.25*(eta*eta) + 4577.34375*(eta*eta)*(xi*xi) + 1276.40625*(xi*xi) - 1594.6875*eta*eta*eta - 728.4375*xi*xi*xi + 442.96875*(eta*eta*eta*eta), -2460.9375*eta*xi + 95.15625*eta + 7884.84375*eta*(xi*xi) - 6103.125*eta*xi*xi*xi + 428.75*xi + 3189.375*xi*(eta*eta) - 885.9375*xi*eta*eta*eta - 18.7109375 + 44.296875*(eta*eta) - 5980.078125*eta*eta*xi*xi - 1957.265625*xi*xi - 442.96875*eta*eta*eta + 3038.4375*(xi*xi*xi) + 332.2265625*(eta*eta*eta*eta) - 1517.578125*xi*xi*xi*xi);
                      case 2:
                        return sign * RealGradient(-1627.5*eta*xi + 3780.0*eta*(xi*xi) - 1575.0*eta*xi*xi*xi + 577.5*xi + 472.5*xi*(eta*eta) + 787.5*xi*(eta*eta*eta) - 52.5 + 840.0*(eta*eta) - 2362.5*eta*eta*xi*xi - 1417.5*xi*xi - 1575.0*eta*eta*eta + 945.0*(xi*xi*xi) + 787.5*(eta*eta*eta*eta), 105.0*eta*xi + 17.5*eta - 2362.5*eta*xi*xi + 3150.0*eta*(xi*xi*xi) - 315.0*xi + 1890.0*xi*(eta*eta) - 1575.0*xi*eta*eta*eta + 10.625 - 131.25*eta*eta - 1181.25*eta*eta*xi*xi + 1863.75*(xi*xi) + 52.5*(eta*eta*eta) - 3465.0*xi*xi*xi + 65.625*(eta*eta*eta*eta) + 1968.75*(xi*xi*xi*xi));
                      case 3:
                        return sign * RealGradient(-521.71875*eta*xi + 115.9375*eta - 708.75*eta*xi*xi + 1214.0625*eta*(xi*xi*xi) - 292.03125*xi + 2392.03125*xi*(eta*eta) - 1525.78125*xi*eta*eta*eta + 22.03125 - 288.75*eta*eta - 935.15625*eta*eta*xi*xi + 908.90625*(xi*xi) - 19.6875*eta*eta*eta - 728.4375*xi*xi*xi + 180.46875*(eta*eta*eta*eta), 584.0625*eta*xi - 27.34375*eta - 2037.65625*eta*xi*xi + 1246.875*eta*(xi*xi*xi) + 148.75*xi - 590.625*xi*eta*eta - 360.9375*xi*eta*eta*eta - 4.3359375 + 18.046875*(eta*eta) + 2288.671875*(eta*eta)*(xi*xi) - 1038.515625*xi*xi + 29.53125*(eta*eta*eta) + 2303.4375*(xi*xi*xi) + 4.1015625*(eta*eta*eta*eta) - 1517.578125*xi*xi*xi*xi);
                      case 4:
                        return sign * RealGradient(-1680.0*eta*xi + 140.0*eta + 5040.0*eta*(xi*xi) - 4200.0*eta*xi*xi*xi + 1260.0*xi - 120.0 - 3360.0*xi*xi + 2520.0*(xi*xi*xi), -700.0*xi + 25.0 + 4200.0*(xi*xi) - 8400.0*xi*xi*xi + 5250.0*(xi*xi*xi*xi));
                      case 5:
                        return sign * RealGradient(-1680.0*eta*xi + 140.0*eta + 5040.0*eta*(xi*xi) - 4200.0*eta*xi*xi*xi, -700.0*xi + 25.0 + 4200.0*(xi*xi) - 8400.0*xi*xi*xi + 5250.0*(xi*xi*xi*xi));
                      case 6:
                        return sign * RealGradient(-1417.5*eta*xi + 183.75*eta + 2657.8125*eta*(xi*xi) - 1328.90625*eta*xi*xi*xi + 4784.0625*xi*(eta*eta) - 2657.8125*xi*eta*eta*eta - 787.5*eta*eta - 5315.625*eta*eta*xi*xi + 767.8125*(eta*eta*eta) - 196.875*eta*eta*eta*eta, 2205.0*eta*xi - 113.75*eta - 7973.4375*eta*xi*xi + 7087.5*eta*(xi*xi*xi) - 446.25*xi - 1949.0625*xi*eta*eta + 393.75*xi*(eta*eta*eta) + 20.0 + 131.25*(eta*eta) + 3986.71875*(eta*eta)*(xi*xi) + 2126.25*(xi*xi) - 45.9375*eta*eta*eta - 3366.5625*xi*xi*xi + 4.1015625*(eta*eta*eta*eta) + 1661.1328125*(xi*xi*xi*xi));
                      case 7:
                        return sign * RealGradient(-840.0*eta*xi + 175.0*eta + 945.0*eta*(xi*xi) - 262.5*eta*xi*xi*xi + 4725.0*xi*(eta*eta) - 4725.0*xi*eta*eta*eta - 1260.0*eta*eta - 3150.0*eta*eta*xi*xi + 2205.0*(eta*eta*eta) - 1050.0*eta*eta*eta*eta, 2520.0*eta*xi - 175.0*eta - 6615.0*eta*xi*xi + 4200.0*eta*(xi*xi*xi) - 245.0*xi - 4725.0*xi*eta*eta + 2100.0*xi*(eta*eta*eta) + 15.0 + 420.0*(eta*eta) + 7087.5*(eta*eta)*(xi*xi) + 840.0*(xi*xi) - 315.0*eta*eta*eta - 945.0*xi*xi*xi + 65.625*(eta*eta*eta*eta) + 328.125*(xi*xi*xi*xi));
                      case 8:
                        return sign * RealGradient(-262.5*eta*xi + 113.75*eta + 137.8125*eta*(xi*xi) - 16.40625*eta*xi*xi*xi + 1949.0625*xi*(eta*eta) - 2657.8125*xi*eta*eta*eta - 1102.5*eta*eta - 590.625*eta*eta*xi*xi + 2657.8125*(eta*eta*eta) - 1771.875*eta*eta*eta*eta, 1575.0*eta*xi - 183.75*eta - 2303.4375*eta*xi*xi + 787.5*eta*(xi*xi*xi) - 96.25*xi - 4784.0625*xi*eta*eta + 3543.75*xi*(eta*eta*eta) + 10.0 + 708.75*(eta*eta) + 3986.71875*(eta*eta)*(xi*xi) + 183.75*(xi*xi) - 885.9375*eta*eta*eta - 111.5625*xi*xi*xi + 332.2265625*(eta*eta*eta*eta) + 20.5078125*(xi*xi*xi*xi));
                      case 9:
                        return sign * RealGradient(0, -140.0*eta + 5.0 + 840.0*(eta*eta) - 1680.0*eta*eta*eta + 1050.0*(eta*eta*eta*eta));
                      case 10:
                        return sign * RealGradient(0, -140.0*eta + 5.0 + 840.0*(eta*eta) - 1680.0*eta*eta*eta + 1050.0*(eta*eta*eta*eta));
                      case 11:
                        return sign * RealGradient(-36.09375*eta*xi + 27.34375*eta - 88.59375*eta*xi*xi - 16.40625*eta*xi*xi*xi + 590.625*xi*(eta*eta) - 1525.78125*xi*eta*eta*eta - 292.03125*eta*eta + 541.40625*(eta*eta)*(xi*xi) + 679.21875*(eta*eta*eta) - 311.71875*eta*eta*eta*eta, 577.5*eta*xi - 115.9375*eta + 59.0625*eta*(xi*xi) - 721.875*eta*xi*xi*xi - 18.59375*xi - 2392.03125*xi*eta*eta + 623.4375*xi*(eta*eta*eta) + 6.4453125 + 260.859375*(eta*eta) + 2288.671875*(eta*eta)*(xi*xi) - 27.890625*xi*xi + 236.25*(eta*eta*eta) + 29.53125*(xi*xi*xi) - 303.515625*eta*eta*eta*eta + 20.5078125*(xi*xi*xi*xi));
                      case 12:
                        return sign * RealGradient(262.5*eta*xi - 17.5*eta - 157.5*eta*xi*xi - 262.5*eta*xi*xi*xi - 1890.0*xi*eta*eta + 787.5*xi*(eta*eta*eta) - 52.5*eta*eta + 2362.5*(eta*eta)*(xi*xi) + 787.5*(eta*eta*eta) - 787.5*eta*eta*eta*eta, -1680.0*eta*xi + 4725.0*eta*(xi*xi) - 3150.0*eta*xi*xi*xi + 87.5*xi - 472.5*xi*eta*eta + 1575.0*xi*(eta*eta*eta) - 6.875 + 813.75*(eta*eta) - 1181.25*eta*eta*xi*xi - 26.25*xi*xi - 1260.0*eta*eta*eta - 367.5*xi*xi*xi + 393.75*(eta*eta*eta*eta) + 328.125*(xi*xi*xi*xi));
                      case 13:
                        return sign * RealGradient(-88.59375*eta*xi - 95.15625*eta + 1328.90625*eta*(xi*xi) - 1328.90625*eta*xi*xi*xi - 3189.375*xi*eta*eta + 3986.71875*xi*(eta*eta*eta) + 1230.46875*(eta*eta) + 1328.90625*(eta*eta)*(xi*xi) - 2628.28125*eta*eta*eta + 1525.78125*(eta*eta*eta*eta), -3622.5*eta*xi + 759.0625*eta + 4784.0625*eta*(xi*xi) - 1771.875*eta*xi*xi*xi - 351.09375*xi + 7530.46875*xi*(eta*eta) - 3051.5625*xi*eta*eta*eta - 5.4296875 - 2127.890625*eta*eta - 5980.078125*eta*eta*xi*xi + 1993.359375*(xi*xi) + 1706.25*(eta*eta*eta) - 3277.96875*xi*xi*xi - 303.515625*eta*eta*eta*eta + 1661.1328125*(xi*xi*xi*xi));
                      case 14:
                        return sign * RealGradient(-4200.0*eta*xi + 700.0*eta + 7560.0*eta*(xi*xi) - 4200.0*eta*xi*xi*xi + 15120.0*xi*(eta*eta) - 12600.0*xi*eta*eta*eta - 4200.0*eta*eta - 12600.0*eta*eta*xi*xi + 7560.0*(eta*eta*eta) - 4200.0*eta*eta*eta*eta, 16800.0*eta*xi - 2800.0*eta - 30240.0*eta*xi*xi + 16800.0*eta*(xi*xi*xi) - 3500.0*xi - 22680.0*xi*eta*eta + 8400.0*xi*(eta*eta*eta) + 375.0 + 6300.0*(eta*eta) + 18900.0*(eta*eta)*(xi*xi) + 10500.0*(xi*xi) - 5040.0*eta*eta*eta - 12600.0*xi*xi*xi + 1050.0*(eta*eta*eta*eta) + 5250.0*(xi*xi*xi*xi));
                      case 15:
                        return RealGradient(70560.0*eta*xi - 19600.0*eta - 70560.0*eta*xi*xi + 16800.0*eta*(xi*xi*xi) - 211680.0*xi*eta*eta + 151200.0*xi*(eta*eta*eta) + 94080.0*(eta*eta) + 100800.0*(eta*eta)*(xi*xi) - 141120.0*eta*eta*eta + 67200.0*(eta*eta*eta*eta), -94080.0*eta*xi + 9800.0*eta + 211680.0*eta*(xi*xi) - 134400.0*eta*xi*xi*xi + 9800.0*xi + 211680.0*xi*(eta*eta) - 134400.0*xi*eta*eta*eta - 700.0 - 35280.0*eta*eta - 226800.0*eta*eta*xi*xi - 35280.0*xi*xi + 47040.0*(eta*eta*eta) + 47040.0*(xi*xi*xi) - 21000.0*eta*eta*eta*eta - 21000.0*xi*xi*xi*xi);
                      case 16:
                        return RealGradient(-70560.0*eta*xi + 9800.0*eta + 141120.0*eta*(xi*xi) - 84000.0*eta*xi*xi*xi + 211680.0*xi*(eta*eta) - 151200.0*xi*eta*eta*eta - 47040.0*eta*eta - 201600.0*eta*eta*xi*xi + 70560.0*(eta*eta*eta) - 33600.0*eta*eta*eta*eta, 188160.0*eta*xi - 19600.0*eta - 423360.0*eta*xi*xi + 268800.0*eta*(xi*xi*xi) - 49000.0*xi - 211680.0*xi*eta*eta + 67200.0*xi*(eta*eta*eta) + 3500.0 + 35280.0*(eta*eta) + 226800.0*(eta*eta)*(xi*xi) + 176400.0*(xi*xi) - 23520.0*eta*eta*eta - 235200.0*xi*xi*xi + 4200.0*(eta*eta*eta*eta) + 105000.0*(xi*xi*xi*xi));
                      case 17:
                        return RealGradient(-60480.0*eta*xi + 6440.0*eta + 131040.0*eta*(xi*xi) - 67200.0*eta*xi*xi*xi + 60480.0*xi*(eta*eta) - 6720.0*eta*eta - 100800.0*eta*eta*xi*xi, 26880.0*eta*xi - 1120.0*eta - 120960.0*eta*xi*xi + 134400.0*eta*(xi*xi*xi) - 14560.0*xi + 560.0 + 80640.0*(xi*xi) - 147840.0*xi*xi*xi + 84000.0*(xi*xi*xi*xi));
                      case 18:
                        return RealGradient(40320.0*eta*xi - 3640.0*eta - 110880.0*eta*xi*xi + 84000.0*eta*(xi*xi*xi) - 30240.0*xi*eta*eta + 3360.0*(eta*eta) + 50400.0*(eta*eta)*(xi*xi), -13440.0*eta*xi + 560.0*eta + 60480.0*eta*(xi*xi) - 67200.0*eta*xi*xi*xi + 18200.0*xi - 700.0 - 100800.0*xi*xi + 184800.0*(xi*xi*xi) - 105000.0*xi*xi*xi*xi);
                      case 19:
                        return RealGradient(560.0*eta - 6720.0*eta*eta + 20160.0*(eta*eta*eta) - 16800.0*eta*eta*eta*eta, 6720.0*eta*xi - 3640.0*eta - 280.0*xi - 30240.0*xi*eta*eta + 33600.0*xi*(eta*eta*eta) + 140.0 + 20160.0*(eta*eta) - 36960.0*eta*eta*eta + 21000.0*(eta*eta*eta*eta));
                      case 20:
                        return RealGradient(-1120.0*eta + 13440.0*(eta*eta) - 40320.0*eta*eta*eta + 33600.0*(eta*eta*eta*eta), -13440.0*eta*xi + 6440.0*eta + 560.0*xi + 60480.0*xi*(eta*eta) - 67200.0*xi*eta*eta*eta - 280.0 - 30240.0*eta*eta + 43680.0*(eta*eta*eta) - 16800.0*eta*eta*eta*eta);
                      case 21:
                        return RealGradient(-35840.0*eta*xi + 5973.33333333333*eta + 52266.6666666667*eta*(xi*xi) - 19911.1111111111*eta*xi*xi*xi + 125440.0*xi*(eta*eta) - 89600.0*xi*eta*eta*eta - 29120.0*eta*eta - 97066.6666666667*eta*eta*xi*xi + 38826.6666666667*(eta*eta*eta) - 15555.5555555556*eta*eta*eta*eta, 56000.0*eta*xi - 3453.33333333333*eta - 170240.0*eta*xi*xi + 129422.222222222*eta*(xi*xi*xi) - 8306.66666666667*xi - 81760.0*xi*eta*eta + 31111.1111111111*xi*(eta*eta*eta) + 420.0 + 6720.0*(eta*eta) + 134400.0*(eta*eta)*(xi*xi) + 35840.0*(xi*xi) - 4355.55555555556*eta*eta*eta - 52764.4444444444*xi*xi*xi + 777.777777777778*(eta*eta*eta*eta) + 24888.8888888889*(xi*xi*xi*xi));
                      case 22:
                        return RealGradient(29120.0*eta*xi - 4293.33333333333*eta - 50773.3333333333*eta*xi*xi + 24888.8888888889*eta*(xi*xi*xi) - 109760.0*xi*eta*eta + 78400.0*xi*(eta*eta*eta) + 22400.0*(eta*eta) + 104533.333333333*(eta*eta)*(xi*xi) - 30613.3333333333*eta*eta*eta + 12444.4444444444*(eta*eta*eta*eta), -62720.0*eta*xi + 4013.33333333333*eta + 185920.0*eta*(xi*xi) - 139377.777777778*eta*xi*xi*xi + 10826.6666666667*xi + 76160.0*xi*(eta*eta) - 24888.8888888889*xi*eta*eta*eta - 560.0 - 6720.0*eta*eta - 117600.0*eta*eta*xi*xi - 45920.0*xi*xi + 3857.77777777778*(eta*eta*eta) + 66702.2222222222*(xi*xi*xi) - 622.222222222222*eta*eta*eta*eta - 31111.1111111111*xi*xi*xi*xi);
                      case 23:
                        return RealGradient(-13440.0*eta*xi + 4013.33333333333*eta + 11573.3333333333*eta*(xi*xi) - 2488.88888888889*eta*xi*xi*xi + 76160.0*xi*(eta*eta) - 78400.0*xi*eta*eta*eta - 31360.0*eta*eta - 37333.3333333333*eta*eta*xi*xi + 61973.3333333333*(eta*eta*eta) - 34844.4444444444*eta*eta*eta*eta, 44800.0*eta*xi - 4293.33333333333*eta - 91840.0*eta*xi*xi + 49777.7777777778*eta*(xi*xi*xi) - 3546.66666666667*xi - 109760.0*xi*eta*eta + 69688.8888888889*xi*(eta*eta*eta) + 280.0 + 14560.0*(eta*eta) + 117600.0*(eta*eta)*(xi*xi) + 10080.0*(xi*xi) - 16924.4444444444*eta*eta*eta - 9955.55555555555*xi*xi*xi + 6222.22222222222*(eta*eta*eta*eta) + 3111.11111111111*(xi*xi*xi*xi));
                      case 24:
                        return RealGradient(13440.0*eta*xi - 3453.33333333333*eta - 13066.6666666667*eta*xi*xi + 3111.11111111111*eta*(xi*xi*xi) - 81760.0*xi*eta*eta + 89600.0*xi*(eta*eta*eta) + 28000.0*(eta*eta) + 46666.6666666667*(eta*eta)*(xi*xi) - 56746.6666666667*eta*eta*eta + 32355.5555555556*(eta*eta*eta*eta), -58240.0*eta*xi + 5973.33333333333*eta + 116480.0*eta*(xi*xi) - 62222.2222222222*eta*xi*xi*xi + 4946.66666666667*xi + 125440.0*xi*(eta*eta) - 64711.1111111111*xi*eta*eta*eta - 420.0 - 17920.0*eta*eta - 134400.0*eta*eta*xi*xi - 13440.0*xi*xi + 17422.2222222222*(eta*eta*eta) + 12817.7777777778*(xi*xi*xi) - 4977.77777777778*eta*eta*eta*eta - 3888.88888888889*xi*xi*xi*xi);
                      case 25:
                        return RealGradient(11200.0*eta*xi - 1431.11111111111*eta - 16426.6666666667*eta*xi*xi + 4977.77777777778*eta*(xi*xi*xi) - 2240.0*xi*eta*eta - 11200.0*xi*eta*eta*eta - 5226.66666666667*eta*eta + 7466.66666666667*(eta*eta)*(xi*xi) + 17173.3333333333*(eta*eta*eta) - 10577.7777777778*eta*eta*eta*eta, 5226.66666666667*eta*xi - 1151.11111111111*eta + 2240.0*eta*(xi*xi) - 9955.55555555555*eta*xi*xi*xi + 715.555555555556*xi - 25760.0*xi*eta*eta + 21155.5555555556*xi*(eta*eta*eta) + 15.5555555555556 + 5600.0*(eta*eta) + 16800.0*(eta*eta)*(xi*xi) - 5600.0*xi*xi - 8337.77777777778*eta*eta*eta + 10951.1111111111*(xi*xi*xi) + 3888.88888888889*(eta*eta*eta*eta) - 6222.22222222222*xi*xi*xi*xi);
                      case 26:
                        return RealGradient(-2240.0*eta*xi - 1057.77777777778*eta + 23893.3333333333*eta*(xi*xi) - 24888.8888888889*eta*xi*xi*xi - 51520.0*xi*eta*eta + 56000.0*xi*(eta*eta*eta) + 14933.3333333333*(eta*eta) + 29866.6666666667*(eta*eta)*(xi*xi) - 27626.6666666667*eta*eta*eta + 13688.8888888889*(eta*eta*eta*eta), -46293.3333333333*eta*xi + 5755.55555555556*eta + 82880.0*eta*(xi*xi) - 39822.2222222222*eta*xi*xi*xi - 5351.11111111111*xi + 80640.0*xi*(eta*eta) - 27377.7777777778*xi*eta*eta*eta + 62.2222222222222 - 13440.0*eta*eta - 84000.0*eta*eta*xi*xi + 32480.0*(xi*xi) + 8835.55555555555*(eta*eta*eta) - 57742.2222222222*xi*xi*xi - 1244.44444444444*eta*eta*eta*eta + 31111.1111111111*(xi*xi*xi*xi));
                      case 27:
                        return RealGradient(-1120.0*eta*xi + 311.111111111111*eta - 373.333333333333*eta*xi*xi + 622.222222222222*eta*(xi*xi*xi) + 12320.0*xi*(eta*eta) - 5600.0*xi*eta*eta*eta - 1493.33333333333*eta*eta - 7466.66666666667*eta*eta*xi*xi - 3733.33333333333*eta*eta*eta + 4977.77777777778*(eta*eta*eta*eta), 1493.33333333333*eta*xi + 591.111111111111*eta - 12320.0*eta*xi*xi + 9955.55555555555*eta*(xi*xi*xi) - 155.555555555556*xi + 5600.0*xi*(eta*eta) - 9955.55555555555*xi*eta*eta*eta - 15.5555555555556 - 3920.0*eta*eta + 8400.0*(eta*eta)*(xi*xi) + 560.0*(xi*xi) + 7217.77777777778*(eta*eta*eta) + 248.888888888889*(xi*xi*xi) - 3888.88888888889*eta*eta*eta*eta - 777.777777777778*xi*xi*xi*xi);
                      case 28:
                        return RealGradient(3360.0*eta*xi + 31.1111111111111*eta - 3733.33333333333*eta*xi*xi - 3111.11111111111*eta*xi*xi*xi - 19040.0*xi*eta*eta - 5600.0*xi*eta*eta*eta - 3733.33333333333*eta*eta + 37333.3333333333*(eta*eta)*(xi*xi) + 17546.6666666667*(eta*eta*eta) - 13688.8888888889*eta*eta*eta*eta, -11946.6666666667*eta*xi - 995.555555555556*eta + 58240.0*eta*(xi*xi) - 49777.7777777778*eta*xi*xi*xi + 964.444444444444*xi - 30240.0*xi*eta*eta + 27377.7777777778*xi*(eta*eta*eta) - 15.5555555555556 + 9520.0*(eta*eta) + 8400.0*(eta*eta)*(xi*xi) - 1680.0*xi*xi - 11075.5555555556*eta*eta*eta - 2737.77777777778*xi*xi*xi + 2644.44444444444*(eta*eta*eta*eta) + 3888.88888888889*(xi*xi*xi*xi));
                      case 29:
                        return RealGradient(-26880.0*eta*xi + 5755.55555555556*eta + 26506.6666666667*eta*(xi*xi) - 4977.77777777778*eta*xi*xi*xi + 80640.0*xi*(eta*eta) - 56000.0*xi*eta*eta*eta - 23146.6666666667*eta*eta - 41066.6666666667*eta*eta*xi*xi + 27626.6666666667*(eta*eta*eta) - 9955.55555555555*eta*eta*eta*eta, 29866.6666666667*eta*xi - 1057.77777777778*eta - 82880.0*eta*xi*xi + 54755.5555555556*eta*(xi*xi*xi) - 3297.77777777778*xi - 51520.0*xi*eta*eta + 19911.1111111111*xi*(eta*eta*eta) + 155.555555555556 - 1120.0*eta*eta + 84000.0*(eta*eta)*(xi*xi) + 12320.0*(xi*xi) + 7964.44444444444*(eta*eta*eta) - 15431.1111111111*xi*xi*xi - 6222.22222222222*eta*eta*eta*eta + 6222.22222222222*(xi*xi*xi*xi));
                      case 30:
                        return RealGradient(11200.0*eta*xi - 1151.11111111111*eta - 25013.3333333333*eta*xi*xi + 15555.5555555556*eta*(xi*xi*xi) - 25760.0*xi*eta*eta + 11200.0*xi*(eta*eta*eta) + 2613.33333333333*(eta*eta) + 31733.3333333333*(eta*eta)*(xi*xi) + 746.666666666667*(eta*eta*eta) - 2488.88888888889*eta*eta*eta*eta, -10453.3333333333*eta*xi - 1431.11111111111*eta + 51520.0*eta*(xi*xi) - 42311.1111111111*eta*xi*xi*xi + 5755.55555555556*xi - 2240.0*xi*eta*eta + 4977.77777777778*xi*(eta*eta*eta) - 77.7777777777778 + 5600.0*(eta*eta) - 16800.0*eta*eta*xi*xi - 28000.0*xi*xi - 5475.55555555556*eta*eta*eta + 41688.8888888889*(xi*xi*xi) + 1244.44444444444*(eta*eta*eta*eta) - 19444.4444444444*xi*xi*xi*xi);
                      case 31:
                        return RealGradient(19040.0*eta*xi - 995.555555555556*eta - 33226.6666666667*eta*xi*xi + 10577.7777777778*eta*(xi*xi*xi) - 30240.0*xi*eta*eta + 5600.0*xi*(eta*eta*eta) - 5973.33333333333*eta*eta + 41066.6666666667*(eta*eta)*(xi*xi) + 19413.3333333333*(eta*eta*eta) - 12444.4444444444*eta*eta*eta*eta, -7466.66666666667*eta*xi + 31.1111111111111*eta + 52640.0*eta*(xi*xi) - 54755.5555555556*eta*xi*xi*xi + 3017.77777777778*xi - 19040.0*xi*eta*eta + 24888.8888888889*xi*(eta*eta*eta) - 108.888888888889 + 1680.0*(eta*eta) - 8400.0*eta*eta*xi*xi - 16240.0*xi*xi - 1244.44444444444*eta*eta*eta + 26631.1111111111*(xi*xi*xi) - 777.777777777778*eta*eta*eta*eta - 13222.2222222222*xi*xi*xi*xi);
                      case 32:
                        return RealGradient(-7840.0*eta*xi + 591.111111111111*eta + 21653.3333333333*eta*(xi*xi) - 15555.5555555556*eta*xi*xi*xi + 5600.0*xi*(eta*eta) + 5600.0*xi*(eta*eta*eta) + 746.666666666667*(eta*eta) - 14933.3333333333*eta*eta*xi*xi - 4106.66666666667*eta*eta*eta + 2488.88888888889*(eta*eta*eta*eta), -2986.66666666667*eta*xi + 311.111111111111*eta - 11200.0*eta*xi*xi + 19911.1111111111*eta*(xi*xi*xi) - 2955.55555555556*xi + 12320.0*xi*(eta*eta) - 4977.77777777778*xi*eta*eta*eta + 77.7777777777778 - 560.0*eta*eta - 8400.0*eta*eta*xi*xi + 19600.0*(xi*xi) - 124.444444444444*eta*eta*eta - 36088.8888888889*xi*xi*xi + 155.555555555556*(eta*eta*eta*eta) + 19444.4444444444*(xi*xi*xi*xi));
                      case 33:
                        return RealGradient(3360.0*eta*xi - 715.555555555556*eta - 4480.0*eta*xi*xi + 1866.66666666667*eta*(xi*xi*xi) - 19040.0*xi*eta*eta + 16800.0*xi*(eta*eta*eta) + 6346.66666666667*(eta*eta) + 11200.0*(eta*eta)*(xi*xi) - 11200.0*eta*eta*eta + 5600.0*(eta*eta*eta*eta), -9706.66666666667*eta*xi + 684.444444444444*eta + 23520.0*eta*(xi*xi) - 14933.3333333333*eta*xi*xi*xi + 1057.77777777778*xi + 20160.0*xi*(eta*eta) - 11200.0*xi*eta*eta*eta - 62.2222222222222 - 1680.0*eta*eta - 25200.0*eta*eta*xi*xi - 3920.0*xi*xi + 1120.0*(eta*eta*eta) + 5226.66666666667*(xi*xi*xi) - 2333.33333333333*xi*xi*xi*xi);
                      case 34:
                        return RealGradient(-3360.0*eta*xi + 684.444444444444*eta + 3360.0*eta*(xi*xi) + 20160.0*xi*(eta*eta) - 16800.0*xi*eta*eta*eta - 4853.33333333333*eta*eta - 16800.0*eta*eta*xi*xi + 7840.0*(eta*eta*eta) - 3733.33333333333*eta*eta*eta*eta, 12693.3333333333*eta*xi - 715.555555555556*eta - 33600.0*eta*xi*xi + 22400.0*eta*(xi*xi*xi) - 622.222222222222*xi - 19040.0*xi*eta*eta + 7466.66666666667*xi*(eta*eta*eta) + 31.1111111111111 + 1680.0*(eta*eta) + 25200.0*(eta*eta)*(xi*xi) + 1680.0*(xi*xi) - 1493.33333333333*eta*eta*eta - 1120.0*xi*xi*xi + 466.666666666667*(eta*eta*eta*eta));
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
                        return sign * RealGradient(-16800.0*eta*xi + 3500.0*eta + 22680.0*eta*(xi*xi) - 8400.0*eta*xi*xi*xi + 2800.0*xi + 30240.0*xi*(eta*eta) - 16800.0*xi*eta*eta*eta - 375.0 - 10500.0*eta*eta - 18900.0*eta*eta*xi*xi - 6300.0*xi*xi + 12600.0*(eta*eta*eta) + 5040.0*(xi*xi*xi) - 5250.0*eta*eta*eta*eta - 1050.0*xi*xi*xi*xi, 4200.0*eta*xi - 15120.0*eta*xi*xi + 12600.0*eta*(xi*xi*xi) - 700.0*xi - 7560.0*xi*eta*eta + 4200.0*xi*(eta*eta*eta) + 12600.0*(eta*eta)*(xi*xi) + 4200.0*(xi*xi) - 7560.0*xi*xi*xi + 4200.0*(xi*xi*xi*xi));
                      case 1:
                        return sign * RealGradient(3622.5*eta*xi + 351.09375*eta - 7530.46875*eta*xi*xi + 3051.5625*eta*(xi*xi*xi) - 759.0625*xi - 4784.0625*xi*eta*eta + 1771.875*xi*(eta*eta*eta) + 5.4296875 - 1993.359375*eta*eta + 5980.078125*(eta*eta)*(xi*xi) + 2127.890625*(xi*xi) + 3277.96875*(eta*eta*eta) - 1706.25*xi*xi*xi - 1661.1328125*eta*eta*eta*eta + 303.515625*(xi*xi*xi*xi), 88.59375*eta*xi + 3189.375*eta*(xi*xi) - 3986.71875*eta*xi*xi*xi + 95.15625*xi - 1328.90625*xi*eta*eta + 1328.90625*xi*(eta*eta*eta) - 1328.90625*eta*eta*xi*xi - 1230.46875*xi*xi + 2628.28125*(xi*xi*xi) - 1525.78125*xi*xi*xi*xi);
                      case 2:
                        return sign * RealGradient(1680.0*eta*xi - 87.5*eta + 472.5*eta*(xi*xi) - 1575.0*eta*xi*xi*xi - 4725.0*xi*eta*eta + 3150.0*xi*(eta*eta*eta) + 6.875 + 26.25*(eta*eta) + 1181.25*(eta*eta)*(xi*xi) - 813.75*xi*xi + 367.5*(eta*eta*eta) + 1260.0*(xi*xi*xi) - 328.125*eta*eta*eta*eta - 393.75*xi*xi*xi*xi, -262.5*eta*xi + 1890.0*eta*(xi*xi) - 787.5*eta*xi*xi*xi + 17.5*xi + 157.5*xi*(eta*eta) + 262.5*xi*(eta*eta*eta) - 2362.5*eta*eta*xi*xi + 52.5*(xi*xi) - 787.5*xi*xi*xi + 787.5*(xi*xi*xi*xi));
                      case 3:
                        return sign * RealGradient(-577.5*eta*xi + 18.59375*eta + 2392.03125*eta*(xi*xi) - 623.4375*eta*xi*xi*xi + 115.9375*xi - 59.0625*xi*eta*eta + 721.875*xi*(eta*eta*eta) - 6.4453125 + 27.890625*(eta*eta) - 2288.671875*eta*eta*xi*xi - 260.859375*xi*xi - 29.53125*eta*eta*eta - 236.25*xi*xi*xi - 20.5078125*eta*eta*eta*eta + 303.515625*(xi*xi*xi*xi), 36.09375*eta*xi - 590.625*eta*xi*xi + 1525.78125*eta*(xi*xi*xi) - 27.34375*xi + 88.59375*xi*(eta*eta) + 16.40625*xi*(eta*eta*eta) - 541.40625*eta*eta*xi*xi + 292.03125*(xi*xi) - 679.21875*xi*xi*xi + 311.71875*(xi*xi*xi*xi));
                      case 4:
                        return sign * RealGradient(140.0*xi - 5.0 - 840.0*xi*xi + 1680.0*(xi*xi*xi) - 1050.0*xi*xi*xi*xi, 0);
                      case 5:
                        return sign * RealGradient(140.0*xi - 5.0 - 840.0*xi*xi + 1680.0*(xi*xi*xi) - 1050.0*xi*xi*xi*xi, 0);
                      case 6:
                        return sign * RealGradient(-1575.0*eta*xi + 96.25*eta + 4784.0625*eta*(xi*xi) - 3543.75*eta*xi*xi*xi + 183.75*xi + 2303.4375*xi*(eta*eta) - 787.5*xi*eta*eta*eta - 10.0 - 183.75*eta*eta - 3986.71875*eta*eta*xi*xi - 708.75*xi*xi + 111.5625*(eta*eta*eta) + 885.9375*(xi*xi*xi) - 20.5078125*eta*eta*eta*eta - 332.2265625*xi*xi*xi*xi, 262.5*eta*xi - 1949.0625*eta*xi*xi + 2657.8125*eta*(xi*xi*xi) - 113.75*xi - 137.8125*xi*eta*eta + 16.40625*xi*(eta*eta*eta) + 590.625*(eta*eta)*(xi*xi) + 1102.5*(xi*xi) - 2657.8125*xi*xi*xi + 1771.875*(xi*xi*xi*xi));
                      case 7:
                        return sign * RealGradient(-2520.0*eta*xi + 245.0*eta + 4725.0*eta*(xi*xi) - 2100.0*eta*xi*xi*xi + 175.0*xi + 6615.0*xi*(eta*eta) - 4200.0*xi*eta*eta*eta - 15.0 - 840.0*eta*eta - 7087.5*eta*eta*xi*xi - 420.0*xi*xi + 945.0*(eta*eta*eta) + 315.0*(xi*xi*xi) - 328.125*eta*eta*eta*eta - 65.625*xi*xi*xi*xi, 840.0*eta*xi - 4725.0*eta*xi*xi + 4725.0*eta*(xi*xi*xi) - 175.0*xi - 945.0*xi*eta*eta + 262.5*xi*(eta*eta*eta) + 3150.0*(eta*eta)*(xi*xi) + 1260.0*(xi*xi) - 2205.0*xi*xi*xi + 1050.0*(xi*xi*xi*xi));
                      case 8:
                        return sign * RealGradient(-2205.0*eta*xi + 446.25*eta + 1949.0625*eta*(xi*xi) - 393.75*eta*xi*xi*xi + 113.75*xi + 7973.4375*xi*(eta*eta) - 7087.5*xi*eta*eta*eta - 20.0 - 2126.25*eta*eta - 3986.71875*eta*eta*xi*xi - 131.25*xi*xi + 3366.5625*(eta*eta*eta) + 45.9375*(xi*xi*xi) - 1661.1328125*eta*eta*eta*eta - 4.1015625*xi*xi*xi*xi, 1417.5*eta*xi - 4784.0625*eta*xi*xi + 2657.8125*eta*(xi*xi*xi) - 183.75*xi - 2657.8125*xi*eta*eta + 1328.90625*xi*(eta*eta*eta) + 5315.625*(eta*eta)*(xi*xi) + 787.5*(xi*xi) - 767.8125*xi*xi*xi + 196.875*(xi*xi*xi*xi));
                      case 9:
                        return sign * RealGradient(700.0*eta - 25.0 - 4200.0*eta*eta + 8400.0*(eta*eta*eta) - 5250.0*eta*eta*eta*eta, 1680.0*eta*xi - 140.0*xi - 5040.0*xi*eta*eta + 4200.0*xi*(eta*eta*eta));
                      case 10:
                        return sign * RealGradient(700.0*eta - 25.0 - 4200.0*eta*eta + 8400.0*(eta*eta*eta) - 5250.0*eta*eta*eta*eta, 1680.0*eta*xi - 1260.0*eta - 140.0*xi - 5040.0*xi*eta*eta + 4200.0*xi*(eta*eta*eta) + 120.0 + 3360.0*(eta*eta) - 2520.0*eta*eta*eta);
                      case 11:
                        return sign * RealGradient(-584.0625*eta*xi - 148.75*eta + 590.625*eta*(xi*xi) + 360.9375*eta*(xi*xi*xi) + 27.34375*xi + 2037.65625*xi*(eta*eta) - 1246.875*xi*eta*eta*eta + 4.3359375 + 1038.515625*(eta*eta) - 2288.671875*eta*eta*xi*xi - 18.046875*xi*xi - 2303.4375*eta*eta*eta - 29.53125*xi*xi*xi + 1517.578125*(eta*eta*eta*eta) - 4.1015625*xi*xi*xi*xi, 521.71875*eta*xi + 292.03125*eta - 2392.03125*eta*xi*xi + 1525.78125*eta*(xi*xi*xi) - 115.9375*xi + 708.75*xi*(eta*eta) - 1214.0625*xi*eta*eta*eta - 22.03125 - 908.90625*eta*eta + 935.15625*(eta*eta)*(xi*xi) + 288.75*(xi*xi) + 728.4375*(eta*eta*eta) + 19.6875*(xi*xi*xi) - 180.46875*xi*xi*xi*xi);
                      case 12:
                        return sign * RealGradient(-105.0*eta*xi + 315.0*eta - 1890.0*eta*xi*xi + 1575.0*eta*(xi*xi*xi) - 17.5*xi + 2362.5*xi*(eta*eta) - 3150.0*xi*eta*eta*eta - 10.625 - 1863.75*eta*eta + 1181.25*(eta*eta)*(xi*xi) + 131.25*(xi*xi) + 3465.0*(eta*eta*eta) - 52.5*xi*xi*xi - 1968.75*eta*eta*eta*eta - 65.625*xi*xi*xi*xi, 1627.5*eta*xi - 577.5*eta - 472.5*eta*xi*xi - 787.5*eta*xi*xi*xi - 3780.0*xi*eta*eta + 1575.0*xi*(eta*eta*eta) + 52.5 + 1417.5*(eta*eta) + 2362.5*(eta*eta)*(xi*xi) - 840.0*xi*xi - 945.0*eta*eta*eta + 1575.0*(xi*xi*xi) - 787.5*xi*xi*xi*xi);
                      case 13:
                        return sign * RealGradient(2460.9375*eta*xi - 428.75*eta - 3189.375*eta*xi*xi + 885.9375*eta*(xi*xi*xi) - 95.15625*xi - 7884.84375*xi*eta*eta + 6103.125*xi*(eta*eta*eta) + 18.7109375 + 1957.265625*(eta*eta) + 5980.078125*(eta*eta)*(xi*xi) - 44.296875*xi*xi - 3038.4375*eta*eta*eta + 442.96875*(xi*xi*xi) + 1517.578125*(eta*eta*eta*eta) - 332.2265625*xi*xi*xi*xi, -4255.78125*eta*xi + 659.53125*eta + 7530.46875*eta*(xi*xi) - 3986.71875*eta*xi*xi*xi + 759.0625*xi + 5118.75*xi*(eta*eta) - 1214.0625*xi*eta*eta*eta - 89.53125 - 1276.40625*eta*eta - 4577.34375*eta*eta*xi*xi - 1811.25*xi*xi + 728.4375*(eta*eta*eta) + 1594.6875*(xi*xi*xi) - 442.96875*xi*xi*xi*xi);
                      case 14:
                        return sign * RealGradient(-8400.0*eta*xi + 1400.0*eta + 15120.0*eta*(xi*xi) - 8400.0*eta*xi*xi*xi + 700.0*xi + 22680.0*xi*(eta*eta) - 16800.0*xi*eta*eta*eta - 75.0 - 6300.0*eta*eta - 18900.0*eta*eta*xi*xi - 2100.0*xi*xi + 10080.0*(eta*eta*eta) + 2520.0*(xi*xi*xi) - 5250.0*eta*eta*eta*eta - 1050.0*xi*xi*xi*xi, 12600.0*eta*xi - 2100.0*eta - 22680.0*eta*xi*xi + 12600.0*eta*(xi*xi*xi) - 2800.0*xi - 15120.0*xi*eta*eta + 4200.0*xi*(eta*eta*eta) + 300.0 + 4200.0*(eta*eta) + 12600.0*(eta*eta)*(xi*xi) + 8400.0*(xi*xi) - 2520.0*eta*eta*eta - 10080.0*xi*xi*xi + 4200.0*(xi*xi*xi*xi));
                      case 15:
                        return RealGradient(188160.0*eta*xi - 49000.0*eta - 211680.0*eta*xi*xi + 67200.0*eta*(xi*xi*xi) - 19600.0*xi - 423360.0*xi*eta*eta + 268800.0*xi*(eta*eta*eta) + 3500.0 + 176400.0*(eta*eta) + 226800.0*(eta*eta)*(xi*xi) + 35280.0*(xi*xi) - 235200.0*eta*eta*eta - 23520.0*xi*xi*xi + 105000.0*(eta*eta*eta*eta) + 4200.0*(xi*xi*xi*xi), -70560.0*eta*xi + 211680.0*eta*(xi*xi) - 151200.0*eta*xi*xi*xi + 9800.0*xi + 141120.0*xi*(eta*eta) - 84000.0*xi*eta*eta*eta - 201600.0*eta*eta*xi*xi - 47040.0*xi*xi + 70560.0*(xi*xi*xi) - 33600.0*xi*xi*xi*xi);
                      case 16:
                        return RealGradient(-94080.0*eta*xi + 9800.0*eta + 211680.0*eta*(xi*xi) - 134400.0*eta*xi*xi*xi + 9800.0*xi + 211680.0*xi*(eta*eta) - 134400.0*xi*eta*eta*eta - 700.0 - 35280.0*eta*eta - 226800.0*eta*eta*xi*xi - 35280.0*xi*xi + 47040.0*(eta*eta*eta) + 47040.0*(xi*xi*xi) - 21000.0*eta*eta*eta*eta - 21000.0*xi*xi*xi*xi, 70560.0*eta*xi - 211680.0*eta*xi*xi + 151200.0*eta*(xi*xi*xi) - 19600.0*xi - 70560.0*xi*eta*eta + 16800.0*xi*(eta*eta*eta) + 100800.0*(eta*eta)*(xi*xi) + 94080.0*(xi*xi) - 141120.0*xi*xi*xi + 67200.0*(xi*xi*xi*xi));
                      case 17:
                        return RealGradient(-13440.0*eta*xi + 560.0*eta + 60480.0*eta*(xi*xi) - 67200.0*eta*xi*xi*xi + 6440.0*xi - 280.0 - 30240.0*xi*xi + 43680.0*(xi*xi*xi) - 16800.0*xi*xi*xi*xi, -1120.0*xi + 13440.0*(xi*xi) - 40320.0*xi*xi*xi + 33600.0*(xi*xi*xi*xi));
                      case 18:
                        return RealGradient(6720.0*eta*xi - 280.0*eta - 30240.0*eta*xi*xi + 33600.0*eta*(xi*xi*xi) - 3640.0*xi + 140.0 + 20160.0*(xi*xi) - 36960.0*xi*xi*xi + 21000.0*(xi*xi*xi*xi), 560.0*xi - 6720.0*xi*xi + 20160.0*(xi*xi*xi) - 16800.0*xi*xi*xi*xi);
                      case 19:
                        return RealGradient(-13440.0*eta*xi + 18200.0*eta + 560.0*xi + 60480.0*xi*(eta*eta) - 67200.0*xi*eta*eta*eta - 700.0 - 100800.0*eta*eta + 184800.0*(eta*eta*eta) - 105000.0*eta*eta*eta*eta, 40320.0*eta*xi - 30240.0*eta*xi*xi - 3640.0*xi - 110880.0*xi*eta*eta + 84000.0*xi*(eta*eta*eta) + 50400.0*(eta*eta)*(xi*xi) + 3360.0*(xi*xi));
                      case 20:
                        return RealGradient(26880.0*eta*xi - 14560.0*eta - 1120.0*xi - 120960.0*xi*eta*eta + 134400.0*xi*(eta*eta*eta) + 560.0 + 80640.0*(eta*eta) - 147840.0*eta*eta*eta + 84000.0*(eta*eta*eta*eta), -60480.0*eta*xi + 60480.0*eta*(xi*xi) + 6440.0*xi + 131040.0*xi*(eta*eta) - 67200.0*xi*eta*eta*eta - 100800.0*eta*eta*xi*xi - 6720.0*xi*xi);
                      case 21:
                        return RealGradient(-58240.0*eta*xi + 4946.66666666667*eta + 125440.0*eta*(xi*xi) - 64711.1111111111*eta*xi*xi*xi + 5973.33333333333*xi + 116480.0*xi*(eta*eta) - 62222.2222222222*xi*eta*eta*eta - 420.0 - 13440.0*eta*eta - 134400.0*eta*eta*xi*xi - 17920.0*xi*xi + 12817.7777777778*(eta*eta*eta) + 17422.2222222222*(xi*xi*xi) - 3888.88888888889*eta*eta*eta*eta - 4977.77777777778*xi*xi*xi*xi, 13440.0*eta*xi - 81760.0*eta*xi*xi + 89600.0*eta*(xi*xi*xi) - 3453.33333333333*xi - 13066.6666666667*xi*eta*eta + 3111.11111111111*xi*(eta*eta*eta) + 46666.6666666667*(eta*eta)*(xi*xi) + 28000.0*(xi*xi) - 56746.6666666667*xi*xi*xi + 32355.5555555556*(xi*xi*xi*xi));
                      case 22:
                        return RealGradient(44800.0*eta*xi - 3546.66666666667*eta - 109760.0*eta*xi*xi + 69688.8888888889*eta*(xi*xi*xi) - 4293.33333333333*xi - 91840.0*xi*eta*eta + 49777.7777777778*xi*(eta*eta*eta) + 280.0 + 10080.0*(eta*eta) + 117600.0*(eta*eta)*(xi*xi) + 14560.0*(xi*xi) - 9955.55555555555*eta*eta*eta - 16924.4444444444*xi*xi*xi + 3111.11111111111*(eta*eta*eta*eta) + 6222.22222222222*(xi*xi*xi*xi), -13440.0*eta*xi + 76160.0*eta*(xi*xi) - 78400.0*eta*xi*xi*xi + 4013.33333333333*xi + 11573.3333333333*xi*(eta*eta) - 2488.88888888889*xi*eta*eta*eta - 37333.3333333333*eta*eta*xi*xi - 31360.0*xi*xi + 61973.3333333333*(xi*xi*xi) - 34844.4444444444*xi*xi*xi*xi);
                      case 23:
                        return RealGradient(-62720.0*eta*xi + 10826.6666666667*eta + 76160.0*eta*(xi*xi) - 24888.8888888889*eta*xi*xi*xi + 4013.33333333333*xi + 185920.0*xi*(eta*eta) - 139377.777777778*xi*eta*eta*eta - 560.0 - 45920.0*eta*eta - 117600.0*eta*eta*xi*xi - 6720.0*xi*xi + 66702.2222222222*(eta*eta*eta) + 3857.77777777778*(xi*xi*xi) - 31111.1111111111*eta*eta*eta*eta - 622.222222222222*xi*xi*xi*xi, 29120.0*eta*xi - 109760.0*eta*xi*xi + 78400.0*eta*(xi*xi*xi) - 4293.33333333333*xi - 50773.3333333333*xi*eta*eta + 24888.8888888889*xi*(eta*eta*eta) + 104533.333333333*(eta*eta)*(xi*xi) + 22400.0*(xi*xi) - 30613.3333333333*xi*xi*xi + 12444.4444444444*(xi*xi*xi*xi));
                      case 24:
                        return RealGradient(56000.0*eta*xi - 8306.66666666667*eta - 81760.0*eta*xi*xi + 31111.1111111111*eta*(xi*xi*xi) - 3453.33333333333*xi - 170240.0*xi*eta*eta + 129422.222222222*xi*(eta*eta*eta) + 420.0 + 35840.0*(eta*eta) + 134400.0*(eta*eta)*(xi*xi) + 6720.0*(xi*xi) - 52764.4444444444*eta*eta*eta - 4355.55555555556*xi*xi*xi + 24888.8888888889*(eta*eta*eta*eta) + 777.777777777778*(xi*xi*xi*xi), -35840.0*eta*xi + 125440.0*eta*(xi*xi) - 89600.0*eta*xi*xi*xi + 5973.33333333333*xi + 52266.6666666667*xi*(eta*eta) - 19911.1111111111*xi*eta*eta*eta - 97066.6666666667*eta*eta*xi*xi - 29120.0*xi*xi + 38826.6666666667*(xi*xi*xi) - 15555.5555555556*xi*xi*xi*xi);
                      case 25:
                        return RealGradient(-10453.3333333333*eta*xi + 5755.55555555556*eta - 2240.0*eta*xi*xi + 4977.77777777778*eta*(xi*xi*xi) - 1431.11111111111*xi + 51520.0*xi*(eta*eta) - 42311.1111111111*xi*eta*eta*eta - 77.7777777777778 - 28000.0*eta*eta - 16800.0*eta*eta*xi*xi + 5600.0*(xi*xi) + 41688.8888888889*(eta*eta*eta) - 5475.55555555556*xi*xi*xi - 19444.4444444444*eta*eta*eta*eta + 1244.44444444444*(xi*xi*xi*xi), 11200.0*eta*xi - 25760.0*eta*xi*xi + 11200.0*eta*(xi*xi*xi) - 1151.11111111111*xi - 25013.3333333333*xi*eta*eta + 15555.5555555556*xi*(eta*eta*eta) + 31733.3333333333*(eta*eta)*(xi*xi) + 2613.33333333333*(xi*xi) + 746.666666666667*(xi*xi*xi) - 2488.88888888889*xi*xi*xi*xi);
                      case 26:
                        return RealGradient(29866.6666666667*eta*xi - 3297.77777777778*eta - 51520.0*eta*xi*xi + 19911.1111111111*eta*(xi*xi*xi) - 1057.77777777778*xi - 82880.0*xi*eta*eta + 54755.5555555556*xi*(eta*eta*eta) + 155.555555555556 + 12320.0*(eta*eta) + 84000.0*(eta*eta)*(xi*xi) - 1120.0*xi*xi - 15431.1111111111*eta*eta*eta + 7964.44444444444*(xi*xi*xi) + 6222.22222222222*(eta*eta*eta*eta) - 6222.22222222222*xi*xi*xi*xi, -26880.0*eta*xi + 80640.0*eta*(xi*xi) - 56000.0*eta*xi*xi*xi + 5755.55555555556*xi + 26506.6666666667*xi*(eta*eta) - 4977.77777777778*xi*eta*eta*eta - 41066.6666666667*eta*eta*xi*xi - 23146.6666666667*xi*xi + 27626.6666666667*(xi*xi*xi) - 9955.55555555555*xi*xi*xi*xi);
                      case 27:
                        return RealGradient(-2986.66666666667*eta*xi - 2955.55555555556*eta + 12320.0*eta*(xi*xi) - 4977.77777777778*eta*xi*xi*xi + 311.111111111111*xi - 11200.0*xi*eta*eta + 19911.1111111111*xi*(eta*eta*eta) + 77.7777777777778 + 19600.0*(eta*eta) - 8400.0*eta*eta*xi*xi - 560.0*xi*xi - 36088.8888888889*eta*eta*eta - 124.444444444444*xi*xi*xi + 19444.4444444444*(eta*eta*eta*eta) + 155.555555555556*(xi*xi*xi*xi), -7840.0*eta*xi + 5600.0*eta*(xi*xi) + 5600.0*eta*(xi*xi*xi) + 591.111111111111*xi + 21653.3333333333*xi*(eta*eta) - 15555.5555555556*xi*eta*eta*eta - 14933.3333333333*eta*eta*xi*xi + 746.666666666667*(xi*xi) - 4106.66666666667*xi*xi*xi + 2488.88888888889*(xi*xi*xi*xi));
                      case 28:
                        return RealGradient(-7466.66666666667*eta*xi + 3017.77777777778*eta - 19040.0*eta*xi*xi + 24888.8888888889*eta*(xi*xi*xi) + 31.1111111111111*xi + 52640.0*xi*(eta*eta) - 54755.5555555556*xi*eta*eta*eta - 108.888888888889 - 16240.0*eta*eta - 8400.0*eta*eta*xi*xi + 1680.0*(xi*xi) + 26631.1111111111*(eta*eta*eta) - 1244.44444444444*xi*xi*xi - 13222.2222222222*eta*eta*eta*eta - 777.777777777778*xi*xi*xi*xi, 19040.0*eta*xi - 30240.0*eta*xi*xi + 5600.0*eta*(xi*xi*xi) - 995.555555555556*xi - 33226.6666666667*xi*eta*eta + 10577.7777777778*xi*(eta*eta*eta) + 41066.6666666667*(eta*eta)*(xi*xi) - 5973.33333333333*xi*xi + 19413.3333333333*(xi*xi*xi) - 12444.4444444444*xi*xi*xi*xi);
                      case 29:
                        return RealGradient(-46293.3333333333*eta*xi - 5351.11111111111*eta + 80640.0*eta*(xi*xi) - 27377.7777777778*eta*xi*xi*xi + 5755.55555555556*xi + 82880.0*xi*(eta*eta) - 39822.2222222222*xi*eta*eta*eta + 62.2222222222222 + 32480.0*(eta*eta) - 84000.0*eta*eta*xi*xi - 13440.0*xi*xi - 57742.2222222222*eta*eta*eta + 8835.55555555555*(xi*xi*xi) + 31111.1111111111*(eta*eta*eta*eta) - 1244.44444444444*xi*xi*xi*xi, -2240.0*eta*xi - 51520.0*eta*xi*xi + 56000.0*eta*(xi*xi*xi) - 1057.77777777778*xi + 23893.3333333333*xi*(eta*eta) - 24888.8888888889*xi*eta*eta*eta + 29866.6666666667*(eta*eta)*(xi*xi) + 14933.3333333333*(xi*xi) - 27626.6666666667*xi*xi*xi + 13688.8888888889*(xi*xi*xi*xi));
                      case 30:
                        return RealGradient(5226.66666666667*eta*xi + 715.555555555556*eta - 25760.0*eta*xi*xi + 21155.5555555556*eta*(xi*xi*xi) - 1151.11111111111*xi + 2240.0*xi*(eta*eta) - 9955.55555555555*xi*eta*eta*eta + 15.5555555555556 - 5600.0*eta*eta + 16800.0*(eta*eta)*(xi*xi) + 5600.0*(xi*xi) + 10951.1111111111*(eta*eta*eta) - 8337.77777777778*xi*xi*xi - 6222.22222222222*eta*eta*eta*eta + 3888.88888888889*(xi*xi*xi*xi), 11200.0*eta*xi - 2240.0*eta*xi*xi - 11200.0*eta*xi*xi*xi - 1431.11111111111*xi - 16426.6666666667*xi*eta*eta + 4977.77777777778*xi*(eta*eta*eta) + 7466.66666666667*(eta*eta)*(xi*xi) - 5226.66666666667*xi*xi + 17173.3333333333*(xi*xi*xi) - 10577.7777777778*xi*xi*xi*xi);
                      case 31:
                        return RealGradient(-11946.6666666667*eta*xi + 964.444444444444*eta - 30240.0*eta*xi*xi + 27377.7777777778*eta*(xi*xi*xi) - 995.555555555556*xi + 58240.0*xi*(eta*eta) - 49777.7777777778*xi*eta*eta*eta - 15.5555555555556 - 1680.0*eta*eta + 8400.0*(eta*eta)*(xi*xi) + 9520.0*(xi*xi) - 2737.77777777778*eta*eta*eta - 11075.5555555556*xi*xi*xi + 3888.88888888889*(eta*eta*eta*eta) + 2644.44444444444*(xi*xi*xi*xi), 3360.0*eta*xi - 19040.0*eta*xi*xi - 5600.0*eta*xi*xi*xi + 31.1111111111111*xi - 3733.33333333333*xi*eta*eta - 3111.11111111111*xi*eta*eta*eta + 37333.3333333333*(eta*eta)*(xi*xi) - 3733.33333333333*xi*xi + 17546.6666666667*(xi*xi*xi) - 13688.8888888889*xi*xi*xi*xi);
                      case 32:
                        return RealGradient(1493.33333333333*eta*xi - 155.555555555556*eta + 5600.0*eta*(xi*xi) - 9955.55555555555*eta*xi*xi*xi + 591.111111111111*xi - 12320.0*xi*eta*eta + 9955.55555555555*xi*(eta*eta*eta) - 15.5555555555556 + 560.0*(eta*eta) + 8400.0*(eta*eta)*(xi*xi) - 3920.0*xi*xi + 248.888888888889*(eta*eta*eta) + 7217.77777777778*(xi*xi*xi) - 777.777777777778*eta*eta*eta*eta - 3888.88888888889*xi*xi*xi*xi, -1120.0*eta*xi + 12320.0*eta*(xi*xi) - 5600.0*eta*xi*xi*xi + 311.111111111111*xi - 373.333333333333*xi*eta*eta + 622.222222222222*xi*(eta*eta*eta) - 7466.66666666667*eta*eta*xi*xi - 1493.33333333333*xi*xi - 3733.33333333333*xi*xi*xi + 4977.77777777778*(xi*xi*xi*xi));
                      case 33:
                        return RealGradient(12693.3333333333*eta*xi - 622.222222222222*eta - 19040.0*eta*xi*xi + 7466.66666666667*eta*(xi*xi*xi) - 715.555555555556*xi - 33600.0*xi*eta*eta + 22400.0*xi*(eta*eta*eta) + 31.1111111111111 + 1680.0*(eta*eta) + 25200.0*(eta*eta)*(xi*xi) + 1680.0*(xi*xi) - 1120.0*eta*eta*eta - 1493.33333333333*xi*xi*xi + 466.666666666667*(xi*xi*xi*xi), -3360.0*eta*xi + 20160.0*eta*(xi*xi) - 16800.0*eta*xi*xi*xi + 684.444444444444*xi + 3360.0*xi*(eta*eta) - 16800.0*eta*eta*xi*xi - 4853.33333333333*xi*xi + 7840.0*(xi*xi*xi) - 3733.33333333333*xi*xi*xi*xi);
                      case 34:
                        return RealGradient(-9706.66666666667*eta*xi + 1057.77777777778*eta + 20160.0*eta*(xi*xi) - 11200.0*eta*xi*xi*xi + 684.444444444444*xi + 23520.0*xi*(eta*eta) - 14933.3333333333*xi*eta*eta*eta - 62.2222222222222 - 3920.0*eta*eta - 25200.0*eta*eta*xi*xi - 1680.0*xi*xi + 5226.66666666667*(eta*eta*eta) + 1120.0*(xi*xi*xi) - 2333.33333333333*eta*eta*eta*eta, 3360.0*eta*xi - 19040.0*eta*xi*xi + 16800.0*eta*(xi*xi*xi) - 715.555555555556*xi - 4480.0*xi*eta*eta + 1866.66666666667*xi*(eta*eta*eta) + 11200.0*(eta*eta)*(xi*xi) + 6346.66666666667*(xi*xi) - 11200.0*xi*xi*xi + 5600.0*(xi*xi*xi*xi));
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
      } // end case FIFTH

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

      // quartic Nedelec (first kind) shape function second derivatives
    case FOURTH:
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
                        return sign * RealGradient(-480.0*eta - 52.5*xi + 420.0*(eta + 1)*(xi + 1) - 787.5*(xi + 1)*(eta + 1)*(eta + 1) + 525.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 114.84375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 472.5 + 900.0*((eta + 1)*(eta + 1)) - 600.0*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 1:
                        return sign * RealGradient(191.111111111111*eta + 21.3888888888889*xi - 171.111111111111*(eta + 1)*(xi + 1) + 320.833333333333*(xi + 1)*((eta + 1)*(eta + 1)) - 213.888888888889*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 46.7881944444444*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 188.611111111111 - 358.333333333333*(eta + 1)*(eta + 1) + 238.888888888889*((eta + 1)*(eta + 1)*(eta + 1)) - 52.2569444444444*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 2:
                        return sign * RealGradient(-151.111111111111*eta - 21.3888888888889*xi + 171.111111111111*(eta + 1)*(xi + 1) - 320.833333333333*(xi + 1)*(eta + 1)*(eta + 1) + 213.888888888889*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 46.7881944444444*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 153.611111111111 + 283.333333333333*((eta + 1)*(eta + 1)) - 188.888888888889*(eta + 1)*(eta + 1)*(eta + 1) + 41.3194444444444*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 3:
                        return sign * RealGradient(360.0*eta + 52.5*xi - 420.0*(eta + 1)*(xi + 1) + 787.5*(xi + 1)*((eta + 1)*(eta + 1)) - 525.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 114.84375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 367.5 - 675.0*(eta + 1)*(eta + 1) + 450.0*((eta + 1)*(eta + 1)*(eta + 1)) - 98.4375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 4:
                        return sign * RealGradient(0, -450.0*eta - 360.0*xi + 1350.0*(eta + 1)*(xi + 1) - 787.5*(eta + 1)*(xi + 1)*(xi + 1) - 1350.0*(xi + 1)*(eta + 1)*(eta + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 690.0 + 450.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 210.0*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 131.25*(eta + 1)*(eta + 1)*(eta + 1));
                      case 5:
                        return sign * RealGradient(0, 158.333333333333*eta + 56.6666666666667*xi - 475.0*(eta + 1)*(xi + 1) + 277.083333333333*(eta + 1)*((xi + 1)*(xi + 1)) + 537.5*(xi + 1)*((eta + 1)*(eta + 1)) - 160.416666666667*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 196.111111111111 - 179.166666666667*(eta + 1)*(eta + 1) - 313.541666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 33.0555555555556*(xi + 1)*(xi + 1) + 93.5763888888889*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 53.4722222222222*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 6:
                        return sign * RealGradient(0, -83.3333333333333*eta - 26.6666666666667*xi + 250.0*(eta + 1)*(xi + 1) - 145.833333333333*(eta + 1)*(xi + 1)*(xi + 1) - 425.0*(xi + 1)*(eta + 1)*(eta + 1) + 160.416666666667*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 101.111111111111 + 141.666666666667*((eta + 1)*(eta + 1)) + 247.916666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 15.5555555555556*((xi + 1)*(xi + 1)) - 93.5763888888889*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 53.4722222222222*(eta + 1)*(eta + 1)*(eta + 1));
                      case 7:
                        return sign * RealGradient(0, 225.0*eta + 90.0*xi - 675.0*(eta + 1)*(xi + 1) + 393.75*(eta + 1)*((xi + 1)*(xi + 1)) + 1012.5*(xi + 1)*((eta + 1)*(eta + 1)) - 393.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 285.0 - 337.5*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 52.5*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 8:
                        return sign * RealGradient(-90.0*eta + 105.0*(eta + 1)*(xi + 1) - 393.75*(xi + 1)*(eta + 1)*(eta + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 114.84375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 90.0 + 337.5*((eta + 1)*(eta + 1)) - 337.5*(eta + 1)*(eta + 1)*(eta + 1) + 98.4375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 9:
                        return sign * RealGradient(37.7777777777778*eta - 42.7777777777778*(eta + 1)*(xi + 1) + 160.416666666667*(xi + 1)*((eta + 1)*(eta + 1)) - 160.416666666667*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 46.7881944444444*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 37.7777777777778 - 141.666666666667*(eta + 1)*(eta + 1) + 141.666666666667*((eta + 1)*(eta + 1)*(eta + 1)) - 41.3194444444444*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 10:
                        return sign * RealGradient(-47.7777777777778*eta + 42.7777777777778*(eta + 1)*(xi + 1) - 160.416666666667*(xi + 1)*(eta + 1)*(eta + 1) + 160.416666666667*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 46.7881944444444*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 47.7777777777778 + 179.166666666667*((eta + 1)*(eta + 1)) - 179.166666666667*(eta + 1)*(eta + 1)*(eta + 1) + 52.2569444444444*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 11:
                        return sign * RealGradient(120.0*eta - 105.0*(eta + 1)*(xi + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)) - 393.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 114.84375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 120.0 - 450.0*(eta + 1)*(eta + 1) + 450.0*((eta + 1)*(eta + 1)*(eta + 1)) - 131.25*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 12:
                        return sign * RealGradient(0, -450.0*eta - 120.0*xi + 900.0*(eta + 1)*(xi + 1) - 393.75*(eta + 1)*(xi + 1)*(xi + 1) - 1350.0*(xi + 1)*(eta + 1)*(eta + 1) + 525.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 510.0 + 675.0*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 52.5*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 262.5*(eta + 1)*(eta + 1)*(eta + 1));
                      case 13:
                        return sign * RealGradient(0, 166.666666666667*eta + 35.5555555555556*xi - 333.333333333333*(eta + 1)*(xi + 1) + 145.833333333333*(eta + 1)*((xi + 1)*(xi + 1)) + 566.666666666667*(xi + 1)*((eta + 1)*(eta + 1)) - 213.888888888889*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 184.444444444444 - 283.333333333333*(eta + 1)*(eta + 1) - 247.916666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 15.5555555555556*(xi + 1)*(xi + 1) + 93.5763888888889*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 106.944444444444*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 14:
                        return sign * RealGradient(0, -316.666666666667*eta - 75.5555555555556*xi + 633.333333333333*(eta + 1)*(xi + 1) - 277.083333333333*(eta + 1)*(xi + 1)*(xi + 1) - 716.666666666667*(xi + 1)*(eta + 1)*(eta + 1) + 213.888888888889*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 354.444444444444 + 358.333333333333*((eta + 1)*(eta + 1)) + 313.541666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 33.0555555555556*((xi + 1)*(xi + 1)) - 93.5763888888889*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 106.944444444444*(eta + 1)*(eta + 1)*(eta + 1));
                      case 15:
                        return sign * RealGradient(0, 900.0*eta + 480.0*xi - 1800.0*(eta + 1)*(xi + 1) + 787.5*(eta + 1)*((xi + 1)*(xi + 1)) + 1800.0*(xi + 1)*((eta + 1)*(eta + 1)) - 525.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1140.0 - 900.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 210.0*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 262.5*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 16:
                        return RealGradient(0, 870.0*eta + 480.0*xi - 1800.0*(eta + 1)*(xi + 1) + 787.5*(eta + 1)*((xi + 1)*(xi + 1)) + 1800.0*(xi + 1)*((eta + 1)*(eta + 1)) - 525.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1118.0 - 870.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 210.0*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 253.75*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 17:
                        return RealGradient(0, 420.0*eta + 360.0*xi - 1350.0*(eta + 1)*(xi + 1) + 787.5*(eta + 1)*((xi + 1)*(xi + 1)) + 1350.0*(xi + 1)*((eta + 1)*(eta + 1)) - 393.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 668.0 - 420.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 210.0*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 122.5*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 18:
                        return RealGradient(0, 60.0*eta + 44.0 - 60.0*(eta + 1)*(eta + 1) + 17.5*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 19:
                        return RealGradient(-390.0*eta + 341.25*(eta + 1)*(xi + 1) - 761.25*(xi + 1)*(eta + 1)*(eta + 1) + 525.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 114.84375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 390.0 + 870.0*((eta + 1)*(eta + 1)) - 600.0*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 20:
                        return RealGradient(-90.0*eta + 78.75*(eta + 1)*(xi + 1) - 367.5*(xi + 1)*(eta + 1)*(eta + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 114.84375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 90.0 + 420.0*((eta + 1)*(eta + 1)) - 450.0*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 21:
                        return RealGradient(-120.0*eta + 105.0*(eta + 1)*(xi + 1) - 52.5*(xi + 1)*(eta + 1)*(eta + 1) - 120.0 + 60.0*((eta + 1)*(eta + 1)), 0);
                      case 22:
                        return RealGradient(292.5*eta - 341.25*(eta + 1)*(xi + 1) + 761.25*(xi + 1)*((eta + 1)*(eta + 1)) - 525.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 114.84375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 292.5 - 652.5*(eta + 1)*(eta + 1) + 450.0*((eta + 1)*(eta + 1)*(eta + 1)) - 98.4375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 23:
                        return RealGradient(67.5*eta - 78.75*(eta + 1)*(xi + 1) + 367.5*(xi + 1)*((eta + 1)*(eta + 1)) - 393.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 114.84375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 67.5 - 315.0*(eta + 1)*(eta + 1) + 337.5*((eta + 1)*(eta + 1)*(eta + 1)) - 98.4375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 24:
                        return RealGradient(90.0*eta - 105.0*(eta + 1)*(xi + 1) + 52.5*(xi + 1)*((eta + 1)*(eta + 1)) + 90.0 - 45.0*(eta + 1)*(eta + 1), 0);
                      case 25:
                        return RealGradient(0, -435.0*eta - 120.0*xi + 900.0*(eta + 1)*(xi + 1) - 393.75*(eta + 1)*(xi + 1)*(xi + 1) - 1350.0*(xi + 1)*(eta + 1)*(eta + 1) + 525.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 497.0 + 652.5*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 52.5*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 253.75*(eta + 1)*(eta + 1)*(eta + 1));
                      case 26:
                        return RealGradient(0, -210.0*eta - 90.0*xi + 675.0*(eta + 1)*(xi + 1) - 393.75*(eta + 1)*(xi + 1)*(xi + 1) - 1012.5*(xi + 1)*(eta + 1)*(eta + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 272.0 + 315.0*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 52.5*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 122.5*(eta + 1)*(eta + 1)*(eta + 1));
                      case 27:
                        return RealGradient(0, -30.0*eta - 26.0 + 45.0*((eta + 1)*(eta + 1)) - 17.5*(eta + 1)*(eta + 1)*(eta + 1));
                      case 28:
                        return RealGradient(0, 0);
                      case 29:
                        return RealGradient(0, 0);
                      case 30:
                        return RealGradient(0, 42.75*eta + 120.0*xi - 90.0*(eta + 1)*(xi + 1) + 39.375*(eta + 1)*((xi + 1)*(xi + 1)) + 105.75 - 52.5*(xi + 1)*(xi + 1));
                      case 31:
                        return RealGradient(0, -42.75*eta - 60.0*xi + 90.0*(eta + 1)*(xi + 1) - 39.375*(eta + 1)*(xi + 1)*(xi + 1) - 74.25 + 26.25*((xi + 1)*(xi + 1)));
                      case 32:
                        return RealGradient(0, 20.25*eta + 90.0*xi - 67.5*(eta + 1)*(xi + 1) + 39.375*(eta + 1)*((xi + 1)*(xi + 1)) + 83.25 - 52.5*(xi + 1)*(xi + 1));
                      case 33:
                        return RealGradient(0, -20.25*eta - 45.0*xi + 67.5*(eta + 1)*(xi + 1) - 39.375*(eta + 1)*(xi + 1)*(xi + 1) - 51.75 + 26.25*((xi + 1)*(xi + 1)));
                      case 34:
                        return RealGradient(0, 0);
                      case 35:
                        return RealGradient(0, 0);
                      case 36:
                        return RealGradient(0, 4.5*eta - 1.5);
                      case 37:
                        return RealGradient(0, 0);
                      case 38:
                        return RealGradient(0, 0);
                      case 39:
                        return RealGradient(0, -4.5*eta - 1.5);
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
                        return sign * RealGradient(-900.0*eta - 480.0*xi + 450.0*(2*eta + 2)*(2*xi + 2) - 393.75*(2*eta + 2)*(xi + 1)*(xi + 1) - 900.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 262.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1140.0 + 900.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 210.0*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 262.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 1:
                        return sign * RealGradient(316.666666666667*eta + 191.111111111111*xi - 179.166666666667*(2*eta + 2)*(2*xi + 2) + 160.416666666667*(2*eta + 2)*((xi + 1)*(xi + 1)) + 358.333333333333*(2*xi + 2)*((eta + 1)*(eta + 1)) - 104.513888888889*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 423.333333333333 - 316.666666666667*(eta + 1)*(eta + 1) - 320.833333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 85.5555555555556*(xi + 1)*(xi + 1) + 93.5763888888889*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 92.3611111111111*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 2:
                        return sign * RealGradient(-166.666666666667*eta - 151.111111111111*xi + 141.666666666667*(2*eta + 2)*(2*xi + 2) - 160.416666666667*(2*eta + 2)*(xi + 1)*(xi + 1) - 283.333333333333*(2*xi + 2)*(eta + 1)*(eta + 1) + 82.6388888888889*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 273.333333333333 + 166.666666666667*((eta + 1)*(eta + 1)) + 320.833333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 85.5555555555556*((xi + 1)*(xi + 1)) - 93.5763888888889*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 48.6111111111111*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 3:
                        return sign * RealGradient(450.0*eta + 360.0*xi - 337.5*(2*eta + 2)*(2*xi + 2) + 393.75*(2*eta + 2)*((xi + 1)*(xi + 1)) + 675.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 196.875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 690.0 - 450.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 210.0*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 4:
                        return sign * RealGradient(0, -120.0*eta - 450.0*xi + 225.0*(2*eta + 2)*(2*xi + 2) - 675.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 262.5*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 196.875*(2*xi + 2)*(eta + 1)*(eta + 1) - 510.0 + 52.5*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 675.0*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(xi + 1));
                      case 5:
                        return sign * RealGradient(0, 47.7777777777778*eta + 158.333333333333*xi - 89.5833333333333*(2*eta + 2)*(2*xi + 2) + 268.75*(2*eta + 2)*((xi + 1)*(xi + 1)) - 104.513888888889*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 80.2083333333333*(2*xi + 2)*((eta + 1)*(eta + 1)) + 185.0 - 21.3888888888889*(eta + 1)*(eta + 1) - 240.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 93.5763888888889*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 237.5*(xi + 1)*(xi + 1) + 92.3611111111111*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 6:
                        return sign * RealGradient(0, -37.7777777777778*eta - 83.3333333333333*xi + 70.8333333333333*(2*eta + 2)*(2*xi + 2) - 212.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 82.6388888888889*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 80.2083333333333*(2*xi + 2)*(eta + 1)*(eta + 1) - 110.0 + 21.3888888888889*((eta + 1)*(eta + 1)) + 240.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.5763888888889*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 125.0*((xi + 1)*(xi + 1)) - 48.6111111111111*(xi + 1)*(xi + 1)*(xi + 1));
                      case 7:
                        return sign * RealGradient(0, 90.0*eta + 225.0*xi - 168.75*(2*eta + 2)*(2*xi + 2) + 506.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 196.875*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*(2*xi + 2)*((eta + 1)*(eta + 1)) + 285.0 - 52.5*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 337.5*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 8:
                        return sign * RealGradient(-225.0*eta - 90.0*xi + 168.75*(2*eta + 2)*(2*xi + 2) - 196.875*(2*eta + 2)*(xi + 1)*(xi + 1) - 506.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 196.875*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 285.0 + 337.5*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 52.5*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 131.25*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 9:
                        return sign * RealGradient(83.3333333333333*eta + 37.7777777777778*xi - 70.8333333333333*(2*eta + 2)*(2*xi + 2) + 80.2083333333333*(2*eta + 2)*((xi + 1)*(xi + 1)) + 212.5*(2*xi + 2)*((eta + 1)*(eta + 1)) - 82.6388888888889*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 110.0 - 125.0*(eta + 1)*(eta + 1) - 240.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 21.3888888888889*(xi + 1)*(xi + 1) + 93.5763888888889*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 48.6111111111111*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 10:
                        return sign * RealGradient(-158.333333333333*eta - 47.7777777777778*xi + 89.5833333333333*(2*eta + 2)*(2*xi + 2) - 80.2083333333333*(2*eta + 2)*(xi + 1)*(xi + 1) - 268.75*(2*xi + 2)*(eta + 1)*(eta + 1) + 104.513888888889*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 185.0 + 237.5*((eta + 1)*(eta + 1)) + 240.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 21.3888888888889*((xi + 1)*(xi + 1)) - 93.5763888888889*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 92.3611111111111*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 11:
                        return sign * RealGradient(450.0*eta + 120.0*xi - 225.0*(2*eta + 2)*(2*xi + 2) + 196.875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 675.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 262.5*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 510.0 - 675.0*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 52.5*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 262.5*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 12:
                        return sign * RealGradient(0, -360.0*eta - 450.0*xi + 337.5*(2*eta + 2)*(2*xi + 2) - 675.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 196.875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 393.75*(2*xi + 2)*(eta + 1)*(eta + 1) - 690.0 + 210.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 450.0*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 13:
                        return sign * RealGradient(0, 151.111111111111*eta + 166.666666666667*xi - 141.666666666667*(2*eta + 2)*(2*xi + 2) + 283.333333333333*(2*eta + 2)*((xi + 1)*(xi + 1)) - 82.6388888888889*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 160.416666666667*(2*xi + 2)*((eta + 1)*(eta + 1)) + 273.333333333333 - 85.5555555555556*(eta + 1)*(eta + 1) - 320.833333333333*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 93.5763888888889*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 166.666666666667*(xi + 1)*(xi + 1) + 48.6111111111111*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 14:
                        return sign * RealGradient(0, -191.111111111111*eta - 316.666666666667*xi + 179.166666666667*(2*eta + 2)*(2*xi + 2) - 358.333333333333*(2*eta + 2)*(xi + 1)*(xi + 1) + 104.513888888889*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 160.416666666667*(2*xi + 2)*(eta + 1)*(eta + 1) - 423.333333333333 + 85.5555555555556*((eta + 1)*(eta + 1)) + 320.833333333333*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.5763888888889*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 316.666666666667*((xi + 1)*(xi + 1)) - 92.3611111111111*(xi + 1)*(xi + 1)*(xi + 1));
                      case 15:
                        return sign * RealGradient(0, 480.0*eta + 900.0*xi - 450.0*(2*eta + 2)*(2*xi + 2) + 900.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 262.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 393.75*(2*xi + 2)*((eta + 1)*(eta + 1)) + 1140.0 - 210.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 900.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 16:
                        return RealGradient(0, 390.0*eta + 870.0*xi - 435.0*(2*eta + 2)*(2*xi + 2) + 900.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 262.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 380.625*(2*xi + 2)*((eta + 1)*(eta + 1)) + 1065.0 - 170.625*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 900.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 17:
                        return RealGradient(0, 90.0*eta + 420.0*xi - 210.0*(2*eta + 2)*(2*xi + 2) + 675.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 262.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 183.75*(2*xi + 2)*((eta + 1)*(eta + 1)) + 465.0 - 39.375*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 675.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 18:
                        return RealGradient(0, 120.0*eta + 60.0*xi - 30.0*(2*eta + 2)*(2*xi + 2) + 26.25*(2*xi + 2)*((eta + 1)*(eta + 1)) + 120.0 - 52.5*(eta + 1)*(eta + 1));
                      case 19:
                        return RealGradient(-870.0*eta - 390.0*xi + 435.0*(2*eta + 2)*(2*xi + 2) - 380.625*(2*eta + 2)*(xi + 1)*(xi + 1) - 900.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 262.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1065.0 + 900.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 170.625*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 262.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 20:
                        return RealGradient(-420.0*eta - 90.0*xi + 210.0*(2*eta + 2)*(2*xi + 2) - 183.75*(2*eta + 2)*(xi + 1)*(xi + 1) - 675.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 262.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 465.0 + 675.0*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 39.375*((xi + 1)*(xi + 1)) - 229.6875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 262.5*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 21:
                        return RealGradient(-60.0*eta - 120.0*xi + 30.0*(2*eta + 2)*(2*xi + 2) - 26.25*(2*eta + 2)*(xi + 1)*(xi + 1) - 120.0 + 52.5*((xi + 1)*(xi + 1)), 0);
                      case 22:
                        return RealGradient(435.0*eta + 292.5*xi - 326.25*(2*eta + 2)*(2*xi + 2) + 380.625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 675.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 196.875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 630.0 - 450.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 170.625*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 23:
                        return RealGradient(210.0*eta + 67.5*xi - 157.5*(2*eta + 2)*(2*xi + 2) + 183.75*(2*eta + 2)*((xi + 1)*(xi + 1)) + 506.25*(2*xi + 2)*((eta + 1)*(eta + 1)) - 196.875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 255.0 - 337.5*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 39.375*(xi + 1)*(xi + 1) + 229.6875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 24:
                        return RealGradient(30.0*eta + 90.0*xi - 22.5*(2*eta + 2)*(2*xi + 2) + 26.25*(2*eta + 2)*((xi + 1)*(xi + 1)) + 90.0 - 52.5*(xi + 1)*(xi + 1), 0);
                      case 25:
                        return RealGradient(0, -292.5*eta - 435.0*xi + 326.25*(2*eta + 2)*(2*xi + 2) - 675.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 196.875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 380.625*(2*xi + 2)*(eta + 1)*(eta + 1) - 630.0 + 170.625*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 450.0*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 26:
                        return RealGradient(0, -67.5*eta - 210.0*xi + 157.5*(2*eta + 2)*(2*xi + 2) - 506.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 196.875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 183.75*(2*xi + 2)*(eta + 1)*(eta + 1) - 255.0 + 39.375*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 337.5*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1));
                      case 27:
                        return RealGradient(0, -90.0*eta - 30.0*xi + 22.5*(2*eta + 2)*(2*xi + 2) - 26.25*(2*xi + 2)*(eta + 1)*(eta + 1) - 90.0 + 52.5*((eta + 1)*(eta + 1)));
                      case 28:
                        return RealGradient(42.75*eta + 33.75 - 45.0*(eta + 1)*(eta + 1) + 13.125*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 29:
                        return RealGradient(-42.75*eta - 33.75 + 45.0*((eta + 1)*(eta + 1)) - 13.125*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 30:
                        return RealGradient(0, 42.75*xi + 33.75 - 45.0*(xi + 1)*(xi + 1) + 13.125*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 31:
                        return RealGradient(0, -42.75*xi - 33.75 + 45.0*((xi + 1)*(xi + 1)) - 13.125*(xi + 1)*(xi + 1)*(xi + 1));
                      case 32:
                        return RealGradient(0, 20.25*xi + 18.75 - 33.75*(xi + 1)*(xi + 1) + 13.125*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 33:
                        return RealGradient(0, -20.25*xi - 18.75 + 33.75*((xi + 1)*(xi + 1)) - 13.125*(xi + 1)*(xi + 1)*(xi + 1));
                      case 34:
                        return RealGradient(20.25*eta + 18.75 - 33.75*(eta + 1)*(eta + 1) + 13.125*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 35:
                        return RealGradient(-20.25*eta - 18.75 + 33.75*((eta + 1)*(eta + 1)) - 13.125*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 36:
                        return RealGradient(0, 4.5*xi);
                      case 37:
                        return RealGradient(-4.5*eta, 0);
                      case 38:
                        return RealGradient(4.5*eta, 0);
                      case 39:
                        return RealGradient(0, -4.5*xi);
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
                        return sign * RealGradient(-480.0*eta - 900.0*xi + 1800.0*(eta + 1)*(xi + 1) - 1800.0*(eta + 1)*(xi + 1)*(xi + 1) + 525.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 787.5*(xi + 1)*(eta + 1)*(eta + 1) - 1140.0 + 210.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 900.0*((xi + 1)*(xi + 1)) - 262.5*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 1:
                        return sign * RealGradient(75.5555555555556*eta + 316.666666666667*xi - 633.333333333333*(eta + 1)*(xi + 1) + 716.666666666667*(eta + 1)*((xi + 1)*(xi + 1)) - 213.888888888889*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 277.083333333333*(xi + 1)*((eta + 1)*(eta + 1)) + 354.444444444444 - 33.0555555555556*(eta + 1)*(eta + 1) - 313.541666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 93.5763888888889*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 358.333333333333*(xi + 1)*(xi + 1) + 106.944444444444*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 2:
                        return sign * RealGradient(-35.5555555555556*eta - 166.666666666667*xi + 333.333333333333*(eta + 1)*(xi + 1) - 566.666666666667*(eta + 1)*(xi + 1)*(xi + 1) + 213.888888888889*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 145.833333333333*(xi + 1)*(eta + 1)*(eta + 1) - 184.444444444444 + 15.5555555555556*((eta + 1)*(eta + 1)) + 247.916666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.5763888888889*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 283.333333333333*((xi + 1)*(xi + 1)) - 106.944444444444*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 3:
                        return sign * RealGradient(120.0*eta + 450.0*xi - 900.0*(eta + 1)*(xi + 1) + 1350.0*(eta + 1)*((xi + 1)*(xi + 1)) - 525.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)) + 510.0 - 52.5*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 675.0*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 4:
                        return sign * RealGradient(0, -120.0*xi + 105.0*(eta + 1)*(xi + 1) - 393.75*(eta + 1)*(xi + 1)*(xi + 1) + 393.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 114.84375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 120.0 + 450.0*((xi + 1)*(xi + 1)) - 450.0*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 5:
                        return sign * RealGradient(0, 47.7777777777778*xi - 42.7777777777778*(eta + 1)*(xi + 1) + 160.416666666667*(eta + 1)*((xi + 1)*(xi + 1)) - 160.416666666667*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 46.7881944444444*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 47.7777777777778 - 179.166666666667*(xi + 1)*(xi + 1) + 179.166666666667*((xi + 1)*(xi + 1)*(xi + 1)) - 52.2569444444444*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 6:
                        return sign * RealGradient(0, -37.7777777777778*xi + 42.7777777777778*(eta + 1)*(xi + 1) - 160.416666666667*(eta + 1)*(xi + 1)*(xi + 1) + 160.416666666667*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 46.7881944444444*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 37.7777777777778 + 141.666666666667*((xi + 1)*(xi + 1)) - 141.666666666667*(xi + 1)*(xi + 1)*(xi + 1) + 41.3194444444444*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 7:
                        return sign * RealGradient(0, 90.0*xi - 105.0*(eta + 1)*(xi + 1) + 393.75*(eta + 1)*((xi + 1)*(xi + 1)) - 393.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 114.84375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 90.0 - 337.5*(xi + 1)*(xi + 1) + 337.5*((xi + 1)*(xi + 1)*(xi + 1)) - 98.4375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 8:
                        return sign * RealGradient(-90.0*eta - 225.0*xi + 675.0*(eta + 1)*(xi + 1) - 1012.5*(eta + 1)*(xi + 1)*(xi + 1) + 393.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 393.75*(xi + 1)*(eta + 1)*(eta + 1) - 285.0 + 52.5*((eta + 1)*(eta + 1)) + 590.625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 337.5*((xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 9:
                        return sign * RealGradient(26.6666666666667*eta + 83.3333333333333*xi - 250.0*(eta + 1)*(xi + 1) + 425.0*(eta + 1)*((xi + 1)*(xi + 1)) - 160.416666666667*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 145.833333333333*(xi + 1)*((eta + 1)*(eta + 1)) + 101.111111111111 - 15.5555555555556*(eta + 1)*(eta + 1) - 247.916666666667*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 93.5763888888889*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 141.666666666667*(xi + 1)*(xi + 1) + 53.4722222222222*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 10:
                        return sign * RealGradient(-56.6666666666667*eta - 158.333333333333*xi + 475.0*(eta + 1)*(xi + 1) - 537.5*(eta + 1)*(xi + 1)*(xi + 1) + 160.416666666667*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 277.083333333333*(xi + 1)*(eta + 1)*(eta + 1) - 196.111111111111 + 33.0555555555556*((eta + 1)*(eta + 1)) + 313.541666666667*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 93.5763888888889*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 179.166666666667*((xi + 1)*(xi + 1)) - 53.4722222222222*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 11:
                        return sign * RealGradient(360.0*eta + 450.0*xi - 1350.0*(eta + 1)*(xi + 1) + 1350.0*(eta + 1)*((xi + 1)*(xi + 1)) - 393.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 787.5*(xi + 1)*((eta + 1)*(eta + 1)) + 690.0 - 210.0*(eta + 1)*(eta + 1) - 787.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 450.0*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 12:
                        return sign * RealGradient(0, -52.5*eta - 360.0*xi + 420.0*(eta + 1)*(xi + 1) - 787.5*(eta + 1)*(xi + 1)*(xi + 1) + 525.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 114.84375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 367.5 + 675.0*((xi + 1)*(xi + 1)) - 450.0*(xi + 1)*(xi + 1)*(xi + 1) + 98.4375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 13:
                        return sign * RealGradient(0, 21.3888888888889*eta + 151.111111111111*xi - 171.111111111111*(eta + 1)*(xi + 1) + 320.833333333333*(eta + 1)*((xi + 1)*(xi + 1)) - 213.888888888889*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 46.7881944444444*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 153.611111111111 - 283.333333333333*(xi + 1)*(xi + 1) + 188.888888888889*((xi + 1)*(xi + 1)*(xi + 1)) - 41.3194444444444*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 14:
                        return sign * RealGradient(0, -21.3888888888889*eta - 191.111111111111*xi + 171.111111111111*(eta + 1)*(xi + 1) - 320.833333333333*(eta + 1)*(xi + 1)*(xi + 1) + 213.888888888889*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 46.7881944444444*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 188.611111111111 + 358.333333333333*((xi + 1)*(xi + 1)) - 238.888888888889*(xi + 1)*(xi + 1)*(xi + 1) + 52.2569444444444*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 15:
                        return sign * RealGradient(0, 52.5*eta + 480.0*xi - 420.0*(eta + 1)*(xi + 1) + 787.5*(eta + 1)*((xi + 1)*(xi + 1)) - 525.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 114.84375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 472.5 - 900.0*(xi + 1)*(xi + 1) + 600.0*((xi + 1)*(xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 16:
                        return RealGradient(0, 390.0*xi - 341.25*(eta + 1)*(xi + 1) + 761.25*(eta + 1)*((xi + 1)*(xi + 1)) - 525.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 114.84375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 390.0 - 870.0*(xi + 1)*(xi + 1) + 600.0*((xi + 1)*(xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 17:
                        return RealGradient(0, 90.0*xi - 78.75*(eta + 1)*(xi + 1) + 367.5*(eta + 1)*((xi + 1)*(xi + 1)) - 393.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 114.84375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 90.0 - 420.0*(xi + 1)*(xi + 1) + 450.0*((xi + 1)*(xi + 1)*(xi + 1)) - 131.25*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 18:
                        return RealGradient(0, 120.0*xi - 105.0*(eta + 1)*(xi + 1) + 52.5*(eta + 1)*((xi + 1)*(xi + 1)) + 120.0 - 60.0*(xi + 1)*(xi + 1));
                      case 19:
                        return RealGradient(-480.0*eta - 870.0*xi + 1800.0*(eta + 1)*(xi + 1) - 1800.0*(eta + 1)*(xi + 1)*(xi + 1) + 525.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 787.5*(xi + 1)*(eta + 1)*(eta + 1) - 1118.0 + 210.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 870.0*((xi + 1)*(xi + 1)) - 253.75*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 20:
                        return RealGradient(-360.0*eta - 420.0*xi + 1350.0*(eta + 1)*(xi + 1) - 1350.0*(eta + 1)*(xi + 1)*(xi + 1) + 393.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 787.5*(xi + 1)*(eta + 1)*(eta + 1) - 668.0 + 210.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 229.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 420.0*((xi + 1)*(xi + 1)) - 122.5*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 21:
                        return RealGradient(-60.0*xi - 44.0 + 60.0*((xi + 1)*(xi + 1)) - 17.5*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 22:
                        return RealGradient(120.0*eta + 435.0*xi - 900.0*(eta + 1)*(xi + 1) + 1350.0*(eta + 1)*((xi + 1)*(xi + 1)) - 525.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)) + 497.0 - 52.5*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 652.5*(xi + 1)*(xi + 1) + 253.75*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 23:
                        return RealGradient(90.0*eta + 210.0*xi - 675.0*(eta + 1)*(xi + 1) + 1012.5*(eta + 1)*((xi + 1)*(xi + 1)) - 393.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 393.75*(xi + 1)*((eta + 1)*(eta + 1)) + 272.0 - 52.5*(eta + 1)*(eta + 1) - 590.625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 229.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 315.0*(xi + 1)*(xi + 1) + 122.5*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 24:
                        return RealGradient(30.0*xi + 26.0 - 45.0*(xi + 1)*(xi + 1) + 17.5*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 25:
                        return RealGradient(0, -292.5*xi + 341.25*(eta + 1)*(xi + 1) - 761.25*(eta + 1)*(xi + 1)*(xi + 1) + 525.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 114.84375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 292.5 + 652.5*((xi + 1)*(xi + 1)) - 450.0*(xi + 1)*(xi + 1)*(xi + 1) + 98.4375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 26:
                        return RealGradient(0, -67.5*xi + 78.75*(eta + 1)*(xi + 1) - 367.5*(eta + 1)*(xi + 1)*(xi + 1) + 393.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 114.84375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 67.5 + 315.0*((xi + 1)*(xi + 1)) - 337.5*(xi + 1)*(xi + 1)*(xi + 1) + 98.4375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 27:
                        return RealGradient(0, -90.0*xi + 105.0*(eta + 1)*(xi + 1) - 52.5*(eta + 1)*(xi + 1)*(xi + 1) - 90.0 + 45.0*((xi + 1)*(xi + 1)));
                      case 28:
                        return RealGradient(120.0*eta + 42.75*xi - 90.0*(eta + 1)*(xi + 1) + 39.375*(xi + 1)*((eta + 1)*(eta + 1)) + 105.75 - 52.5*(eta + 1)*(eta + 1), 0);
                      case 29:
                        return RealGradient(-60.0*eta - 42.75*xi + 90.0*(eta + 1)*(xi + 1) - 39.375*(xi + 1)*(eta + 1)*(eta + 1) - 74.25 + 26.25*((eta + 1)*(eta + 1)), 0);
                      case 30:
                        return RealGradient(0, 0);
                      case 31:
                        return RealGradient(0, 0);
                      case 32:
                        return RealGradient(0, 0);
                      case 33:
                        return RealGradient(0, 0);
                      case 34:
                        return RealGradient(90.0*eta + 20.25*xi - 67.5*(eta + 1)*(xi + 1) + 39.375*(xi + 1)*((eta + 1)*(eta + 1)) + 83.25 - 52.5*(eta + 1)*(eta + 1), 0);
                      case 35:
                        return RealGradient(-45.0*eta - 20.25*xi + 67.5*(eta + 1)*(xi + 1) - 39.375*(xi + 1)*(eta + 1)*(eta + 1) - 51.75 + 26.25*((eta + 1)*(eta + 1)), 0);
                      case 36:
                        return RealGradient(0, 0);
                      case 37:
                        return RealGradient(1.5 - 4.5*xi, 0);
                      case 38:
                        return RealGradient(4.5*xi + 1.5, 0);
                      case 39:
                        return RealGradient(0, 0);
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
                        return sign * RealGradient(1344.0*eta*xi - 1680.0*eta - 840.0*xi + 480.0 + 1344.0*(eta*eta), -4032.0*eta*xi + 1680.0*eta + 2520.0*xi - 480.0 - 1344.0*eta*eta - 2688.0*xi*xi);
                      case 1:
                        return sign * RealGradient(-547.555555555556*eta*xi + 640.888888888889*eta + 342.222222222222*xi - 191.111111111111 - 497.777777777778*eta*eta, 1493.33333333333*eta*xi - 472.888888888889*eta - 1008.0*xi + 177.777777777778 + 199.111111111111*(eta*eta) + 1095.11111111111*(xi*xi));
                      case 2:
                        return sign * RealGradient(547.555555555556*eta*xi - 248.888888888889*eta - 342.222222222222*xi + 151.111111111111 + 49.7777777777778*(eta*eta), -149.333333333333*eta*xi - 87.1111111111111*eta + 840.0*xi - 113.777777777778 + 248.888888888889*(eta*eta) - 1095.11111111111*xi*xi);
                      case 3:
                        return sign * RealGradient(-1344.0*eta*xi + 504.0*eta + 840.0*xi - 360.0, -2016.0*xi + 288.0 + 2688.0*(xi*xi));
                      case 4:
                        return sign * RealGradient(168.0*eta*(3 - 8*xi), -2016.0*xi + 288.0 + 2688.0*(xi*xi));
                      case 5:
                        return sign * RealGradient(99.5555555555556*eta*(-6*eta - 4*xi + 3), 1792.0*eta*xi - 597.333333333333*eta - 821.333333333333*xi + 160.0 + 298.666666666667*(eta*eta) + 796.444444444444*(xi*xi));
                      case 6:
                        return sign * RealGradient(6.22222222222222*eta*(-48*eta - 8*xi + 15), 896.0*eta*xi - 522.666666666667*eta - 186.666666666667*xi + 64.0 + 597.333333333333*(eta*eta) + 99.5555555555556*(xi*xi));
                      case 7:
                        return sign * RealGradient(0, 0);
                      case 8:
                        return sign * RealGradient(0, 0);
                      case 9:
                        return sign * RealGradient(6.22222222222222*eta*(-40*eta + 8*xi + 7), 746.666666666667*eta*xi - 317.333333333333*eta + 12.4444444444444*xi + 23.1111111111111 - 49.7777777777778*eta*eta - 99.5555555555555*xi*xi);
                      case 10:
                        return sign * RealGradient(99.5555555555556*eta*(-2*eta + 4*xi - 1), 597.333333333333*eta*xi - 522.666666666667*eta + 771.555555555556*xi - 135.111111111111 + 497.777777777778*(eta*eta) - 796.444444444444*xi*xi);
                      case 11:
                        return sign * RealGradient(168.0*eta*(8*eta + 8*xi - 5), -4032.0*eta*xi + 2520.0*eta + 3360.0*xi - 960.0 - 1344.0*eta*eta - 2688.0*xi*xi);
                      case 12:
                        return RealGradient(2016.0*eta*(-4*eta - 2*xi + 3), 24192.0*eta*xi - 12096.0*eta - 9072.0*xi + 2160.0 + 12096.0*(eta*eta) + 8064.0*(xi*xi));
                      case 13:
                        return RealGradient(1008.0*eta*(12*eta + 16*xi - 9), -36288.0*eta*xi + 18144.0*eta + 36288.0*xi - 8640.0 - 8064.0*eta*eta - 32256.0*xi*xi);
                      case 14:
                        return RealGradient(1008.0*eta*(-8*eta - 12*xi + 9), 24192.0*eta*xi - 6048.0*eta - 21168.0*xi + 3456.0 + 24192.0*(xi*xi));
                      case 15:
                        return RealGradient(1008.0*eta*(4*eta + 16*xi - 7), -12096.0*eta*xi + 3024.0*eta + 28224.0*xi - 4608.0 - 32256.0*xi*xi);
                      case 16:
                        return RealGradient(0, -2016.0*eta + 144.0 + 4032.0*(eta*eta));
                      case 17:
                        return RealGradient(0, 4032.0*eta - 288.0 - 8064.0*eta*eta);
                      case 18:
                        return RealGradient(252.0*eta*(-24*eta - 12*xi + 13), 18144.0*eta*xi - 8064.0*eta - 6804.0*xi + 1548.0 + 7056.0*(eta*eta) + 6048.0*(xi*xi));
                      case 19:
                        return RealGradient(252.0*eta*(28*eta + 16*xi - 13), -21168.0*eta*xi + 9576.0*eta + 9324.0*xi - 2196.0 - 6048.0*eta*eta - 8064.0*xi*xi);
                      case 20:
                        return RealGradient(1008.0*eta*(1 - xi), 1008.0*eta - 1512.0*xi + 144.0 - 2016.0*eta*eta + 2016.0*(xi*xi));
                      case 21:
                        return RealGradient(252.0*eta*(-12*eta + 16*xi - 3), 9072.0*eta*xi - 5544.0*eta + 6804.0*xi - 936.0 + 4032.0*(eta*eta) - 8064.0*xi*xi);
                      case 22:
                        return RealGradient(1008.0*eta*(4*eta + 2*xi - 3), -12096.0*eta*xi + 4536.0*eta + 4536.0*xi - 972.0 - 3024.0*eta*eta - 4032.0*xi*xi);
                      case 23:
                        return RealGradient(2016.0*eta*(-eta - 2*xi + 1), 6048.0*eta*xi - 1512.0*eta - 8064.0*xi + 1440.0 + 8064.0*(xi*xi));
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
                        return sign * RealGradient(2688.0*eta*xi - 2520.0*eta - 1680.0*xi + 720.0 + 2016.0*(eta*eta) + 672.0*(xi*xi), -2688.0*eta*xi + 840.0*eta + 1680.0*xi - 240.0 - 672.0*eta*eta - 2016.0*xi*xi);
                      case 1:
                        return sign * RealGradient(-995.555555555556*eta*xi + 522.666666666667*eta + 640.888888888889*xi - 200.0 - 298.666666666667*eta*eta - 273.777777777778*xi*xi, 398.222222222222*eta*xi + 99.5555555555556*eta - 472.888888888889*xi + 19.5555555555556 - 199.111111111111*eta*eta + 746.666666666667*(xi*xi));
                      case 2:
                        return sign * RealGradient(99.5555555555556*eta*xi + 317.333333333333*eta - 248.888888888889*xi - 16.0 - 373.333333333333*eta*eta + 273.777777777778*(xi*xi), 497.777777777778*eta*xi - 43.5555555555556*eta - 87.1111111111111*xi + 12.4444444444444 - 24.8888888888889*eta*eta - 74.6666666666667*xi*xi);
                      case 3:
                        return sign * RealGradient(504.0*xi - 72.0 - 672.0*xi*xi, 0);
                      case 4:
                        return sign * RealGradient(504.0*xi - 72.0 - 672.0*xi*xi, 0);
                      case 5:
                        return sign * RealGradient(-1194.66666666667*eta*xi + 522.666666666667*eta + 298.666666666667*xi - 80.0 - 448.0*eta*eta - 199.111111111111*xi*xi, 597.333333333333*eta*xi - 93.3333333333333*eta - 597.333333333333*xi + 56.0 + 24.8888888888889*(eta*eta) + 896.0*(xi*xi));
                      case 6:
                        return sign * RealGradient(-597.333333333333*eta*xi + 597.333333333333*eta + 93.3333333333333*xi - 56.0 - 896.0*eta*eta - 24.8888888888889*xi*xi, 1194.66666666667*eta*xi - 298.666666666667*eta - 522.666666666667*xi + 80.0 + 199.111111111111*(eta*eta) + 448.0*(xi*xi));
                      case 7:
                        return sign * RealGradient(0, -504.0*eta + 72.0 + 672.0*(eta*eta));
                      case 8:
                        return sign * RealGradient(0, -504.0*eta + 72.0 + 672.0*(eta*eta));
                      case 9:
                        return sign * RealGradient(-497.777777777778*eta*xi + 87.1111111111111*eta + 43.5555555555556*xi - 12.4444444444444 + 74.6666666666667*(eta*eta) + 24.8888888888889*(xi*xi), -99.5555555555556*eta*xi + 248.888888888889*eta - 317.333333333333*xi + 16.0 - 273.777777777778*eta*eta + 373.333333333333*(xi*xi));
                      case 10:
                        return sign * RealGradient(-398.222222222222*eta*xi + 472.888888888889*eta - 99.5555555555556*xi - 19.5555555555556 - 746.666666666667*eta*eta + 199.111111111111*(xi*xi), 995.555555555556*eta*xi - 640.888888888889*eta - 522.666666666667*xi + 200.0 + 273.777777777778*(eta*eta) + 298.666666666667*(xi*xi));
                      case 11:
                        return sign * RealGradient(2688.0*eta*xi - 1680.0*eta - 840.0*xi + 240.0 + 2016.0*(eta*eta) + 672.0*(xi*xi), -2688.0*eta*xi + 1680.0*eta + 2520.0*xi - 720.0 - 672.0*eta*eta - 2016.0*xi*xi);
                      case 12:
                        return RealGradient(-16128.0*eta*xi + 18144.0*eta + 6048.0*xi - 3240.0 - 18144.0*eta*eta - 2016.0*xi*xi, 24192.0*eta*xi - 9072.0*eta - 12096.0*xi + 2160.0 + 8064.0*(eta*eta) + 12096.0*(xi*xi));
                      case 13:
                        return RealGradient(24192.0*eta*xi - 12096.0*eta - 9072.0*xi + 2160.0 + 12096.0*(eta*eta) + 8064.0*(xi*xi), -16128.0*eta*xi + 6048.0*eta + 18144.0*xi - 3240.0 - 2016.0*eta*eta - 18144.0*xi*xi);
                      case 14:
                        return RealGradient(-16128.0*eta*xi + 4032.0*eta + 9072.0*xi - 1944.0 - 6048.0*xi*xi, -6048.0*xi + 432.0 + 12096.0*(xi*xi));
                      case 15:
                        return RealGradient(8064.0*eta*xi - 2016.0*eta - 7056.0*xi + 1152.0 + 8064.0*(xi*xi), 3024.0*xi - 216.0 - 6048.0*xi*xi);
                      case 16:
                        return RealGradient(3024.0*eta - 216.0 - 6048.0*eta*eta, 8064.0*eta*xi - 7056.0*eta - 2016.0*xi + 1152.0 + 8064.0*(eta*eta));
                      case 17:
                        return RealGradient(-6048.0*eta + 432.0 + 12096.0*(eta*eta), -16128.0*eta*xi + 9072.0*eta + 4032.0*xi - 1944.0 - 6048.0*eta*eta);
                      case 18:
                        return RealGradient(-12096.0*eta*xi + 9576.0*eta + 3276.0*xi - 1332.0 - 10584.0*eta*eta - 1512.0*xi*xi, 14112.0*eta*xi - 3276.0*eta - 8064.0*xi + 1044.0 + 2016.0*(eta*eta) + 9072.0*(xi*xi));
                      case 19:
                        return RealGradient(14112.0*eta*xi - 8064.0*eta - 3276.0*xi + 1044.0 + 9072.0*(eta*eta) + 2016.0*(xi*xi), -12096.0*eta*xi + 3276.0*eta + 9576.0*xi - 1332.0 - 1512.0*eta*eta - 10584.0*xi*xi);
                      case 20:
                        return RealGradient(-1512.0*eta + 1008.0*xi - 216.0 + 3024.0*(eta*eta) - 504.0*xi*xi, -4032.0*eta*xi + 2016.0*eta + 1008.0*xi - 360.0 - 2016.0*eta*eta);
                      case 21:
                        return RealGradient(-6048.0*eta*xi + 4536.0*eta - 756.0*xi - 216.0 - 6048.0*eta*eta + 2016.0*(xi*xi), 8064.0*eta*xi - 3024.0*eta - 5544.0*xi + 1188.0 + 1008.0*(eta*eta) + 4536.0*(xi*xi));
                      case 22:
                        return RealGradient(8064.0*eta*xi - 5544.0*eta - 3024.0*xi + 1188.0 + 4536.0*(eta*eta) + 1008.0*(xi*xi), -6048.0*eta*xi - 756.0*eta + 4536.0*xi - 216.0 + 2016.0*(eta*eta) - 6048.0*xi*xi);
                      case 23:
                        return RealGradient(-4032.0*eta*xi + 1008.0*eta + 2016.0*xi - 360.0 - 2016.0*xi*xi, 1008.0*eta - 1512.0*xi - 216.0 - 504.0*eta*eta + 3024.0*(xi*xi));
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
                        return sign * RealGradient(4032.0*eta*xi - 3360.0*eta - 2520.0*xi + 960.0 + 2688.0*(eta*eta) + 1344.0*(xi*xi), 168.0*xi*(-8*eta - 8*xi + 5));
                      case 1:
                        return sign * RealGradient(-597.333333333333*eta*xi - 771.555555555556*eta + 522.666666666667*xi + 135.111111111111 + 796.444444444444*(eta*eta) - 497.777777777778*xi*xi, 99.5555555555556*xi*(-4*eta + 2*xi + 1));
                      case 2:
                        return sign * RealGradient(-746.666666666667*eta*xi - 12.4444444444444*eta + 317.333333333333*xi - 23.1111111111111 + 99.5555555555555*(eta*eta) + 49.7777777777778*(xi*xi), 6.22222222222222*xi*(-8*eta + 40*xi - 7));
                      case 3:
                        return sign * RealGradient(0, 0);
                      case 4:
                        return sign * RealGradient(0, 0);
                      case 5:
                        return sign * RealGradient(-896.0*eta*xi + 186.666666666667*eta + 522.666666666667*xi - 64.0 - 99.5555555555556*eta*eta - 597.333333333333*xi*xi, 6.22222222222222*xi*(8*eta + 48*xi - 15));
                      case 6:
                        return sign * RealGradient(-1792.0*eta*xi + 821.333333333333*eta + 597.333333333333*xi - 160.0 - 796.444444444444*eta*eta - 298.666666666667*xi*xi, 99.5555555555556*xi*(4*eta + 6*xi - 3));
                      case 7:
                        return sign * RealGradient(2016.0*eta - 288.0 - 2688.0*eta*eta, 168.0*xi*(8*eta - 3));
                      case 8:
                        return sign * RealGradient(2016.0*eta - 288.0 - 2688.0*eta*eta, 1344.0*eta*xi - 840.0*eta - 504.0*xi + 360.0);
                      case 9:
                        return sign * RealGradient(149.333333333333*eta*xi - 840.0*eta + 87.1111111111111*xi + 113.777777777778 + 1095.11111111111*(eta*eta) - 248.888888888889*xi*xi, -547.555555555556*eta*xi + 342.222222222222*eta + 248.888888888889*xi - 151.111111111111 - 49.7777777777778*xi*xi);
                      case 10:
                        return sign * RealGradient(-1493.33333333333*eta*xi + 1008.0*eta + 472.888888888889*xi - 177.777777777778 - 1095.11111111111*eta*eta - 199.111111111111*xi*xi, 547.555555555556*eta*xi - 342.222222222222*eta - 640.888888888889*xi + 191.111111111111 + 497.777777777778*(xi*xi));
                      case 11:
                        return sign * RealGradient(4032.0*eta*xi - 2520.0*eta - 1680.0*xi + 480.0 + 2688.0*(eta*eta) + 1344.0*(xi*xi), -1344.0*eta*xi + 840.0*eta + 1680.0*xi - 480.0 - 1344.0*xi*xi);
                      case 12:
                        return RealGradient(-36288.0*eta*xi + 36288.0*eta + 18144.0*xi - 8640.0 - 32256.0*eta*eta - 8064.0*xi*xi, 1008.0*xi*(16*eta + 12*xi - 9));
                      case 13:
                        return RealGradient(24192.0*eta*xi - 9072.0*eta - 12096.0*xi + 2160.0 + 8064.0*(eta*eta) + 12096.0*(xi*xi), 2016.0*xi*(-2*eta - 4*xi + 3));
                      case 14:
                        return RealGradient(4032.0*xi - 288.0 - 8064.0*xi*xi, 0);
                      case 15:
                        return RealGradient(-2016.0*xi + 144.0 + 4032.0*(xi*xi), 0);
                      case 16:
                        return RealGradient(-12096.0*eta*xi + 28224.0*eta + 3024.0*xi - 4608.0 - 32256.0*eta*eta, 1008.0*xi*(16*eta + 4*xi - 7));
                      case 17:
                        return RealGradient(24192.0*eta*xi - 21168.0*eta - 6048.0*xi + 3456.0 + 24192.0*(eta*eta), 1008.0*xi*(-12*eta - 8*xi + 9));
                      case 18:
                        return RealGradient(-21168.0*eta*xi + 9324.0*eta + 9576.0*xi - 2196.0 - 8064.0*eta*eta - 6048.0*xi*xi, 252.0*xi*(16*eta + 28*xi - 13));
                      case 19:
                        return RealGradient(18144.0*eta*xi - 6804.0*eta - 8064.0*xi + 1548.0 + 6048.0*(eta*eta) + 7056.0*(xi*xi), 252.0*xi*(-12*eta - 24*xi + 13));
                      case 20:
                        return RealGradient(6048.0*eta*xi - 8064.0*eta - 1512.0*xi + 1440.0 + 8064.0*(eta*eta), 2016.0*xi*(-2*eta - xi + 1));
                      case 21:
                        return RealGradient(-12096.0*eta*xi + 4536.0*eta + 4536.0*xi - 972.0 - 4032.0*eta*eta - 3024.0*xi*xi, 1008.0*xi*(2*eta + 4*xi - 3));
                      case 22:
                        return RealGradient(9072.0*eta*xi + 6804.0*eta - 5544.0*xi - 936.0 - 8064.0*eta*eta + 4032.0*(xi*xi), 252.0*xi*(16*eta - 12*xi - 3));
                      case 23:
                        return RealGradient(-1512.0*eta + 1008.0*xi + 144.0 + 2016.0*(eta*eta) - 2016.0*xi*xi, 1008.0*xi*(1 - eta));
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
      } // end case FOURTH

      // quintic Nedelec (first kind) shape function second derivatives
    case FIFTH:
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
                        return sign * RealGradient(-3281.25*eta - 525.0*xi + 6562.5*(eta + 1)*(xi + 1) - 2953.125*(eta + 1)*(xi + 1)*(xi + 1) - 19687.5*(xi + 1)*(eta + 1)*(eta + 1) + 22968.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 11484.375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 3543.75 + 9843.75*((eta + 1)*(eta + 1)) + 8859.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 236.25*((xi + 1)*(xi + 1)) - 10335.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 930.234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 11484.375*(eta + 1)*(eta + 1)*(eta + 1) + 5742.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1033.59375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 1:
                        return sign * RealGradient(1030.517578125*eta + 159.55078125*xi - 1994.384765625*(eta + 1)*(xi + 1) + 853.6376953125*(eta + 1)*((xi + 1)*(xi + 1)) + 5983.154296875*(xi + 1)*((eta + 1)*(eta + 1)) - 6980.3466796875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3490.17333984375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 628.231201171875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1107.626953125 - 3091.552734375*(eta + 1)*(eta + 1) - 2560.9130859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 68.291015625*(xi + 1)*(xi + 1) + 2987.73193359375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1493.86596679688*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 268.895874023438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 3606.8115234375*((eta + 1)*(eta + 1)*(eta + 1)) - 1803.40576171875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 324.613037109375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 2:
                        return sign * RealGradient(-902.34375*eta - 177.1875*xi + 2214.84375*(eta + 1)*(xi + 1) - 1107.421875*(eta + 1)*(xi + 1)*(xi + 1) - 6644.53125*(xi + 1)*(eta + 1)*(eta + 1) + 7751.953125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 3875.9765625*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 697.67578125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1007.34375 + 2707.03125*((eta + 1)*(eta + 1)) + 3322.265625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 88.59375*((xi + 1)*(xi + 1)) - 3875.9765625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1937.98828125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 348.837890625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3158.203125*(eta + 1)*(eta + 1)*(eta + 1) + 1579.1015625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 284.23828125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 3:
                        return sign * RealGradient(456.298828125*eta + 113.61328125*xi - 1420.166015625*(eta + 1)*(xi + 1) + 853.6376953125*(eta + 1)*((xi + 1)*(xi + 1)) + 4260.498046875*(xi + 1)*((eta + 1)*(eta + 1)) - 4970.5810546875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2485.29052734375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 447.352294921875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 533.408203125 - 1368.896484375*(eta + 1)*(eta + 1) - 2560.9130859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 68.291015625*(xi + 1)*(xi + 1) + 2987.73193359375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1493.86596679688*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 268.895874023438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1597.0458984375*((eta + 1)*(eta + 1)*(eta + 1)) - 798.52294921875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 143.734130859375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 4:
                        return sign * RealGradient(-1968.75*eta - 420.0*xi + 5250.0*(eta + 1)*(xi + 1) - 2953.125*(eta + 1)*(xi + 1)*(xi + 1) - 15750.0*(xi + 1)*(eta + 1)*(eta + 1) + 18375.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 9187.5*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1653.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2231.25 + 5906.25*((eta + 1)*(eta + 1)) + 8859.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 236.25*((xi + 1)*(xi + 1)) - 10335.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 930.234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1) + 3445.3125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 620.15625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 5:
                        return sign * RealGradient(0, 2250.0*eta + 1968.75*xi - 11812.5*(eta + 1)*(xi + 1) + 15750.0*(eta + 1)*((xi + 1)*(xi + 1)) - 5906.25*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 20671.875*(xi + 1)*((eta + 1)*(eta + 1)) - 13781.25*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3100.78125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 3843.75 - 3937.5*(eta + 1)*(eta + 1) - 27562.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 2625.0*(xi + 1)*(xi + 1) + 18375.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4134.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2625.0*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 984.375*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 590.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 6:
                        return sign * RealGradient(0, -671.484375*eta - 298.388671875*xi + 3525.29296875*(eta + 1)*(xi + 1) - 4700.390625*(eta + 1)*(xi + 1)*(xi + 1) + 1762.646484375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 6492.2607421875*(xi + 1)*(eta + 1)*(eta + 1) + 4188.2080078125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 896.319580078125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 913.037109375 + 1236.62109375*((eta + 1)*(eta + 1)) + 8656.34765625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3246.13037109375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 397.8515625*((xi + 1)*(xi + 1)) - 5584.27734375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1195.0927734375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 797.75390625*(eta + 1)*(eta + 1)*(eta + 1) + 2094.10400390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 149.1943359375*(xi + 1)*(xi + 1)*(xi + 1) - 448.159790039063*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 170.7275390625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 7:
                        return sign * RealGradient(0, 393.75*eta + 147.65625*xi - 2067.1875*(eta + 1)*(xi + 1) + 2756.25*(eta + 1)*((xi + 1)*(xi + 1)) - 1033.59375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5684.765625*(xi + 1)*((eta + 1)*(eta + 1)) - 4651.171875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1162.79296875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 513.28125 - 1082.8125*(eta + 1)*(eta + 1) - 7579.6875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2842.3828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1) + 6201.5625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 885.9375*((eta + 1)*(eta + 1)*(eta + 1)) - 2325.5859375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*((xi + 1)*(xi + 1)*(xi + 1)) + 581.396484375*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 221.484375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 8:
                        return sign * RealGradient(0, -165.234375*eta - 52.294921875*xi + 867.48046875*(eta + 1)*(xi + 1) - 1156.640625*(eta + 1)*(xi + 1)*(xi + 1) + 433.740234375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 2874.6826171875*(xi + 1)*(eta + 1)*(eta + 1) + 2982.3486328125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 896.319580078125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 207.568359375 + 547.55859375*((eta + 1)*(eta + 1)) + 3832.91015625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1437.34130859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 69.7265625*((xi + 1)*(xi + 1)) - 3976.46484375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1195.0927734375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 568.06640625*(eta + 1)*(eta + 1)*(eta + 1) + 1491.17431640625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 26.1474609375*(xi + 1)*(xi + 1)*(xi + 1) - 448.159790039063*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 170.7275390625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 9:
                        return sign * RealGradient(0, 900.0*eta + 393.75*xi - 4725.0*(eta + 1)*(xi + 1) + 6300.0*(eta + 1)*((xi + 1)*(xi + 1)) - 2362.5*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 12403.125*(xi + 1)*((eta + 1)*(eta + 1)) - 11025.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3100.78125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1218.75 - 2362.5*(eta + 1)*(eta + 1) - 16537.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 525.0*(xi + 1)*(xi + 1) + 14700.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4134.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2100.0*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 590.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 10:
                        return sign * RealGradient(-393.75*eta + 1050.0*(eta + 1)*(xi + 1) - 590.625*(eta + 1)*(xi + 1)*(xi + 1) - 6300.0*(xi + 1)*(eta + 1)*(eta + 1) + 11025.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 7350.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1653.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 393.75 + 2362.5*((eta + 1)*(eta + 1)) + 3543.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6201.5625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 930.234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 4134.375*(eta + 1)*(eta + 1)*(eta + 1) + 2756.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 620.15625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 11:
                        return sign * RealGradient(91.259765625*eta - 284.033203125*(eta + 1)*(xi + 1) + 170.7275390625*(eta + 1)*((xi + 1)*(xi + 1)) + 1704.19921875*(xi + 1)*((eta + 1)*(eta + 1)) - 2982.3486328125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1988.232421875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 447.352294921875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 91.259765625 - 547.55859375*(eta + 1)*(eta + 1) - 1024.365234375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1792.63916015625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1195.0927734375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 268.895874023438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 958.2275390625*((eta + 1)*(eta + 1)*(eta + 1)) - 638.818359375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 143.734130859375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 12:
                        return sign * RealGradient(-180.46875*eta + 442.96875*(eta + 1)*(xi + 1) - 221.484375*(eta + 1)*(xi + 1)*(xi + 1) - 2657.8125*(xi + 1)*(eta + 1)*(eta + 1) + 4651.171875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 3100.78125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 697.67578125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 180.46875 + 1082.8125*((eta + 1)*(eta + 1)) + 1328.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2325.5859375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1550.390625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 348.837890625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1894.921875*(eta + 1)*(eta + 1)*(eta + 1) + 1263.28125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 284.23828125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 13:
                        return sign * RealGradient(206.103515625*eta - 398.876953125*(eta + 1)*(xi + 1) + 170.7275390625*(eta + 1)*((xi + 1)*(xi + 1)) + 2393.26171875*(xi + 1)*((eta + 1)*(eta + 1)) - 4188.2080078125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2792.138671875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 628.231201171875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 206.103515625 - 1236.62109375*(eta + 1)*(eta + 1) - 1024.365234375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1792.63916015625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1195.0927734375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 268.895874023438*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2164.0869140625*((eta + 1)*(eta + 1)*(eta + 1)) - 1442.724609375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 324.613037109375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 14:
                        return sign * RealGradient(-656.25*eta + 1312.5*(eta + 1)*(xi + 1) - 590.625*(eta + 1)*(xi + 1)*(xi + 1) - 7875.0*(xi + 1)*(eta + 1)*(eta + 1) + 13781.25*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 9187.5*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 656.25 + 3937.5*((eta + 1)*(eta + 1)) + 3543.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6201.5625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 930.234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1) + 4593.75*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1033.59375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 15:
                        return sign * RealGradient(0, 2250.0*eta + 656.25*xi - 7875.0*(eta + 1)*(xi + 1) + 7875.0*(eta + 1)*((xi + 1)*(xi + 1)) - 2362.5*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 20671.875*(xi + 1)*((eta + 1)*(eta + 1)) - 18375.0*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2718.75 - 5906.25*(eta + 1)*(eta + 1) - 20671.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 656.25*(xi + 1)*(xi + 1) + 18375.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5167.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5250.0*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1476.5625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 16:
                        return sign * RealGradient(0, -413.0859375*eta - 87.158203125*xi + 1445.80078125*(eta + 1)*(xi + 1) - 1445.80078125*(eta + 1)*(xi + 1)*(xi + 1) + 433.740234375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 4791.1376953125*(xi + 1)*(eta + 1)*(eta + 1) + 4970.5810546875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 1493.86596679688*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 475.341796875 + 1368.896484375*((eta + 1)*(eta + 1)) + 4791.1376953125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1437.34130859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 87.158203125*((xi + 1)*(xi + 1)) - 4970.5810546875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1493.86596679688*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1420.166015625*(eta + 1)*(eta + 1)*(eta + 1) + 1491.17431640625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 26.1474609375*(xi + 1)*(xi + 1)*(xi + 1) - 448.159790039063*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 426.81884765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 17:
                        return sign * RealGradient(0, 984.375*eta + 246.09375*xi - 3445.3125*(eta + 1)*(xi + 1) + 3445.3125*(eta + 1)*((xi + 1)*(xi + 1)) - 1033.59375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9474.609375*(xi + 1)*((eta + 1)*(eta + 1)) - 7751.953125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1937.98828125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1160.15625 - 2707.03125*(eta + 1)*(eta + 1) - 9474.609375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2842.3828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 246.09375*(xi + 1)*(xi + 1) + 7751.953125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1937.98828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2214.84375*((eta + 1)*(eta + 1)*(eta + 1)) - 2325.5859375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*((xi + 1)*(xi + 1)*(xi + 1)) + 581.396484375*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 553.7109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 18:
                        return sign * RealGradient(0, -1678.7109375*eta - 497.314453125*xi + 5875.48828125*(eta + 1)*(xi + 1) - 5875.48828125*(eta + 1)*(xi + 1)*(xi + 1) + 1762.646484375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 10820.4345703125*(xi + 1)*(eta + 1)*(eta + 1) + 6980.3466796875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 1493.86596679688*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 2033.935546875 + 3091.552734375*((eta + 1)*(eta + 1)) + 10820.4345703125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3246.13037109375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 497.314453125*((xi + 1)*(xi + 1)) - 6980.3466796875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1493.86596679688*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1994.384765625*(eta + 1)*(eta + 1)*(eta + 1) + 2094.10400390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 149.1943359375*(xi + 1)*(xi + 1)*(xi + 1) - 448.159790039063*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 426.81884765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 19:
                        return sign * RealGradient(0, 5625.0*eta + 3281.25*xi - 19687.5*(eta + 1)*(xi + 1) + 19687.5*(eta + 1)*((xi + 1)*(xi + 1)) - 5906.25*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 34453.125*(xi + 1)*((eta + 1)*(eta + 1)) - 22968.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 7968.75 - 9843.75*(eta + 1)*(eta + 1) - 34453.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 3281.25*(xi + 1)*(xi + 1) + 22968.75*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5167.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 6562.5*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 984.375*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1476.5625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 20:
                        return RealGradient(0, 5287.5*eta + 3239.0625*xi - 19434.375*(eta + 1)*(xi + 1) + 19687.5*(eta + 1)*((xi + 1)*(xi + 1)) - 5906.25*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 34010.15625*(xi + 1)*((eta + 1)*(eta + 1)) - 22673.4375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5101.5234375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 7645.3125 - 9253.125*(eta + 1)*(eta + 1) - 34453.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 3281.25*(xi + 1)*(xi + 1) + 22968.75*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5167.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 6168.75*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 984.375*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1387.96875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 21:
                        return RealGradient(0, -2081.25*eta - 1926.5625*xi + 11559.375*(eta + 1)*(xi + 1) - 15750.0*(eta + 1)*(xi + 1)*(xi + 1) + 5906.25*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 20228.90625*(xi + 1)*(eta + 1)*(eta + 1) + 13485.9375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 3034.3359375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3660.9375 + 3642.1875*((eta + 1)*(eta + 1)) + 27562.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10335.9375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2625.0*((xi + 1)*(xi + 1)) - 18375.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2428.125*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 984.375*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 546.328125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 22:
                        return RealGradient(0, 843.75*eta + 126.5625*xi - 759.375*(eta + 1)*(xi + 1) + 1328.90625*(xi + 1)*((eta + 1)*(eta + 1)) - 885.9375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 199.3359375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 829.6875 - 1476.5625*(eta + 1)*(eta + 1) + 984.375*((eta + 1)*(eta + 1)*(eta + 1)) - 221.484375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 23:
                        return RealGradient(0, -675.0*eta - 126.5625*xi + 759.375*(eta + 1)*(xi + 1) - 1328.90625*(xi + 1)*(eta + 1)*(eta + 1) + 885.9375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 199.3359375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 689.0625 + 1181.25*((eta + 1)*(eta + 1)) - 787.5*(eta + 1)*(eta + 1)*(eta + 1) + 177.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 24:
                        return RealGradient(-2559.375*eta + 5118.75*(eta + 1)*(xi + 1) - 2303.4375*(eta + 1)*(xi + 1)*(xi + 1) - 18506.25*(xi + 1)*(eta + 1)*(eta + 1) + 22673.4375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 11484.375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2559.375 + 9253.125*((eta + 1)*(eta + 1)) + 8327.8125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10203.046875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 930.234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 11336.71875*(eta + 1)*(eta + 1)*(eta + 1) + 5742.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1033.59375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 25:
                        return RealGradient(525.0*eta - 1050.0*(eta + 1)*(xi + 1) + 472.5*(eta + 1)*((xi + 1)*(xi + 1)) + 7284.375*(xi + 1)*((eta + 1)*(eta + 1)) - 13485.9375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 9187.5*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2067.1875*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 525.0 - 3642.1875*(eta + 1)*(eta + 1) - 3277.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6068.671875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4134.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 930.234375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 6742.96875*((eta + 1)*(eta + 1)*(eta + 1)) - 4593.75*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1033.59375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 26:
                        return RealGradient(-1181.25*eta + 2362.5*(eta + 1)*(xi + 1) - 1063.125*(eta + 1)*(xi + 1)*(xi + 1) - 2953.125*(xi + 1)*(eta + 1)*(eta + 1) + 885.9375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 1181.25 + 1476.5625*((eta + 1)*(eta + 1)) + 1328.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 398.671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 442.96875*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 27:
                        return RealGradient(590.625*eta - 1181.25*(eta + 1)*(xi + 1) + 531.5625*(eta + 1)*((xi + 1)*(xi + 1)) + 2362.5*(xi + 1)*((eta + 1)*(eta + 1)) - 885.9375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625 - 1181.25*(eta + 1)*(eta + 1) - 1063.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.671875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 442.96875*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 28:
                        return RealGradient(-1535.625*eta + 4095.0*(eta + 1)*(xi + 1) - 2303.4375*(eta + 1)*(xi + 1)*(xi + 1) - 14805.0*(xi + 1)*(eta + 1)*(eta + 1) + 18138.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 9187.5*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1653.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1535.625 + 5551.875*((eta + 1)*(eta + 1)) + 8327.8125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10203.046875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 930.234375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 6802.03125*(eta + 1)*(eta + 1)*(eta + 1) + 3445.3125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 620.15625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 29:
                        return RealGradient(315.0*eta - 840.0*(eta + 1)*(xi + 1) + 472.5*(eta + 1)*((xi + 1)*(xi + 1)) + 5827.5*(xi + 1)*((eta + 1)*(eta + 1)) - 10788.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 7350.0*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1653.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 315.0 - 2185.3125*(eta + 1)*(eta + 1) - 3277.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6068.671875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4134.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 930.234375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 4045.78125*((eta + 1)*(eta + 1)*(eta + 1)) - 2756.25*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 620.15625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 30:
                        return RealGradient(-708.75*eta + 1890.0*(eta + 1)*(xi + 1) - 1063.125*(eta + 1)*(xi + 1)*(xi + 1) - 2362.5*(xi + 1)*(eta + 1)*(eta + 1) + 708.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 708.75 + 885.9375*((eta + 1)*(eta + 1)) + 1328.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 398.671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 265.78125*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 31:
                        return RealGradient(354.375*eta - 945.0*(eta + 1)*(xi + 1) + 531.5625*(eta + 1)*((xi + 1)*(xi + 1)) + 1890.0*(xi + 1)*((eta + 1)*(eta + 1)) - 708.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 354.375 - 708.75*(eta + 1)*(eta + 1) - 1063.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.671875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 265.78125*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 32:
                        return RealGradient(0, 2115.0*eta + 647.8125*xi - 7773.75*(eta + 1)*(xi + 1) + 7875.0*(eta + 1)*((xi + 1)*(xi + 1)) - 2362.5*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 20406.09375*(xi + 1)*((eta + 1)*(eta + 1)) - 18138.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5101.5234375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2586.5625 - 5551.875*(eta + 1)*(eta + 1) - 20671.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 656.25*(xi + 1)*(xi + 1) + 18375.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5167.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4935.0*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1387.96875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 33:
                        return RealGradient(0, -832.5*eta - 385.3125*xi + 4623.75*(eta + 1)*(xi + 1) - 6300.0*(eta + 1)*(xi + 1)*(xi + 1) + 2362.5*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 12137.34375*(xi + 1)*(eta + 1)*(eta + 1) + 10788.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 3034.3359375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1148.4375 + 2185.3125*((eta + 1)*(eta + 1)) + 16537.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6201.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 525.0*((xi + 1)*(xi + 1)) - 14700.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1942.5*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 546.328125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 34:
                        return RealGradient(0, 337.5*eta + 25.3125*xi - 303.75*(eta + 1)*(xi + 1) + 797.34375*(xi + 1)*((eta + 1)*(eta + 1)) - 708.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 199.3359375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 334.6875 - 885.9375*(eta + 1)*(eta + 1) + 787.5*((eta + 1)*(eta + 1)*(eta + 1)) - 221.484375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1));
                      case 35:
                        return RealGradient(0, -270.0*eta - 25.3125*xi + 303.75*(eta + 1)*(xi + 1) - 797.34375*(xi + 1)*(eta + 1)*(eta + 1) + 708.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 199.3359375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 272.8125 + 708.75*((eta + 1)*(eta + 1)) - 630.0*(eta + 1)*(eta + 1)*(eta + 1) + 177.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)));
                      case 36:
                        return RealGradient(60.0*eta + 60.0 - 247.5*(eta + 1)*(eta + 1) + 318.75*((eta + 1)*(eta + 1)*(eta + 1)) - 164.0625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 29.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 37:
                        return RealGradient(60.0*eta + 60.0 - 247.5*(eta + 1)*(eta + 1) + 318.75*((eta + 1)*(eta + 1)*(eta + 1)) - 164.0625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 29.53125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 38:
                        return RealGradient(-30.0*eta - 30.0 + 123.75*((eta + 1)*(eta + 1)) - 159.375*(eta + 1)*(eta + 1)*(eta + 1) + 82.03125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 14.765625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 39:
                        return RealGradient(0, 594.0*eta + 1147.5*xi - 2295.0*(eta + 1)*(xi + 1) + 2362.5*(eta + 1)*((xi + 1)*(xi + 1)) - 708.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 956.25*(xi + 1)*((eta + 1)*(eta + 1)) + 1444.5 - 247.5*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 295.3125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1181.25*(xi + 1)*(xi + 1) + 354.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 40:
                        return RealGradient(0, 396.0*eta + 382.5*xi - 1530.0*(eta + 1)*(xi + 1) + 1575.0*(eta + 1)*((xi + 1)*(xi + 1)) - 472.5*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 956.25*(xi + 1)*((eta + 1)*(eta + 1)) + 679.5 - 247.5*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 295.3125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 393.75*(xi + 1)*(xi + 1) + 118.125*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 41:
                        return RealGradient(0, -247.5*eta - 191.25*xi + 956.25*(eta + 1)*(xi + 1) - 984.375*(eta + 1)*(xi + 1)*(xi + 1) + 295.3125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 478.125*(xi + 1)*(eta + 1)*(eta + 1) - 389.25 + 123.75*((eta + 1)*(eta + 1)) + 492.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 147.65625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 196.875*((xi + 1)*(xi + 1)) - 59.0625*(xi + 1)*(xi + 1)*(xi + 1));
                      case 42:
                        return RealGradient(0, -216.0*eta - 675.0*xi + 1350.0*(eta + 1)*(xi + 1) - 1890.0*(eta + 1)*(xi + 1)*(xi + 1) + 708.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 562.5*(xi + 1)*(eta + 1)*(eta + 1) - 783.0 + 90.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 295.3125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 945.0*((xi + 1)*(xi + 1)) - 354.375*(xi + 1)*(xi + 1)*(xi + 1));
                      case 43:
                        return RealGradient(0, -144.0*eta - 225.0*xi + 900.0*(eta + 1)*(xi + 1) - 1260.0*(eta + 1)*(xi + 1)*(xi + 1) + 472.5*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 562.5*(xi + 1)*(eta + 1)*(eta + 1) - 333.0 + 90.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 295.3125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 315.0*((xi + 1)*(xi + 1)) - 118.125*(xi + 1)*(xi + 1)*(xi + 1));
                      case 44:
                        return RealGradient(0, 90.0*eta + 112.5*xi - 562.5*(eta + 1)*(xi + 1) + 787.5*(eta + 1)*((xi + 1)*(xi + 1)) - 295.3125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 281.25*(xi + 1)*((eta + 1)*(eta + 1)) + 184.5 - 45.0*(eta + 1)*(eta + 1) - 393.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 147.65625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 157.5*(xi + 1)*(xi + 1) + 59.0625*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 45:
                        return RealGradient(-7.5*eta - 7.5 + 90.0*((eta + 1)*(eta + 1)) - 187.5*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 29.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 46:
                        return RealGradient(-7.5*eta - 7.5 + 90.0*((eta + 1)*(eta + 1)) - 187.5*(eta + 1)*(eta + 1)*(eta + 1) + 131.25*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 29.53125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 47:
                        return RealGradient(3.75*eta + 3.75 - 45.0*(eta + 1)*(eta + 1) + 93.75*((eta + 1)*(eta + 1)*(eta + 1)) - 65.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 14.765625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 48:
                        return RealGradient(0, 81.0*eta + 33.75*xi - 67.5*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 74.25 - 33.75*(eta + 1)*(eta + 1));
                      case 49:
                        return RealGradient(0, -54.0*eta - 33.75*xi + 67.5*(eta + 1)*(xi + 1) - 28.125*(xi + 1)*(eta + 1)*(eta + 1) - 60.75 + 22.5*((eta + 1)*(eta + 1)));
                      case 50:
                        return RealGradient(-30.0*eta - 30.0 + 33.75*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 51:
                        return RealGradient(7.5*eta + 7.5 - 22.5*(eta + 1)*(eta + 1) + 9.375*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 52:
                        return RealGradient(-30.0*eta - 30.0 + 33.75*((eta + 1)*(eta + 1)) - 9.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 53:
                        return RealGradient(7.5*eta + 7.5 - 22.5*(eta + 1)*(eta + 1) + 9.375*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 54:
                        return RealGradient(0, 54.0*eta + 11.25*xi - 45.0*(eta + 1)*(xi + 1) + 28.125*(xi + 1)*((eta + 1)*(eta + 1)) + 51.75 - 33.75*(eta + 1)*(eta + 1));
                      case 55:
                        return RealGradient(0, -36.0*eta - 11.25*xi + 45.0*(eta + 1)*(xi + 1) - 28.125*(xi + 1)*(eta + 1)*(eta + 1) - 38.25 + 22.5*((eta + 1)*(eta + 1)));
                      case 56:
                        return RealGradient(0, 0);
                      case 57:
                        return RealGradient(0, 3.75*xi - 0.75);
                      case 58:
                        return RealGradient(0, -3.75*xi - 0.75);
                      case 59:
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
                        return sign * RealGradient(-5625.0*eta - 3281.25*xi + 4921.875*(2*eta + 2)*(2*xi + 2) - 9843.75*(2*eta + 2)*(xi + 1)*(xi + 1) + 2953.125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 17226.5625*(2*xi + 2)*(eta + 1)*(eta + 1) + 11484.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 7968.75 + 9843.75*((eta + 1)*(eta + 1)) + 34453.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10335.9375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3281.25*((xi + 1)*(xi + 1)) - 22968.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 6562.5*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 984.375*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1476.5625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 1:
                        return sign * RealGradient(1678.7109375*eta + 1030.517578125*xi - 1545.7763671875*(2*eta + 2)*(2*xi + 2) + 2991.5771484375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 853.6376953125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 5410.21728515625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 3606.8115234375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 811.532592773438*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2429.443359375 - 2937.744140625*(eta + 1)*(eta + 1) - 10470.5200195313*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2987.73193359375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 997.1923828125*(xi + 1)*(xi + 1) + 6980.3466796875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1570.57800292969*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1958.49609375*((eta + 1)*(eta + 1)*(eta + 1)) - 1991.8212890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 284.5458984375*((xi + 1)*(xi + 1)*(xi + 1)) + 448.159790039063*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 440.66162109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 2:
                        return sign * RealGradient(-984.375*eta - 902.34375*xi + 1353.515625*(2*eta + 2)*(2*xi + 2) - 3322.265625*(2*eta + 2)*(xi + 1)*(xi + 1) + 1107.421875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 4737.3046875*(2*xi + 2)*(eta + 1)*(eta + 1) + 3158.203125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 710.595703125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1722.65625 + 1722.65625*((eta + 1)*(eta + 1)) + 11627.9296875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3875.9765625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1107.421875*((xi + 1)*(xi + 1)) - 7751.953125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1744.189453125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1148.4375*(eta + 1)*(eta + 1)*(eta + 1) + 2583.984375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 369.140625*(xi + 1)*(xi + 1)*(xi + 1) - 581.396484375*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 258.3984375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 3:
                        return sign * RealGradient(413.0859375*eta + 456.298828125*xi - 684.4482421875*(2*eta + 2)*(2*xi + 2) + 2130.2490234375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 853.6376953125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 2395.56884765625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1597.0458984375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 359.335327148438*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 800.537109375 - 722.900390625*(eta + 1)*(eta + 1) - 7455.87158203125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2987.73193359375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 710.0830078125*(xi + 1)*(xi + 1) + 4970.5810546875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1118.38073730469*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 481.93359375*((eta + 1)*(eta + 1)*(eta + 1)) - 1991.8212890625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 284.5458984375*((xi + 1)*(xi + 1)*(xi + 1)) + 448.159790039063*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 108.43505859375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 4:
                        return sign * RealGradient(-2250.0*eta - 1968.75*xi + 2953.125*(2*eta + 2)*(2*xi + 2) - 7875.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 2953.125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 10335.9375*(2*xi + 2)*(eta + 1)*(eta + 1) + 6890.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3843.75 + 3937.5*((eta + 1)*(eta + 1)) + 27562.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10335.9375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2625.0*((xi + 1)*(xi + 1)) - 18375.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2625.0*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 984.375*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 5:
                        return sign * RealGradient(0, 656.25*eta + 2250.0*xi - 1968.75*(2*eta + 2)*(2*xi + 2) + 10335.9375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 9187.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 2583.984375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 3937.5*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1181.25*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 2718.75 - 656.25*(eta + 1)*(eta + 1) - 20671.875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 18375.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5906.25*(xi + 1)*(xi + 1) + 6201.5625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 196.875*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 5250.0*((xi + 1)*(xi + 1)*(xi + 1)) - 1476.5625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 6:
                        return sign * RealGradient(0, -206.103515625*eta - 671.484375*xi + 618.310546875*(2*eta + 2)*(2*xi + 2) - 3246.13037109375*(2*eta + 2)*(xi + 1)*(xi + 1) + 2885.44921875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 811.532592773438*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1196.630859375*(2*xi + 2)*(eta + 1)*(eta + 1) + 341.455078125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 821.630859375 + 199.4384765625*((eta + 1)*(eta + 1)) + 6282.31201171875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 5584.27734375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1570.57800292969*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1762.646484375*((xi + 1)*(xi + 1)) - 1792.63916015625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 56.9091796875*(eta + 1)*(eta + 1)*(eta + 1) + 1593.45703125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 448.159790039063*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1566.796875*(xi + 1)*(xi + 1)*(xi + 1) + 440.66162109375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 7:
                        return sign * RealGradient(0, 180.46875*eta + 393.75*xi - 541.40625*(2*eta + 2)*(2*xi + 2) + 2842.3828125*(2*eta + 2)*((xi + 1)*(xi + 1)) - 2526.5625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 710.595703125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1328.90625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 442.96875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 541.40625 - 221.484375*(eta + 1)*(eta + 1) - 6976.7578125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1744.189453125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1033.59375*(xi + 1)*(xi + 1) + 2325.5859375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 73.828125*((eta + 1)*(eta + 1)*(eta + 1)) - 2067.1875*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 581.396484375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 918.75*((xi + 1)*(xi + 1)*(xi + 1)) - 258.3984375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 8:
                        return sign * RealGradient(0, -91.259765625*eta - 165.234375*xi + 273.779296875*(2*eta + 2)*(2*xi + 2) - 1437.34130859375*(2*eta + 2)*(xi + 1)*(xi + 1) + 1277.63671875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 359.335327148438*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 852.099609375*(2*xi + 2)*(eta + 1)*(eta + 1) + 341.455078125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 242.724609375 + 142.0166015625*((eta + 1)*(eta + 1)) + 4473.52294921875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 3976.46484375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1118.38073730469*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 433.740234375*((xi + 1)*(xi + 1)) - 1792.63916015625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 56.9091796875*(eta + 1)*(eta + 1)*(eta + 1) + 1593.45703125*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 448.159790039063*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 385.546875*(xi + 1)*(xi + 1)*(xi + 1) + 108.43505859375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 9:
                        return sign * RealGradient(0, 393.75*eta + 900.0*xi - 1181.25*(2*eta + 2)*(2*xi + 2) + 6201.5625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 5512.5*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 3150.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1181.25*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1218.75 - 525.0*(eta + 1)*(eta + 1) - 16537.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 14700.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2362.5*(xi + 1)*(xi + 1) + 6201.5625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 196.875*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2100.0*((xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 10:
                        return sign * RealGradient(-900.0*eta - 393.75*xi + 1181.25*(2*eta + 2)*(2*xi + 2) - 3150.0*(2*eta + 2)*(xi + 1)*(xi + 1) + 1181.25*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 6201.5625*(2*xi + 2)*(eta + 1)*(eta + 1) + 5512.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 1218.75 + 2362.5*((eta + 1)*(eta + 1)) + 16537.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6201.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 525.0*((xi + 1)*(xi + 1)) - 14700.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2100.0*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 11:
                        return sign * RealGradient(165.234375*eta + 91.259765625*xi - 273.779296875*(2*eta + 2)*(2*xi + 2) + 852.099609375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 341.455078125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1437.34130859375*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1277.63671875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 359.335327148438*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 242.724609375 - 433.740234375*(eta + 1)*(eta + 1) - 4473.52294921875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1792.63916015625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 142.0166015625*(xi + 1)*(xi + 1) + 3976.46484375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1118.38073730469*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 385.546875*((eta + 1)*(eta + 1)*(eta + 1)) - 1593.45703125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 56.9091796875*((xi + 1)*(xi + 1)*(xi + 1)) + 448.159790039063*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 108.43505859375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 12:
                        return sign * RealGradient(-393.75*eta - 180.46875*xi + 541.40625*(2*eta + 2)*(2*xi + 2) - 1328.90625*(2*eta + 2)*(xi + 1)*(xi + 1) + 442.96875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 2842.3828125*(2*xi + 2)*(eta + 1)*(eta + 1) + 2526.5625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 710.595703125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 541.40625 + 1033.59375*((eta + 1)*(eta + 1)) + 6976.7578125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2325.5859375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 221.484375*((xi + 1)*(xi + 1)) - 6201.5625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1744.189453125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 918.75*(eta + 1)*(eta + 1)*(eta + 1) + 2067.1875*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 73.828125*(xi + 1)*(xi + 1)*(xi + 1) - 581.396484375*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 258.3984375*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 13:
                        return sign * RealGradient(671.484375*eta + 206.103515625*xi - 618.310546875*(2*eta + 2)*(2*xi + 2) + 1196.630859375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 341.455078125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 3246.13037109375*(2*xi + 2)*((eta + 1)*(eta + 1)) - 2885.44921875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 811.532592773438*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 821.630859375 - 1762.646484375*(eta + 1)*(eta + 1) - 6282.31201171875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 1792.63916015625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 199.4384765625*(xi + 1)*(xi + 1) + 5584.27734375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 1570.57800292969*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1566.796875*((eta + 1)*(eta + 1)*(eta + 1)) - 1593.45703125*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 56.9091796875*((xi + 1)*(xi + 1)*(xi + 1)) + 448.159790039063*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 440.66162109375*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 14:
                        return sign * RealGradient(-2250.0*eta - 656.25*xi + 1968.75*(2*eta + 2)*(2*xi + 2) - 3937.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 1181.25*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 10335.9375*(2*xi + 2)*(eta + 1)*(eta + 1) + 9187.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 2718.75 + 5906.25*((eta + 1)*(eta + 1)) + 20671.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6201.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 656.25*((xi + 1)*(xi + 1)) - 18375.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 5250.0*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 196.875*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1476.5625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 15:
                        return sign * RealGradient(0, 1968.75*eta + 2250.0*xi - 2953.125*(2*eta + 2)*(2*xi + 2) + 10335.9375*(2*eta + 2)*((xi + 1)*(xi + 1)) - 6890.625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 7875.0*(2*xi + 2)*((eta + 1)*(eta + 1)) - 2953.125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 3843.75 - 2625.0*(eta + 1)*(eta + 1) - 27562.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 18375.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3937.5*(xi + 1)*(xi + 1) + 10335.9375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 984.375*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2625.0*((xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 16:
                        return sign * RealGradient(0, -456.298828125*eta - 413.0859375*xi + 684.4482421875*(2*eta + 2)*(2*xi + 2) - 2395.56884765625*(2*eta + 2)*(xi + 1)*(xi + 1) + 1597.0458984375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 359.335327148438*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2130.2490234375*(2*xi + 2)*(eta + 1)*(eta + 1) + 853.6376953125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 800.537109375 + 710.0830078125*((eta + 1)*(eta + 1)) + 7455.87158203125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 4970.5810546875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1118.38073730469*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 722.900390625*((xi + 1)*(xi + 1)) - 2987.73193359375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 284.5458984375*(eta + 1)*(eta + 1)*(eta + 1) + 1991.8212890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 448.159790039063*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 481.93359375*(xi + 1)*(xi + 1)*(xi + 1) + 108.43505859375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 17:
                        return sign * RealGradient(0, 902.34375*eta + 984.375*xi - 1353.515625*(2*eta + 2)*(2*xi + 2) + 4737.3046875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 3158.203125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 710.595703125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 3322.265625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 1107.421875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1722.65625 - 1107.421875*(eta + 1)*(eta + 1) - 11627.9296875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 7751.953125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1744.189453125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1722.65625*(xi + 1)*(xi + 1) + 3875.9765625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 369.140625*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 581.396484375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1148.4375*((xi + 1)*(xi + 1)*(xi + 1)) - 258.3984375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 18:
                        return sign * RealGradient(0, -1030.517578125*eta - 1678.7109375*xi + 1545.7763671875*(2*eta + 2)*(2*xi + 2) - 5410.21728515625*(2*eta + 2)*(xi + 1)*(xi + 1) + 3606.8115234375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 811.532592773438*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2991.5771484375*(2*xi + 2)*(eta + 1)*(eta + 1) + 853.6376953125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 2429.443359375 + 997.1923828125*((eta + 1)*(eta + 1)) + 10470.5200195313*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6980.3466796875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1570.57800292969*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2937.744140625*((xi + 1)*(xi + 1)) - 2987.73193359375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 284.5458984375*(eta + 1)*(eta + 1)*(eta + 1) + 1991.8212890625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 448.159790039063*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1958.49609375*(xi + 1)*(xi + 1)*(xi + 1) + 440.66162109375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 19:
                        return sign * RealGradient(0, 3281.25*eta + 5625.0*xi - 4921.875*(2*eta + 2)*(2*xi + 2) + 17226.5625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 11484.375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 2583.984375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 9843.75*(2*xi + 2)*((eta + 1)*(eta + 1)) - 2953.125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 7968.75 - 3281.25*(eta + 1)*(eta + 1) - 34453.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 22968.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 9843.75*(xi + 1)*(xi + 1) + 10335.9375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 984.375*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 6562.5*((xi + 1)*(xi + 1)*(xi + 1)) - 1476.5625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 20:
                        return RealGradient(0, 2559.375*eta + 5287.5*xi - 4626.5625*(2*eta + 2)*(2*xi + 2) + 17005.078125*(2*eta + 2)*((xi + 1)*(xi + 1)) - 11484.375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 2583.984375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 9253.125*(2*xi + 2)*((eta + 1)*(eta + 1)) - 2775.9375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 7115.625 - 2559.375*(eta + 1)*(eta + 1) - 34010.15625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 22968.75*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 9717.1875*(xi + 1)*(xi + 1) + 10203.046875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 767.8125*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 6562.5*((xi + 1)*(xi + 1)*(xi + 1)) - 1476.5625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 21:
                        return RealGradient(0, -525.0*eta - 2081.25*xi + 1821.09375*(2*eta + 2)*(2*xi + 2) - 10114.453125*(2*eta + 2)*(xi + 1)*(xi + 1) + 9187.5*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 2583.984375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3642.1875*(2*xi + 2)*(eta + 1)*(eta + 1) + 1092.65625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 2456.25 + 525.0*((eta + 1)*(eta + 1)) + 20228.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 18375.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5167.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 5779.6875*((xi + 1)*(xi + 1)) - 6068.671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 157.5*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5250.0*(xi + 1)*(xi + 1)*(xi + 1) + 1476.5625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 22:
                        return RealGradient(0, 1181.25*eta + 843.75*xi - 738.28125*(2*eta + 2)*(2*xi + 2) + 664.453125*(2*eta + 2)*((xi + 1)*(xi + 1)) + 1476.5625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 442.96875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1687.5 - 1181.25*(eta + 1)*(eta + 1) - 1328.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 379.6875*(xi + 1)*(xi + 1) + 398.671875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 354.375*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 23:
                        return RealGradient(0, -590.625*eta - 675.0*xi + 590.625*(2*eta + 2)*(2*xi + 2) - 664.453125*(2*eta + 2)*(xi + 1)*(xi + 1) - 1181.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 354.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1096.875 + 590.625*((eta + 1)*(eta + 1)) + 1328.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 379.6875*((xi + 1)*(xi + 1)) - 398.671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 177.1875*(eta + 1)*(eta + 1)*(eta + 1));
                      case 24:
                        return RealGradient(-5287.5*eta - 2559.375*xi + 4626.5625*(2*eta + 2)*(2*xi + 2) - 9253.125*(2*eta + 2)*(xi + 1)*(xi + 1) + 2775.9375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 17005.078125*(2*xi + 2)*(eta + 1)*(eta + 1) + 11484.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 2583.984375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 7115.625 + 9717.1875*((eta + 1)*(eta + 1)) + 34010.15625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10203.046875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2559.375*((xi + 1)*(xi + 1)) - 22968.75*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5167.96875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 6562.5*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 767.8125*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1476.5625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 25:
                        return RealGradient(2081.25*eta + 525.0*xi - 1821.09375*(2*eta + 2)*(2*xi + 2) + 3642.1875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1092.65625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 10114.453125*(2*xi + 2)*((eta + 1)*(eta + 1)) - 9187.5*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 2583.984375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 2456.25 - 5779.6875*(eta + 1)*(eta + 1) - 20228.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6068.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 525.0*(xi + 1)*(xi + 1) + 18375.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 5167.96875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 5250.0*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 157.5*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 1476.5625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 26:
                        return RealGradient(-843.75*eta - 1181.25*xi + 738.28125*(2*eta + 2)*(2*xi + 2) - 1476.5625*(2*eta + 2)*(xi + 1)*(xi + 1) + 442.96875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 664.453125*(2*xi + 2)*(eta + 1)*(eta + 1) - 1687.5 + 379.6875*((eta + 1)*(eta + 1)) + 1328.90625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 398.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1181.25*((xi + 1)*(xi + 1)) - 354.375*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 27:
                        return RealGradient(675.0*eta + 590.625*xi - 590.625*(2*eta + 2)*(2*xi + 2) + 1181.25*(2*eta + 2)*((xi + 1)*(xi + 1)) - 354.375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 664.453125*(2*xi + 2)*((eta + 1)*(eta + 1)) + 1096.875 - 379.6875*(eta + 1)*(eta + 1) - 1328.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(xi + 1)*(xi + 1) + 177.1875*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 28:
                        return RealGradient(-2115.0*eta - 1535.625*xi + 2775.9375*(2*eta + 2)*(2*xi + 2) - 7402.5*(2*eta + 2)*(xi + 1)*(xi + 1) + 2775.9375*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 10203.046875*(2*xi + 2)*(eta + 1)*(eta + 1) + 6890.625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1550.390625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 3358.125 + 3886.875*((eta + 1)*(eta + 1)) + 27208.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 10203.046875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2047.5*((xi + 1)*(xi + 1)) - 18375.0*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 4134.375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 2625.0*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 767.8125*(xi + 1)*(xi + 1)*(xi + 1) - 1550.390625*(xi + 1)*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 590.625*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 29:
                        return RealGradient(832.5*eta + 315.0*xi - 1092.65625*(2*eta + 2)*(2*xi + 2) + 2913.75*(2*eta + 2)*((xi + 1)*(xi + 1)) - 1092.65625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 6068.671875*(2*xi + 2)*((eta + 1)*(eta + 1)) - 5512.5*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 1550.390625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 1087.5 - 2311.875*(eta + 1)*(eta + 1) - 16183.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6068.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 420.0*(xi + 1)*(xi + 1) + 14700.0*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) - 4134.375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2100.0*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 157.5*((xi + 1)*(xi + 1)*(xi + 1)) + 1550.390625*((xi + 1)*(xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) - 590.625*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 30:
                        return RealGradient(-337.5*eta - 708.75*xi + 442.96875*(2*eta + 2)*(2*xi + 2) - 1181.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 442.96875*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 398.671875*(2*xi + 2)*(eta + 1)*(eta + 1) - 911.25 + 151.875*((eta + 1)*(eta + 1)) + 1063.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 398.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 945.0*((xi + 1)*(xi + 1)) - 354.375*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 31:
                        return RealGradient(270.0*eta + 354.375*xi - 354.375*(2*eta + 2)*(2*xi + 2) + 945.0*(2*eta + 2)*((xi + 1)*(xi + 1)) - 354.375*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 398.671875*(2*xi + 2)*((eta + 1)*(eta + 1)) + 556.875 - 151.875*(eta + 1)*(eta + 1) - 1063.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 472.5*(xi + 1)*(xi + 1) + 177.1875*((xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 32:
                        return RealGradient(0, 1535.625*eta + 2115.0*xi - 2775.9375*(2*eta + 2)*(2*xi + 2) + 10203.046875*(2*eta + 2)*((xi + 1)*(xi + 1)) - 6890.625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 7402.5*(2*xi + 2)*((eta + 1)*(eta + 1)) - 2775.9375*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 3358.125 - 2047.5*(eta + 1)*(eta + 1) - 27208.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 18375.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3886.875*(xi + 1)*(xi + 1) + 10203.046875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 767.8125*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2625.0*((xi + 1)*(xi + 1)*(xi + 1)) - 590.625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 33:
                        return RealGradient(0, -315.0*eta - 832.5*xi + 1092.65625*(2*eta + 2)*(2*xi + 2) - 6068.671875*(2*eta + 2)*(xi + 1)*(xi + 1) + 5512.5*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2913.75*(2*xi + 2)*(eta + 1)*(eta + 1) + 1092.65625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 1087.5 + 420.0*((eta + 1)*(eta + 1)) + 16183.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 14700.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4134.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2311.875*((xi + 1)*(xi + 1)) - 6068.671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 157.5*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2100.0*(xi + 1)*(xi + 1)*(xi + 1) + 590.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 34:
                        return RealGradient(0, 708.75*eta + 337.5*xi - 442.96875*(2*eta + 2)*(2*xi + 2) + 398.671875*(2*eta + 2)*((xi + 1)*(xi + 1)) + 1181.25*(2*xi + 2)*((eta + 1)*(eta + 1)) - 442.96875*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 911.25 - 945.0*(eta + 1)*(eta + 1) - 1063.125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 151.875*(xi + 1)*(xi + 1) + 398.671875*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 354.375*((eta + 1)*(eta + 1)*(eta + 1)));
                      case 35:
                        return RealGradient(0, -354.375*eta - 270.0*xi + 354.375*(2*eta + 2)*(2*xi + 2) - 398.671875*(2*eta + 2)*(xi + 1)*(xi + 1) - 945.0*(2*xi + 2)*(eta + 1)*(eta + 1) + 354.375*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 556.875 + 472.5*((eta + 1)*(eta + 1)) + 1063.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 151.875*((xi + 1)*(xi + 1)) - 398.671875*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 177.1875*(eta + 1)*(eta + 1)*(eta + 1));
                      case 36:
                        return RealGradient(594.0*eta + 60.0*xi - 123.75*(2*eta + 2)*(2*xi + 2) + 478.125*(2*xi + 2)*((eta + 1)*(eta + 1)) - 328.125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 73.828125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 582.0 - 1147.5*(eta + 1)*(eta + 1) + 787.5*((eta + 1)*(eta + 1)*(eta + 1)) - 177.1875*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 37:
                        return RealGradient(396.0*eta + 60.0*xi - 123.75*(2*eta + 2)*(2*xi + 2) + 478.125*(2*xi + 2)*((eta + 1)*(eta + 1)) - 328.125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 73.828125*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 408.0 - 765.0*(eta + 1)*(eta + 1) + 525.0*((eta + 1)*(eta + 1)*(eta + 1)) - 118.125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 38:
                        return RealGradient(-247.5*eta - 30.0*xi + 61.875*(2*eta + 2)*(2*xi + 2) - 239.0625*(2*xi + 2)*(eta + 1)*(eta + 1) + 164.0625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 36.9140625*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 247.5 + 478.125*((eta + 1)*(eta + 1)) - 328.125*(eta + 1)*(eta + 1)*(eta + 1) + 73.828125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 39:
                        return RealGradient(0, 60.0*eta + 594.0*xi - 123.75*(2*eta + 2)*(2*xi + 2) + 478.125*(2*eta + 2)*((xi + 1)*(xi + 1)) - 328.125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 582.0 - 1147.5*(xi + 1)*(xi + 1) + 787.5*((xi + 1)*(xi + 1)*(xi + 1)) - 177.1875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 40:
                        return RealGradient(0, 60.0*eta + 396.0*xi - 123.75*(2*eta + 2)*(2*xi + 2) + 478.125*(2*eta + 2)*((xi + 1)*(xi + 1)) - 328.125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 408.0 - 765.0*(xi + 1)*(xi + 1) + 525.0*((xi + 1)*(xi + 1)*(xi + 1)) - 118.125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 41:
                        return RealGradient(0, -30.0*eta - 247.5*xi + 61.875*(2*eta + 2)*(2*xi + 2) - 239.0625*(2*eta + 2)*(xi + 1)*(xi + 1) + 164.0625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 36.9140625*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 247.5 + 478.125*((xi + 1)*(xi + 1)) - 328.125*(xi + 1)*(xi + 1)*(xi + 1) + 73.828125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 42:
                        return RealGradient(0, -7.5*eta - 216.0*xi + 45.0*(2*eta + 2)*(2*xi + 2) - 281.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 262.5*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 73.828125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 214.5 + 675.0*((xi + 1)*(xi + 1)) - 630.0*(xi + 1)*(xi + 1)*(xi + 1) + 177.1875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 43:
                        return RealGradient(0, -7.5*eta - 144.0*xi + 45.0*(2*eta + 2)*(2*xi + 2) - 281.25*(2*eta + 2)*(xi + 1)*(xi + 1) + 262.5*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)) - 73.828125*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 145.5 + 450.0*((xi + 1)*(xi + 1)) - 420.0*(xi + 1)*(xi + 1)*(xi + 1) + 118.125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 44:
                        return RealGradient(0, 3.75*eta + 90.0*xi - 22.5*(2*eta + 2)*(2*xi + 2) + 140.625*(2*eta + 2)*((xi + 1)*(xi + 1)) - 131.25*(2*eta + 2)*(xi + 1)*(xi + 1)*(xi + 1) + 36.9140625*(2*eta + 2)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 90.0 - 281.25*(xi + 1)*(xi + 1) + 262.5*((xi + 1)*(xi + 1)*(xi + 1)) - 73.828125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 45:
                        return RealGradient(-216.0*eta - 7.5*xi + 45.0*(2*eta + 2)*(2*xi + 2) - 281.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 262.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 73.828125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 214.5 + 675.0*((eta + 1)*(eta + 1)) - 630.0*(eta + 1)*(eta + 1)*(eta + 1) + 177.1875*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 46:
                        return RealGradient(-144.0*eta - 7.5*xi + 45.0*(2*eta + 2)*(2*xi + 2) - 281.25*(2*xi + 2)*(eta + 1)*(eta + 1) + 262.5*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)) - 73.828125*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 145.5 + 450.0*((eta + 1)*(eta + 1)) - 420.0*(eta + 1)*(eta + 1)*(eta + 1) + 118.125*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 47:
                        return RealGradient(90.0*eta + 3.75*xi - 22.5*(2*eta + 2)*(2*xi + 2) + 140.625*(2*xi + 2)*((eta + 1)*(eta + 1)) - 131.25*(2*xi + 2)*(eta + 1)*(eta + 1)*(eta + 1) + 36.9140625*(2*xi + 2)*((eta + 1)*(eta + 1)*(eta + 1)*(eta + 1)) + 90.0 - 281.25*(eta + 1)*(eta + 1) + 262.5*((eta + 1)*(eta + 1)*(eta + 1)) - 73.828125*(eta + 1)*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 48:
                        return RealGradient(0, 30.0*eta + 81.0*xi - 16.875*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 75.0 - 33.75*(xi + 1)*(xi + 1));
                      case 49:
                        return RealGradient(0, -7.5*eta - 54.0*xi + 11.25*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*eta + 2)*(xi + 1)*(xi + 1) - 52.5 + 33.75*((xi + 1)*(xi + 1)));
                      case 50:
                        return RealGradient(-81.0*eta - 30.0*xi + 16.875*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 75.0 + 33.75*((eta + 1)*(eta + 1)), 0);
                      case 51:
                        return RealGradient(54.0*eta + 7.5*xi - 11.25*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*xi + 2)*((eta + 1)*(eta + 1)) + 52.5 - 33.75*(eta + 1)*(eta + 1), 0);
                      case 52:
                        return RealGradient(-54.0*eta - 30.0*xi + 16.875*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*xi + 2)*(eta + 1)*(eta + 1) - 60.0 + 22.5*((eta + 1)*(eta + 1)), 0);
                      case 53:
                        return RealGradient(36.0*eta + 7.5*xi - 11.25*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*xi + 2)*((eta + 1)*(eta + 1)) + 37.5 - 22.5*(eta + 1)*(eta + 1), 0);
                      case 54:
                        return RealGradient(0, 30.0*eta + 54.0*xi - 16.875*(2*eta + 2)*(2*xi + 2) + 14.0625*(2*eta + 2)*((xi + 1)*(xi + 1)) + 60.0 - 22.5*(xi + 1)*(xi + 1));
                      case 55:
                        return RealGradient(0, -7.5*eta - 36.0*xi + 11.25*(2*eta + 2)*(2*xi + 2) - 14.0625*(2*eta + 2)*(xi + 1)*(xi + 1) - 37.5 + 22.5*((xi + 1)*(xi + 1)));
                      case 56:
                        return RealGradient(0, 0);
                      case 57:
                        return RealGradient(0, 0);
                      case 58:
                        return RealGradient(0, 0);
                      case 59:
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
                        return sign * RealGradient(-3281.25*eta - 5625.0*xi + 19687.5*(eta + 1)*(xi + 1) - 34453.125*(eta + 1)*(xi + 1)*(xi + 1) + 22968.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 19687.5*(xi + 1)*(eta + 1)*(eta + 1) + 5906.25*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 7968.75 + 3281.25*((eta + 1)*(eta + 1)) + 34453.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 22968.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5167.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 9843.75*((xi + 1)*(xi + 1)) - 10335.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6562.5*(xi + 1)*(xi + 1)*(xi + 1) + 1476.5625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 1:
                        return sign * RealGradient(497.314453125*eta + 1678.7109375*xi - 5875.48828125*(eta + 1)*(xi + 1) + 10820.4345703125*(eta + 1)*((xi + 1)*(xi + 1)) - 6980.3466796875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1493.86596679688*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 5875.48828125*(xi + 1)*((eta + 1)*(eta + 1)) - 1762.646484375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 2033.935546875 - 497.314453125*(eta + 1)*(eta + 1) - 10820.4345703125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6980.3466796875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1493.86596679688*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3091.552734375*(xi + 1)*(xi + 1) + 3246.13037109375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 149.1943359375*((eta + 1)*(eta + 1)*(eta + 1)) - 2094.10400390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 448.159790039063*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1994.384765625*((xi + 1)*(xi + 1)*(xi + 1)) - 426.81884765625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 2:
                        return sign * RealGradient(-246.09375*eta - 984.375*xi + 3445.3125*(eta + 1)*(xi + 1) - 9474.609375*(eta + 1)*(xi + 1)*(xi + 1) + 7751.953125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 1937.98828125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3445.3125*(xi + 1)*(eta + 1)*(eta + 1) + 1033.59375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 1160.15625 + 246.09375*((eta + 1)*(eta + 1)) + 9474.609375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 7751.953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1937.98828125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2707.03125*((xi + 1)*(xi + 1)) - 2842.3828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 73.828125*(eta + 1)*(eta + 1)*(eta + 1) + 2325.5859375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 581.396484375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2214.84375*(xi + 1)*(xi + 1)*(xi + 1) + 553.7109375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 3:
                        return sign * RealGradient(87.158203125*eta + 413.0859375*xi - 1445.80078125*(eta + 1)*(xi + 1) + 4791.1376953125*(eta + 1)*((xi + 1)*(xi + 1)) - 4970.5810546875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1493.86596679688*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1445.80078125*(xi + 1)*((eta + 1)*(eta + 1)) - 433.740234375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 475.341796875 - 87.158203125*(eta + 1)*(eta + 1) - 4791.1376953125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 4970.5810546875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1493.86596679688*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1368.896484375*(xi + 1)*(xi + 1) + 1437.34130859375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 26.1474609375*((eta + 1)*(eta + 1)*(eta + 1)) - 1491.17431640625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 448.159790039063*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1420.166015625*((xi + 1)*(xi + 1)*(xi + 1)) - 426.81884765625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 4:
                        return sign * RealGradient(-656.25*eta - 2250.0*xi + 7875.0*(eta + 1)*(xi + 1) - 20671.875*(eta + 1)*(xi + 1)*(xi + 1) + 18375.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 7875.0*(xi + 1)*(eta + 1)*(eta + 1) + 2362.5*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 2718.75 + 656.25*((eta + 1)*(eta + 1)) + 20671.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 18375.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5167.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 5906.25*((xi + 1)*(xi + 1)) - 6201.5625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 196.875*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 5250.0*(xi + 1)*(xi + 1)*(xi + 1) + 1476.5625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 5:
                        return sign * RealGradient(0, 656.25*xi - 1312.5*(eta + 1)*(xi + 1) + 7875.0*(eta + 1)*((xi + 1)*(xi + 1)) - 13781.25*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9187.5*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2067.1875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 590.625*(xi + 1)*((eta + 1)*(eta + 1)) + 656.25 - 3543.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 930.234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 3937.5*(xi + 1)*(xi + 1) + 6890.625*((xi + 1)*(xi + 1)*(xi + 1)) - 4593.75*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1033.59375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 6:
                        return sign * RealGradient(0, -206.103515625*xi + 398.876953125*(eta + 1)*(xi + 1) - 2393.26171875*(eta + 1)*(xi + 1)*(xi + 1) + 4188.2080078125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 2792.138671875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 628.231201171875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 170.7275390625*(xi + 1)*(eta + 1)*(eta + 1) - 206.103515625 + 1024.365234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1792.63916015625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1195.0927734375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 268.895874023438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1236.62109375*((xi + 1)*(xi + 1)) - 2164.0869140625*(xi + 1)*(xi + 1)*(xi + 1) + 1442.724609375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 324.613037109375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 7:
                        return sign * RealGradient(0, 180.46875*xi - 442.96875*(eta + 1)*(xi + 1) + 2657.8125*(eta + 1)*((xi + 1)*(xi + 1)) - 4651.171875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3100.78125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 697.67578125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 221.484375*(xi + 1)*((eta + 1)*(eta + 1)) + 180.46875 - 1328.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 2325.5859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 348.837890625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1082.8125*(xi + 1)*(xi + 1) + 1894.921875*((xi + 1)*(xi + 1)*(xi + 1)) - 1263.28125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 284.23828125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 8:
                        return sign * RealGradient(0, -91.259765625*xi + 284.033203125*(eta + 1)*(xi + 1) - 1704.19921875*(eta + 1)*(xi + 1)*(xi + 1) + 2982.3486328125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 1988.232421875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 447.352294921875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 170.7275390625*(xi + 1)*(eta + 1)*(eta + 1) - 91.259765625 + 1024.365234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 1792.63916015625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1195.0927734375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 268.895874023438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 547.55859375*((xi + 1)*(xi + 1)) - 958.2275390625*(xi + 1)*(xi + 1)*(xi + 1) + 638.818359375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 143.734130859375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 9:
                        return sign * RealGradient(0, 393.75*xi - 1050.0*(eta + 1)*(xi + 1) + 6300.0*(eta + 1)*((xi + 1)*(xi + 1)) - 11025.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 7350.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1653.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 590.625*(xi + 1)*((eta + 1)*(eta + 1)) + 393.75 - 3543.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 6201.5625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 930.234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2362.5*(xi + 1)*(xi + 1) + 4134.375*((xi + 1)*(xi + 1)*(xi + 1)) - 2756.25*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 620.15625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 10:
                        return sign * RealGradient(-393.75*eta - 900.0*xi + 4725.0*(eta + 1)*(xi + 1) - 12403.125*(eta + 1)*(xi + 1)*(xi + 1) + 11025.0*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 3100.78125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6300.0*(xi + 1)*(eta + 1)*(eta + 1) + 2362.5*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 1218.75 + 525.0*((eta + 1)*(eta + 1)) + 16537.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 14700.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4134.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2362.5*((xi + 1)*(xi + 1)) - 6201.5625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 196.875*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2100.0*(xi + 1)*(xi + 1)*(xi + 1) + 590.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 11:
                        return sign * RealGradient(52.294921875*eta + 165.234375*xi - 867.48046875*(eta + 1)*(xi + 1) + 2874.6826171875*(eta + 1)*((xi + 1)*(xi + 1)) - 2982.3486328125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 896.319580078125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1156.640625*(xi + 1)*((eta + 1)*(eta + 1)) - 433.740234375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 207.568359375 - 69.7265625*(eta + 1)*(eta + 1) - 3832.91015625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3976.46484375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1195.0927734375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 547.55859375*(xi + 1)*(xi + 1) + 1437.34130859375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 26.1474609375*((eta + 1)*(eta + 1)*(eta + 1)) - 1491.17431640625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 448.159790039063*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 568.06640625*((xi + 1)*(xi + 1)*(xi + 1)) - 170.7275390625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 12:
                        return sign * RealGradient(-147.65625*eta - 393.75*xi + 2067.1875*(eta + 1)*(xi + 1) - 5684.765625*(eta + 1)*(xi + 1)*(xi + 1) + 4651.171875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 1162.79296875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2756.25*(xi + 1)*(eta + 1)*(eta + 1) + 1033.59375*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 513.28125 + 196.875*((eta + 1)*(eta + 1)) + 7579.6875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6201.5625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1082.8125*((xi + 1)*(xi + 1)) - 2842.3828125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 73.828125*(eta + 1)*(eta + 1)*(eta + 1) + 2325.5859375*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 581.396484375*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 885.9375*(xi + 1)*(xi + 1)*(xi + 1) + 221.484375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 13:
                        return sign * RealGradient(298.388671875*eta + 671.484375*xi - 3525.29296875*(eta + 1)*(xi + 1) + 6492.2607421875*(eta + 1)*((xi + 1)*(xi + 1)) - 4188.2080078125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 896.319580078125*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 4700.390625*(xi + 1)*((eta + 1)*(eta + 1)) - 1762.646484375*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 913.037109375 - 397.8515625*(eta + 1)*(eta + 1) - 8656.34765625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 5584.27734375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1195.0927734375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 1236.62109375*(xi + 1)*(xi + 1) + 3246.13037109375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 149.1943359375*((eta + 1)*(eta + 1)*(eta + 1)) - 2094.10400390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 448.159790039063*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 797.75390625*((xi + 1)*(xi + 1)*(xi + 1)) - 170.7275390625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 14:
                        return sign * RealGradient(-1968.75*eta - 2250.0*xi + 11812.5*(eta + 1)*(xi + 1) - 20671.875*(eta + 1)*(xi + 1)*(xi + 1) + 13781.25*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 3100.78125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 15750.0*(xi + 1)*(eta + 1)*(eta + 1) + 5906.25*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 3843.75 + 2625.0*((eta + 1)*(eta + 1)) + 27562.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 18375.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4134.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 3937.5*((xi + 1)*(xi + 1)) - 10335.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2625.0*(xi + 1)*(xi + 1)*(xi + 1) + 590.625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 15:
                        return sign * RealGradient(0, 420.0*eta + 1968.75*xi - 5250.0*(eta + 1)*(xi + 1) + 15750.0*(eta + 1)*((xi + 1)*(xi + 1)) - 18375.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9187.5*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1653.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2953.125*(xi + 1)*((eta + 1)*(eta + 1)) + 2231.25 - 236.25*(eta + 1)*(eta + 1) - 8859.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 930.234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 5906.25*(xi + 1)*(xi + 1) + 6890.625*((xi + 1)*(xi + 1)*(xi + 1)) - 3445.3125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 620.15625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 16:
                        return sign * RealGradient(0, -113.61328125*eta - 456.298828125*xi + 1420.166015625*(eta + 1)*(xi + 1) - 4260.498046875*(eta + 1)*(xi + 1)*(xi + 1) + 4970.5810546875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 2485.29052734375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 447.352294921875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 853.6376953125*(xi + 1)*(eta + 1)*(eta + 1) - 533.408203125 + 68.291015625*((eta + 1)*(eta + 1)) + 2560.9130859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2987.73193359375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1493.86596679688*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 268.895874023438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1368.896484375*((xi + 1)*(xi + 1)) - 1597.0458984375*(xi + 1)*(xi + 1)*(xi + 1) + 798.52294921875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 143.734130859375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 17:
                        return sign * RealGradient(0, 177.1875*eta + 902.34375*xi - 2214.84375*(eta + 1)*(xi + 1) + 6644.53125*(eta + 1)*((xi + 1)*(xi + 1)) - 7751.953125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3875.9765625*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 697.67578125*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1107.421875*(xi + 1)*((eta + 1)*(eta + 1)) + 1007.34375 - 88.59375*(eta + 1)*(eta + 1) - 3322.265625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 3875.9765625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1937.98828125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 348.837890625*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2707.03125*(xi + 1)*(xi + 1) + 3158.203125*((xi + 1)*(xi + 1)*(xi + 1)) - 1579.1015625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 284.23828125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 18:
                        return sign * RealGradient(0, -159.55078125*eta - 1030.517578125*xi + 1994.384765625*(eta + 1)*(xi + 1) - 5983.154296875*(eta + 1)*(xi + 1)*(xi + 1) + 6980.3466796875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 3490.17333984375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 628.231201171875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 853.6376953125*(xi + 1)*(eta + 1)*(eta + 1) - 1107.626953125 + 68.291015625*((eta + 1)*(eta + 1)) + 2560.9130859375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 2987.73193359375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1493.86596679688*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 268.895874023438*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3091.552734375*((xi + 1)*(xi + 1)) - 3606.8115234375*(xi + 1)*(xi + 1)*(xi + 1) + 1803.40576171875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 324.613037109375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 19:
                        return sign * RealGradient(0, 525.0*eta + 3281.25*xi - 6562.5*(eta + 1)*(xi + 1) + 19687.5*(eta + 1)*((xi + 1)*(xi + 1)) - 22968.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 11484.375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2067.1875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2953.125*(xi + 1)*((eta + 1)*(eta + 1)) + 3543.75 - 236.25*(eta + 1)*(eta + 1) - 8859.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10335.9375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 930.234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 9843.75*(xi + 1)*(xi + 1) + 11484.375*((xi + 1)*(xi + 1)*(xi + 1)) - 5742.1875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1033.59375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 20:
                        return RealGradient(0, 2559.375*xi - 5118.75*(eta + 1)*(xi + 1) + 18506.25*(eta + 1)*((xi + 1)*(xi + 1)) - 22673.4375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 11484.375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 2067.1875*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2303.4375*(xi + 1)*((eta + 1)*(eta + 1)) + 2559.375 - 8327.8125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10203.046875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 930.234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 9253.125*(xi + 1)*(xi + 1) + 11336.71875*((xi + 1)*(xi + 1)*(xi + 1)) - 5742.1875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1033.59375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 21:
                        return RealGradient(0, -525.0*xi + 1050.0*(eta + 1)*(xi + 1) - 7284.375*(eta + 1)*(xi + 1)*(xi + 1) + 13485.9375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 9187.5*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2067.1875*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 472.5*(xi + 1)*(eta + 1)*(eta + 1) - 525.0 + 3277.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6068.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4134.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 930.234375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3642.1875*((xi + 1)*(xi + 1)) - 6742.96875*(xi + 1)*(xi + 1)*(xi + 1) + 4593.75*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1033.59375*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 22:
                        return RealGradient(0, 1181.25*xi - 2362.5*(eta + 1)*(xi + 1) + 2953.125*(eta + 1)*((xi + 1)*(xi + 1)) - 885.9375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1063.125*(xi + 1)*((eta + 1)*(eta + 1)) + 1181.25 - 1328.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1476.5625*(xi + 1)*(xi + 1) + 442.96875*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 23:
                        return RealGradient(0, -590.625*xi + 1181.25*(eta + 1)*(xi + 1) - 2362.5*(eta + 1)*(xi + 1)*(xi + 1) + 885.9375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 531.5625*(xi + 1)*(eta + 1)*(eta + 1) - 590.625 + 1063.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 398.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1181.25*((xi + 1)*(xi + 1)) - 442.96875*(xi + 1)*(xi + 1)*(xi + 1));
                      case 24:
                        return RealGradient(-3239.0625*eta - 5287.5*xi + 19434.375*(eta + 1)*(xi + 1) - 34010.15625*(eta + 1)*(xi + 1)*(xi + 1) + 22673.4375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 5101.5234375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 19687.5*(xi + 1)*(eta + 1)*(eta + 1) + 5906.25*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 7645.3125 + 3281.25*((eta + 1)*(eta + 1)) + 34453.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 22968.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5167.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 9253.125*((xi + 1)*(xi + 1)) - 10335.9375*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(eta + 1) + 6890.625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 6168.75*(xi + 1)*(xi + 1)*(xi + 1) + 1387.96875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 25:
                        return RealGradient(1926.5625*eta + 2081.25*xi - 11559.375*(eta + 1)*(xi + 1) + 20228.90625*(eta + 1)*((xi + 1)*(xi + 1)) - 13485.9375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3034.3359375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 15750.0*(xi + 1)*((eta + 1)*(eta + 1)) - 5906.25*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 3660.9375 - 2625.0*(eta + 1)*(eta + 1) - 27562.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 18375.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 3642.1875*(xi + 1)*(xi + 1) + 10335.9375*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 984.375*((eta + 1)*(eta + 1)*(eta + 1)) - 6890.625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 2428.125*((xi + 1)*(xi + 1)*(xi + 1)) - 546.328125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 26:
                        return RealGradient(-126.5625*eta - 843.75*xi + 759.375*(eta + 1)*(xi + 1) - 1328.90625*(eta + 1)*(xi + 1)*(xi + 1) + 885.9375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 199.3359375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 829.6875 + 1476.5625*((xi + 1)*(xi + 1)) - 984.375*(xi + 1)*(xi + 1)*(xi + 1) + 221.484375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 27:
                        return RealGradient(126.5625*eta + 675.0*xi - 759.375*(eta + 1)*(xi + 1) + 1328.90625*(eta + 1)*((xi + 1)*(xi + 1)) - 885.9375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 199.3359375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 689.0625 - 1181.25*(xi + 1)*(xi + 1) + 787.5*((xi + 1)*(xi + 1)*(xi + 1)) - 177.1875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 28:
                        return RealGradient(-647.8125*eta - 2115.0*xi + 7773.75*(eta + 1)*(xi + 1) - 20406.09375*(eta + 1)*(xi + 1)*(xi + 1) + 18138.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 5101.5234375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 7875.0*(xi + 1)*(eta + 1)*(eta + 1) + 2362.5*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 2586.5625 + 656.25*((eta + 1)*(eta + 1)) + 20671.875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 18375.0*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 5167.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 5551.875*((xi + 1)*(xi + 1)) - 6201.5625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 196.875*(eta + 1)*(eta + 1)*(eta + 1) + 5512.5*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 1550.390625*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 4935.0*(xi + 1)*(xi + 1)*(xi + 1) + 1387.96875*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 29:
                        return RealGradient(385.3125*eta + 832.5*xi - 4623.75*(eta + 1)*(xi + 1) + 12137.34375*(eta + 1)*((xi + 1)*(xi + 1)) - 10788.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 3034.3359375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 6300.0*(xi + 1)*((eta + 1)*(eta + 1)) - 2362.5*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1148.4375 - 525.0*(eta + 1)*(eta + 1) - 16537.5*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 14700.0*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 4134.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 2185.3125*(xi + 1)*(xi + 1) + 6201.5625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 196.875*((eta + 1)*(eta + 1)*(eta + 1)) - 5512.5*(eta + 1)*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1550.390625*((eta + 1)*(eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 1942.5*((xi + 1)*(xi + 1)*(xi + 1)) - 546.328125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 30:
                        return RealGradient(-25.3125*eta - 337.5*xi + 303.75*(eta + 1)*(xi + 1) - 797.34375*(eta + 1)*(xi + 1)*(xi + 1) + 708.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 199.3359375*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) - 334.6875 + 885.9375*((xi + 1)*(xi + 1)) - 787.5*(xi + 1)*(xi + 1)*(xi + 1) + 221.484375*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)), 0);
                      case 31:
                        return RealGradient(25.3125*eta + 270.0*xi - 303.75*(eta + 1)*(xi + 1) + 797.34375*(eta + 1)*((xi + 1)*(xi + 1)) - 708.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 199.3359375*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) + 272.8125 - 708.75*(xi + 1)*(xi + 1) + 630.0*((xi + 1)*(xi + 1)*(xi + 1)) - 177.1875*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1), 0);
                      case 32:
                        return RealGradient(0, 1535.625*xi - 4095.0*(eta + 1)*(xi + 1) + 14805.0*(eta + 1)*((xi + 1)*(xi + 1)) - 18138.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 9187.5*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 1653.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2303.4375*(xi + 1)*((eta + 1)*(eta + 1)) + 1535.625 - 8327.8125*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 10203.046875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 5167.96875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 930.234375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 5551.875*(xi + 1)*(xi + 1) + 6802.03125*((xi + 1)*(xi + 1)*(xi + 1)) - 3445.3125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 620.15625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 33:
                        return RealGradient(0, -315.0*xi + 840.0*(eta + 1)*(xi + 1) - 5827.5*(eta + 1)*(xi + 1)*(xi + 1) + 10788.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 7350.0*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1653.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 472.5*(xi + 1)*(eta + 1)*(eta + 1) - 315.0 + 3277.96875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 6068.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 4134.375*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 930.234375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 2185.3125*((xi + 1)*(xi + 1)) - 4045.78125*(xi + 1)*(xi + 1)*(xi + 1) + 2756.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 620.15625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 34:
                        return RealGradient(0, 708.75*xi - 1890.0*(eta + 1)*(xi + 1) + 2362.5*(eta + 1)*((xi + 1)*(xi + 1)) - 708.75*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 1063.125*(xi + 1)*((eta + 1)*(eta + 1)) + 708.75 - 1328.90625*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) + 398.671875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)*(xi + 1)) - 885.9375*(xi + 1)*(xi + 1) + 265.78125*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 35:
                        return RealGradient(0, -354.375*xi + 945.0*(eta + 1)*(xi + 1) - 1890.0*(eta + 1)*(xi + 1)*(xi + 1) + 708.75*(eta + 1)*((xi + 1)*(xi + 1)*(xi + 1)) - 531.5625*(xi + 1)*(eta + 1)*(eta + 1) - 354.375 + 1063.125*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) - 398.671875*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 708.75*((xi + 1)*(xi + 1)) - 265.78125*(xi + 1)*(xi + 1)*(xi + 1));
                      case 36:
                        return RealGradient(1147.5*eta + 594.0*xi - 2295.0*(eta + 1)*(xi + 1) + 956.25*(eta + 1)*((xi + 1)*(xi + 1)) + 2362.5*(xi + 1)*((eta + 1)*(eta + 1)) - 708.75*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 1444.5 - 1181.25*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 247.5*(xi + 1)*(xi + 1) + 295.3125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 354.375*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 37:
                        return RealGradient(382.5*eta + 396.0*xi - 1530.0*(eta + 1)*(xi + 1) + 956.25*(eta + 1)*((xi + 1)*(xi + 1)) + 1575.0*(xi + 1)*((eta + 1)*(eta + 1)) - 472.5*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 679.5 - 393.75*(eta + 1)*(eta + 1) - 984.375*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 247.5*(xi + 1)*(xi + 1) + 295.3125*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 118.125*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 38:
                        return RealGradient(-191.25*eta - 247.5*xi + 956.25*(eta + 1)*(xi + 1) - 478.125*(eta + 1)*(xi + 1)*(xi + 1) - 984.375*(xi + 1)*(eta + 1)*(eta + 1) + 295.3125*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 389.25 + 196.875*((eta + 1)*(eta + 1)) + 492.1875*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 123.75*((xi + 1)*(xi + 1)) - 147.65625*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 59.0625*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 39:
                        return RealGradient(0, 60.0*xi + 60.0 - 247.5*(xi + 1)*(xi + 1) + 318.75*((xi + 1)*(xi + 1)*(xi + 1)) - 164.0625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 29.53125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 40:
                        return RealGradient(0, 60.0*xi + 60.0 - 247.5*(xi + 1)*(xi + 1) + 318.75*((xi + 1)*(xi + 1)*(xi + 1)) - 164.0625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 29.53125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 41:
                        return RealGradient(0, -30.0*xi - 30.0 + 123.75*((xi + 1)*(xi + 1)) - 159.375*(xi + 1)*(xi + 1)*(xi + 1) + 82.03125*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 14.765625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 42:
                        return RealGradient(0, -7.5*xi - 7.5 + 90.0*((xi + 1)*(xi + 1)) - 187.5*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 29.53125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 43:
                        return RealGradient(0, -7.5*xi - 7.5 + 90.0*((xi + 1)*(xi + 1)) - 187.5*(xi + 1)*(xi + 1)*(xi + 1) + 131.25*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)) - 29.53125*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1));
                      case 44:
                        return RealGradient(0, 3.75*xi + 3.75 - 45.0*(xi + 1)*(xi + 1) + 93.75*((xi + 1)*(xi + 1)*(xi + 1)) - 65.625*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1) + 14.765625*((xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)*(xi + 1)));
                      case 45:
                        return RealGradient(-675.0*eta - 216.0*xi + 1350.0*(eta + 1)*(xi + 1) - 562.5*(eta + 1)*(xi + 1)*(xi + 1) - 1890.0*(xi + 1)*(eta + 1)*(eta + 1) + 708.75*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 783.0 + 945.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 90.0*((xi + 1)*(xi + 1)) - 295.3125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 354.375*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 46:
                        return RealGradient(-225.0*eta - 144.0*xi + 900.0*(eta + 1)*(xi + 1) - 562.5*(eta + 1)*(xi + 1)*(xi + 1) - 1260.0*(xi + 1)*(eta + 1)*(eta + 1) + 472.5*(xi + 1)*((eta + 1)*(eta + 1)*(eta + 1)) - 333.0 + 315.0*((eta + 1)*(eta + 1)) + 787.5*((eta + 1)*(eta + 1))*((xi + 1)*(xi + 1)) + 90.0*((xi + 1)*(xi + 1)) - 295.3125*(xi + 1)*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) - 118.125*(eta + 1)*(eta + 1)*(eta + 1), 0);
                      case 47:
                        return RealGradient(112.5*eta + 90.0*xi - 562.5*(eta + 1)*(xi + 1) + 281.25*(eta + 1)*((xi + 1)*(xi + 1)) + 787.5*(xi + 1)*((eta + 1)*(eta + 1)) - 295.3125*(xi + 1)*(eta + 1)*(eta + 1)*(eta + 1) + 184.5 - 157.5*(eta + 1)*(eta + 1) - 393.75*(eta + 1)*(eta + 1)*(xi + 1)*(xi + 1) - 45.0*(xi + 1)*(xi + 1) + 147.65625*((xi + 1)*(xi + 1))*((eta + 1)*(eta + 1)*(eta + 1)) + 59.0625*((eta + 1)*(eta + 1)*(eta + 1)), 0);
                      case 48:
                        return RealGradient(0, 30.0*xi + 30.0 - 33.75*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 49:
                        return RealGradient(0, -7.5*xi - 7.5 + 22.5*((xi + 1)*(xi + 1)) - 9.375*(xi + 1)*(xi + 1)*(xi + 1));
                      case 50:
                        return RealGradient(-33.75*eta - 81.0*xi + 67.5*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 74.25 + 33.75*((xi + 1)*(xi + 1)), 0);
                      case 51:
                        return RealGradient(33.75*eta + 54.0*xi - 67.5*(eta + 1)*(xi + 1) + 28.125*(eta + 1)*((xi + 1)*(xi + 1)) + 60.75 - 22.5*(xi + 1)*(xi + 1), 0);
                      case 52:
                        return RealGradient(-11.25*eta - 54.0*xi + 45.0*(eta + 1)*(xi + 1) - 28.125*(eta + 1)*(xi + 1)*(xi + 1) - 51.75 + 33.75*((xi + 1)*(xi + 1)), 0);
                      case 53:
                        return RealGradient(11.25*eta + 36.0*xi - 45.0*(eta + 1)*(xi + 1) + 28.125*(eta + 1)*((xi + 1)*(xi + 1)) + 38.25 - 22.5*(xi + 1)*(xi + 1), 0);
                      case 54:
                        return RealGradient(0, 30.0*xi + 30.0 - 33.75*(xi + 1)*(xi + 1) + 9.375*((xi + 1)*(xi + 1)*(xi + 1)));
                      case 55:
                        return RealGradient(0, -7.5*xi - 7.5 + 22.5*((xi + 1)*(xi + 1)) - 9.375*(xi + 1)*(xi + 1)*(xi + 1));
                      case 56:
                        return RealGradient(3.75*eta - 0.75, 0);
                      case 57:
                        return RealGradient(0, 0);
                      case 58:
                        return RealGradient(0, 0);
                      case 59:
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
                        return sign * RealGradient(30240.0*eta*xi - 12600.0*eta - 12600.0*eta*xi*xi - 8400.0*xi - 25200.0*xi*eta*eta + 2100.0 + 22680.0*(eta*eta) + 7560.0*(xi*xi) - 12600.0*eta*eta*eta, -45360.0*eta*xi + 8400.0*eta + 50400.0*eta*(xi*xi) + 12600.0*xi + 37800.0*xi*(eta*eta) - 1400.0 - 15120.0*eta*eta - 30240.0*xi*xi + 8400.0*(eta*eta*eta) + 21000.0*(xi*xi*xi));
                      case 1:
                        return sign * RealGradient(-10237.5*eta*xi + 4255.78125*eta + 3642.1875*eta*(xi*xi) + 2552.8125*xi + 9154.6875*xi*(eta*eta) - 659.53125 - 7530.46875*eta*eta - 2185.3125*xi*xi + 3986.71875*(eta*eta*eta), 15769.6875*eta*xi - 2460.9375*eta - 18309.375*eta*xi*xi - 3914.53125*xi - 11960.15625*xi*eta*eta + 428.75 + 3189.375*(eta*eta) + 9115.3125*(xi*xi) - 885.9375*eta*eta*eta - 6070.3125*xi*xi*xi);
                      case 2:
                        return sign * RealGradient(7560.0*eta*xi - 1627.5*eta - 4725.0*eta*xi*xi - 2835.0*xi - 4725.0*xi*eta*eta + 577.5 + 472.5*(eta*eta) + 2835.0*(xi*xi) + 787.5*(eta*eta*eta), -4725.0*eta*xi + 105.0*eta + 9450.0*eta*(xi*xi) + 3727.5*xi - 2362.5*xi*eta*eta - 315.0 + 1890.0*(eta*eta) - 10395.0*xi*xi - 1575.0*eta*eta*eta + 7875.0*(xi*xi*xi));
                      case 3:
                        return sign * RealGradient(-1417.5*eta*xi - 521.71875*eta + 3642.1875*eta*(xi*xi) + 1817.8125*xi - 1870.3125*xi*eta*eta - 292.03125 + 2392.03125*(eta*eta) - 2185.3125*xi*xi - 1525.78125*eta*eta*eta, -4075.3125*eta*xi + 584.0625*eta + 3740.625*eta*(xi*xi) - 2077.03125*xi + 4577.34375*xi*(eta*eta) + 148.75 - 590.625*eta*eta + 6910.3125*(xi*xi) - 360.9375*eta*eta*eta - 6070.3125*xi*xi*xi);
                      case 4:
                        return sign * RealGradient(10080.0*eta*xi - 1680.0*eta - 12600.0*eta*xi*xi - 6720.0*xi + 1260.0 + 7560.0*(xi*xi), 8400.0*xi - 700.0 - 25200.0*xi*xi + 21000.0*(xi*xi*xi));
                      case 5:
                        return sign * RealGradient(840.0*eta*(12*xi - 2 - 15*xi*xi), 8400.0*xi - 700.0 - 25200.0*xi*xi + 21000.0*(xi*xi*xi));
                      case 6:
                        return sign * RealGradient(88.59375*eta*(-120*eta*xi + 54*eta + 60*xi - 16 - 30*eta*eta - 45*xi*xi), -15946.875*eta*xi + 2205.0*eta + 21262.5*eta*(xi*xi) + 4252.5*xi + 7973.4375*xi*(eta*eta) - 446.25 - 1949.0625*eta*eta - 10099.6875*xi*xi + 393.75*(eta*eta*eta) + 6644.53125*(xi*xi*xi));
                      case 7:
                        return sign * RealGradient(105.0*eta*(-60*eta*xi + 45*eta + 18*xi - 8 - 45*eta*eta - 15*xi*xi/2), -13230.0*eta*xi + 2520.0*eta + 12600.0*eta*(xi*xi) + 1680.0*xi + 14175.0*xi*(eta*eta) - 245.0 - 4725.0*eta*eta - 2835.0*xi*xi + 2100.0*(eta*eta*eta) + 1312.5*(xi*xi*xi));
                      case 8:
                        return sign * RealGradient(3.28125*eta*(-360*eta*xi + 594*eta + 84*xi - 80 - 810*eta*eta - 15*xi*xi), -4606.875*eta*xi + 1575.0*eta + 2362.5*eta*(xi*xi) + 367.5*xi + 7973.4375*xi*(eta*eta) - 96.25 - 4784.0625*eta*eta - 334.6875*xi*xi + 3543.75*(eta*eta*eta) + 82.03125*(xi*xi*xi));
                      case 9:
                        return sign * RealGradient(0, 0);
                      case 10:
                        return sign * RealGradient(0, 0);
                      case 11:
                        return sign * RealGradient(3.28125*eta*(330*eta*xi + 180*eta - 54*xi - 11 - 465*eta*eta - 15*xi*xi), 118.125*eta*xi + 577.5*eta - 2165.625*eta*xi*xi - 55.78125*xi + 4577.34375*xi*(eta*eta) - 18.59375 - 2392.03125*eta*eta + 88.59375*(xi*xi) + 623.4375*(eta*eta*eta) + 82.03125*(xi*xi*xi));
                      case 12:
                        return sign * RealGradient(105.0*eta*(45*eta*xi - 18*eta - 3*xi + 5/2 + 15*(eta*eta)/2 - 15*xi*xi/2), 9450.0*eta*xi - 1680.0*eta - 9450.0*eta*xi*xi - 52.5*xi - 2362.5*xi*eta*eta + 87.5 - 472.5*eta*eta - 1102.5*xi*xi + 1575.0*(eta*eta*eta) + 1312.5*(xi*xi*xi));
                      case 13:
                        return sign * RealGradient(88.59375*eta*(30*eta*xi - 36*eta + 30*xi - 1 + 45*(eta*eta) - 45*xi*xi), 9568.125*eta*xi - 3622.5*eta - 5315.625*eta*xi*xi + 3986.71875*xi - 11960.15625*xi*eta*eta - 351.09375 + 7530.46875*(eta*eta) - 9833.90625*xi*xi - 3051.5625*eta*eta*eta + 6644.53125*(xi*xi*xi));
                      case 14:
                        return sign * RealGradient(840.0*eta*(-30*eta*xi + 18*eta + 18*xi - 5 - 15*eta*eta - 15*xi*xi), -60480.0*eta*xi + 16800.0*eta + 50400.0*eta*(xi*xi) + 21000.0*xi + 37800.0*xi*(eta*eta) - 3500.0 - 22680.0*eta*eta - 37800.0*xi*xi + 8400.0*(eta*eta*eta) + 21000.0*(xi*xi*xi));
                      case 15:
                        return RealGradient(10080.0*eta*(20*eta*xi - 21*eta - 14*xi + 7 + 15*(eta*eta) + 5*(xi*xi)), 423360.0*eta*xi - 94080.0*eta - 403200.0*eta*xi*xi - 70560.0*xi - 453600.0*xi*eta*eta + 9800.0 + 211680.0*(eta*eta) + 141120.0*(xi*xi) - 134400.0*eta*eta*eta - 84000.0*xi*xi*xi);
                      case 16:
                        return RealGradient(10080.0*eta*(-40*eta*xi + 21*eta + 28*xi - 7 - 15*eta*eta - 25*xi*xi), -846720.0*eta*xi + 188160.0*eta + 806400.0*eta*(xi*xi) + 352800.0*xi + 453600.0*xi*(eta*eta) - 49000.0 - 211680.0*eta*eta - 705600.0*xi*xi + 67200.0*(eta*eta*eta) + 420000.0*(xi*xi*xi));
                      case 17:
                        return RealGradient(20160.0*eta*(-10*eta*xi + 3*eta + 13*xi - 3 - 10*xi*xi), -241920.0*eta*xi + 26880.0*eta + 403200.0*eta*(xi*xi) + 161280.0*xi - 14560.0 - 443520.0*xi*xi + 336000.0*(xi*xi*xi));
                      case 18:
                        return RealGradient(10080.0*eta*(10*eta*xi - 3*eta - 22*xi + 4 + 25*(xi*xi)), 120960.0*eta*xi - 13440.0*eta - 201600.0*eta*xi*xi - 201600.0*xi + 18200.0 + 554400.0*(xi*xi) - 420000.0*xi*xi*xi);
                      case 19:
                        return RealGradient(0, 6720.0*eta - 280.0 - 30240.0*eta*eta + 33600.0*(eta*eta*eta));
                      case 20:
                        return RealGradient(0, -13440.0*eta + 560.0 + 60480.0*(eta*eta) - 67200.0*eta*eta*eta);
                      case 21:
                        return RealGradient(8960.0*eta*(-65*eta*xi/3 + 14*eta + 35*xi/3 - 4 - 10*eta*eta - 20*xi*xi/3), -340480.0*eta*xi + 56000.0*eta + 388266.666666667*eta*(xi*xi) + 71680.0*xi + 268800.0*xi*(eta*eta) - 8306.66666666667 - 81760.0*eta*eta - 158293.333333333*xi*xi + 31111.1111111111*(eta*eta*eta) + 99555.5555555555*(xi*xi*xi));
                      case 22:
                        return RealGradient(2240.0*eta*(280*eta*xi/3 - 49*eta - 136*xi/3 + 13 + 35*(eta*eta) + 100*(xi*xi)/3), 371840.0*eta*xi - 62720.0*eta - 418133.333333333*eta*xi*xi - 91840.0*xi - 235200.0*xi*eta*eta + 10826.6666666667 + 76160.0*(eta*eta) + 200106.666666667*(xi*xi) - 24888.8888888889*eta*eta*eta - 124444.444444444*xi*xi*xi);
                      case 23:
                        return RealGradient(2240.0*eta*(-100*eta*xi/3 + 34*eta + 31*xi/3 - 6 - 35*eta*eta - 10*xi*xi/3), -183680.0*eta*xi + 44800.0*eta + 149333.333333333*eta*(xi*xi) + 20160.0*xi + 235200.0*xi*(eta*eta) - 3546.66666666667 - 109760.0*eta*eta - 29866.6666666667*xi*xi + 69688.8888888889*(eta*eta*eta) + 12444.4444444444*(xi*xi*xi));
                      case 24:
                        return RealGradient(1120.0*eta*(250*eta*xi/3 - 73*eta - 70*xi/3 + 12 + 80*(eta*eta) + 25*(xi*xi)/3), 232960.0*eta*xi - 58240.0*eta - 186666.666666667*eta*xi*xi - 26880.0*xi - 268800.0*xi*eta*eta + 4946.66666666667 + 125440.0*(eta*eta) + 38453.3333333333*(xi*xi) - 64711.1111111111*eta*eta*eta - 15555.5555555556*xi*xi*xi);
                      case 25:
                        return RealGradient(2240.0*eta*(20*eta*xi/3 - eta - 44*xi/3 + 5 - 5*eta*eta + 20*(xi*xi)/3), 4480.0*eta*xi + 5226.66666666667*eta - 29866.6666666667*eta*xi*xi - 11200.0*xi + 33600.0*xi*(eta*eta) + 715.555555555556 - 25760.0*eta*eta + 32853.3333333333*(xi*xi) + 21155.5555555556*(eta*eta*eta) - 24888.8888888889*xi*xi*xi);
                      case 26:
                        return RealGradient(2240.0*eta*(80*eta*xi/3 - 23*eta + 64*xi/3 - 1 + 25*(eta*eta) - 100*xi*xi/3), 165760.0*eta*xi - 46293.3333333333*eta - 119466.666666667*eta*xi*xi + 64960.0*xi - 168000.0*xi*eta*eta - 5351.11111111111 + 80640.0*(eta*eta) - 173226.666666667*xi*xi - 27377.7777777778*eta*eta*eta + 124444.444444444*(xi*xi*xi));
                      case 27:
                        return RealGradient(1120.0*eta*(-40*eta*xi/3 + 11*eta - 2*xi/3 - 1 - 5*eta*eta + 5*(xi*xi)/3), -24640.0*eta*xi + 1493.33333333333*eta + 29866.6666666667*eta*(xi*xi) + 1120.0*xi + 16800.0*xi*(eta*eta) - 155.555555555556 + 5600.0*(eta*eta) + 746.666666666667*(xi*xi) - 9955.55555555556*eta*eta*eta - 3111.11111111111*xi*xi*xi);
                      case 28:
                        return RealGradient(1120.0*eta*(200*eta*xi/3 - 17*eta - 20*xi/3 + 3 - 5*eta*eta - 25*xi*xi/3), 116480.0*eta*xi - 11946.6666666667*eta - 149333.333333333*eta*xi*xi - 3360.0*xi + 16800.0*xi*(eta*eta) + 964.444444444444 - 30240.0*eta*eta - 8213.33333333333*xi*xi + 27377.7777777778*(eta*eta*eta) + 15555.5555555556*(xi*xi*xi));
                      case 29:
                        return RealGradient(2240.0*eta*(-110*eta*xi/3 + 36*eta + 71*xi/3 - 12 - 25*eta*eta - 20*xi*xi/3), -165760.0*eta*xi + 29866.6666666667*eta + 164266.666666667*eta*(xi*xi) + 24640.0*xi + 168000.0*xi*(eta*eta) - 3297.77777777778 - 51520.0*eta*eta - 46293.3333333333*xi*xi + 19911.1111111111*(eta*eta*eta) + 24888.8888888889*(xi*xi*xi));
                      case 30:
                        return RealGradient(1120.0*eta*(170*eta*xi/3 - 23*eta - 134*xi/3 + 10 + 10*(eta*eta) + 125*(xi*xi)/3), 103040.0*eta*xi - 10453.3333333333*eta - 126933.333333333*eta*xi*xi - 56000.0*xi - 33600.0*xi*eta*eta + 5755.55555555556 - 2240.0*eta*eta + 125066.666666667*(xi*xi) + 4977.77777777778*(eta*eta*eta) - 77777.7777777778*xi*xi*xi);
                      case 31:
                        return RealGradient(1120.0*eta*(220*eta*xi/3 - 27*eta - 178*xi/3 + 17 + 5*(eta*eta) + 85*(xi*xi)/3), 105280.0*eta*xi - 7466.66666666667*eta - 164266.666666667*eta*xi*xi - 32480.0*xi - 16800.0*xi*eta*eta + 3017.77777777778 - 19040.0*eta*eta + 79893.3333333333*(xi*xi) + 24888.8888888889*(eta*eta*eta) - 52888.8888888889*xi*xi*xi);
                      case 32:
                        return RealGradient(1120.0*eta*(-80*eta*xi/3 + 5*eta + 116*xi/3 - 7 + 5*(eta*eta) - 125*xi*xi/3), -22400.0*eta*xi - 2986.66666666667*eta + 59733.3333333333*eta*(xi*xi) + 39200.0*xi - 16800.0*xi*eta*eta - 2955.55555555556 + 12320.0*(eta*eta) - 108266.666666667*xi*xi - 4977.77777777778*eta*eta*eta + 77777.7777777778*(xi*xi*xi));
                      case 33:
                        return RealGradient(1120.0*eta*(20*eta*xi - 17*eta - 8*xi + 3 + 15*(eta*eta) + 5*(xi*xi)), 47040.0*eta*xi - 9706.66666666667*eta - 44800.0*eta*xi*xi - 7840.0*xi - 50400.0*xi*eta*eta + 1057.77777777778 + 20160.0*(eta*eta) + 15680.0*(xi*xi) - 11200.0*eta*eta*eta - 9333.33333333333*xi*xi*xi);
                      case 34:
                        return RealGradient(3360.0*eta*(-10*eta*xi + 6*eta + 2*xi - 1 - 5*eta*eta), -67200.0*eta*xi + 12693.3333333333*eta + 67200.0*eta*(xi*xi) + 3360.0*xi + 50400.0*xi*(eta*eta) - 622.222222222222 - 19040.0*eta*eta - 3360.0*xi*xi + 7466.66666666667*(eta*eta*eta));
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
                        return sign * RealGradient(45360.0*eta*xi - 16800.0*eta - 25200.0*eta*xi*xi - 12600.0*xi - 37800.0*xi*eta*eta + 2800.0 + 30240.0*(eta*eta) + 15120.0*(xi*xi) - 16800.0*eta*eta*eta - 4200.0*xi*xi*xi, -30240.0*eta*xi + 4200.0*eta + 37800.0*eta*(xi*xi) + 8400.0*xi + 25200.0*xi*(eta*eta) - 700.0 - 7560.0*eta*eta - 22680.0*xi*xi + 4200.0*(eta*eta*eta) + 16800.0*(xi*xi*xi));
                      case 1:
                        return sign * RealGradient(-15060.9375*eta*xi + 3622.5*eta + 9154.6875*eta*(xi*xi) + 4255.78125*xi + 11960.15625*xi*(eta*eta) - 759.0625 - 4784.0625*eta*eta - 5118.75*xi*xi + 1771.875*(eta*eta*eta) + 1214.0625*(xi*xi*xi), 6378.75*eta*xi + 88.59375*eta - 11960.15625*eta*xi*xi - 2460.9375*xi - 2657.8125*xi*eta*eta + 95.15625 - 1328.90625*eta*eta + 7884.84375*(xi*xi) + 1328.90625*(eta*eta*eta) - 6103.125*xi*xi*xi);
                      case 2:
                        return sign * RealGradient(945.0*eta*xi + 1680.0*eta - 4725.0*eta*xi*xi - 1627.5*xi + 2362.5*xi*(eta*eta) - 4725.0*eta*eta + 3780.0*(xi*xi) + 3150.0*(eta*eta*eta) - 1575.0*xi*xi*xi, 3780.0*eta*xi - 262.5*eta - 2362.5*eta*xi*xi + 105.0*xi - 4725.0*xi*eta*eta + 17.5 + 157.5*(eta*eta) - 2362.5*xi*xi + 262.5*(eta*eta*eta) + 3150.0*(xi*xi*xi));
                      case 3:
                        return sign * RealGradient(4784.0625*eta*xi - 577.5*eta - 1870.3125*eta*xi*xi - 521.71875*xi - 4577.34375*xi*eta*eta + 115.9375 - 59.0625*eta*eta - 708.75*xi*xi + 721.875*(eta*eta*eta) + 1214.0625*(xi*xi*xi), -1181.25*eta*xi + 36.09375*eta + 4577.34375*eta*(xi*xi) + 584.0625*xi - 1082.8125*xi*eta*eta - 27.34375 + 88.59375*(eta*eta) - 2037.65625*xi*xi + 16.40625*(eta*eta*eta) + 1246.875*(xi*xi*xi));
                      case 4:
                        return sign * RealGradient(-1680.0*xi + 140.0 + 5040.0*(xi*xi) - 4200.0*xi*xi*xi, 0);
                      case 5:
                        return sign * RealGradient(-1680.0*xi + 140.0 + 5040.0*(xi*xi) - 4200.0*xi*xi*xi, 0);
                      case 6:
                        return sign * RealGradient(9568.125*eta*xi - 1575.0*eta - 10631.25*eta*xi*xi - 1417.5*xi - 7973.4375*xi*eta*eta + 183.75 + 2303.4375*(eta*eta) + 2657.8125*(xi*xi) - 787.5*eta*eta*eta - 1328.90625*xi*xi*xi, -3898.125*eta*xi + 262.5*eta + 7973.4375*eta*(xi*xi) + 2205.0*xi + 1181.25*xi*(eta*eta) - 113.75 - 137.8125*eta*eta - 7973.4375*xi*xi + 16.40625*(eta*eta*eta) + 7087.5*(xi*xi*xi));
                      case 7:
                        return sign * RealGradient(9450.0*eta*xi - 2520.0*eta - 6300.0*eta*xi*xi - 840.0*xi - 14175.0*xi*eta*eta + 175.0 + 6615.0*(eta*eta) + 945.0*(xi*xi) - 4200.0*eta*eta*eta - 262.5*xi*xi*xi, -9450.0*eta*xi + 840.0*eta + 14175.0*eta*(xi*xi) + 2520.0*xi + 6300.0*xi*(eta*eta) - 175.0 - 945.0*eta*eta - 6615.0*xi*xi + 262.5*(eta*eta*eta) + 4200.0*(xi*xi*xi));
                      case 8:
                        return sign * RealGradient(3898.125*eta*xi - 2205.0*eta - 1181.25*eta*xi*xi - 262.5*xi - 7973.4375*xi*eta*eta + 113.75 + 7973.4375*(eta*eta) + 137.8125*(xi*xi) - 7087.5*eta*eta*eta - 16.40625*xi*xi*xi, -9568.125*eta*xi + 1417.5*eta + 7973.4375*eta*(xi*xi) + 1575.0*xi + 10631.25*xi*(eta*eta) - 183.75 - 2657.8125*eta*eta - 2303.4375*xi*xi + 1328.90625*(eta*eta*eta) + 787.5*(xi*xi*xi));
                      case 9:
                        return sign * RealGradient(0, 1680.0*eta - 140.0 - 5040.0*eta*eta + 4200.0*(eta*eta*eta));
                      case 10:
                        return sign * RealGradient(0, 1680.0*eta - 140.0 - 5040.0*eta*eta + 4200.0*(eta*eta*eta));
                      case 11:
                        return sign * RealGradient(1181.25*eta*xi - 584.0625*eta + 1082.8125*eta*(xi*xi) - 36.09375*xi - 4577.34375*xi*eta*eta + 27.34375 + 2037.65625*(eta*eta) - 88.59375*xi*xi - 1246.875*eta*eta*eta - 16.40625*xi*xi*xi, -4784.0625*eta*xi + 521.71875*eta + 4577.34375*eta*(xi*xi) + 577.5*xi + 1870.3125*xi*(eta*eta) - 115.9375 + 708.75*(eta*eta) + 59.0625*(xi*xi) - 1214.0625*eta*eta*eta - 721.875*xi*xi*xi);
                      case 12:
                        return sign * RealGradient(-3780.0*eta*xi - 105.0*eta + 4725.0*eta*(xi*xi) + 262.5*xi + 2362.5*xi*(eta*eta) - 17.5 + 2362.5*(eta*eta) - 157.5*xi*xi - 3150.0*eta*eta*eta - 262.5*xi*xi*xi, -945.0*eta*xi + 1627.5*eta - 2362.5*eta*xi*xi - 1680.0*xi + 4725.0*xi*(eta*eta) - 3780.0*eta*eta + 4725.0*(xi*xi) + 1575.0*(eta*eta*eta) - 3150.0*xi*xi*xi);
                      case 13:
                        return sign * RealGradient(-6378.75*eta*xi + 2460.9375*eta + 2657.8125*eta*(xi*xi) - 88.59375*xi + 11960.15625*xi*(eta*eta) - 95.15625 - 7884.84375*eta*eta + 1328.90625*(xi*xi) + 6103.125*(eta*eta*eta) - 1328.90625*xi*xi*xi, 15060.9375*eta*xi - 4255.78125*eta - 11960.15625*eta*xi*xi - 3622.5*xi - 9154.6875*xi*eta*eta + 759.0625 + 5118.75*(eta*eta) + 4784.0625*(xi*xi) - 1214.0625*eta*eta*eta - 1771.875*xi*xi*xi);
                      case 14:
                        return sign * RealGradient(30240.0*eta*xi - 8400.0*eta - 25200.0*eta*xi*xi - 4200.0*xi - 37800.0*xi*eta*eta + 700.0 + 22680.0*(eta*eta) + 7560.0*(xi*xi) - 16800.0*eta*eta*eta - 4200.0*xi*xi*xi, -45360.0*eta*xi + 12600.0*eta + 37800.0*eta*(xi*xi) + 16800.0*xi + 25200.0*xi*(eta*eta) - 2800.0 - 15120.0*eta*eta - 30240.0*xi*xi + 4200.0*(eta*eta*eta) + 16800.0*(xi*xi*xi));
                      case 15:
                        return RealGradient(-423360.0*eta*xi + 188160.0*eta + 201600.0*eta*(xi*xi) + 70560.0*xi + 453600.0*xi*(eta*eta) - 19600.0 - 423360.0*eta*eta - 70560.0*xi*xi + 268800.0*(eta*eta*eta) + 16800.0*(xi*xi*xi), 423360.0*eta*xi - 70560.0*eta - 453600.0*eta*xi*xi - 94080.0*xi - 403200.0*xi*eta*eta + 9800.0 + 141120.0*(eta*eta) + 211680.0*(xi*xi) - 84000.0*eta*eta*eta - 134400.0*xi*xi*xi);
                      case 16:
                        return RealGradient(423360.0*eta*xi - 94080.0*eta - 403200.0*eta*xi*xi - 70560.0*xi - 453600.0*xi*eta*eta + 9800.0 + 211680.0*(eta*eta) + 141120.0*(xi*xi) - 134400.0*eta*eta*eta - 84000.0*xi*xi*xi, -423360.0*eta*xi + 70560.0*eta + 453600.0*eta*(xi*xi) + 188160.0*xi + 201600.0*xi*(eta*eta) - 19600.0 - 70560.0*eta*eta - 423360.0*xi*xi + 16800.0*(eta*eta*eta) + 268800.0*(xi*xi*xi));
                      case 17:
                        return RealGradient(120960.0*eta*xi - 13440.0*eta - 201600.0*eta*xi*xi - 60480.0*xi + 6440.0 + 131040.0*(xi*xi) - 67200.0*xi*xi*xi, 26880.0*xi - 1120.0 - 120960.0*xi*xi + 134400.0*(xi*xi*xi));
                      case 18:
                        return RealGradient(-60480.0*eta*xi + 6720.0*eta + 100800.0*eta*(xi*xi) + 40320.0*xi - 3640.0 - 110880.0*xi*xi + 84000.0*(xi*xi*xi), -13440.0*xi + 560.0 + 60480.0*(xi*xi) - 67200.0*xi*xi*xi);
                      case 19:
                        return RealGradient(-13440.0*eta + 560.0 + 60480.0*(eta*eta) - 67200.0*eta*eta*eta, -60480.0*eta*xi + 40320.0*eta + 6720.0*xi + 100800.0*xi*(eta*eta) - 3640.0 - 110880.0*eta*eta + 84000.0*(eta*eta*eta));
                      case 20:
                        return RealGradient(26880.0*eta - 1120.0 - 120960.0*eta*eta + 134400.0*(eta*eta*eta), 120960.0*eta*xi - 60480.0*eta - 13440.0*xi - 201600.0*xi*eta*eta + 6440.0 + 131040.0*(eta*eta) - 67200.0*eta*eta*eta);
                      case 21:
                        return RealGradient(250880.0*eta*xi - 58240.0*eta - 194133.333333333*eta*xi*xi - 35840.0*xi - 268800.0*xi*eta*eta + 5973.33333333333 + 116480.0*(eta*eta) + 52266.6666666667*(xi*xi) - 62222.2222222222*eta*eta*eta - 19911.1111111111*xi*xi*xi, -163520.0*eta*xi + 13440.0*eta + 268800.0*eta*(xi*xi) + 56000.0*xi + 93333.3333333333*xi*(eta*eta) - 3453.33333333333 - 13066.6666666667*eta*eta - 170240.0*xi*xi + 3111.11111111111*(eta*eta*eta) + 129422.222222222*(xi*xi*xi));
                      case 22:
                        return RealGradient(-219520.0*eta*xi + 44800.0*eta + 209066.666666667*eta*(xi*xi) + 29120.0*xi + 235200.0*xi*(eta*eta) - 4293.33333333333 - 91840.0*eta*eta - 50773.3333333333*xi*xi + 49777.7777777778*(eta*eta*eta) + 24888.8888888889*(xi*xi*xi), 152320.0*eta*xi - 13440.0*eta - 235200.0*eta*xi*xi - 62720.0*xi - 74666.6666666667*xi*eta*eta + 4013.33333333333 + 11573.3333333333*(eta*eta) + 185920.0*(xi*xi) - 2488.88888888889*eta*eta*eta - 139377.777777778*xi*xi*xi);
                      case 23:
                        return RealGradient(152320.0*eta*xi - 62720.0*eta - 74666.6666666667*eta*xi*xi - 13440.0*xi - 235200.0*xi*eta*eta + 4013.33333333333 + 185920.0*(eta*eta) + 11573.3333333333*(xi*xi) - 139377.777777778*eta*eta*eta - 2488.88888888889*xi*xi*xi, -219520.0*eta*xi + 29120.0*eta + 235200.0*eta*(xi*xi) + 44800.0*xi + 209066.666666667*xi*(eta*eta) - 4293.33333333333 - 50773.3333333333*eta*eta - 91840.0*xi*xi + 24888.8888888889*(eta*eta*eta) + 49777.7777777778*(xi*xi*xi));
                      case 24:
                        return RealGradient(-163520.0*eta*xi + 56000.0*eta + 93333.3333333333*eta*(xi*xi) + 13440.0*xi + 268800.0*xi*(eta*eta) - 3453.33333333333 - 170240.0*eta*eta - 13066.6666666667*xi*xi + 129422.222222222*(eta*eta*eta) + 3111.11111111111*(xi*xi*xi), 250880.0*eta*xi - 35840.0*eta - 268800.0*eta*xi*xi - 58240.0*xi - 194133.333333333*xi*eta*eta + 5973.33333333333 + 52266.6666666667*(eta*eta) + 116480.0*(xi*xi) - 19911.1111111111*eta*eta*eta - 62222.2222222222*xi*xi*xi);
                      case 25:
                        return RealGradient(-4480.0*eta*xi - 10453.3333333333*eta + 14933.3333333333*eta*(xi*xi) + 11200.0*xi - 33600.0*xi*eta*eta - 1431.11111111111 + 51520.0*(eta*eta) - 16426.6666666667*xi*xi - 42311.1111111111*eta*eta*eta + 4977.77777777778*(xi*xi*xi), -51520.0*eta*xi + 11200.0*eta + 33600.0*eta*(xi*xi) + 5226.66666666667*xi + 63466.6666666667*xi*(eta*eta) - 1151.11111111111 - 25013.3333333333*eta*eta + 2240.0*(xi*xi) + 15555.5555555556*(eta*eta*eta) - 9955.55555555555*xi*xi*xi);
                      case 26:
                        return RealGradient(-103040.0*eta*xi + 29866.6666666667*eta + 59733.3333333333*eta*(xi*xi) - 2240.0*xi + 168000.0*xi*(eta*eta) - 1057.77777777778 - 82880.0*eta*eta + 23893.3333333333*(xi*xi) + 54755.5555555556*(eta*eta*eta) - 24888.8888888889*xi*xi*xi, 161280.0*eta*xi - 26880.0*eta - 168000.0*eta*xi*xi - 46293.3333333333*xi - 82133.3333333333*xi*eta*eta + 5755.55555555556 + 26506.6666666667*(eta*eta) + 82880.0*(xi*xi) - 4977.77777777778*eta*eta*eta - 39822.2222222222*xi*xi*xi);
                      case 27:
                        return RealGradient(24640.0*eta*xi - 2986.66666666667*eta - 14933.3333333333*eta*xi*xi - 1120.0*xi - 16800.0*xi*eta*eta + 311.111111111111 - 11200.0*eta*eta - 373.333333333333*xi*xi + 19911.1111111111*(eta*eta*eta) + 622.222222222222*(xi*xi*xi), 11200.0*eta*xi - 7840.0*eta + 16800.0*eta*(xi*xi) + 1493.33333333333*xi - 29866.6666666667*xi*eta*eta + 591.111111111111 + 21653.3333333333*(eta*eta) - 12320.0*xi*xi - 15555.5555555556*eta*eta*eta + 9955.55555555555*(xi*xi*xi));
                      case 28:
                        return RealGradient(-38080.0*eta*xi - 7466.66666666667*eta + 74666.6666666667*eta*(xi*xi) + 3360.0*xi - 16800.0*xi*eta*eta + 31.1111111111111 + 52640.0*(eta*eta) - 3733.33333333333*xi*xi - 54755.5555555556*eta*eta*eta - 3111.11111111111*xi*xi*xi, -60480.0*eta*xi + 19040.0*eta + 16800.0*eta*(xi*xi) - 11946.6666666667*xi + 82133.3333333333*xi*(eta*eta) - 995.555555555556 - 33226.6666666667*eta*eta + 58240.0*(xi*xi) + 10577.7777777778*(eta*eta*eta) - 49777.7777777778*xi*xi*xi);
                      case 29:
                        return RealGradient(161280.0*eta*xi - 46293.3333333333*eta - 82133.3333333333*eta*xi*xi - 26880.0*xi - 168000.0*xi*eta*eta + 5755.55555555556 + 82880.0*(eta*eta) + 26506.6666666667*(xi*xi) - 39822.2222222222*eta*eta*eta - 4977.77777777778*xi*xi*xi, -103040.0*eta*xi - 2240.0*eta + 168000.0*eta*(xi*xi) + 29866.6666666667*xi + 59733.3333333333*xi*(eta*eta) - 1057.77777777778 + 23893.3333333333*(eta*eta) - 82880.0*xi*xi - 24888.8888888889*eta*eta*eta + 54755.5555555556*(xi*xi*xi));
                      case 30:
                        return RealGradient(-51520.0*eta*xi + 5226.66666666667*eta + 63466.6666666667*eta*(xi*xi) + 11200.0*xi + 33600.0*xi*(eta*eta) - 1151.11111111111 + 2240.0*(eta*eta) - 25013.3333333333*xi*xi - 9955.55555555555*eta*eta*eta + 15555.5555555556*(xi*xi*xi), -4480.0*eta*xi + 11200.0*eta - 33600.0*eta*xi*xi - 10453.3333333333*xi + 14933.3333333333*xi*(eta*eta) - 1431.11111111111 - 16426.6666666667*eta*eta + 51520.0*(xi*xi) + 4977.77777777778*(eta*eta*eta) - 42311.1111111111*xi*xi*xi);
                      case 31:
                        return RealGradient(-60480.0*eta*xi - 11946.6666666667*eta + 82133.3333333333*eta*(xi*xi) + 19040.0*xi + 16800.0*xi*(eta*eta) - 995.555555555556 + 58240.0*(eta*eta) - 33226.6666666667*xi*xi - 49777.7777777778*eta*eta*eta + 10577.7777777778*(xi*xi*xi), -38080.0*eta*xi + 3360.0*eta - 16800.0*eta*xi*xi - 7466.66666666667*xi + 74666.6666666667*xi*(eta*eta) + 31.1111111111111 - 3733.33333333333*eta*eta + 52640.0*(xi*xi) - 3111.11111111111*eta*eta*eta - 54755.5555555556*xi*xi*xi);
                      case 32:
                        return RealGradient(11200.0*eta*xi + 1493.33333333333*eta - 29866.6666666667*eta*xi*xi - 7840.0*xi + 16800.0*xi*(eta*eta) + 591.111111111111 - 12320.0*eta*eta + 21653.3333333333*(xi*xi) + 9955.55555555555*(eta*eta*eta) - 15555.5555555556*xi*xi*xi, 24640.0*eta*xi - 1120.0*eta - 16800.0*eta*xi*xi - 2986.66666666667*xi - 14933.3333333333*xi*eta*eta + 311.111111111111 - 373.333333333333*eta*eta - 11200.0*xi*xi + 622.222222222222*(eta*eta*eta) + 19911.1111111111*(xi*xi*xi));
                      case 33:
                        return RealGradient(-38080.0*eta*xi + 12693.3333333333*eta + 22400.0*eta*(xi*xi) + 3360.0*xi + 50400.0*xi*(eta*eta) - 715.555555555556 - 33600.0*eta*eta - 4480.0*xi*xi + 22400.0*(eta*eta*eta) + 1866.66666666667*(xi*xi*xi), 40320.0*eta*xi - 3360.0*eta - 50400.0*eta*xi*xi - 9706.66666666667*xi - 33600.0*xi*eta*eta + 684.444444444444 + 3360.0*(eta*eta) + 23520.0*(xi*xi) - 14933.3333333333*xi*xi*xi);
                      case 34:
                        return RealGradient(40320.0*eta*xi - 9706.66666666667*eta - 33600.0*eta*xi*xi - 3360.0*xi - 50400.0*xi*eta*eta + 684.444444444444 + 23520.0*(eta*eta) + 3360.0*(xi*xi) - 14933.3333333333*eta*eta*eta, -38080.0*eta*xi + 3360.0*eta + 50400.0*eta*(xi*xi) + 12693.3333333333*xi + 22400.0*xi*(eta*eta) - 715.555555555556 - 4480.0*eta*eta - 33600.0*xi*xi + 1866.66666666667*(eta*eta*eta) + 22400.0*(xi*xi*xi));
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
                        return sign * RealGradient(60480.0*eta*xi - 21000.0*eta - 37800.0*eta*xi*xi - 16800.0*xi - 50400.0*xi*eta*eta + 3500.0 + 37800.0*(eta*eta) + 22680.0*(xi*xi) - 21000.0*eta*eta*eta - 8400.0*xi*xi*xi, 840.0*xi*(30*eta*xi - 18*eta - 18*xi + 5 + 15*(eta*eta) + 15*(xi*xi)));
                      case 1:
                        return sign * RealGradient(-9568.125*eta*xi - 3986.71875*eta + 11960.15625*eta*(xi*xi) + 3622.5*xi + 5315.625*xi*(eta*eta) + 351.09375 + 9833.90625*(eta*eta) - 7530.46875*xi*xi - 6644.53125*eta*eta*eta + 3051.5625*(xi*xi*xi), 88.59375*xi*(-30*eta*xi - 30*eta + 36*xi + 1 + 45*(eta*eta) - 45*xi*xi));
                      case 2:
                        return sign * RealGradient(-9450.0*eta*xi + 52.5*eta + 2362.5*eta*(xi*xi) + 1680.0*xi + 9450.0*xi*(eta*eta) - 87.5 + 1102.5*(eta*eta) + 472.5*(xi*xi) - 1312.5*eta*eta*eta - 1575.0*xi*xi*xi, 105.0*xi*(-45*eta*xi + 3*eta + 18*xi - 5/2 + 15*(eta*eta)/2 - 15*xi*xi/2));
                      case 3:
                        return sign * RealGradient(-118.125*eta*xi + 55.78125*eta - 4577.34375*eta*xi*xi - 577.5*xi + 2165.625*xi*(eta*eta) + 18.59375 - 88.59375*eta*eta + 2392.03125*(xi*xi) - 82.03125*eta*eta*eta - 623.4375*xi*xi*xi, 3.28125*xi*(-330*eta*xi + 54*eta - 180*xi + 11 + 15*(eta*eta) + 465*(xi*xi)));
                      case 4:
                        return sign * RealGradient(0, 0);
                      case 5:
                        return sign * RealGradient(0, 0);
                      case 6:
                        return sign * RealGradient(4606.875*eta*xi - 367.5*eta - 7973.4375*eta*xi*xi - 1575.0*xi - 2362.5*xi*eta*eta + 96.25 + 334.6875*(eta*eta) + 4784.0625*(xi*xi) - 82.03125*eta*eta*eta - 3543.75*xi*xi*xi, 3.28125*xi*(360*eta*xi - 84*eta - 594*xi + 80 + 15*(eta*eta) + 810*(xi*xi)));
                      case 7:
                        return sign * RealGradient(13230.0*eta*xi - 1680.0*eta - 14175.0*eta*xi*xi - 2520.0*xi - 12600.0*xi*eta*eta + 245.0 + 2835.0*(eta*eta) + 4725.0*(xi*xi) - 1312.5*eta*eta*eta - 2100.0*xi*xi*xi, 105.0*xi*(60*eta*xi - 18*eta - 45*xi + 8 + 15*(eta*eta)/2 + 45*(xi*xi)));
                      case 8:
                        return sign * RealGradient(15946.875*eta*xi - 4252.5*eta - 7973.4375*eta*xi*xi - 2205.0*xi - 21262.5*xi*eta*eta + 446.25 + 10099.6875*(eta*eta) + 1949.0625*(xi*xi) - 6644.53125*eta*eta*eta - 393.75*xi*xi*xi, 88.59375*xi*(120*eta*xi - 60*eta - 54*xi + 16 + 45*(eta*eta) + 30*(xi*xi)));
                      case 9:
                        return sign * RealGradient(-8400.0*eta + 700.0 + 25200.0*(eta*eta) - 21000.0*eta*eta*eta, 840.0*xi*(-12*eta + 2 + 15*(eta*eta)));
                      case 10:
                        return sign * RealGradient(-8400.0*eta + 700.0 + 25200.0*(eta*eta) - 21000.0*eta*eta*eta, -10080.0*eta*xi + 6720.0*eta + 1680.0*xi + 12600.0*xi*(eta*eta) - 1260.0 - 7560.0*eta*eta);
                      case 11:
                        return sign * RealGradient(4075.3125*eta*xi + 2077.03125*eta - 4577.34375*eta*xi*xi - 584.0625*xi - 3740.625*xi*eta*eta - 148.75 - 6910.3125*eta*eta + 590.625*(xi*xi) + 6070.3125*(eta*eta*eta) + 360.9375*(xi*xi*xi), 1417.5*eta*xi - 1817.8125*eta + 1870.3125*eta*(xi*xi) + 521.71875*xi - 3642.1875*xi*eta*eta + 292.03125 + 2185.3125*(eta*eta) - 2392.03125*xi*xi + 1525.78125*(xi*xi*xi));
                      case 12:
                        return sign * RealGradient(4725.0*eta*xi - 3727.5*eta + 2362.5*eta*(xi*xi) - 105.0*xi - 9450.0*xi*eta*eta + 315.0 + 10395.0*(eta*eta) - 1890.0*xi*xi - 7875.0*eta*eta*eta + 1575.0*(xi*xi*xi), -7560.0*eta*xi + 2835.0*eta + 4725.0*eta*(xi*xi) + 1627.5*xi + 4725.0*xi*(eta*eta) - 577.5 - 2835.0*eta*eta - 472.5*xi*xi - 787.5*xi*xi*xi);
                      case 13:
                        return sign * RealGradient(-15769.6875*eta*xi + 3914.53125*eta + 11960.15625*eta*(xi*xi) + 2460.9375*xi + 18309.375*xi*(eta*eta) - 428.75 - 9115.3125*eta*eta - 3189.375*xi*xi + 6070.3125*(eta*eta*eta) + 885.9375*(xi*xi*xi), 10237.5*eta*xi - 2552.8125*eta - 9154.6875*eta*xi*xi - 4255.78125*xi - 3642.1875*xi*eta*eta + 659.53125 + 2185.3125*(eta*eta) + 7530.46875*(xi*xi) - 3986.71875*xi*xi*xi);
                      case 14:
                        return sign * RealGradient(45360.0*eta*xi - 12600.0*eta - 37800.0*eta*xi*xi - 8400.0*xi - 50400.0*xi*eta*eta + 1400.0 + 30240.0*(eta*eta) + 15120.0*(xi*xi) - 21000.0*eta*eta*eta - 8400.0*xi*xi*xi, -30240.0*eta*xi + 8400.0*eta + 25200.0*eta*(xi*xi) + 12600.0*xi + 12600.0*xi*(eta*eta) - 2100.0 - 7560.0*eta*eta - 22680.0*xi*xi + 12600.0*(xi*xi*xi));
                      case 15:
                        return RealGradient(-846720.0*eta*xi + 352800.0*eta + 453600.0*eta*(xi*xi) + 188160.0*xi + 806400.0*xi*(eta*eta) - 49000.0 - 705600.0*eta*eta - 211680.0*xi*xi + 420000.0*(eta*eta*eta) + 67200.0*(xi*xi*xi), 10080.0*xi*(-40*eta*xi + 28*eta + 21*xi - 7 - 25*eta*eta - 15*xi*xi));
                      case 16:
                        return RealGradient(423360.0*eta*xi - 70560.0*eta - 453600.0*eta*xi*xi - 94080.0*xi - 403200.0*xi*eta*eta + 9800.0 + 141120.0*(eta*eta) + 211680.0*(xi*xi) - 84000.0*eta*eta*eta - 134400.0*xi*xi*xi, 10080.0*xi*(20*eta*xi - 14*eta - 21*xi + 7 + 5*(eta*eta) + 15*(xi*xi)));
                      case 17:
                        return RealGradient(-13440.0*xi + 560.0 + 60480.0*(xi*xi) - 67200.0*xi*xi*xi, 0);
                      case 18:
                        return RealGradient(6720.0*xi - 280.0 - 30240.0*xi*xi + 33600.0*(xi*xi*xi), 0);
                      case 19:
                        return RealGradient(120960.0*eta*xi - 201600.0*eta - 13440.0*xi - 201600.0*xi*eta*eta + 18200.0 + 554400.0*(eta*eta) - 420000.0*eta*eta*eta, 10080.0*xi*(10*eta*xi - 22*eta - 3*xi + 4 + 25*(eta*eta)));
                      case 20:
                        return RealGradient(-241920.0*eta*xi + 161280.0*eta + 26880.0*xi + 403200.0*xi*(eta*eta) - 14560.0 - 443520.0*eta*eta + 336000.0*(eta*eta*eta), 20160.0*xi*(-10*eta*xi + 13*eta + 3*xi - 3 - 10*eta*eta));
                      case 21:
                        return RealGradient(232960.0*eta*xi - 26880.0*eta - 268800.0*eta*xi*xi - 58240.0*xi - 186666.666666667*xi*eta*eta + 4946.66666666667 + 38453.3333333333*(eta*eta) + 125440.0*(xi*xi) - 15555.5555555556*eta*eta*eta - 64711.1111111111*xi*xi*xi, 1120.0*xi*(250*eta*xi/3 - 70*eta/3 - 73*xi + 12 + 25*(eta*eta)/3 + 80*(xi*xi)));
                      case 22:
                        return RealGradient(-183680.0*eta*xi + 20160.0*eta + 235200.0*eta*(xi*xi) + 44800.0*xi + 149333.333333333*xi*(eta*eta) - 3546.66666666667 - 29866.6666666667*eta*eta - 109760.0*xi*xi + 12444.4444444444*(eta*eta*eta) + 69688.8888888889*(xi*xi*xi), 2240.0*xi*(-100*eta*xi/3 + 31*eta/3 + 34*xi - 6 - 10*eta*eta/3 - 35*xi*xi));
                      case 23:
                        return RealGradient(371840.0*eta*xi - 91840.0*eta - 235200.0*eta*xi*xi - 62720.0*xi - 418133.333333333*xi*eta*eta + 10826.6666666667 + 200106.666666667*(eta*eta) + 76160.0*(xi*xi) - 124444.444444444*eta*eta*eta - 24888.8888888889*xi*xi*xi, 2240.0*xi*(280*eta*xi/3 - 136*eta/3 - 49*xi + 13 + 100*(eta*eta)/3 + 35*(xi*xi)));
                      case 24:
                        return RealGradient(-340480.0*eta*xi + 71680.0*eta + 268800.0*eta*(xi*xi) + 56000.0*xi + 388266.666666667*xi*(eta*eta) - 8306.66666666667 - 158293.333333333*eta*eta - 81760.0*xi*xi + 99555.5555555555*(eta*eta*eta) + 31111.1111111111*(xi*xi*xi), 8960.0*xi*(-65*eta*xi/3 + 35*eta/3 + 14*xi - 4 - 20*eta*eta/3 - 10*xi*xi));
                      case 25:
                        return RealGradient(103040.0*eta*xi - 56000.0*eta - 33600.0*eta*xi*xi - 10453.3333333333*xi - 126933.333333333*xi*eta*eta + 5755.55555555556 + 125066.666666667*(eta*eta) - 2240.0*xi*xi - 77777.7777777778*eta*eta*eta + 4977.77777777778*(xi*xi*xi), 1120.0*xi*(170*eta*xi/3 - 134*eta/3 - 23*xi + 10 + 125*(eta*eta)/3 + 10*(xi*xi)));
                      case 26:
                        return RealGradient(-165760.0*eta*xi + 24640.0*eta + 168000.0*eta*(xi*xi) + 29866.6666666667*xi + 164266.666666667*xi*(eta*eta) - 3297.77777777778 - 46293.3333333333*eta*eta - 51520.0*xi*xi + 24888.8888888889*(eta*eta*eta) + 19911.1111111111*(xi*xi*xi), 2240.0*xi*(-110*eta*xi/3 + 71*eta/3 + 36*xi - 12 - 20*eta*eta/3 - 25*xi*xi));
                      case 27:
                        return RealGradient(-22400.0*eta*xi + 39200.0*eta - 16800.0*eta*xi*xi - 2986.66666666667*xi + 59733.3333333333*xi*(eta*eta) - 2955.55555555556 - 108266.666666667*eta*eta + 12320.0*(xi*xi) + 77777.7777777778*(eta*eta*eta) - 4977.77777777778*xi*xi*xi, 1120.0*xi*(-80*eta*xi/3 + 116*eta/3 + 5*xi - 7 - 125*eta*eta/3 + 5*(xi*xi)));
                      case 28:
                        return RealGradient(105280.0*eta*xi - 32480.0*eta - 16800.0*eta*xi*xi - 7466.66666666667*xi - 164266.666666667*xi*eta*eta + 3017.77777777778 + 79893.3333333333*(eta*eta) - 19040.0*xi*xi - 52888.8888888889*eta*eta*eta + 24888.8888888889*(xi*xi*xi), 1120.0*xi*(220*eta*xi/3 - 178*eta/3 - 27*xi + 17 + 85*(eta*eta)/3 + 5*(xi*xi)));
                      case 29:
                        return RealGradient(165760.0*eta*xi + 64960.0*eta - 168000.0*eta*xi*xi - 46293.3333333333*xi - 119466.666666667*xi*eta*eta - 5351.11111111111 - 173226.666666667*eta*eta + 80640.0*(xi*xi) + 124444.444444444*(eta*eta*eta) - 27377.7777777778*xi*xi*xi, 2240.0*xi*(80*eta*xi/3 + 64*eta/3 - 23*xi - 1 - 100*eta*eta/3 + 25*(xi*xi)));
                      case 30:
                        return RealGradient(4480.0*eta*xi - 11200.0*eta + 33600.0*eta*(xi*xi) + 5226.66666666667*xi - 29866.6666666667*xi*eta*eta + 715.555555555556 + 32853.3333333333*(eta*eta) - 25760.0*xi*xi - 24888.8888888889*eta*eta*eta + 21155.5555555556*(xi*xi*xi), 2240.0*xi*(20*eta*xi/3 - 44*eta/3 - xi + 5 + 20*(eta*eta)/3 - 5*xi*xi));
                      case 31:
                        return RealGradient(116480.0*eta*xi - 3360.0*eta + 16800.0*eta*(xi*xi) - 11946.6666666667*xi - 149333.333333333*xi*eta*eta + 964.444444444444 - 8213.33333333333*eta*eta - 30240.0*xi*xi + 15555.5555555556*(eta*eta*eta) + 27377.7777777778*(xi*xi*xi), 1120.0*xi*(200*eta*xi/3 - 20*eta/3 - 17*xi + 3 - 25*eta*eta/3 - 5*xi*xi));
                      case 32:
                        return RealGradient(-24640.0*eta*xi + 1120.0*eta + 16800.0*eta*(xi*xi) + 1493.33333333333*xi + 29866.6666666667*xi*(eta*eta) - 155.555555555556 + 746.666666666667*(eta*eta) + 5600.0*(xi*xi) - 3111.11111111111*eta*eta*eta - 9955.55555555556*xi*xi*xi, 1120.0*xi*(-40*eta*xi/3 - 2*eta/3 + 11*xi - 1 + 5*(eta*eta)/3 - 5*xi*xi));
                      case 33:
                        return RealGradient(-67200.0*eta*xi + 3360.0*eta + 50400.0*eta*(xi*xi) + 12693.3333333333*xi + 67200.0*xi*(eta*eta) - 622.222222222222 - 3360.0*eta*eta - 19040.0*xi*xi + 7466.66666666667*(xi*xi*xi), 3360.0*xi*(-10*eta*xi + 2*eta + 6*xi - 1 - 5*xi*xi));
                      case 34:
                        return RealGradient(47040.0*eta*xi - 7840.0*eta - 50400.0*eta*xi*xi - 9706.66666666667*xi - 44800.0*xi*eta*eta + 1057.77777777778 + 15680.0*(eta*eta) + 20160.0*(xi*xi) - 9333.33333333333*eta*eta*eta - 11200.0*xi*xi*xi, 1120.0*xi*(20*eta*xi - 8*eta - 17*xi + 3 + 5*(eta*eta) + 15*(xi*xi)));
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
      } // end case FIFTH

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
