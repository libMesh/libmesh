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
// https://web.archive.org/web/20240719125223/https://www.dealii.org/reports/nedelec/nedelec.pdf
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

  const Order totalorder = order + add_p_level*elem->p_level();
  libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

  const char sign = i >= totalorder * elem->n_edges() || elem->edge_orientation(i / totalorder) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / totalorder * 2 + 1) * totalorder - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (totalorder)
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
                  return sign * RealGradient(-81.*eta/4. - 9.*xi + 162.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 135.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 324.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99./4. + 81.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 270.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 45.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 1:
                  return sign * RealGradient(27.*eta/8. + 15.*xi/4. - 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 135.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 135.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 75.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 51./8. - 27.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 15.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 75.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 15.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2., 0.);
                case 2:
                  return sign * RealGradient(-27.*eta/4. - 6.*xi + 108.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 135.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 216.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 45./4. + 27.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 270.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 15.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 3:
                  return sign * RealGradient(0., 27.*xi/4. - 54.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 216.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 45.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 27./4. - 180.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 54.*(xi + 1.)*(xi + 1.)/(2.*2.) + 45.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 4:
                  return sign * RealGradient(0., -9.*xi/8. + 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 90.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 75.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 45.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 9./8. + 90.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 75.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 9.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 15.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2.);
                case 5:
                  return sign * RealGradient(0., 9.*xi/4. - 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 144.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 45.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 9./4. - 180.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 18.*(xi + 1.)*(xi + 1.)/(2.*2.) + 15.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 6:
                  return sign * RealGradient(-9.*eta/4. + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 45.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 144.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 9./4. + 18.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 180.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 15.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 7:
                  return sign * RealGradient(9.*eta/8. - 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 45.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 90.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 75.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 9./8. - 9.*(eta + 1.)*(eta + 1.)/(2.*2.) - 90.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 75.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 15.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2., 0.);
                case 8:
                  return sign * RealGradient(-27.*eta/4. + 54.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 45.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 216.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 27./4. + 54.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 180.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 45.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 9:
                  return sign * RealGradient(0., 6.*eta + 27.*xi/4. - 108.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 216.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 135.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 45./4. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 270.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 27.*(xi + 1.)*(xi + 1.)/(2.*2.) + 15.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 10:
                  return sign * RealGradient(0., -15.*eta/4. - 27.*xi/8. + 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 135.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 75.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 135.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 51./8. + 15.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 75.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 27.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 15.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2.);
                case 11:
                  return sign * RealGradient(0., 9.*eta + 81.*xi/4. - 162.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 324.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 135.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 99./4. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 270.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 81.*(xi + 1.)*(xi + 1.)/(2.*2.) + 45.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 12:
                  return RealGradient(0., 18.*xi - 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 324.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 18. - 270.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 81.*(xi + 1.)*(xi + 1.)/(2.*2.) + 45.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 13:
                  return RealGradient(0., -9.*xi/2. + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 216.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 30.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 9./2. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 54.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 45.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                case 14:
                  return RealGradient(-18.*eta + 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 324.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18. + 81.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 270.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 45.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 15:
                  return RealGradient(9.*eta/2. - 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 30.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 216.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 9./2. - 54.*(eta + 1.)*(eta + 1.)/(2.*2.) - 180.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 45.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                case 16:
                  return RealGradient(-6.*eta + 96.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 216.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 6. + 27.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 270.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 15.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 17:
                  return RealGradient(3.*eta/2. - 24.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 30.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 144.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3./2. - 18.*(eta + 1.)*(eta + 1.)/(2.*2.) - 180.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 15.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                case 18:
                  return RealGradient(0., 6.*xi - 96.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 216.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 6. - 270.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 27.*(xi + 1.)*(xi + 1.)/(2.*2.) + 15.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 19:
                  return RealGradient(0., -3.*xi/2. + 24.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 144.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 120.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 30.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 3./2. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 15.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                case 20:
                  return RealGradient(2.*eta + 2. - 9.*(eta + 1.)*(eta + 1.)/(2.*2.) + 5.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                case 21:
                  return RealGradient(0., 2.*xi + 2. - 9.*(xi + 1.)*(xi + 1.)/(2.*2.) + 5.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 22:
                  return RealGradient(0., -xi/2. - 1./2. + 6.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                case 23:
                  return RealGradient(-eta/2. - 1./2. + 6.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 5.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
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
                  return sign * RealGradient(120.*eta*xi - 54.*eta - 45.*eta*xi*xi - 36.*xi - 90.*xi*eta*eta + 9. + 90.*(eta*eta) + 30.*(xi*xi) - 45.*eta*eta*eta, -60.*eta*xi + 90.*eta*(xi*xi) + 18.*xi + 45.*xi*(eta*eta) - 60.*xi*xi + 45.*(xi*xi*xi));
                case 1:
                  return sign * RealGradient(-75.*eta*xi/2. - 9.*eta/4. + 45.*eta*(xi*xi)/2. + 15.*xi + 45.*xi*(eta*eta)/2. - 3./2. + 15.*(eta*eta) - 15.*xi*xi - 45.*eta*eta*eta/4., -45.*eta*xi*xi/2. - 21.*xi/4. + 45.*xi*(eta*eta)/4. + 105.*(xi*xi)/4. - 45.*xi*xi*xi/2.);
                case 2:
                  return sign * RealGradient(30.*eta*xi - 3.*eta - 45.*eta*xi*xi - 24.*xi + 3. + 30.*(xi*xi), 9.*xi - 45.*xi*xi + 45.*(xi*xi*xi));
                case 3:
                  return sign * RealGradient(30.*eta*xi - 3.*eta - 45.*eta*xi*xi, 9.*xi - 45.*xi*xi + 45.*(xi*xi*xi));
                case 4:
                  return sign * RealGradient(45.*eta*xi/2. - 6.*eta - 45.*eta*xi*xi/4. - 45.*xi*eta*eta + 75.*(eta*eta)/4. - 45.*eta*eta*eta/4., -45.*eta*xi/2. + 45.*eta*(xi*xi) + 6.*xi + 45.*xi*(eta*eta)/4. - 75.*xi*xi/4. + 45.*(xi*xi*xi)/4.);
                case 5:
                  return sign * RealGradient(-9.*eta + 45.*(eta*eta) - 45.*eta*eta*eta, -30.*eta*xi + 3.*xi + 45.*xi*(eta*eta));
                case 6:
                  return sign * RealGradient(-9.*eta + 45.*(eta*eta) - 45.*eta*eta*eta, -30.*eta*xi + 24.*eta + 3.*xi + 45.*xi*(eta*eta) - 3. - 30.*eta*eta);
                case 7:
                  return sign * RealGradient(21.*eta/4. - 45.*eta*xi*xi/4. + 45.*xi*(eta*eta)/2. - 105.*eta*eta/4. + 45.*(eta*eta*eta)/2., 75.*eta*xi/2. - 15.*eta - 45.*eta*xi*xi/2. + 9.*xi/4. - 45.*xi*eta*eta/2. + 3./2. + 15.*(eta*eta) - 15.*xi*xi + 45.*(xi*xi*xi)/4.);
                case 8:
                  return sign * RealGradient(60.*eta*xi - 18.*eta - 45.*eta*xi*xi - 90.*xi*eta*eta + 60.*(eta*eta) - 45.*eta*eta*eta, -120.*eta*xi + 36.*eta + 90.*eta*(xi*xi) + 54.*xi + 45.*xi*(eta*eta) - 9. - 30.*eta*eta - 90.*xi*xi + 45.*(xi*xi*xi));
                case 9:
                  return RealGradient(-300.*eta*xi + 180.*eta + 90.*eta*(xi*xi) + 360.*xi*(eta*eta) - 450.*eta*eta + 270.*(eta*eta*eta), 300.*eta*xi - 360.*eta*xi*xi - 60.*xi - 270.*xi*eta*eta + 150.*(xi*xi) - 90.*xi*xi*xi);
                case 10:
                  return RealGradient(300.*eta*xi - 60.*eta - 270.*eta*xi*xi - 360.*xi*eta*eta + 150.*(eta*eta) - 90.*eta*eta*eta, -300.*eta*xi + 360.*eta*(xi*xi) + 180.*xi + 90.*xi*(eta*eta) - 450.*xi*xi + 270.*(xi*xi*xi));
                case 11:
                  return RealGradient(360.*eta*xi - 60.*eta - 180.*eta*xi*xi - 360.*xi*eta*eta + 60.*(eta*eta), -120.*eta*xi + 360.*eta*(xi*xi) + 60.*xi - 240.*xi*xi + 180.*(xi*xi*xi));
                case 12:
                  return RealGradient(-240.*eta*xi + 30.*eta + 270.*eta*(xi*xi) + 180.*xi*(eta*eta) - 30.*eta*eta, 60.*eta*xi - 180.*eta*xi*xi - 90.*xi + 360.*(xi*xi) - 270.*xi*xi*xi);
                case 13:
                  return RealGradient(60.*eta*xi - 90.*eta - 180.*xi*eta*eta + 360.*(eta*eta) - 270.*eta*eta*eta, -240.*eta*xi + 180.*eta*(xi*xi) + 30.*xi + 270.*xi*(eta*eta) - 30.*xi*xi);
                case 14:
                  return RealGradient(-120.*eta*xi + 60.*eta + 360.*xi*(eta*eta) - 240.*eta*eta + 180.*(eta*eta*eta), 360.*eta*xi - 360.*eta*xi*xi - 60.*xi - 180.*xi*eta*eta + 60.*(xi*xi));
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
                  return sign * RealGradient(-64.*eta - 30.*xi + 960.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1920.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1120.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 4800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 86. + 480.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 7200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 120.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 9600.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 640.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 70.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 2450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 280.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                case 1:
                  return sign * RealGradient(272.*eta/27. + 95.*xi/9. - 3040.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 6880.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 12320.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 3800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 15200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 6650.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. + 523./27. - 680.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 8600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 15400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 430.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 34400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 15050.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. + 2720.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 61600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 770.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. + 26950.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/27. - 1190.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/27., 0.);
                case 2:
                  return sign * RealGradient(-128.*eta/27. - 50.*xi/9. + 1600.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 5440.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 12320.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 2000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 8000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 3500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. - 262./27. + 320.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 6800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 15400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 340.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 27200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 11900.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. - 1280.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 61600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 770.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. - 26950.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/27. + 560.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/27., 0.);
                case 3:
                  return sign * RealGradient(16.*eta + 15.*xi - 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1440.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1050.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 29. - 120.*(eta + 1.)*(eta + 1.)/(2.*2.) - 5400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 90.*(xi + 1.)*(xi + 1.)/(2.*2.) + 7200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 160.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 70.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 2450.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 70.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                case 4:
                  return sign * RealGradient(0., -16.*xi + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 480.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 280.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 16. + 3600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 240.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 480.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 280.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                case 5:
                  return sign * RealGradient(0., 68.*xi/27. - 760.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 1900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 3800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/3. + 6650.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 1720.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/9. - 3080.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 68./27. - 4300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 8600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/3. - 15050.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 340.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 7700.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 15400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 26950.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/27. + 680.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 1190.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/27.);
                case 6:
                  return sign * RealGradient(0., -32.*xi/27. + 400.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 1000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 2000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/3. - 3500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 1360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/9. + 3080.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 32./27. + 3400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 6800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/3. + 11900.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 160.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 7700.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 15400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 26950.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/27. - 320.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 560.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/27.);
                case 7:
                  return sign * RealGradient(0., 4.*xi - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1050.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 280.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4. - 2700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 120.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 70.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                case 8:
                  return sign * RealGradient(-4.*eta + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 280.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 1800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1050.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 4. + 60.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 2700.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 5400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3150.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 120.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 70.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                case 9:
                  return sign * RealGradient(32.*eta/27. - 400.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 1360.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 3080.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 1000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 2000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/3. + 3500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. + 32./27. - 160.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 3400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7700.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. + 6800.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/3. - 11900.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. + 320.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 15400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 26950.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/27. - 560.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/27., 0.);
                case 10:
                  return sign * RealGradient(-68.*eta/27. + 760.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 1720.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 3080.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 1900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 3800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/3. - 6650.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. - 68./27. + 340.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 4300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 7700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. - 8600.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/3. + 15050.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. - 680.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 15400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 26950.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/27. + 1190.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/27., 0.);
                case 11:
                  return sign * RealGradient(16.*eta - 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 480.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 280.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 3600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 16. - 240.*(eta + 1.)*(eta + 1.)/(2.*2.) - 3600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 7200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 480.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 280.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                case 12:
                  return sign * RealGradient(0., -15.*eta - 16.*xi + 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1050.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1440.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 1120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 29. + 90.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 5400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 120.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 70.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 160.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 70.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                case 13:
                  return sign * RealGradient(0., 50.*eta/9. + 128.*xi/27. - 1600.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 2000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 8000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 3500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 5440.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/9. - 12320.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 262./27. - 340.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 6800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 27200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 11900.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 320.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 15400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. + 770.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 61600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 26950.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/27. + 1280.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 560.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/27.);
                case 14:
                  return sign * RealGradient(0., -95.*eta/9. - 272.*xi/27. + 3040.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 3800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 15200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 6650.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 6880.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/9. + 12320.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 523./27. + 430.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 8600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 34400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 15050.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 680.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 15400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. - 770.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 61600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 26950.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/27. - 2720.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 1190.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/27.);
                case 15:
                  return sign * RealGradient(0., 30.*eta + 64.*xi - 960.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1920.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 1120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 86. - 120.*(eta + 1.)*(eta + 1.)/(2.*2.) - 7200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 9600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 480.*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 70.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 640.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 280.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                case 16:
                  return RealGradient(0., 52.*xi - 780.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3480.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1560.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 910.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 52. - 6960.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 9600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 464.*(xi + 1.)*(xi + 1.)/(2.*2.) + 4060.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 640.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 280.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                case 17:
                  return RealGradient(0., 12.*xi - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1680.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 210.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 12. - 3360.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 224.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1960.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 480.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 280.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                case 18:
                  return RealGradient(0., 16.*xi - 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 240.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 480.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 280.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 16. - 480.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 32.*(xi + 1.)*(xi + 1.)/(2.*2.) + 280.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)));
                case 19:
                  return RealGradient(-52.*eta + 780.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1560.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 910.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3480.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 4800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 52. + 464.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 6960.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4060.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 9600.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 640.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 280.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                case 20:
                  return RealGradient(-12.*eta + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 210.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1680.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 3600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 12. + 224.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 3360.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1960.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 7200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 480.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 280.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                case 21:
                  return RealGradient(-16.*eta + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 480.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 280.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 240.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 16. + 32.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 480.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 280.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                case 22:
                  return RealGradient(13.*eta - 390.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1170.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 910.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1740.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1050.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 13. - 116.*(eta + 1.)*(eta + 1.)/(2.*2.) - 5220.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4060.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 7200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 160.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 70.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                case 23:
                  return RealGradient(3.*eta - 90.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 210.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 840.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 1800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1050.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 3. - 56.*(eta + 1.)*(eta + 1.)/(2.*2.) - 2520.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1960.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 5400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 120.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2450.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 70.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                case 24:
                  return RealGradient(4.*eta - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 360.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 280.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 4. - 8.*(eta + 1.)*(eta + 1.)/(2.*2.) - 360.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 280.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                case 25:
                  return RealGradient(0., -13.*xi + 390.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1740.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1050.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1170.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 910.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 13. + 5220.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 116.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4060.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 160.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 70.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                case 26:
                  return RealGradient(0., -3.*xi + 90.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 840.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1050.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 210.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3. + 2520.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 56.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1960.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 120.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 70.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                case 27:
                  return RealGradient(0., -4.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 280.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4. + 360.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 8.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 280.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.));
                case 28:
                  return RealGradient(12.*eta - 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 171.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 240.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 105.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 12. - 114.*(eta + 1.)*(eta + 1.)/(2.*2.) + 160.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 70.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                case 29:
                  return RealGradient(-6.*eta + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 171.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 240.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 105.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 6. + 57.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 80.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 35.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                case 30:
                  return RealGradient(0., 12.*xi - 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 171.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 240.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 105.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 12. - 114.*(xi + 1.)*(xi + 1.)/(2.*2.) + 160.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 70.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                case 31:
                  return RealGradient(0., -6.*xi + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 171.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 240.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 105.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 6. + 57.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 80.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 35.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                case 32:
                  return RealGradient(0., 2.*xi - 6.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 81.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 105.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 2. - 54.*(xi + 1.)*(xi + 1.)/(2.*2.) + 120.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 70.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                case 33:
                  return RealGradient(0., -xi + 6.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 81.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 105.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1. + 27.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 60.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 35.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                case 34:
                  return RealGradient(2.*eta - 6.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 81.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 105.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 2. - 54.*(eta + 1.)*(eta + 1.)/(2.*2.) + 120.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 70.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                case 35:
                  return RealGradient(-eta + 6.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 81.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 105.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 1. + 27.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 60.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 35.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                case 36:
                  return RealGradient(0., 6.*xi - 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 18.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 6. - 12.*(xi + 1.)*(xi + 1.)/(2.*2.));
                case 37:
                  return RealGradient(-6.*eta + 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 18.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 6. + 12.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                case 38:
                  return RealGradient(3.*eta - 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 18.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 3. - 6.*(eta + 1.)*(eta + 1.)/(2.*2.), 0.);
                case 39:
                  return RealGradient(0., -3.*xi + 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 18.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 3. + 6.*((xi + 1.)*(xi + 1.)/(2.*2.)));
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
                  return sign * RealGradient(720.*eta*xi - 160.*eta - 840.*eta*xi*xi + 224.*eta*(xi*xi*xi) - 120.*xi - 1260.*xi*eta*eta + 672.*xi*(eta*eta*eta) + 16. + 480.*(eta*eta) + 672.*(eta*eta)*(xi*xi) + 240.*(xi*xi) - 560.*eta*eta*eta - 140.*xi*xi*xi + 224.*(eta*eta*eta*eta), -240.*eta*xi + 840.*eta*(xi*xi) - 672.*eta*xi*xi*xi + 40.*xi + 420.*xi*(eta*eta) - 224.*xi*eta*eta*eta - 672.*eta*eta*xi*xi - 240.*xi*xi + 420.*(xi*xi*xi) - 224.*xi*xi*xi*xi);
                case 1:
                  return sign * RealGradient(-200.*eta*xi - 76.*eta/27. + 2884.*eta*(xi*xi)/9. - 2464.*eta*xi*xi*xi/27. + 380.*xi/9. + 784.*xi*(eta*eta)/3. - 896.*xi*eta*eta*eta/9. - 68./27. + 608.*(eta*eta)/9. - 2240.*eta*eta*xi*xi/9. - 860.*xi*xi/9. - 3472.*eta*eta*eta/27. + 1540.*(xi*xi*xi)/27. + 1792.*(eta*eta*eta*eta)/27., 176.*eta*xi/9. - 2128.*eta*xi*xi/9. + 2240.*eta*(xi*xi*xi)/9. - 296.*xi/27. + 448.*xi*(eta*eta)/9. - 1792.*xi*eta*eta*eta/27. + 896.*(eta*eta)*(xi*xi)/9. + 800.*(xi*xi)/9. - 168.*xi*xi*xi + 2464.*(xi*xi*xi*xi)/27.);
                case 2:
                  return sign * RealGradient(-16.*eta*xi + 112.*eta/27. - 1120.*eta*xi*xi/9. + 2464.*eta*(xi*xi*xi)/27. - 200.*xi/9. + 476.*xi*(eta*eta)/3. - 1120.*xi*eta*eta*eta/9. + 32./27. - 104.*eta*eta/9. + 224.*(eta*eta)*(xi*xi)/9. + 680.*(xi*xi)/9. - 56.*eta*eta*eta/27. - 1540.*xi*xi*xi/27. + 224.*(eta*eta*eta*eta)/27., 112.*eta*xi/9. - 392.*eta*xi*xi/9. - 224.*eta*xi*xi*xi/9. + 152.*xi/27. - 196.*xi*eta*eta/9. - 224.*xi*eta*eta*eta/27. + 1120.*(eta*eta)*(xi*xi)/9. - 512.*xi*xi/9. + 140.*(xi*xi*xi) - 2464.*xi*xi*xi*xi/27.);
                case 3:
                  return sign * RealGradient(-72.*eta*xi + 4.*eta + 252.*eta*(xi*xi) - 224.*eta*xi*xi*xi + 60.*xi - 4. - 180.*xi*xi + 140.*(xi*xi*xi), -16.*xi + 144.*(xi*xi) - 336.*xi*xi*xi + 224.*(xi*xi*xi*xi));
                case 4:
                  return sign * RealGradient(-72.*eta*xi + 4.*eta + 252.*eta*(xi*xi) - 224.*eta*xi*xi*xi, -16.*xi + 144.*(xi*xi) - 336.*xi*xi*xi + 224.*(xi*xi*xi*xi));
                case 5:
                  return sign * RealGradient(-80.*eta*xi + 8.*eta + 448.*eta*(xi*xi)/3. - 1792.*eta*xi*xi*xi/27. + 784.*xi*(eta*eta)/3. - 448.*xi*eta*eta*eta/3. - 32.*eta*eta - 896.*eta*eta*xi*xi/3. + 280.*(eta*eta*eta)/9. - 224.*eta*eta*eta*eta/27., 56.*eta*xi - 896.*eta*xi*xi/3. + 896.*eta*(xi*xi*xi)/3. - 12.*xi - 140.*xi*eta*eta/3. + 224.*xi*(eta*eta*eta)/27. + 448.*(eta*eta)*(xi*xi)/3. + 80.*(xi*xi) - 1232.*xi*xi*xi/9. + 1792.*(xi*xi*xi*xi)/27.);
                case 6:
                  return sign * RealGradient(-56.*eta*xi + 12.*eta + 140.*eta*(xi*xi)/3. - 224.*eta*xi*xi*xi/27. + 896.*xi*(eta*eta)/3. - 896.*xi*eta*eta*eta/3. - 80.*eta*eta - 448.*eta*eta*xi*xi/3. + 1232.*(eta*eta*eta)/9. - 1792.*eta*eta*eta*eta/27., 80.*eta*xi - 784.*eta*xi*xi/3. + 448.*eta*(xi*xi*xi)/3. - 8.*xi - 448.*xi*eta*eta/3. + 1792.*xi*(eta*eta*eta)/27. + 896.*(eta*eta)*(xi*xi)/3. + 32.*(xi*xi) - 280.*xi*xi*xi/9. + 224.*(xi*xi*xi*xi)/27.);
                case 7:
                  return sign * RealGradient(16.*eta - 144.*eta*eta + 336.*(eta*eta*eta) - 224.*eta*eta*eta*eta, 72.*eta*xi - 4.*xi - 252.*xi*eta*eta + 224.*xi*(eta*eta*eta));
                case 8:
                  return sign * RealGradient(16.*eta - 144.*eta*eta + 336.*(eta*eta*eta) - 224.*eta*eta*eta*eta, 72.*eta*xi - 60.*eta - 4.*xi - 252.*xi*eta*eta + 224.*xi*(eta*eta*eta) + 4. + 180.*(eta*eta) - 140.*eta*eta*eta);
                case 9:
                  return sign * RealGradient(-112.*eta*xi/9. - 152.*eta/27. + 196.*eta*(xi*xi)/9. + 224.*eta*(xi*xi*xi)/27. + 392.*xi*(eta*eta)/9. + 224.*xi*(eta*eta*eta)/9. + 512.*(eta*eta)/9. - 1120.*eta*eta*xi*xi/9. - 140.*eta*eta*eta + 2464.*(eta*eta*eta*eta)/27., 16.*eta*xi + 200.*eta/9. - 476.*eta*xi*xi/3. + 1120.*eta*(xi*xi*xi)/9. - 112.*xi/27. + 1120.*xi*(eta*eta)/9. - 2464.*xi*eta*eta*eta/27. - 32./27. - 680.*eta*eta/9. - 224.*eta*eta*xi*xi/9. + 104.*(xi*xi)/9. + 1540.*(eta*eta*eta)/27. + 56.*(xi*xi*xi)/27. - 224.*xi*xi*xi*xi/27.);
                case 10:
                  return sign * RealGradient(-176.*eta*xi/9. + 296.*eta/27. - 448.*eta*xi*xi/9. + 1792.*eta*(xi*xi*xi)/27. + 2128.*xi*(eta*eta)/9. - 2240.*xi*eta*eta*eta/9. - 800.*eta*eta/9. - 896.*eta*eta*xi*xi/9. + 168.*(eta*eta*eta) - 2464.*eta*eta*eta*eta/27., 200.*eta*xi - 380.*eta/9. - 784.*eta*xi*xi/3. + 896.*eta*(xi*xi*xi)/9. + 76.*xi/27. - 2884.*xi*eta*eta/9. + 2464.*xi*(eta*eta*eta)/27. + 68./27. + 860.*(eta*eta)/9. + 2240.*(eta*eta)*(xi*xi)/9. - 608.*xi*xi/9. - 1540.*eta*eta*eta/27. + 3472.*(xi*xi*xi)/27. - 1792.*xi*xi*xi*xi/27.);
                case 11:
                  return sign * RealGradient(240.*eta*xi - 40.*eta - 420.*eta*xi*xi + 224.*eta*(xi*xi*xi) - 840.*xi*eta*eta + 672.*xi*(eta*eta*eta) + 240.*(eta*eta) + 672.*(eta*eta)*(xi*xi) - 420.*eta*eta*eta + 224.*(eta*eta*eta*eta), -720.*eta*xi + 120.*eta + 1260.*eta*(xi*xi) - 672.*eta*xi*xi*xi + 160.*xi + 840.*xi*(eta*eta) - 224.*xi*eta*eta*eta - 16. - 240.*eta*eta - 672.*eta*eta*xi*xi - 480.*xi*xi + 140.*(eta*eta*eta) + 560.*(xi*xi*xi) - 224.*xi*xi*xi*xi);
                case 12:
                  return RealGradient(-3240.*eta*xi + 960.*eta + 3024.*eta*(xi*xi) - 672.*eta*xi*xi*xi + 9072.*xi*(eta*eta) - 6048.*xi*eta*eta*eta - 4320.*eta*eta - 4032.*eta*eta*xi*xi + 6048.*(eta*eta*eta) - 2688.*eta*eta*eta*eta, 2160.*eta*xi - 6048.*eta*xi*xi + 4032.*eta*(xi*xi*xi) - 240.*xi - 4536.*xi*eta*eta + 2688.*xi*(eta*eta*eta) + 6048.*(eta*eta)*(xi*xi) + 1080.*(xi*xi) - 1512.*xi*xi*xi + 672.*(xi*xi*xi*xi));
                case 13:
                  return RealGradient(2160.*eta*xi - 240.*eta - 4536.*eta*xi*xi + 2688.*eta*(xi*xi*xi) - 6048.*xi*eta*eta + 4032.*xi*(eta*eta*eta) + 1080.*(eta*eta) + 6048.*(eta*eta)*(xi*xi) - 1512.*eta*eta*eta + 672.*(eta*eta*eta*eta), -3240.*eta*xi + 9072.*eta*(xi*xi) - 6048.*eta*xi*xi*xi + 960.*xi + 3024.*xi*(eta*eta) - 672.*xi*eta*eta*eta - 4032.*eta*eta*xi*xi - 4320.*xi*xi + 6048.*(xi*xi*xi) - 2688.*xi*xi*xi*xi);
                case 14:
                  return RealGradient(-1944.*eta*xi + 144.*eta + 4536.*eta*(xi*xi) - 2016.*eta*xi*xi*xi + 2016.*xi*(eta*eta) - 144.*eta*eta - 4032.*eta*eta*xi*xi, 432.*eta*xi - 3024.*eta*xi*xi + 4032.*eta*(xi*xi*xi) - 216.*xi + 1728.*(xi*xi) - 3528.*xi*xi*xi + 2016.*(xi*xi*xi*xi));
                case 15:
                  return RealGradient(1152.*eta*xi - 72.*eta - 3528.*eta*xi*xi + 2688.*eta*(xi*xi*xi) - 1008.*xi*eta*eta + 72.*(eta*eta) + 2016.*(eta*eta)*(xi*xi), -216.*eta*xi + 1512.*eta*(xi*xi) - 2016.*eta*xi*xi*xi + 288.*xi - 2304.*xi*xi + 4704.*(xi*xi*xi) - 2688.*xi*xi*xi*xi);
                case 16:
                  return RealGradient(-216.*eta*xi + 288.*eta + 1512.*xi*(eta*eta) - 2016.*xi*eta*eta*eta - 2304.*eta*eta + 4704.*(eta*eta*eta) - 2688.*eta*eta*eta*eta, 1152.*eta*xi - 1008.*eta*xi*xi - 72.*xi - 3528.*xi*eta*eta + 2688.*xi*(eta*eta*eta) + 2016.*(eta*eta)*(xi*xi) + 72.*(xi*xi));
                case 17:
                  return RealGradient(432.*eta*xi - 216.*eta - 3024.*xi*eta*eta + 4032.*xi*(eta*eta*eta) + 1728.*(eta*eta) - 3528.*eta*eta*eta + 2016.*(eta*eta*eta*eta), -1944.*eta*xi + 2016.*eta*(xi*xi) + 144.*xi + 4536.*xi*(eta*eta) - 2016.*xi*eta*eta*eta - 4032.*eta*eta*xi*xi - 144.*xi*xi);
                case 18:
                  return RealGradient(-1332.*eta*xi + 216.*eta + 1638.*eta*(xi*xi) - 504.*eta*xi*xi*xi + 4788.*xi*(eta*eta) - 3528.*xi*eta*eta*eta - 1098.*eta*eta - 3024.*eta*eta*xi*xi + 1554.*(eta*eta*eta) - 672.*eta*eta*eta*eta, 1044.*eta*xi - 4032.*eta*xi*xi + 3024.*eta*(xi*xi*xi) - 144.*xi - 1638.*xi*eta*eta + 672.*xi*(eta*eta*eta) + 3528.*(eta*eta)*(xi*xi) + 774.*(xi*xi) - 1134.*xi*xi*xi + 504.*(xi*xi*xi*xi));
                case 19:
                  return RealGradient(1044.*eta*xi - 144.*eta - 1638.*eta*xi*xi + 672.*eta*(xi*xi*xi) - 4032.*xi*eta*eta + 3024.*xi*(eta*eta*eta) + 774.*(eta*eta) + 3528.*(eta*eta)*(xi*xi) - 1134.*eta*eta*eta + 504.*(eta*eta*eta*eta), -1332.*eta*xi + 4788.*eta*(xi*xi) - 3528.*eta*xi*xi*xi + 216.*xi + 1638.*xi*(eta*eta) - 504.*xi*eta*eta*eta - 3024.*eta*eta*xi*xi - 1098.*xi*xi + 1554.*(xi*xi*xi) - 672.*xi*xi*xi*xi);
                case 20:
                  return RealGradient(-216.*eta*xi - 48.*eta + 504.*eta*(xi*xi) - 168.*eta*xi*xi*xi - 756.*xi*eta*eta + 1008.*xi*(eta*eta*eta) + 720.*(eta*eta) - 1344.*eta*eta*eta + 672.*(eta*eta*eta*eta), -360.*eta*xi + 504.*eta*(xi*xi) + 12.*xi + 1008.*xi*(eta*eta) - 672.*xi*eta*eta*eta - 1008.*eta*eta*xi*xi + 72.*(xi*xi) - 252.*xi*xi*xi + 168.*(xi*xi*xi*xi));
                case 21:
                  return RealGradient(-216.*eta*xi + 66.*eta - 378.*eta*xi*xi + 672.*eta*(xi*xi*xi) + 2268.*xi*(eta*eta) - 2016.*xi*eta*eta*eta - 486.*eta*eta - 1512.*eta*eta*xi*xi + 756.*(eta*eta*eta) - 336.*eta*eta*eta*eta, 1188.*eta*xi - 2772.*eta*xi*xi + 1512.*eta*(xi*xi*xi) + 6.*xi - 1512.*xi*eta*eta + 336.*xi*(eta*eta*eta) + 2016.*(eta*eta)*(xi*xi) - 468.*xi*xi + 1134.*(xi*xi*xi) - 672.*xi*xi*xi*xi);
                case 22:
                  return RealGradient(1188.*eta*xi + 6.*eta - 1512.*eta*xi*xi + 336.*eta*(xi*xi*xi) - 2772.*xi*eta*eta + 1512.*xi*(eta*eta*eta) - 468.*eta*eta + 2016.*(eta*eta)*(xi*xi) + 1134.*(eta*eta*eta) - 672.*eta*eta*eta*eta, -216.*eta*xi + 2268.*eta*(xi*xi) - 2016.*eta*xi*xi*xi + 66.*xi - 378.*xi*eta*eta + 672.*xi*(eta*eta*eta) - 1512.*eta*eta*xi*xi - 486.*xi*xi + 756.*(xi*xi*xi) - 336.*xi*xi*xi*xi);
                case 23:
                  return RealGradient(-360.*eta*xi + 12.*eta + 1008.*eta*(xi*xi) - 672.*eta*xi*xi*xi + 504.*xi*(eta*eta) + 72.*(eta*eta) - 1008.*eta*eta*xi*xi - 252.*eta*eta*eta + 168.*(eta*eta*eta*eta), -216.*eta*xi - 756.*eta*xi*xi + 1008.*eta*(xi*xi*xi) - 48.*xi + 504.*xi*(eta*eta) - 168.*xi*eta*eta*eta + 720.*(xi*xi) - 1344.*xi*xi*xi + 672.*(xi*xi*xi*xi));
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
                  return sign * RealGradient(-625.*eta/4. - 75.*xi + 3750.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 13125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 17500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 7875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 22500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 52500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 52500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 875./4. + 1875.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 78750.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 105000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 47250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 525.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 183750.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 183750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 66150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 4375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 245000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 110250.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 700.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 245000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 88200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 4375.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 110250.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 315.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 1:
                  return sign * RealGradient(12125.*eta/512. + 2865.*xi/128. - 71625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. + 527625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 340375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 291375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 214875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/32. - 501375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 501375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. - 180495.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. + 22615./512. - 36375.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 1582875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 1021125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 874125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 21105.*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 3693375.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 3693375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 1329615.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 84875.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 2382625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 13615.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. + 2382625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. - 857745.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. - 84875.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 11655.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 30555.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/128., 0.);
                case 2:
                  return sign * RealGradient(-375.*eta/32. - 105.*xi/8. + 2625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/4. - 28875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 23625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 23625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 7875.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 18375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 18375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 6615.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/2. - 765./32. + 1125.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 86625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 70875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 70875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 1155.*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 202125.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 202125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 72765.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/4. - 2625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. - 165375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 59535.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/2. + 2625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. - 59535.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/4. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/8., 0.);
                case 3:
                  return sign * RealGradient(2125.*eta/512. + 705.*xi/128. - 17625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. + 233625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 242375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 291375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 52875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/32. - 123375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 123375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. - 44415.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. + 4775./512. - 6375.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 700875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 727125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 874125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 9345.*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 1635375.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 1635375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 588735.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 14875.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 1696625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 9695.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. + 1696625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. - 610785.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. - 14875.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 11655.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 5355.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/128., 0.);
                case 4:
                  return sign * RealGradient(-125.*eta/4. - 30.*xi + 1500.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 7875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 14000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 7875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 9000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 21000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 21000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7560.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 235./4. + 375.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 47250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 84000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 47250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 315.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 110250.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 110250.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 196000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 110250.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 560.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 196000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 70560.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 110250.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 315.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 315.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 5:
                  return sign * RealGradient(0., 125.*xi/4. - 750.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 9000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 31500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 42000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 2625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 3500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1575.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 125./4. - 31500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 110250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 147000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 66150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 750.*(xi + 1.)*(xi + 1.)/(2.*2.) + 42000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18900.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 147000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 196000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 2625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 66150.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 3500.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 6:
                  return sign * RealGradient(0., -2425.*xi/512. + 14325.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. - 42975.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/16. + 300825.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 100275.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. + 180495.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 105525.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/128. + 68075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 58275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. - 2425./512. + 316575.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 2216025.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 738675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16. - 1329615.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 7275.*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 204225.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/16. + 174825.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. + 1429575.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 476525.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. + 857745.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 50925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 1223775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 407925.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 16975.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/32. - 30555.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/128.);
                case 7:
                  return sign * RealGradient(0., 75.*xi/32. - 525.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/4. + 1575.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 11025.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 7350.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 6615.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/2. + 5775.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/8. - 4725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. + 75./32. - 17325.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 121275.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 40425.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 72765.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/4. - 225.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 14175.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 66150.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 59535.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/2. + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 33075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 59535.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/4. - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/8.);
                case 8:
                  return sign * RealGradient(0., -425.*xi/512. + 3525.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. - 10575.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/16. + 74025.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 24675.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. + 44415.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 46725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/128. + 48475.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 58275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. - 425./512. + 140175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 981225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 327075.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16. - 588735.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 1275.*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 145425.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/16. + 174825.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. + 1017975.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 339325.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. + 610785.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 8925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 1223775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 407925.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 2975.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/32. - 5355.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/128.);
                case 9:
                  return sign * RealGradient(0., 25.*xi/4. - 300.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 12600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 16800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 7560.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 1575.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1575.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 25./4. - 18900.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 66150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 88200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.) + 33600.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18900.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 117600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 156800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 70560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 66150.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 700.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 315.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 10:
                  return sign * RealGradient(-25.*eta/4. + 300.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1575.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1575.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 3600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 12600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 16800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7560.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 25./4. + 150.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 18900.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 33600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18900.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 66150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 88200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 117600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 66150.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 156800.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 70560.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 700.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 88200.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 315.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 11:
                  return sign * RealGradient(425.*eta/512. - 3525.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. + 46725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 48475.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 58275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 10575.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/16. - 74025.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 24675.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 44415.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. + 425./512. - 1275.*(eta + 1.)*(eta + 1.)/(2.*2.)/64. - 140175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 145425.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 174825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. + 981225.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 327075.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/16. + 588735.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 8925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 1017975.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 1223775.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 339325.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 610785.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. - 2975.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/32. - 407925.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 5355.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/128., 0.);
                case 12:
                  return sign * RealGradient(-75.*eta/32. + 525.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/4. - 5775.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 4725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 1575.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 11025.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 7350.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 6615.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/2. - 75./32. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 17325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 14175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. - 121275.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 40425.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 72765.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/4. - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 66150.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 59535.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/2. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 33075.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 59535.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/4. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/8., 0.);
                case 13:
                  return sign * RealGradient(2425.*eta/512. - 14325.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. + 105525.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 68075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 58275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 42975.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/16. - 300825.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 100275.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 180495.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. + 2425./512. - 7275.*(eta + 1.)*(eta + 1.)/(2.*2.)/64. - 316575.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 204225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 174825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. + 2216025.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 738675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/16. + 1329615.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 50925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 1429575.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 1223775.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 476525.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 857745.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/32. - 16975.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/32. - 407925.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. + 30555.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/128., 0.);
                case 14:
                  return sign * RealGradient(-125.*eta/4. + 750.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 3500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1575.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 9000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 31500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 42000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 125./4. + 750.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 31500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 42000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18900.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 110250.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 147000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 66150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 2625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 147000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 66150.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 196000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 88200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 3500.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 88200.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 15:
                  return sign * RealGradient(0., 30.*eta + 125.*xi/4. - 1500.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 9000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 21000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 21000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 7560.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 7875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 14000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 7875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 235./4. - 315.*(eta + 1.)*(eta + 1.)/(2.*2.) - 47250.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 110250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 110250.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 375.*(xi + 1.)*(xi + 1.)/(2.*2.) + 84000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 47250.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 196000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 196000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 70560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 875.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 110250.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 315.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 110250.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 315.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 16:
                  return sign * RealGradient(0., -705.*eta/128. - 2125.*xi/512. + 17625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. - 52875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 123375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 123375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. + 44415.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 233625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/128. + 242375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 291375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. - 4775./512. + 9345.*((eta + 1.)*(eta + 1.)/(2.*2.))/128. + 700875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 1635375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 1635375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. - 588735.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 6375.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 727125.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 874125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 9695.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 1696625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 1696625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. + 610785.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 14875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 2039625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 11655.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 14875.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. - 5355.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/128.);
                case 17:
                  return sign * RealGradient(0., 105.*eta/8. + 375.*xi/32. - 2625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/4. + 7875.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 18375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 18375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. - 6615.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/2. + 28875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/8. - 23625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 23625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. + 765./32. - 1155.*(eta + 1.)*(eta + 1.)/(2.*2.)/8. - 86625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 202125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 202125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. + 72765.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/4. - 1125.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 70875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 70875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. - 59535.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/2. + 2625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. + 165375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. + 59535.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/4. - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/8.);
                case 18:
                  return sign * RealGradient(0., -2865.*eta/128. - 12125.*xi/512. + 71625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. - 214875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 501375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 501375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. + 180495.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 527625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/128. + 340375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 291375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. - 22615./512. + 21105.*((eta + 1.)*(eta + 1.)/(2.*2.))/128. + 1582875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 3693375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 3693375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. - 1329615.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 36375.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 1021125.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 874125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 13615.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 2382625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 2382625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. + 857745.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/32. - 84875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 2039625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 11655.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 84875.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. - 30555.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/128.);
                case 19:
                  return sign * RealGradient(0., 75.*eta + 625.*xi/4. - 3750.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 22500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 52500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 52500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 13125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 17500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 7875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 875./4. - 525.*(eta + 1.)*(eta + 1.)/(2.*2.) - 78750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 183750.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 183750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 66150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 1875.*(xi + 1.)*(xi + 1.)/(2.*2.) + 105000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 47250.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 700.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 245000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 245000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 4375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 110250.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 315.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 110250.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 4375.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 20:
                  return RealGradient(0., 975.*xi/8. - 2925.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 21150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 51825.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 52500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 20475.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 13650.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 12285.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 975./8. - 74025.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 362775.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 183750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 66150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 3525.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 98700.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 44415.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 241850.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 245000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 17275.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. + 217665.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. - 110250.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 4375.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 21:
                  return RealGradient(0., -25.*xi + 600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 8325.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 30825.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 42000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 18900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1260.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 25. + 58275.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 215775.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 147000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 66150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 2775.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 38850.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 34965.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 143850.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 196000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 88200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 10275.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. - 129465.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 88200.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 3500.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                case 22:
                  return RealGradient(0., 225.*xi/4. - 1350.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2025.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 6300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2835.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 225./4. - 23625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 1125.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 15750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 9450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 675.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. + 8505.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2.);
                case 23:
                  return RealGradient(0., -225.*xi/8. + 675.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2700.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2025.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 3150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2835.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 225./8. + 9450.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 14175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 225.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 12600.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5670.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 675.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. - 8505.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2.);
                case 24:
                  return RealGradient(-975.*eta/8. + 2925.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 20475.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 13650.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 12285.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 21150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 51825.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 52500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 975./8. + 3525.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 74025.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 98700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 44415.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 362775.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 183750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 66150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 17275.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 241850.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 217665.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 245000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 88200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 4375.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 110250.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 25:
                  return RealGradient(25.*eta - 600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1260.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 8325.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 30825.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 42000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 18900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 25. - 2775.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 58275.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 38850.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 34965.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. + 215775.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 147000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 66150.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 10275.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 143850.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 129465.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 196000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 88200.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 3500.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                case 26:
                  return RealGradient(-225.*eta/4. + 1350.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 3375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2025.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 225./4. + 1125.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 23625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 15750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 675.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                case 27:
                  return RealGradient(225.*eta/8. - 675.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 4725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 3150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 2700.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2025.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 225./8. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.) - 9450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 12600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5670.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 14175.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. + 675.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 9450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2., 0.);
                case 28:
                  return RealGradient(-195.*eta/8. + 1170.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 12285.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 10920.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 12285.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 8460.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 20730.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 21000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7560.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 195./8. + 705.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 44415.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 78960.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 44415.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 217665.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 110250.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 3455.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 193480.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 217665.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 196000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 70560.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 110250.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 315.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 29:
                  return RealGradient(5.*eta - 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1260.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2240.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1260.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 3330.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 12330.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 16800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 7560.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 5. - 555.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 34965.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 31080.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 34965.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. + 129465.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 88200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 39690.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 2055.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 115080.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 129465.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 156800.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 70560.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 700.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 88200.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 315.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                case 30:
                  return RealGradient(-45.*eta/4. + 540.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2835.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5040.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1350.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 810.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 45./4. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 12600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. - 8505.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 135.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 7560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                case 31:
                  return RealGradient(45.*eta/8. - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2835.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 2520.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 1080.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 810.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 45./8. - 45.*(eta + 1.)*(eta + 1.)/(2.*2.) - 5670.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 10080.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5670.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 8505.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. + 135.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 7560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2., 0.);
                case 32:
                  return RealGradient(0., 195.*xi/8. - 1170.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 8460.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 20730.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 21000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 7560.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 12285.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 10920.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 12285.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 195./8. - 44415.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 217665.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 110250.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 705.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 78960.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 44415.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 193480.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 196000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 70560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 3455.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. + 217665.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. - 110250.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 315.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 33:
                  return RealGradient(0., -5.*xi + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3330.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 12330.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 16800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 7560.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 1260.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2240.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1260.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 5. + 34965.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 129465.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 88200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 555.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 31080.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 34965.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 115080.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 156800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 70560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 2055.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. - 129465.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 88200.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 700.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 315.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                case 34:
                  return RealGradient(0., 45.*xi/4. - 540.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1350.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 810.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 5040.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2835.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 45./4. - 14175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 8505.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 225.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 12600.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 7560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 135.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. + 8505.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2.);
                case 35:
                  return RealGradient(0., -45.*xi/8. + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1080.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 810.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 2520.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2835.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 45./8. + 5670.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8505.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 45.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10080.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5670.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 7560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 135.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. - 8505.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2.);
                case 36:
                  return RealGradient(36.*eta - 288.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 240.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 2376.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 6120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 6300.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 2268.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 36. - 594.*(eta + 1.)*(eta + 1.)/(2.*2.) - 1980.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5100.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5250.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1890.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 1530.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 567.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                case 37:
                  return RealGradient(12.*eta - 192.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 240.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 1584.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4080.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1512.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 12. - 198.*(eta + 1.)*(eta + 1.)/(2.*2.) - 1980.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5100.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5250.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1890.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 510.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 189.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                case 38:
                  return RealGradient(-6.*eta + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 990.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2550.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 945.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 6. + 99.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 990.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2550.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2625.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 945.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 255.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. - 189.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/2., 0.);
                case 39:
                  return RealGradient(0., 36.*xi - 288.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2376.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 6300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 2268.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 240.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 36. - 1980.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5100.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5250.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1890.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 594.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1530.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 567.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 40:
                  return RealGradient(0., 12.*xi - 192.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1584.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4080.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1512.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 240.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 12. - 1980.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5100.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5250.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1890.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 198.*(xi + 1.)*(xi + 1.)/(2.*2.) + 510.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 189.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                case 41:
                  return RealGradient(0., -6.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 990.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2550.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 945.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 6. + 990.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2550.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 945.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 99.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 255.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. - 189.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/2.);
                case 42:
                  return RealGradient(0., -9.*xi/2. + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 864.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5040.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 2268.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 30.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 9./2. + 720.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1890.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 216.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 900.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1260.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 567.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                case 43:
                  return RealGradient(0., -3.*xi/2. + 24.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 576.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1512.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 30.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 3./2. + 720.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1890.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 72.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 300.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 420.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 189.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                case 44:
                  return RealGradient(0., 3.*xi/4. - 15.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 360.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 945.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 15.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 3./4. - 360.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2100.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 945.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 36.*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 210.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 189.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/2.);
                case 45:
                  return RealGradient(-9.*eta/2. + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 30.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 864.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 3600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5040.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 2268.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 9./2. + 216.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 720.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1890.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1260.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 567.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 46:
                  return RealGradient(-3.*eta/2. + 24.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 30.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 576.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1512.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 3./2. + 72.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 720.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1890.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 300.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 420.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 189.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                case 47:
                  return RealGradient(3.*eta/4. - 15.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 15.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 1500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 945.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 3./4. - 36.*(eta + 1.)*(eta + 1.)/(2.*2.) - 360.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2100.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 945.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 150.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 210.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 189.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/2., 0.);
                case 48:
                  return RealGradient(0., 18.*xi - 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 324.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 18. - 270.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 81.*(xi + 1.)*(xi + 1.)/(2.*2.) + 45.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 49:
                  return RealGradient(0., -9.*xi/2. + 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 216.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 30.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 9./2. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 54.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 45.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                case 50:
                  return RealGradient(-18.*eta + 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 324.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18. + 81.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 270.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 45.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 51:
                  return RealGradient(9.*eta/2. - 36.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 30.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 216.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 9./2. - 54.*(eta + 1.)*(eta + 1.)/(2.*2.) - 180.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 45.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                case 52:
                  return RealGradient(-6.*eta + 96.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 216.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 6. + 27.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 270.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 15.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                case 53:
                  return RealGradient(3.*eta/2. - 24.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 30.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 144.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3./2. - 18.*(eta + 1.)*(eta + 1.)/(2.*2.) - 180.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 15.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                case 54:
                  return RealGradient(0., 6.*xi - 96.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 216.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 120.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 6. - 270.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 27.*(xi + 1.)*(xi + 1.)/(2.*2.) + 15.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 55:
                  return RealGradient(0., -3.*xi/2. + 24.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 144.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 120.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 30.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 3./2. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 15.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                case 56:
                  return RealGradient(2.*eta + 2. - 9.*(eta + 1.)*(eta + 1.)/(2.*2.) + 5.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                case 57:
                  return RealGradient(0., 2.*xi + 2. - 9.*(xi + 1.)*(xi + 1.)/(2.*2.) + 5.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                case 58:
                  return RealGradient(0., -xi/2. - 1./2. + 6.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                case 59:
                  return RealGradient(-eta/2. - 1./2. + 6.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 5.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
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
                  return sign * RealGradient(2800.*eta*xi - 375.*eta - 6300.*eta*xi*xi + 5040.*eta*(xi*xi*xi) - 1050.*eta*xi*xi*xi*xi - 300.*xi - 8400.*xi*eta*eta + 10080.*xi*(eta*eta*eta) - 4200.*xi*eta*eta*eta*eta + 25. + 1750.*(eta*eta) + 11340.*(eta*eta)*(xi*xi) - 4200.*eta*eta*xi*xi*xi + 1050.*(xi*xi) - 6300.*xi*xi*eta*eta*eta - 3500.*eta*eta*eta - 1400.*xi*xi*xi + 3150.*(eta*eta*eta*eta) + 630.*(xi*xi*xi*xi) - 1050.*eta*eta*eta*eta*eta, -700.*eta*xi + 4200.*eta*(xi*xi) - 7560.*eta*xi*xi*xi + 4200.*eta*(xi*xi*xi*xi) + 75.*xi + 2100.*xi*(eta*eta) - 2520.*xi*eta*eta*eta + 1050.*xi*(eta*eta*eta*eta) - 7560.*eta*eta*xi*xi + 6300.*(eta*eta)*(xi*xi*xi) - 700.*xi*xi + 4200.*(xi*xi)*(eta*eta*eta) + 2100.*(xi*xi*xi) - 2520.*xi*xi*xi*xi + 1050.*(xi*xi*xi*xi*xi));
                case 1:
                  return sign * RealGradient(-12145.*eta*xi/16. + 695.*eta/128. + 136185.*eta*(xi*xi)/64. - 6825.*eta*xi*xi*xi/4. + 19425.*eta*(xi*xi*xi*xi)/64. + 2865.*xi/32. + 7245.*xi*(eta*eta)/4. - 25515.*xi*eta*eta*eta/16. + 14175.*xi*(eta*eta*eta*eta)/32. - 485./128. + 11235.*(eta*eta)/64. - 240975.*eta*eta*xi*xi/64. + 48825.*(eta*eta)*(xi*xi*xi)/32. - 21105.*xi*xi/64. + 127575.*(xi*xi)*(eta*eta*eta)/64. - 42525.*eta*eta*eta/64. + 13615.*(xi*xi*xi)/32. + 104895.*(eta*eta*eta*eta)/128. - 11655.*xi*xi*xi*xi/64. - 42525.*eta*eta*eta*eta*eta/128., 3045.*eta*xi/32. - 39375.*eta*xi*xi/32. + 84105.*eta*(xi*xi*xi)/32. - 48825.*eta*xi*xi*xi*xi/32. - 2395.*xi/128. + 2835.*xi*(eta*eta)/64. - 14175.*xi*eta*eta*eta/32. + 42525.*xi*(eta*eta*eta*eta)/128. + 25515.*(eta*eta)*(xi*xi)/16. - 127575.*eta*eta*xi*xi*xi/64. + 1715.*(xi*xi)/8. - 14175.*xi*xi*eta*eta*eta/32. - 41755.*xi*xi*xi/64. + 48615.*(xi*xi*xi*xi)/64. - 19425.*xi*xi*xi*xi*xi/64.);
                case 2:
                  return sign * RealGradient(55.*eta/8. - 3255.*eta*xi*xi/4. + 1260.*eta*(xi*xi*xi) - 1575.*eta*xi*xi*xi*xi/4. - 105.*xi/2. + 840.*xi*(eta*eta) - 1575.*xi*eta*eta*eta + 1575.*xi*(eta*eta*eta*eta)/2. + 15./8. - 175.*eta*eta/4. + 945.*(eta*eta)*(xi*xi)/4. - 1575.*eta*eta*xi*xi*xi/2. + 1155.*(xi*xi)/4. + 1575.*(xi*xi)*(eta*eta*eta)/4. + 35.*(eta*eta*eta)/4. - 945.*xi*xi*xi/2. + 735.*(eta*eta*eta*eta)/8. + 945.*(xi*xi*xi*xi)/4. - 525.*eta*eta*eta*eta*eta/8., 35.*eta*xi/2. + 105.*eta*(xi*xi)/2. - 1575.*eta*xi*xi*xi/2. + 1575.*eta*(xi*xi*xi*xi)/2. + 85.*xi/8. - 525.*xi*eta*eta/4. + 105.*xi*(eta*eta*eta)/2. + 525.*xi*(eta*eta*eta*eta)/8. + 945.*(eta*eta)*(xi*xi) - 1575.*eta*eta*xi*xi*xi/4. - 315.*xi*xi/2. - 1575.*xi*xi*eta*eta*eta/2. + 2485.*(xi*xi*xi)/4. - 3465.*xi*xi*xi*xi/4. + 1575.*(xi*xi*xi*xi*xi)/4.);
                case 3:
                  return sign * RealGradient(1855.*eta*xi/16. - 825.*eta/128. - 16695.*eta*xi*xi/64. - 945.*eta*xi*xi*xi/4. + 19425.*eta*(xi*xi*xi*xi)/64. + 705.*xi/32. - 1155.*xi*eta*eta/4. - 315.*xi*eta*eta*eta/16. + 5775.*xi*(eta*eta*eta*eta)/32. - 85./128. + 595.*(eta*eta)/64. + 76545.*(eta*eta)*(xi*xi)/64. - 9975.*eta*eta*xi*xi*xi/32. - 9345.*xi*xi/64. - 48825.*xi*xi*eta*eta*eta/64. + 595.*(eta*eta*eta)/64. + 9695.*(xi*xi*xi)/32. - 945.*eta*eta*eta*eta/128. - 11655.*xi*xi*xi*xi/64. - 525.*eta*eta*eta*eta*eta/128., -875.*eta*xi/32. + 9345.*eta*(xi*xi)/32. - 21735.*eta*xi*xi*xi/32. + 9975.*eta*(xi*xi*xi*xi)/32. - 555.*xi/128. + 1155.*xi*(eta*eta)/64. + 945.*xi*(eta*eta*eta)/32. + 525.*xi*(eta*eta*eta*eta)/128. - 4725.*eta*eta*xi*xi/16. + 48825.*(eta*eta)*(xi*xi*xi)/64. + 595.*(xi*xi)/8. - 5775.*xi*xi*eta*eta*eta/32. - 22155.*xi*xi*xi/64. + 36855.*(xi*xi*xi*xi)/64. - 19425.*xi*xi*xi*xi*xi/64.);
                case 4:
                  return sign * RealGradient(140.*eta*xi - 5.*eta - 840.*eta*xi*xi + 1680.*eta*(xi*xi*xi) - 1050.*eta*xi*xi*xi*xi - 120.*xi + 5. + 630.*(xi*xi) - 1120.*xi*xi*xi + 630.*(xi*xi*xi*xi), 25.*xi - 350.*xi*xi + 1400.*(xi*xi*xi) - 2100.*xi*xi*xi*xi + 1050.*(xi*xi*xi*xi*xi));
                case 5:
                  return sign * RealGradient(140.*eta*xi - 5.*eta - 840.*eta*xi*xi + 1680.*eta*(xi*xi*xi) - 1050.*eta*xi*xi*xi*xi, 25.*xi - 350.*xi*xi + 1400.*(xi*xi*xi) - 2100.*xi*xi*xi*xi + 1050.*(xi*xi*xi*xi*xi));
                case 6:
                  return sign * RealGradient(735.*eta*xi/4. - 10.*eta - 2835.*eta*xi*xi/4. + 14175.*eta*(xi*xi*xi)/16. - 42525.*eta*xi*xi*xi*xi/128. - 1575.*xi*eta*eta/2. + 12285.*xi*(eta*eta*eta)/16. - 1575.*xi*eta*eta*eta*eta/8. + 385.*(eta*eta)/8. + 76545.*(eta*eta)*(xi*xi)/32. - 14175.*eta*eta*xi*xi*xi/8. - 42525.*xi*xi*eta*eta*eta/32. - 245.*eta*eta*eta/4. + 1785.*(eta*eta*eta*eta)/64. - 525.*eta*eta*eta*eta*eta/128., -455.*eta*xi/4. + 2205.*eta*(xi*xi)/2. - 42525.*eta*xi*xi*xi/16. + 14175.*eta*(xi*xi*xi*xi)/8. + 20.*xi + 525.*xi*(eta*eta)/4. - 735.*xi*eta*eta*eta/16. + 525.*xi*(eta*eta*eta*eta)/128. - 31185.*eta*eta*xi*xi/32. + 42525.*(eta*eta)*(xi*xi*xi)/32. - 1785.*xi*xi/8. + 1575.*(xi*xi)*(eta*eta*eta)/8. + 2835.*(xi*xi*xi)/4. - 53865.*xi*xi*xi*xi/64. + 42525.*(xi*xi*xi*xi*xi)/128.);
                case 7:
                  return sign * RealGradient(175.*eta*xi - 15.*eta - 420.*eta*xi*xi + 315.*eta*(xi*xi*xi) - 525.*eta*xi*xi*xi*xi/8. - 1260.*xi*eta*eta + 2205.*xi*(eta*eta*eta) - 1050.*xi*eta*eta*eta*eta + 245.*(eta*eta)/2. + 4725.*(eta*eta)*(xi*xi)/2. - 1050.*eta*eta*xi*xi*xi - 4725.*xi*xi*eta*eta*eta/2. - 280.*eta*eta*eta + 945.*(eta*eta*eta*eta)/4. - 525.*eta*eta*eta*eta*eta/8., -175.*eta*xi + 1260.*eta*(xi*xi) - 2205.*eta*xi*xi*xi + 1050.*eta*(xi*xi*xi*xi) + 15.*xi + 420.*xi*(eta*eta) - 315.*xi*eta*eta*eta + 525.*xi*(eta*eta*eta*eta)/8. - 4725.*eta*eta*xi*xi/2. + 4725.*(eta*eta)*(xi*xi*xi)/2. - 245.*xi*xi/2. + 1050.*(xi*xi)*(eta*eta*eta) + 280.*(xi*xi*xi) - 945.*xi*xi*xi*xi/4. + 525.*(xi*xi*xi*xi*xi)/8.);
                case 8:
                  return sign * RealGradient(455.*eta*xi/4. - 20.*eta - 525.*eta*xi*xi/4. + 735.*eta*(xi*xi*xi)/16. - 525.*eta*xi*xi*xi*xi/128. - 2205.*xi*eta*eta/2. + 42525.*xi*(eta*eta*eta)/16. - 14175.*xi*eta*eta*eta*eta/8. + 1785.*(eta*eta)/8. + 31185.*(eta*eta)*(xi*xi)/32. - 1575.*eta*eta*xi*xi*xi/8. - 42525.*xi*xi*eta*eta*eta/32. - 2835.*eta*eta*eta/4. + 53865.*(eta*eta*eta*eta)/64. - 42525.*eta*eta*eta*eta*eta/128., -735.*eta*xi/4. + 1575.*eta*(xi*xi)/2. - 12285.*eta*xi*xi*xi/16. + 1575.*eta*(xi*xi*xi*xi)/8. + 10.*xi + 2835.*xi*(eta*eta)/4. - 14175.*xi*eta*eta*eta/16. + 42525.*xi*(eta*eta*eta*eta)/128. - 76545.*eta*eta*xi*xi/32. + 42525.*(eta*eta)*(xi*xi*xi)/32. - 385.*xi*xi/8. + 14175.*(xi*xi)*(eta*eta*eta)/8. + 245.*(xi*xi*xi)/4. - 1785.*xi*xi*xi*xi/64. + 525.*(xi*xi*xi*xi*xi)/128.);
                case 9:
                  return sign * RealGradient(-25.*eta + 350.*(eta*eta) - 1400.*eta*eta*eta + 2100.*(eta*eta*eta*eta) - 1050.*eta*eta*eta*eta*eta, -140.*eta*xi + 5.*xi + 840.*xi*(eta*eta) - 1680.*xi*eta*eta*eta + 1050.*xi*(eta*eta*eta*eta));
                case 10:
                  return sign * RealGradient(-25.*eta + 350.*(eta*eta) - 1400.*eta*eta*eta + 2100.*(eta*eta*eta*eta) - 1050.*eta*eta*eta*eta*eta, -140.*eta*xi + 120.*eta + 5.*xi + 840.*xi*(eta*eta) - 1680.*xi*eta*eta*eta + 1050.*xi*(eta*eta*eta*eta) - 5. - 630.*eta*eta + 1120.*(eta*eta*eta) - 630.*eta*eta*eta*eta);
                case 11:
                  return sign * RealGradient(875.*eta*xi/32. + 555.*eta/128. - 1155.*eta*xi*xi/64. - 945.*eta*xi*xi*xi/32. - 525.*eta*xi*xi*xi*xi/128. - 9345.*xi*eta*eta/32. + 21735.*xi*(eta*eta*eta)/32. - 9975.*xi*eta*eta*eta*eta/32. - 595.*eta*eta/8. + 4725.*(eta*eta)*(xi*xi)/16. + 5775.*(eta*eta)*(xi*xi*xi)/32. - 48825.*xi*xi*eta*eta*eta/64. + 22155.*(eta*eta*eta)/64. - 36855.*eta*eta*eta*eta/64. + 19425.*(eta*eta*eta*eta*eta)/64., -1855.*eta*xi/16. - 705.*eta/32. + 1155.*eta*(xi*xi)/4. + 315.*eta*(xi*xi*xi)/16. - 5775.*eta*xi*xi*xi*xi/32. + 825.*xi/128. + 16695.*xi*(eta*eta)/64. + 945.*xi*(eta*eta*eta)/4. - 19425.*xi*eta*eta*eta*eta/64. + 85./128. + 9345.*(eta*eta)/64. - 76545.*eta*eta*xi*xi/64. + 48825.*(eta*eta)*(xi*xi*xi)/64. - 595.*xi*xi/64. + 9975.*(xi*xi)*(eta*eta*eta)/32. - 9695.*eta*eta*eta/32. - 595.*xi*xi*xi/64. + 11655.*(eta*eta*eta*eta)/64. + 945.*(xi*xi*xi*xi)/128. + 525.*(xi*xi*xi*xi*xi)/128.);
                case 12:
                  return sign * RealGradient(-35.*eta*xi/2. - 85.*eta/8. + 525.*eta*(xi*xi)/4. - 105.*eta*xi*xi*xi/2. - 525.*eta*xi*xi*xi*xi/8. - 105.*xi*eta*eta/2. + 1575.*xi*(eta*eta*eta)/2. - 1575.*xi*eta*eta*eta*eta/2. + 315.*(eta*eta)/2. - 945.*eta*eta*xi*xi + 1575.*(eta*eta)*(xi*xi*xi)/2. + 1575.*(xi*xi)*(eta*eta*eta)/4. - 2485.*eta*eta*eta/4. + 3465.*(eta*eta*eta*eta)/4. - 1575.*eta*eta*eta*eta*eta/4., 105.*eta/2. - 840.*eta*xi*xi + 1575.*eta*(xi*xi*xi) - 1575.*eta*xi*xi*xi*xi/2. - 55.*xi/8. + 3255.*xi*(eta*eta)/4. - 1260.*xi*eta*eta*eta + 1575.*xi*(eta*eta*eta*eta)/4. - 15./8. - 1155.*eta*eta/4. - 945.*eta*eta*xi*xi/4. - 1575.*eta*eta*xi*xi*xi/4. + 175.*(xi*xi)/4. + 1575.*(xi*xi)*(eta*eta*eta)/2. + 945.*(eta*eta*eta)/2. - 35.*xi*xi*xi/4. - 945.*eta*eta*eta*eta/4. - 735.*xi*xi*xi*xi/8. + 525.*(xi*xi*xi*xi*xi)/8.);
                case 13:
                  return sign * RealGradient(-3045.*eta*xi/32. + 2395.*eta/128. - 2835.*eta*xi*xi/64. + 14175.*eta*(xi*xi*xi)/32. - 42525.*eta*xi*xi*xi*xi/128. + 39375.*xi*(eta*eta)/32. - 84105.*xi*eta*eta*eta/32. + 48825.*xi*(eta*eta*eta*eta)/32. - 1715.*eta*eta/8. - 25515.*eta*eta*xi*xi/16. + 14175.*(eta*eta)*(xi*xi*xi)/32. + 127575.*(xi*xi)*(eta*eta*eta)/64. + 41755.*(eta*eta*eta)/64. - 48615.*eta*eta*eta*eta/64. + 19425.*(eta*eta*eta*eta*eta)/64., 12145.*eta*xi/16. - 2865.*eta/32. - 7245.*eta*xi*xi/4. + 25515.*eta*(xi*xi*xi)/16. - 14175.*eta*xi*xi*xi*xi/32. - 695.*xi/128. - 136185.*xi*eta*eta/64. + 6825.*xi*(eta*eta*eta)/4. - 19425.*xi*eta*eta*eta*eta/64. + 485./128. + 21105.*(eta*eta)/64. + 240975.*(eta*eta)*(xi*xi)/64. - 127575.*eta*eta*xi*xi*xi/64. - 11235.*xi*xi/64. - 48825.*xi*xi*eta*eta*eta/32. - 13615.*eta*eta*eta/32. + 42525.*(xi*xi*xi)/64. + 11655.*(eta*eta*eta*eta)/64. - 104895.*xi*xi*xi*xi/128. + 42525.*(xi*xi*xi*xi*xi)/128.);
                case 14:
                  return sign * RealGradient(700.*eta*xi - 75.*eta - 2100.*eta*xi*xi + 2520.*eta*(xi*xi*xi) - 1050.*eta*xi*xi*xi*xi - 4200.*xi*eta*eta + 7560.*xi*(eta*eta*eta) - 4200.*xi*eta*eta*eta*eta + 700.*(eta*eta) + 7560.*(eta*eta)*(xi*xi) - 4200.*eta*eta*xi*xi*xi - 6300.*xi*xi*eta*eta*eta - 2100.*eta*eta*eta + 2520.*(eta*eta*eta*eta) - 1050.*eta*eta*eta*eta*eta, -2800.*eta*xi + 300.*eta + 8400.*eta*(xi*xi) - 10080.*eta*xi*xi*xi + 4200.*eta*(xi*xi*xi*xi) + 375.*xi + 6300.*xi*(eta*eta) - 5040.*xi*eta*eta*eta + 1050.*xi*(eta*eta*eta*eta) - 25. - 1050.*eta*eta - 11340.*eta*eta*xi*xi + 6300.*(eta*eta)*(xi*xi*xi) - 1750.*xi*xi + 4200.*(xi*xi)*(eta*eta*eta) + 1400.*(eta*eta*eta) + 3500.*(xi*xi*xi) - 630.*eta*eta*eta*eta - 3150.*xi*xi*xi*xi + 1050.*(xi*xi*xi*xi*xi));
                case 15:
                  return RealGradient(-19600.*eta*xi + 3500.*eta + 35280.*eta*(xi*xi) - 23520.*eta*xi*xi*xi + 4200.*eta*(xi*xi*xi*xi) + 94080.*xi*(eta*eta) - 141120.*xi*eta*eta*eta + 67200.*xi*(eta*eta*eta*eta) - 24500.*eta*eta - 105840.*eta*eta*xi*xi + 33600.*(eta*eta)*(xi*xi*xi) + 75600.*(xi*xi)*(eta*eta*eta) + 58800.*(eta*eta*eta) - 58800.*eta*eta*eta*eta + 21000.*(eta*eta*eta*eta*eta), 9800.*eta*xi - 47040.*eta*xi*xi + 70560.*eta*(xi*xi*xi) - 33600.*eta*xi*xi*xi*xi - 700.*xi - 35280.*xi*eta*eta + 47040.*xi*(eta*eta*eta) - 21000.*xi*eta*eta*eta*eta + 105840.*(eta*eta)*(xi*xi) - 75600.*eta*eta*xi*xi*xi + 4900.*(xi*xi) - 67200.*xi*xi*eta*eta*eta - 11760.*xi*xi*xi + 11760.*(xi*xi*xi*xi) - 4200.*xi*xi*xi*xi*xi);
                case 16:
                  return RealGradient(9800.*eta*xi - 700.*eta - 35280.*eta*xi*xi + 47040.*eta*(xi*xi*xi) - 21000.*eta*xi*xi*xi*xi - 47040.*xi*eta*eta + 70560.*xi*(eta*eta*eta) - 33600.*xi*eta*eta*eta*eta + 4900.*(eta*eta) + 105840.*(eta*eta)*(xi*xi) - 67200.*eta*eta*xi*xi*xi - 75600.*xi*xi*eta*eta*eta - 11760.*eta*eta*eta + 11760.*(eta*eta*eta*eta) - 4200.*eta*eta*eta*eta*eta, -19600.*eta*xi + 94080.*eta*(xi*xi) - 141120.*eta*xi*xi*xi + 67200.*eta*(xi*xi*xi*xi) + 3500.*xi + 35280.*xi*(eta*eta) - 23520.*xi*eta*eta*eta + 4200.*xi*(eta*eta*eta*eta) - 105840.*eta*eta*xi*xi + 75600.*(eta*eta)*(xi*xi*xi) - 24500.*xi*xi + 33600.*(xi*xi)*(eta*eta*eta) + 58800.*(xi*xi*xi) - 58800.*xi*xi*xi*xi + 21000.*(xi*xi*xi*xi*xi));
                case 17:
                  return RealGradient(6440.*eta*xi - 280.*eta - 30240.*eta*xi*xi + 43680.*eta*(xi*xi*xi) - 16800.*eta*xi*xi*xi*xi - 6720.*xi*eta*eta + 280.*(eta*eta) + 30240.*(eta*eta)*(xi*xi) - 33600.*eta*eta*xi*xi*xi, -1120.*eta*xi + 13440.*eta*(xi*xi) - 40320.*eta*xi*xi*xi + 33600.*eta*(xi*xi*xi*xi) + 560.*xi - 7280.*xi*xi + 26880.*(xi*xi*xi) - 36960.*xi*xi*xi*xi + 16800.*(xi*xi*xi*xi*xi));
                case 18:
                  return RealGradient(-3640.*eta*xi + 140.*eta + 20160.*eta*(xi*xi) - 36960.*eta*xi*xi*xi + 21000.*eta*(xi*xi*xi*xi) + 3360.*xi*(eta*eta) - 140.*eta*eta - 15120.*eta*eta*xi*xi + 16800.*(eta*eta)*(xi*xi*xi), 560.*eta*xi - 6720.*eta*xi*xi + 20160.*eta*(xi*xi*xi) - 16800.*eta*xi*xi*xi*xi - 700.*xi + 9100.*(xi*xi) - 33600.*xi*xi*xi + 46200.*(xi*xi*xi*xi) - 21000.*xi*xi*xi*xi*xi);
                case 19:
                  return RealGradient(560.*eta*xi - 700.*eta - 6720.*xi*eta*eta + 20160.*xi*(eta*eta*eta) - 16800.*xi*eta*eta*eta*eta + 9100.*(eta*eta) - 33600.*eta*eta*eta + 46200.*(eta*eta*eta*eta) - 21000.*eta*eta*eta*eta*eta, -3640.*eta*xi + 3360.*eta*(xi*xi) + 140.*xi + 20160.*xi*(eta*eta) - 36960.*xi*eta*eta*eta + 21000.*xi*(eta*eta*eta*eta) - 15120.*eta*eta*xi*xi - 140.*xi*xi + 16800.*(xi*xi)*(eta*eta*eta));
                case 20:
                  return RealGradient(-1120.*eta*xi + 560.*eta + 13440.*xi*(eta*eta) - 40320.*xi*eta*eta*eta + 33600.*xi*(eta*eta*eta*eta) - 7280.*eta*eta + 26880.*(eta*eta*eta) - 36960.*eta*eta*eta*eta + 16800.*(eta*eta*eta*eta*eta), 6440.*eta*xi - 6720.*eta*xi*xi - 280.*xi - 30240.*xi*eta*eta + 43680.*xi*(eta*eta*eta) - 16800.*xi*eta*eta*eta*eta + 30240.*(eta*eta)*(xi*xi) + 280.*(xi*xi) - 33600.*xi*xi*eta*eta*eta);
                case 21:
                  return RealGradient(17920.*eta*xi/3. - 420.*eta - 17920.*eta*xi*xi + 156800.*eta*(xi*xi*xi)/9. - 44800.*eta*xi*xi*xi*xi/9. - 29120.*xi*eta*eta + 116480.*xi*(eta*eta*eta)/3. - 140000.*xi*eta*eta*eta*eta/9. + 7420.*(eta*eta)/3. + 62720.*(eta*eta)*(xi*xi) - 291200.*eta*eta*xi*xi*xi/9. - 44800.*xi*xi*eta*eta*eta - 4480.*eta*eta*eta + 28840.*(eta*eta*eta*eta)/9. - 7000.*eta*eta*eta*eta*eta/9., -10360.*eta*xi/3. + 28000.*eta*(xi*xi) - 170240.*eta*xi*xi*xi/3. + 291200.*eta*(xi*xi*xi*xi)/9. + 420.*xi + 6720.*xi*(eta*eta) - 39200.*xi*eta*eta*eta/9. + 7000.*xi*(eta*eta*eta*eta)/9. - 40880.*eta*eta*xi*xi + 44800.*(eta*eta)*(xi*xi*xi) - 12460.*xi*xi/3. + 140000.*(xi*xi)*(eta*eta*eta)/9. + 35840.*(xi*xi*xi)/3. - 118720.*xi*xi*xi*xi/9. + 44800.*(xi*xi*xi*xi*xi)/9.);
                case 22:
                  return RealGradient(-12880.*eta*xi/3. + 280.*eta + 14560.*eta*(xi*xi) - 152320.*eta*xi*xi*xi/9. + 56000.*eta*(xi*xi*xi*xi)/9. + 22400.*xi*(eta*eta) - 91840.*xi*eta*eta*eta/3. + 112000.*xi*(eta*eta*eta*eta)/9. - 5320.*eta*eta/3. - 54880.*eta*eta*xi*xi + 313600.*(eta*eta)*(xi*xi*xi)/9. + 39200.*(xi*xi)*(eta*eta*eta) + 3360.*(eta*eta*eta) - 22400.*eta*eta*eta*eta/9. + 5600.*(eta*eta*eta*eta*eta)/9., 12040.*eta*xi/3. - 31360.*eta*xi*xi + 185920.*eta*(xi*xi*xi)/3. - 313600.*eta*xi*xi*xi*xi/9. - 560.*xi - 6720.*xi*eta*eta + 34720.*xi*(eta*eta*eta)/9. - 5600.*xi*eta*eta*eta*eta/9. + 38080.*(eta*eta)*(xi*xi) - 39200.*eta*eta*xi*xi*xi + 16240.*(xi*xi)/3. - 112000.*xi*xi*eta*eta*eta/9. - 45920.*xi*xi*xi/3. + 150080.*(xi*xi*xi*xi)/9. - 56000.*xi*xi*xi*xi*xi/9.);
                case 23:
                  return RealGradient(12040.*eta*xi/3. - 560.*eta - 6720.*eta*xi*xi + 34720.*eta*(xi*xi*xi)/9. - 5600.*eta*xi*xi*xi*xi/9. - 31360.*xi*eta*eta + 185920.*xi*(eta*eta*eta)/3. - 313600.*xi*eta*eta*eta*eta/9. + 16240.*(eta*eta)/3. + 38080.*(eta*eta)*(xi*xi) - 112000.*eta*eta*xi*xi*xi/9. - 39200.*xi*xi*eta*eta*eta - 45920.*eta*eta*eta/3. + 150080.*(eta*eta*eta*eta)/9. - 56000.*eta*eta*eta*eta*eta/9., -12880.*eta*xi/3. + 22400.*eta*(xi*xi) - 91840.*eta*xi*xi*xi/3. + 112000.*eta*(xi*xi*xi*xi)/9. + 280.*xi + 14560.*xi*(eta*eta) - 152320.*xi*eta*eta*eta/9. + 56000.*xi*(eta*eta*eta*eta)/9. - 54880.*eta*eta*xi*xi + 39200.*(eta*eta)*(xi*xi*xi) - 5320.*xi*xi/3. + 313600.*(xi*xi)*(eta*eta*eta)/9. + 3360.*(xi*xi*xi) - 22400.*xi*xi*xi*xi/9. + 5600.*(xi*xi*xi*xi*xi)/9.);
                case 24:
                  return RealGradient(-10360.*eta*xi/3. + 420.*eta + 6720.*eta*(xi*xi) - 39200.*eta*xi*xi*xi/9. + 7000.*eta*(xi*xi*xi*xi)/9. + 28000.*xi*(eta*eta) - 170240.*xi*eta*eta*eta/3. + 291200.*xi*(eta*eta*eta*eta)/9. - 12460.*eta*eta/3. - 40880.*eta*eta*xi*xi + 140000.*(eta*eta)*(xi*xi*xi)/9. + 44800.*(xi*xi)*(eta*eta*eta) + 35840.*(eta*eta*eta)/3. - 118720.*eta*eta*eta*eta/9. + 44800.*(eta*eta*eta*eta*eta)/9., 17920.*eta*xi/3. - 29120.*eta*xi*xi + 116480.*eta*(xi*xi*xi)/3. - 140000.*eta*xi*xi*xi*xi/9. - 420.*xi - 17920.*xi*eta*eta + 156800.*xi*(eta*eta*eta)/9. - 44800.*xi*eta*eta*eta*eta/9. + 62720.*(eta*eta)*(xi*xi) - 44800.*eta*eta*xi*xi*xi + 7420.*(xi*xi)/3. - 291200.*xi*xi*eta*eta*eta/9. - 4480.*xi*xi*xi + 28840.*(xi*xi*xi*xi)/9. - 7000.*xi*xi*xi*xi*xi/9.);
                case 25:
                  return RealGradient(-12880.*eta*xi/9. - 700.*eta/9. + 5600.*eta*(xi*xi) - 49280.*eta*xi*xi*xi/9. + 11200.*eta*(xi*xi*xi*xi)/9. - 15680.*xi*eta*eta/3. + 51520.*xi*(eta*eta*eta)/3. - 95200.*xi*eta*eta*eta*eta/9. + 25900.*(eta*eta)/9. - 1120.*eta*eta*xi*xi + 22400.*(eta*eta)*(xi*xi*xi)/9. - 5600.*xi*xi*eta*eta*eta - 28000.*eta*eta*eta/3. + 93800.*(eta*eta*eta*eta)/9. - 35000.*eta*eta*eta*eta*eta/9., -10360.*eta*xi/9. + 7840.*eta*(xi*xi)/3. + 2240.*eta*(xi*xi*xi)/3. - 22400.*eta*xi*xi*xi*xi/9. + 140.*xi/9. + 5600.*xi*(eta*eta) - 75040.*xi*eta*eta*eta/9. + 35000.*xi*(eta*eta*eta*eta)/9. - 12880.*eta*eta*xi*xi + 5600.*(eta*eta)*(xi*xi*xi) + 3220.*(xi*xi)/9. + 95200.*(xi*xi)*(eta*eta*eta)/9. - 5600.*xi*xi*xi/3. + 24640.*(xi*xi*xi*xi)/9. - 11200.*xi*xi*xi*xi*xi/9.);
                case 26:
                  return RealGradient(-9520.*eta*xi/9. + 1400.*eta/9. - 1120.*eta*xi*xi + 71680.*eta*(xi*xi*xi)/9. - 56000.*eta*xi*xi*xi*xi/9. + 44800.*xi*(eta*eta)/3. - 82880.*xi*eta*eta*eta/3. + 123200.*xi*(eta*eta*eta*eta)/9. - 14840.*eta*eta/9. - 25760.*eta*eta*xi*xi + 89600.*(eta*eta)*(xi*xi*xi)/9. + 28000.*(xi*xi)*(eta*eta*eta) + 12320.*(eta*eta*eta)/3. - 34720.*eta*eta*eta*eta/9. + 11200.*(eta*eta*eta*eta*eta)/9., 51800.*eta*xi/9. - 69440.*eta*xi*xi/3. + 82880.*eta*(xi*xi*xi)/3. - 89600.*eta*xi*xi*xi*xi/9. + 560.*xi/9. - 13440.*xi*eta*eta + 79520.*xi*(eta*eta*eta)/9. - 11200.*xi*eta*eta*eta*eta/9. + 40320.*(eta*eta)*(xi*xi) - 28000.*eta*eta*xi*xi*xi - 24080.*xi*xi/9. - 123200.*xi*xi*eta*eta*eta/9. + 32480.*(xi*xi*xi)/3. - 129920.*xi*xi*xi*xi/9. + 56000.*(xi*xi*xi*xi*xi)/9.);
                case 27:
                  return RealGradient(2800.*eta*xi/9. + 700.*eta/9. - 560.*eta*xi*xi - 1120.*eta*xi*xi*xi/9. + 1400.*eta*(xi*xi*xi*xi)/9. - 4480.*xi*eta*eta/3. - 11200.*xi*eta*eta*eta/3. + 44800.*xi*(eta*eta*eta*eta)/9. - 13300.*eta*eta/9. + 6160.*(eta*eta)*(xi*xi) - 22400.*eta*eta*xi*xi*xi/9. - 2800.*xi*xi*eta*eta*eta + 19600.*(eta*eta*eta)/3. - 81200.*eta*eta*eta*eta/9. + 35000.*(eta*eta*eta*eta*eta)/9., 5320.*eta*xi/9. + 2240.*eta*(xi*xi)/3. - 12320.*eta*xi*xi*xi/3. + 22400.*eta*(xi*xi*xi*xi)/9. - 140.*xi/9. - 3920.*xi*eta*eta + 64960.*xi*(eta*eta*eta)/9. - 35000.*xi*eta*eta*eta*eta/9. + 2800.*(eta*eta)*(xi*xi) + 2800.*(eta*eta)*(xi*xi*xi) - 700.*xi*xi/9. - 44800.*xi*xi*eta*eta*eta/9. + 560.*(xi*xi*xi)/3. + 560.*(xi*xi*xi*xi)/9. - 1400.*xi*xi*xi*xi*xi/9.);
                case 28:
                  return RealGradient(280.*eta*xi/9. - 980.*eta/9. + 1680.*eta*(xi*xi) - 11200.*eta*xi*xi*xi/9. - 7000.*eta*xi*xi*xi*xi/9. - 11200.*xi*eta*eta/3. + 52640.*xi*(eta*eta*eta)/3. - 123200.*xi*eta*eta*eta*eta/9. + 13580.*(eta*eta)/9. - 9520.*eta*eta*xi*xi + 112000.*(eta*eta)*(xi*xi*xi)/9. - 2800.*xi*xi*eta*eta*eta - 16240.*eta*eta*eta/3. + 59920.*(eta*eta*eta*eta)/9. - 23800.*eta*eta*eta*eta*eta/9., -8960.*eta*xi/9. - 17920.*eta*xi*xi/3. + 58240.*eta*(xi*xi*xi)/3. - 112000.*eta*xi*xi*xi*xi/9. - 140.*xi/9. + 9520.*xi*(eta*eta) - 99680.*xi*eta*eta*eta/9. + 23800.*xi*(eta*eta*eta*eta)/9. - 15120.*eta*eta*xi*xi + 2800.*(eta*eta)*(xi*xi*xi) + 4340.*(xi*xi)/9. + 123200.*(xi*xi)*(eta*eta*eta)/9. - 560.*xi*xi*xi - 6160.*xi*xi*xi*xi/9. + 7000.*(xi*xi*xi*xi*xi)/9.);
                case 29:
                  return RealGradient(51800.*eta*xi/9. + 560.*eta/9. - 13440.*eta*xi*xi + 79520.*eta*(xi*xi*xi)/9. - 11200.*eta*xi*xi*xi*xi/9. - 69440.*xi*eta*eta/3. + 82880.*xi*(eta*eta*eta)/3. - 89600.*xi*eta*eta*eta*eta/9. - 24080.*eta*eta/9. + 40320.*(eta*eta)*(xi*xi) - 123200.*eta*eta*xi*xi*xi/9. - 28000.*xi*xi*eta*eta*eta + 32480.*(eta*eta*eta)/3. - 129920.*eta*eta*eta*eta/9. + 56000.*(eta*eta*eta*eta*eta)/9., -9520.*eta*xi/9. + 44800.*eta*(xi*xi)/3. - 82880.*eta*xi*xi*xi/3. + 123200.*eta*(xi*xi*xi*xi)/9. + 1400.*xi/9. - 1120.*xi*eta*eta + 71680.*xi*(eta*eta*eta)/9. - 56000.*xi*eta*eta*eta*eta/9. - 25760.*eta*eta*xi*xi + 28000.*(eta*eta)*(xi*xi*xi) - 14840.*xi*xi/9. + 89600.*(xi*xi)*(eta*eta*eta)/9. + 12320.*(xi*xi*xi)/3. - 34720.*xi*xi*xi*xi/9. + 11200.*(xi*xi*xi*xi*xi)/9.);
                case 30:
                  return RealGradient(-10360.*eta*xi/9. + 140.*eta/9. + 5600.*eta*(xi*xi) - 75040.*eta*xi*xi*xi/9. + 35000.*eta*(xi*xi*xi*xi)/9. + 7840.*xi*(eta*eta)/3. + 2240.*xi*(eta*eta*eta)/3. - 22400.*xi*eta*eta*eta*eta/9. + 3220.*(eta*eta)/9. - 12880.*eta*eta*xi*xi + 95200.*(eta*eta)*(xi*xi*xi)/9. + 5600.*(xi*xi)*(eta*eta*eta) - 5600.*eta*eta*eta/3. + 24640.*(eta*eta*eta*eta)/9. - 11200.*eta*eta*eta*eta*eta/9., -12880.*eta*xi/9. - 15680.*eta*xi*xi/3. + 51520.*eta*(xi*xi*xi)/3. - 95200.*eta*xi*xi*xi*xi/9. - 700.*xi/9. + 5600.*xi*(eta*eta) - 49280.*xi*eta*eta*eta/9. + 11200.*xi*(eta*eta*eta*eta)/9. - 1120.*eta*eta*xi*xi - 5600.*eta*eta*xi*xi*xi + 25900.*(xi*xi)/9. + 22400.*(xi*xi)*(eta*eta*eta)/9. - 28000.*xi*xi*xi/3. + 93800.*(xi*xi*xi*xi)/9. - 35000.*xi*xi*xi*xi*xi/9.);
                case 31:
                  return RealGradient(-8960.*eta*xi/9. - 140.*eta/9. + 9520.*eta*(xi*xi) - 99680.*eta*xi*xi*xi/9. + 23800.*eta*(xi*xi*xi*xi)/9. - 17920.*xi*eta*eta/3. + 58240.*xi*(eta*eta*eta)/3. - 112000.*xi*eta*eta*eta*eta/9. + 4340.*(eta*eta)/9. - 15120.*eta*eta*xi*xi + 123200.*(eta*eta)*(xi*xi*xi)/9. + 2800.*(xi*xi)*(eta*eta*eta) - 560.*eta*eta*eta - 6160.*eta*eta*eta*eta/9. + 7000.*(eta*eta*eta*eta*eta)/9., 280.*eta*xi/9. - 11200.*eta*xi*xi/3. + 52640.*eta*(xi*xi*xi)/3. - 123200.*eta*xi*xi*xi*xi/9. - 980.*xi/9. + 1680.*xi*(eta*eta) - 11200.*xi*eta*eta*eta/9. - 7000.*xi*eta*eta*eta*eta/9. - 9520.*eta*eta*xi*xi - 2800.*eta*eta*xi*xi*xi + 13580.*(xi*xi)/9. + 112000.*(xi*xi)*(eta*eta*eta)/9. - 16240.*xi*xi*xi/3. + 59920.*(xi*xi*xi*xi)/9. - 23800.*xi*xi*xi*xi*xi/9.);
                case 32:
                  return RealGradient(5320.*eta*xi/9. - 140.*eta/9. - 3920.*eta*xi*xi + 64960.*eta*(xi*xi*xi)/9. - 35000.*eta*xi*xi*xi*xi/9. + 2240.*xi*(eta*eta)/3. - 12320.*xi*eta*eta*eta/3. + 22400.*xi*(eta*eta*eta*eta)/9. - 700.*eta*eta/9. + 2800.*(eta*eta)*(xi*xi) - 44800.*eta*eta*xi*xi*xi/9. + 2800.*(xi*xi)*(eta*eta*eta) + 560.*(eta*eta*eta)/3. + 560.*(eta*eta*eta*eta)/9. - 1400.*eta*eta*eta*eta*eta/9., 2800.*eta*xi/9. - 4480.*eta*xi*xi/3. - 11200.*eta*xi*xi*xi/3. + 44800.*eta*(xi*xi*xi*xi)/9. + 700.*xi/9. - 560.*xi*eta*eta - 1120.*xi*eta*eta*eta/9. + 1400.*xi*(eta*eta*eta*eta)/9. + 6160.*(eta*eta)*(xi*xi) - 2800.*eta*eta*xi*xi*xi - 13300.*xi*xi/9. - 22400.*xi*xi*eta*eta*eta/9. + 19600.*(xi*xi*xi)/3. - 81200.*xi*xi*xi*xi/9. + 35000.*(xi*xi*xi*xi*xi)/9.);
                case 33:
                  return RealGradient(-6440.*eta*xi/9. + 280.*eta/9. + 1680.*eta*(xi*xi) - 4480.*eta*xi*xi*xi/3. + 1400.*eta*(xi*xi*xi*xi)/3. + 19040.*xi*(eta*eta)/3. - 11200.*xi*eta*eta*eta + 5600.*xi*(eta*eta*eta*eta) - 2800.*eta*eta/9. - 9520.*eta*eta*xi*xi + 11200.*(eta*eta)*(xi*xi*xi)/3. + 8400.*(xi*xi)*(eta*eta*eta) + 560.*(eta*eta*eta) - 280.*eta*eta*eta*eta, 6160.*eta*xi/9. - 14560.*eta*xi*xi/3. + 7840.*eta*(xi*xi*xi) - 11200.*eta*xi*xi*xi*xi/3. - 560.*xi/9. - 1680.*xi*eta*eta + 1120.*xi*(eta*eta*eta) + 10080.*(eta*eta)*(xi*xi) - 8400.*eta*eta*xi*xi*xi + 4760.*(xi*xi)/9. - 5600.*xi*xi*eta*eta*eta - 3920.*xi*xi*xi/3. + 3920.*(xi*xi*xi*xi)/3. - 1400.*xi*xi*xi*xi*xi/3.);
                case 34:
                  return RealGradient(6160.*eta*xi/9. - 560.*eta/9. - 1680.*eta*xi*xi + 1120.*eta*(xi*xi*xi) - 14560.*xi*eta*eta/3. + 7840.*xi*(eta*eta*eta) - 11200.*xi*eta*eta*eta*eta/3. + 4760.*(eta*eta)/9. + 10080.*(eta*eta)*(xi*xi) - 5600.*eta*eta*xi*xi*xi - 8400.*xi*xi*eta*eta*eta - 3920.*eta*eta*eta/3. + 3920.*(eta*eta*eta*eta)/3. - 1400.*eta*eta*eta*eta*eta/3., -6440.*eta*xi/9. + 19040.*eta*(xi*xi)/3. - 11200.*eta*xi*xi*xi + 5600.*eta*(xi*xi*xi*xi) + 280.*xi/9. + 1680.*xi*(eta*eta) - 4480.*xi*eta*eta*eta/3. + 1400.*xi*(eta*eta*eta*eta)/3. - 9520.*eta*eta*xi*xi + 8400.*(eta*eta)*(xi*xi*xi) - 2800.*xi*xi/9. + 11200.*(xi*xi)*(eta*eta*eta)/3. + 560.*(xi*xi*xi) - 280.*xi*xi*xi*xi);
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
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << totalorder);
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

  const Order totalorder = order + add_p_level*elem->p_level();
  libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

  const char sign = i >= totalorder * elem->n_edges() || elem->edge_orientation(i / totalorder) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / totalorder * 2 + 1) * totalorder - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (totalorder)
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
                        return sign * RealGradient(81.*eta/2. + 15.*xi/2. - 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 39. - 162.*(eta + 1.)*(eta + 1.)/(2.*2.) + 90.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 1:
                        return sign * RealGradient(-135.*eta/8. - 15.*xi/4. + 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 135.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 75.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 135./8. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2., 0.);
                      case 2:
                        return sign * RealGradient(27.*eta + 15.*xi/2. - 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 57./2. - 108.*(eta + 1.)*(eta + 1.)/(2.*2.) + 60.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 3:
                        return sign * RealGradient(0., -27.*eta/2. - 27.*xi + 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 135./4. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 135.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 4:
                        return sign * RealGradient(0., 45.*eta/8. + 9.*xi/2. - 90.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 90.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 9. - 45.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. - 45.*(xi + 1.)*(xi + 1.)/(2.*2.)/4.);
                      case 5:
                        return sign * RealGradient(0., -9.*eta - 9.*xi + 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 63./4. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 45.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 6:
                        return sign * RealGradient(9.*eta - 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 9. - 72.*(eta + 1.)*(eta + 1.)/(2.*2.) + 60.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 7:
                        return sign * RealGradient(-45.*eta/8. + 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 90.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 75.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 45./8. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2., 0.);
                      case 8:
                        return sign * RealGradient(27.*eta/2. - 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 27./2. - 108.*(eta + 1.)*(eta + 1.)/(2.*2.) + 90.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 9:
                        return sign * RealGradient(0., -27.*eta - 27.*xi/2. + 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 135./4. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 45.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 10:
                        return sign * RealGradient(0., 135.*eta/8. + 27.*xi/4. - 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 135.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 81./4. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. - 45.*(xi + 1.)*(xi + 1.)/(2.*2.)/4.);
                      case 11:
                        return sign * RealGradient(0., -81.*eta/2. - 81.*xi/2. + 324.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 243./4. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 135.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 12:
                        return RealGradient(0., -36.*eta - 81.*xi/2. + 324.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 117./2. + 60.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 135.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 13:
                        return RealGradient(0., 9.*eta + 27.*xi - 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 63./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 135.*(xi + 1.)*(xi + 1.)/(2.*2.)/2.);
                      case 14:
                        return RealGradient(36.*eta - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 36. - 162.*(eta + 1.)*(eta + 1.)/(2.*2.) + 90.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 15:
                        return RealGradient(-9.*eta + 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 9. + 108.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 90.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 16:
                        return RealGradient(24.*eta - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 24. - 108.*(eta + 1.)*(eta + 1.)/(2.*2.) + 60.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 17:
                        return RealGradient(-6.*eta + 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 6. + 72.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 60.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 18:
                        return RealGradient(0., -24.*eta - 27.*xi/2. + 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 63./2. + 60.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 45.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 19:
                        return RealGradient(0., 6.*eta + 9.*xi - 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 27./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 45.*(xi + 1.)*(xi + 1.)/(2.*2.)/2.);
                      case 20:
                        return RealGradient(0., 0.);
                      case 21:
                        return RealGradient(0., -9.*xi/2. - 5./2. + 15.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 22:
                        return RealGradient(0., 3.*xi + 5./2. - 15.*(xi + 1.)*(xi + 1.)/(2.*2.)/2.);
                      case 23:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(81.*eta/2. + 81.*xi/2. - 324.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 243./4. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 135.*(xi + 1.)*(xi + 1.)/(2.*2.)/2., 0.);
                      case 1:
                        return sign * RealGradient(-27.*eta/4. - 135.*xi/8. + 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 135.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 81./4. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 135.*((xi + 1.)*(xi + 1.)/(2.*2.))/4., 0.);
                      case 2:
                        return sign * RealGradient(27.*eta/2. + 27.*xi - 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 135./4. - 45.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 135.*(xi + 1.)*(xi + 1.)/(2.*2.)/2., 0.);
                      case 3:
                        return sign * RealGradient(0., -27.*xi/2. + 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 27./2. + 108.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 90.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 4:
                        return sign * RealGradient(0., 45.*xi/8. - 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 90.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 75.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 45./8. - 45.*(xi + 1.)*(xi + 1.)/(2.*2.) + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2.);
                      case 5:
                        return sign * RealGradient(0., -9.*xi + 45.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 9. + 72.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 60.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 6:
                        return sign * RealGradient(9.*eta + 9.*xi - 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 63./4. - 45.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 45.*(xi + 1.)*(xi + 1.)/(2.*2.)/2., 0.);
                      case 7:
                        return sign * RealGradient(-9.*eta/2. - 45.*xi/8. + 90.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 90.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 9. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 45.*((xi + 1.)*(xi + 1.)/(2.*2.))/4., 0.);
                      case 8:
                        return sign * RealGradient(27.*eta + 27.*xi/2. - 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 135./4. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 45.*(xi + 1.)*(xi + 1.)/(2.*2.)/2., 0.);
                      case 9:
                        return sign * RealGradient(0., -15.*eta/2. - 27.*xi + 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 57./2. + 108.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 60.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 10:
                        return sign * RealGradient(0., 15.*eta/4. + 135.*xi/8. - 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 135.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 75.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 135./8. - 135.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2.);
                      case 11:
                        return sign * RealGradient(0., -15.*eta/2. - 81.*xi/2. + 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 39. + 162.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 90.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 12:
                        return RealGradient(0., -36.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 36. + 162.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 90.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 13:
                        return RealGradient(0., 9.*xi - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 9. - 108.*(xi + 1.)*(xi + 1.)/(2.*2.) + 90.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 14:
                        return RealGradient(81.*eta/2. + 36.*xi - 324.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 117./2. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.), 0.);
                      case 15:
                        return RealGradient(-27.*eta - 9.*xi + 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 63./2. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.)), 0.);
                      case 16:
                        return RealGradient(27.*eta/2. + 24.*xi - 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 63./2. - 45.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.), 0.);
                      case 17:
                        return RealGradient(-9.*eta - 6.*xi + 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 27./2. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.)), 0.);
                      case 18:
                        return RealGradient(0., -24.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 24. + 108.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 60.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 19:
                        return RealGradient(0., 6.*xi - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 6. - 72.*(xi + 1.)*(xi + 1.)/(2.*2.) + 60.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 20:
                        return RealGradient(-9.*eta/2. - 5./2. + 15.*((eta + 1.)*(eta + 1.)/(2.*2.))/2., 0.);
                      case 21:
                        return RealGradient(0., 0.);
                      case 22:
                        return RealGradient(0., 0.);
                      case 23:
                        return RealGradient(3.*eta + 5./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.)/2., 0.);
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
                        return sign * RealGradient(-90.*eta*xi + 120.*eta + 60.*xi - 36. - 90.*eta*eta, 180.*eta*xi - 60.*eta - 120.*xi + 18. + 45.*(eta*eta) + 135.*(xi*xi));
                      case 1:
                        return sign * RealGradient(45.*eta*xi - 75.*eta/2. - 30.*xi + 15. + 45.*(eta*eta)/2., -45.*eta*xi + 105.*xi/2. - 21./4. + 45.*(eta*eta)/4. - 135.*xi*xi/2.);
                      case 2:
                        return sign * RealGradient(-90.*eta*xi + 30.*eta + 60.*xi - 24., -90.*xi + 9. + 135.*(xi*xi));
                      case 3:
                        return sign * RealGradient(-90.*eta*xi + 30.*eta, -90.*xi + 9. + 135.*(xi*xi));
                      case 4:
                        return sign * RealGradient(-45.*eta*xi/2. + 45.*eta/2. - 45.*eta*eta, 90.*eta*xi - 45.*eta/2. - 75.*xi/2. + 6. + 45.*(eta*eta)/4. + 135.*(xi*xi)/4.);
                      case 5:
                        return sign * RealGradient(0., -30.*eta + 3. + 45.*(eta*eta));
                      case 6:
                        return sign * RealGradient(0., -30.*eta + 3. + 45.*(eta*eta));
                      case 7:
                        return sign * RealGradient(-45.*eta*xi/2. + 45.*(eta*eta)/2., -45.*eta*xi + 75.*eta/2. - 30.*xi + 9./4. - 45.*eta*eta/2. + 135.*(xi*xi)/4.);
                      case 8:
                        return sign * RealGradient(-90.*eta*xi + 60.*eta - 90.*eta*eta, 180.*eta*xi - 120.*eta - 180.*xi + 54. + 45.*(eta*eta) + 135.*(xi*xi));
                      case 9:
                        return RealGradient(180.*eta*xi - 300.*eta + 360.*(eta*eta), -720.*eta*xi + 300.*eta + 300.*xi - 60. - 270.*eta*eta - 270.*xi*xi);
                      case 10:
                        return RealGradient(-540.*eta*xi + 300.*eta - 360.*eta*eta, 720.*eta*xi - 300.*eta - 900.*xi + 180. + 90.*(eta*eta) + 810.*(xi*xi));
                      case 11:
                        return RealGradient(-360.*eta*xi + 360.*eta - 360.*eta*eta, 720.*eta*xi - 120.*eta - 480.*xi + 60. + 540.*(xi*xi));
                      case 12:
                        return RealGradient(540.*eta*xi - 240.*eta + 180.*(eta*eta), -360.*eta*xi + 60.*eta + 720.*xi - 90. - 810.*xi*xi);
                      case 13:
                        return RealGradient(60.*eta - 180.*eta*eta, 360.*eta*xi - 240.*eta - 60.*xi + 30. + 270.*(eta*eta));
                      case 14:
                        return RealGradient(-120.*eta + 360.*(eta*eta), -720.*eta*xi + 360.*eta + 120.*xi - 60. - 180.*eta*eta);
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
                        return sign * RealGradient(-180.*eta*xi + 180.*eta + 120.*xi - 54. - 135.*eta*eta - 45.*xi*xi, 90.*eta*xi - 60.*xi + 90.*(xi*xi));
                      case 1:
                        return sign * RealGradient(45.*eta*xi + 30.*eta - 75.*xi/2. - 9./4. - 135.*eta*eta/4. + 45.*(xi*xi)/2., 45.*eta*xi/2. - 45.*xi*xi/2.);
                      case 2:
                        return sign * RealGradient(30.*xi - 3. - 45.*xi*xi, 0.);
                      case 3:
                        return sign * RealGradient(30.*xi - 3. - 45.*xi*xi, 0.);
                      case 4:
                        return sign * RealGradient(-90.*eta*xi + 75.*eta/2. + 45.*xi/2. - 6. - 135.*eta*eta/4. - 45.*xi*xi/4., 45.*eta*xi/2. - 45.*xi/2. + 45.*(xi*xi));
                      case 5:
                        return sign * RealGradient(90.*eta - 9. - 135.*eta*eta, 90.*eta*xi - 30.*xi);
                      case 6:
                        return sign * RealGradient(90.*eta - 9. - 135.*eta*eta, 90.*eta*xi - 60.*eta - 30.*xi + 24.);
                      case 7:
                        return sign * RealGradient(45.*eta*xi - 105.*eta/2. + 21./4. + 135.*(eta*eta)/2. - 45.*xi*xi/4., -45.*eta*xi + 30.*eta + 75.*xi/2. - 15. - 45.*xi*xi/2.);
                      case 8:
                        return sign * RealGradient(-180.*eta*xi + 120.*eta + 60.*xi - 18. - 135.*eta*eta - 45.*xi*xi, 90.*eta*xi - 60.*eta - 120.*xi + 36. + 90.*(xi*xi));
                      case 9:
                        return RealGradient(720.*eta*xi - 900.*eta - 300.*xi + 180. + 810.*(eta*eta) + 90.*(xi*xi), -540.*eta*xi + 300.*xi - 360.*xi*xi);
                      case 10:
                        return RealGradient(-720.*eta*xi + 300.*eta + 300.*xi - 60. - 270.*eta*eta - 270.*xi*xi, 180.*eta*xi - 300.*xi + 360.*(xi*xi));
                      case 11:
                        return RealGradient(-720.*eta*xi + 120.*eta + 360.*xi - 60. - 180.*xi*xi, -120.*xi + 360.*(xi*xi));
                      case 12:
                        return RealGradient(360.*eta*xi - 60.*eta - 240.*xi + 30. + 270.*(xi*xi), 60.*xi - 180.*xi*xi);
                      case 13:
                        return RealGradient(-360.*eta*xi + 720.*eta + 60.*xi - 90. - 810.*eta*eta, 540.*eta*xi - 240.*xi + 180.*(xi*xi));
                      case 14:
                        return RealGradient(720.*eta*xi - 480.*eta - 120.*xi + 60. + 540.*(eta*eta), -360.*eta*xi + 360.*xi - 360.*xi*xi);
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
                        return sign * RealGradient(240.*eta + 60.*xi - 1920.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1680.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 7200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 9600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 270. - 1800.*(eta + 1.)*(eta + 1.)/(2.*2.) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 105.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 2400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1050.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 1:
                        return sign * RealGradient(-760.*eta/9. - 215.*xi/9. + 6880.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 6160.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/9. - 8600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 34400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 15050.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. - 880./9. + 1900.*((eta + 1.)*(eta + 1.)/(2.*2.))/3. + 7700.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. + 385.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 30800.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 13475.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. - 7600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 3325.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9., 0.);
                      case 2:
                        return sign * RealGradient(400.*eta/9. + 170.*xi/9. - 5440.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 6160.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/9. + 6800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 27200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 11900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. + 520./9. - 1000.*(eta + 1.)*(eta + 1.)/(2.*2.)/3. - 7700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. - 385.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 30800.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 13475.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. + 4000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 1750.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9., 0.);
                      case 3:
                        return sign * RealGradient(-120.*eta - 45.*xi + 1440.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1680.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 5400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 7200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 150. + 900.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 105.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 4:
                        return sign * RealGradient(0., 60.*eta + 120.*xi - 1800.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 5400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 164. - 240.*(eta + 1.)*(eta + 1.)/(2.*2.) - 10800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 720.*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 140.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 560.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 5:
                        return sign * RealGradient(0., -190.*eta/9. - 170.*xi/9. + 1900.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 1900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 13300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 4300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 7700.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 1012./27. + 860.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 4300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 30100.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 340.*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 7700.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/3. - 1540.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 53900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 2380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27.);
                      case 6:
                        return sign * RealGradient(0., 100.*eta/9. + 80.*xi/9. - 1000.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 1000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 3400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 7700.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 508./27. - 680.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 3400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 23800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 160.*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7700.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/3. + 1540.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 53900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 1120.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27.);
                      case 7:
                        return sign * RealGradient(0., -30.*eta - 30.*xi + 900.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2700.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2700.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 56. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 8100.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 180.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 140.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 8:
                        return sign * RealGradient(30.*eta - 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 420.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 2700.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 5400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 30. - 450.*(eta + 1.)*(eta + 1.)/(2.*2.) - 3150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 9:
                        return sign * RealGradient(-100.*eta/9. + 1360.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 1540.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/9. - 3400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 6800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/3. - 11900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. - 100./9. + 500.*((eta + 1.)*(eta + 1.)/(2.*2.))/3. + 3850.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 7700.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/3. + 13475.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. - 1000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/3. + 1750.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9., 0.);
                      case 10:
                        return sign * RealGradient(190.*eta/9. - 1720.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 1540.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/9. + 4300.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 8600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/3. + 15050.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/9. + 190./9. - 950.*(eta + 1.)*(eta + 1.)/(2.*2.)/3. - 3850.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7700.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/3. - 13475.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9. + 1900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/3. - 3325.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/9., 0.);
                      case 11:
                        return sign * RealGradient(-60.*eta + 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 420.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 3600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 7200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 60. + 900.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 3150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1050.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 12:
                        return sign * RealGradient(0., 120.*eta + 60.*xi - 1800.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 5400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 164. - 720.*(eta + 1.)*(eta + 1.)/(2.*2.) - 10800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 240.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 140.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 13:
                        return sign * RealGradient(0., -400.*eta/9. - 160.*xi/9. + 2000.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 4000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 6800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 15400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 1552./27. + 2720.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 13600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 23800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 640.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 30800.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. - 6160.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 53900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 1120.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27.);
                      case 14:
                        return sign * RealGradient(0., 760.*eta/9. + 340.*xi/9. - 3800.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 7600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 13300.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 8600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 15400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 3028./27. - 3440.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 17200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 30100.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 1360.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 30800.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. + 6160.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 53900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 2380.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27.);
                      case 15:
                        return sign * RealGradient(0., -240.*eta - 240.*xi + 3600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 7200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 7200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 416. + 960.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 14400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 960.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 560.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 16:
                        return RealGradient(0., -195.*eta - 232.*xi + 3480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 7200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 6960.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 4060.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 375. + 780.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 14400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 960.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 455.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 560.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 17:
                        return RealGradient(0., -45.*eta - 112.*xi + 1680.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 5400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 1960.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 145. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 10800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 720.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 105.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 560.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 18:
                        return RealGradient(0., -60.*eta - 16.*xi + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 480.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 280.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 60. + 240.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.));
                      case 19:
                        return RealGradient(195.*eta - 1560.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1365.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 6960.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 9600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 195. - 1740.*(eta + 1.)*(eta + 1.)/(2.*2.) - 6090.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 2400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1050.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 20:
                        return RealGradient(45.*eta - 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 315.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 3360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 7200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 45. - 840.*(eta + 1.)*(eta + 1.)/(2.*2.) - 2940.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1050.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 21:
                        return RealGradient(60.*eta - 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 420.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 480.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 60. - 120.*(eta + 1.)*(eta + 1.)/(2.*2.) - 420.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.), 0.);
                      case 22:
                        return RealGradient(-195.*eta/2. + 1170.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1365.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 5220.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 7200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 195./2. + 870.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 6090.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1200.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 23:
                        return RealGradient(-45.*eta/2. + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 315.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 2520.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 5400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 45./2. + 420.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 2940.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 24:
                        return RealGradient(-30.*eta + 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 420.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 30. + 60.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 420.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)), 0.);
                      case 25:
                        return RealGradient(0., 195.*eta/2. + 58.*xi - 1740.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 5220.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4060.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 285./2. - 585.*(eta + 1.)*(eta + 1.)/(2.*2.) - 10800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 240.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 455.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 140.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 26:
                        return RealGradient(0., 45.*eta/2. + 28.*xi - 840.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2700.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2520.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 1960.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 95./2. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.) - 8100.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 180.*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 105.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 140.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 27:
                        return RealGradient(0., 30.*eta + 4.*xi - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 280.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 30. - 180.*(eta + 1.)*(eta + 1.)/(2.*2.) + 140.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)));
                      case 28:
                        return RealGradient(-9.*eta - 9. + 171.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 120.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 105.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2., 0.);
                      case 29:
                        return RealGradient(9.*eta + 9. - 171.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 120.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 105.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 30:
                        return RealGradient(0., -9.*eta - 57.*xi + 171.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 210.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 54. + 240.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 140.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 31:
                        return RealGradient(0., 9.*eta + 57.*xi/2. - 171.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 360.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 210.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 63./2. - 120.*(xi + 1.)*(xi + 1.)/(2.*2.) + 70.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 32:
                        return RealGradient(0., -3.*eta/2. - 27.*xi + 81.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 210.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 53./2. + 180.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 140.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 33:
                        return RealGradient(0., 3.*eta/2. + 27.*xi/2. - 81.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 210.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14. - 90.*(xi + 1.)*(xi + 1.)/(2.*2.) + 70.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 34:
                        return RealGradient(-3.*eta/2. - 3./2. + 81.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 90.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 105.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2., 0.);
                      case 35:
                        return RealGradient(3.*eta/2. + 3./2. - 81.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 90.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 105.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 36:
                        return RealGradient(0., -9.*eta/2. - 6.*xi + 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 9./2.);
                      case 37:
                        return RealGradient(9.*eta/2. + 9./2. - 9.*(eta + 1.)*(eta + 1.)/(2.*2.), 0.);
                      case 38:
                        return RealGradient(-9.*eta/2. - 9./2. + 9.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 39:
                        return RealGradient(0., 9.*eta/2. + 3.*xi - 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 9./2.);
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
                        return sign * RealGradient(240.*eta + 240.*xi - 3600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 7200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 7200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 416. - 960.*(eta + 1.)*(eta + 1.)/(2.*2.) - 14400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 960.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 560.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 1:
                        return sign * RealGradient(-340.*eta/9. - 760.*xi/9. + 3800.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 8600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 15400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 7600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 13300.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 3028./27. + 1360.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 17200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 30800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 3440.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 30100.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. - 2380.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 53900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 6160.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27., 0.);
                      case 2:
                        return sign * RealGradient(160.*eta/9. + 400.*xi/9. - 2000.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 6800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 15400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 4000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 7000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 1552./27. - 640.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 13600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 30800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 2720.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 23800.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. + 1120.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 53900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 6160.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27., 0.);
                      case 3:
                        return sign * RealGradient(-60.*eta - 120.*xi + 1800.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 5400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 164. + 240.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 10800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 720.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 560.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 4:
                        return sign * RealGradient(0., 60.*xi - 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 420.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 60. - 3150.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 900.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1800.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1050.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 5:
                        return sign * RealGradient(0., -190.*xi/9. + 1720.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 4300.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 8600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/3. - 15050.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 1540.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 190./9. + 3850.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 7700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/3. + 13475.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 950.*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 1900.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/3. + 3325.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9.);
                      case 6:
                        return sign * RealGradient(0., 100.*xi/9. - 1360.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 3400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 6800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/3. + 11900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 1540.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 100./9. - 3850.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7700.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/3. - 13475.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 500.*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 1000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/3. - 1750.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9.);
                      case 7:
                        return sign * RealGradient(0., -30.*xi + 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2700.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 420.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 30. + 3150.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 450.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 900.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 8:
                        return sign * RealGradient(30.*eta + 30.*xi - 900.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2700.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2700.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 56. - 180.*(eta + 1.)*(eta + 1.)/(2.*2.) - 8100.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 180.*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 140.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 140.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 9:
                        return sign * RealGradient(-80.*eta/9. - 100.*xi/9. + 1000.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 3400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7700.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 1000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 7000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 508./27. + 160.*((eta + 1.)*(eta + 1.)/(2.*2.))/3. + 3400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/3. + 680.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 23800.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. - 1120.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/27. + 53900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27. - 1540.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27., 0.);
                      case 10:
                        return sign * RealGradient(170.*eta/9. + 190.*xi/9. - 1900.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 4300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 7700.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 1900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 13300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 1012./27. - 340.*(eta + 1.)*(eta + 1.)/(2.*2.)/3. - 4300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7700.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/3. - 860.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 30100.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. + 2380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/27. - 53900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/27. + 1540.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/27., 0.);
                      case 11:
                        return sign * RealGradient(-120.*eta - 60.*xi + 1800.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 164. + 720.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 10800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 240.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 560.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 140.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 12:
                        return sign * RealGradient(0., 45.*eta + 120.*xi - 1440.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 5400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1680.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 150. - 105.*(eta + 1.)*(eta + 1.)/(2.*2.) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 900.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 13:
                        return sign * RealGradient(0., -170.*eta/9. - 400.*xi/9. + 5440.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. - 6800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 27200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 11900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 6160.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 520./9. + 385.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 7700.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 30800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 13475.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 1000.*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 4000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 1750.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9.);
                      case 14:
                        return sign * RealGradient(0., 215.*eta/9. + 760.*xi/9. - 6880.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/9. + 8600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 34400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 15050.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/9. + 6160.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 880./9. - 385.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 7700.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 30800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 13475.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9. - 1900.*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 7600.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 3325.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/9.);
                      case 15:
                        return sign * RealGradient(0., -60.*eta - 240.*xi + 1920.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 7200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 9600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1680.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 270. + 105.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1800.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1050.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 16:
                        return RealGradient(0., -195.*xi + 1560.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 6960.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 9600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1365.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 195. + 6090.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1740.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1050.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 17:
                        return RealGradient(0., -45.*xi + 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 315.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 45. + 2940.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 840.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1800.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1050.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 18:
                        return RealGradient(0., -60.*xi + 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 480.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 420.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 60. + 420.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 120.*((xi + 1.)*(xi + 1.)/(2.*2.)));
                      case 19:
                        return RealGradient(232.*eta + 195.*xi - 3480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 6960.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4060.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 7200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 375. - 960.*(eta + 1.)*(eta + 1.)/(2.*2.) - 14400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 780.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 455.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 20:
                        return RealGradient(112.*eta + 45.*xi - 1680.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3360.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1960.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 5400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 145. - 720.*(eta + 1.)*(eta + 1.)/(2.*2.) - 10800.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 180.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 560.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4900.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 105.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 21:
                        return RealGradient(16.*eta + 60.*xi - 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 480.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 280.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 60. - 240.*(xi + 1.)*(xi + 1.)/(2.*2.) + 140.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 22:
                        return RealGradient(-58.*eta - 195.*xi/2. + 1740.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 5220.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4060.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 285./2. + 240.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 10800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 585.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 455.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 23:
                        return RealGradient(-28.*eta - 45.*xi/2. + 840.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2520.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1960.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2700.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 2100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 95./2. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 8100.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 135.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4900.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 105.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 24:
                        return RealGradient(-4.*eta - 30.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 280.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 30. + 180.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 140.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 25:
                        return RealGradient(0., 195.*xi/2. - 1170.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 5220.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1365.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 195./2. - 6090.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 870.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 26:
                        return RealGradient(0., 45.*xi/2. - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2520.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 315.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 45./2. - 2940.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 420.*(xi + 1.)*(xi + 1.)/(2.*2.) + 900.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 27:
                        return RealGradient(0., 30.*xi - 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 360.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 420.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 30. - 420.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 28:
                        return RealGradient(-57.*eta - 9.*xi + 171.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 210.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 54. + 240.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 29:
                        return RealGradient(57.*eta/2. + 9.*xi - 171.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 210.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 63./2. - 120.*(eta + 1.)*(eta + 1.)/(2.*2.) + 70.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 30:
                        return RealGradient(0., -9.*xi - 9. + 171.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 120.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 105.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 31:
                        return RealGradient(0., 9.*xi + 9. - 171.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 120.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 105.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2.);
                      case 32:
                        return RealGradient(0., -3.*xi/2. - 3./2. + 81.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 90.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 105.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 33:
                        return RealGradient(0., 3.*xi/2. + 3./2. - 81.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 90.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 105.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2.);
                      case 34:
                        return RealGradient(-27.*eta - 3.*xi/2. + 81.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 210.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 53./2. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 140.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 35:
                        return RealGradient(27.*eta/2. + 3.*xi/2. - 81.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 210.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 14. - 90.*(eta + 1.)*(eta + 1.)/(2.*2.) + 70.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 36:
                        return RealGradient(0., -9.*xi/2. - 9./2. + 9.*((xi + 1.)*(xi + 1.)/(2.*2.)));
                      case 37:
                        return RealGradient(6.*eta + 9.*xi/2. - 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 9./2., 0.);
                      case 38:
                        return RealGradient(-3.*eta - 9.*xi/2. + 18.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 9./2., 0.);
                      case 39:
                        return RealGradient(0., 9.*xi/2. + 9./2. - 9.*(xi + 1.)*(xi + 1.)/(2.*2.));
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
                        return sign * RealGradient(-1680.*eta*xi + 720.*eta + 672.*eta*(xi*xi) + 480.*xi + 1344.*xi*(eta*eta) - 120. - 1260.*eta*eta - 420.*xi*xi + 672.*(eta*eta*eta), 1680.*eta*xi - 240.*eta - 2016.*eta*xi*xi - 480.*xi - 1344.*xi*eta*eta + 40. + 420.*(eta*eta) + 1260.*(xi*xi) - 224.*eta*eta*eta - 896.*xi*xi*xi);
                      case 1:
                        return sign * RealGradient(5768.*eta*xi/9. - 200.*eta - 2464.*eta*xi*xi/9. - 1720.*xi/9. - 4480.*xi*eta*eta/9. + 380./9. + 784.*(eta*eta)/3. + 1540.*(xi*xi)/9. - 896.*eta*eta*eta/9., -4256.*eta*xi/9. + 176.*eta/9. + 2240.*eta*(xi*xi)/3. + 1600.*xi/9. + 1792.*xi*(eta*eta)/9. - 296./27. + 448.*(eta*eta)/9. - 504.*xi*xi - 1792.*eta*eta*eta/27. + 9856.*(xi*xi*xi)/27.);
                      case 2:
                        return sign * RealGradient(-2240.*eta*xi/9. - 16.*eta + 2464.*eta*(xi*xi)/9. + 1360.*xi/9. + 448.*xi*(eta*eta)/9. - 200./9. + 476.*(eta*eta)/3. - 1540.*xi*xi/9. - 1120.*eta*eta*eta/9., -784.*eta*xi/9. + 112.*eta/9. - 224.*eta*xi*xi/3. - 1024.*xi/9. + 2240.*xi*(eta*eta)/9. + 152./27. - 196.*eta*eta/9. + 420.*(xi*xi) - 224.*eta*eta*eta/27. - 9856.*xi*xi*xi/27.);
                      case 3:
                        return sign * RealGradient(504.*eta*xi - 72.*eta - 672.*eta*xi*xi - 360.*xi + 60. + 420.*(xi*xi), 288.*xi - 16. - 1008.*xi*xi + 896.*(xi*xi*xi));
                      case 4:
                        return sign * RealGradient(504.*eta*xi - 72.*eta - 672.*eta*xi*xi, 288.*xi - 16. - 1008.*xi*xi + 896.*(xi*xi*xi));
                      case 5:
                        return sign * RealGradient(896.*eta*xi/3. - 80.*eta - 1792.*eta*xi*xi/9. - 1792.*xi*eta*eta/3. + 784.*(eta*eta)/3. - 448.*eta*eta*eta/3., -1792.*eta*xi/3. + 56.*eta + 896.*eta*(xi*xi) + 160.*xi + 896.*xi*(eta*eta)/3. - 12. - 140.*eta*eta/3. - 1232.*xi*xi/3. + 224.*(eta*eta*eta)/27. + 7168.*(xi*xi*xi)/27.);
                      case 6:
                        return sign * RealGradient(280.*eta*xi/3. - 56.*eta - 224.*eta*xi*xi/9. - 896.*xi*eta*eta/3. + 896.*(eta*eta)/3. - 896.*eta*eta*eta/3., -1568.*eta*xi/3. + 80.*eta + 448.*eta*(xi*xi) + 64.*xi + 1792.*xi*(eta*eta)/3. - 8. - 448.*eta*eta/3. - 280.*xi*xi/3. + 1792.*(eta*eta*eta)/27. + 896.*(xi*xi*xi)/27.);
                      case 7:
                        return sign * RealGradient(0., 72.*eta - 4. - 252.*eta*eta + 224.*(eta*eta*eta));
                      case 8:
                        return sign * RealGradient(0., 72.*eta - 4. - 252.*eta*eta + 224.*(eta*eta*eta));
                      case 9:
                        return sign * RealGradient(392.*eta*xi/9. - 112.*eta/9. + 224.*eta*(xi*xi)/9. - 2240.*xi*eta*eta/9. + 392.*(eta*eta)/9. + 224.*(eta*eta*eta)/9., -952.*eta*xi/3. + 16.*eta + 1120.*eta*(xi*xi)/3. + 208.*xi/9. - 448.*xi*eta*eta/9. - 112./27. + 1120.*(eta*eta)/9. + 56.*(xi*xi)/9. - 2464.*eta*eta*eta/27. - 896.*xi*xi*xi/27.);
                      case 10:
                        return sign * RealGradient(-896.*eta*xi/9. - 176.*eta/9. + 1792.*eta*(xi*xi)/9. - 1792.*xi*eta*eta/9. + 2128.*(eta*eta)/9. - 2240.*eta*eta*eta/9., -1568.*eta*xi/3. + 200.*eta + 896.*eta*(xi*xi)/3. - 1216.*xi/9. + 4480.*xi*(eta*eta)/9. + 76./27. - 2884.*eta*eta/9. + 3472.*(xi*xi)/9. + 2464.*(eta*eta*eta)/27. - 7168.*xi*xi*xi/27.);
                      case 11:
                        return sign * RealGradient(-840.*eta*xi + 240.*eta + 672.*eta*(xi*xi) + 1344.*xi*(eta*eta) - 840.*eta*eta + 672.*(eta*eta*eta), 2520.*eta*xi - 720.*eta - 2016.*eta*xi*xi - 960.*xi - 1344.*xi*eta*eta + 160. + 840.*(eta*eta) + 1680.*(xi*xi) - 224.*eta*eta*eta - 896.*xi*xi*xi);
                      case 12:
                        return RealGradient(6048.*eta*xi - 3240.*eta - 2016.*eta*xi*xi - 8064.*xi*eta*eta + 9072.*(eta*eta) - 6048.*eta*eta*eta, -12096.*eta*xi + 2160.*eta + 12096.*eta*(xi*xi) + 2160.*xi + 12096.*xi*(eta*eta) - 240. - 4536.*eta*eta - 4536.*xi*xi + 2688.*(eta*eta*eta) + 2688.*(xi*xi*xi));
                      case 13:
                        return RealGradient(-9072.*eta*xi + 2160.*eta + 8064.*eta*(xi*xi) + 12096.*xi*(eta*eta) - 6048.*eta*eta + 4032.*(eta*eta*eta), 18144.*eta*xi - 3240.*eta - 18144.*eta*xi*xi - 8640.*xi - 8064.*xi*eta*eta + 960. + 3024.*(eta*eta) + 18144.*(xi*xi) - 672.*eta*eta*eta - 10752.*xi*xi*xi);
                      case 14:
                        return RealGradient(9072.*eta*xi - 1944.*eta - 6048.*eta*xi*xi - 8064.*xi*eta*eta + 2016.*(eta*eta), -6048.*eta*xi + 432.*eta + 12096.*eta*(xi*xi) + 3456.*xi - 216. - 10584.*xi*xi + 8064.*(xi*xi*xi));
                      case 15:
                        return RealGradient(-7056.*eta*xi + 1152.*eta + 8064.*eta*(xi*xi) + 4032.*xi*(eta*eta) - 1008.*eta*eta, 3024.*eta*xi - 216.*eta - 6048.*eta*xi*xi - 4608.*xi + 288. + 14112.*(xi*xi) - 10752.*xi*xi*xi);
                      case 16:
                        return RealGradient(-216.*eta + 1512.*(eta*eta) - 2016.*eta*eta*eta, -2016.*eta*xi + 1152.*eta + 144.*xi + 4032.*xi*(eta*eta) - 72. - 3528.*eta*eta + 2688.*(eta*eta*eta));
                      case 17:
                        return RealGradient(432.*eta - 3024.*eta*eta + 4032.*(eta*eta*eta), 4032.*eta*xi - 1944.*eta - 288.*xi - 8064.*xi*eta*eta + 144. + 4536.*(eta*eta) - 2016.*eta*eta*eta);
                      case 18:
                        return RealGradient(3276.*eta*xi - 1332.*eta - 1512.*eta*xi*xi - 6048.*xi*eta*eta + 4788.*(eta*eta) - 3528.*eta*eta*eta, -8064.*eta*xi + 1044.*eta + 9072.*eta*(xi*xi) + 1548.*xi + 7056.*xi*(eta*eta) - 144. - 1638.*eta*eta - 3402.*xi*xi + 672.*(eta*eta*eta) + 2016.*(xi*xi*xi));
                      case 19:
                        return RealGradient(-3276.*eta*xi + 1044.*eta + 2016.*eta*(xi*xi) + 7056.*xi*(eta*eta) - 4032.*eta*eta + 3024.*(eta*eta*eta), 9576.*eta*xi - 1332.*eta - 10584.*eta*xi*xi - 2196.*xi - 6048.*xi*eta*eta + 216. + 1638.*(eta*eta) + 4662.*(xi*xi) - 504.*eta*eta*eta - 2688.*xi*xi*xi);
                      case 20:
                        return RealGradient(1008.*eta*xi - 216.*eta - 504.*eta*xi*xi - 756.*eta*eta + 1008.*(eta*eta*eta), 1008.*eta*xi - 360.*eta + 144.*xi - 2016.*xi*eta*eta + 12. + 1008.*(eta*eta) - 756.*xi*xi - 672.*eta*eta*eta + 672.*(xi*xi*xi));
                      case 21:
                        return RealGradient(-756.*eta*xi - 216.*eta + 2016.*eta*(xi*xi) - 3024.*xi*eta*eta + 2268.*(eta*eta) - 2016.*eta*eta*eta, -5544.*eta*xi + 1188.*eta + 4536.*eta*(xi*xi) - 936.*xi + 4032.*xi*(eta*eta) + 6. - 1512.*eta*eta + 3402.*(xi*xi) + 336.*(eta*eta*eta) - 2688.*xi*xi*xi);
                      case 22:
                        return RealGradient(-3024.*eta*xi + 1188.*eta + 1008.*eta*(xi*xi) + 4032.*xi*(eta*eta) - 2772.*eta*eta + 1512.*(eta*eta*eta), 4536.*eta*xi - 216.*eta - 6048.*eta*xi*xi - 972.*xi - 3024.*xi*eta*eta + 66. - 378.*eta*eta + 2268.*(xi*xi) + 672.*(eta*eta*eta) - 1344.*xi*xi*xi);
                      case 23:
                        return RealGradient(2016.*eta*xi - 360.*eta - 2016.*eta*xi*xi - 2016.*xi*eta*eta + 504.*(eta*eta), -1512.*eta*xi - 216.*eta + 3024.*eta*(xi*xi) + 1440.*xi - 48. + 504.*(eta*eta) - 4032.*xi*xi - 168.*eta*eta*eta + 2688.*(xi*xi*xi));
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
                        return sign * RealGradient(-2520.*eta*xi + 960.*eta + 1344.*eta*(xi*xi) + 720.*xi + 2016.*xi*(eta*eta) - 160. - 1680.*eta*eta - 840.*xi*xi + 896.*(eta*eta*eta) + 224.*(xi*xi*xi), 840.*eta*xi - 1344.*eta*xi*xi - 240.*xi - 672.*xi*eta*eta + 840.*(xi*xi) - 672.*xi*xi*xi);
                      case 1:
                        return sign * RealGradient(1568.*eta*xi/3. + 1216.*eta/9. - 4480.*eta*xi*xi/9. - 200.*xi - 896.*xi*eta*eta/3. - 76./27. - 3472.*eta*eta/9. + 2884.*(xi*xi)/9. + 7168.*(eta*eta*eta)/27. - 2464.*xi*xi*xi/27., 896.*eta*xi/9. + 1792.*eta*(xi*xi)/9. + 176.*xi/9. - 1792.*xi*eta*eta/9. - 2128.*xi*xi/9. + 2240.*(xi*xi*xi)/9.);
                      case 2:
                        return sign * RealGradient(952.*eta*xi/3. - 208.*eta/9. + 448.*eta*(xi*xi)/9. - 16.*xi - 1120.*xi*eta*eta/3. + 112./27. - 56.*eta*eta/9. - 1120.*xi*xi/9. + 896.*(eta*eta*eta)/27. + 2464.*(xi*xi*xi)/27., -392.*eta*xi/9. + 2240.*eta*(xi*xi)/9. + 112.*xi/9. - 224.*xi*eta*eta/9. - 392.*xi*xi/9. - 224.*xi*xi*xi/9.);
                      case 3:
                        return sign * RealGradient(-72.*xi + 4. + 252.*(xi*xi) - 224.*xi*xi*xi, 0.);
                      case 4:
                        return sign * RealGradient(-72.*xi + 4. + 252.*(xi*xi) - 224.*xi*xi*xi, 0.);
                      case 5:
                        return sign * RealGradient(1568.*eta*xi/3. - 64.*eta - 1792.*eta*xi*xi/3. - 80.*xi - 448.*xi*eta*eta + 8. + 280.*(eta*eta)/3. + 448.*(xi*xi)/3. - 896.*eta*eta*eta/27. - 1792.*xi*xi*xi/27., -280.*eta*xi/3. + 896.*eta*(xi*xi)/3. + 56.*xi + 224.*xi*(eta*eta)/9. - 896.*xi*xi/3. + 896.*(xi*xi*xi)/3.);
                      case 6:
                        return sign * RealGradient(1792.*eta*xi/3. - 160.*eta - 896.*eta*xi*xi/3. - 56.*xi - 896.*xi*eta*eta + 12. + 1232.*(eta*eta)/3. + 140.*(xi*xi)/3. - 7168.*eta*eta*eta/27. - 224.*xi*xi*xi/27., -896.*eta*xi/3. + 1792.*eta*(xi*xi)/3. + 80.*xi + 1792.*xi*(eta*eta)/9. - 784.*xi*xi/3. + 448.*(xi*xi*xi)/3.);
                      case 7:
                        return sign * RealGradient(-288.*eta + 16. + 1008.*(eta*eta) - 896.*eta*eta*eta, -504.*eta*xi + 72.*xi + 672.*xi*(eta*eta));
                      case 8:
                        return sign * RealGradient(-288.*eta + 16. + 1008.*(eta*eta) - 896.*eta*eta*eta, -504.*eta*xi + 360.*eta + 72.*xi + 672.*xi*(eta*eta) - 60. - 420.*eta*eta);
                      case 9:
                        return sign * RealGradient(784.*eta*xi/9. + 1024.*eta/9. - 2240.*eta*xi*xi/9. - 112.*xi/9. + 224.*xi*(eta*eta)/3. - 152./27. - 420.*eta*eta + 196.*(xi*xi)/9. + 9856.*(eta*eta*eta)/27. + 224.*(xi*xi*xi)/27., 2240.*eta*xi/9. - 1360.*eta/9. - 448.*eta*xi*xi/9. + 16.*xi - 2464.*xi*eta*eta/9. + 200./9. + 1540.*(eta*eta)/9. - 476.*xi*xi/3. + 1120.*(xi*xi*xi)/9.);
                      case 10:
                        return sign * RealGradient(4256.*eta*xi/9. - 1600.*eta/9. - 1792.*eta*xi*xi/9. - 176.*xi/9. - 2240.*xi*eta*eta/3. + 296./27. + 504.*(eta*eta) - 448.*xi*xi/9. - 9856.*eta*eta*eta/27. + 1792.*(xi*xi*xi)/27., -5768.*eta*xi/9. + 1720.*eta/9. + 4480.*eta*(xi*xi)/9. + 200.*xi + 2464.*xi*(eta*eta)/9. - 380./9. - 1540.*eta*eta/9. - 784.*xi*xi/3. + 896.*(xi*xi*xi)/9.);
                      case 11:
                        return sign * RealGradient(-1680.*eta*xi + 480.*eta + 1344.*eta*(xi*xi) + 240.*xi + 2016.*xi*(eta*eta) - 40. - 1260.*eta*eta - 420.*xi*xi + 896.*(eta*eta*eta) + 224.*(xi*xi*xi), 1680.*eta*xi - 480.*eta - 1344.*eta*xi*xi - 720.*xi - 672.*xi*eta*eta + 120. + 420.*(eta*eta) + 1260.*(xi*xi) - 672.*xi*xi*xi);
                      case 12:
                        return RealGradient(18144.*eta*xi - 8640.*eta - 8064.*eta*xi*xi - 3240.*xi - 18144.*xi*eta*eta + 960. + 18144.*(eta*eta) + 3024.*(xi*xi) - 10752.*eta*eta*eta - 672.*xi*xi*xi, -9072.*eta*xi + 12096.*eta*(xi*xi) + 2160.*xi + 8064.*xi*(eta*eta) - 6048.*xi*xi + 4032.*(xi*xi*xi));
                      case 13:
                        return RealGradient(-12096.*eta*xi + 2160.*eta + 12096.*eta*(xi*xi) + 2160.*xi + 12096.*xi*(eta*eta) - 240. - 4536.*eta*eta - 4536.*xi*xi + 2688.*(eta*eta*eta) + 2688.*(xi*xi*xi), 6048.*eta*xi - 8064.*eta*xi*xi - 3240.*xi - 2016.*xi*eta*eta + 9072.*(xi*xi) - 6048.*xi*xi*xi);
                      case 14:
                        return RealGradient(4032.*eta*xi - 288.*eta - 8064.*eta*xi*xi - 1944.*xi + 144. + 4536.*(xi*xi) - 2016.*xi*xi*xi, 432.*xi - 3024.*xi*xi + 4032.*(xi*xi*xi));
                      case 15:
                        return RealGradient(-2016.*eta*xi + 144.*eta + 4032.*eta*(xi*xi) + 1152.*xi - 72. - 3528.*xi*xi + 2688.*(xi*xi*xi), -216.*xi + 1512.*(xi*xi) - 2016.*xi*xi*xi);
                      case 16:
                        return RealGradient(3024.*eta*xi - 4608.*eta - 216.*xi - 6048.*xi*eta*eta + 288. + 14112.*(eta*eta) - 10752.*eta*eta*eta, -7056.*eta*xi + 4032.*eta*(xi*xi) + 1152.*xi + 8064.*xi*(eta*eta) - 1008.*xi*xi);
                      case 17:
                        return RealGradient(-6048.*eta*xi + 3456.*eta + 432.*xi + 12096.*xi*(eta*eta) - 216. - 10584.*eta*eta + 8064.*(eta*eta*eta), 9072.*eta*xi - 8064.*eta*xi*xi - 1944.*xi - 6048.*xi*eta*eta + 2016.*(xi*xi));
                      case 18:
                        return RealGradient(9576.*eta*xi - 2196.*eta - 6048.*eta*xi*xi - 1332.*xi - 10584.*xi*eta*eta + 216. + 4662.*(eta*eta) + 1638.*(xi*xi) - 2688.*eta*eta*eta - 504.*xi*xi*xi, -3276.*eta*xi + 7056.*eta*(xi*xi) + 1044.*xi + 2016.*xi*(eta*eta) - 4032.*xi*xi + 3024.*(xi*xi*xi));
                      case 19:
                        return RealGradient(-8064.*eta*xi + 1548.*eta + 7056.*eta*(xi*xi) + 1044.*xi + 9072.*xi*(eta*eta) - 144. - 3402.*eta*eta - 1638.*xi*xi + 2016.*(eta*eta*eta) + 672.*(xi*xi*xi), 3276.*eta*xi - 6048.*eta*xi*xi - 1332.*xi - 1512.*xi*eta*eta + 4788.*(xi*xi) - 3528.*xi*xi*xi);
                      case 20:
                        return RealGradient(-1512.*eta*xi + 1440.*eta - 216.*xi + 3024.*xi*(eta*eta) - 48. - 4032.*eta*eta + 504.*(xi*xi) + 2688.*(eta*eta*eta) - 168.*xi*xi*xi, 2016.*eta*xi - 2016.*eta*xi*xi - 360.*xi - 2016.*xi*eta*eta + 504.*(xi*xi));
                      case 21:
                        return RealGradient(4536.*eta*xi - 972.*eta - 3024.*eta*xi*xi - 216.*xi - 6048.*xi*eta*eta + 66. + 2268.*(eta*eta) - 378.*xi*xi - 1344.*eta*eta*eta + 672.*(xi*xi*xi), -3024.*eta*xi + 4032.*eta*(xi*xi) + 1188.*xi + 1008.*xi*(eta*eta) - 2772.*xi*xi + 1512.*(xi*xi*xi));
                      case 22:
                        return RealGradient(-5544.*eta*xi - 936.*eta + 4032.*eta*(xi*xi) + 1188.*xi + 4536.*xi*(eta*eta) + 6. + 3402.*(eta*eta) - 1512.*xi*xi - 2688.*eta*eta*eta + 336.*(xi*xi*xi), -756.*eta*xi - 3024.*eta*xi*xi - 216.*xi + 2016.*xi*(eta*eta) + 2268.*(xi*xi) - 2016.*xi*xi*xi);
                      case 23:
                        return RealGradient(1008.*eta*xi + 144.*eta - 2016.*eta*xi*xi - 360.*xi + 12. - 756.*eta*eta + 1008.*(xi*xi) + 672.*(eta*eta*eta) - 672.*xi*xi*xi, 1008.*eta*xi - 216.*xi - 504.*xi*eta*eta - 756.*xi*xi + 1008.*(xi*xi*xi));
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
                        return sign * RealGradient(1875.*eta/2. + 525.*xi/2. - 13125.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 26250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 15750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 78750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 183750.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 183750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 66150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 1125. - 11250.*(eta + 1.)*(eta + 1.)/(2.*2.) - 157500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 94500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1050.*(xi + 1.)*(xi + 1.)/(2.*2.) + 367500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 367500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 132300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 26250.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 630.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 79380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 26250.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 1:
                        return sign * RealGradient(-71625.*eta/256. - 21105.*xi/256. + 527625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. - 1021125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 291375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 1582875.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 3693375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 3693375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 1329615.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. - 10875./32. + 214875.*((eta + 1.)*(eta + 1.)/(2.*2.))/64. + 3063375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 874125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 40845.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 7147875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 7147875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 2573235.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64. - 501375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 11655.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. - 2039625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/32. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/32. + 501375.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 180495.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64., 0.);
                      case 2:
                        return sign * RealGradient(2625.*eta/16. + 1155.*xi/16. - 28875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/8. + 70875.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 23625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 86625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. - 202125.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 202125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 72765.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/4. + 1785./8. - 7875.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 212625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 70875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 2835.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 496125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 496125.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 178605.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/4. + 18375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. + 165375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. - 59535.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/2. - 18375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 6615.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/4., 0.);
                      case 3:
                        return sign * RealGradient(-17625.*eta/256. - 9345.*xi/256. + 233625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. - 727125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 291375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 700875.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 1635375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 1635375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 588735.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. - 3195./32. + 52875.*((eta + 1.)*(eta + 1.)/(2.*2.))/64. + 2181375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 874125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 29085.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 5089875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 5089875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 1832355.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64. - 123375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 11655.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. - 2039625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/32. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/32. + 123375.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 44415.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64., 0.);
                      case 4:
                        return sign * RealGradient(375.*eta + 315.*xi/2. - 7875.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 21000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 15750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 47250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 110250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 110250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 39690.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 1005./2. - 4500.*(eta + 1.)*(eta + 1.)/(2.*2.) - 126000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 94500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 840.*(xi + 1.)*(xi + 1.)/(2.*2.) + 294000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 294000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 105840.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 10500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 630.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 79380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 10500.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 3780.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 5:
                        return sign * RealGradient(0., -375.*eta/2. - 375.*xi + 9000.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 47250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 84000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 47250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 31500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 42000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 2125./4. + 2625.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 165375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 294000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 7875.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 220500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 99225.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1750.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 392000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 7000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 176400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 7875.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 6:
                        return sign * RealGradient(0., 14325.*eta/256. + 7275.*xi/128. - 42975.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/16. + 902475.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 100275.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 902475.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 316575.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/32. - 204225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/16. + 174825.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. + 55325./512. - 105525.*(eta + 1.)*(eta + 1.)/(2.*2.)/256. - 6648075.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 738675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. - 6648075.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 152775.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 4288725.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 3671325.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 68075.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 476525.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 4288725.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 16975.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. + 407925.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 58275.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/256. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 152775.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/256.);
                      case 7:
                        return sign * RealGradient(0., -525.*eta/16. - 225.*xi/8. + 1575.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 33075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 14700.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 33075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 17325.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 14175.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 1875./32. + 5775.*((eta + 1.)*(eta + 1.)/(2.*2.))/16. + 363825.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 80850.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 363825.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. + 4725.*((xi + 1.)*(xi + 1.)/(2.*2.))/16. - 297675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 297675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 132300.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 66150.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 4725.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/16. + 297675.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16.);
                      case 8:
                        return sign * RealGradient(0., 3525.*eta/256. + 1275.*xi/128. - 10575.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/16. + 222075.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 24675.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 222075.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 140175.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/32. - 145425.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/16. + 174825.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. + 11725./512. - 46725.*(eta + 1.)*(eta + 1.)/(2.*2.)/256. - 2943675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 327075.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. - 2943675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 26775.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 3053925.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 3671325.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 48475.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 339325.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 3053925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 2975.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. + 407925.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 58275.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/256. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 26775.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/256.);
                      case 9:
                        return sign * RealGradient(0., -75.*eta - 75.*xi + 3600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 33600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 18900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 33600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 575./4. + 1575.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 99225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1575.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 176400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 99225.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 313600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 176400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 10:
                        return sign * RealGradient(75.*eta - 1575.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 66150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 88200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 39690.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 75. - 1800.*(eta + 1.)*(eta + 1.)/(2.*2.) - 50400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 37800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 176400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 235200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 105840.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 6300.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 132300.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 79380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 8400.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 3780.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 11:
                        return sign * RealGradient(-3525.*eta/256. + 46725.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. - 145425.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 58275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 140175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/32. + 981225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 327075.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/16. + 588735.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. - 3525./256. + 10575.*((eta + 1.)*(eta + 1.)/(2.*2.))/32. + 436275.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 174825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. - 3053925.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 1017975.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/16. - 1832355.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64. - 74025.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 1223775.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 407925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/32. + 24675.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/16. - 44415.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64., 0.);
                      case 12:
                        return sign * RealGradient(525.*eta/16. - 5775.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/8. + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 17325.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 121275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 40425.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 72765.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/4. + 525./16. - 1575.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 42525.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 297675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 99225.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 178605.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/4. + 11025.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 66150.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 59535.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/2. - 3675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 6615.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/4., 0.);
                      case 13:
                        return sign * RealGradient(-14325.*eta/256. + 105525.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. - 204225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 58275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 316575.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/32. + 2216025.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 738675.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/16. + 1329615.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/64. - 14325./256. + 42975.*((eta + 1.)*(eta + 1.)/(2.*2.))/32. + 612675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 174825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. - 4288725.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 1429575.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/16. - 2573235.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64. - 300825.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 1223775.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 407925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. + 734265.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/32. + 100275.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/16. - 180495.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/64., 0.);
                      case 14:
                        return sign * RealGradient(375.*eta/2. - 2625.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 5250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 31500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 110250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 147000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 66150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 375./2. - 4500.*(eta + 1.)*(eta + 1.)/(2.*2.) - 63000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 37800.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 220500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 294000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 132300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 15750.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 132300.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 79380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 21000.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 15:
                        return sign * RealGradient(0., -375.*eta - 375.*xi/2. + 9000.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 31500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 42000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 47250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 84000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 47250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 2125./4. + 7875.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 165375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 99225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 2625.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 294000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 165375.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 7000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 392000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1750.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 220500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 16:
                        return sign * RealGradient(0., 17625.*eta/256. + 6375.*xi/256. - 52875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. + 370125.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 123375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 222075.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 700875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/64. - 727125.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 874125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. + 45875./512. - 233625.*(eta + 1.)*(eta + 1.)/(2.*2.)/256. - 4906125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 1635375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 2943675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 44625.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 5089875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 6118875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 242375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 1696625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 3053925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 14875.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. + 2039625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. - 291375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/256. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 26775.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/256.);
                      case 17:
                        return sign * RealGradient(0., -2625.*eta/16. - 1125.*xi/16. + 7875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 55125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 18375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 33075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 86625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/4. + 70875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 70875.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. - 7125./32. + 28875.*((eta + 1.)*(eta + 1.)/(2.*2.))/16. + 606375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 202125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 363825.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. + 7875.*((xi + 1.)*(xi + 1.)/(2.*2.))/16. - 496125.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 496125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 23625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. - 165375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 23625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/16. + 297675.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16.);
                      case 18:
                        return sign * RealGradient(0., 71625.*eta/256. + 36375.*xi/256. - 214875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. + 1504125.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 501375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 902475.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 1582875.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/64. - 1021125.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 874125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. + 203875./512. - 527625.*(eta + 1.)*(eta + 1.)/(2.*2.)/256. - 11080125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 3693375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 6648075.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 254625.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 7147875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 6118875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 340375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/128. - 2382625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 4288725.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. + 84875.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. + 2039625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/32. - 291375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/256. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 152775.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/256.);
                      case 19:
                        return sign * RealGradient(0., -1875.*eta/2. - 1875.*xi/2. + 22500.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 78750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 105000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 47250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 78750.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 105000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 47250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 6875./4. + 13125.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 275625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 367500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 13125.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 367500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 165375.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 8750.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 490000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 8750.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 220500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 7875.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 20:
                        return RealGradient(0., -2925.*eta/4. - 3525.*xi/4. + 21150.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 155475.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 105000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 47250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 74025.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 98700.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 44415.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 11925./8. + 20475.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 1088325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 367500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 51825.*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 362775.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 652995.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 6825.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 490000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 8750.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 220500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 12285.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 7875.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 21:
                        return RealGradient(0., 150.*eta + 2775.*xi/8. - 8325.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 92475.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 84000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 47250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 58275.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 38850.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 34965.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 3775./8. - 1050.*(eta + 1.)*(eta + 1.)/(2.*2.) - 647325.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 294000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 165375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 30825.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 215775.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 388395.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 1400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 392000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 7000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 630.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2.);
                      case 22:
                        return RealGradient(0., -675.*eta/2. - 1125.*xi/8. + 3375.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 6075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. - 23625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 15750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 3375./8. + 4725.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 42525.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. + 2025.*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 25515.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 3150.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2.);
                      case 23:
                        return RealGradient(0., 675.*eta/4. + 225.*xi/2. - 2700.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 6075.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 9450.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 12600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5670.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 2025./8. - 4725.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 42525.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. - 2025.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 14175.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 25515.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4.);
                      case 24:
                        return RealGradient(2925.*eta/4. - 20475.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 20475.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 12285.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 74025.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 362775.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 183750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 66150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 2925./4. - 10575.*(eta + 1.)*(eta + 1.)/(2.*2.) - 148050.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 88830.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 362775.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 367500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 132300.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 51825.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 217665.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 79380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 26250.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 25:
                        return RealGradient(-150.*eta + 2100.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2520.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 58275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 215775.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 147000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 66150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 150. + 8325.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 58275.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 34965.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 215775.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 294000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 132300.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 30825.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 129465.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 79380.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 21000.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 9450.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                      case 26:
                        return RealGradient(675.*eta/2. - 4725.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 9450.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5670.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 23625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 675./2. - 3375.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 23625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 14175.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 2025.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 27:
                        return RealGradient(-675.*eta/4. + 4725.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2835.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 9450.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 14175.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 675./4. + 1350.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 18900.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 11340.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 2025.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 28:
                        return RealGradient(585.*eta/2. - 12285.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 16380.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 12285.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 44415.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 217665.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 110250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 39690.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 585./2. - 4230.*(eta + 1.)*(eta + 1.)/(2.*2.) - 118440.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 88830.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 290220.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 294000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 105840.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 10365.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 217665.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 79380.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 10500.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 3780.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 29:
                        return RealGradient(-60.*eta + 1260.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2520.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 34965.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 129465.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 88200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 39690.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 60. + 1665.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 46620.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 34965.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 172620.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 235200.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 105840.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) - 6165.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 129465.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 79380.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) + 8400.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 3780.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                      case 30:
                        return RealGradient(135.*eta - 2835.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 7560.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5670.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 8505.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 135. - 675.*(eta + 1.)*(eta + 1.)/(2.*2.) - 18900.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 11340.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 405.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 31:
                        return RealGradient(-135.*eta/2. + 2835.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 3780.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 2835.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5670.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8505.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 135./2. + 540.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 15120.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 11340.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 11340.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 405.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)), 0.);
                      case 32:
                        return RealGradient(0., -585.*eta/2. - 705.*xi/4. + 8460.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 31095.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 42000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 44415.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 78960.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 44415.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 3555./8. + 12285.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 652995.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 220500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 99225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 10365.*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 290220.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 652995.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 5460.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 392000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1750.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 220500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 12285.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 33:
                        return RealGradient(0., 60.*eta + 555.*xi/8. - 3330.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 18495.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 33600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 34965.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 31080.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 34965.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 995./8. - 630.*(eta + 1.)*(eta + 1.)/(2.*2.) - 388395.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 176400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 99225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 6165.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 172620.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 388395.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 1120.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 313600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 630.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2.);
                      case 34:
                        return RealGradient(0., -135.*eta - 225.*xi/8. + 1350.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1215.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 12600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 1215./8. + 2835.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 25515.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. + 405.*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 11340.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 25515.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 2520.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2.);
                      case 35:
                        return RealGradient(0., 135.*eta/2. + 45.*xi/2. - 1080.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1215.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 5670.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 10080.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5670.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 675./8. - 2835.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 25515.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. - 405.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 11340.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 25515.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 1260.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4.);
                      case 36:
                        return RealGradient(-72.*eta + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1980.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 5100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1890.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 72. + 1188.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 3060.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3150.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1134.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                      case 37:
                        return RealGradient(-48.*eta + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1980.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 5100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1890.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 48. + 792.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2040.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2100.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 756.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.), 0.);
                      case 38:
                        return RealGradient(30.*eta - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 990.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2550.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 945.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 30. - 495.*(eta + 1.)*(eta + 1.)/(2.*2.) + 1275.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.))/2., 0.);
                      case 39:
                        return RealGradient(0., -72.*eta - 297.*xi + 2376.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 9180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 12600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5670.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1980.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 333. + 120.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 7650.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 2295.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3150.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 40:
                        return RealGradient(0., -48.*eta - 99.*xi + 1584.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 6120.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3780.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1980.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 135. + 120.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 7650.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 765.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1050.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2.);
                      case 41:
                        return RealGradient(0., 30.*eta + 99.*xi/2. - 990.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3825.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 990.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 147./2. - 60.*(eta + 1.)*(eta + 1.)/(2.*2.) - 3825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 765.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4.);
                      case 42:
                        return RealGradient(0., 9.*eta + 108.*xi - 864.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 5400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10080.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 5670.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 720.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 225./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 4500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1350.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2520.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2.);
                      case 43:
                        return RealGradient(0., 6.*eta + 36.*xi - 576.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6720.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3780.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 720.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 81./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 4500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 450.*(xi + 1.)*(xi + 1.)/(2.*2.) + 840.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2.);
                      case 44:
                        return RealGradient(0., -15.*eta/4. - 18.*xi + 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 21. + 15.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 2250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 225.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 420.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4.);
                      case 45:
                        return RealGradient(9.*eta - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 720.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 3000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1890.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 9. - 432.*(eta + 1.)*(eta + 1.)/(2.*2.) + 1800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2520.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1134.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 46:
                        return RealGradient(6.*eta - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 720.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 3000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1890.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.) + 6. - 288.*(eta + 1.)*(eta + 1.)/(2.*2.) + 1200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1680.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 756.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)), 0.);
                      case 47:
                        return RealGradient(-15.*eta/4. + 15.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 1500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 945.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)) - 15./4. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 750.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1050.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.*2.)/2., 0.);
                      case 48:
                        return RealGradient(0., -36.*eta - 81.*xi/2. + 324.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 117./2. + 60.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 135.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 49:
                        return RealGradient(0., 9.*eta + 27.*xi - 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 63./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 135.*(xi + 1.)*(xi + 1.)/(2.*2.)/2.);
                      case 50:
                        return RealGradient(36.*eta - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 36. - 162.*(eta + 1.)*(eta + 1.)/(2.*2.) + 90.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 51:
                        return RealGradient(-9.*eta + 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 9. + 108.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 90.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 52:
                        return RealGradient(24.*eta - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 24. - 108.*(eta + 1.)*(eta + 1.)/(2.*2.) + 60.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 53:
                        return RealGradient(-6.*eta + 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 6. + 72.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 60.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 54:
                        return RealGradient(0., -24.*eta - 27.*xi/2. + 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 63./2. + 60.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 45.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 55:
                        return RealGradient(0., 6.*eta + 9.*xi - 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 27./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.) - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 45.*(xi + 1.)*(xi + 1.)/(2.*2.)/2.);
                      case 56:
                        return RealGradient(0., 0.);
                      case 57:
                        return RealGradient(0., -9.*xi/2. - 5./2. + 15.*((xi + 1.)*(xi + 1.)/(2.*2.))/2.);
                      case 58:
                        return RealGradient(0., 3.*xi + 5./2. - 15.*(xi + 1.)*(xi + 1.)/(2.*2.)/2.);
                      case 59:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(1875.*eta/2. + 1875.*xi/2. - 22500.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 78750.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 105000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 47250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 78750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 105000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 47250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 6875./4. - 13125.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 275625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 367500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 165375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 13125.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 367500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 165375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 8750.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 490000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 8750.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 7875.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 1:
                        return sign * RealGradient(-36375.*eta/256. - 71625.*xi/256. + 214875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. - 1582875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 1021125.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 874125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 1504125.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 501375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/16. - 902475.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. - 203875./512. + 254625.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 11080125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 7147875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 6118875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 527625.*((xi + 1.)*(xi + 1.)/(2.*2.))/256. - 3693375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 6648075.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. - 84875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 2382625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. - 340375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 4288725.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 152775.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/256. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 291375.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/256., 0.);
                      case 2:
                        return sign * RealGradient(1125.*eta/16. + 2625.*xi/16. - 7875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 86625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 70875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 70875.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 55125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. - 18375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 33075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. + 7125./32. - 7875.*(eta + 1.)*(eta + 1.)/(2.*2.)/16. - 606375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 496125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 496125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 28875.*(xi + 1.)*(xi + 1.)/(2.*2.)/16. + 202125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 363825.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. + 2625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 23625.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. + 297675.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/16. - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 23625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16., 0.);
                      case 3:
                        return sign * RealGradient(-6375.*eta/256. - 17625.*xi/256. + 52875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. - 700875.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 727125.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 874125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 370125.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 123375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/16. - 222075.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. - 45875./512. + 44625.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 4906125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 5089875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 6118875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 233625.*((xi + 1.)*(xi + 1.)/(2.*2.))/256. - 1635375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 2943675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. - 14875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 1696625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. - 242375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 3053925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 26775.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/256. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 291375.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/256., 0.);
                      case 4:
                        return sign * RealGradient(375.*eta/2. + 375.*xi - 9000.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 47250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 84000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 47250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 31500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 42000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 2125./4. - 2625.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 165375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 294000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 165375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 7875.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 220500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99225.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1750.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 392000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 7000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 5:
                        return sign * RealGradient(0., -375.*xi/2. + 2625.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 31500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 110250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 147000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 66150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 5250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 3150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 375./2. + 63000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 294000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 132300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 4500.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 37800.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 132300.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 79380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 15750.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 21000.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 9450.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 6:
                        return sign * RealGradient(0., 14325.*xi/256. - 105525.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. + 316575.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 2216025.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 738675.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16. - 1329615.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 204225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 58275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 14325./256. - 612675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 4288725.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 1429575.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16. + 2573235.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64. - 42975.*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 174825.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/16. - 1223775.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 407925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/32. + 300825.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 100275.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16. + 180495.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64.);
                      case 7:
                        return sign * RealGradient(0., -525.*xi/16. + 5775.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/8. - 17325.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 121275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 40425.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 72765.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/4. - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/8. + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 525./16. + 42525.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 297675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 99225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 178605.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/4. + 1575.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 66150.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 59535.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/2. - 11025.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 3675.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 6615.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/4.);
                      case 8:
                        return sign * RealGradient(0., 3525.*xi/256. - 46725.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. + 140175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 981225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 327075.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/16. - 588735.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 145425.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 58275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 3525./256. - 436275.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 3053925.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 1017975.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16. + 1832355.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64. - 10575.*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 174825.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/16. - 1223775.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 407925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/32. + 74025.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 24675.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16. + 44415.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64.);
                      case 9:
                        return sign * RealGradient(0., -75.*xi + 1575.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 66150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 88200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 3150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 75. + 50400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 235200.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 105840.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 1800.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 37800.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 132300.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 79380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 6300.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 8400.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 3780.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 10:
                        return sign * RealGradient(75.*eta + 75.*xi - 3600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 18900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 33600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 33600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 575./4. - 1575.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 176400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 99225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1575.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 176400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99225.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 313600.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 11:
                        return sign * RealGradient(-1275.*eta/128. - 3525.*xi/256. + 10575.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/16. - 140175.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 145425.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 174825.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. - 222075.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 24675.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 222075.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. - 11725./512. + 26775.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 2943675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 3053925.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 3671325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 46725.*((xi + 1.)*(xi + 1.)/(2.*2.))/256. - 327075.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 2943675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. - 2975.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/16. + 339325.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 407925.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 48475.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 3053925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 26775.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/256. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 58275.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/256., 0.);
                      case 12:
                        return sign * RealGradient(225.*eta/8. + 525.*xi/16. - 1575.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 17325.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 14175.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 33075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. - 14700.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 33075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. + 1875./32. - 4725.*(eta + 1.)*(eta + 1.)/(2.*2.)/16. - 363825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 297675.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 297675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 5775.*(xi + 1.)*(xi + 1.)/(2.*2.)/16. + 80850.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 363825.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 132300.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 66150.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. + 297675.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/16. - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 4725.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/16., 0.);
                      case 13:
                        return sign * RealGradient(-7275.*eta/128. - 14325.*xi/256. + 42975.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/16. - 316575.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/32. + 204225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 174825.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/32. - 902475.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 100275.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 902475.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. - 55325./512. + 152775.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 6648075.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 4288725.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 3671325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 105525.*((xi + 1.)*(xi + 1.)/(2.*2.))/256. - 738675.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 6648075.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. - 16975.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/16. + 476525.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 407925.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 68075.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/128. - 4288725.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/64. + 152775.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/256. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 58275.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/256., 0.);
                      case 14:
                        return sign * RealGradient(375.*eta + 375.*xi/2. - 9000.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 31500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 42000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 18900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 47250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 84000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 47250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 2125./4. - 7875.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 165375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 220500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 99225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 2625.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 294000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 165375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 392000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 1750.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 7875.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 15:
                        return sign * RealGradient(0., -315.*eta/2. - 375.*xi + 7875.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 47250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 110250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 110250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 21000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 15750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1005./2. + 840.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 126000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 294000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 294000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 105840.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 4500.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 94500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 630.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 79380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 10500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 10500.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 3780.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 16:
                        return sign * RealGradient(0., 9345.*eta/256. + 17625.*xi/256. - 233625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. + 700875.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 1635375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 1635375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. - 588735.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 727125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 291375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 3195./32. - 29085.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 2181375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 5089875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 5089875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. + 1832355.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64. - 52875.*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 874125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. + 11655.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/32. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/32. + 123375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 123375.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. + 44415.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64.);
                      case 17:
                        return sign * RealGradient(0., -1155.*eta/16. - 2625.*xi/16. + 28875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/8. - 86625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 202125.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 202125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. + 72765.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/4. - 70875.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/8. + 23625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. - 1785./8. + 2835.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 212625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 496125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 496125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. - 178605.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/4. + 7875.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 70875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/4. + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. + 59535.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/2. - 18375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 18375.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. - 6615.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/4.);
                      case 18:
                        return sign * RealGradient(0., 21105.*eta/256. + 71625.*xi/256. - 527625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/128. + 1582875.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 3693375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/64. + 3693375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/64. - 1329615.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/64. + 1021125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 291375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. + 10875./32. - 40845.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 3063375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 7147875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 7147875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. + 2573235.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64. - 214875.*(xi + 1.)*(xi + 1.)/(2.*2.)/64. + 874125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. + 11655.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/64. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/32. - 734265.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/32. + 501375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 501375.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. + 180495.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/64.);
                      case 19:
                        return sign * RealGradient(0., -525.*eta/2. - 1875.*xi/2. + 13125.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 78750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 183750.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 183750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 66150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 26250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 15750.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1125. + 1050.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 157500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 367500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 367500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 132300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 11250.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 94500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 630.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 79380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 26250.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 26250.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 9450.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 20:
                        return RealGradient(0., -2925.*xi/4. + 20475.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 74025.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 362775.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 183750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 66150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 20475.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 12285.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2925./4. + 148050.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 362775.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 367500.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 132300.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 10575.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 88830.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 217665.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 79380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 51825.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 26250.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 9450.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 21:
                        return RealGradient(0., 150.*xi - 2100.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 58275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 215775.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 147000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 66150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2520.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 150. - 58275.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 215775.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 294000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 132300.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 8325.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 34965.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 129465.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 79380.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 30825.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 21000.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 9450.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                      case 22:
                        return RealGradient(0., -675.*xi/2. + 4725.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 23625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 9450.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 5670.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 675./2. + 23625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 14175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3375.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2025.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2.);
                      case 23:
                        return RealGradient(0., 675.*xi/4. - 4725.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 9450.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 14175.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2835.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 675./4. - 18900.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1350.*(xi + 1.)*(xi + 1.)/(2.*2.) + 11340.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2025.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2.);
                      case 24:
                        return RealGradient(3525.*eta/4. + 2925.*xi/4. - 21150.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 74025.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 98700.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 44415.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 155475.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 105000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 47250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 11925./8. - 51825.*(eta + 1.)*(eta + 1.)/(2.*2.)/8. - 1088325.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 362775.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 652995.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 20475.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 367500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 165375.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 8750.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 490000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 6825.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 220500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 7875.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 12285.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4., 0.);
                      case 25:
                        return RealGradient(-2775.*eta/8. - 150.*xi + 8325.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 58275.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 38850.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 34965.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 92475.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. + 84000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 47250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 3775./8. + 30825.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 647325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 215775.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 388395.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 1050.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 294000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 165375.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 7000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 392000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 220500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 630.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)), 0.);
                      case 26:
                        return RealGradient(1125.*eta/8. + 675.*xi/2. - 3375.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 23625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 15750.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 6075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 3375./8. - 2025.*(eta + 1.)*(eta + 1.)/(2.*2.)/8. - 42525.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 14175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 25515.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 4725.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 3150.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 27:
                        return RealGradient(-225.*eta/2. - 675.*xi/4. + 2700.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 9450.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 12600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5670.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 6075.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 2025./8. + 2025.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 42525.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 14175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 25515.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 4725.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4., 0.);
                      case 28:
                        return RealGradient(705.*eta/4. + 585.*xi/2. - 8460.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 44415.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 78960.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 44415.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 31095.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 42000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 3555./8. - 10365.*(eta + 1.)*(eta + 1.)/(2.*2.)/8. - 652995.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 290220.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 652995.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 12285.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 220500.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99225.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1750.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 392000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 220500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 5460.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 176400.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 12285.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4., 0.);
                      case 29:
                        return RealGradient(-555.*eta/8. - 60.*xi + 3330.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 34965.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 31080.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 34965.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 18495.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 33600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 18900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 995./8. + 6165.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 388395.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 172620.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 388395.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 630.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 176400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 99225.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 313600.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 176400.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1120.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 176400.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 630.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)), 0.);
                      case 30:
                        return RealGradient(225.*eta/8. + 135.*xi - 1350.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 12600.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 1215.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 1215./8. - 405.*(eta + 1.)*(eta + 1.)/(2.*2.)/8. - 25515.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 11340.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 25515.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4. - 2835.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 2520.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 31:
                        return RealGradient(-45.*eta/2. - 135.*xi/2. + 1080.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 5670.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 10080.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5670.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 1215.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 675./8. + 405.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 25515.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 11340.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 25515.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 2835.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 1260.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4., 0.);
                      case 32:
                        return RealGradient(0., -585.*xi/2. + 12285.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 44415.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 217665.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 110250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 39690.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 16380.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 12285.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 585./2. + 118440.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 290220.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 294000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 105840.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 4230.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 88830.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 217665.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 220500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 79380.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 10365.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 10500.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 3780.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 33:
                        return RealGradient(0., 60.*xi - 1260.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 34965.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 129465.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 88200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 39690.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 3360.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2520.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 60. - 46620.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 172620.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 235200.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 105840.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 1665.*(xi + 1.)*(xi + 1.)/(2.*2.) + 34965.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 129465.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 176400.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 79380.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 6165.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 8400.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 3780.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                      case 34:
                        return RealGradient(0., -135.*xi + 2835.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 14175.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 8505.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 7560.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 5670.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 135. + 18900.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 11340.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 675.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 14175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 405.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 35:
                        return RealGradient(0., 135.*xi/2. - 2835.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 5670.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8505.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 3780.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 2835.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 135./2. - 15120.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 11340.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 540.*(xi + 1.)*(xi + 1.)/(2.*2.) + 11340.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 405.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 36:
                        return RealGradient(-297.*eta - 72.*xi + 2376.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1980.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 9180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 12600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 5670.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 333. + 2295.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 7650.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 120.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 3150.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2., 0.);
                      case 37:
                        return RealGradient(-99.*eta - 48.*xi + 1584.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1980.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 6120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 3780.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 135. + 765.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 7650.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 120.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10500.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 1050.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2., 0.);
                      case 38:
                        return RealGradient(99.*eta/2. + 30.*xi - 990.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 990.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 3825.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 5250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 147./2. - 765.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 3825.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.) + 5250.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4., 0.);
                      case 39:
                        return RealGradient(0., -72.*xi + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1980.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1890.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 72. + 1188.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3060.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1134.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 40:
                        return RealGradient(0., -48.*xi + 240.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 1980.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 5250.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1890.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 48. + 792.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2040.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2100.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 756.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.));
                      case 41:
                        return RealGradient(0., 30.*xi - 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 990.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2550.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 945.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 30. - 495.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1275.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.))/2.);
                      case 42:
                        return RealGradient(0., 9.*xi - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 720.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1890.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 9. - 432.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1800.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2520.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 1134.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                      case 43:
                        return RealGradient(0., 6.*xi - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 720.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4200.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 1890.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.) + 6. - 288.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1680.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 756.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)));
                      case 44:
                        return RealGradient(0., -15.*xi/4. + 15.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 1500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) + 945.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)) - 15./4. + 180.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 750.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1050.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.*2.)/2.);
                      case 45:
                        return RealGradient(108.*eta + 9.*xi - 864.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 720.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 5400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 10080.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 5670.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 225./2. - 1350.*(eta + 1.)*(eta + 1.)/(2.*2.) - 4500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 15.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 2520.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 46:
                        return RealGradient(36.*eta + 6.*xi - 576.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 720.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 3600.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 6720.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3780.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 81./2. - 450.*(eta + 1.)*(eta + 1.)/(2.*2.) - 4500.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 15.*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 840.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2., 0.);
                      case 47:
                        return RealGradient(-18.*eta - 15.*xi/4. + 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 360.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 2250.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 4200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 21. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 2250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 4200.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. - 420.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4., 0.);
                      case 48:
                        return RealGradient(0., -36.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 36. + 162.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 90.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 49:
                        return RealGradient(0., 9.*xi - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 9. - 108.*(xi + 1.)*(xi + 1.)/(2.*2.) + 90.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 50:
                        return RealGradient(81.*eta/2. + 36.*xi - 324.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 270.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 117./2. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.), 0.);
                      case 51:
                        return RealGradient(-27.*eta - 9.*xi + 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 270.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 63./2. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.)), 0.);
                      case 52:
                        return RealGradient(27.*eta/2. + 24.*xi - 216.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 270.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 180.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 63./2. - 45.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 60.*(xi + 1.)*(xi + 1.)/(2.*2.), 0.);
                      case 53:
                        return RealGradient(-9.*eta - 6.*xi + 144.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 180.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 180.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 27./2. + 45.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 225.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 15.*((xi + 1.)*(xi + 1.)/(2.*2.)), 0.);
                      case 54:
                        return RealGradient(0., -24.*xi + 120.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 270.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 24. + 108.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 60.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 55:
                        return RealGradient(0., 6.*xi - 30.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 180.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 6. - 72.*(xi + 1.)*(xi + 1.)/(2.*2.) + 60.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 56:
                        return RealGradient(-9.*eta/2. - 5./2. + 15.*((eta + 1.)*(eta + 1.)/(2.*2.))/2., 0.);
                      case 57:
                        return RealGradient(0., 0.);
                      case 58:
                        return RealGradient(0., 0.);
                      case 59:
                        return RealGradient(3.*eta + 5./2. - 15.*(eta + 1.)*(eta + 1.)/(2.*2.)/2., 0.);
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
                        return sign * RealGradient(-12600.*eta*xi + 2800.*eta + 15120.*eta*(xi*xi) - 4200.*eta*xi*xi*xi + 2100.*xi + 22680.*xi*(eta*eta) - 12600.*xi*eta*eta*eta - 300. - 8400.*eta*eta - 12600.*eta*eta*xi*xi - 4200.*xi*xi + 10080.*(eta*eta*eta) + 2520.*(xi*xi*xi) - 4200.*eta*eta*eta*eta, 8400.*eta*xi - 700.*eta - 22680.*eta*xi*xi + 16800.*eta*(xi*xi*xi) - 1400.*xi - 15120.*xi*eta*eta + 8400.*xi*(eta*eta*eta) + 75. + 2100.*(eta*eta) + 18900.*(eta*eta)*(xi*xi) + 6300.*(xi*xi) - 2520.*eta*eta*eta - 10080.*xi*xi*xi + 1050.*(eta*eta*eta*eta) + 5250.*(xi*xi*xi*xi));
                      case 1:
                        return sign * RealGradient(136185.*eta*xi/32. - 12145.*eta/16. - 20475.*eta*xi*xi/4. + 19425.*eta*(xi*xi*xi)/16. - 21105.*xi/32. - 240975.*xi*eta*eta/32. + 127575.*xi*(eta*eta*eta)/32. + 2865./32. + 7245.*(eta*eta)/4. + 146475.*(eta*eta)*(xi*xi)/32. + 40845.*(xi*xi)/32. - 25515.*eta*eta*eta/16. - 11655.*xi*xi*xi/16. + 14175.*(eta*eta*eta*eta)/32., -39375.*eta*xi/16. + 3045.*eta/32. + 252315.*eta*(xi*xi)/32. - 48825.*eta*xi*xi*xi/8. + 1715.*xi/4. + 25515.*xi*(eta*eta)/8. - 14175.*xi*eta*eta*eta/16. - 2395./128. + 2835.*(eta*eta)/64. - 382725.*eta*eta*xi*xi/64. - 125265.*xi*xi/64. - 14175.*eta*eta*eta/32. + 48615.*(xi*xi*xi)/16. + 42525.*(eta*eta*eta*eta)/128. - 97125.*xi*xi*xi*xi/64.);
                      case 2:
                        return sign * RealGradient(-3255.*eta*xi/2. + 3780.*eta*(xi*xi) - 1575.*eta*xi*xi*xi + 1155.*xi/2. + 945.*xi*(eta*eta)/2. + 1575.*xi*(eta*eta*eta)/2. - 105./2. + 840.*(eta*eta) - 4725.*eta*eta*xi*xi/2. - 2835.*xi*xi/2. - 1575.*eta*eta*eta + 945.*(xi*xi*xi) + 1575.*(eta*eta*eta*eta)/2., 105.*eta*xi + 35.*eta/2. - 4725.*eta*xi*xi/2. + 3150.*eta*(xi*xi*xi) - 315.*xi + 1890.*xi*(eta*eta) - 1575.*xi*eta*eta*eta + 85./8. - 525.*eta*eta/4. - 4725.*eta*eta*xi*xi/4. + 7455.*(xi*xi)/4. + 105.*(eta*eta*eta)/2. - 3465.*xi*xi*xi + 525.*(eta*eta*eta*eta)/8. + 7875.*(xi*xi*xi*xi)/4.);
                      case 3:
                        return sign * RealGradient(-16695.*eta*xi/32. + 1855.*eta/16. - 2835.*eta*xi*xi/4. + 19425.*eta*(xi*xi*xi)/16. - 9345.*xi/32. + 76545.*xi*(eta*eta)/32. - 48825.*xi*eta*eta*eta/32. + 705./32. - 1155.*eta*eta/4. - 29925.*eta*eta*xi*xi/32. + 29085.*(xi*xi)/32. - 315.*eta*eta*eta/16. - 11655.*xi*xi*xi/16. + 5775.*(eta*eta*eta*eta)/32., 9345.*eta*xi/16. - 875.*eta/32. - 65205.*eta*xi*xi/32. + 9975.*eta*(xi*xi*xi)/8. + 595.*xi/4. - 4725.*xi*eta*eta/8. - 5775.*xi*eta*eta*eta/16. - 555./128. + 1155.*(eta*eta)/64. + 146475.*(eta*eta)*(xi*xi)/64. - 66465.*xi*xi/64. + 945.*(eta*eta*eta)/32. + 36855.*(xi*xi*xi)/16. + 525.*(eta*eta*eta*eta)/128. - 97125.*xi*xi*xi*xi/64.);
                      case 4:
                        return sign * RealGradient(-1680.*eta*xi + 140.*eta + 5040.*eta*(xi*xi) - 4200.*eta*xi*xi*xi + 1260.*xi - 120. - 3360.*xi*xi + 2520.*(xi*xi*xi), -700.*xi + 25. + 4200.*(xi*xi) - 8400.*xi*xi*xi + 5250.*(xi*xi*xi*xi));
                      case 5:
                        return sign * RealGradient(-1680.*eta*xi + 140.*eta + 5040.*eta*(xi*xi) - 4200.*eta*xi*xi*xi, -700.*xi + 25. + 4200.*(xi*xi) - 8400.*xi*xi*xi + 5250.*(xi*xi*xi*xi));
                      case 6:
                        return sign * RealGradient(-2835.*eta*xi/2. + 735.*eta/4. + 42525.*eta*(xi*xi)/16. - 42525.*eta*xi*xi*xi/32. + 76545.*xi*(eta*eta)/16. - 42525.*xi*eta*eta*eta/16. - 1575.*eta*eta/2. - 42525.*eta*eta*xi*xi/8. + 12285.*(eta*eta*eta)/16. - 1575.*eta*eta*eta*eta/8., 2205.*eta*xi - 455.*eta/4. - 127575.*eta*xi*xi/16. + 14175.*eta*(xi*xi*xi)/2. - 1785.*xi/4. - 31185.*xi*eta*eta/16. + 1575.*xi*(eta*eta*eta)/4. + 20. + 525.*(eta*eta)/4. + 127575.*(eta*eta)*(xi*xi)/32. + 8505.*(xi*xi)/4. - 735.*eta*eta*eta/16. - 53865.*xi*xi*xi/16. + 525.*(eta*eta*eta*eta)/128. + 212625.*(xi*xi*xi*xi)/128.);
                      case 7:
                        return sign * RealGradient(-840.*eta*xi + 175.*eta + 945.*eta*(xi*xi) - 525.*eta*xi*xi*xi/2. + 4725.*xi*(eta*eta) - 4725.*xi*eta*eta*eta - 1260.*eta*eta - 3150.*eta*eta*xi*xi + 2205.*(eta*eta*eta) - 1050.*eta*eta*eta*eta, 2520.*eta*xi - 175.*eta - 6615.*eta*xi*xi + 4200.*eta*(xi*xi*xi) - 245.*xi - 4725.*xi*eta*eta + 2100.*xi*(eta*eta*eta) + 15. + 420.*(eta*eta) + 14175.*(eta*eta)*(xi*xi)/2. + 840.*(xi*xi) - 315.*eta*eta*eta - 945.*xi*xi*xi + 525.*(eta*eta*eta*eta)/8. + 2625.*(xi*xi*xi*xi)/8.);
                      case 8:
                        return sign * RealGradient(-525.*eta*xi/2. + 455.*eta/4. + 2205.*eta*(xi*xi)/16. - 525.*eta*xi*xi*xi/32. + 31185.*xi*(eta*eta)/16. - 42525.*xi*eta*eta*eta/16. - 2205.*eta*eta/2. - 4725.*eta*eta*xi*xi/8. + 42525.*(eta*eta*eta)/16. - 14175.*eta*eta*eta*eta/8., 1575.*eta*xi - 735.*eta/4. - 36855.*eta*xi*xi/16. + 1575.*eta*(xi*xi*xi)/2. - 385.*xi/4. - 76545.*xi*eta*eta/16. + 14175.*xi*(eta*eta*eta)/4. + 10. + 2835.*(eta*eta)/4. + 127575.*(eta*eta)*(xi*xi)/32. + 735.*(xi*xi)/4. - 14175.*eta*eta*eta/16. - 1785.*xi*xi*xi/16. + 42525.*(eta*eta*eta*eta)/128. + 2625.*(xi*xi*xi*xi)/128.);
                      case 9:
                        return sign * RealGradient(0., -140.*eta + 5. + 840.*(eta*eta) - 1680.*eta*eta*eta + 1050.*(eta*eta*eta*eta));
                      case 10:
                        return sign * RealGradient(0., -140.*eta + 5. + 840.*(eta*eta) - 1680.*eta*eta*eta + 1050.*(eta*eta*eta*eta));
                      case 11:
                        return sign * RealGradient(-1155.*eta*xi/32. + 875.*eta/32. - 2835.*eta*xi*xi/32. - 525.*eta*xi*xi*xi/32. + 4725.*xi*(eta*eta)/8. - 48825.*xi*eta*eta*eta/32. - 9345.*eta*eta/32. + 17325.*(eta*eta)*(xi*xi)/32. + 21735.*(eta*eta*eta)/32. - 9975.*eta*eta*eta*eta/32., 1155.*eta*xi/2. - 1855.*eta/16. + 945.*eta*(xi*xi)/16. - 5775.*eta*xi*xi*xi/8. - 595.*xi/32. - 76545.*xi*eta*eta/32. + 9975.*xi*(eta*eta*eta)/16. + 825./128. + 16695.*(eta*eta)/64. + 146475.*(eta*eta)*(xi*xi)/64. - 1785.*xi*xi/64. + 945.*(eta*eta*eta)/4. + 945.*(xi*xi*xi)/32. - 19425.*eta*eta*eta*eta/64. + 2625.*(xi*xi*xi*xi)/128.);
                      case 12:
                        return sign * RealGradient(525.*eta*xi/2. - 35.*eta/2. - 315.*eta*xi*xi/2. - 525.*eta*xi*xi*xi/2. - 1890.*xi*eta*eta + 1575.*xi*(eta*eta*eta)/2. - 105.*eta*eta/2. + 4725.*(eta*eta)*(xi*xi)/2. + 1575.*(eta*eta*eta)/2. - 1575.*eta*eta*eta*eta/2., -1680.*eta*xi + 4725.*eta*(xi*xi) - 3150.*eta*xi*xi*xi + 175.*xi/2. - 945.*xi*eta*eta/2. + 1575.*xi*(eta*eta*eta) - 55./8. + 3255.*(eta*eta)/4. - 4725.*eta*eta*xi*xi/4. - 105.*xi*xi/4. - 1260.*eta*eta*eta - 735.*xi*xi*xi/2. + 1575.*(eta*eta*eta*eta)/4. + 2625.*(xi*xi*xi*xi)/8.);
                      case 13:
                        return sign * RealGradient(-2835.*eta*xi/32. - 3045.*eta/32. + 42525.*eta*(xi*xi)/32. - 42525.*eta*xi*xi*xi/32. - 25515.*xi*eta*eta/8. + 127575.*xi*(eta*eta*eta)/32. + 39375.*(eta*eta)/32. + 42525.*(eta*eta)*(xi*xi)/32. - 84105.*eta*eta*eta/32. + 48825.*(eta*eta*eta*eta)/32., -7245.*eta*xi/2. + 12145.*eta/16. + 76545.*eta*(xi*xi)/16. - 14175.*eta*xi*xi*xi/8. - 11235.*xi/32. + 240975.*xi*(eta*eta)/32. - 48825.*xi*eta*eta*eta/16. - 695./128. - 136185.*eta*eta/64. - 382725.*eta*eta*xi*xi/64. + 127575.*(xi*xi)/64. + 6825.*(eta*eta*eta)/4. - 104895.*xi*xi*xi/32. - 19425.*eta*eta*eta*eta/64. + 212625.*(xi*xi*xi*xi)/128.);
                      case 14:
                        return sign * RealGradient(-4200.*eta*xi + 700.*eta + 7560.*eta*(xi*xi) - 4200.*eta*xi*xi*xi + 15120.*xi*(eta*eta) - 12600.*xi*eta*eta*eta - 4200.*eta*eta - 12600.*eta*eta*xi*xi + 7560.*(eta*eta*eta) - 4200.*eta*eta*eta*eta, 16800.*eta*xi - 2800.*eta - 30240.*eta*xi*xi + 16800.*eta*(xi*xi*xi) - 3500.*xi - 22680.*xi*eta*eta + 8400.*xi*(eta*eta*eta) + 375. + 6300.*(eta*eta) + 18900.*(eta*eta)*(xi*xi) + 10500.*(xi*xi) - 5040.*eta*eta*eta - 12600.*xi*xi*xi + 1050.*(eta*eta*eta*eta) + 5250.*(xi*xi*xi*xi));
                      case 15:
                        return RealGradient(70560.*eta*xi - 19600.*eta - 70560.*eta*xi*xi + 16800.*eta*(xi*xi*xi) - 211680.*xi*eta*eta + 151200.*xi*(eta*eta*eta) + 94080.*(eta*eta) + 100800.*(eta*eta)*(xi*xi) - 141120.*eta*eta*eta + 67200.*(eta*eta*eta*eta), -94080.*eta*xi + 9800.*eta + 211680.*eta*(xi*xi) - 134400.*eta*xi*xi*xi + 9800.*xi + 211680.*xi*(eta*eta) - 134400.*xi*eta*eta*eta - 700. - 35280.*eta*eta - 226800.*eta*eta*xi*xi - 35280.*xi*xi + 47040.*(eta*eta*eta) + 47040.*(xi*xi*xi) - 21000.*eta*eta*eta*eta - 21000.*xi*xi*xi*xi);
                      case 16:
                        return RealGradient(-70560.*eta*xi + 9800.*eta + 141120.*eta*(xi*xi) - 84000.*eta*xi*xi*xi + 211680.*xi*(eta*eta) - 151200.*xi*eta*eta*eta - 47040.*eta*eta - 201600.*eta*eta*xi*xi + 70560.*(eta*eta*eta) - 33600.*eta*eta*eta*eta, 188160.*eta*xi - 19600.*eta - 423360.*eta*xi*xi + 268800.*eta*(xi*xi*xi) - 49000.*xi - 211680.*xi*eta*eta + 67200.*xi*(eta*eta*eta) + 3500. + 35280.*(eta*eta) + 226800.*(eta*eta)*(xi*xi) + 176400.*(xi*xi) - 23520.*eta*eta*eta - 235200.*xi*xi*xi + 4200.*(eta*eta*eta*eta) + 105000.*(xi*xi*xi*xi));
                      case 17:
                        return RealGradient(-60480.*eta*xi + 6440.*eta + 131040.*eta*(xi*xi) - 67200.*eta*xi*xi*xi + 60480.*xi*(eta*eta) - 6720.*eta*eta - 100800.*eta*eta*xi*xi, 26880.*eta*xi - 1120.*eta - 120960.*eta*xi*xi + 134400.*eta*(xi*xi*xi) - 14560.*xi + 560. + 80640.*(xi*xi) - 147840.*xi*xi*xi + 84000.*(xi*xi*xi*xi));
                      case 18:
                        return RealGradient(40320.*eta*xi - 3640.*eta - 110880.*eta*xi*xi + 84000.*eta*(xi*xi*xi) - 30240.*xi*eta*eta + 3360.*(eta*eta) + 50400.*(eta*eta)*(xi*xi), -13440.*eta*xi + 560.*eta + 60480.*eta*(xi*xi) - 67200.*eta*xi*xi*xi + 18200.*xi - 700. - 100800.*xi*xi + 184800.*(xi*xi*xi) - 105000.*xi*xi*xi*xi);
                      case 19:
                        return RealGradient(560.*eta - 6720.*eta*eta + 20160.*(eta*eta*eta) - 16800.*eta*eta*eta*eta, 6720.*eta*xi - 3640.*eta - 280.*xi - 30240.*xi*eta*eta + 33600.*xi*(eta*eta*eta) + 140. + 20160.*(eta*eta) - 36960.*eta*eta*eta + 21000.*(eta*eta*eta*eta));
                      case 20:
                        return RealGradient(-1120.*eta + 13440.*(eta*eta) - 40320.*eta*eta*eta + 33600.*(eta*eta*eta*eta), -13440.*eta*xi + 6440.*eta + 560.*xi + 60480.*xi*(eta*eta) - 67200.*xi*eta*eta*eta - 280. - 30240.*eta*eta + 43680.*(eta*eta*eta) - 16800.*eta*eta*eta*eta);
                      case 21:
                        return RealGradient(-35840.*eta*xi + 17920.*eta/3. + 156800.*eta*(xi*xi)/3. - 179200.*eta*xi*xi*xi/9. + 125440.*xi*(eta*eta) - 89600.*xi*eta*eta*eta - 29120.*eta*eta - 291200.*eta*eta*xi*xi/3. + 116480.*(eta*eta*eta)/3. - 140000.*eta*eta*eta*eta/9., 56000.*eta*xi - 10360.*eta/3. - 170240.*eta*xi*xi + 1164800.*eta*(xi*xi*xi)/9. - 24920.*xi/3. - 81760.*xi*eta*eta + 280000.*xi*(eta*eta*eta)/9. + 420. + 6720.*(eta*eta) + 134400.*(eta*eta)*(xi*xi) + 35840.*(xi*xi) - 39200.*eta*eta*eta/9. - 474880.*xi*xi*xi/9. + 7000.*(eta*eta*eta*eta)/9. + 224000.*(xi*xi*xi*xi)/9.);
                      case 22:
                        return RealGradient(29120.*eta*xi - 12880.*eta/3. - 152320.*eta*xi*xi/3. + 224000.*eta*(xi*xi*xi)/9. - 109760.*xi*eta*eta + 78400.*xi*(eta*eta*eta) + 22400.*(eta*eta) + 313600.*(eta*eta)*(xi*xi)/3. - 91840.*eta*eta*eta/3. + 112000.*(eta*eta*eta*eta)/9., -62720.*eta*xi + 12040.*eta/3. + 185920.*eta*(xi*xi) - 1254400.*eta*xi*xi*xi/9. + 32480.*xi/3. + 76160.*xi*(eta*eta) - 224000.*xi*eta*eta*eta/9. - 560. - 6720.*eta*eta - 117600.*eta*eta*xi*xi - 45920.*xi*xi + 34720.*(eta*eta*eta)/9. + 600320.*(xi*xi*xi)/9. - 5600.*eta*eta*eta*eta/9. - 280000.*xi*xi*xi*xi/9.);
                      case 23:
                        return RealGradient(-13440.*eta*xi + 12040.*eta/3. + 34720.*eta*(xi*xi)/3. - 22400.*eta*xi*xi*xi/9. + 76160.*xi*(eta*eta) - 78400.*xi*eta*eta*eta - 31360.*eta*eta - 112000.*eta*eta*xi*xi/3. + 185920.*(eta*eta*eta)/3. - 313600.*eta*eta*eta*eta/9., 44800.*eta*xi - 12880.*eta/3. - 91840.*eta*xi*xi + 448000.*eta*(xi*xi*xi)/9. - 10640.*xi/3. - 109760.*xi*eta*eta + 627200.*xi*(eta*eta*eta)/9. + 280. + 14560.*(eta*eta) + 117600.*(eta*eta)*(xi*xi) + 10080.*(xi*xi) - 152320.*eta*eta*eta/9. - 89600.*xi*xi*xi/9. + 56000.*(eta*eta*eta*eta)/9. + 28000.*(xi*xi*xi*xi)/9.);
                      case 24:
                        return RealGradient(13440.*eta*xi - 10360.*eta/3. - 39200.*eta*xi*xi/3. + 28000.*eta*(xi*xi*xi)/9. - 81760.*xi*eta*eta + 89600.*xi*(eta*eta*eta) + 28000.*(eta*eta) + 140000.*(eta*eta)*(xi*xi)/3. - 170240.*eta*eta*eta/3. + 291200.*(eta*eta*eta*eta)/9., -58240.*eta*xi + 17920.*eta/3. + 116480.*eta*(xi*xi) - 560000.*eta*xi*xi*xi/9. + 14840.*xi/3. + 125440.*xi*(eta*eta) - 582400.*xi*eta*eta*eta/9. - 420. - 17920.*eta*eta - 134400.*eta*eta*xi*xi - 13440.*xi*xi + 156800.*(eta*eta*eta)/9. + 115360.*(xi*xi*xi)/9. - 44800.*eta*eta*eta*eta/9. - 35000.*xi*xi*xi*xi/9.);
                      case 25:
                        return RealGradient(11200.*eta*xi - 12880.*eta/9. - 49280.*eta*xi*xi/3. + 44800.*eta*(xi*xi*xi)/9. - 2240.*xi*eta*eta - 11200.*xi*eta*eta*eta - 15680.*eta*eta/3. + 22400.*(eta*eta)*(xi*xi)/3. + 51520.*(eta*eta*eta)/3. - 95200.*eta*eta*eta*eta/9., 15680.*eta*xi/3. - 10360.*eta/9. + 2240.*eta*(xi*xi) - 89600.*eta*xi*xi*xi/9. + 6440.*xi/9. - 25760.*xi*eta*eta + 190400.*xi*(eta*eta*eta)/9. + 140./9. + 5600.*(eta*eta) + 16800.*(eta*eta)*(xi*xi) - 5600.*xi*xi - 75040.*eta*eta*eta/9. + 98560.*(xi*xi*xi)/9. + 35000.*(eta*eta*eta*eta)/9. - 56000.*xi*xi*xi*xi/9.);
                      case 26:
                        return RealGradient(-2240.*eta*xi - 9520.*eta/9. + 71680.*eta*(xi*xi)/3. - 224000.*eta*xi*xi*xi/9. - 51520.*xi*eta*eta + 56000.*xi*(eta*eta*eta) + 44800.*(eta*eta)/3. + 89600.*(eta*eta)*(xi*xi)/3. - 82880.*eta*eta*eta/3. + 123200.*(eta*eta*eta*eta)/9., -138880.*eta*xi/3. + 51800.*eta/9. + 82880.*eta*(xi*xi) - 358400.*eta*xi*xi*xi/9. - 48160.*xi/9. + 80640.*xi*(eta*eta) - 246400.*xi*eta*eta*eta/9. + 560./9. - 13440.*eta*eta - 84000.*eta*eta*xi*xi + 32480.*(xi*xi) + 79520.*(eta*eta*eta)/9. - 519680.*xi*xi*xi/9. - 11200.*eta*eta*eta*eta/9. + 280000.*(xi*xi*xi*xi)/9.);
                      case 27:
                        return RealGradient(-1120.*eta*xi + 2800.*eta/9. - 1120.*eta*xi*xi/3. + 5600.*eta*(xi*xi*xi)/9. + 12320.*xi*(eta*eta) - 5600.*xi*eta*eta*eta - 4480.*eta*eta/3. - 22400.*eta*eta*xi*xi/3. - 11200.*eta*eta*eta/3. + 44800.*(eta*eta*eta*eta)/9., 4480.*eta*xi/3. + 5320.*eta/9. - 12320.*eta*xi*xi + 89600.*eta*(xi*xi*xi)/9. - 1400.*xi/9. + 5600.*xi*(eta*eta) - 89600.*xi*eta*eta*eta/9. - 140./9. - 3920.*eta*eta + 8400.*(eta*eta)*(xi*xi) + 560.*(xi*xi) + 64960.*(eta*eta*eta)/9. + 2240.*(xi*xi*xi)/9. - 35000.*eta*eta*eta*eta/9. - 7000.*xi*xi*xi*xi/9.);
                      case 28:
                        return RealGradient(3360.*eta*xi + 280.*eta/9. - 11200.*eta*xi*xi/3. - 28000.*eta*xi*xi*xi/9. - 19040.*xi*eta*eta - 5600.*xi*eta*eta*eta - 11200.*eta*eta/3. + 112000.*(eta*eta)*(xi*xi)/3. + 52640.*(eta*eta*eta)/3. - 123200.*eta*eta*eta*eta/9., -35840.*eta*xi/3. - 8960.*eta/9. + 58240.*eta*(xi*xi) - 448000.*eta*xi*xi*xi/9. + 8680.*xi/9. - 30240.*xi*eta*eta + 246400.*xi*(eta*eta*eta)/9. - 140./9. + 9520.*(eta*eta) + 8400.*(eta*eta)*(xi*xi) - 1680.*xi*xi - 99680.*eta*eta*eta/9. - 24640.*xi*xi*xi/9. + 23800.*(eta*eta*eta*eta)/9. + 35000.*(xi*xi*xi*xi)/9.);
                      case 29:
                        return RealGradient(-26880.*eta*xi + 51800.*eta/9. + 79520.*eta*(xi*xi)/3. - 44800.*eta*xi*xi*xi/9. + 80640.*xi*(eta*eta) - 56000.*xi*eta*eta*eta - 69440.*eta*eta/3. - 123200.*eta*eta*xi*xi/3. + 82880.*(eta*eta*eta)/3. - 89600.*eta*eta*eta*eta/9., 89600.*eta*xi/3. - 9520.*eta/9. - 82880.*eta*xi*xi + 492800.*eta*(xi*xi*xi)/9. - 29680.*xi/9. - 51520.*xi*eta*eta + 179200.*xi*(eta*eta*eta)/9. + 1400./9. - 1120.*eta*eta + 84000.*(eta*eta)*(xi*xi) + 12320.*(xi*xi) + 71680.*(eta*eta*eta)/9. - 138880.*xi*xi*xi/9. - 56000.*eta*eta*eta*eta/9. + 56000.*(xi*xi*xi*xi)/9.);
                      case 30:
                        return RealGradient(11200.*eta*xi - 10360.*eta/9. - 75040.*eta*xi*xi/3. + 140000.*eta*(xi*xi*xi)/9. - 25760.*xi*eta*eta + 11200.*xi*(eta*eta*eta) + 7840.*(eta*eta)/3. + 95200.*(eta*eta)*(xi*xi)/3. + 2240.*(eta*eta*eta)/3. - 22400.*eta*eta*eta*eta/9., -31360.*eta*xi/3. - 12880.*eta/9. + 51520.*eta*(xi*xi) - 380800.*eta*xi*xi*xi/9. + 51800.*xi/9. - 2240.*xi*eta*eta + 44800.*xi*(eta*eta*eta)/9. - 700./9. + 5600.*(eta*eta) - 16800.*eta*eta*xi*xi - 28000.*xi*xi - 49280.*eta*eta*eta/9. + 375200.*(xi*xi*xi)/9. + 11200.*(eta*eta*eta*eta)/9. - 175000.*xi*xi*xi*xi/9.);
                      case 31:
                        return RealGradient(19040.*eta*xi - 8960.*eta/9. - 99680.*eta*xi*xi/3. + 95200.*eta*(xi*xi*xi)/9. - 30240.*xi*eta*eta + 5600.*xi*(eta*eta*eta) - 17920.*eta*eta/3. + 123200.*(eta*eta)*(xi*xi)/3. + 58240.*(eta*eta*eta)/3. - 112000.*eta*eta*eta*eta/9., -22400.*eta*xi/3. + 280.*eta/9. + 52640.*eta*(xi*xi) - 492800.*eta*xi*xi*xi/9. + 27160.*xi/9. - 19040.*xi*eta*eta + 224000.*xi*(eta*eta*eta)/9. - 980./9. + 1680.*(eta*eta) - 8400.*eta*eta*xi*xi - 16240.*xi*xi - 11200.*eta*eta*eta/9. + 239680.*(xi*xi*xi)/9. - 7000.*eta*eta*eta*eta/9. - 119000.*xi*xi*xi*xi/9.);
                      case 32:
                        return RealGradient(-7840.*eta*xi + 5320.*eta/9. + 64960.*eta*(xi*xi)/3. - 140000.*eta*xi*xi*xi/9. + 5600.*xi*(eta*eta) + 5600.*xi*(eta*eta*eta) + 2240.*(eta*eta)/3. - 44800.*eta*eta*xi*xi/3. - 12320.*eta*eta*eta/3. + 22400.*(eta*eta*eta*eta)/9., -8960.*eta*xi/3. + 2800.*eta/9. - 11200.*eta*xi*xi + 179200.*eta*(xi*xi*xi)/9. - 26600.*xi/9. + 12320.*xi*(eta*eta) - 44800.*xi*eta*eta*eta/9. + 700./9. - 560.*eta*eta - 8400.*eta*eta*xi*xi + 19600.*(xi*xi) - 1120.*eta*eta*eta/9. - 324800.*xi*xi*xi/9. + 1400.*(eta*eta*eta*eta)/9. + 175000.*(xi*xi*xi*xi)/9.);
                      case 33:
                        return RealGradient(3360.*eta*xi - 6440.*eta/9. - 4480.*eta*xi*xi + 5600.*eta*(xi*xi*xi)/3. - 19040.*xi*eta*eta + 16800.*xi*(eta*eta*eta) + 19040.*(eta*eta)/3. + 11200.*(eta*eta)*(xi*xi) - 11200.*eta*eta*eta + 5600.*(eta*eta*eta*eta), -29120.*eta*xi/3. + 6160.*eta/9. + 23520.*eta*(xi*xi) - 44800.*eta*xi*xi*xi/3. + 9520.*xi/9. + 20160.*xi*(eta*eta) - 11200.*xi*eta*eta*eta - 560./9. - 1680.*eta*eta - 25200.*eta*eta*xi*xi - 3920.*xi*xi + 1120.*(eta*eta*eta) + 15680.*(xi*xi*xi)/3. - 7000.*xi*xi*xi*xi/3.);
                      case 34:
                        return RealGradient(-3360.*eta*xi + 6160.*eta/9. + 3360.*eta*(xi*xi) + 20160.*xi*(eta*eta) - 16800.*xi*eta*eta*eta - 14560.*eta*eta/3. - 16800.*eta*eta*xi*xi + 7840.*(eta*eta*eta) - 11200.*eta*eta*eta*eta/3., 38080.*eta*xi/3. - 6440.*eta/9. - 33600.*eta*xi*xi + 22400.*eta*(xi*xi*xi) - 5600.*xi/9. - 19040.*xi*eta*eta + 22400.*xi*(eta*eta*eta)/3. + 280./9. + 1680.*(eta*eta) + 25200.*(eta*eta)*(xi*xi) + 1680.*(xi*xi) - 4480.*eta*eta*eta/3. - 1120.*xi*xi*xi + 1400.*(eta*eta*eta*eta)/3.);
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
                        return sign * RealGradient(-16800.*eta*xi + 3500.*eta + 22680.*eta*(xi*xi) - 8400.*eta*xi*xi*xi + 2800.*xi + 30240.*xi*(eta*eta) - 16800.*xi*eta*eta*eta - 375. - 10500.*eta*eta - 18900.*eta*eta*xi*xi - 6300.*xi*xi + 12600.*(eta*eta*eta) + 5040.*(xi*xi*xi) - 5250.*eta*eta*eta*eta - 1050.*xi*xi*xi*xi, 4200.*eta*xi - 15120.*eta*xi*xi + 12600.*eta*(xi*xi*xi) - 700.*xi - 7560.*xi*eta*eta + 4200.*xi*(eta*eta*eta) + 12600.*(eta*eta)*(xi*xi) + 4200.*(xi*xi) - 7560.*xi*xi*xi + 4200.*(xi*xi*xi*xi));
                      case 1:
                        return sign * RealGradient(7245.*eta*xi/2. + 11235.*eta/32. - 240975.*eta*xi*xi/32. + 48825.*eta*(xi*xi*xi)/16. - 12145.*xi/16. - 76545.*xi*eta*eta/16. + 14175.*xi*(eta*eta*eta)/8. + 695./128. - 127575.*eta*eta/64. + 382725.*(eta*eta)*(xi*xi)/64. + 136185.*(xi*xi)/64. + 104895.*(eta*eta*eta)/32. - 6825.*xi*xi*xi/4. - 212625.*eta*eta*eta*eta/128. + 19425.*(xi*xi*xi*xi)/64., 2835.*eta*xi/32. + 25515.*eta*(xi*xi)/8. - 127575.*eta*xi*xi*xi/32. + 3045.*xi/32. - 42525.*xi*eta*eta/32. + 42525.*xi*(eta*eta*eta)/32. - 42525.*eta*eta*xi*xi/32. - 39375.*xi*xi/32. + 84105.*(xi*xi*xi)/32. - 48825.*xi*xi*xi*xi/32.);
                      case 2:
                        return sign * RealGradient(1680.*eta*xi - 175.*eta/2. + 945.*eta*(xi*xi)/2. - 1575.*eta*xi*xi*xi - 4725.*xi*eta*eta + 3150.*xi*(eta*eta*eta) + 55./8. + 105.*(eta*eta)/4. + 4725.*(eta*eta)*(xi*xi)/4. - 3255.*xi*xi/4. + 735.*(eta*eta*eta)/2. + 1260.*(xi*xi*xi) - 2625.*eta*eta*eta*eta/8. - 1575.*xi*xi*xi*xi/4., -525.*eta*xi/2. + 1890.*eta*(xi*xi) - 1575.*eta*xi*xi*xi/2. + 35.*xi/2. + 315.*xi*(eta*eta)/2. + 525.*xi*(eta*eta*eta)/2. - 4725.*eta*eta*xi*xi/2. + 105.*(xi*xi)/2. - 1575.*xi*xi*xi/2. + 1575.*(xi*xi*xi*xi)/2.);
                      case 3:
                        return sign * RealGradient(-1155.*eta*xi/2. + 595.*eta/32. + 76545.*eta*(xi*xi)/32. - 9975.*eta*xi*xi*xi/16. + 1855.*xi/16. - 945.*xi*eta*eta/16. + 5775.*xi*(eta*eta*eta)/8. - 825./128. + 1785.*(eta*eta)/64. - 146475.*eta*eta*xi*xi/64. - 16695.*xi*xi/64. - 945.*eta*eta*eta/32. - 945.*xi*xi*xi/4. - 2625.*eta*eta*eta*eta/128. + 19425.*(xi*xi*xi*xi)/64., 1155.*eta*xi/32. - 4725.*eta*xi*xi/8. + 48825.*eta*(xi*xi*xi)/32. - 875.*xi/32. + 2835.*xi*(eta*eta)/32. + 525.*xi*(eta*eta*eta)/32. - 17325.*eta*eta*xi*xi/32. + 9345.*(xi*xi)/32. - 21735.*xi*xi*xi/32. + 9975.*(xi*xi*xi*xi)/32.);
                      case 4:
                        return sign * RealGradient(140.*xi - 5. - 840.*xi*xi + 1680.*(xi*xi*xi) - 1050.*xi*xi*xi*xi, 0.);
                      case 5:
                        return sign * RealGradient(140.*xi - 5. - 840.*xi*xi + 1680.*(xi*xi*xi) - 1050.*xi*xi*xi*xi, 0.);
                      case 6:
                        return sign * RealGradient(-1575.*eta*xi + 385.*eta/4. + 76545.*eta*(xi*xi)/16. - 14175.*eta*xi*xi*xi/4. + 735.*xi/4. + 36855.*xi*(eta*eta)/16. - 1575.*xi*eta*eta*eta/2. - 10. - 735.*eta*eta/4. - 127575.*eta*eta*xi*xi/32. - 2835.*xi*xi/4. + 1785.*(eta*eta*eta)/16. + 14175.*(xi*xi*xi)/16. - 2625.*eta*eta*eta*eta/128. - 42525.*xi*xi*xi*xi/128., 525.*eta*xi/2. - 31185.*eta*xi*xi/16. + 42525.*eta*(xi*xi*xi)/16. - 455.*xi/4. - 2205.*xi*eta*eta/16. + 525.*xi*(eta*eta*eta)/32. + 4725.*(eta*eta)*(xi*xi)/8. + 2205.*(xi*xi)/2. - 42525.*xi*xi*xi/16. + 14175.*(xi*xi*xi*xi)/8.);
                      case 7:
                        return sign * RealGradient(-2520.*eta*xi + 245.*eta + 4725.*eta*(xi*xi) - 2100.*eta*xi*xi*xi + 175.*xi + 6615.*xi*(eta*eta) - 4200.*xi*eta*eta*eta - 15. - 840.*eta*eta - 14175.*eta*eta*xi*xi/2. - 420.*xi*xi + 945.*(eta*eta*eta) + 315.*(xi*xi*xi) - 2625.*eta*eta*eta*eta/8. - 525.*xi*xi*xi*xi/8., 840.*eta*xi - 4725.*eta*xi*xi + 4725.*eta*(xi*xi*xi) - 175.*xi - 945.*xi*eta*eta + 525.*xi*(eta*eta*eta)/2. + 3150.*(eta*eta)*(xi*xi) + 1260.*(xi*xi) - 2205.*xi*xi*xi + 1050.*(xi*xi*xi*xi));
                      case 8:
                        return sign * RealGradient(-2205.*eta*xi + 1785.*eta/4. + 31185.*eta*(xi*xi)/16. - 1575.*eta*xi*xi*xi/4. + 455.*xi/4. + 127575.*xi*(eta*eta)/16. - 14175.*xi*eta*eta*eta/2. - 20. - 8505.*eta*eta/4. - 127575.*eta*eta*xi*xi/32. - 525.*xi*xi/4. + 53865.*(eta*eta*eta)/16. + 735.*(xi*xi*xi)/16. - 212625.*eta*eta*eta*eta/128. - 525.*xi*xi*xi*xi/128., 2835.*eta*xi/2. - 76545.*eta*xi*xi/16. + 42525.*eta*(xi*xi*xi)/16. - 735.*xi/4. - 42525.*xi*eta*eta/16. + 42525.*xi*(eta*eta*eta)/32. + 42525.*(eta*eta)*(xi*xi)/8. + 1575.*(xi*xi)/2. - 12285.*xi*xi*xi/16. + 1575.*(xi*xi*xi*xi)/8.);
                      case 9:
                        return sign * RealGradient(700.*eta - 25. - 4200.*eta*eta + 8400.*(eta*eta*eta) - 5250.*eta*eta*eta*eta, 1680.*eta*xi - 140.*xi - 5040.*xi*eta*eta + 4200.*xi*(eta*eta*eta));
                      case 10:
                        return sign * RealGradient(700.*eta - 25. - 4200.*eta*eta + 8400.*(eta*eta*eta) - 5250.*eta*eta*eta*eta, 1680.*eta*xi - 1260.*eta - 140.*xi - 5040.*xi*eta*eta + 4200.*xi*(eta*eta*eta) + 120. + 3360.*(eta*eta) - 2520.*eta*eta*eta);
                      case 11:
                        return sign * RealGradient(-9345.*eta*xi/16. - 595.*eta/4. + 4725.*eta*(xi*xi)/8. + 5775.*eta*(xi*xi*xi)/16. + 875.*xi/32. + 65205.*xi*(eta*eta)/32. - 9975.*xi*eta*eta*eta/8. + 555./128. + 66465.*(eta*eta)/64. - 146475.*eta*eta*xi*xi/64. - 1155.*xi*xi/64. - 36855.*eta*eta*eta/16. - 945.*xi*xi*xi/32. + 97125.*(eta*eta*eta*eta)/64. - 525.*xi*xi*xi*xi/128., 16695.*eta*xi/32. + 9345.*eta/32. - 76545.*eta*xi*xi/32. + 48825.*eta*(xi*xi*xi)/32. - 1855.*xi/16. + 2835.*xi*(eta*eta)/4. - 19425.*xi*eta*eta*eta/16. - 705./32. - 29085.*eta*eta/32. + 29925.*(eta*eta)*(xi*xi)/32. + 1155.*(xi*xi)/4. + 11655.*(eta*eta*eta)/16. + 315.*(xi*xi*xi)/16. - 5775.*xi*xi*xi*xi/32.);
                      case 12:
                        return sign * RealGradient(-105.*eta*xi + 315.*eta - 1890.*eta*xi*xi + 1575.*eta*(xi*xi*xi) - 35.*xi/2. + 4725.*xi*(eta*eta)/2. - 3150.*xi*eta*eta*eta - 85./8. - 7455.*eta*eta/4. + 4725.*(eta*eta)*(xi*xi)/4. + 525.*(xi*xi)/4. + 3465.*(eta*eta*eta) - 105.*xi*xi*xi/2. - 7875.*eta*eta*eta*eta/4. - 525.*xi*xi*xi*xi/8., 3255.*eta*xi/2. - 1155.*eta/2. - 945.*eta*xi*xi/2. - 1575.*eta*xi*xi*xi/2. - 3780.*xi*eta*eta + 1575.*xi*(eta*eta*eta) + 105./2. + 2835.*(eta*eta)/2. + 4725.*(eta*eta)*(xi*xi)/2. - 840.*xi*xi - 945.*eta*eta*eta + 1575.*(xi*xi*xi) - 1575.*xi*xi*xi*xi/2.);
                      case 13:
                        return sign * RealGradient(39375.*eta*xi/16. - 1715.*eta/4. - 25515.*eta*xi*xi/8. + 14175.*eta*(xi*xi*xi)/16. - 3045.*xi/32. - 252315.*xi*eta*eta/32. + 48825.*xi*(eta*eta*eta)/8. + 2395./128. + 125265.*(eta*eta)/64. + 382725.*(eta*eta)*(xi*xi)/64. - 2835.*xi*xi/64. - 48615.*eta*eta*eta/16. + 14175.*(xi*xi*xi)/32. + 97125.*(eta*eta*eta*eta)/64. - 42525.*xi*xi*xi*xi/128., -136185.*eta*xi/32. + 21105.*eta/32. + 240975.*eta*(xi*xi)/32. - 127575.*eta*xi*xi*xi/32. + 12145.*xi/16. + 20475.*xi*(eta*eta)/4. - 19425.*xi*eta*eta*eta/16. - 2865./32. - 40845.*eta*eta/32. - 146475.*eta*eta*xi*xi/32. - 7245.*xi*xi/4. + 11655.*(eta*eta*eta)/16. + 25515.*(xi*xi*xi)/16. - 14175.*xi*xi*xi*xi/32.);
                      case 14:
                        return sign * RealGradient(-8400.*eta*xi + 1400.*eta + 15120.*eta*(xi*xi) - 8400.*eta*xi*xi*xi + 700.*xi + 22680.*xi*(eta*eta) - 16800.*xi*eta*eta*eta - 75. - 6300.*eta*eta - 18900.*eta*eta*xi*xi - 2100.*xi*xi + 10080.*(eta*eta*eta) + 2520.*(xi*xi*xi) - 5250.*eta*eta*eta*eta - 1050.*xi*xi*xi*xi, 12600.*eta*xi - 2100.*eta - 22680.*eta*xi*xi + 12600.*eta*(xi*xi*xi) - 2800.*xi - 15120.*xi*eta*eta + 4200.*xi*(eta*eta*eta) + 300. + 4200.*(eta*eta) + 12600.*(eta*eta)*(xi*xi) + 8400.*(xi*xi) - 2520.*eta*eta*eta - 10080.*xi*xi*xi + 4200.*(xi*xi*xi*xi));
                      case 15:
                        return RealGradient(188160.*eta*xi - 49000.*eta - 211680.*eta*xi*xi + 67200.*eta*(xi*xi*xi) - 19600.*xi - 423360.*xi*eta*eta + 268800.*xi*(eta*eta*eta) + 3500. + 176400.*(eta*eta) + 226800.*(eta*eta)*(xi*xi) + 35280.*(xi*xi) - 235200.*eta*eta*eta - 23520.*xi*xi*xi + 105000.*(eta*eta*eta*eta) + 4200.*(xi*xi*xi*xi), -70560.*eta*xi + 211680.*eta*(xi*xi) - 151200.*eta*xi*xi*xi + 9800.*xi + 141120.*xi*(eta*eta) - 84000.*xi*eta*eta*eta - 201600.*eta*eta*xi*xi - 47040.*xi*xi + 70560.*(xi*xi*xi) - 33600.*xi*xi*xi*xi);
                      case 16:
                        return RealGradient(-94080.*eta*xi + 9800.*eta + 211680.*eta*(xi*xi) - 134400.*eta*xi*xi*xi + 9800.*xi + 211680.*xi*(eta*eta) - 134400.*xi*eta*eta*eta - 700. - 35280.*eta*eta - 226800.*eta*eta*xi*xi - 35280.*xi*xi + 47040.*(eta*eta*eta) + 47040.*(xi*xi*xi) - 21000.*eta*eta*eta*eta - 21000.*xi*xi*xi*xi, 70560.*eta*xi - 211680.*eta*xi*xi + 151200.*eta*(xi*xi*xi) - 19600.*xi - 70560.*xi*eta*eta + 16800.*xi*(eta*eta*eta) + 100800.*(eta*eta)*(xi*xi) + 94080.*(xi*xi) - 141120.*xi*xi*xi + 67200.*(xi*xi*xi*xi));
                      case 17:
                        return RealGradient(-13440.*eta*xi + 560.*eta + 60480.*eta*(xi*xi) - 67200.*eta*xi*xi*xi + 6440.*xi - 280. - 30240.*xi*xi + 43680.*(xi*xi*xi) - 16800.*xi*xi*xi*xi, -1120.*xi + 13440.*(xi*xi) - 40320.*xi*xi*xi + 33600.*(xi*xi*xi*xi));
                      case 18:
                        return RealGradient(6720.*eta*xi - 280.*eta - 30240.*eta*xi*xi + 33600.*eta*(xi*xi*xi) - 3640.*xi + 140. + 20160.*(xi*xi) - 36960.*xi*xi*xi + 21000.*(xi*xi*xi*xi), 560.*xi - 6720.*xi*xi + 20160.*(xi*xi*xi) - 16800.*xi*xi*xi*xi);
                      case 19:
                        return RealGradient(-13440.*eta*xi + 18200.*eta + 560.*xi + 60480.*xi*(eta*eta) - 67200.*xi*eta*eta*eta - 700. - 100800.*eta*eta + 184800.*(eta*eta*eta) - 105000.*eta*eta*eta*eta, 40320.*eta*xi - 30240.*eta*xi*xi - 3640.*xi - 110880.*xi*eta*eta + 84000.*xi*(eta*eta*eta) + 50400.*(eta*eta)*(xi*xi) + 3360.*(xi*xi));
                      case 20:
                        return RealGradient(26880.*eta*xi - 14560.*eta - 1120.*xi - 120960.*xi*eta*eta + 134400.*xi*(eta*eta*eta) + 560. + 80640.*(eta*eta) - 147840.*eta*eta*eta + 84000.*(eta*eta*eta*eta), -60480.*eta*xi + 60480.*eta*(xi*xi) + 6440.*xi + 131040.*xi*(eta*eta) - 67200.*xi*eta*eta*eta - 100800.*eta*eta*xi*xi - 6720.*xi*xi);
                      case 21:
                        return RealGradient(-58240.*eta*xi + 14840.*eta/3. + 125440.*eta*(xi*xi) - 582400.*eta*xi*xi*xi/9. + 17920.*xi/3. + 116480.*xi*(eta*eta) - 560000.*xi*eta*eta*eta/9. - 420. - 13440.*eta*eta - 134400.*eta*eta*xi*xi - 17920.*xi*xi + 115360.*(eta*eta*eta)/9. + 156800.*(xi*xi*xi)/9. - 35000.*eta*eta*eta*eta/9. - 44800.*xi*xi*xi*xi/9., 13440.*eta*xi - 81760.*eta*xi*xi + 89600.*eta*(xi*xi*xi) - 10360.*xi/3. - 39200.*xi*eta*eta/3. + 28000.*xi*(eta*eta*eta)/9. + 140000.*(eta*eta)*(xi*xi)/3. + 28000.*(xi*xi) - 170240.*xi*xi*xi/3. + 291200.*(xi*xi*xi*xi)/9.);
                      case 22:
                        return RealGradient(44800.*eta*xi - 10640.*eta/3. - 109760.*eta*xi*xi + 627200.*eta*(xi*xi*xi)/9. - 12880.*xi/3. - 91840.*xi*eta*eta + 448000.*xi*(eta*eta*eta)/9. + 280. + 10080.*(eta*eta) + 117600.*(eta*eta)*(xi*xi) + 14560.*(xi*xi) - 89600.*eta*eta*eta/9. - 152320.*xi*xi*xi/9. + 28000.*(eta*eta*eta*eta)/9. + 56000.*(xi*xi*xi*xi)/9., -13440.*eta*xi + 76160.*eta*(xi*xi) - 78400.*eta*xi*xi*xi + 12040.*xi/3. + 34720.*xi*(eta*eta)/3. - 22400.*xi*eta*eta*eta/9. - 112000.*eta*eta*xi*xi/3. - 31360.*xi*xi + 185920.*(xi*xi*xi)/3. - 313600.*xi*xi*xi*xi/9.);
                      case 23:
                        return RealGradient(-62720.*eta*xi + 32480.*eta/3. + 76160.*eta*(xi*xi) - 224000.*eta*xi*xi*xi/9. + 12040.*xi/3. + 185920.*xi*(eta*eta) - 1254400.*xi*eta*eta*eta/9. - 560. - 45920.*eta*eta - 117600.*eta*eta*xi*xi - 6720.*xi*xi + 600320.*(eta*eta*eta)/9. + 34720.*(xi*xi*xi)/9. - 280000.*eta*eta*eta*eta/9. - 5600.*xi*xi*xi*xi/9., 29120.*eta*xi - 109760.*eta*xi*xi + 78400.*eta*(xi*xi*xi) - 12880.*xi/3. - 152320.*xi*eta*eta/3. + 224000.*xi*(eta*eta*eta)/9. + 313600.*(eta*eta)*(xi*xi)/3. + 22400.*(xi*xi) - 91840.*xi*xi*xi/3. + 112000.*(xi*xi*xi*xi)/9.);
                      case 24:
                        return RealGradient(56000.*eta*xi - 24920.*eta/3. - 81760.*eta*xi*xi + 280000.*eta*(xi*xi*xi)/9. - 10360.*xi/3. - 170240.*xi*eta*eta + 1164800.*xi*(eta*eta*eta)/9. + 420. + 35840.*(eta*eta) + 134400.*(eta*eta)*(xi*xi) + 6720.*(xi*xi) - 474880.*eta*eta*eta/9. - 39200.*xi*xi*xi/9. + 224000.*(eta*eta*eta*eta)/9. + 7000.*(xi*xi*xi*xi)/9., -35840.*eta*xi + 125440.*eta*(xi*xi) - 89600.*eta*xi*xi*xi + 17920.*xi/3. + 156800.*xi*(eta*eta)/3. - 179200.*xi*eta*eta*eta/9. - 291200.*eta*eta*xi*xi/3. - 29120.*xi*xi + 116480.*(xi*xi*xi)/3. - 140000.*xi*xi*xi*xi/9.);
                      case 25:
                        return RealGradient(-31360.*eta*xi/3. + 51800.*eta/9. - 2240.*eta*xi*xi + 44800.*eta*(xi*xi*xi)/9. - 12880.*xi/9. + 51520.*xi*(eta*eta) - 380800.*xi*eta*eta*eta/9. - 700./9. - 28000.*eta*eta - 16800.*eta*eta*xi*xi + 5600.*(xi*xi) + 375200.*(eta*eta*eta)/9. - 49280.*xi*xi*xi/9. - 175000.*eta*eta*eta*eta/9. + 11200.*(xi*xi*xi*xi)/9., 11200.*eta*xi - 25760.*eta*xi*xi + 11200.*eta*(xi*xi*xi) - 10360.*xi/9. - 75040.*xi*eta*eta/3. + 140000.*xi*(eta*eta*eta)/9. + 95200.*(eta*eta)*(xi*xi)/3. + 7840.*(xi*xi)/3. + 2240.*(xi*xi*xi)/3. - 22400.*xi*xi*xi*xi/9.);
                      case 26:
                        return RealGradient(89600.*eta*xi/3. - 29680.*eta/9. - 51520.*eta*xi*xi + 179200.*eta*(xi*xi*xi)/9. - 9520.*xi/9. - 82880.*xi*eta*eta + 492800.*xi*(eta*eta*eta)/9. + 1400./9. + 12320.*(eta*eta) + 84000.*(eta*eta)*(xi*xi) - 1120.*xi*xi - 138880.*eta*eta*eta/9. + 71680.*(xi*xi*xi)/9. + 56000.*(eta*eta*eta*eta)/9. - 56000.*xi*xi*xi*xi/9., -26880.*eta*xi + 80640.*eta*(xi*xi) - 56000.*eta*xi*xi*xi + 51800.*xi/9. + 79520.*xi*(eta*eta)/3. - 44800.*xi*eta*eta*eta/9. - 123200.*eta*eta*xi*xi/3. - 69440.*xi*xi/3. + 82880.*(xi*xi*xi)/3. - 89600.*xi*xi*xi*xi/9.);
                      case 27:
                        return RealGradient(-8960.*eta*xi/3. - 26600.*eta/9. + 12320.*eta*(xi*xi) - 44800.*eta*xi*xi*xi/9. + 2800.*xi/9. - 11200.*xi*eta*eta + 179200.*xi*(eta*eta*eta)/9. + 700./9. + 19600.*(eta*eta) - 8400.*eta*eta*xi*xi - 560.*xi*xi - 324800.*eta*eta*eta/9. - 1120.*xi*xi*xi/9. + 175000.*(eta*eta*eta*eta)/9. + 1400.*(xi*xi*xi*xi)/9., -7840.*eta*xi + 5600.*eta*(xi*xi) + 5600.*eta*(xi*xi*xi) + 5320.*xi/9. + 64960.*xi*(eta*eta)/3. - 140000.*xi*eta*eta*eta/9. - 44800.*eta*eta*xi*xi/3. + 2240.*(xi*xi)/3. - 12320.*xi*xi*xi/3. + 22400.*(xi*xi*xi*xi)/9.);
                      case 28:
                        return RealGradient(-22400.*eta*xi/3. + 27160.*eta/9. - 19040.*eta*xi*xi + 224000.*eta*(xi*xi*xi)/9. + 280.*xi/9. + 52640.*xi*(eta*eta) - 492800.*xi*eta*eta*eta/9. - 980./9. - 16240.*eta*eta - 8400.*eta*eta*xi*xi + 1680.*(xi*xi) + 239680.*(eta*eta*eta)/9. - 11200.*xi*xi*xi/9. - 119000.*eta*eta*eta*eta/9. - 7000.*xi*xi*xi*xi/9., 19040.*eta*xi - 30240.*eta*xi*xi + 5600.*eta*(xi*xi*xi) - 8960.*xi/9. - 99680.*xi*eta*eta/3. + 95200.*xi*(eta*eta*eta)/9. + 123200.*(eta*eta)*(xi*xi)/3. - 17920.*xi*xi/3. + 58240.*(xi*xi*xi)/3. - 112000.*xi*xi*xi*xi/9.);
                      case 29:
                        return RealGradient(-138880.*eta*xi/3. - 48160.*eta/9. + 80640.*eta*(xi*xi) - 246400.*eta*xi*xi*xi/9. + 51800.*xi/9. + 82880.*xi*(eta*eta) - 358400.*xi*eta*eta*eta/9. + 560./9. + 32480.*(eta*eta) - 84000.*eta*eta*xi*xi - 13440.*xi*xi - 519680.*eta*eta*eta/9. + 79520.*(xi*xi*xi)/9. + 280000.*(eta*eta*eta*eta)/9. - 11200.*xi*xi*xi*xi/9., -2240.*eta*xi - 51520.*eta*xi*xi + 56000.*eta*(xi*xi*xi) - 9520.*xi/9. + 71680.*xi*(eta*eta)/3. - 224000.*xi*eta*eta*eta/9. + 89600.*(eta*eta)*(xi*xi)/3. + 44800.*(xi*xi)/3. - 82880.*xi*xi*xi/3. + 123200.*(xi*xi*xi*xi)/9.);
                      case 30:
                        return RealGradient(15680.*eta*xi/3. + 6440.*eta/9. - 25760.*eta*xi*xi + 190400.*eta*(xi*xi*xi)/9. - 10360.*xi/9. + 2240.*xi*(eta*eta) - 89600.*xi*eta*eta*eta/9. + 140./9. - 5600.*eta*eta + 16800.*(eta*eta)*(xi*xi) + 5600.*(xi*xi) + 98560.*(eta*eta*eta)/9. - 75040.*xi*xi*xi/9. - 56000.*eta*eta*eta*eta/9. + 35000.*(xi*xi*xi*xi)/9., 11200.*eta*xi - 2240.*eta*xi*xi - 11200.*eta*xi*xi*xi - 12880.*xi/9. - 49280.*xi*eta*eta/3. + 44800.*xi*(eta*eta*eta)/9. + 22400.*(eta*eta)*(xi*xi)/3. - 15680.*xi*xi/3. + 51520.*(xi*xi*xi)/3. - 95200.*xi*xi*xi*xi/9.);
                      case 31:
                        return RealGradient(-35840.*eta*xi/3. + 8680.*eta/9. - 30240.*eta*xi*xi + 246400.*eta*(xi*xi*xi)/9. - 8960.*xi/9. + 58240.*xi*(eta*eta) - 448000.*xi*eta*eta*eta/9. - 140./9. - 1680.*eta*eta + 8400.*(eta*eta)*(xi*xi) + 9520.*(xi*xi) - 24640.*eta*eta*eta/9. - 99680.*xi*xi*xi/9. + 35000.*(eta*eta*eta*eta)/9. + 23800.*(xi*xi*xi*xi)/9., 3360.*eta*xi - 19040.*eta*xi*xi - 5600.*eta*xi*xi*xi + 280.*xi/9. - 11200.*xi*eta*eta/3. - 28000.*xi*eta*eta*eta/9. + 112000.*(eta*eta)*(xi*xi)/3. - 11200.*xi*xi/3. + 52640.*(xi*xi*xi)/3. - 123200.*xi*xi*xi*xi/9.);
                      case 32:
                        return RealGradient(4480.*eta*xi/3. - 1400.*eta/9. + 5600.*eta*(xi*xi) - 89600.*eta*xi*xi*xi/9. + 5320.*xi/9. - 12320.*xi*eta*eta + 89600.*xi*(eta*eta*eta)/9. - 140./9. + 560.*(eta*eta) + 8400.*(eta*eta)*(xi*xi) - 3920.*xi*xi + 2240.*(eta*eta*eta)/9. + 64960.*(xi*xi*xi)/9. - 7000.*eta*eta*eta*eta/9. - 35000.*xi*xi*xi*xi/9., -1120.*eta*xi + 12320.*eta*(xi*xi) - 5600.*eta*xi*xi*xi + 2800.*xi/9. - 1120.*xi*eta*eta/3. + 5600.*xi*(eta*eta*eta)/9. - 22400.*eta*eta*xi*xi/3. - 4480.*xi*xi/3. - 11200.*xi*xi*xi/3. + 44800.*(xi*xi*xi*xi)/9.);
                      case 33:
                        return RealGradient(38080.*eta*xi/3. - 5600.*eta/9. - 19040.*eta*xi*xi + 22400.*eta*(xi*xi*xi)/3. - 6440.*xi/9. - 33600.*xi*eta*eta + 22400.*xi*(eta*eta*eta) + 280./9. + 1680.*(eta*eta) + 25200.*(eta*eta)*(xi*xi) + 1680.*(xi*xi) - 1120.*eta*eta*eta - 4480.*xi*xi*xi/3. + 1400.*(xi*xi*xi*xi)/3., -3360.*eta*xi + 20160.*eta*(xi*xi) - 16800.*eta*xi*xi*xi + 6160.*xi/9. + 3360.*xi*(eta*eta) - 16800.*eta*eta*xi*xi - 14560.*xi*xi/3. + 7840.*(xi*xi*xi) - 11200.*xi*xi*xi*xi/3.);
                      case 34:
                        return RealGradient(-29120.*eta*xi/3. + 9520.*eta/9. + 20160.*eta*(xi*xi) - 11200.*eta*xi*xi*xi + 6160.*xi/9. + 23520.*xi*(eta*eta) - 44800.*xi*eta*eta*eta/3. - 560./9. - 3920.*eta*eta - 25200.*eta*eta*xi*xi - 1680.*xi*xi + 15680.*(eta*eta*eta)/3. + 1120.*(xi*xi*xi) - 7000.*eta*eta*eta*eta/3., 3360.*eta*xi - 19040.*eta*xi*xi + 16800.*eta*(xi*xi*xi) - 6440.*xi/9. - 4480.*xi*eta*eta + 5600.*xi*(eta*eta*eta)/3. + 11200.*(eta*eta)*(xi*xi) + 19040.*(xi*xi)/3. - 11200.*xi*xi*xi + 5600.*(xi*xi*xi*xi));
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
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << totalorder);
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

  const Order totalorder = order + add_p_level*elem->p_level();
  libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

  const char sign = i >= totalorder * elem->n_edges() || elem->edge_orientation(i / totalorder) ? 1 : -1;
  const unsigned int ii = sign > 0 ? i : (i / totalorder * 2 + 1) * totalorder - 1 - i;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (totalorder)
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
                        return sign * RealGradient(-135.*eta/4. - 105./4. + 135.*((eta + 1.)*(eta + 1.))/4. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 1:
                        return sign * RealGradient(135.*eta/8. + 105./8. - 135.*(eta + 1.)*(eta + 1.)/8. + 75.*((eta + 1.)*(eta + 1.)*(eta + 1.))/16., 0.);
                      case 2:
                        return sign * RealGradient(-135.*eta/4. - 105./4. + 135.*((eta + 1.)*(eta + 1.))/4. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 3:
                        return sign * RealGradient(0., 54.*eta + 135.*xi/4. - 135.*(eta + 1.)*(xi + 1.)/2. + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 243./4. - 45.*(eta + 1.)*(eta + 1.)/2.);
                      case 4:
                        return sign * RealGradient(0., -45.*eta/2. - 45.*xi/8. + 225.*(eta + 1.)*(xi + 1.)/8. - 225.*(xi + 1.)*(eta + 1.)*(eta + 1.)/16. - 189./8. + 45.*((eta + 1.)*(eta + 1.))/4.);
                      case 5:
                        return sign * RealGradient(0., 36.*eta + 45.*xi/4. - 45.*(eta + 1.)*(xi + 1.) + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 153./4. - 45.*(eta + 1.)*(eta + 1.)/2.);
                      case 6:
                        return sign * RealGradient(-45.*eta/4. - 45./4. + 45.*((eta + 1.)*(eta + 1.))/2. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 7:
                        return sign * RealGradient(45.*eta/8. + 45./8. - 45.*(eta + 1.)*(eta + 1.)/4. + 75.*((eta + 1.)*(eta + 1.)*(eta + 1.))/16., 0.);
                      case 8:
                        return sign * RealGradient(-45.*eta/4. - 45./4. + 45.*((eta + 1.)*(eta + 1.))/2. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 9:
                        return sign * RealGradient(0., 54.*eta + 45.*xi/4. - 45.*(eta + 1.)*(xi + 1.) + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 207./4. - 135.*(eta + 1.)*(eta + 1.)/4.);
                      case 10:
                        return sign * RealGradient(0., -135.*eta/4. - 45.*xi/8. + 225.*(eta + 1.)*(xi + 1.)/8. - 225.*(xi + 1.)*(eta + 1.)*(eta + 1.)/16. - 261./8. + 135.*((eta + 1.)*(eta + 1.))/8.);
                      case 11:
                        return sign * RealGradient(0., 81.*eta + 135.*xi/4. - 135.*(eta + 1.)*(xi + 1.)/2. + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 297./4. - 135.*(eta + 1.)*(eta + 1.)/4.);
                      case 12:
                        return RealGradient(0., 81.*eta + 135.*xi/4. - 135.*(eta + 1.)*(xi + 1.)/2. + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 297./4. - 135.*(eta + 1.)*(eta + 1.)/4.);
                      case 13:
                        return RealGradient(0., -54.*eta - 135.*xi/4. + 135.*(eta + 1.)*(xi + 1.)/2. - 225.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 243./4. + 45.*((eta + 1.)*(eta + 1.))/2.);
                      case 14:
                        return RealGradient(-30.*eta - 30. + 135.*((eta + 1.)*(eta + 1.))/4. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 15:
                        return RealGradient(15.*eta/2. + 15./2. - 45.*(eta + 1.)*(eta + 1.)/2. + 75.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8., 0.);
                      case 16:
                        return RealGradient(-30.*eta - 30. + 135.*((eta + 1.)*(eta + 1.))/4. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 17:
                        return RealGradient(15.*eta/2. + 15./2. - 45.*(eta + 1.)*(eta + 1.)/2. + 75.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8., 0.);
                      case 18:
                        return RealGradient(0., 54.*eta + 45.*xi/4. - 45.*(eta + 1.)*(xi + 1.) + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 207./4. - 135.*(eta + 1.)*(eta + 1.)/4.);
                      case 19:
                        return RealGradient(0., -36.*eta - 45.*xi/4. + 45.*(eta + 1.)*(xi + 1.) - 225.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 153./4. + 45.*((eta + 1.)*(eta + 1.))/2.);
                      case 20:
                        return RealGradient(0., 0.);
                      case 21:
                        return RealGradient(0., 15.*xi/4. - 3./4.);
                      case 22:
                        return RealGradient(0., -15.*xi/4. - 3./4.);
                      case 23:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(-81.*eta - 135.*xi/4. + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 297./4. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 1:
                        return sign * RealGradient(135.*eta/4. + 135.*xi/8. - 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 135./4. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)/4., 0.);
                      case 2:
                        return sign * RealGradient(-54.*eta - 135.*xi/4. + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 243./4. + 90.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 3:
                        return sign * RealGradient(0., 45.*eta/4. + 54.*xi - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 207./4. - 135.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 4:
                        return sign * RealGradient(0., -45.*eta/8. - 45.*xi/2. + 90.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. - 45./2. + 225.*((xi + 1.)*(xi + 1.)/(2.*2.))/4.);
                      case 5:
                        return sign * RealGradient(0., 45.*eta/4. + 36.*xi - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 153./4. - 90.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 6:
                        return sign * RealGradient(-36.*eta - 45.*xi/4. + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 153./4. + 90.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 7:
                        return sign * RealGradient(45.*eta/2. + 45.*xi/8. - 90.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 45./2. - 225.*(eta + 1.)*(eta + 1.)/(2.*2.)/4., 0.);
                      case 8:
                        return sign * RealGradient(-54.*eta - 45.*xi/4. + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 207./4. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 9:
                        return sign * RealGradient(0., 135.*eta/4. + 54.*xi - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 243./4. - 90.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 10:
                        return sign * RealGradient(0., -135.*eta/8. - 135.*xi/4. + 135.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. - 135./4. + 225.*((xi + 1.)*(xi + 1.)/(2.*2.))/4.);
                      case 11:
                        return sign * RealGradient(0., 135.*eta/4. + 81.*xi - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 297./4. - 135.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 12:
                        return RealGradient(0., 30.*eta + 81.*xi - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 75. - 135.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 13:
                        return RealGradient(0., -15.*eta/2. - 54.*xi + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 105./2. + 135.*((xi + 1.)*(xi + 1.)/(2.*2.)));
                      case 14:
                        return RealGradient(-81.*eta - 30.*xi + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 75. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 15:
                        return RealGradient(54.*eta + 15.*xi/2. - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 105./2. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.), 0.);
                      case 16:
                        return RealGradient(-54.*eta - 30.*xi + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 60. + 90.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 17:
                        return RealGradient(36.*eta + 15.*xi/2. - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 75./2. - 90.*(eta + 1.)*(eta + 1.)/(2.*2.), 0.);
                      case 18:
                        return RealGradient(0., 30.*eta + 54.*xi - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 60. - 90.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 19:
                        return RealGradient(0., -15.*eta/2. - 36.*xi + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 75./2. + 90.*((xi + 1.)*(xi + 1.)/(2.*2.)));
                      case 20:
                        return RealGradient(0., 0.);
                      case 21:
                        return RealGradient(0., 0.);
                      case 22:
                        return RealGradient(0., 0.);
                      case 23:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(-135.*eta/4. - 81.*xi + 135.*(eta + 1.)*(xi + 1.)/2. - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 297./4. + 135.*((xi + 1.)*(xi + 1.))/4., 0.);
                      case 1:
                        return sign * RealGradient(45.*eta/8. + 135.*xi/4. - 225.*(eta + 1.)*(xi + 1.)/8. + 225.*(eta + 1.)*((xi + 1.)*(xi + 1.))/16. + 261./8. - 135.*(xi + 1.)*(xi + 1.)/8., 0.);
                      case 2:
                        return sign * RealGradient(-45.*eta/4. - 54.*xi + 45.*(eta + 1.)*(xi + 1.) - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 207./4. + 135.*((xi + 1.)*(xi + 1.))/4., 0.);
                      case 3:
                        return sign * RealGradient(0., 45.*xi/4. + 45./4. - 45.*(xi + 1.)*(xi + 1.)/2. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 4:
                        return sign * RealGradient(0., -45.*xi/8. - 45./8. + 45.*((xi + 1.)*(xi + 1.))/4. - 75.*(xi + 1.)*(xi + 1.)*(xi + 1.)/16.);
                      case 5:
                        return sign * RealGradient(0., 45.*xi/4. + 45./4. - 45.*(xi + 1.)*(xi + 1.)/2. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 6:
                        return sign * RealGradient(-45.*eta/4. - 36.*xi + 45.*(eta + 1.)*(xi + 1.) - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 153./4. + 45.*((xi + 1.)*(xi + 1.))/2., 0.);
                      case 7:
                        return sign * RealGradient(45.*eta/8. + 45.*xi/2. - 225.*(eta + 1.)*(xi + 1.)/8. + 225.*(eta + 1.)*((xi + 1.)*(xi + 1.))/16. + 189./8. - 45.*(xi + 1.)*(xi + 1.)/4., 0.);
                      case 8:
                        return sign * RealGradient(-135.*eta/4. - 54.*xi + 135.*(eta + 1.)*(xi + 1.)/2. - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 243./4. + 45.*((xi + 1.)*(xi + 1.))/2., 0.);
                      case 9:
                        return sign * RealGradient(0., 135.*xi/4. + 105./4. - 135.*(xi + 1.)*(xi + 1.)/4. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 10:
                        return sign * RealGradient(0., -135.*xi/8. - 105./8. + 135.*((xi + 1.)*(xi + 1.))/8. - 75.*(xi + 1.)*(xi + 1.)*(xi + 1.)/16.);
                      case 11:
                        return sign * RealGradient(0., 135.*xi/4. + 105./4. - 135.*(xi + 1.)*(xi + 1.)/4. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 12:
                        return RealGradient(0., 30.*xi + 30. - 135.*(xi + 1.)*(xi + 1.)/4. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 13:
                        return RealGradient(0., -15.*xi/2. - 15./2. + 45.*((xi + 1.)*(xi + 1.))/2. - 75.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8.);
                      case 14:
                        return RealGradient(-135.*eta/4. - 81.*xi + 135.*(eta + 1.)*(xi + 1.)/2. - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 297./4. + 135.*((xi + 1.)*(xi + 1.))/4., 0.);
                      case 15:
                        return RealGradient(135.*eta/4. + 54.*xi - 135.*(eta + 1.)*(xi + 1.)/2. + 225.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. + 243./4. - 45.*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 16:
                        return RealGradient(-45.*eta/4. - 54.*xi + 45.*(eta + 1.)*(xi + 1.) - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 207./4. + 135.*((xi + 1.)*(xi + 1.))/4., 0.);
                      case 17:
                        return RealGradient(45.*eta/4. + 36.*xi - 45.*(eta + 1.)*(xi + 1.) + 225.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. + 153./4. - 45.*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 18:
                        return RealGradient(0., 30.*xi + 30. - 135.*(xi + 1.)*(xi + 1.)/4. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 19:
                        return RealGradient(0., -15.*xi/2. - 15./2. + 45.*((xi + 1.)*(xi + 1.))/2. - 75.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8.);
                      case 20:
                        return RealGradient(15.*eta/4. - 3./4., 0.);
                      case 21:
                        return RealGradient(0., 0.);
                      case 22:
                        return RealGradient(0., 0.);
                      case 23:
                        return RealGradient(-15.*eta/4. - 3./4., 0.);
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
                        return sign * RealGradient(60. - 90.*eta, 180.*eta + 270.*xi - 120.);
                      case 1:
                        return sign * RealGradient(45.*eta - 30., -45.*eta - 135.*xi + 105./2.);
                      case 2:
                        return sign * RealGradient(60. - 90.*eta, 270.*xi - 90.);
                      case 3:
                        return sign * RealGradient(-90.*eta, 270.*xi - 90.);
                      case 4:
                        return sign * RealGradient(-45.*eta/2., 90.*eta + 135.*xi/2. - 75./2.);
                      case 5:
                        return sign * RealGradient(0., 0.);
                      case 6:
                        return sign * RealGradient(0., 0.);
                      case 7:
                        return sign * RealGradient(-45.*eta/2., -45.*eta + 135.*xi/2. - 30.);
                      case 8:
                        return sign * RealGradient(-90.*eta, 180.*eta + 270.*xi - 180.);
                      case 9:
                        return RealGradient(180.*eta, -720.*eta - 540.*xi + 300.);
                      case 10:
                        return RealGradient(-540.*eta, 720.*eta + 1620.*xi - 900.);
                      case 11:
                        return RealGradient(-360.*eta, 720.*eta + 1080.*xi - 480.);
                      case 12:
                        return RealGradient(540.*eta, -360.*eta - 1620.*xi + 720.);
                      case 13:
                        return RealGradient(0., 360.*eta - 60.);
                      case 14:
                        return RealGradient(0., 120. - 720.*eta);
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
                        return sign * RealGradient(-180.*eta - 90.*xi + 120., 90.*eta + 180.*xi - 60.);
                      case 1:
                        return sign * RealGradient(45.*eta + 45.*xi - 75./2., 45.*eta/2. - 45.*xi);
                      case 2:
                        return sign * RealGradient(30. - 90.*xi, 0.);
                      case 3:
                        return sign * RealGradient(30. - 90.*xi, 0.);
                      case 4:
                        return sign * RealGradient(-90.*eta - 45.*xi/2. + 45./2., 45.*eta/2. + 90.*xi - 45./2.);
                      case 5:
                        return sign * RealGradient(0., 90.*eta - 30.);
                      case 6:
                        return sign * RealGradient(0., 90.*eta - 30.);
                      case 7:
                        return sign * RealGradient(45.*eta - 45.*xi/2., -45.*eta - 45.*xi + 75./2.);
                      case 8:
                        return sign * RealGradient(-180.*eta - 90.*xi + 60., 90.*eta + 180.*xi - 120.);
                      case 9:
                        return RealGradient(720.*eta + 180.*xi - 300., -540.*eta - 720.*xi + 300.);
                      case 10:
                        return RealGradient(-720.*eta - 540.*xi + 300., 180.*eta + 720.*xi - 300.);
                      case 11:
                        return RealGradient(-720.*eta - 360.*xi + 360., 720.*xi - 120.);
                      case 12:
                        return RealGradient(360.*eta + 540.*xi - 240., 60. - 360.*xi);
                      case 13:
                        return RealGradient(60. - 360.*eta, 540.*eta + 360.*xi - 240.);
                      case 14:
                        return RealGradient(720.*eta - 120., -360.*eta - 720.*xi + 360.);
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
                        return sign * RealGradient(-270.*eta - 180.*xi + 180., 90.*xi);
                      case 1:
                        return sign * RealGradient(-135.*eta/2. + 45.*xi + 30., 45.*xi/2.);
                      case 2:
                        return sign * RealGradient(0., 0.);
                      case 3:
                        return sign * RealGradient(0., 0.);
                      case 4:
                        return sign * RealGradient(-135.*eta/2. - 90.*xi + 75./2., 45.*xi/2.);
                      case 5:
                        return sign * RealGradient(90. - 270.*eta, 90.*xi);
                      case 6:
                        return sign * RealGradient(90. - 270.*eta, 90.*xi - 60.);
                      case 7:
                        return sign * RealGradient(135.*eta + 45.*xi - 105./2., 30. - 45.*xi);
                      case 8:
                        return sign * RealGradient(-270.*eta - 180.*xi + 120., 90.*xi - 60.);
                      case 9:
                        return RealGradient(1620.*eta + 720.*xi - 900., -540.*xi);
                      case 10:
                        return RealGradient(-540.*eta - 720.*xi + 300., 180.*xi);
                      case 11:
                        return RealGradient(120. - 720.*xi, 0.);
                      case 12:
                        return RealGradient(360.*xi - 60., 0.);
                      case 13:
                        return RealGradient(-1620.*eta - 360.*xi + 720., 540.*xi);
                      case 14:
                        return RealGradient(1080.*eta + 720.*xi - 480., -360.*xi);
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
                        return sign * RealGradient(-480.*eta - 105.*xi/2. + 420.*(eta + 1.)*(xi + 1.) - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 525.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 3675.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. - 945./2. + 900.*((eta + 1.)*(eta + 1.)) - 600.*(eta + 1.)*(eta + 1.)*(eta + 1.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4., 0.);
                      case 1:
                        return sign * RealGradient(1720.*eta/9. + 385.*xi/18. - 1540.*(eta + 1.)*(xi + 1.)/9. + 1925.*(xi + 1.)*((eta + 1.)*(eta + 1.))/6. - 1925.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/9. + 13475.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/288. + 3395./18. - 1075.*(eta + 1.)*(eta + 1.)/3. + 2150.*((eta + 1.)*(eta + 1.)*(eta + 1.))/9. - 7525.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/144., 0.);
                      case 2:
                        return sign * RealGradient(-1360.*eta/9. - 385.*xi/18. + 1540.*(eta + 1.)*(xi + 1.)/9. - 1925.*(xi + 1.)*(eta + 1.)*(eta + 1.)/6. + 1925.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/9. - 13475.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/288. - 2765./18. + 850.*((eta + 1.)*(eta + 1.))/3. - 1700.*(eta + 1.)*(eta + 1.)*(eta + 1.)/9. + 2975.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/72., 0.);
                      case 3:
                        return sign * RealGradient(360.*eta + 105.*xi/2. - 420.*(eta + 1.)*(xi + 1.) + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 525.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 3675.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 735./2. - 675.*(eta + 1.)*(eta + 1.) + 450.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16., 0.);
                      case 4:
                        return sign * RealGradient(0., -450.*eta - 360.*xi + 1350.*(eta + 1.)*(xi + 1.) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 1350.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 690. + 450.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. + 210.*((xi + 1.)*(xi + 1.)) - 3675.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)/4.);
                      case 5:
                        return sign * RealGradient(0., 475.*eta/3. + 170.*xi/3. - 475.*(eta + 1.)*(xi + 1.) + 3325.*(eta + 1.)*((xi + 1.)*(xi + 1.))/12. + 1075.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 1925.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/12. + 1765./9. - 1075.*(eta + 1.)*(eta + 1.)/6. - 7525.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/24. - 595.*(xi + 1.)*(xi + 1.)/18. + 13475.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/144. + 1925.*((eta + 1.)*(eta + 1.)*(eta + 1.))/36.);
                      case 6:
                        return sign * RealGradient(0., -250.*eta/3. - 80.*xi/3. + 250.*(eta + 1.)*(xi + 1.) - 875.*(eta + 1.)*(xi + 1.)*(xi + 1.)/6. - 425.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 1925.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/12. - 910./9. + 425.*((eta + 1.)*(eta + 1.))/3. + 2975.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/12. + 140.*((xi + 1.)*(xi + 1.))/9. - 13475.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/144. - 1925.*(eta + 1.)*(eta + 1.)*(eta + 1.)/36.);
                      case 7:
                        return sign * RealGradient(0., 225.*eta + 90.*xi - 675.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. + 2025.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 285. - 675.*(eta + 1.)*(eta + 1.)/2. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 105.*(xi + 1.)*(xi + 1.)/2. + 3675.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.))/4.);
                      case 8:
                        return sign * RealGradient(-90.*eta + 105.*(eta + 1.)*(xi + 1.) - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)/4. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 3675.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. - 90. + 675.*((eta + 1.)*(eta + 1.))/2. - 675.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16., 0.);
                      case 9:
                        return sign * RealGradient(340.*eta/9. - 385.*(eta + 1.)*(xi + 1.)/9. + 1925.*(xi + 1.)*((eta + 1.)*(eta + 1.))/12. - 1925.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/12. + 13475.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/288. + 340./9. - 425.*(eta + 1.)*(eta + 1.)/3. + 425.*((eta + 1.)*(eta + 1.)*(eta + 1.))/3. - 2975.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/72., 0.);
                      case 10:
                        return sign * RealGradient(-430.*eta/9. + 385.*(eta + 1.)*(xi + 1.)/9. - 1925.*(xi + 1.)*(eta + 1.)*(eta + 1.)/12. + 1925.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/12. - 13475.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/288. - 430./9. + 1075.*((eta + 1.)*(eta + 1.))/6. - 1075.*(eta + 1.)*(eta + 1.)*(eta + 1.)/6. + 7525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/144., 0.);
                      case 11:
                        return sign * RealGradient(120.*eta - 105.*(eta + 1.)*(xi + 1.) + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 3675.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 120. - 450.*(eta + 1.)*(eta + 1.) + 450.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4., 0.);
                      case 12:
                        return sign * RealGradient(0., -450.*eta - 120.*xi + 900.*(eta + 1.)*(xi + 1.) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. - 1350.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 525.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 510. + 675.*((eta + 1.)*(eta + 1.)) + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. + 105.*((xi + 1.)*(xi + 1.))/2. - 3675.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2.);
                      case 13:
                        return sign * RealGradient(0., 500.*eta/3. + 320.*xi/9. - 1000.*(eta + 1.)*(xi + 1.)/3. + 875.*(eta + 1.)*((xi + 1.)*(xi + 1.))/6. + 1700.*(xi + 1.)*((eta + 1.)*(eta + 1.))/3. - 1925.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/9. + 1660./9. - 850.*(eta + 1.)*(eta + 1.)/3. - 2975.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/12. - 140.*(xi + 1.)*(xi + 1.)/9. + 13475.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/144. + 1925.*((eta + 1.)*(eta + 1.)*(eta + 1.))/18.);
                      case 14:
                        return sign * RealGradient(0., -950.*eta/3. - 680.*xi/9. + 1900.*(eta + 1.)*(xi + 1.)/3. - 3325.*(eta + 1.)*(xi + 1.)*(xi + 1.)/12. - 2150.*(xi + 1.)*(eta + 1.)*(eta + 1.)/3. + 1925.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/9. - 3190./9. + 1075.*((eta + 1.)*(eta + 1.))/3. + 7525.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/24. + 595.*((xi + 1.)*(xi + 1.))/18. - 13475.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/144. - 1925.*(eta + 1.)*(eta + 1.)*(eta + 1.)/18.);
                      case 15:
                        return sign * RealGradient(0., 900.*eta + 480.*xi - 1800.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. + 1800.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 525.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 1140. - 900.*(eta + 1.)*(eta + 1.) - 1575.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 210.*(xi + 1.)*(xi + 1.) + 3675.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.))/2.);
                      case 16:
                        return RealGradient(0., 870.*eta + 480.*xi - 1800.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. + 1800.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 525.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 1118. - 870.*(eta + 1.)*(eta + 1.) - 1575.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 210.*(xi + 1.)*(xi + 1.) + 3675.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 1015.*((eta + 1.)*(eta + 1.)*(eta + 1.))/4.);
                      case 17:
                        return RealGradient(0., 420.*eta + 360.*xi - 1350.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. + 1350.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 668. - 420.*(eta + 1.)*(eta + 1.) - 1575.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 210.*(xi + 1.)*(xi + 1.) + 3675.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 245.*((eta + 1.)*(eta + 1.)*(eta + 1.))/2.);
                      case 18:
                        return RealGradient(0., 60.*eta + 44. - 60.*(eta + 1.)*(eta + 1.) + 35.*((eta + 1.)*(eta + 1.)*(eta + 1.))/2.);
                      case 19:
                        return RealGradient(-390.*eta + 1365.*(eta + 1.)*(xi + 1.)/4. - 3045.*(xi + 1.)*(eta + 1.)*(eta + 1.)/4. + 525.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 3675.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. - 390. + 870.*((eta + 1.)*(eta + 1.)) - 600.*(eta + 1.)*(eta + 1.)*(eta + 1.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4., 0.);
                      case 20:
                        return RealGradient(-90.*eta + 315.*(eta + 1.)*(xi + 1.)/4. - 735.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 3675.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. - 90. + 420.*((eta + 1.)*(eta + 1.)) - 450.*(eta + 1.)*(eta + 1.)*(eta + 1.) + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4., 0.);
                      case 21:
                        return RealGradient(-120.*eta + 105.*(eta + 1.)*(xi + 1.) - 105.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 120. + 60.*((eta + 1.)*(eta + 1.)), 0.);
                      case 22:
                        return RealGradient(585.*eta/2. - 1365.*(eta + 1.)*(xi + 1.)/4. + 3045.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. - 525.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 3675.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 585./2. - 1305.*(eta + 1.)*(eta + 1.)/2. + 450.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16., 0.);
                      case 23:
                        return RealGradient(135.*eta/2. - 315.*(eta + 1.)*(xi + 1.)/4. + 735.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 3675.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 135./2. - 315.*(eta + 1.)*(eta + 1.) + 675.*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16., 0.);
                      case 24:
                        return RealGradient(90.*eta - 105.*(eta + 1.)*(xi + 1.) + 105.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. + 90. - 45.*(eta + 1.)*(eta + 1.), 0.);
                      case 25:
                        return RealGradient(0., -435.*eta - 120.*xi + 900.*(eta + 1.)*(xi + 1.) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. - 1350.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 525.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 497. + 1305.*((eta + 1.)*(eta + 1.))/2. + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. + 105.*((xi + 1.)*(xi + 1.))/2. - 3675.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 1015.*(eta + 1.)*(eta + 1.)*(eta + 1.)/4.);
                      case 26:
                        return RealGradient(0., -210.*eta - 90.*xi + 675.*(eta + 1.)*(xi + 1.) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. - 2025.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 272. + 315.*((eta + 1.)*(eta + 1.)) + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. + 105.*((xi + 1.)*(xi + 1.))/2. - 3675.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 245.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2.);
                      case 27:
                        return RealGradient(0., -30.*eta - 26. + 45.*((eta + 1.)*(eta + 1.)) - 35.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2.);
                      case 28:
                        return RealGradient(0., 0.);
                      case 29:
                        return RealGradient(0., 0.);
                      case 30:
                        return RealGradient(0., 171.*eta/4. + 120.*xi - 90.*(eta + 1.)*(xi + 1.) + 315.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. + 423./4. - 105.*(xi + 1.)*(xi + 1.)/2.);
                      case 31:
                        return RealGradient(0., -171.*eta/4. - 60.*xi + 90.*(eta + 1.)*(xi + 1.) - 315.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 297./4. + 105.*((xi + 1.)*(xi + 1.))/4.);
                      case 32:
                        return RealGradient(0., 81.*eta/4. + 90.*xi - 135.*(eta + 1.)*(xi + 1.)/2. + 315.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. + 333./4. - 105.*(xi + 1.)*(xi + 1.)/2.);
                      case 33:
                        return RealGradient(0., -81.*eta/4. - 45.*xi + 135.*(eta + 1.)*(xi + 1.)/2. - 315.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 207./4. + 105.*((xi + 1.)*(xi + 1.))/4.);
                      case 34:
                        return RealGradient(0., 0.);
                      case 35:
                        return RealGradient(0., 0.);
                      case 36:
                        return RealGradient(0., 9.*eta/2. - 3./2.);
                      case 37:
                        return RealGradient(0., 0.);
                      case 38:
                        return RealGradient(0., 0.);
                      case 39:
                        return RealGradient(0., -9.*eta/2. - 3./2.);
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
                        return sign * RealGradient(-900.*eta - 480.*xi + 7200.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 6300.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 14400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1140. + 3600.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 12600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 840.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7350.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 2100.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 1:
                        return sign * RealGradient(950.*eta/3. + 1720.*xi/9. - 8600.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 7700.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. + 17200.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. - 30100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 1270./3. - 3800.*(eta + 1.)*(eta + 1.)/(2.*2.)/3. - 15400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. - 3080.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 26950.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. + 6650.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9., 0.);
                      case 2:
                        return sign * RealGradient(-500.*eta/3. - 1360.*xi/9. + 6800.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 7700.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. - 13600.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. + 23800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 820./3. + 2000.*((eta + 1.)*(eta + 1.)/(2.*2.))/3. + 15400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. + 3080.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 26950.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. - 3500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9., 0.);
                      case 3:
                        return sign * RealGradient(450.*eta + 360.*xi - 5400.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 6300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 10800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 6300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 690. - 1800.*(eta + 1.)*(eta + 1.)/(2.*2.) - 12600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 840.*(xi + 1.)*(xi + 1.)/(2.*2.) + 7350.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 1050.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 4:
                        return sign * RealGradient(0., -120.*eta - 450.*xi + 3600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 10800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3150.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 510. + 210.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 9450.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7350.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2700.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2100.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 5:
                        return sign * RealGradient(0., 430.*eta/9. + 475.*xi/3. - 4300.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 4300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 30100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 3850.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. + 185. - 770.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 3850.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 26950.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 950.*(xi + 1.)*(xi + 1.)/(2.*2.) + 6650.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9.);
                      case 6:
                        return sign * RealGradient(0., -340.*eta/9. - 250.*xi/3. + 3400.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 3400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 23800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 3850.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. - 110. + 770.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 3850.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 26950.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 500.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3500.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9.);
                      case 7:
                        return sign * RealGradient(0., 90.*eta + 225.*xi - 2700.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 8100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 6300.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 3150.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 285. - 210.*(eta + 1.)*(eta + 1.)/(2.*2.) - 9450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7350.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1350.*(xi + 1.)*(xi + 1.)/(2.*2.) + 1050.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 8:
                        return sign * RealGradient(-225.*eta - 90.*xi + 2700.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3150.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 8100.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 6300.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 285. + 1350.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 9450.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 210.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7350.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 1050.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 9:
                        return sign * RealGradient(250.*eta/3. + 340.*xi/9. - 3400.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 3850.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. + 3400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 23800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. + 110. - 500.*(eta + 1.)*(eta + 1.)/(2.*2.) - 3850.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 770.*(xi + 1.)*(xi + 1.)/(2.*2.)/9. + 26950.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. + 3500.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9., 0.);
                      case 10:
                        return sign * RealGradient(-475.*eta/3. - 430.*xi/9. + 4300.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 3850.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. - 4300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 30100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/9. - 185. + 950.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 3850.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 770.*((xi + 1.)*(xi + 1.)/(2.*2.))/9. - 26950.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9. - 6650.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/9., 0.);
                      case 11:
                        return sign * RealGradient(450.*eta + 120.*xi - 3600.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 3150.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 10800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 8400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 510. - 2700.*(eta + 1.)*(eta + 1.)/(2.*2.) - 9450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 210.*(xi + 1.)*(xi + 1.)/(2.*2.) + 7350.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 2100.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 12:
                        return sign * RealGradient(0., -360.*eta - 450.*xi + 5400.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 10800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 6300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 690. + 840.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 12600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7350.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1800.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1050.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 13:
                        return sign * RealGradient(0., 1360.*eta/9. + 500.*xi/3. - 6800.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. + 13600.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 23800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 7700.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/3. + 820./3. - 3080.*(eta + 1.)*(eta + 1.)/(2.*2.)/9. - 15400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 26950.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 2000.*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 3500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9.);
                      case 14:
                        return sign * RealGradient(0., -1720.*eta/9. - 950.*xi/3. + 8600.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/3. - 17200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/3. + 30100.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/9. - 7700.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/3. - 1270./3. + 3080.*((eta + 1.)*(eta + 1.)/(2.*2.))/9. + 15400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 26950.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9. + 3800.*((xi + 1.)*(xi + 1.)/(2.*2.))/3. - 6650.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/9.);
                      case 15:
                        return sign * RealGradient(0., 480.*eta + 900.*xi - 7200.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 14400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 6300.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 1140. - 840.*(eta + 1.)*(eta + 1.)/(2.*2.) - 12600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7350.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3600.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 16:
                        return RealGradient(0., 390.*eta + 870.*xi - 6960.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 14400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 6090.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 1065. - 1365.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 12600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7350.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 3600.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 17:
                        return RealGradient(0., 90.*eta + 420.*xi - 3360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 10800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 8400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2940.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 465. - 315.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 9450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 7350.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2700.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 18:
                        return RealGradient(0., 120.*eta + 60.*xi - 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 420.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 120. - 210.*(eta + 1.)*(eta + 1.)/(2.*2.));
                      case 19:
                        return RealGradient(-870.*eta - 390.*xi + 6960.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 6090.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 14400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1065. + 3600.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 12600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 1365.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 7350.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 2100.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 20:
                        return RealGradient(-420.*eta - 90.*xi + 3360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 2940.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 10800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 465. + 2700.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 9450.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 315.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 7350.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) - 2100.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 21:
                        return RealGradient(-60.*eta - 120.*xi + 480.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 420.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 120. + 210.*((xi + 1.)*(xi + 1.)/(2.*2.)), 0.);
                      case 22:
                        return RealGradient(435.*eta + 585.*xi/2. - 5220.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 6090.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 10800.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 6300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 630. - 1800.*(eta + 1.)*(eta + 1.)/(2.*2.) - 12600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 1365.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 7350.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 1050.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 23:
                        return RealGradient(210.*eta + 135.*xi/2. - 2520.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2940.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 8100.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 6300.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 255. - 1350.*(eta + 1.)*(eta + 1.)/(2.*2.) - 9450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 315.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 7350.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 1050.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 24:
                        return RealGradient(30.*eta + 90.*xi - 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 420.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 90. - 210.*(xi + 1.)*(xi + 1.)/(2.*2.), 0.);
                      case 25:
                        return RealGradient(0., -585.*eta/2. - 435.*xi + 5220.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 10800.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 6090.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 630. + 1365.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 12600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7350.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1800.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1050.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 26:
                        return RealGradient(0., -135.*eta/2. - 210.*xi + 2520.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 8100.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2940.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 255. + 315.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 9450.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 7350.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1350.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 1050.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 27:
                        return RealGradient(0., -90.*eta - 30.*xi + 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 420.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 90. + 210.*((eta + 1.)*(eta + 1.)/(2.*2.)));
                      case 28:
                        return RealGradient(171.*eta/4. + 135./4. - 180.*(eta + 1.)*(eta + 1.)/(2.*2.) + 105.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 29:
                        return RealGradient(-171.*eta/4. - 135./4. + 180.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 105.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 30:
                        return RealGradient(0., 171.*xi/4. + 135./4. - 180.*(xi + 1.)*(xi + 1.)/(2.*2.) + 105.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 31:
                        return RealGradient(0., -171.*xi/4. - 135./4. + 180.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 105.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 32:
                        return RealGradient(0., 81.*xi/4. + 75./4. - 135.*(xi + 1.)*(xi + 1.)/(2.*2.) + 105.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)));
                      case 33:
                        return RealGradient(0., -81.*xi/4. - 75./4. + 135.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 105.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.));
                      case 34:
                        return RealGradient(81.*eta/4. + 75./4. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.) + 105.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)), 0.);
                      case 35:
                        return RealGradient(-81.*eta/4. - 75./4. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 105.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.), 0.);
                      case 36:
                        return RealGradient(0., 9.*xi/2.);
                      case 37:
                        return RealGradient(-9.*eta/2., 0.);
                      case 38:
                        return RealGradient(9.*eta/2., 0.);
                      case 39:
                        return RealGradient(0., -9.*xi/2.);
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
                        return sign * RealGradient(-480.*eta - 900.*xi + 1800.*(eta + 1.)*(xi + 1.) - 1800.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 525.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 1140. + 210.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 3675.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 900.*((xi + 1.)*(xi + 1.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 1:
                        return sign * RealGradient(680.*eta/9. + 950.*xi/3. - 1900.*(eta + 1.)*(xi + 1.)/3. + 2150.*(eta + 1.)*((xi + 1.)*(xi + 1.))/3. - 1925.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/9. + 3325.*(xi + 1.)*((eta + 1.)*(eta + 1.))/12. + 3190./9. - 595.*(eta + 1.)*(eta + 1.)/18. - 7525.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/24. + 13475.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/144. - 1075.*(xi + 1.)*(xi + 1.)/3. + 1925.*((xi + 1.)*(xi + 1.)*(xi + 1.))/18., 0.);
                      case 2:
                        return sign * RealGradient(-320.*eta/9. - 500.*xi/3. + 1000.*(eta + 1.)*(xi + 1.)/3. - 1700.*(eta + 1.)*(xi + 1.)*(xi + 1.)/3. + 1925.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/9. - 875.*(xi + 1.)*(eta + 1.)*(eta + 1.)/6. - 1660./9. + 140.*((eta + 1.)*(eta + 1.))/9. + 2975.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/12. - 13475.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/144. + 850.*((xi + 1.)*(xi + 1.))/3. - 1925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/18., 0.);
                      case 3:
                        return sign * RealGradient(120.*eta + 450.*xi - 900.*(eta + 1.)*(xi + 1.) + 1350.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 525.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. + 510. - 105.*(eta + 1.)*(eta + 1.)/2. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 3675.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 675.*(xi + 1.)*(xi + 1.) + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.))/2., 0.);
                      case 4:
                        return sign * RealGradient(0., -120.*xi + 105.*(eta + 1.)*(xi + 1.) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 3675.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 120. + 450.*((xi + 1.)*(xi + 1.)) - 450.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4.);
                      case 5:
                        return sign * RealGradient(0., 430.*xi/9. - 385.*(eta + 1.)*(xi + 1.)/9. + 1925.*(eta + 1.)*((xi + 1.)*(xi + 1.))/12. - 1925.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/12. + 13475.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/288. + 430./9. - 1075.*(xi + 1.)*(xi + 1.)/6. + 1075.*((xi + 1.)*(xi + 1.)*(xi + 1.))/6. - 7525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/144.);
                      case 6:
                        return sign * RealGradient(0., -340.*xi/9. + 385.*(eta + 1.)*(xi + 1.)/9. - 1925.*(eta + 1.)*(xi + 1.)*(xi + 1.)/12. + 1925.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/12. - 13475.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/288. - 340./9. + 425.*((xi + 1.)*(xi + 1.))/3. - 425.*(xi + 1.)*(xi + 1.)*(xi + 1.)/3. + 2975.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/72.);
                      case 7:
                        return sign * RealGradient(0., 90.*xi - 105.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 3675.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 90. - 675.*(xi + 1.)*(xi + 1.)/2. + 675.*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16.);
                      case 8:
                        return sign * RealGradient(-90.*eta - 225.*xi + 675.*(eta + 1.)*(xi + 1.) - 2025.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)/4. - 285. + 105.*((eta + 1.)*(eta + 1.))/2. + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 3675.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 675.*((xi + 1.)*(xi + 1.))/2. - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)/4., 0.);
                      case 9:
                        return sign * RealGradient(80.*eta/3. + 250.*xi/3. - 250.*(eta + 1.)*(xi + 1.) + 425.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 1925.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/12. + 875.*(xi + 1.)*((eta + 1.)*(eta + 1.))/6. + 910./9. - 140.*(eta + 1.)*(eta + 1.)/9. - 2975.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/12. + 13475.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/144. - 425.*(xi + 1.)*(xi + 1.)/3. + 1925.*((xi + 1.)*(xi + 1.)*(xi + 1.))/36., 0.);
                      case 10:
                        return sign * RealGradient(-170.*eta/3. - 475.*xi/3. + 475.*(eta + 1.)*(xi + 1.) - 1075.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 1925.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/12. - 3325.*(xi + 1.)*(eta + 1.)*(eta + 1.)/12. - 1765./9. + 595.*((eta + 1.)*(eta + 1.))/18. + 7525.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/24. - 13475.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/144. + 1075.*((xi + 1.)*(xi + 1.))/6. - 1925.*(xi + 1.)*(xi + 1.)*(xi + 1.)/36., 0.);
                      case 11:
                        return sign * RealGradient(360.*eta + 450.*xi - 1350.*(eta + 1.)*(xi + 1.) + 1350.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. + 690. - 210.*(eta + 1.)*(eta + 1.) - 1575.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 3675.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 450.*(xi + 1.)*(xi + 1.) + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.))/4., 0.);
                      case 12:
                        return sign * RealGradient(0., -105.*eta/2. - 360.*xi + 420.*(eta + 1.)*(xi + 1.) - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 525.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 3675.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 735./2. + 675.*((xi + 1.)*(xi + 1.)) - 450.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/16.);
                      case 13:
                        return sign * RealGradient(0., 385.*eta/18. + 1360.*xi/9. - 1540.*(eta + 1.)*(xi + 1.)/9. + 1925.*(eta + 1.)*((xi + 1.)*(xi + 1.))/6. - 1925.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/9. + 13475.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/288. + 2765./18. - 850.*(xi + 1.)*(xi + 1.)/3. + 1700.*((xi + 1.)*(xi + 1.)*(xi + 1.))/9. - 2975.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/72.);
                      case 14:
                        return sign * RealGradient(0., -385.*eta/18. - 1720.*xi/9. + 1540.*(eta + 1.)*(xi + 1.)/9. - 1925.*(eta + 1.)*(xi + 1.)*(xi + 1.)/6. + 1925.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/9. - 13475.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/288. - 3395./18. + 1075.*((xi + 1.)*(xi + 1.))/3. - 2150.*(xi + 1.)*(xi + 1.)*(xi + 1.)/9. + 7525.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/144.);
                      case 15:
                        return sign * RealGradient(0., 105.*eta/2. + 480.*xi - 420.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 525.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 3675.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 945./2. - 900.*(xi + 1.)*(xi + 1.) + 600.*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4.);
                      case 16:
                        return RealGradient(0., 390.*xi - 1365.*(eta + 1.)*(xi + 1.)/4. + 3045.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. - 525.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 3675.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 390. - 870.*(xi + 1.)*(xi + 1.) + 600.*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4.);
                      case 17:
                        return RealGradient(0., 90.*xi - 315.*(eta + 1.)*(xi + 1.)/4. + 735.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 3675.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 90. - 420.*(xi + 1.)*(xi + 1.) + 450.*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4.);
                      case 18:
                        return RealGradient(0., 120.*xi - 105.*(eta + 1.)*(xi + 1.) + 105.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. + 120. - 60.*(xi + 1.)*(xi + 1.));
                      case 19:
                        return RealGradient(-480.*eta - 870.*xi + 1800.*(eta + 1.)*(xi + 1.) - 1800.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 525.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 1118. + 210.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 3675.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 870.*((xi + 1.)*(xi + 1.)) - 1015.*(xi + 1.)*(xi + 1.)*(xi + 1.)/4., 0.);
                      case 20:
                        return RealGradient(-360.*eta - 420.*xi + 1350.*(eta + 1.)*(xi + 1.) - 1350.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 1575.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 668. + 210.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 3675.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 420.*((xi + 1.)*(xi + 1.)) - 245.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 21:
                        return RealGradient(-60.*xi - 44. + 60.*((xi + 1.)*(xi + 1.)) - 35.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 22:
                        return RealGradient(120.*eta + 435.*xi - 900.*(eta + 1.)*(xi + 1.) + 1350.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 525.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. + 497. - 105.*(eta + 1.)*(eta + 1.)/2. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 3675.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 1305.*(xi + 1.)*(xi + 1.)/2. + 1015.*((xi + 1.)*(xi + 1.)*(xi + 1.))/4., 0.);
                      case 23:
                        return RealGradient(90.*eta + 210.*xi - 675.*(eta + 1.)*(xi + 1.) + 2025.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 1575.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. + 272. - 105.*(eta + 1.)*(eta + 1.)/2. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 3675.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 315.*(xi + 1.)*(xi + 1.) + 245.*((xi + 1.)*(xi + 1.)*(xi + 1.))/2., 0.);
                      case 24:
                        return RealGradient(30.*xi + 26. - 45.*(xi + 1.)*(xi + 1.) + 35.*((xi + 1.)*(xi + 1.)*(xi + 1.))/2., 0.);
                      case 25:
                        return RealGradient(0., -585.*xi/2. + 1365.*(eta + 1.)*(xi + 1.)/4. - 3045.*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. + 525.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 3675.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 585./2. + 1305.*((xi + 1.)*(xi + 1.))/2. - 450.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/16.);
                      case 26:
                        return RealGradient(0., -135.*xi/2. + 315.*(eta + 1.)*(xi + 1.)/4. - 735.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 3675.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 135./2. + 315.*((xi + 1.)*(xi + 1.)) - 675.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/16.);
                      case 27:
                        return RealGradient(0., -90.*xi + 105.*(eta + 1.)*(xi + 1.) - 105.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 90. + 45.*((xi + 1.)*(xi + 1.)));
                      case 28:
                        return RealGradient(120.*eta + 171.*xi/4. - 90.*(eta + 1.)*(xi + 1.) + 315.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 423./4. - 105.*(eta + 1.)*(eta + 1.)/2., 0.);
                      case 29:
                        return RealGradient(-60.*eta - 171.*xi/4. + 90.*(eta + 1.)*(xi + 1.) - 315.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 297./4. + 105.*((eta + 1.)*(eta + 1.))/4., 0.);
                      case 30:
                        return RealGradient(0., 0.);
                      case 31:
                        return RealGradient(0., 0.);
                      case 32:
                        return RealGradient(0., 0.);
                      case 33:
                        return RealGradient(0., 0.);
                      case 34:
                        return RealGradient(90.*eta + 81.*xi/4. - 135.*(eta + 1.)*(xi + 1.)/2. + 315.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 333./4. - 105.*(eta + 1.)*(eta + 1.)/2., 0.);
                      case 35:
                        return RealGradient(-45.*eta - 81.*xi/4. + 135.*(eta + 1.)*(xi + 1.)/2. - 315.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 207./4. + 105.*((eta + 1.)*(eta + 1.))/4., 0.);
                      case 36:
                        return RealGradient(0., 0.);
                      case 37:
                        return RealGradient(3./2. - 9.*xi/2., 0.);
                      case 38:
                        return RealGradient(9.*xi/2. + 3./2., 0.);
                      case 39:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(1344.*eta*xi - 1680.*eta - 840.*xi + 480. + 1344.*(eta*eta), -4032.*eta*xi + 1680.*eta + 2520.*xi - 480. - 1344.*eta*eta - 2688.*xi*xi);
                      case 1:
                        return sign * RealGradient(-4928.*eta*xi/9. + 5768.*eta/9. + 3080.*xi/9. - 1720./9. - 4480.*eta*eta/9., 4480.*eta*xi/3. - 4256.*eta/9. - 1008.*xi + 1600./9. + 1792.*(eta*eta)/9. + 9856.*(xi*xi)/9.);
                      case 2:
                        return sign * RealGradient(4928.*eta*xi/9. - 2240.*eta/9. - 3080.*xi/9. + 1360./9. + 448.*(eta*eta)/9., -448.*eta*xi/3. - 784.*eta/9. + 840.*xi - 1024./9. + 2240.*(eta*eta)/9. - 9856.*xi*xi/9.);
                      case 3:
                        return sign * RealGradient(-1344.*eta*xi + 504.*eta + 840.*xi - 360., -2016.*xi + 288. + 2688.*(xi*xi));
                      case 4:
                        return sign * RealGradient(168.*eta*(3. - 8.*xi), -2016.*xi + 288. + 2688.*(xi*xi));
                      case 5:
                        return sign * RealGradient(896.*eta*(-6.*eta - 4.*xi + 3.)/9., 1792.*eta*xi - 1792.*eta/3. - 2464.*xi/3. + 160. + 896.*(eta*eta)/3. + 7168.*(xi*xi)/9.);
                      case 6:
                        return sign * RealGradient(56.*eta*(-48.*eta - 8.*xi + 15.)/9., 896.*eta*xi - 1568.*eta/3. - 560.*xi/3. + 64. + 1792.*(eta*eta)/3. + 896.*(xi*xi)/9.);
                      case 7:
                        return sign * RealGradient(0., 0.);
                      case 8:
                        return sign * RealGradient(0., 0.);
                      case 9:
                        return sign * RealGradient(56.*eta*(-40.*eta + 8.*xi + 7.)/9., 2240.*eta*xi/3. - 952.*eta/3. + 112.*xi/9. + 208./9. - 448.*eta*eta/9. - 896.*xi*xi/9.);
                      case 10:
                        return sign * RealGradient(896.*eta*(-2.*eta + 4.*xi - 1.)/9., 1792.*eta*xi/3. - 1568.*eta/3. + 6944.*xi/9. - 1216./9. + 4480.*(eta*eta)/9. - 7168.*xi*xi/9.);
                      case 11:
                        return sign * RealGradient(168.*eta*(8.*eta + 8.*xi - 5.), -4032.*eta*xi + 2520.*eta + 3360.*xi - 960. - 1344.*eta*eta - 2688.*xi*xi);
                      case 12:
                        return RealGradient(2016.*eta*(-4.*eta - 2.*xi + 3.), 24192.*eta*xi - 12096.*eta - 9072.*xi + 2160. + 12096.*(eta*eta) + 8064.*(xi*xi));
                      case 13:
                        return RealGradient(1008.*eta*(12.*eta + 16.*xi - 9.), -36288.*eta*xi + 18144.*eta + 36288.*xi - 8640. - 8064.*eta*eta - 32256.*xi*xi);
                      case 14:
                        return RealGradient(1008.*eta*(-8.*eta - 12.*xi + 9.), 24192.*eta*xi - 6048.*eta - 21168.*xi + 3456. + 24192.*(xi*xi));
                      case 15:
                        return RealGradient(1008.*eta*(4.*eta + 16.*xi - 7.), -12096.*eta*xi + 3024.*eta + 28224.*xi - 4608. - 32256.*xi*xi);
                      case 16:
                        return RealGradient(0., -2016.*eta + 144. + 4032.*(eta*eta));
                      case 17:
                        return RealGradient(0., 4032.*eta - 288. - 8064.*eta*eta);
                      case 18:
                        return RealGradient(252.*eta*(-24.*eta - 12.*xi + 13.), 18144.*eta*xi - 8064.*eta - 6804.*xi + 1548. + 7056.*(eta*eta) + 6048.*(xi*xi));
                      case 19:
                        return RealGradient(252.*eta*(28.*eta + 16.*xi - 13.), -21168.*eta*xi + 9576.*eta + 9324.*xi - 2196. - 6048.*eta*eta - 8064.*xi*xi);
                      case 20:
                        return RealGradient(1008.*eta*(1. - xi), 1008.*eta - 1512.*xi + 144. - 2016.*eta*eta + 2016.*(xi*xi));
                      case 21:
                        return RealGradient(252.*eta*(-12.*eta + 16.*xi - 3.), 9072.*eta*xi - 5544.*eta + 6804.*xi - 936. + 4032.*(eta*eta) - 8064.*xi*xi);
                      case 22:
                        return RealGradient(1008.*eta*(4.*eta + 2.*xi - 3.), -12096.*eta*xi + 4536.*eta + 4536.*xi - 972. - 3024.*eta*eta - 4032.*xi*xi);
                      case 23:
                        return RealGradient(2016.*eta*(-eta - 2.*xi + 1.), 6048.*eta*xi - 1512.*eta - 8064.*xi + 1440. + 8064.*(xi*xi));
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
                        return sign * RealGradient(2688.*eta*xi - 2520.*eta - 1680.*xi + 720. + 2016.*(eta*eta) + 672.*(xi*xi), -2688.*eta*xi + 840.*eta + 1680.*xi - 240. - 672.*eta*eta - 2016.*xi*xi);
                      case 1:
                        return sign * RealGradient(-8960.*eta*xi/9. + 1568.*eta/3. + 5768.*xi/9. - 200. - 896.*eta*eta/3. - 2464.*xi*xi/9., 3584.*eta*xi/9. + 896.*eta/9. - 4256.*xi/9. + 176./9. - 1792.*eta*eta/9. + 2240.*(xi*xi)/3.);
                      case 2:
                        return sign * RealGradient(896.*eta*xi/9. + 952.*eta/3. - 2240.*xi/9. - 16. - 1120.*eta*eta/3. + 2464.*(xi*xi)/9., 4480.*eta*xi/9. - 392.*eta/9. - 784.*xi/9. + 112./9. - 224.*eta*eta/9. - 224.*xi*xi/3.);
                      case 3:
                        return sign * RealGradient(504.*xi - 72. - 672.*xi*xi, 0.);
                      case 4:
                        return sign * RealGradient(504.*xi - 72. - 672.*xi*xi, 0.);
                      case 5:
                        return sign * RealGradient(-3584.*eta*xi/3. + 1568.*eta/3. + 896.*xi/3. - 80. - 448.*eta*eta - 1792.*xi*xi/9., 1792.*eta*xi/3. - 280.*eta/3. - 1792.*xi/3. + 56. + 224.*(eta*eta)/9. + 896.*(xi*xi));
                      case 6:
                        return sign * RealGradient(-1792.*eta*xi/3. + 1792.*eta/3. + 280.*xi/3. - 56. - 896.*eta*eta - 224.*xi*xi/9., 3584.*eta*xi/3. - 896.*eta/3. - 1568.*xi/3. + 80. + 1792.*(eta*eta)/9. + 448.*(xi*xi));
                      case 7:
                        return sign * RealGradient(0., -504.*eta + 72. + 672.*(eta*eta));
                      case 8:
                        return sign * RealGradient(0., -504.*eta + 72. + 672.*(eta*eta));
                      case 9:
                        return sign * RealGradient(-4480.*eta*xi/9. + 784.*eta/9. + 392.*xi/9. - 112./9. + 224.*(eta*eta)/3. + 224.*(xi*xi)/9., -896.*eta*xi/9. + 2240.*eta/9. - 952.*xi/3. + 16. - 2464.*eta*eta/9. + 1120.*(xi*xi)/3.);
                      case 10:
                        return sign * RealGradient(-3584.*eta*xi/9. + 4256.*eta/9. - 896.*xi/9. - 176./9. - 2240.*eta*eta/3. + 1792.*(xi*xi)/9., 8960.*eta*xi/9. - 5768.*eta/9. - 1568.*xi/3. + 200. + 2464.*(eta*eta)/9. + 896.*(xi*xi)/3.);
                      case 11:
                        return sign * RealGradient(2688.*eta*xi - 1680.*eta - 840.*xi + 240. + 2016.*(eta*eta) + 672.*(xi*xi), -2688.*eta*xi + 1680.*eta + 2520.*xi - 720. - 672.*eta*eta - 2016.*xi*xi);
                      case 12:
                        return RealGradient(-16128.*eta*xi + 18144.*eta + 6048.*xi - 3240. - 18144.*eta*eta - 2016.*xi*xi, 24192.*eta*xi - 9072.*eta - 12096.*xi + 2160. + 8064.*(eta*eta) + 12096.*(xi*xi));
                      case 13:
                        return RealGradient(24192.*eta*xi - 12096.*eta - 9072.*xi + 2160. + 12096.*(eta*eta) + 8064.*(xi*xi), -16128.*eta*xi + 6048.*eta + 18144.*xi - 3240. - 2016.*eta*eta - 18144.*xi*xi);
                      case 14:
                        return RealGradient(-16128.*eta*xi + 4032.*eta + 9072.*xi - 1944. - 6048.*xi*xi, -6048.*xi + 432. + 12096.*(xi*xi));
                      case 15:
                        return RealGradient(8064.*eta*xi - 2016.*eta - 7056.*xi + 1152. + 8064.*(xi*xi), 3024.*xi - 216. - 6048.*xi*xi);
                      case 16:
                        return RealGradient(3024.*eta - 216. - 6048.*eta*eta, 8064.*eta*xi - 7056.*eta - 2016.*xi + 1152. + 8064.*(eta*eta));
                      case 17:
                        return RealGradient(-6048.*eta + 432. + 12096.*(eta*eta), -16128.*eta*xi + 9072.*eta + 4032.*xi - 1944. - 6048.*eta*eta);
                      case 18:
                        return RealGradient(-12096.*eta*xi + 9576.*eta + 3276.*xi - 1332. - 10584.*eta*eta - 1512.*xi*xi, 14112.*eta*xi - 3276.*eta - 8064.*xi + 1044. + 2016.*(eta*eta) + 9072.*(xi*xi));
                      case 19:
                        return RealGradient(14112.*eta*xi - 8064.*eta - 3276.*xi + 1044. + 9072.*(eta*eta) + 2016.*(xi*xi), -12096.*eta*xi + 3276.*eta + 9576.*xi - 1332. - 1512.*eta*eta - 10584.*xi*xi);
                      case 20:
                        return RealGradient(-1512.*eta + 1008.*xi - 216. + 3024.*(eta*eta) - 504.*xi*xi, -4032.*eta*xi + 2016.*eta + 1008.*xi - 360. - 2016.*eta*eta);
                      case 21:
                        return RealGradient(-6048.*eta*xi + 4536.*eta - 756.*xi - 216. - 6048.*eta*eta + 2016.*(xi*xi), 8064.*eta*xi - 3024.*eta - 5544.*xi + 1188. + 1008.*(eta*eta) + 4536.*(xi*xi));
                      case 22:
                        return RealGradient(8064.*eta*xi - 5544.*eta - 3024.*xi + 1188. + 4536.*(eta*eta) + 1008.*(xi*xi), -6048.*eta*xi - 756.*eta + 4536.*xi - 216. + 2016.*(eta*eta) - 6048.*xi*xi);
                      case 23:
                        return RealGradient(-4032.*eta*xi + 1008.*eta + 2016.*xi - 360. - 2016.*xi*xi, 1008.*eta - 1512.*xi - 216. - 504.*eta*eta + 3024.*(xi*xi));
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
                        return sign * RealGradient(4032.*eta*xi - 3360.*eta - 2520.*xi + 960. + 2688.*(eta*eta) + 1344.*(xi*xi), 168.*xi*(-8.*eta - 8.*xi + 5.));
                      case 1:
                        return sign * RealGradient(-1792.*eta*xi/3. - 6944.*eta/9. + 1568.*xi/3. + 1216./9. + 7168.*(eta*eta)/9. - 4480.*xi*xi/9., 896.*xi*(-4.*eta + 2.*xi + 1.)/9.);
                      case 2:
                        return sign * RealGradient(-2240.*eta*xi/3. - 112.*eta/9. + 952.*xi/3. - 208./9. + 896.*(eta*eta)/9. + 448.*(xi*xi)/9., 56.*xi*(-8.*eta + 40.*xi - 7.)/9.);
                      case 3:
                        return sign * RealGradient(0., 0.);
                      case 4:
                        return sign * RealGradient(0., 0.);
                      case 5:
                        return sign * RealGradient(-896.*eta*xi + 560.*eta/3. + 1568.*xi/3. - 64. - 896.*eta*eta/9. - 1792.*xi*xi/3., 56.*xi*(8.*eta + 48.*xi - 15.)/9.);
                      case 6:
                        return sign * RealGradient(-1792.*eta*xi + 2464.*eta/3. + 1792.*xi/3. - 160. - 7168.*eta*eta/9. - 896.*xi*xi/3., 896.*xi*(4.*eta + 6.*xi - 3.)/9.);
                      case 7:
                        return sign * RealGradient(2016.*eta - 288. - 2688.*eta*eta, 168.*xi*(8.*eta - 3.));
                      case 8:
                        return sign * RealGradient(2016.*eta - 288. - 2688.*eta*eta, 1344.*eta*xi - 840.*eta - 504.*xi + 360.);
                      case 9:
                        return sign * RealGradient(448.*eta*xi/3. - 840.*eta + 784.*xi/9. + 1024./9. + 9856.*(eta*eta)/9. - 2240.*xi*xi/9., -4928.*eta*xi/9. + 3080.*eta/9. + 2240.*xi/9. - 1360./9. - 448.*xi*xi/9.);
                      case 10:
                        return sign * RealGradient(-4480.*eta*xi/3. + 1008.*eta + 4256.*xi/9. - 1600./9. - 9856.*eta*eta/9. - 1792.*xi*xi/9., 4928.*eta*xi/9. - 3080.*eta/9. - 5768.*xi/9. + 1720./9. + 4480.*(xi*xi)/9.);
                      case 11:
                        return sign * RealGradient(4032.*eta*xi - 2520.*eta - 1680.*xi + 480. + 2688.*(eta*eta) + 1344.*(xi*xi), -1344.*eta*xi + 840.*eta + 1680.*xi - 480. - 1344.*xi*xi);
                      case 12:
                        return RealGradient(-36288.*eta*xi + 36288.*eta + 18144.*xi - 8640. - 32256.*eta*eta - 8064.*xi*xi, 1008.*xi*(16.*eta + 12.*xi - 9.));
                      case 13:
                        return RealGradient(24192.*eta*xi - 9072.*eta - 12096.*xi + 2160. + 8064.*(eta*eta) + 12096.*(xi*xi), 2016.*xi*(-2.*eta - 4.*xi + 3.));
                      case 14:
                        return RealGradient(4032.*xi - 288. - 8064.*xi*xi, 0.);
                      case 15:
                        return RealGradient(-2016.*xi + 144. + 4032.*(xi*xi), 0.);
                      case 16:
                        return RealGradient(-12096.*eta*xi + 28224.*eta + 3024.*xi - 4608. - 32256.*eta*eta, 1008.*xi*(16.*eta + 4.*xi - 7.));
                      case 17:
                        return RealGradient(24192.*eta*xi - 21168.*eta - 6048.*xi + 3456. + 24192.*(eta*eta), 1008.*xi*(-12.*eta - 8.*xi + 9.));
                      case 18:
                        return RealGradient(-21168.*eta*xi + 9324.*eta + 9576.*xi - 2196. - 8064.*eta*eta - 6048.*xi*xi, 252.*xi*(16.*eta + 28.*xi - 13.));
                      case 19:
                        return RealGradient(18144.*eta*xi - 6804.*eta - 8064.*xi + 1548. + 6048.*(eta*eta) + 7056.*(xi*xi), 252.*xi*(-12.*eta - 24.*xi + 13.));
                      case 20:
                        return RealGradient(6048.*eta*xi - 8064.*eta - 1512.*xi + 1440. + 8064.*(eta*eta), 2016.*xi*(-2.*eta - xi + 1.));
                      case 21:
                        return RealGradient(-12096.*eta*xi + 4536.*eta + 4536.*xi - 972. - 4032.*eta*eta - 3024.*xi*xi, 1008.*xi*(2.*eta + 4.*xi - 3.));
                      case 22:
                        return RealGradient(9072.*eta*xi + 6804.*eta - 5544.*xi - 936. - 8064.*eta*eta + 4032.*(xi*xi), 252.*xi*(16.*eta - 12.*xi - 3.));
                      case 23:
                        return RealGradient(-1512.*eta + 1008.*xi + 144. + 2016.*(eta*eta) - 2016.*xi*xi, 1008.*xi*(1. - eta));
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
                        return sign * RealGradient(-13125.*eta/4. - 525.*xi + 13125.*(eta + 1.)*(xi + 1.)/2. - 23625.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 39375.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 91875.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 91875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 33075.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 14175./4. + 39375.*((eta + 1.)*(eta + 1.))/4. + 70875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. + 945.*((xi + 1.)*(xi + 1.))/4. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 165375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 59535.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 91875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 91875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 33075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 1:
                        return sign * RealGradient(527625.*eta/512. + 40845.*xi/256. - 1021125.*(eta + 1.)*(xi + 1.)/512. + 874125.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. + 3063375.*(xi + 1.)*((eta + 1.)*(eta + 1.))/512. - 7147875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 7147875.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/2048. - 2573235.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. + 567105./512. - 1582875.*(eta + 1.)*(eta + 1.)/512. - 2622375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/1024. - 34965.*(xi + 1.)*(xi + 1.)/512. + 6118875.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. - 6118875.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. + 2202795.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8192. + 3693375.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 3693375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2048. + 1329615.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4096., 0.);
                      case 2:
                        return sign * RealGradient(-28875.*eta/32. - 2835.*xi/16. + 70875.*(eta + 1.)*(xi + 1.)/32. - 70875.*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. - 212625.*(xi + 1.)*(eta + 1.)*(eta + 1.)/32. + 496125.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 496125.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. + 178605.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/256. - 32235./32. + 86625.*((eta + 1.)*(eta + 1.))/32. + 212625.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/64. + 2835.*((xi + 1.)*(xi + 1.))/32. - 496125.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. + 496125.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/256. - 178605.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. - 202125.*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 202125.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/128. - 72765.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/256., 0.);
                      case 3:
                        return sign * RealGradient(233625.*eta/512. + 29085.*xi/256. - 727125.*(eta + 1.)*(xi + 1.)/512. + 874125.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. + 2181375.*(xi + 1.)*((eta + 1.)*(eta + 1.))/512. - 5089875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 5089875.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/2048. - 1832355.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. + 273105./512. - 700875.*(eta + 1.)*(eta + 1.)/512. - 2622375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/1024. - 34965.*(xi + 1.)*(xi + 1.)/512. + 6118875.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. - 6118875.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. + 2202795.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8192. + 1635375.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 1635375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2048. + 588735.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4096., 0.);
                      case 4:
                        return sign * RealGradient(-7875.*eta/4. - 420.*xi + 5250.*(eta + 1.)*(xi + 1.) - 23625.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 15750.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 18375.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 18375.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 6615.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 8925./4. + 23625.*((eta + 1.)*(eta + 1.))/4. + 70875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. + 945.*((xi + 1.)*(xi + 1.))/4. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 165375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 59535.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 55125.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 55125.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 19845.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 5:
                        return sign * RealGradient(0., 2250.*eta + 7875.*xi/4. - 23625.*(eta + 1.)*(xi + 1.)/2. + 15750.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 23625.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 165375.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. - 55125.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 99225.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 15375./4. - 7875.*(eta + 1.)*(eta + 1.)/2. - 55125.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 2625.*(xi + 1.)*(xi + 1.) + 18375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 33075.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 2625.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 55125.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 7875.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8.);
                      case 6:
                        return sign * RealGradient(0., -42975.*eta/64. - 152775.*xi/512. + 902475.*(eta + 1.)*(xi + 1.)/256. - 300825.*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. + 902475.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 6648075.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. + 4288725.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 3671325.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. - 467475./512. + 316575.*((eta + 1.)*(eta + 1.))/256. + 2216025.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/256. - 6648075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 50925.*((xi + 1.)*(xi + 1.))/128. - 1429575.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/256. + 1223775.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 204225.*(eta + 1.)*(eta + 1.)*(eta + 1.)/256. + 4288725.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2048. - 152775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. - 3671325.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8192. + 174825.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/1024.);
                      case 7:
                        return sign * RealGradient(0., 1575.*eta/4. + 4725.*xi/32. - 33075.*(eta + 1.)*(xi + 1.)/16. + 11025.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. - 33075.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 363825.*(xi + 1.)*((eta + 1.)*(eta + 1.))/64. - 297675.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 297675.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/256. + 16425./32. - 17325.*(eta + 1.)*(eta + 1.)/16. - 121275.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/16. + 363825.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 1575.*(xi + 1.)*(xi + 1.)/8. + 99225.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 99225.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 14175.*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. + 297675.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/512. - 14175.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64.);
                      case 8:
                        return sign * RealGradient(0., -10575.*eta/64. - 26775.*xi/512. + 222075.*(eta + 1.)*(xi + 1.)/256. - 74025.*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. + 222075.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 2943675.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. + 3053925.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 3671325.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. - 106275./512. + 140175.*((eta + 1.)*(eta + 1.))/256. + 981225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/256. - 2943675.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 8925.*((xi + 1.)*(xi + 1.))/128. - 1017975.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/256. + 1223775.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 145425.*(eta + 1.)*(eta + 1.)*(eta + 1.)/256. + 3053925.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2048. - 26775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. - 3671325.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8192. + 174825.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/1024.);
                      case 9:
                        return sign * RealGradient(0., 900.*eta + 1575.*xi/4. - 4725.*(eta + 1.)*(xi + 1.) + 6300.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 99225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. - 11025.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 99225.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 4875./4. - 4725.*(eta + 1.)*(eta + 1.)/2. - 33075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 99225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 525.*(xi + 1.)*(xi + 1.) + 14700.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 33075.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 2100.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 11025.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8.);
                      case 10:
                        return sign * RealGradient(-1575.*eta/4. + 1050.*(eta + 1.)*(xi + 1.) - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 6300.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 11025.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 7350.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 6615.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 1575./4. + 4725.*((eta + 1.)*(eta + 1.))/2. + 14175.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/4. - 99225.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 33075.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 59535.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 33075.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 11025.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 19845.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 11:
                        return sign * RealGradient(46725.*eta/512. - 145425.*(eta + 1.)*(xi + 1.)/512. + 174825.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. + 436275.*(xi + 1.)*((eta + 1.)*(eta + 1.))/256. - 3053925.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 1017975.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/512. - 1832355.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. + 46725./512. - 140175.*(eta + 1.)*(eta + 1.)/256. - 524475.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/512. + 3671325.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. - 1223775.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 2202795.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8192. + 981225.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 327075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 588735.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4096., 0.);
                      case 12:
                        return sign * RealGradient(-5775.*eta/32. + 14175.*(eta + 1.)*(xi + 1.)/32. - 14175.*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. - 42525.*(xi + 1.)*(eta + 1.)*(eta + 1.)/16. + 297675.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 99225.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 178605.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/256. - 5775./32. + 17325.*((eta + 1.)*(eta + 1.))/16. + 42525.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/32. - 297675.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. + 99225.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 178605.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. - 121275.*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 40425.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 72765.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/256., 0.);
                      case 13:
                        return sign * RealGradient(105525.*eta/512. - 204225.*(eta + 1.)*(xi + 1.)/512. + 174825.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. + 612675.*(xi + 1.)*((eta + 1.)*(eta + 1.))/256. - 4288725.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 1429575.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/512. - 2573235.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. + 105525./512. - 316575.*(eta + 1.)*(eta + 1.)/256. - 524475.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/512. + 3671325.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. - 1223775.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 2202795.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8192. + 2216025.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 738675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 1329615.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4096., 0.);
                      case 14:
                        return sign * RealGradient(-2625.*eta/4. + 2625.*(eta + 1.)*(xi + 1.)/2. - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 7875.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 55125.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 18375.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 33075.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 2625./4. + 7875.*((eta + 1.)*(eta + 1.))/2. + 14175.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/4. - 99225.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 33075.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 59535.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 55125.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 18375.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 33075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 15:
                        return sign * RealGradient(0., 2250.*eta + 2625.*xi/4. - 7875.*(eta + 1.)*(xi + 1.) + 7875.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 165375.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. - 18375.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 165375.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 10875./4. - 23625.*(eta + 1.)*(eta + 1.)/4. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 99225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 2625.*(xi + 1.)*(xi + 1.)/4. + 18375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 5250.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 11025.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 23625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16.);
                      case 16:
                        return sign * RealGradient(0., -52875.*eta/128. - 44625.*xi/512. + 370125.*(eta + 1.)*(xi + 1.)/256. - 370125.*(eta + 1.)*(xi + 1.)*(xi + 1.)/256. + 222075.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 4906125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. + 5089875.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 6118875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. - 243375./512. + 700875.*((eta + 1.)*(eta + 1.))/512. + 4906125.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/1024. - 2943675.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 44625.*((xi + 1.)*(xi + 1.))/512. - 5089875.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 6118875.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4096. - 727125.*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 3053925.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2048. - 26775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. - 3671325.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8192. + 874125.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/2048.);
                      case 17:
                        return sign * RealGradient(0., 7875.*eta/8. + 7875.*xi/32. - 55125.*(eta + 1.)*(xi + 1.)/16. + 55125.*(eta + 1.)*((xi + 1.)*(xi + 1.))/16. - 33075.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 606375.*(xi + 1.)*((eta + 1.)*(eta + 1.))/64. - 496125.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 496125.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/256. + 37125./32. - 86625.*(eta + 1.)*(eta + 1.)/32. - 606375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. + 363825.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 7875.*(xi + 1.)*(xi + 1.)/32. + 496125.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 496125.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/256. + 70875.*((eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. + 297675.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/512. - 70875.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128.);
                      case 18:
                        return sign * RealGradient(0., -214875.*eta/128. - 254625.*xi/512. + 1504125.*(eta + 1.)*(xi + 1.)/256. - 1504125.*(eta + 1.)*(xi + 1.)*(xi + 1.)/256. + 902475.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 11080125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. + 7147875.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 6118875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4096. - 1041375./512. + 1582875.*((eta + 1.)*(eta + 1.))/512. + 11080125.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/1024. - 6648075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 254625.*((xi + 1.)*(xi + 1.))/512. - 7147875.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/1024. + 6118875.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4096. - 1021125.*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 4288725.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2048. - 152775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. - 3671325.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8192. + 874125.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/2048.);
                      case 19:
                        return sign * RealGradient(0., 5625.*eta + 13125.*xi/4. - 39375.*(eta + 1.)*(xi + 1.)/2. + 39375.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 23625.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 275625.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. - 91875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 165375.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 31875./4. - 39375.*(eta + 1.)*(eta + 1.)/4. - 275625.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 13125.*(xi + 1.)*(xi + 1.)/4. + 91875.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 13125.*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 55125.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 7875.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 23625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16.);
                      case 20:
                        return RealGradient(0., 10575.*eta/2. + 51825.*xi/16. - 155475.*(eta + 1.)*(xi + 1.)/8. + 39375.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 23625.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 1088325.*(xi + 1.)*((eta + 1.)*(eta + 1.))/32. - 362775.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 652995.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/128. + 122325./16. - 74025.*(eta + 1.)*(eta + 1.)/8. - 275625.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 13125.*(xi + 1.)*(xi + 1.)/4. + 91875.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 24675.*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 55125.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 7875.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 44415.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32.);
                      case 21:
                        return RealGradient(0., -8325.*eta/4. - 30825.*xi/16. + 92475.*(eta + 1.)*(xi + 1.)/8. - 15750.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 23625.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 647325.*(xi + 1.)*(eta + 1.)*(eta + 1.)/32. + 215775.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 388395.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. - 58575./16. + 58275.*((eta + 1.)*(eta + 1.))/16. + 55125.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 2625.*((xi + 1.)*(xi + 1.)) - 18375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 33075.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 19425.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 55125.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. - 99225.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 34965.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64.);
                      case 22:
                        return RealGradient(0., 3375.*eta/4. + 2025.*xi/16. - 6075.*(eta + 1.)*(xi + 1.)/8. + 42525.*(xi + 1.)*((eta + 1.)*(eta + 1.))/32. - 14175.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 25515.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/128. + 13275./16. - 23625.*(eta + 1.)*(eta + 1.)/16. + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 14175.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64.);
                      case 23:
                        return RealGradient(0., -675.*eta - 2025.*xi/16. + 6075.*(eta + 1.)*(xi + 1.)/8. - 42525.*(xi + 1.)*(eta + 1.)*(eta + 1.)/32. + 14175.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 25515.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. - 11025./16. + 4725.*((eta + 1.)*(eta + 1.))/4. - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16.);
                      case 24:
                        return RealGradient(-20475.*eta/8. + 20475.*(eta + 1.)*(xi + 1.)/4. - 36855.*(eta + 1.)*(xi + 1.)*(xi + 1.)/16. - 74025.*(xi + 1.)*(eta + 1.)*(eta + 1.)/4. + 362775.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 91875.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 33075.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 20475./8. + 74025.*((eta + 1.)*(eta + 1.))/8. + 133245.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/16. - 652995.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 165375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 59535.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 362775.*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 91875.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 33075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 25:
                        return RealGradient(525.*eta - 1050.*(eta + 1.)*(xi + 1.) + 945.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. + 58275.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. - 215775.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 18375.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 33075.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 525. - 58275.*(eta + 1.)*(eta + 1.)/16. - 104895.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 388395.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 33075.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 59535.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. + 215775.*((eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 18375.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 33075.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32., 0.);
                      case 26:
                        return RealGradient(-4725.*eta/4. + 4725.*(eta + 1.)*(xi + 1.)/2. - 8505.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 23625.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. + 14175.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 4725./4. + 23625.*((eta + 1.)*(eta + 1.))/16. + 42525.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/32. - 25515.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 14175.*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 27:
                        return RealGradient(4725.*eta/8. - 4725.*(eta + 1.)*(xi + 1.)/4. + 8505.*(eta + 1.)*((xi + 1.)*(xi + 1.))/16. + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 14175.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 4725./8. - 4725.*(eta + 1.)*(eta + 1.)/4. - 8505.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 25515.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. + 14175.*((eta + 1.)*(eta + 1.)*(eta + 1.))/32., 0.);
                      case 28:
                        return RealGradient(-12285.*eta/8. + 4095.*(eta + 1.)*(xi + 1.) - 36855.*(eta + 1.)*(xi + 1.)*(xi + 1.)/16. - 14805.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 72555.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 18375.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 6615.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 12285./8. + 44415.*((eta + 1.)*(eta + 1.))/8. + 133245.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/16. - 652995.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 165375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 59535.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 217665.*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 55125.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 19845.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 29:
                        return RealGradient(315.*eta - 840.*(eta + 1.)*(xi + 1.) + 945.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. + 11655.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 43155.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 7350.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)) - 6615.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 315. - 34965.*(eta + 1.)*(eta + 1.)/16. - 104895.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 388395.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 33075.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 59535.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. + 129465.*((eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 11025.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 19845.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32., 0.);
                      case 30:
                        return RealGradient(-2835.*eta/4. + 1890.*(eta + 1.)*(xi + 1.) - 8505.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 4725.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 2835.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 2835./4. + 14175.*((eta + 1.)*(eta + 1.))/16. + 42525.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/32. - 25515.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. - 8505.*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 31:
                        return RealGradient(2835.*eta/8. - 945.*(eta + 1.)*(xi + 1.) + 8505.*(eta + 1.)*((xi + 1.)*(xi + 1.))/16. + 1890.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 2835.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 2835./8. - 2835.*(eta + 1.)*(eta + 1.)/4. - 8505.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 25515.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/64. + 8505.*((eta + 1.)*(eta + 1.)*(eta + 1.))/32., 0.);
                      case 32:
                        return RealGradient(0., 2115.*eta + 10365.*xi/16. - 31095.*(eta + 1.)*(xi + 1.)/4. + 7875.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 652995.*(xi + 1.)*((eta + 1.)*(eta + 1.))/32. - 72555.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 652995.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/128. + 41385./16. - 44415.*(eta + 1.)*(eta + 1.)/8. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 99225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 2625.*(xi + 1.)*(xi + 1.)/4. + 18375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. + 4935.*((eta + 1.)*(eta + 1.)*(eta + 1.)) - 11025.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 99225.*((xi + 1.)*(xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64. - 44415.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32.);
                      case 33:
                        return RealGradient(0., -1665.*eta/2. - 6165.*xi/16. + 18495.*(eta + 1.)*(xi + 1.)/4. - 6300.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 4725.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 388395.*(xi + 1.)*(eta + 1.)*(eta + 1.)/32. + 43155.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 388395.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. - 18375./16. + 34965.*((eta + 1.)*(eta + 1.))/16. + 33075.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 99225.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 525.*((xi + 1.)*(xi + 1.)) - 14700.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.) + 33075.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 3885.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 11025.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. - 99225.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 34965.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64.);
                      case 34:
                        return RealGradient(0., 675.*eta/2. + 405.*xi/16. - 1215.*(eta + 1.)*(xi + 1.)/4. + 25515.*(xi + 1.)*((eta + 1.)*(eta + 1.))/32. - 2835.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 25515.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/128. + 5355./16. - 14175.*(eta + 1.)*(eta + 1.)/16. + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 14175.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64.);
                      case 35:
                        return RealGradient(0., -270.*eta - 405.*xi/16. + 1215.*(eta + 1.)*(xi + 1.)/4. - 25515.*(xi + 1.)*(eta + 1.)*(eta + 1.)/32. + 2835.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 25515.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. - 4365./16. + 2835.*((eta + 1.)*(eta + 1.))/4. - 630.*(eta + 1.)*(eta + 1.)*(eta + 1.) + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/16.);
                      case 36:
                        return RealGradient(60.*eta + 60. - 495.*(eta + 1.)*(eta + 1.)/2. + 1275.*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 2625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32., 0.);
                      case 37:
                        return RealGradient(60.*eta + 60. - 495.*(eta + 1.)*(eta + 1.)/2. + 1275.*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 2625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32., 0.);
                      case 38:
                        return RealGradient(-30.*eta - 30. + 495.*((eta + 1.)*(eta + 1.))/4. - 1275.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 2625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/64., 0.);
                      case 39:
                        return RealGradient(0., 594.*eta + 2295.*xi/2. - 2295.*(eta + 1.)*(xi + 1.) + 4725.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 2835.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 3825.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. + 2889./2. - 495.*(eta + 1.)*(eta + 1.)/2. - 7875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 4725.*(xi + 1.)*(xi + 1.)/4. + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 40:
                        return RealGradient(0., 396.*eta + 765.*xi/2. - 1530.*(eta + 1.)*(xi + 1.) + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 945.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 3825.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. + 1359./2. - 495.*(eta + 1.)*(eta + 1.)/2. - 7875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 1575.*(xi + 1.)*(xi + 1.)/4. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 41:
                        return RealGradient(0., -495.*eta/2. - 765.*xi/4. + 3825.*(eta + 1.)*(xi + 1.)/4. - 7875.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 4725.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 3825.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 1557./4. + 495.*((eta + 1.)*(eta + 1.))/4. + 7875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/16. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 1575.*((xi + 1.)*(xi + 1.))/8. - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)/16.);
                      case 42:
                        return RealGradient(0., -216.*eta - 675.*xi + 1350.*(eta + 1.)*(xi + 1.) - 1890.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 2835.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 1125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 783. + 90.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 945.*((xi + 1.)*(xi + 1.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8.);
                      case 43:
                        return RealGradient(0., -144.*eta - 225.*xi + 900.*(eta + 1.)*(xi + 1.) - 1260.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 945.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 1125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 333. + 90.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 4725.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 315.*((xi + 1.)*(xi + 1.)) - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8.);
                      case 44:
                        return RealGradient(0., 90.*eta + 225.*xi/2. - 1125.*(eta + 1.)*(xi + 1.)/2. + 1575.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 1125.*(xi + 1.)*((eta + 1.)*(eta + 1.))/4. + 369./2. - 45.*(eta + 1.)*(eta + 1.) - 1575.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. + 4725.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/32. - 315.*(xi + 1.)*(xi + 1.)/2. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.))/16.);
                      case 45:
                        return RealGradient(-15.*eta/2. - 15./2. + 90.*((eta + 1.)*(eta + 1.)) - 375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 46:
                        return RealGradient(-15.*eta/2. - 15./2. + 90.*((eta + 1.)*(eta + 1.)) - 375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 525.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32., 0.);
                      case 47:
                        return RealGradient(15.*eta/4. + 15./4. - 45.*(eta + 1.)*(eta + 1.) + 375.*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 525.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.))/64., 0.);
                      case 48:
                        return RealGradient(0., 81.*eta + 135.*xi/4. - 135.*(eta + 1.)*(xi + 1.)/2. + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 297./4. - 135.*(eta + 1.)*(eta + 1.)/4.);
                      case 49:
                        return RealGradient(0., -54.*eta - 135.*xi/4. + 135.*(eta + 1.)*(xi + 1.)/2. - 225.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 243./4. + 45.*((eta + 1.)*(eta + 1.))/2.);
                      case 50:
                        return RealGradient(-30.*eta - 30. + 135.*((eta + 1.)*(eta + 1.))/4. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 51:
                        return RealGradient(15.*eta/2. + 15./2. - 45.*(eta + 1.)*(eta + 1.)/2. + 75.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8., 0.);
                      case 52:
                        return RealGradient(-30.*eta - 30. + 135.*((eta + 1.)*(eta + 1.))/4. - 75.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 53:
                        return RealGradient(15.*eta/2. + 15./2. - 45.*(eta + 1.)*(eta + 1.)/2. + 75.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8., 0.);
                      case 54:
                        return RealGradient(0., 54.*eta + 45.*xi/4. - 45.*(eta + 1.)*(xi + 1.) + 225.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 207./4. - 135.*(eta + 1.)*(eta + 1.)/4.);
                      case 55:
                        return RealGradient(0., -36.*eta - 45.*xi/4. + 45.*(eta + 1.)*(xi + 1.) - 225.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. - 153./4. + 45.*((eta + 1.)*(eta + 1.))/2.);
                      case 56:
                        return RealGradient(0., 0.);
                      case 57:
                        return RealGradient(0., 15.*xi/4. - 3./4.);
                      case 58:
                        return RealGradient(0., -15.*xi/4. - 3./4.);
                      case 59:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(-5625.*eta - 13125.*xi/4. + 78750.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 157500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 94500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 275625.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 367500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 165375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 31875./4. + 39375.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 551250.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 330750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 13125.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 735000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 330750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 52500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 441000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 198450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 23625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 1:
                        return sign * RealGradient(214875.*eta/128. + 527625.*xi/512. - 1582875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. + 3063375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 874125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 11080125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 3693375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 6648075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. + 1243875./512. - 1504125.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 21443625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 6118875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 1021125.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 7147875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. - 12866175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 501375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 291375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/128. + 3671325.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 902475.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128., 0.);
                      case 2:
                        return sign * RealGradient(-7875.*eta/8. - 28875.*xi/32. + 86625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/4. - 212625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 70875.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 606375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/8. + 202125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 363825.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. - 55125./32. + 55125.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 1488375.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 496125.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 70875.*((xi + 1.)*(xi + 1.)/(2.*2.))/16. - 496125.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 893025.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 18375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 165375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 23625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/8. - 297675.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 33075.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8., 0.);
                      case 3:
                        return sign * RealGradient(52875.*eta/128. + 233625.*xi/512. - 700875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. + 2181375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/64. - 874125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 4906125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 1635375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/32. + 2943675.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. + 409875./512. - 370125.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 15269625.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 6118875.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 727125.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 5089875.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. - 9161775.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 123375.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. - 2039625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 291375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/128. + 3671325.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 222075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128., 0.);
                      case 4:
                        return sign * RealGradient(-2250.*eta - 7875.*xi/4. + 47250.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 126000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 94500.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 165375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 220500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 15375./4. + 15750.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 441000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 330750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 10500.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 588000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 264600.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 21000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 441000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 198450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 5:
                        return sign * RealGradient(0., 2625.*eta/4. + 2250.*xi - 31500.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 165375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 294000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 63000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 37800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 10875./4. - 2625.*(eta + 1.)*(eta + 1.)/(2.*2.) - 330750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 588000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 330750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 23625.*(xi + 1.)*(xi + 1.)/(2.*2.) + 198450.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 352800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 198450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 42000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 23625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 6:
                        return sign * RealGradient(0., -105525.*eta/512. - 42975.*xi/64. + 316575.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. - 6648075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 738675.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. - 6648075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 612675.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/32. + 174825.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/16. - 420675./512. + 204225.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 12866175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 1429575.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/8. + 12866175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 902475.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 3671325.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. - 58275.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/128. + 407925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 100275.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/8. + 902475.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128.);
                      case 7:
                        return sign * RealGradient(0., 5775.*eta/32. + 1575.*xi/4. - 17325.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 363825.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 80850.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 363825.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. + 42525.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 17325./32. - 14175.*(eta + 1.)*(eta + 1.)/(2.*2.)/16. - 893025.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 198450.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 893025.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 33075.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 297675.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. + 4725.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/8. - 132300.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 297675.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 7350.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 33075.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8.);
                      case 8:
                        return sign * RealGradient(0., -46725.*eta/512. - 10575.*xi/64. + 140175.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. - 2943675.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 327075.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/8. - 2943675.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 436275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/32. + 174825.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/16. - 124275./512. + 145425.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 9161775.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 1017975.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/8. + 9161775.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 222075.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 3671325.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. - 58275.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/128. + 407925.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/4. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 24675.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/8. + 222075.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128.);
                      case 9:
                        return sign * RealGradient(0., 1575.*eta/4. + 900.*xi - 18900.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 99225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 176400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 99225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 50400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 37800.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4875./4. - 2100.*(eta + 1.)*(eta + 1.)/(2.*2.) - 264600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 470400.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 264600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 9450.*(xi + 1.)*(xi + 1.)/(2.*2.) + 198450.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 352800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 198450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 16800.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 9450.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 10:
                        return sign * RealGradient(-900.*eta - 1575.*xi/4. + 18900.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 50400.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 37800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 99225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 176400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 4875./4. + 9450.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 264600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 198450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2100.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 470400.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 264600.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 16800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 352800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 198450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 11:
                        return sign * RealGradient(10575.*eta/64. + 46725.*xi/512. - 140175.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. + 436275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 174825.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 2943675.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 327075.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 2943675.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. + 124275./512. - 222075.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 9161775.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 3671325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 145425.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 1017975.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/8. - 9161775.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 24675.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/8. - 407925.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 58275.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/128. + 3671325.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 222075.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128., 0.);
                      case 12:
                        return sign * RealGradient(-1575.*eta/4. - 5775.*xi/32. + 17325.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 42525.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 363825.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/8. + 80850.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 363825.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/8. - 17325./32. + 33075.*((eta + 1.)*(eta + 1.)/(2.*2.))/8. + 893025.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 297675.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 14175.*((xi + 1.)*(xi + 1.)/(2.*2.))/16. - 198450.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 893025.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8. - 7350.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 132300.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/8. - 297675.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4. + 33075.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/8., 0.);
                      case 13:
                        return sign * RealGradient(42975.*eta/64. + 105525.*xi/512. - 316575.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/32. + 612675.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/32. - 174825.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/16. + 6648075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/128. - 738675.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/8. + 6648075.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/128. + 420675./512. - 902475.*(eta + 1.)*(eta + 1.)/(2.*2.)/128. - 12866175.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 3671325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/64. - 204225.*(xi + 1.)*(xi + 1.)/(2.*2.)/256. + 1429575.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/8. - 12866175.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128. + 100275.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/8. - 407925.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/4. + 58275.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/128. + 3671325.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/64. - 902475.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/128., 0.);
                      case 14:
                        return sign * RealGradient(-2250.*eta - 2625.*xi/4. + 31500.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 63000.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 37800.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 165375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 294000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 165375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 10875./4. + 23625.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 330750.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 198450.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2625.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 588000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 330750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 42000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 352800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) - 198450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 23625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 15:
                        return sign * RealGradient(0., 7875.*eta/4. + 2250.*xi - 47250.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 165375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 220500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 99225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 126000.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 94500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 15375./4. - 10500.*(eta + 1.)*(eta + 1.)/(2.*2.) - 441000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 588000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 264600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 15750.*(xi + 1.)*(xi + 1.)/(2.*2.) + 330750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 441000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 198450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 21000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 9450.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 16:
                        return sign * RealGradient(0., -233625.*eta/512. - 52875.*xi/128. + 700875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. - 4906125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 1635375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 2943675.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 2181375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 874125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. - 409875./512. + 727125.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 15269625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 5089875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 9161775.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 370125.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 6118875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. - 291375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/128. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 123375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 222075.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128.);
                      case 17:
                        return sign * RealGradient(0., 28875.*eta/32. + 7875.*xi/8. - 86625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/4. + 606375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/8. - 202125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 363825.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/8. + 212625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. - 70875.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. + 55125./32. - 70875.*(eta + 1.)*(eta + 1.)/(2.*2.)/16. - 1488375.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 496125.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 893025.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8. - 55125.*(xi + 1.)*(xi + 1.)/(2.*2.)/8. + 496125.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/4. + 23625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/8. - 165375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 297675.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4. + 18375.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 33075.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/8.);
                      case 18:
                        return sign * RealGradient(0., -527625.*eta/512. - 214875.*xi/128. + 1582875.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/64. - 11080125.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/128. + 3693375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/32. - 6648075.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/128. - 3063375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/64. + 874125.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/32. - 1243875./512. + 1021125.*((eta + 1.)*(eta + 1.)/(2.*2.))/256. + 21443625.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 7147875.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 12866175.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128. + 1504125.*((xi + 1.)*(xi + 1.)/(2.*2.))/128. - 6118875.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/64. - 291375.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/128. + 2039625.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/16. - 3671325.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/64. - 501375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/32. + 902475.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/128.);
                      case 19:
                        return sign * RealGradient(0., 13125.*eta/4. + 5625.*xi - 78750.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 275625.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 367500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 157500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 94500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 31875./4. - 13125.*(eta + 1.)*(eta + 1.)/(2.*2.) - 551250.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 735000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 330750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 39375.*(xi + 1.)*(xi + 1.)/(2.*2.) + 330750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 441000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 198450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 52500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 23625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 20:
                        return RealGradient(0., 20475.*eta/8. + 10575.*xi/2. - 74025.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 1088325.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 367500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 165375.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 148050.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 88830.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 56925./8. - 20475.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 1088325.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 735000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 330750.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 155475.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 652995.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. + 12285.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 441000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 198450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 52500.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 23625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 21:
                        return RealGradient(0., -525.*eta - 8325.*xi/4. + 58275.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 647325.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 294000.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 165375.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 58275.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 34965.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 9825./4. + 2100.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 647325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 588000.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 330750.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 92475.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 388395.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 1260.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 352800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 198450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 42000.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 23625.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 22:
                        return RealGradient(0., 4725.*eta/4. + 3375.*xi/4. - 23625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 42525.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/4. + 23625.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3375./2. - 4725.*(eta + 1.)*(eta + 1.)/(2.*2.) - 42525.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. - 6075.*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 25515.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)));
                      case 23:
                        return RealGradient(0., -4725.*eta/8. - 675.*xi + 9450.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 42525.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. - 18900.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 11340.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 8775./8. + 4725.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 42525.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. + 6075.*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 25515.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2.);
                      case 24:
                        return RealGradient(-10575.*eta/2. - 20475.*xi/8. + 74025.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 148050.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 88830.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1088325.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/4. + 367500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 165375.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 56925./8. + 155475.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 1088325.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 652995.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 20475.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 735000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 330750.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 52500.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 441000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 12285.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. - 198450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 23625.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 25:
                        return RealGradient(8325.*eta/4. + 525.*xi - 58275.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 58275.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 34965.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 647325.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. - 294000.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 165375.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 9825./4. - 92475.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 647325.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 388395.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 2100.*(xi + 1.)*(xi + 1.)/(2.*2.) + 588000.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 330750.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 42000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 352800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1260.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 198450.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 23625.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 26:
                        return RealGradient(-3375.*eta/4. - 4725.*xi/4. + 23625.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 23625.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 42525.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 3375./2. + 6075.*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 42525.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 25515.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 4725.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 27:
                        return RealGradient(675.*eta + 4725.*xi/8. - 9450.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 18900.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 11340.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 42525.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 8775./8. - 6075.*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 42525.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 25515.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 4725.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2., 0.);
                      case 28:
                        return RealGradient(-2115.*eta - 12285.*xi/8. + 44415.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 118440.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 88830.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 652995.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/4. + 220500.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 99225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 26865./8. + 31095.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 435330.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 652995.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 8190.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 588000.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 264600.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 21000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 441000.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 12285.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. - 198450.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 9450.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 29:
                        return RealGradient(1665.*eta/2. + 315.*xi - 34965.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 46620.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 34965.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 388395.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. - 176400.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 99225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 2175./2. - 18495.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 258930.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 388395.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 1680.*(xi + 1.)*(xi + 1.)/(2.*2.) + 470400.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 264600.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) + 16800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 352800.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1260.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) + 198450.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) - 9450.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 30:
                        return RealGradient(-675.*eta/2. - 2835.*xi/4. + 14175.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 18900.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 14175.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 25515.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.)/4. - 3645./4. + 1215.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. + 17010.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 25515.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)/2. + 3780.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.), 0.);
                      case 31:
                        return RealGradient(270.*eta + 2835.*xi/8. - 5670.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 15120.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 11340.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 25515.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.))/4. + 4455./8. - 1215.*(eta + 1.)*(eta + 1.)/(2.*2.)/2. - 17010.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 25515.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2. - 1890.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.))/2., 0.);
                      case 32:
                        return RealGradient(0., 12285.*eta/8. + 2115.*xi - 44415.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 652995.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/4. - 220500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 99225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 118440.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 88830.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 26865./8. - 8190.*(eta + 1.)*(eta + 1.)/(2.*2.) - 435330.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 588000.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 264600.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 31095.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 652995.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. + 12285.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. - 441000.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 198450.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 21000.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 9450.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 33:
                        return RealGradient(0., -315.*eta - 1665.*xi/2. + 34965.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. - 388395.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. + 176400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 99225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 46620.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 34965.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2175./2. + 1680.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 258930.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) - 470400.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 264600.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 18495.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 388395.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 1260.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 352800.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 198450.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 16800.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 9450.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 34:
                        return RealGradient(0., 2835.*eta/4. + 675.*xi/2. - 14175.*(eta/2. + 1./2.)*(xi/2. + 1./2.)/2. + 25515.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.))/4. + 18900.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 14175.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 3645./4. - 3780.*(eta + 1.)*(eta + 1.)/(2.*2.) - 17010.*(eta + 1.)*(eta + 1.)/(2.*2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 1215.*(xi + 1.)*(xi + 1.)/(2.*2.)/2. + 25515.*((xi + 1.)*(xi + 1.)/(2.*2.))*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.))/2. + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)));
                      case 35:
                        return RealGradient(0., -2835.*eta/8. - 270.*xi + 5670.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 25515.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.)/4. - 15120.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 11340.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4455./8. + 1890.*((eta + 1.)*(eta + 1.)/(2.*2.)) + 17010.*((eta + 1.)*(eta + 1.)/(2.*2.))*((xi + 1.)*(xi + 1.)/(2.*2.)) + 1215.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 25515.*(xi + 1.)*(xi + 1.)/(2.*2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2. - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)/2.);
                      case 36:
                        return RealGradient(594.*eta + 60.*xi - 1980.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 7650.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 10500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 582. - 4590.*(eta + 1.)*(eta + 1.)/(2.*2.) + 6300.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 37:
                        return RealGradient(396.*eta + 60.*xi - 1980.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 7650.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 10500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)) + 408. - 3060.*(eta + 1.)*(eta + 1.)/(2.*2.) + 4200.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 1890.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.), 0.);
                      case 38:
                        return RealGradient(-495.*eta/2. - 30.*xi + 990.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3825.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 5250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/2. - 495./2. + 3825.*((eta + 1.)*(eta + 1.)/(2.*2.))/2. - 2625.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/4., 0.);
                      case 39:
                        return RealGradient(0., 60.*eta + 594.*xi - 1980.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 7650.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 582. - 4590.*(xi + 1.)*(xi + 1.)/(2.*2.) + 6300.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 40:
                        return RealGradient(0., 60.*eta + 396.*xi - 1980.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 7650.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 10500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)) + 408. - 3060.*(xi + 1.)*(xi + 1.)/(2.*2.) + 4200.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 1890.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.));
                      case 41:
                        return RealGradient(0., -30.*eta - 495.*xi/2. + 990.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 3825.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 5250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/2. - 495./2. + 3825.*((xi + 1.)*(xi + 1.)/(2.*2.))/2. - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/4.);
                      case 42:
                        return RealGradient(0., -15.*eta/2. - 216.*xi + 720.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 4500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 429./2. + 2700.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 5040.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 2835.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 43:
                        return RealGradient(0., -15.*eta/2. - 144.*xi + 720.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 4500.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) + 8400.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.) - 291./2. + 1800.*((xi + 1.)*(xi + 1.)/(2.*2.)) - 3360.*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 1890.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)));
                      case 44:
                        return RealGradient(0., 15.*eta/4. + 90.*xi - 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2250.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) - 4200.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.) + 4725.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.))/2. + 90. - 1125.*(xi + 1.)*(xi + 1.)/(2.*2.) + 2100.*((xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.)) - 4725.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/(2.*2.*2.*2.)/4.);
                      case 45:
                        return RealGradient(-216.*eta - 15.*xi/2. + 720.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 4500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 429./2. + 2700.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 5040.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 46:
                        return RealGradient(-144.*eta - 15.*xi/2. + 720.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 4500.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) + 8400.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.) - 291./2. + 1800.*((eta + 1.)*(eta + 1.)/(2.*2.)) - 3360.*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 1890.*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)), 0.);
                      case 47:
                        return RealGradient(90.*eta + 15.*xi/4. - 360.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 2250.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) - 4200.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.) + 4725.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.))/2. + 90. - 1125.*(eta + 1.)*(eta + 1.)/(2.*2.) + 2100.*((eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.)) - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/(2.*2.*2.*2.)/4., 0.);
                      case 48:
                        return RealGradient(0., 30.*eta + 81.*xi - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 75. - 135.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 49:
                        return RealGradient(0., -15.*eta/2. - 54.*xi + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 105./2. + 135.*((xi + 1.)*(xi + 1.)/(2.*2.)));
                      case 50:
                        return RealGradient(-81.*eta - 30.*xi + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 75. + 135.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 51:
                        return RealGradient(54.*eta + 15.*xi/2. - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 105./2. - 135.*(eta + 1.)*(eta + 1.)/(2.*2.), 0.);
                      case 52:
                        return RealGradient(-54.*eta - 30.*xi + 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(xi/2. + 1./2.)*(eta + 1.)*(eta + 1.)/(2.*2.) - 60. + 90.*((eta + 1.)*(eta + 1.)/(2.*2.)), 0.);
                      case 53:
                        return RealGradient(36.*eta + 15.*xi/2. - 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(xi/2. + 1./2.)*((eta + 1.)*(eta + 1.)/(2.*2.)) + 75./2. - 90.*(eta + 1.)*(eta + 1.)/(2.*2.), 0.);
                      case 54:
                        return RealGradient(0., 30.*eta + 54.*xi - 270.*(eta/2. + 1./2.)*(xi/2. + 1./2.) + 225.*(eta/2. + 1./2.)*((xi + 1.)*(xi + 1.)/(2.*2.)) + 60. - 90.*(xi + 1.)*(xi + 1.)/(2.*2.));
                      case 55:
                        return RealGradient(0., -15.*eta/2. - 36.*xi + 180.*(eta/2. + 1./2.)*(xi/2. + 1./2.) - 225.*(eta/2. + 1./2.)*(xi + 1.)*(xi + 1.)/(2.*2.) - 75./2. + 90.*((xi + 1.)*(xi + 1.)/(2.*2.)));
                      case 56:
                        return RealGradient(0., 0.);
                      case 57:
                        return RealGradient(0., 0.);
                      case 58:
                        return RealGradient(0., 0.);
                      case 59:
                        return RealGradient(0., 0.);
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
                        return sign * RealGradient(-13125.*eta/4. - 5625.*xi + 39375.*(eta + 1.)*(xi + 1.)/2. - 275625.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 91875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 165375.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 39375.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 23625.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 31875./4. + 13125.*((eta + 1.)*(eta + 1.))/4. + 275625.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 91875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 39375.*((xi + 1.)*(xi + 1.))/4. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 7875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 55125.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. - 13125.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 23625.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/16., 0.);
                      case 1:
                        return sign * RealGradient(254625.*eta/512. + 214875.*xi/128. - 1504125.*(eta + 1.)*(xi + 1.)/256. + 11080125.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. - 7147875.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 6118875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. + 1504125.*(xi + 1.)*((eta + 1.)*(eta + 1.))/256. - 902475.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 1041375./512. - 254625.*(eta + 1.)*(eta + 1.)/512. - 11080125.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/1024. + 7147875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 6118875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4096. - 1582875.*(xi + 1.)*(xi + 1.)/512. + 6648075.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. + 152775.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 4288725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8192. + 1021125.*((xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 874125.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048., 0.);
                      case 2:
                        return sign * RealGradient(-7875.*eta/32. - 7875.*xi/8. + 55125.*(eta + 1.)*(xi + 1.)/16. - 606375.*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. + 496125.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 496125.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/256. - 55125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/16. + 33075.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 37125./32. + 7875.*((eta + 1.)*(eta + 1.))/32. + 606375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/64. - 496125.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 496125.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/256. + 86625.*((xi + 1.)*(xi + 1.))/32. - 363825.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 297675.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/512. - 70875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 70875.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/128., 0.);
                      case 3:
                        return sign * RealGradient(44625.*eta/512. + 52875.*xi/128. - 370125.*(eta + 1.)*(xi + 1.)/256. + 4906125.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. - 5089875.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 6118875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. + 370125.*(xi + 1.)*((eta + 1.)*(eta + 1.))/256. - 222075.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 243375./512. - 44625.*(eta + 1.)*(eta + 1.)/512. - 4906125.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/1024. + 5089875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 6118875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4096. - 700875.*(xi + 1.)*(xi + 1.)/512. + 2943675.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. + 26775.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 3053925.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8192. + 727125.*((xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 874125.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048., 0.);
                      case 4:
                        return sign * RealGradient(-2625.*eta/4. - 2250.*xi + 7875.*(eta + 1.)*(xi + 1.) - 165375.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 18375.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 165375.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 7875.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 10875./4. + 2625.*((eta + 1.)*(eta + 1.))/4. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 18375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 23625.*((xi + 1.)*(xi + 1.))/4. - 99225.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 11025.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. - 5250.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 23625.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/16., 0.);
                      case 5:
                        return sign * RealGradient(0., 2625.*xi/4. - 2625.*(eta + 1.)*(xi + 1.)/2. + 7875.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 55125.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 18375.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 33075.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 2625./4. - 14175.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. + 99225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 33075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 59535.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 7875.*(xi + 1.)*(xi + 1.)/2. + 55125.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 18375.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 33075.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 6:
                        return sign * RealGradient(0., -105525.*xi/512. + 204225.*(eta + 1.)*(xi + 1.)/512. - 612675.*(eta + 1.)*(xi + 1.)*(xi + 1.)/256. + 4288725.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 1429575.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/512. + 2573235.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. - 174825.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. - 105525./512. + 524475.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/512. - 3671325.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 1223775.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 2202795.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8192. + 316575.*((xi + 1.)*(xi + 1.))/256. - 2216025.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 738675.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 1329615.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4096.);
                      case 7:
                        return sign * RealGradient(0., 5775.*xi/32. - 14175.*(eta + 1.)*(xi + 1.)/32. + 42525.*(eta + 1.)*((xi + 1.)*(xi + 1.))/16. - 297675.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 99225.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. - 178605.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/256. + 14175.*(xi + 1.)*((eta + 1.)*(eta + 1.))/64. + 5775./32. - 42525.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 297675.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 99225.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 178605.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 17325.*(xi + 1.)*(xi + 1.)/16. + 121275.*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 40425.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 72765.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/256.);
                      case 8:
                        return sign * RealGradient(0., -46725.*xi/512. + 145425.*(eta + 1.)*(xi + 1.)/512. - 436275.*(eta + 1.)*(xi + 1.)*(xi + 1.)/256. + 3053925.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 1017975.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/512. + 1832355.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. - 174825.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. - 46725./512. + 524475.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/512. - 3671325.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 1223775.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 2202795.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8192. + 140175.*((xi + 1.)*(xi + 1.))/256. - 981225.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 327075.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 588735.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4096.);
                      case 9:
                        return sign * RealGradient(0., 1575.*xi/4. - 1050.*(eta + 1.)*(xi + 1.) + 6300.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 11025.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 7350.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)) - 6615.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 1575./4. - 14175.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. + 99225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 33075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 59535.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 4725.*(xi + 1.)*(xi + 1.)/2. + 33075.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 11025.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 19845.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 10:
                        return sign * RealGradient(-1575.*eta/4. - 900.*xi + 4725.*(eta + 1.)*(xi + 1.) - 99225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 11025.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 99225.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 6300.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 4875./4. + 525.*((eta + 1.)*(eta + 1.)) + 33075.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 14700.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 33075.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 4725.*((xi + 1.)*(xi + 1.))/2. - 99225.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 11025.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. - 2100.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8., 0.);
                      case 11:
                        return sign * RealGradient(26775.*eta/512. + 10575.*xi/64. - 222075.*(eta + 1.)*(xi + 1.)/256. + 2943675.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. - 3053925.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 3671325.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. + 74025.*(xi + 1.)*((eta + 1.)*(eta + 1.))/64. - 222075.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 106275./512. - 8925.*(eta + 1.)*(eta + 1.)/128. - 981225.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/256. + 1017975.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/256. - 1223775.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. - 140175.*(xi + 1.)*(xi + 1.)/256. + 2943675.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. + 26775.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 3053925.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8192. + 145425.*((xi + 1.)*(xi + 1.)*(xi + 1.))/256. - 174825.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024., 0.);
                      case 12:
                        return sign * RealGradient(-4725.*eta/32. - 1575.*xi/4. + 33075.*(eta + 1.)*(xi + 1.)/16. - 363825.*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. + 297675.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 297675.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/256. - 11025.*(xi + 1.)*(eta + 1.)*(eta + 1.)/4. + 33075.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/32. - 16425./32. + 1575.*((eta + 1.)*(eta + 1.))/8. + 121275.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/16. - 99225.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 99225.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. + 17325.*((xi + 1.)*(xi + 1.))/16. - 363825.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/128. - 4725.*(eta + 1.)*(eta + 1.)*(eta + 1.)/64. + 297675.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 297675.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/512. - 14175.*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 14175.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64., 0.);
                      case 13:
                        return sign * RealGradient(152775.*eta/512. + 42975.*xi/64. - 902475.*(eta + 1.)*(xi + 1.)/256. + 6648075.*(eta + 1.)*((xi + 1.)*(xi + 1.))/1024. - 4288725.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 3671325.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. + 300825.*(xi + 1.)*((eta + 1.)*(eta + 1.))/64. - 902475.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/512. + 467475./512. - 50925.*(eta + 1.)*(eta + 1.)/128. - 2216025.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/256. + 1429575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/256. - 1223775.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. - 316575.*(xi + 1.)*(xi + 1.)/256. + 6648075.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/2048. + 152775.*((eta + 1.)*(eta + 1.)*(eta + 1.))/1024. - 4288725.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 3671325.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8192. + 204225.*((xi + 1.)*(xi + 1.)*(xi + 1.))/256. - 174825.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024., 0.);
                      case 14:
                        return sign * RealGradient(-7875.*eta/4. - 2250.*xi + 23625.*(eta + 1.)*(xi + 1.)/2. - 165375.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 55125.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 99225.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. - 15750.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 23625.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 15375./4. + 2625.*((eta + 1.)*(eta + 1.)) + 55125.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. - 18375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 33075.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8. + 7875.*((xi + 1.)*(xi + 1.))/2. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 7875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 55125.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 4725.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8., 0.);
                      case 15:
                        return sign * RealGradient(0., 420.*eta + 7875.*xi/4. - 5250.*(eta + 1.)*(xi + 1.) + 15750.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 18375.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 18375.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 6615.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 23625.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 8925./4. - 945.*(eta + 1.)*(eta + 1.)/4. - 70875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 59535.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 23625.*(xi + 1.)*(xi + 1.)/4. + 55125.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 55125.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 19845.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 16:
                        return sign * RealGradient(0., -29085.*eta/256. - 233625.*xi/512. + 727125.*(eta + 1.)*(xi + 1.)/512. - 2181375.*(eta + 1.)*(xi + 1.)*(xi + 1.)/512. + 5089875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 5089875.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 1832355.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. - 874125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. - 273105./512. + 34965.*((eta + 1.)*(eta + 1.))/512. + 2622375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/1024. - 6118875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 6118875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. - 2202795.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8192. + 700875.*((xi + 1.)*(xi + 1.))/512. - 1635375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 1635375.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/2048. - 588735.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4096.);
                      case 17:
                        return sign * RealGradient(0., 2835.*eta/16. + 28875.*xi/32. - 70875.*(eta + 1.)*(xi + 1.)/32. + 212625.*(eta + 1.)*((xi + 1.)*(xi + 1.))/32. - 496125.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 496125.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 178605.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/256. + 70875.*(xi + 1.)*((eta + 1.)*(eta + 1.))/64. + 32235./32. - 2835.*(eta + 1.)*(eta + 1.)/32. - 212625.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/64. + 496125.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/128. - 496125.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/256. + 178605.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/512. - 86625.*(xi + 1.)*(xi + 1.)/32. + 202125.*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 202125.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. + 72765.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/256.);
                      case 18:
                        return sign * RealGradient(0., -40845.*eta/256. - 527625.*xi/512. + 1021125.*(eta + 1.)*(xi + 1.)/512. - 3063375.*(eta + 1.)*(xi + 1.)*(xi + 1.)/512. + 7147875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/1024. - 7147875.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 2573235.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. - 874125.*(xi + 1.)*(eta + 1.)*(eta + 1.)/1024. - 567105./512. + 34965.*((eta + 1.)*(eta + 1.))/512. + 2622375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/1024. - 6118875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2048. + 6118875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4096. - 2202795.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8192. + 1582875.*((xi + 1.)*(xi + 1.))/512. - 3693375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/1024. + 3693375.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/2048. - 1329615.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4096.);
                      case 19:
                        return sign * RealGradient(0., 525.*eta + 13125.*xi/4. - 13125.*(eta + 1.)*(xi + 1.)/2. + 39375.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 91875.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 91875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 33075.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 23625.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 14175./4. - 945.*(eta + 1.)*(eta + 1.)/4. - 70875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 59535.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 39375.*(xi + 1.)*(xi + 1.)/4. + 91875.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 91875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 33075.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 20:
                        return RealGradient(0., 20475.*xi/8. - 20475.*(eta + 1.)*(xi + 1.)/4. + 74025.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. - 362775.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 91875.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 33075.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 36855.*(xi + 1.)*((eta + 1.)*(eta + 1.))/16. + 20475./8. - 133245.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/16. + 652995.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 59535.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 74025.*(xi + 1.)*(xi + 1.)/8. + 362775.*((xi + 1.)*(xi + 1.)*(xi + 1.))/32. - 91875.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 33075.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 21:
                        return RealGradient(0., -525.*xi + 1050.*(eta + 1.)*(xi + 1.) - 58275.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. + 215775.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 18375.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 33075.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 945.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 525. + 104895.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/32. - 388395.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 33075.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 59535.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 58275.*((xi + 1.)*(xi + 1.))/16. - 215775.*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 18375.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 33075.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32.);
                      case 22:
                        return RealGradient(0., 4725.*xi/4. - 4725.*(eta + 1.)*(xi + 1.)/2. + 23625.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. - 14175.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 8505.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 4725./4. - 42525.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 25515.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 23625.*(xi + 1.)*(xi + 1.)/16. + 14175.*((xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 23:
                        return RealGradient(0., -4725.*xi/8. + 4725.*(eta + 1.)*(xi + 1.)/4. - 4725.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 14175.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 8505.*(xi + 1.)*(eta + 1.)*(eta + 1.)/16. - 4725./8. + 8505.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 25515.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 4725.*((xi + 1.)*(xi + 1.))/4. - 14175.*(xi + 1.)*(xi + 1.)*(xi + 1.)/32.);
                      case 24:
                        return RealGradient(-51825.*eta/16. - 10575.*xi/2. + 155475.*(eta + 1.)*(xi + 1.)/8. - 1088325.*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 362775.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 652995.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. - 39375.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. + 23625.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 122325./16. + 13125.*((eta + 1.)*(eta + 1.))/4. + 275625.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 91875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 74025.*((xi + 1.)*(xi + 1.))/8. - 165375.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 7875.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 55125.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. - 24675.*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 44415.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32., 0.);
                      case 25:
                        return RealGradient(30825.*eta/16. + 8325.*xi/4. - 92475.*(eta + 1.)*(xi + 1.)/8. + 647325.*(eta + 1.)*((xi + 1.)*(xi + 1.))/32. - 215775.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 388395.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/128. + 15750.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 23625.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 58575./16. - 2625.*(eta + 1.)*(eta + 1.) - 55125.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 18375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 33075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. - 58275.*(xi + 1.)*(xi + 1.)/16. + 165375.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 7875.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 55125.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. + 19425.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 34965.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64., 0.);
                      case 26:
                        return RealGradient(-2025.*eta/16. - 3375.*xi/4. + 6075.*(eta + 1.)*(xi + 1.)/8. - 42525.*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 14175.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/16. - 25515.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. - 13275./16. + 23625.*((xi + 1.)*(xi + 1.))/16. - 7875.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 14175.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64., 0.);
                      case 27:
                        return RealGradient(2025.*eta/16. + 675.*xi - 6075.*(eta + 1.)*(xi + 1.)/8. + 42525.*(eta + 1.)*((xi + 1.)*(xi + 1.))/32. - 14175.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 25515.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/128. + 11025./16. - 4725.*(xi + 1.)*(xi + 1.)/4. + 1575.*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16., 0.);
                      case 28:
                        return RealGradient(-10365.*eta/16. - 2115.*xi + 31095.*(eta + 1.)*(xi + 1.)/4. - 652995.*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 72555.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 652995.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. - 7875.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 41385./16. + 2625.*((eta + 1.)*(eta + 1.))/4. + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 18375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 165375.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. + 44415.*((xi + 1.)*(xi + 1.))/8. - 99225.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 1575.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8. + 11025.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 99225.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. - 4935.*(xi + 1.)*(xi + 1.)*(xi + 1.) + 44415.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32., 0.);
                      case 29:
                        return RealGradient(6165.*eta/16. + 1665.*xi/2. - 18495.*(eta + 1.)*(xi + 1.)/4. + 388395.*(eta + 1.)*((xi + 1.)*(xi + 1.))/32. - 43155.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 388395.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/128. + 6300.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 4725.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 18375./16. - 525.*(eta + 1.)*(eta + 1.) - 33075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 14700.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 33075.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. - 34965.*(xi + 1.)*(xi + 1.)/16. + 99225.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 1575.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8. - 11025.*(eta + 1.)*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 99225.*((eta + 1.)*(eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. + 3885.*((xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 34965.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64., 0.);
                      case 30:
                        return RealGradient(-405.*eta/16. - 675.*xi/2. + 1215.*(eta + 1.)*(xi + 1.)/4. - 25515.*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 2835.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 25515.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/128. - 5355./16. + 14175.*((xi + 1.)*(xi + 1.))/16. - 1575.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 14175.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64., 0.);
                      case 31:
                        return RealGradient(405.*eta/16. + 270.*xi - 1215.*(eta + 1.)*(xi + 1.)/4. + 25515.*(eta + 1.)*((xi + 1.)*(xi + 1.))/32. - 2835.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 25515.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/128. + 4365./16. - 2835.*(xi + 1.)*(xi + 1.)/4. + 630.*((xi + 1.)*(xi + 1.)*(xi + 1.)) - 2835.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16., 0.);
                      case 32:
                        return RealGradient(0., 12285.*xi/8. - 4095.*(eta + 1.)*(xi + 1.) + 14805.*(eta + 1.)*((xi + 1.)*(xi + 1.)) - 72555.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 18375.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/2. - 6615.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 36855.*(xi + 1.)*((eta + 1.)*(eta + 1.))/16. + 12285./8. - 133245.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/16. + 652995.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 165375.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 59535.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 44415.*(xi + 1.)*(xi + 1.)/8. + 217665.*((xi + 1.)*(xi + 1.)*(xi + 1.))/32. - 55125.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 19845.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 33:
                        return RealGradient(0., -315.*xi + 840.*(eta + 1.)*(xi + 1.) - 11655.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. + 43155.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 7350.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.) + 6615.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 945.*(xi + 1.)*(eta + 1.)*(eta + 1.)/2. - 315. + 104895.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/32. - 388395.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 33075.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/8. - 59535.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 34965.*((xi + 1.)*(xi + 1.))/16. - 129465.*(xi + 1.)*(xi + 1.)*(xi + 1.)/32. + 11025.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 19845.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32.);
                      case 34:
                        return RealGradient(0., 2835.*xi/4. - 1890.*(eta + 1.)*(xi + 1.) + 4725.*(eta + 1.)*((xi + 1.)*(xi + 1.))/2. - 2835.*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/4. + 8505.*(xi + 1.)*((eta + 1.)*(eta + 1.))/8. + 2835./4. - 42525.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/32. + 25515.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.)*(xi + 1.))/64. - 14175.*(xi + 1.)*(xi + 1.)/16. + 8505.*((xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 35:
                        return RealGradient(0., -2835.*xi/8. + 945.*(eta + 1.)*(xi + 1.) - 1890.*(eta + 1.)*(xi + 1.)*(xi + 1.) + 2835.*(eta + 1.)*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 8505.*(xi + 1.)*(eta + 1.)*(eta + 1.)/16. - 2835./8. + 8505.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/8. - 25515.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64. + 2835.*((xi + 1.)*(xi + 1.))/4. - 8505.*(xi + 1.)*(xi + 1.)*(xi + 1.)/32.);
                      case 36:
                        return RealGradient(2295.*eta/2. + 594.*xi - 2295.*(eta + 1.)*(xi + 1.) + 3825.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 2835.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/4. + 2889./2. - 4725.*(eta + 1.)*(eta + 1.)/4. - 7875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 495.*(xi + 1.)*(xi + 1.)/2. + 4725.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 2835.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8., 0.);
                      case 37:
                        return RealGradient(765.*eta/2. + 396.*xi - 1530.*(eta + 1.)*(xi + 1.) + 3825.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.)) - 945.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/2. + 1359./2. - 1575.*(eta + 1.)*(eta + 1.)/4. - 7875.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 495.*(xi + 1.)*(xi + 1.)/2. + 4725.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.))/8., 0.);
                      case 38:
                        return RealGradient(-765.*eta/4. - 495.*xi/2. + 3825.*(eta + 1.)*(xi + 1.)/4. - 3825.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 7875.*(xi + 1.)*(eta + 1.)*(eta + 1.)/8. + 4725.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/16. - 1557./4. + 1575.*((eta + 1.)*(eta + 1.))/8. + 7875.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/16. + 495.*((xi + 1.)*(xi + 1.))/4. - 4725.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/32. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)/16., 0.);
                      case 39:
                        return RealGradient(0., 60.*xi + 60. - 495.*(xi + 1.)*(xi + 1.)/2. + 1275.*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 40:
                        return RealGradient(0., 60.*xi + 60. - 495.*(xi + 1.)*(xi + 1.)/2. + 1275.*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 2625.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/16. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32.);
                      case 41:
                        return RealGradient(0., -30.*xi - 30. + 495.*((xi + 1.)*(xi + 1.))/4. - 1275.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 2625.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/32. - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/64.);
                      case 42:
                        return RealGradient(0., -15.*xi/2. - 15./2. + 90.*((xi + 1.)*(xi + 1.)) - 375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32.);
                      case 43:
                        return RealGradient(0., -15.*xi/2. - 15./2. + 90.*((xi + 1.)*(xi + 1.)) - 375.*(xi + 1.)*(xi + 1.)*(xi + 1.)/2. + 525.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 945.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/32.);
                      case 44:
                        return RealGradient(0., 15.*xi/4. + 15./4. - 45.*(xi + 1.)*(xi + 1.) + 375.*((xi + 1.)*(xi + 1.)*(xi + 1.))/4. - 525.*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)/8. + 945.*((xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.)*(xi + 1.))/64.);
                      case 45:
                        return RealGradient(-675.*eta - 216.*xi + 1350.*(eta + 1.)*(xi + 1.) - 1125.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 1890.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 2835.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/4. - 783. + 945.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. + 90.*((xi + 1.)*(xi + 1.)) - 4725.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 2835.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 46:
                        return RealGradient(-225.*eta - 144.*xi + 900.*(eta + 1.)*(xi + 1.) - 1125.*(eta + 1.)*(xi + 1.)*(xi + 1.)/2. - 1260.*(xi + 1.)*(eta + 1.)*(eta + 1.) + 945.*(xi + 1.)*((eta + 1.)*(eta + 1.)*(eta + 1.))/2. - 333. + 315.*((eta + 1.)*(eta + 1.)) + 1575.*((eta + 1.)*(eta + 1.))*((xi + 1.)*(xi + 1.))/2. + 90.*((xi + 1.)*(xi + 1.)) - 4725.*(xi + 1.)*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. - 945.*(eta + 1.)*(eta + 1.)*(eta + 1.)/8., 0.);
                      case 47:
                        return RealGradient(225.*eta/2. + 90.*xi - 1125.*(eta + 1.)*(xi + 1.)/2. + 1125.*(eta + 1.)*((xi + 1.)*(xi + 1.))/4. + 1575.*(xi + 1.)*((eta + 1.)*(eta + 1.))/2. - 4725.*(xi + 1.)*(eta + 1.)*(eta + 1.)*(eta + 1.)/16. + 369./2. - 315.*(eta + 1.)*(eta + 1.)/2. - 1575.*(eta + 1.)*(eta + 1.)*(xi + 1.)*(xi + 1.)/4. - 45.*(xi + 1.)*(xi + 1.) + 4725.*((xi + 1.)*(xi + 1.))*((eta + 1.)*(eta + 1.)*(eta + 1.))/32. + 945.*((eta + 1.)*(eta + 1.)*(eta + 1.))/16., 0.);
                      case 48:
                        return RealGradient(0., 30.*xi + 30. - 135.*(xi + 1.)*(xi + 1.)/4. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 49:
                        return RealGradient(0., -15.*xi/2. - 15./2. + 45.*((xi + 1.)*(xi + 1.))/2. - 75.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8.);
                      case 50:
                        return RealGradient(-135.*eta/4. - 81.*xi + 135.*(eta + 1.)*(xi + 1.)/2. - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 297./4. + 135.*((xi + 1.)*(xi + 1.))/4., 0.);
                      case 51:
                        return RealGradient(135.*eta/4. + 54.*xi - 135.*(eta + 1.)*(xi + 1.)/2. + 225.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. + 243./4. - 45.*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 52:
                        return RealGradient(-45.*eta/4. - 54.*xi + 45.*(eta + 1.)*(xi + 1.) - 225.*(eta + 1.)*(xi + 1.)*(xi + 1.)/8. - 207./4. + 135.*((xi + 1.)*(xi + 1.))/4., 0.);
                      case 53:
                        return RealGradient(45.*eta/4. + 36.*xi - 45.*(eta + 1.)*(xi + 1.) + 225.*(eta + 1.)*((xi + 1.)*(xi + 1.))/8. + 153./4. - 45.*(xi + 1.)*(xi + 1.)/2., 0.);
                      case 54:
                        return RealGradient(0., 30.*xi + 30. - 135.*(xi + 1.)*(xi + 1.)/4. + 75.*((xi + 1.)*(xi + 1.)*(xi + 1.))/8.);
                      case 55:
                        return RealGradient(0., -15.*xi/2. - 15./2. + 45.*((xi + 1.)*(xi + 1.))/2. - 75.*(xi + 1.)*(xi + 1.)*(xi + 1.)/8.);
                      case 56:
                        return RealGradient(15.*eta/4. - 3./4., 0.);
                      case 57:
                        return RealGradient(0., 0.);
                      case 58:
                        return RealGradient(0., 0.);
                      case 59:
                        return RealGradient(-15.*eta/4. - 3./4., 0.);
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
                        return sign * RealGradient(30240.*eta*xi - 12600.*eta - 12600.*eta*xi*xi - 8400.*xi - 25200.*xi*eta*eta + 2100. + 22680.*(eta*eta) + 7560.*(xi*xi) - 12600.*eta*eta*eta, -45360.*eta*xi + 8400.*eta + 50400.*eta*(xi*xi) + 12600.*xi + 37800.*xi*(eta*eta) - 1400. - 15120.*eta*eta - 30240.*xi*xi + 8400.*(eta*eta*eta) + 21000.*(xi*xi*xi));
                      case 1:
                        return sign * RealGradient(-20475.*eta*xi/2. + 136185.*eta/32. + 58275.*eta*(xi*xi)/16. + 40845.*xi/16. + 146475.*xi*(eta*eta)/16. - 21105./32. - 240975.*eta*eta/32. - 34965.*xi*xi/16. + 127575.*(eta*eta*eta)/32., 252315.*eta*xi/16. - 39375.*eta/16. - 146475.*eta*xi*xi/8. - 125265.*xi/32. - 382725.*xi*eta*eta/32. + 1715./4. + 25515.*(eta*eta)/8. + 145845.*(xi*xi)/16. - 14175.*eta*eta*eta/16. - 97125.*xi*xi*xi/16.);
                      case 2:
                        return sign * RealGradient(7560.*eta*xi - 3255.*eta/2. - 4725.*eta*xi*xi - 2835.*xi - 4725.*xi*eta*eta + 1155./2. + 945.*(eta*eta)/2. + 2835.*(xi*xi) + 1575.*(eta*eta*eta)/2., -4725.*eta*xi + 105.*eta + 9450.*eta*(xi*xi) + 7455.*xi/2. - 4725.*xi*eta*eta/2. - 315. + 1890.*(eta*eta) - 10395.*xi*xi - 1575.*eta*eta*eta + 7875.*(xi*xi*xi));
                      case 3:
                        return sign * RealGradient(-2835.*eta*xi/2. - 16695.*eta/32. + 58275.*eta*(xi*xi)/16. + 29085.*xi/16. - 29925.*xi*eta*eta/16. - 9345./32. + 76545.*(eta*eta)/32. - 34965.*xi*xi/16. - 48825.*eta*eta*eta/32., -65205.*eta*xi/16. + 9345.*eta/16. + 29925.*eta*(xi*xi)/8. - 66465.*xi/32. + 146475.*xi*(eta*eta)/32. + 595./4. - 4725.*eta*eta/8. + 110565.*(xi*xi)/16. - 5775.*eta*eta*eta/16. - 97125.*xi*xi*xi/16.);
                      case 4:
                        return sign * RealGradient(10080.*eta*xi - 1680.*eta - 12600.*eta*xi*xi - 6720.*xi + 1260. + 7560.*(xi*xi), 8400.*xi - 700. - 25200.*xi*xi + 21000.*(xi*xi*xi));
                      case 5:
                        return sign * RealGradient(840.*eta*(12.*xi - 2. - 15.*xi*xi), 8400.*xi - 700. - 25200.*xi*xi + 21000.*(xi*xi*xi));
                      case 6:
                        return sign * RealGradient(2835.*eta*(-120.*eta*xi + 54.*eta + 60.*xi - 16. - 30.*eta*eta - 45.*xi*xi)/32., -127575.*eta*xi/8. + 2205.*eta + 42525.*eta*(xi*xi)/2. + 8505.*xi/2. + 127575.*xi*(eta*eta)/16. - 1785./4. - 31185.*eta*eta/16. - 161595.*xi*xi/16. + 1575.*(eta*eta*eta)/4. + 212625.*(xi*xi*xi)/32.);
                      case 7:
                        return sign * RealGradient(105.*eta*(-60.*eta*xi + 45.*eta + 18.*xi - 8. - 45.*eta*eta - 15.*xi*xi/2.), -13230.*eta*xi + 2520.*eta + 12600.*eta*(xi*xi) + 1680.*xi + 14175.*xi*(eta*eta) - 245. - 4725.*eta*eta - 2835.*xi*xi + 2100.*(eta*eta*eta) + 2625.*(xi*xi*xi)/2.);
                      case 8:
                        return sign * RealGradient(105.*eta*(-360.*eta*xi + 594.*eta + 84.*xi - 80. - 810.*eta*eta - 15.*xi*xi)/32., -36855.*eta*xi/8. + 1575.*eta + 4725.*eta*(xi*xi)/2. + 735.*xi/2. + 127575.*xi*(eta*eta)/16. - 385./4. - 76545.*eta*eta/16. - 5355.*xi*xi/16. + 14175.*(eta*eta*eta)/4. + 2625.*(xi*xi*xi)/32.);
                      case 9:
                        return sign * RealGradient(0., 0.);
                      case 10:
                        return sign * RealGradient(0., 0.);
                      case 11:
                        return sign * RealGradient(105.*eta*(330.*eta*xi + 180.*eta - 54.*xi - 11. - 465.*eta*eta - 15.*xi*xi)/32., 945.*eta*xi/8. + 1155.*eta/2. - 17325.*eta*xi*xi/8. - 1785.*xi/32. + 146475.*xi*(eta*eta)/32. - 595./32. - 76545.*eta*eta/32. + 2835.*(xi*xi)/32. + 9975.*(eta*eta*eta)/16. + 2625.*(xi*xi*xi)/32.);
                      case 12:
                        return sign * RealGradient(105.*eta*(45.*eta*xi - 18.*eta - 3.*xi + 5./2. + 15.*(eta*eta)/2. - 15.*xi*xi/2.), 9450.*eta*xi - 1680.*eta - 9450.*eta*xi*xi - 105.*xi/2. - 4725.*xi*eta*eta/2. + 175./2. - 945.*eta*eta/2. - 2205.*xi*xi/2. + 1575.*(eta*eta*eta) + 2625.*(xi*xi*xi)/2.);
                      case 13:
                        return sign * RealGradient(2835.*eta*(30.*eta*xi - 36.*eta + 30.*xi - 1. + 45.*(eta*eta) - 45.*xi*xi)/32., 76545.*eta*xi/8. - 7245.*eta/2. - 42525.*eta*xi*xi/8. + 127575.*xi/32. - 382725.*xi*eta*eta/32. - 11235./32. + 240975.*(eta*eta)/32. - 314685.*xi*xi/32. - 48825.*eta*eta*eta/16. + 212625.*(xi*xi*xi)/32.);
                      case 14:
                        return sign * RealGradient(840.*eta*(-30.*eta*xi + 18.*eta + 18.*xi - 5. - 15.*eta*eta - 15.*xi*xi), -60480.*eta*xi + 16800.*eta + 50400.*eta*(xi*xi) + 21000.*xi + 37800.*xi*(eta*eta) - 3500. - 22680.*eta*eta - 37800.*xi*xi + 8400.*(eta*eta*eta) + 21000.*(xi*xi*xi));
                      case 15:
                        return RealGradient(10080.*eta*(20.*eta*xi - 21.*eta - 14.*xi + 7. + 15.*(eta*eta) + 5.*(xi*xi)), 423360.*eta*xi - 94080.*eta - 403200.*eta*xi*xi - 70560.*xi - 453600.*xi*eta*eta + 9800. + 211680.*(eta*eta) + 141120.*(xi*xi) - 134400.*eta*eta*eta - 84000.*xi*xi*xi);
                      case 16:
                        return RealGradient(10080.*eta*(-40.*eta*xi + 21.*eta + 28.*xi - 7. - 15.*eta*eta - 25.*xi*xi), -846720.*eta*xi + 188160.*eta + 806400.*eta*(xi*xi) + 352800.*xi + 453600.*xi*(eta*eta) - 49000. - 211680.*eta*eta - 705600.*xi*xi + 67200.*(eta*eta*eta) + 420000.*(xi*xi*xi));
                      case 17:
                        return RealGradient(20160.*eta*(-10.*eta*xi + 3.*eta + 13.*xi - 3. - 10.*xi*xi), -241920.*eta*xi + 26880.*eta + 403200.*eta*(xi*xi) + 161280.*xi - 14560. - 443520.*xi*xi + 336000.*(xi*xi*xi));
                      case 18:
                        return RealGradient(10080.*eta*(10.*eta*xi - 3.*eta - 22.*xi + 4. + 25.*(xi*xi)), 120960.*eta*xi - 13440.*eta - 201600.*eta*xi*xi - 201600.*xi + 18200. + 554400.*(xi*xi) - 420000.*xi*xi*xi);
                      case 19:
                        return RealGradient(0., 6720.*eta - 280. - 30240.*eta*eta + 33600.*(eta*eta*eta));
                      case 20:
                        return RealGradient(0., -13440.*eta + 560. + 60480.*(eta*eta) - 67200.*eta*eta*eta);
                      case 21:
                        return RealGradient(8960.*eta*(-65.*eta*xi/3. + 14.*eta + 35.*xi/3. - 4. - 10.*eta*eta - 20.*xi*xi/3.), -340480.*eta*xi + 56000.*eta + 1164800.*eta*(xi*xi)/3. + 71680.*xi + 268800.*xi*(eta*eta) - 24920./3. - 81760.*eta*eta - 474880.*xi*xi/3. + 280000.*(eta*eta*eta)/9. + 896000.*(xi*xi*xi)/9.);
                      case 22:
                        return RealGradient(2240.*eta*(280.*eta*xi/3. - 49.*eta - 136.*xi/3. + 13. + 35.*(eta*eta) + 100.*(xi*xi)/3.), 371840.*eta*xi - 62720.*eta - 1254400.*eta*xi*xi/3. - 91840.*xi - 235200.*xi*eta*eta + 32480./3. + 76160.*(eta*eta) + 600320.*(xi*xi)/3. - 224000.*eta*eta*eta/9. - 1120000.*xi*xi*xi/9.);
                      case 23:
                        return RealGradient(2240.*eta*(-100.*eta*xi/3. + 34.*eta + 31.*xi/3. - 6. - 35.*eta*eta - 10.*xi*xi/3.), -183680.*eta*xi + 44800.*eta + 448000.*eta*(xi*xi)/3. + 20160.*xi + 235200.*xi*(eta*eta) - 10640./3. - 109760.*eta*eta - 89600.*xi*xi/3. + 627200.*(eta*eta*eta)/9. + 112000.*(xi*xi*xi)/9.);
                      case 24:
                        return RealGradient(1120.*eta*(250.*eta*xi/3. - 73.*eta - 70.*xi/3. + 12. + 80.*(eta*eta) + 25.*(xi*xi)/3.), 232960.*eta*xi - 58240.*eta - 560000.*eta*xi*xi/3. - 26880.*xi - 268800.*xi*eta*eta + 14840./3. + 125440.*(eta*eta) + 115360.*(xi*xi)/3. - 582400.*eta*eta*eta/9. - 140000.*xi*xi*xi/9.);
                      case 25:
                        return RealGradient(2240.*eta*(20.*eta*xi/3. - eta - 44.*xi/3. + 5. - 5.*eta*eta + 20.*(xi*xi)/3.), 4480.*eta*xi + 15680.*eta/3. - 89600.*eta*xi*xi/3. - 11200.*xi + 33600.*xi*(eta*eta) + 6440./9. - 25760.*eta*eta + 98560.*(xi*xi)/3. + 190400.*(eta*eta*eta)/9. - 224000.*xi*xi*xi/9.);
                      case 26:
                        return RealGradient(2240.*eta*(80.*eta*xi/3. - 23.*eta + 64.*xi/3. - 1. + 25.*(eta*eta) - 100.*xi*xi/3.), 165760.*eta*xi - 138880.*eta/3. - 358400.*eta*xi*xi/3. + 64960.*xi - 168000.*xi*eta*eta - 48160./9. + 80640.*(eta*eta) - 519680.*xi*xi/3. - 246400.*eta*eta*eta/9. + 1120000.*(xi*xi*xi)/9.);
                      case 27:
                        return RealGradient(1120.*eta*(-40.*eta*xi/3. + 11.*eta - 2.*xi/3. - 1. - 5.*eta*eta + 5.*(xi*xi)/3.), -24640.*eta*xi + 4480.*eta/3. + 89600.*eta*(xi*xi)/3. + 1120.*xi + 16800.*xi*(eta*eta) - 1400./9. + 5600.*(eta*eta) + 2240.*(xi*xi)/3. - 89600.*eta*eta*eta/9. - 28000.*xi*xi*xi/9.);
                      case 28:
                        return RealGradient(1120.*eta*(200.*eta*xi/3. - 17.*eta - 20.*xi/3. + 3. - 5.*eta*eta - 25.*xi*xi/3.), 116480.*eta*xi - 35840.*eta/3. - 448000.*eta*xi*xi/3. - 3360.*xi + 16800.*xi*(eta*eta) + 8680./9. - 30240.*eta*eta - 24640.*xi*xi/3. + 246400.*(eta*eta*eta)/9. + 140000.*(xi*xi*xi)/9.);
                      case 29:
                        return RealGradient(2240.*eta*(-110.*eta*xi/3. + 36.*eta + 71.*xi/3. - 12. - 25.*eta*eta - 20.*xi*xi/3.), -165760.*eta*xi + 89600.*eta/3. + 492800.*eta*(xi*xi)/3. + 24640.*xi + 168000.*xi*(eta*eta) - 29680./9. - 51520.*eta*eta - 138880.*xi*xi/3. + 179200.*(eta*eta*eta)/9. + 224000.*(xi*xi*xi)/9.);
                      case 30:
                        return RealGradient(1120.*eta*(170.*eta*xi/3. - 23.*eta - 134.*xi/3. + 10. + 10.*(eta*eta) + 125.*(xi*xi)/3.), 103040.*eta*xi - 31360.*eta/3. - 380800.*eta*xi*xi/3. - 56000.*xi - 33600.*xi*eta*eta + 51800./9. - 2240.*eta*eta + 375200.*(xi*xi)/3. + 44800.*(eta*eta*eta)/9. - 700000.*xi*xi*xi/9.);
                      case 31:
                        return RealGradient(1120.*eta*(220.*eta*xi/3. - 27.*eta - 178.*xi/3. + 17. + 5.*(eta*eta) + 85.*(xi*xi)/3.), 105280.*eta*xi - 22400.*eta/3. - 492800.*eta*xi*xi/3. - 32480.*xi - 16800.*xi*eta*eta + 27160./9. - 19040.*eta*eta + 239680.*(xi*xi)/3. + 224000.*(eta*eta*eta)/9. - 476000.*xi*xi*xi/9.);
                      case 32:
                        return RealGradient(1120.*eta*(-80.*eta*xi/3. + 5.*eta + 116.*xi/3. - 7. + 5.*(eta*eta) - 125.*xi*xi/3.), -22400.*eta*xi - 8960.*eta/3. + 179200.*eta*(xi*xi)/3. + 39200.*xi - 16800.*xi*eta*eta - 26600./9. + 12320.*(eta*eta) - 324800.*xi*xi/3. - 44800.*eta*eta*eta/9. + 700000.*(xi*xi*xi)/9.);
                      case 33:
                        return RealGradient(1120.*eta*(20.*eta*xi - 17.*eta - 8.*xi + 3. + 15.*(eta*eta) + 5.*(xi*xi)), 47040.*eta*xi - 29120.*eta/3. - 44800.*eta*xi*xi - 7840.*xi - 50400.*xi*eta*eta + 9520./9. + 20160.*(eta*eta) + 15680.*(xi*xi) - 11200.*eta*eta*eta - 28000.*xi*xi*xi/3.);
                      case 34:
                        return RealGradient(3360.*eta*(-10.*eta*xi + 6.*eta + 2.*xi - 1. - 5.*eta*eta), -67200.*eta*xi + 38080.*eta/3. + 67200.*eta*(xi*xi) + 3360.*xi + 50400.*xi*(eta*eta) - 5600./9. - 19040.*eta*eta - 3360.*xi*xi + 22400.*(eta*eta*eta)/3.);
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
                        return sign * RealGradient(45360.*eta*xi - 16800.*eta - 25200.*eta*xi*xi - 12600.*xi - 37800.*xi*eta*eta + 2800. + 30240.*(eta*eta) + 15120.*(xi*xi) - 16800.*eta*eta*eta - 4200.*xi*xi*xi, -30240.*eta*xi + 4200.*eta + 37800.*eta*(xi*xi) + 8400.*xi + 25200.*xi*(eta*eta) - 700. - 7560.*eta*eta - 22680.*xi*xi + 4200.*(eta*eta*eta) + 16800.*(xi*xi*xi));
                      case 1:
                        return sign * RealGradient(-240975.*eta*xi/16. + 7245.*eta/2. + 146475.*eta*(xi*xi)/16. + 136185.*xi/32. + 382725.*xi*(eta*eta)/32. - 12145./16. - 76545.*eta*eta/16. - 20475.*xi*xi/4. + 14175.*(eta*eta*eta)/8. + 19425.*(xi*xi*xi)/16., 25515.*eta*xi/4. + 2835.*eta/32. - 382725.*eta*xi*xi/32. - 39375.*xi/16. - 42525.*xi*eta*eta/16. + 3045./32. - 42525.*eta*eta/32. + 252315.*(xi*xi)/32. + 42525.*(eta*eta*eta)/32. - 48825.*xi*xi*xi/8.);
                      case 2:
                        return sign * RealGradient(945.*eta*xi + 1680.*eta - 4725.*eta*xi*xi - 3255.*xi/2. + 4725.*xi*(eta*eta)/2. - 4725.*eta*eta + 3780.*(xi*xi) + 3150.*(eta*eta*eta) - 1575.*xi*xi*xi, 3780.*eta*xi - 525.*eta/2. - 4725.*eta*xi*xi/2. + 105.*xi - 4725.*xi*eta*eta + 35./2. + 315.*(eta*eta)/2. - 4725.*xi*xi/2. + 525.*(eta*eta*eta)/2. + 3150.*(xi*xi*xi));
                      case 3:
                        return sign * RealGradient(76545.*eta*xi/16. - 1155.*eta/2. - 29925.*eta*xi*xi/16. - 16695.*xi/32. - 146475.*xi*eta*eta/32. + 1855./16. - 945.*eta*eta/16. - 2835.*xi*xi/4. + 5775.*(eta*eta*eta)/8. + 19425.*(xi*xi*xi)/16., -4725.*eta*xi/4. + 1155.*eta/32. + 146475.*eta*(xi*xi)/32. + 9345.*xi/16. - 17325.*xi*eta*eta/16. - 875./32. + 2835.*(eta*eta)/32. - 65205.*xi*xi/32. + 525.*(eta*eta*eta)/32. + 9975.*(xi*xi*xi)/8.);
                      case 4:
                        return sign * RealGradient(-1680.*xi + 140. + 5040.*(xi*xi) - 4200.*xi*xi*xi, 0.);
                      case 5:
                        return sign * RealGradient(-1680.*xi + 140. + 5040.*(xi*xi) - 4200.*xi*xi*xi, 0.);
                      case 6:
                        return sign * RealGradient(76545.*eta*xi/8. - 1575.*eta - 42525.*eta*xi*xi/4. - 2835.*xi/2. - 127575.*xi*eta*eta/16. + 735./4. + 36855.*(eta*eta)/16. + 42525.*(xi*xi)/16. - 1575.*eta*eta*eta/2. - 42525.*xi*xi*xi/32., -31185.*eta*xi/8. + 525.*eta/2. + 127575.*eta*(xi*xi)/16. + 2205.*xi + 4725.*xi*(eta*eta)/4. - 455./4. - 2205.*eta*eta/16. - 127575.*xi*xi/16. + 525.*(eta*eta*eta)/32. + 14175.*(xi*xi*xi)/2.);
                      case 7:
                        return sign * RealGradient(9450.*eta*xi - 2520.*eta - 6300.*eta*xi*xi - 840.*xi - 14175.*xi*eta*eta + 175. + 6615.*(eta*eta) + 945.*(xi*xi) - 4200.*eta*eta*eta - 525.*xi*xi*xi/2., -9450.*eta*xi + 840.*eta + 14175.*eta*(xi*xi) + 2520.*xi + 6300.*xi*(eta*eta) - 175. - 945.*eta*eta - 6615.*xi*xi + 525.*(eta*eta*eta)/2. + 4200.*(xi*xi*xi));
                      case 8:
                        return sign * RealGradient(31185.*eta*xi/8. - 2205.*eta - 4725.*eta*xi*xi/4. - 525.*xi/2. - 127575.*xi*eta*eta/16. + 455./4. + 127575.*(eta*eta)/16. + 2205.*(xi*xi)/16. - 14175.*eta*eta*eta/2. - 525.*xi*xi*xi/32., -76545.*eta*xi/8. + 2835.*eta/2. + 127575.*eta*(xi*xi)/16. + 1575.*xi + 42525.*xi*(eta*eta)/4. - 735./4. - 42525.*eta*eta/16. - 36855.*xi*xi/16. + 42525.*(eta*eta*eta)/32. + 1575.*(xi*xi*xi)/2.);
                      case 9:
                        return sign * RealGradient(0., 1680.*eta - 140. - 5040.*eta*eta + 4200.*(eta*eta*eta));
                      case 10:
                        return sign * RealGradient(0., 1680.*eta - 140. - 5040.*eta*eta + 4200.*(eta*eta*eta));
                      case 11:
                        return sign * RealGradient(4725.*eta*xi/4. - 9345.*eta/16. + 17325.*eta*(xi*xi)/16. - 1155.*xi/32. - 146475.*xi*eta*eta/32. + 875./32. + 65205.*(eta*eta)/32. - 2835.*xi*xi/32. - 9975.*eta*eta*eta/8. - 525.*xi*xi*xi/32., -76545.*eta*xi/16. + 16695.*eta/32. + 146475.*eta*(xi*xi)/32. + 1155.*xi/2. + 29925.*xi*(eta*eta)/16. - 1855./16. + 2835.*(eta*eta)/4. + 945.*(xi*xi)/16. - 19425.*eta*eta*eta/16. - 5775.*xi*xi*xi/8.);
                      case 12:
                        return sign * RealGradient(-3780.*eta*xi - 105.*eta + 4725.*eta*(xi*xi) + 525.*xi/2. + 4725.*xi*(eta*eta)/2. - 35./2. + 4725.*(eta*eta)/2. - 315.*xi*xi/2. - 3150.*eta*eta*eta - 525.*xi*xi*xi/2., -945.*eta*xi + 3255.*eta/2. - 4725.*eta*xi*xi/2. - 1680.*xi + 4725.*xi*(eta*eta) - 3780.*eta*eta + 4725.*(xi*xi) + 1575.*(eta*eta*eta) - 3150.*xi*xi*xi);
                      case 13:
                        return sign * RealGradient(-25515.*eta*xi/4. + 39375.*eta/16. + 42525.*eta*(xi*xi)/16. - 2835.*xi/32. + 382725.*xi*(eta*eta)/32. - 3045./32. - 252315.*eta*eta/32. + 42525.*(xi*xi)/32. + 48825.*(eta*eta*eta)/8. - 42525.*xi*xi*xi/32., 240975.*eta*xi/16. - 136185.*eta/32. - 382725.*eta*xi*xi/32. - 7245.*xi/2. - 146475.*xi*eta*eta/16. + 12145./16. + 20475.*(eta*eta)/4. + 76545.*(xi*xi)/16. - 19425.*eta*eta*eta/16. - 14175.*xi*xi*xi/8.);
                      case 14:
                        return sign * RealGradient(30240.*eta*xi - 8400.*eta - 25200.*eta*xi*xi - 4200.*xi - 37800.*xi*eta*eta + 700. + 22680.*(eta*eta) + 7560.*(xi*xi) - 16800.*eta*eta*eta - 4200.*xi*xi*xi, -45360.*eta*xi + 12600.*eta + 37800.*eta*(xi*xi) + 16800.*xi + 25200.*xi*(eta*eta) - 2800. - 15120.*eta*eta - 30240.*xi*xi + 4200.*(eta*eta*eta) + 16800.*(xi*xi*xi));
                      case 15:
                        return RealGradient(-423360.*eta*xi + 188160.*eta + 201600.*eta*(xi*xi) + 70560.*xi + 453600.*xi*(eta*eta) - 19600. - 423360.*eta*eta - 70560.*xi*xi + 268800.*(eta*eta*eta) + 16800.*(xi*xi*xi), 423360.*eta*xi - 70560.*eta - 453600.*eta*xi*xi - 94080.*xi - 403200.*xi*eta*eta + 9800. + 141120.*(eta*eta) + 211680.*(xi*xi) - 84000.*eta*eta*eta - 134400.*xi*xi*xi);
                      case 16:
                        return RealGradient(423360.*eta*xi - 94080.*eta - 403200.*eta*xi*xi - 70560.*xi - 453600.*xi*eta*eta + 9800. + 211680.*(eta*eta) + 141120.*(xi*xi) - 134400.*eta*eta*eta - 84000.*xi*xi*xi, -423360.*eta*xi + 70560.*eta + 453600.*eta*(xi*xi) + 188160.*xi + 201600.*xi*(eta*eta) - 19600. - 70560.*eta*eta - 423360.*xi*xi + 16800.*(eta*eta*eta) + 268800.*(xi*xi*xi));
                      case 17:
                        return RealGradient(120960.*eta*xi - 13440.*eta - 201600.*eta*xi*xi - 60480.*xi + 6440. + 131040.*(xi*xi) - 67200.*xi*xi*xi, 26880.*xi - 1120. - 120960.*xi*xi + 134400.*(xi*xi*xi));
                      case 18:
                        return RealGradient(-60480.*eta*xi + 6720.*eta + 100800.*eta*(xi*xi) + 40320.*xi - 3640. - 110880.*xi*xi + 84000.*(xi*xi*xi), -13440.*xi + 560. + 60480.*(xi*xi) - 67200.*xi*xi*xi);
                      case 19:
                        return RealGradient(-13440.*eta + 560. + 60480.*(eta*eta) - 67200.*eta*eta*eta, -60480.*eta*xi + 40320.*eta + 6720.*xi + 100800.*xi*(eta*eta) - 3640. - 110880.*eta*eta + 84000.*(eta*eta*eta));
                      case 20:
                        return RealGradient(26880.*eta - 1120. - 120960.*eta*eta + 134400.*(eta*eta*eta), 120960.*eta*xi - 60480.*eta - 13440.*xi - 201600.*xi*eta*eta + 6440. + 131040.*(eta*eta) - 67200.*eta*eta*eta);
                      case 21:
                        return RealGradient(250880.*eta*xi - 58240.*eta - 582400.*eta*xi*xi/3. - 35840.*xi - 268800.*xi*eta*eta + 17920./3. + 116480.*(eta*eta) + 156800.*(xi*xi)/3. - 560000.*eta*eta*eta/9. - 179200.*xi*xi*xi/9., -163520.*eta*xi + 13440.*eta + 268800.*eta*(xi*xi) + 56000.*xi + 280000.*xi*(eta*eta)/3. - 10360./3. - 39200.*eta*eta/3. - 170240.*xi*xi + 28000.*(eta*eta*eta)/9. + 1164800.*(xi*xi*xi)/9.);
                      case 22:
                        return RealGradient(-219520.*eta*xi + 44800.*eta + 627200.*eta*(xi*xi)/3. + 29120.*xi + 235200.*xi*(eta*eta) - 12880./3. - 91840.*eta*eta - 152320.*xi*xi/3. + 448000.*(eta*eta*eta)/9. + 224000.*(xi*xi*xi)/9., 152320.*eta*xi - 13440.*eta - 235200.*eta*xi*xi - 62720.*xi - 224000.*xi*eta*eta/3. + 12040./3. + 34720.*(eta*eta)/3. + 185920.*(xi*xi) - 22400.*eta*eta*eta/9. - 1254400.*xi*xi*xi/9.);
                      case 23:
                        return RealGradient(152320.*eta*xi - 62720.*eta - 224000.*eta*xi*xi/3. - 13440.*xi - 235200.*xi*eta*eta + 12040./3. + 185920.*(eta*eta) + 34720.*(xi*xi)/3. - 1254400.*eta*eta*eta/9. - 22400.*xi*xi*xi/9., -219520.*eta*xi + 29120.*eta + 235200.*eta*(xi*xi) + 44800.*xi + 627200.*xi*(eta*eta)/3. - 12880./3. - 152320.*eta*eta/3. - 91840.*xi*xi + 224000.*(eta*eta*eta)/9. + 448000.*(xi*xi*xi)/9.);
                      case 24:
                        return RealGradient(-163520.*eta*xi + 56000.*eta + 280000.*eta*(xi*xi)/3. + 13440.*xi + 268800.*xi*(eta*eta) - 10360./3. - 170240.*eta*eta - 39200.*xi*xi/3. + 1164800.*(eta*eta*eta)/9. + 28000.*(xi*xi*xi)/9., 250880.*eta*xi - 35840.*eta - 268800.*eta*xi*xi - 58240.*xi - 582400.*xi*eta*eta/3. + 17920./3. + 156800.*(eta*eta)/3. + 116480.*(xi*xi) - 179200.*eta*eta*eta/9. - 560000.*xi*xi*xi/9.);
                      case 25:
                        return RealGradient(-4480.*eta*xi - 31360.*eta/3. + 44800.*eta*(xi*xi)/3. + 11200.*xi - 33600.*xi*eta*eta - 12880./9. + 51520.*(eta*eta) - 49280.*xi*xi/3. - 380800.*eta*eta*eta/9. + 44800.*(xi*xi*xi)/9., -51520.*eta*xi + 11200.*eta + 33600.*eta*(xi*xi) + 15680.*xi/3. + 190400.*xi*(eta*eta)/3. - 10360./9. - 75040.*eta*eta/3. + 2240.*(xi*xi) + 140000.*(eta*eta*eta)/9. - 89600.*xi*xi*xi/9.);
                      case 26:
                        return RealGradient(-103040.*eta*xi + 89600.*eta/3. + 179200.*eta*(xi*xi)/3. - 2240.*xi + 168000.*xi*(eta*eta) - 9520./9. - 82880.*eta*eta + 71680.*(xi*xi)/3. + 492800.*(eta*eta*eta)/9. - 224000.*xi*xi*xi/9., 161280.*eta*xi - 26880.*eta - 168000.*eta*xi*xi - 138880.*xi/3. - 246400.*xi*eta*eta/3. + 51800./9. + 79520.*(eta*eta)/3. + 82880.*(xi*xi) - 44800.*eta*eta*eta/9. - 358400.*xi*xi*xi/9.);
                      case 27:
                        return RealGradient(24640.*eta*xi - 8960.*eta/3. - 44800.*eta*xi*xi/3. - 1120.*xi - 16800.*xi*eta*eta + 2800./9. - 11200.*eta*eta - 1120.*xi*xi/3. + 179200.*(eta*eta*eta)/9. + 5600.*(xi*xi*xi)/9., 11200.*eta*xi - 7840.*eta + 16800.*eta*(xi*xi) + 4480.*xi/3. - 89600.*xi*eta*eta/3. + 5320./9. + 64960.*(eta*eta)/3. - 12320.*xi*xi - 140000.*eta*eta*eta/9. + 89600.*(xi*xi*xi)/9.);
                      case 28:
                        return RealGradient(-38080.*eta*xi - 22400.*eta/3. + 224000.*eta*(xi*xi)/3. + 3360.*xi - 16800.*xi*eta*eta + 280./9. + 52640.*(eta*eta) - 11200.*xi*xi/3. - 492800.*eta*eta*eta/9. - 28000.*xi*xi*xi/9., -60480.*eta*xi + 19040.*eta + 16800.*eta*(xi*xi) - 35840.*xi/3. + 246400.*xi*(eta*eta)/3. - 8960./9. - 99680.*eta*eta/3. + 58240.*(xi*xi) + 95200.*(eta*eta*eta)/9. - 448000.*xi*xi*xi/9.);
                      case 29:
                        return RealGradient(161280.*eta*xi - 138880.*eta/3. - 246400.*eta*xi*xi/3. - 26880.*xi - 168000.*xi*eta*eta + 51800./9. + 82880.*(eta*eta) + 79520.*(xi*xi)/3. - 358400.*eta*eta*eta/9. - 44800.*xi*xi*xi/9., -103040.*eta*xi - 2240.*eta + 168000.*eta*(xi*xi) + 89600.*xi/3. + 179200.*xi*(eta*eta)/3. - 9520./9. + 71680.*(eta*eta)/3. - 82880.*xi*xi - 224000.*eta*eta*eta/9. + 492800.*(xi*xi*xi)/9.);
                      case 30:
                        return RealGradient(-51520.*eta*xi + 15680.*eta/3. + 190400.*eta*(xi*xi)/3. + 11200.*xi + 33600.*xi*(eta*eta) - 10360./9. + 2240.*(eta*eta) - 75040.*xi*xi/3. - 89600.*eta*eta*eta/9. + 140000.*(xi*xi*xi)/9., -4480.*eta*xi + 11200.*eta - 33600.*eta*xi*xi - 31360.*xi/3. + 44800.*xi*(eta*eta)/3. - 12880./9. - 49280.*eta*eta/3. + 51520.*(xi*xi) + 44800.*(eta*eta*eta)/9. - 380800.*xi*xi*xi/9.);
                      case 31:
                        return RealGradient(-60480.*eta*xi - 35840.*eta/3. + 246400.*eta*(xi*xi)/3. + 19040.*xi + 16800.*xi*(eta*eta) - 8960./9. + 58240.*(eta*eta) - 99680.*xi*xi/3. - 448000.*eta*eta*eta/9. + 95200.*(xi*xi*xi)/9., -38080.*eta*xi + 3360.*eta - 16800.*eta*xi*xi - 22400.*xi/3. + 224000.*xi*(eta*eta)/3. + 280./9. - 11200.*eta*eta/3. + 52640.*(xi*xi) - 28000.*eta*eta*eta/9. - 492800.*xi*xi*xi/9.);
                      case 32:
                        return RealGradient(11200.*eta*xi + 4480.*eta/3. - 89600.*eta*xi*xi/3. - 7840.*xi + 16800.*xi*(eta*eta) + 5320./9. - 12320.*eta*eta + 64960.*(xi*xi)/3. + 89600.*(eta*eta*eta)/9. - 140000.*xi*xi*xi/9., 24640.*eta*xi - 1120.*eta - 16800.*eta*xi*xi - 8960.*xi/3. - 44800.*xi*eta*eta/3. + 2800./9. - 1120.*eta*eta/3. - 11200.*xi*xi + 5600.*(eta*eta*eta)/9. + 179200.*(xi*xi*xi)/9.);
                      case 33:
                        return RealGradient(-38080.*eta*xi + 38080.*eta/3. + 22400.*eta*(xi*xi) + 3360.*xi + 50400.*xi*(eta*eta) - 6440./9. - 33600.*eta*eta - 4480.*xi*xi + 22400.*(eta*eta*eta) + 5600.*(xi*xi*xi)/3., 40320.*eta*xi - 3360.*eta - 50400.*eta*xi*xi - 29120.*xi/3. - 33600.*xi*eta*eta + 6160./9. + 3360.*(eta*eta) + 23520.*(xi*xi) - 44800.*xi*xi*xi/3.);
                      case 34:
                        return RealGradient(40320.*eta*xi - 29120.*eta/3. - 33600.*eta*xi*xi - 3360.*xi - 50400.*xi*eta*eta + 6160./9. + 23520.*(eta*eta) + 3360.*(xi*xi) - 44800.*eta*eta*eta/3., -38080.*eta*xi + 3360.*eta + 50400.*eta*(xi*xi) + 38080.*xi/3. + 22400.*xi*(eta*eta) - 6440./9. - 4480.*eta*eta - 33600.*xi*xi + 5600.*(eta*eta*eta)/3. + 22400.*(xi*xi*xi));
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
                        return sign * RealGradient(60480.*eta*xi - 21000.*eta - 37800.*eta*xi*xi - 16800.*xi - 50400.*xi*eta*eta + 3500. + 37800.*(eta*eta) + 22680.*(xi*xi) - 21000.*eta*eta*eta - 8400.*xi*xi*xi, 840.*xi*(30.*eta*xi - 18.*eta - 18.*xi + 5. + 15.*(eta*eta) + 15.*(xi*xi)));
                      case 1:
                        return sign * RealGradient(-76545.*eta*xi/8. - 127575.*eta/32. + 382725.*eta*(xi*xi)/32. + 7245.*xi/2. + 42525.*xi*(eta*eta)/8. + 11235./32. + 314685.*(eta*eta)/32. - 240975.*xi*xi/32. - 212625.*eta*eta*eta/32. + 48825.*(xi*xi*xi)/16., 2835.*xi*(-30.*eta*xi - 30.*eta + 36.*xi + 1. + 45.*(eta*eta) - 45.*xi*xi)/32.);
                      case 2:
                        return sign * RealGradient(-9450.*eta*xi + 105.*eta/2. + 4725.*eta*(xi*xi)/2. + 1680.*xi + 9450.*xi*(eta*eta) - 175./2. + 2205.*(eta*eta)/2. + 945.*(xi*xi)/2. - 2625.*eta*eta*eta/2. - 1575.*xi*xi*xi, 105.*xi*(-45.*eta*xi + 3.*eta + 18.*xi - 5./2. + 15.*(eta*eta)/2. - 15.*xi*xi/2.));
                      case 3:
                        return sign * RealGradient(-945.*eta*xi/8. + 1785.*eta/32. - 146475.*eta*xi*xi/32. - 1155.*xi/2. + 17325.*xi*(eta*eta)/8. + 595./32. - 2835.*eta*eta/32. + 76545.*(xi*xi)/32. - 2625.*eta*eta*eta/32. - 9975.*xi*xi*xi/16., 105.*xi*(-330.*eta*xi + 54.*eta - 180.*xi + 11. + 15.*(eta*eta) + 465.*(xi*xi))/32.);
                      case 4:
                        return sign * RealGradient(0., 0.);
                      case 5:
                        return sign * RealGradient(0., 0.);
                      case 6:
                        return sign * RealGradient(36855.*eta*xi/8. - 735.*eta/2. - 127575.*eta*xi*xi/16. - 1575.*xi - 4725.*xi*eta*eta/2. + 385./4. + 5355.*(eta*eta)/16. + 76545.*(xi*xi)/16. - 2625.*eta*eta*eta/32. - 14175.*xi*xi*xi/4., 105.*xi*(360.*eta*xi - 84.*eta - 594.*xi + 80. + 15.*(eta*eta) + 810.*(xi*xi))/32.);
                      case 7:
                        return sign * RealGradient(13230.*eta*xi - 1680.*eta - 14175.*eta*xi*xi - 2520.*xi - 12600.*xi*eta*eta + 245. + 2835.*(eta*eta) + 4725.*(xi*xi) - 2625.*eta*eta*eta/2. - 2100.*xi*xi*xi, 105.*xi*(60.*eta*xi - 18.*eta - 45.*xi + 8. + 15.*(eta*eta)/2. + 45.*(xi*xi)));
                      case 8:
                        return sign * RealGradient(127575.*eta*xi/8. - 8505.*eta/2. - 127575.*eta*xi*xi/16. - 2205.*xi - 42525.*xi*eta*eta/2. + 1785./4. + 161595.*(eta*eta)/16. + 31185.*(xi*xi)/16. - 212625.*eta*eta*eta/32. - 1575.*xi*xi*xi/4., 2835.*xi*(120.*eta*xi - 60.*eta - 54.*xi + 16. + 45.*(eta*eta) + 30.*(xi*xi))/32.);
                      case 9:
                        return sign * RealGradient(-8400.*eta + 700. + 25200.*(eta*eta) - 21000.*eta*eta*eta, 840.*xi*(-12.*eta + 2. + 15.*(eta*eta)));
                      case 10:
                        return sign * RealGradient(-8400.*eta + 700. + 25200.*(eta*eta) - 21000.*eta*eta*eta, -10080.*eta*xi + 6720.*eta + 1680.*xi + 12600.*xi*(eta*eta) - 1260. - 7560.*eta*eta);
                      case 11:
                        return sign * RealGradient(65205.*eta*xi/16. + 66465.*eta/32. - 146475.*eta*xi*xi/32. - 9345.*xi/16. - 29925.*xi*eta*eta/8. - 595./4. - 110565.*eta*eta/16. + 4725.*(xi*xi)/8. + 97125.*(eta*eta*eta)/16. + 5775.*(xi*xi*xi)/16., 2835.*eta*xi/2. - 29085.*eta/16. + 29925.*eta*(xi*xi)/16. + 16695.*xi/32. - 58275.*xi*eta*eta/16. + 9345./32. + 34965.*(eta*eta)/16. - 76545.*xi*xi/32. + 48825.*(xi*xi*xi)/32.);
                      case 12:
                        return sign * RealGradient(4725.*eta*xi - 7455.*eta/2. + 4725.*eta*(xi*xi)/2. - 105.*xi - 9450.*xi*eta*eta + 315. + 10395.*(eta*eta) - 1890.*xi*xi - 7875.*eta*eta*eta + 1575.*(xi*xi*xi), -7560.*eta*xi + 2835.*eta + 4725.*eta*(xi*xi) + 3255.*xi/2. + 4725.*xi*(eta*eta) - 1155./2. - 2835.*eta*eta - 945.*xi*xi/2. - 1575.*xi*xi*xi/2.);
                      case 13:
                        return sign * RealGradient(-252315.*eta*xi/16. + 125265.*eta/32. + 382725.*eta*(xi*xi)/32. + 39375.*xi/16. + 146475.*xi*(eta*eta)/8. - 1715./4. - 145845.*eta*eta/16. - 25515.*xi*xi/8. + 97125.*(eta*eta*eta)/16. + 14175.*(xi*xi*xi)/16., 20475.*eta*xi/2. - 40845.*eta/16. - 146475.*eta*xi*xi/16. - 136185.*xi/32. - 58275.*xi*eta*eta/16. + 21105./32. + 34965.*(eta*eta)/16. + 240975.*(xi*xi)/32. - 127575.*xi*xi*xi/32.);
                      case 14:
                        return sign * RealGradient(45360.*eta*xi - 12600.*eta - 37800.*eta*xi*xi - 8400.*xi - 50400.*xi*eta*eta + 1400. + 30240.*(eta*eta) + 15120.*(xi*xi) - 21000.*eta*eta*eta - 8400.*xi*xi*xi, -30240.*eta*xi + 8400.*eta + 25200.*eta*(xi*xi) + 12600.*xi + 12600.*xi*(eta*eta) - 2100. - 7560.*eta*eta - 22680.*xi*xi + 12600.*(xi*xi*xi));
                      case 15:
                        return RealGradient(-846720.*eta*xi + 352800.*eta + 453600.*eta*(xi*xi) + 188160.*xi + 806400.*xi*(eta*eta) - 49000. - 705600.*eta*eta - 211680.*xi*xi + 420000.*(eta*eta*eta) + 67200.*(xi*xi*xi), 10080.*xi*(-40.*eta*xi + 28.*eta + 21.*xi - 7. - 25.*eta*eta - 15.*xi*xi));
                      case 16:
                        return RealGradient(423360.*eta*xi - 70560.*eta - 453600.*eta*xi*xi - 94080.*xi - 403200.*xi*eta*eta + 9800. + 141120.*(eta*eta) + 211680.*(xi*xi) - 84000.*eta*eta*eta - 134400.*xi*xi*xi, 10080.*xi*(20.*eta*xi - 14.*eta - 21.*xi + 7. + 5.*(eta*eta) + 15.*(xi*xi)));
                      case 17:
                        return RealGradient(-13440.*xi + 560. + 60480.*(xi*xi) - 67200.*xi*xi*xi, 0.);
                      case 18:
                        return RealGradient(6720.*xi - 280. - 30240.*xi*xi + 33600.*(xi*xi*xi), 0.);
                      case 19:
                        return RealGradient(120960.*eta*xi - 201600.*eta - 13440.*xi - 201600.*xi*eta*eta + 18200. + 554400.*(eta*eta) - 420000.*eta*eta*eta, 10080.*xi*(10.*eta*xi - 22.*eta - 3.*xi + 4. + 25.*(eta*eta)));
                      case 20:
                        return RealGradient(-241920.*eta*xi + 161280.*eta + 26880.*xi + 403200.*xi*(eta*eta) - 14560. - 443520.*eta*eta + 336000.*(eta*eta*eta), 20160.*xi*(-10.*eta*xi + 13.*eta + 3.*xi - 3. - 10.*eta*eta));
                      case 21:
                        return RealGradient(232960.*eta*xi - 26880.*eta - 268800.*eta*xi*xi - 58240.*xi - 560000.*xi*eta*eta/3. + 14840./3. + 115360.*(eta*eta)/3. + 125440.*(xi*xi) - 140000.*eta*eta*eta/9. - 582400.*xi*xi*xi/9., 1120.*xi*(250.*eta*xi/3. - 70.*eta/3. - 73.*xi + 12. + 25.*(eta*eta)/3. + 80.*(xi*xi)));
                      case 22:
                        return RealGradient(-183680.*eta*xi + 20160.*eta + 235200.*eta*(xi*xi) + 44800.*xi + 448000.*xi*(eta*eta)/3. - 10640./3. - 89600.*eta*eta/3. - 109760.*xi*xi + 112000.*(eta*eta*eta)/9. + 627200.*(xi*xi*xi)/9., 2240.*xi*(-100.*eta*xi/3. + 31.*eta/3. + 34.*xi - 6. - 10.*eta*eta/3. - 35.*xi*xi));
                      case 23:
                        return RealGradient(371840.*eta*xi - 91840.*eta - 235200.*eta*xi*xi - 62720.*xi - 1254400.*xi*eta*eta/3. + 32480./3. + 600320.*(eta*eta)/3. + 76160.*(xi*xi) - 1120000.*eta*eta*eta/9. - 224000.*xi*xi*xi/9., 2240.*xi*(280.*eta*xi/3. - 136.*eta/3. - 49.*xi + 13. + 100.*(eta*eta)/3. + 35.*(xi*xi)));
                      case 24:
                        return RealGradient(-340480.*eta*xi + 71680.*eta + 268800.*eta*(xi*xi) + 56000.*xi + 1164800.*xi*(eta*eta)/3. - 24920./3. - 474880.*eta*eta/3. - 81760.*xi*xi + 896000.*(eta*eta*eta)/9. + 280000.*(xi*xi*xi)/9., 8960.*xi*(-65.*eta*xi/3. + 35.*eta/3. + 14.*xi - 4. - 20.*eta*eta/3. - 10.*xi*xi));
                      case 25:
                        return RealGradient(103040.*eta*xi - 56000.*eta - 33600.*eta*xi*xi - 31360.*xi/3. - 380800.*xi*eta*eta/3. + 51800./9. + 375200.*(eta*eta)/3. - 2240.*xi*xi - 700000.*eta*eta*eta/9. + 44800.*(xi*xi*xi)/9., 1120.*xi*(170.*eta*xi/3. - 134.*eta/3. - 23.*xi + 10. + 125.*(eta*eta)/3. + 10.*(xi*xi)));
                      case 26:
                        return RealGradient(-165760.*eta*xi + 24640.*eta + 168000.*eta*(xi*xi) + 89600.*xi/3. + 492800.*xi*(eta*eta)/3. - 29680./9. - 138880.*eta*eta/3. - 51520.*xi*xi + 224000.*(eta*eta*eta)/9. + 179200.*(xi*xi*xi)/9., 2240.*xi*(-110.*eta*xi/3. + 71.*eta/3. + 36.*xi - 12. - 20.*eta*eta/3. - 25.*xi*xi));
                      case 27:
                        return RealGradient(-22400.*eta*xi + 39200.*eta - 16800.*eta*xi*xi - 8960.*xi/3. + 179200.*xi*(eta*eta)/3. - 26600./9. - 324800.*eta*eta/3. + 12320.*(xi*xi) + 700000.*(eta*eta*eta)/9. - 44800.*xi*xi*xi/9., 1120.*xi*(-80.*eta*xi/3. + 116.*eta/3. + 5.*xi - 7. - 125.*eta*eta/3. + 5.*(xi*xi)));
                      case 28:
                        return RealGradient(105280.*eta*xi - 32480.*eta - 16800.*eta*xi*xi - 22400.*xi/3. - 492800.*xi*eta*eta/3. + 27160./9. + 239680.*(eta*eta)/3. - 19040.*xi*xi - 476000.*eta*eta*eta/9. + 224000.*(xi*xi*xi)/9., 1120.*xi*(220.*eta*xi/3. - 178.*eta/3. - 27.*xi + 17. + 85.*(eta*eta)/3. + 5.*(xi*xi)));
                      case 29:
                        return RealGradient(165760.*eta*xi + 64960.*eta - 168000.*eta*xi*xi - 138880.*xi/3. - 358400.*xi*eta*eta/3. - 48160./9. - 519680.*eta*eta/3. + 80640.*(xi*xi) + 1120000.*(eta*eta*eta)/9. - 246400.*xi*xi*xi/9., 2240.*xi*(80.*eta*xi/3. + 64.*eta/3. - 23.*xi - 1. - 100.*eta*eta/3. + 25.*(xi*xi)));
                      case 30:
                        return RealGradient(4480.*eta*xi - 11200.*eta + 33600.*eta*(xi*xi) + 15680.*xi/3. - 89600.*xi*eta*eta/3. + 6440./9. + 98560.*(eta*eta)/3. - 25760.*xi*xi - 224000.*eta*eta*eta/9. + 190400.*(xi*xi*xi)/9., 2240.*xi*(20.*eta*xi/3. - 44.*eta/3. - xi + 5. + 20.*(eta*eta)/3. - 5.*xi*xi));
                      case 31:
                        return RealGradient(116480.*eta*xi - 3360.*eta + 16800.*eta*(xi*xi) - 35840.*xi/3. - 448000.*xi*eta*eta/3. + 8680./9. - 24640.*eta*eta/3. - 30240.*xi*xi + 140000.*(eta*eta*eta)/9. + 246400.*(xi*xi*xi)/9., 1120.*xi*(200.*eta*xi/3. - 20.*eta/3. - 17.*xi + 3. - 25.*eta*eta/3. - 5.*xi*xi));
                      case 32:
                        return RealGradient(-24640.*eta*xi + 1120.*eta + 16800.*eta*(xi*xi) + 4480.*xi/3. + 89600.*xi*(eta*eta)/3. - 1400./9. + 2240.*(eta*eta)/3. + 5600.*(xi*xi) - 28000.*eta*eta*eta/9. - 89600.*xi*xi*xi/9., 1120.*xi*(-40.*eta*xi/3. - 2.*eta/3. + 11.*xi - 1. + 5.*(eta*eta)/3. - 5.*xi*xi));
                      case 33:
                        return RealGradient(-67200.*eta*xi + 3360.*eta + 50400.*eta*(xi*xi) + 38080.*xi/3. + 67200.*xi*(eta*eta) - 5600./9. - 3360.*eta*eta - 19040.*xi*xi + 22400.*(xi*xi*xi)/3., 3360.*xi*(-10.*eta*xi + 2.*eta + 6.*xi - 1. - 5.*xi*xi));
                      case 34:
                        return RealGradient(47040.*eta*xi - 7840.*eta - 50400.*eta*xi*xi - 29120.*xi/3. - 44800.*xi*eta*eta + 9520./9. + 15680.*(eta*eta) + 20160.*(xi*xi) - 28000.*eta*eta*eta/3. - 11200.*xi*xi*xi, 1120.*xi*(20.*eta*xi - 8.*eta - 17.*xi + 3. + 5.*(eta*eta) + 15.*(xi*xi)));
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
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << totalorder);

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
