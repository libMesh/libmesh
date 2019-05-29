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



// Local includes
#include "libmesh/quadrature_gm.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

void QGrundmann_Moller::init_3D(const ElemType, unsigned int)
{
  // Nearly all GM rules contain negative weights, so if you are not
  // allowing rules with negative weights, we cannot continue!
  if (!allow_rules_with_negative_weights)
    libmesh_error_msg("You requested a Grundmann-Moller rule but\n"        \
                      << "are not allowing rules with negative weights!\n" \
                      << "Either select a different quadrature class or\n" \
                      << "set allow_rules_with_negative_weights==true.");

  switch (_type)
    {
    case TET4:
    case TET10:
      {
        switch(get_order())
          {
            // We hard-code the first few orders based on output from
            // the mp-quadrature library:
            //
            // https://code.google.com/p/mp-quadrature
            //
            // The points are ordered in such a way that the amount of
            // round-off error in the quadrature calculations is
            // (hopefully) minimized.  These orderings were obtained
            // via a simple permutation optimization strategy designed
            // to produce the smallest overall average error while
            // integrating all polynomials of degree <= d.
          case CONSTANT:
          case FIRST:
            {
              _points.resize(1);
              _weights.resize(1);

              _points[0](0) = Real(1)/4;
              _points[0](1) = Real(1)/4;
              _points[0](2) = Real(1)/4;

              _weights[0] = Real(1)/6;
              return;
            }

          case SECOND:
          case THIRD:
            {
              _points.resize(5);
              _weights.resize(5);

              _points[0](0) = Real(1)/6;  _points[0](1) = Real(1)/6;  _points[0](2) = Real(1)/2;  _weights[0] =  Real(3)/40;
              _points[1](0) = Real(1)/6;  _points[1](1) = Real(1)/6;  _points[1](2) = Real(1)/6;  _weights[1] =  Real(3)/40;
              _points[2](0) = Real(1)/4;  _points[2](1) = Real(1)/4;  _points[2](2) = Real(1)/4;  _weights[2] = -Real(2)/15;
              _points[3](0) = Real(1)/2;  _points[3](1) = Real(1)/6;  _points[3](2) = Real(1)/6;  _weights[3] =  Real(3)/40;
              _points[4](0) = Real(1)/6;  _points[4](1) = Real(1)/2;  _points[4](2) = Real(1)/6;  _weights[4] =  Real(3)/40;
              return;
            }

          case FOURTH:
          case FIFTH:
            {
              _points.resize(15);
              _weights.resize(15);

              _points[ 0](0) = Real(1)/8;  _points[ 0](1) = Real(3)/8;  _points[ 0](2) = Real(1)/8;  _weights[ 0] =  Real(16)/315;
              _points[ 1](0) = Real(1)/8;  _points[ 1](1) = Real(1)/8;  _points[ 1](2) = Real(5)/8;  _weights[ 1] =  Real(16)/315;
              _points[ 2](0) = Real(3)/8;  _points[ 2](1) = Real(1)/8;  _points[ 2](2) = Real(1)/8;  _weights[ 2] =  Real(16)/315;
              _points[ 3](0) = Real(1)/6;  _points[ 3](1) = Real(1)/2;  _points[ 3](2) = Real(1)/6;  _weights[ 3] = -Real(27)/280;
              _points[ 4](0) = Real(3)/8;  _points[ 4](1) = Real(1)/8;  _points[ 4](2) = Real(3)/8;  _weights[ 4] =  Real(16)/315;
              _points[ 5](0) = Real(1)/8;  _points[ 5](1) = Real(3)/8;  _points[ 5](2) = Real(3)/8;  _weights[ 5] =  Real(16)/315;
              _points[ 6](0) = Real(1)/6;  _points[ 6](1) = Real(1)/6;  _points[ 6](2) = Real(1)/6;  _weights[ 6] = -Real(27)/280;
              _points[ 7](0) = Real(1)/6;  _points[ 7](1) = Real(1)/6;  _points[ 7](2) = Real(1)/2;  _weights[ 7] = -Real(27)/280;
              _points[ 8](0) = Real(1)/8;  _points[ 8](1) = Real(1)/8;  _points[ 8](2) = Real(1)/8;  _weights[ 8] =  Real(16)/315;
              _points[ 9](0) = Real(1)/4;  _points[ 9](1) = Real(1)/4;  _points[ 9](2) = Real(1)/4;  _weights[ 9] =    Real(2)/45;
              _points[10](0) = Real(1)/8;  _points[10](1) = Real(5)/8;  _points[10](2) = Real(1)/8;  _weights[10] =  Real(16)/315;
              _points[11](0) = Real(1)/2;  _points[11](1) = Real(1)/6;  _points[11](2) = Real(1)/6;  _weights[11] = -Real(27)/280;
              _points[12](0) = Real(1)/8;  _points[12](1) = Real(1)/8;  _points[12](2) = Real(3)/8;  _weights[12] =  Real(16)/315;
              _points[13](0) = Real(5)/8;  _points[13](1) = Real(1)/8;  _points[13](2) = Real(1)/8;  _weights[13] =  Real(16)/315;
              _points[14](0) = Real(3)/8;  _points[14](1) = Real(3)/8;  _points[14](2) = Real(1)/8;  _weights[14] =  Real(16)/315;
              return;
            }

          case SIXTH:
          case SEVENTH:
            {
              _points.resize(35);
              _weights.resize(35);

              _points[ 0](0) = Real(3)/10;  _points[ 0](1) = Real(1)/10;  _points[ 0](2) = Real(3)/10;  _weights[ 0] = Real(3125)/72576;
              _points[ 1](0) =  Real(1)/6;  _points[ 1](1) =  Real(1)/2;  _points[ 1](2) =  Real(1)/6;  _weights[ 1] =   Real(243)/4480;
              _points[ 2](0) =  Real(1)/6;  _points[ 2](1) =  Real(1)/6;  _points[ 2](2) =  Real(1)/2;  _weights[ 2] =   Real(243)/4480;
              _points[ 3](0) =  Real(1)/2;  _points[ 3](1) = Real(1)/10;  _points[ 3](2) = Real(1)/10;  _weights[ 3] = Real(3125)/72576;
              _points[ 4](0) = Real(3)/10;  _points[ 4](1) = Real(1)/10;  _points[ 4](2) = Real(1)/10;  _weights[ 4] = Real(3125)/72576;
              _points[ 5](0) = Real(3)/10;  _points[ 5](1) = Real(3)/10;  _points[ 5](2) = Real(1)/10;  _weights[ 5] = Real(3125)/72576;
              _points[ 6](0) =  Real(1)/2;  _points[ 6](1) =  Real(1)/6;  _points[ 6](2) =  Real(1)/6;  _weights[ 6] =   Real(243)/4480;
              _points[ 7](0) = Real(3)/10;  _points[ 7](1) = Real(1)/10;  _points[ 7](2) =  Real(1)/2;  _weights[ 7] = Real(3125)/72576;
              _points[ 8](0) = Real(7)/10;  _points[ 8](1) = Real(1)/10;  _points[ 8](2) = Real(1)/10;  _weights[ 8] = Real(3125)/72576;
              _points[ 9](0) =  Real(3)/8;  _points[ 9](1) =  Real(1)/8;  _points[ 9](2) =  Real(1)/8;  _weights[ 9] =  -Real(256)/2835;
              _points[10](0) = Real(3)/10;  _points[10](1) = Real(3)/10;  _points[10](2) = Real(3)/10;  _weights[10] = Real(3125)/72576;
              _points[11](0) = Real(1)/10;  _points[11](1) =  Real(1)/2;  _points[11](2) = Real(3)/10;  _weights[11] = Real(3125)/72576;
              _points[12](0) = Real(1)/10;  _points[12](1) = Real(1)/10;  _points[12](2) = Real(7)/10;  _weights[12] = Real(3125)/72576;
              _points[13](0) =  Real(1)/2;  _points[13](1) = Real(1)/10;  _points[13](2) = Real(3)/10;  _weights[13] = Real(3125)/72576;
              _points[14](0) = Real(1)/10;  _points[14](1) = Real(7)/10;  _points[14](2) = Real(1)/10;  _weights[14] = Real(3125)/72576;
              _points[15](0) = Real(1)/10;  _points[15](1) =  Real(1)/2;  _points[15](2) = Real(1)/10;  _weights[15] = Real(3125)/72576;
              _points[16](0) =  Real(1)/6;  _points[16](1) =  Real(1)/6;  _points[16](2) =  Real(1)/6;  _weights[16] =   Real(243)/4480;
              _points[17](0) =  Real(3)/8;  _points[17](1) =  Real(1)/8;  _points[17](2) =  Real(3)/8;  _weights[17] =  -Real(256)/2835;
              _points[18](0) =  Real(1)/8;  _points[18](1) =  Real(1)/8;  _points[18](2) =  Real(5)/8;  _weights[18] =  -Real(256)/2835;
              _points[19](0) = Real(1)/10;  _points[19](1) = Real(1)/10;  _points[19](2) = Real(3)/10;  _weights[19] = Real(3125)/72576;
              _points[20](0) =  Real(1)/8;  _points[20](1) =  Real(3)/8;  _points[20](2) =  Real(3)/8;  _weights[20] =  -Real(256)/2835;
              _points[21](0) =  Real(5)/8;  _points[21](1) =  Real(1)/8;  _points[21](2) =  Real(1)/8;  _weights[21] =  -Real(256)/2835;
              _points[22](0) =  Real(1)/8;  _points[22](1) =  Real(5)/8;  _points[22](2) =  Real(1)/8;  _weights[22] =  -Real(256)/2835;
              _points[23](0) = Real(1)/10;  _points[23](1) = Real(3)/10;  _points[23](2) = Real(1)/10;  _weights[23] = Real(3125)/72576;
              _points[24](0) =  Real(1)/4;  _points[24](1) =  Real(1)/4;  _points[24](2) =  Real(1)/4;  _weights[24] =     -Real(8)/945;
              _points[25](0) =  Real(1)/8;  _points[25](1) =  Real(1)/8;  _points[25](2) =  Real(3)/8;  _weights[25] =  -Real(256)/2835;
              _points[26](0) =  Real(3)/8;  _points[26](1) =  Real(3)/8;  _points[26](2) =  Real(1)/8;  _weights[26] =  -Real(256)/2835;
              _points[27](0) =  Real(1)/8;  _points[27](1) =  Real(3)/8;  _points[27](2) =  Real(1)/8;  _weights[27] =  -Real(256)/2835;
              _points[28](0) = Real(1)/10;  _points[28](1) = Real(3)/10;  _points[28](2) =  Real(1)/2;  _weights[28] = Real(3125)/72576;
              _points[29](0) = Real(3)/10;  _points[29](1) =  Real(1)/2;  _points[29](2) = Real(1)/10;  _weights[29] = Real(3125)/72576;
              _points[30](0) = Real(1)/10;  _points[30](1) = Real(1)/10;  _points[30](2) =  Real(1)/2;  _weights[30] = Real(3125)/72576;
              _points[31](0) =  Real(1)/2;  _points[31](1) = Real(3)/10;  _points[31](2) = Real(1)/10;  _weights[31] = Real(3125)/72576;
              _points[32](0) =  Real(1)/8;  _points[32](1) =  Real(1)/8;  _points[32](2) =  Real(1)/8;  _weights[32] =  -Real(256)/2835;
              _points[33](0) = Real(1)/10;  _points[33](1) = Real(3)/10;  _points[33](2) = Real(3)/10;  _weights[33] = Real(3125)/72576;
              _points[34](0) = Real(1)/10;  _points[34](1) = Real(1)/10;  _points[34](2) = Real(1)/10;  _weights[34] = Real(3125)/72576;
              return;
            }

          default:
            {
              // Untested above _order=23 but should work...
              gm_rule(get_order()/2, /*dim=*/3);
              return;
            }
          } // end switch (order)
      } // end case TET4, TET10

    default:
      libmesh_error_msg("ERROR: Unsupported element type: " << Utility::enum_to_string(_type));
    } // end switch (_type)
}

} // namespace libMesh
