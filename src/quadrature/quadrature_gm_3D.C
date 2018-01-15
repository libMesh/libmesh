// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{



void QGrundmann_Moller::init_3D(const ElemType type_in,
                                unsigned int p)
{
  // Nearly all GM rules contain negative weights, so if you are not
  // allowing rules with negative weights, we cannot continue!
  if (!allow_rules_with_negative_weights)
    libmesh_error_msg("You requested a Grundmann-Moller rule but\n"        \
                      << "are not allowing rules with negative weights!\n" \
                      << "Either select a different quadrature class or\n" \
                      << "set allow_rules_with_negative_weights==true.");

  switch (type_in)
    {
    case TET4:
    case TET10:
      {
        switch(_order + 2*p)
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

              _points[0](0) = 1.L/4.L;
              _points[0](1) = 1.L/4.L;
              _points[0](2) = 1.L/4.L;

              _weights[0] = 1.L/6.L;
              return;
            }

          case SECOND:
          case THIRD:
            {
              _points.resize(5);
              _weights.resize(5);

              _points[0](0) = 1.L/6.L;  _points[0](1) = 1.L/6.L;  _points[0](2) = 1.L/2.L;  _weights[0] =  3.L/40.L;
              _points[1](0) = 1.L/6.L;  _points[1](1) = 1.L/6.L;  _points[1](2) = 1.L/6.L;  _weights[1] =  3.L/40.L;
              _points[2](0) = 1.L/4.L;  _points[2](1) = 1.L/4.L;  _points[2](2) = 1.L/4.L;  _weights[2] = -2.L/15.L;
              _points[3](0) = 1.L/2.L;  _points[3](1) = 1.L/6.L;  _points[3](2) = 1.L/6.L;  _weights[3] =  3.L/40.L;
              _points[4](0) = 1.L/6.L;  _points[4](1) = 1.L/2.L;  _points[4](2) = 1.L/6.L;  _weights[4] =  3.L/40.L;
              return;
            }

          case FOURTH:
          case FIFTH:
            {
              _points.resize(15);
              _weights.resize(15);

              _points[ 0](0) = 1.L/8.L;  _points[ 0](1) = 3.L/8.L;  _points[ 0](2) = 1.L/8.L;  _weights[ 0] =  16.L/315.L;
              _points[ 1](0) = 1.L/8.L;  _points[ 1](1) = 1.L/8.L;  _points[ 1](2) = 5.L/8.L;  _weights[ 1] =  16.L/315.L;
              _points[ 2](0) = 3.L/8.L;  _points[ 2](1) = 1.L/8.L;  _points[ 2](2) = 1.L/8.L;  _weights[ 2] =  16.L/315.L;
              _points[ 3](0) = 1.L/6.L;  _points[ 3](1) = 1.L/2.L;  _points[ 3](2) = 1.L/6.L;  _weights[ 3] = -27.L/280.L;
              _points[ 4](0) = 3.L/8.L;  _points[ 4](1) = 1.L/8.L;  _points[ 4](2) = 3.L/8.L;  _weights[ 4] =  16.L/315.L;
              _points[ 5](0) = 1.L/8.L;  _points[ 5](1) = 3.L/8.L;  _points[ 5](2) = 3.L/8.L;  _weights[ 5] =  16.L/315.L;
              _points[ 6](0) = 1.L/6.L;  _points[ 6](1) = 1.L/6.L;  _points[ 6](2) = 1.L/6.L;  _weights[ 6] = -27.L/280.L;
              _points[ 7](0) = 1.L/6.L;  _points[ 7](1) = 1.L/6.L;  _points[ 7](2) = 1.L/2.L;  _weights[ 7] = -27.L/280.L;
              _points[ 8](0) = 1.L/8.L;  _points[ 8](1) = 1.L/8.L;  _points[ 8](2) = 1.L/8.L;  _weights[ 8] =  16.L/315.L;
              _points[ 9](0) = 1.L/4.L;  _points[ 9](1) = 1.L/4.L;  _points[ 9](2) = 1.L/4.L;  _weights[ 9] =    2.L/45.L;
              _points[10](0) = 1.L/8.L;  _points[10](1) = 5.L/8.L;  _points[10](2) = 1.L/8.L;  _weights[10] =  16.L/315.L;
              _points[11](0) = 1.L/2.L;  _points[11](1) = 1.L/6.L;  _points[11](2) = 1.L/6.L;  _weights[11] = -27.L/280.L;
              _points[12](0) = 1.L/8.L;  _points[12](1) = 1.L/8.L;  _points[12](2) = 3.L/8.L;  _weights[12] =  16.L/315.L;
              _points[13](0) = 5.L/8.L;  _points[13](1) = 1.L/8.L;  _points[13](2) = 1.L/8.L;  _weights[13] =  16.L/315.L;
              _points[14](0) = 3.L/8.L;  _points[14](1) = 3.L/8.L;  _points[14](2) = 1.L/8.L;  _weights[14] =  16.L/315.L;
              return;
            }

          case SIXTH:
          case SEVENTH:
            {
              _points.resize(35);
              _weights.resize(35);

              _points[ 0](0) = 3.L/10.L;  _points[ 0](1) = 1.L/10.L;  _points[ 0](2) = 3.L/10.L;  _weights[ 0] = 3125.L/72576.L;
              _points[ 1](0) =  1.L/6.L;  _points[ 1](1) =  1.L/2.L;  _points[ 1](2) =  1.L/6.L;  _weights[ 1] =   243.L/4480.L;
              _points[ 2](0) =  1.L/6.L;  _points[ 2](1) =  1.L/6.L;  _points[ 2](2) =  1.L/2.L;  _weights[ 2] =   243.L/4480.L;
              _points[ 3](0) =  1.L/2.L;  _points[ 3](1) = 1.L/10.L;  _points[ 3](2) = 1.L/10.L;  _weights[ 3] = 3125.L/72576.L;
              _points[ 4](0) = 3.L/10.L;  _points[ 4](1) = 1.L/10.L;  _points[ 4](2) = 1.L/10.L;  _weights[ 4] = 3125.L/72576.L;
              _points[ 5](0) = 3.L/10.L;  _points[ 5](1) = 3.L/10.L;  _points[ 5](2) = 1.L/10.L;  _weights[ 5] = 3125.L/72576.L;
              _points[ 6](0) =  1.L/2.L;  _points[ 6](1) =  1.L/6.L;  _points[ 6](2) =  1.L/6.L;  _weights[ 6] =   243.L/4480.L;
              _points[ 7](0) = 3.L/10.L;  _points[ 7](1) = 1.L/10.L;  _points[ 7](2) =  1.L/2.L;  _weights[ 7] = 3125.L/72576.L;
              _points[ 8](0) = 7.L/10.L;  _points[ 8](1) = 1.L/10.L;  _points[ 8](2) = 1.L/10.L;  _weights[ 8] = 3125.L/72576.L;
              _points[ 9](0) =  3.L/8.L;  _points[ 9](1) =  1.L/8.L;  _points[ 9](2) =  1.L/8.L;  _weights[ 9] =  -256.L/2835.L;
              _points[10](0) = 3.L/10.L;  _points[10](1) = 3.L/10.L;  _points[10](2) = 3.L/10.L;  _weights[10] = 3125.L/72576.L;
              _points[11](0) = 1.L/10.L;  _points[11](1) =  1.L/2.L;  _points[11](2) = 3.L/10.L;  _weights[11] = 3125.L/72576.L;
              _points[12](0) = 1.L/10.L;  _points[12](1) = 1.L/10.L;  _points[12](2) = 7.L/10.L;  _weights[12] = 3125.L/72576.L;
              _points[13](0) =  1.L/2.L;  _points[13](1) = 1.L/10.L;  _points[13](2) = 3.L/10.L;  _weights[13] = 3125.L/72576.L;
              _points[14](0) = 1.L/10.L;  _points[14](1) = 7.L/10.L;  _points[14](2) = 1.L/10.L;  _weights[14] = 3125.L/72576.L;
              _points[15](0) = 1.L/10.L;  _points[15](1) =  1.L/2.L;  _points[15](2) = 1.L/10.L;  _weights[15] = 3125.L/72576.L;
              _points[16](0) =  1.L/6.L;  _points[16](1) =  1.L/6.L;  _points[16](2) =  1.L/6.L;  _weights[16] =   243.L/4480.L;
              _points[17](0) =  3.L/8.L;  _points[17](1) =  1.L/8.L;  _points[17](2) =  3.L/8.L;  _weights[17] =  -256.L/2835.L;
              _points[18](0) =  1.L/8.L;  _points[18](1) =  1.L/8.L;  _points[18](2) =  5.L/8.L;  _weights[18] =  -256.L/2835.L;
              _points[19](0) = 1.L/10.L;  _points[19](1) = 1.L/10.L;  _points[19](2) = 3.L/10.L;  _weights[19] = 3125.L/72576.L;
              _points[20](0) =  1.L/8.L;  _points[20](1) =  3.L/8.L;  _points[20](2) =  3.L/8.L;  _weights[20] =  -256.L/2835.L;
              _points[21](0) =  5.L/8.L;  _points[21](1) =  1.L/8.L;  _points[21](2) =  1.L/8.L;  _weights[21] =  -256.L/2835.L;
              _points[22](0) =  1.L/8.L;  _points[22](1) =  5.L/8.L;  _points[22](2) =  1.L/8.L;  _weights[22] =  -256.L/2835.L;
              _points[23](0) = 1.L/10.L;  _points[23](1) = 3.L/10.L;  _points[23](2) = 1.L/10.L;  _weights[23] = 3125.L/72576.L;
              _points[24](0) =  1.L/4.L;  _points[24](1) =  1.L/4.L;  _points[24](2) =  1.L/4.L;  _weights[24] =     -8.L/945.L;
              _points[25](0) =  1.L/8.L;  _points[25](1) =  1.L/8.L;  _points[25](2) =  3.L/8.L;  _weights[25] =  -256.L/2835.L;
              _points[26](0) =  3.L/8.L;  _points[26](1) =  3.L/8.L;  _points[26](2) =  1.L/8.L;  _weights[26] =  -256.L/2835.L;
              _points[27](0) =  1.L/8.L;  _points[27](1) =  3.L/8.L;  _points[27](2) =  1.L/8.L;  _weights[27] =  -256.L/2835.L;
              _points[28](0) = 1.L/10.L;  _points[28](1) = 3.L/10.L;  _points[28](2) =  1.L/2.L;  _weights[28] = 3125.L/72576.L;
              _points[29](0) = 3.L/10.L;  _points[29](1) =  1.L/2.L;  _points[29](2) = 1.L/10.L;  _weights[29] = 3125.L/72576.L;
              _points[30](0) = 1.L/10.L;  _points[30](1) = 1.L/10.L;  _points[30](2) =  1.L/2.L;  _weights[30] = 3125.L/72576.L;
              _points[31](0) =  1.L/2.L;  _points[31](1) = 3.L/10.L;  _points[31](2) = 1.L/10.L;  _weights[31] = 3125.L/72576.L;
              _points[32](0) =  1.L/8.L;  _points[32](1) =  1.L/8.L;  _points[32](2) =  1.L/8.L;  _weights[32] =  -256.L/2835.L;
              _points[33](0) = 1.L/10.L;  _points[33](1) = 3.L/10.L;  _points[33](2) = 3.L/10.L;  _weights[33] = 3125.L/72576.L;
              _points[34](0) = 1.L/10.L;  _points[34](1) = 1.L/10.L;  _points[34](2) = 1.L/10.L;  _weights[34] = 3125.L/72576.L;
              return;
            }

          default:
            {
              // Untested above _order=23 but should work...
              gm_rule((_order + 2*p)/2, /*dim=*/3);
              return;
            }
          } // end switch (order)
      } // end case TET4, TET10

    default:
      libmesh_error_msg("ERROR: Unsupported element type: " << type_in);
    } // end switch (type_in)
}

} // namespace libMesh
