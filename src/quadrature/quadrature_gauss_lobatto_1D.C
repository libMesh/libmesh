// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_gauss_lobatto.h"

namespace libMesh
{

void QGaussLobatto::init_1D(const ElemType,
                            unsigned int p)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules
  switch(_order + 2*p)
    {
      // Since Gauss-Lobatto rules must include the endpoints of the
      // domain, there is no 1-point rule.  The two-point
      // Gauss-Lobatto rule is equivalent to the trapezoidal rule.
    case CONSTANT:
    case FIRST:
      {
        _points.resize (2);
        _weights.resize(2);

        _points[0](0) = -1.0L;
        _points[1]    = -_points[0];

        _weights[0]   = 1.;
        _weights[1]   = _weights[0];

        return;
      }

      // The three-point Gauss-Lobatto rule is equivalent to Simpsons' rule.
      // It can integrate cubic polynomials exactly.
    case SECOND:
    case THIRD:
      {
        _points.resize (3);
        _weights.resize(3);

        _points[0](0) = -1.0L;
        _points[1]    = 0.0L;
        _points[2]    = -_points[0];

        _weights[0]   = 1.0L / 3.0L;
        _weights[1]   = 4.0L / 3.0L;
        _weights[2]   = _weights[0];
        return;
      }

      // The four-point Gauss-Lobatto rule can integrate 2*4-3 =
      // 5th-order polynomials exactly.
    case FOURTH:
    case FIFTH:
      {
        _points.resize (4);
        _weights.resize(4);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -std::sqrt(1.0L/5.0L);
        _points[ 2]    = -_points[1];
        _points[ 3]    = -_points[0];

        _weights[ 0]   = 1.0L/6.0L;
        _weights[ 1]   = 5.0L/6.0L;
        _weights[ 2]   = _weights[1];
        _weights[ 3]   = _weights[0];

        return;
      }

      // The five-point Gauss-Lobatto rule can integrate 2*5-3 =
      // 7th-order polynomials exactly.
    case SIXTH:
    case SEVENTH:
      {
        _points.resize (5);
        _weights.resize(5);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -std::sqrt(3.0L/7.0L);
        _points[ 2](0) = 0.;
        _points[ 3]    = -_points[1];
        _points[ 4]    = -_points[0];

        _weights[ 0]   = 1.0L/10.0L;
        _weights[ 1]   = 49.0L/90.0L;
        _weights[ 2]   = 32.0L/45.0L;
        _weights[ 3]   = _weights[1];
        _weights[ 4]   = _weights[0];

        return;
      }

    default:
      libmesh_error_msg("Quadrature rule " << _order << " not supported!");
    }
}

} // namespace libMesh
