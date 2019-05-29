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
#include "libmesh/quadrature_trap.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

void QTrap::init_3D(const ElemType, unsigned int)
{
#if LIBMESH_DIM == 3

  //-----------------------------------------------------------------------
  // 3D quadrature rules
  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
        // We compute the 3D quadrature rule as a tensor
        // product of the 1D quadrature rule.
        QTrap q1D(1);
        q1D.init(EDGE2);

        tensor_product_hex( q1D );

        return;
      }



      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
      {
        _points.resize(4);
        _weights.resize(4);

        _points[0](0) = 0.;
        _points[0](1) = 0.;
        _points[0](2) = 0.;

        _points[1](0) = 1.;
        _points[1](1) = 0.;
        _points[1](2) = 0.;

        _points[2](0) = 0.;
        _points[2](1) = 1.;
        _points[2](2) = 0.;

        _points[3](0) = 0.;
        _points[3](1) = 0.;
        _points[3](2) = 1.;



        _weights[0] = .0416666666666666666666666666666666666666666667;
        _weights[1] = _weights[0];
        _weights[2] = _weights[0];
        _weights[3] = _weights[0];

        return;
      }



      //---------------------------------------------
      // Prism quadrature rules
    case PRISM6:
    case PRISM15:
    case PRISM18:
      {
        // We compute the 3D quadrature rule as a tensor
        // product of the 1D quadrature rule and a 2D
        // triangle quadrature rule

        QTrap q1D(1);
        QTrap q2D(2);

        // Initialize
        q1D.init(EDGE2);
        q2D.init(TRI3);

        tensor_product_prism(q1D, q2D);

        return;
      }


      //---------------------------------------------
      // Unsupported type
    default:
      libmesh_error_msg("ERROR: Unsupported type: " << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
