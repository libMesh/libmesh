// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_simpson.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

void QSimpson::init_3D(const ElemType, unsigned int)
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
        QSimpson q1D(1);
        q1D.init(EDGE2);

        tensor_product_hex( q1D );

        return;
      }



      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
    case TET14:
      {
        // This rule is created by combining 8 subtets
        // which use the trapezoidal rule.  The weights
        // may seem a bit odd, but they are correct,
        // and should add up to 1/6, the volume of the
        // reference tet.  The points of this rule are
        // at the nodal points of the TET10, allowing
        // you to generate diagonal element stiffness
        // matrices when using quadratic elements.
        // It should be able to integrate something
        // better than linears, but I'm not sure how
        // high.

        _points.resize(10);
        _weights.resize(10);

        _points[0](0) = 0.;   _points[5](0) = .5;
        _points[0](1) = 0.;   _points[5](1) = .5;
        _points[0](2) = 0.;   _points[5](2) = 0.;

        _points[1](0) = 1.;   _points[6](0) = 0.;
        _points[1](1) = 0.;   _points[6](1) = .5;
        _points[1](2) = 0.;   _points[6](2) = 0.;

        _points[2](0) = 0.;   _points[7](0) = 0.;
        _points[2](1) = 1.;   _points[7](1) = 0.;
        _points[2](2) = 0.;   _points[7](2) = .5;

        _points[3](0) = 0.;   _points[8](0) = .5;
        _points[3](1) = 0.;   _points[8](1) = 0.;
        _points[3](2) = 1.;   _points[8](2) = .5;

        _points[4](0) = .5;   _points[9](0) = 0.;
        _points[4](1) = 0.;   _points[9](1) = .5;
        _points[4](2) = 0.;   _points[9](2) = .5;


        _weights[0] = Real(1)/192;
        _weights[1] = _weights[0];
        _weights[2] = _weights[0];
        _weights[3] = _weights[0];

        _weights[4] = Real(14)/576;
        _weights[5] = _weights[4];
        _weights[6] = _weights[4];
        _weights[7] = _weights[4];
        _weights[8] = _weights[4];
        _weights[9] = _weights[4];

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

        QSimpson q1D(1);
        QSimpson q2D(2);

        // Initialize
        q1D.init(EDGE2);
        q2D.init(TRI3);

        tensor_product_prism(q1D, q2D);

        return;
      }


      //---------------------------------------------
      // Pyramid quadrature rules
    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
      {
        _points.resize(14);
        _weights.resize(14);

        _points[0](0) = -1.;
        _points[0](1) = -1.;
        _points[0](2) = 0.;

        _points[1](0) = 1.;
        _points[1](1) = -1.;
        _points[1](2) = 0.;

        _points[2](0) = 1.;
        _points[2](1) = 1.;
        _points[2](2) = 0.;

        _points[3](0) = -1.;
        _points[3](1) = 1.;
        _points[3](2) = 0.;

        _points[4](0) = 0.;
        _points[4](1) = 0.;
        _points[4](2) = 1.;

        _points[5](0) = 0.;
        _points[5](1) = -1.;
        _points[5](2) = 0.;

        _points[6](0) = 1.;
        _points[6](1) = 0.;
        _points[6](2) = 0.;

        _points[7](0) = 0.;
        _points[7](1) = 1.;
        _points[7](2) = 0.;

        _points[8](0) = -1.;
        _points[8](1) = 0.;
        _points[8](2) = 0.;

        _points[9](0) = 0.;
        _points[9](1) = -0.5;
        _points[9](2) = 0.5;

        _points[10](0) = 0.5;
        _points[10](1) = 0.;
        _points[10](2) = 0.5;

        _points[11](0) = 0.;
        _points[11](1) = 0.5;
        _points[11](2) = 0.5;

        _points[12](0) = -0.5;
        _points[12](1) = 0.;
        _points[12](2) = 0.5;

        _points[13](0) = 0.;
        _points[13](1) = 0.;
        _points[13](2) = 0.;

        // These are of dubious value since we can't integrate on the
        // vertex where the mapping Jacobian is ill-defined, and even
        // if we could it looks like we'd need negative weight at that
        // vertex to give us exact integrals of both z and z^2.  So I
        // punt and just use QTrap weights.
        _weights[0] = 1/Real(4);
        _weights[1] = _weights[0];
        _weights[2] = _weights[0];
        _weights[3] = _weights[0];
        _weights[4] = 1/Real(3);
        _weights[5] = 0;
        _weights[6] = 0;
        _weights[7] = 0;
        _weights[8] = 0;
        _weights[9] = 0;
        _weights[10] = 0;
        _weights[11] = 0;
        _weights[12] = 0;
        _weights[13] = 0;

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
