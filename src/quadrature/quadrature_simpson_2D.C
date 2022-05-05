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

void QSimpson::init_2D(const ElemType, unsigned int)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules
  switch (_type)
    {


      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        // We compute the 2D quadrature rule as a tensor
        // product of the 1D quadrature rule.
        QSimpson q1D(1);
        q1D.init(EDGE2);
        tensor_product_quad( q1D );
        return;
      }


      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRISHELL3:
    case TRI3SUBDIVISION:
    case TRI6:
    case TRI7:
      {
        // The integral of the vertex Lagrange basis functions of the
        // TRI6 is equal to 0, while the integral of the edge Lagrange
        // basis functions is 1/6, so attempting to derive the weights
        // of a nodal quadrature rule using this approach shows that
        // we only require three quadrature points to get a rule which
        // is exact for quadratics.
        //
        // Unfortunately, it is not possible to derive a nodal
        // quadrature rule on the TRI6 which is exact for cubics, so
        // we instead choose the weights such that they are all
        // strictly positive while still integrating linears
        // exactly. This avoids issues with this nodal quadrature rule
        // producing a singular elemental mass matrix.  We use the
        // "optimization" approach described in quadrature_nodal_2D.C
        // to choose the weights.

        _points.resize(6);
        _weights.resize(6);

        _points[0](0) = 0.;
        _points[0](1) = 0.;

        _points[1](0) = 1.;
        _points[1](1) = 0.;

        _points[2](0) = 0.;
        _points[2](1) = 1.;

        _points[3](0) = 0.5;
        _points[3](1) = 0.;

        _points[4](0) = 0.5;
        _points[4](1) = 0.5;

        _points[5](0) = 0.;
        _points[5](1) = 0.5;

        _weights[0] = Real(1)/38; // 0.0263157894736842
        _weights[1] = Real(1)/38;
        _weights[2] = Real(1)/38;
        _weights[3] = Real(8)/57; // 0.140350877192982
        _weights[4] = Real(8)/57;
        _weights[5] = Real(8)/57;

        return;
      }


      //---------------------------------------------
      // Unsupported type
    default:
      libmesh_error_msg("Element type not supported!:" << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
