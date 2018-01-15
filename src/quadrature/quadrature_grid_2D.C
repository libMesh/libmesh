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
#include "libmesh/quadrature_grid.h"

namespace libMesh
{


void QGrid::init_2D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules

  // We ignore p - the grid rule is just for experimentation

  switch (type_in)
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
        QGrid q1D(1,_order);
        q1D.init(EDGE2);
        tensor_product_quad( q1D );
        return;
      }


      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRISHELL3:
    case TRI6:
      {
        const unsigned int np = (_order + 1)*(_order + 2)/2;
        const Real weight = Real(0.5)/np;
        const Real dx = Real(1)/(_order+1);
        _points.resize(np);
        _weights.resize(np);

        unsigned int pt = 0;
        for (int i = 0; i != _order + 1; ++i)
          {
            for (int j = 0; j != _order + 1 - i; ++j)
              {
                _points[pt](0) = (i+0.5)*dx;
                _points[pt](1) = (j+0.5)*dx;
                _weights[pt] = weight;
                pt++;
              }
          }
        return;
      }

      //---------------------------------------------
      // Unsupported type
    default:
      libmesh_error_msg("Element type not supported!:" << type_in);
    }
#endif
}

} // namespace libMesh
