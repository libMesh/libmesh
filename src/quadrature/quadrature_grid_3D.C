// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


void QGrid::init_3D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM == 3

  //-----------------------------------------------------------------------
  // 3D quadrature rules

  // We ignore p - the grid rule is just for experimentation
  switch (type_in)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QGrid q1D(1,_order);
	q1D.init(EDGE2);

	tensor_product_hex( q1D );

	return;
      }



      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
      {
        _points.resize((_order+1)*(_order+2)*(_order+3)/6);
        _weights.resize((_order+1)*(_order+2)*(_order+3)/6);

        unsigned int pt = 0;
        for (int i = 0; i != _order + 1; ++i)
          {
            for (int j = 0; j != _order + 1 - i; ++j)
              {
                for (int k = 0; k != _order + 1 - i - j; ++k)
                  {
                    _points[pt](0) = (double)i / (double)_order;
                    _points[pt](1) = (double)j / (double)_order;
                    _points[pt](2) = (double)k / (double)_order;
                    _weights[pt] = 1.0 / (double)(_order+1) /
                      (double)(_order+2) / (double)(_order+3);
                    pt++;
                  }
              }
          }
	return;
      }


      // Prism quadrature rules
    case PRISM6:
    case PRISM15:
    case PRISM18:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule and a 2D
	// triangle quadrature rule

	QGrid q1D(1,_order);
	QGrid q2D(2,_order);

	// Initialize
	q1D.init(EDGE2);
	q2D.init(TRI3);

	tensor_product_prism(q1D, q2D);

	return;
      }



      //---------------------------------------------
      // Pyramid
    case PYRAMID5:
    case PYRAMID14:
      {
        _points.resize((_order+1)*(_order+2)*(_order+3)/6);
        _weights.resize((_order+1)*(_order+2)*(_order+3)/6);

        unsigned int pt = 0;
        for (int k = 0; k != _order + 1; ++k)
          {
            for (int i = 0; i != _order + 1 - k; ++i)
              {
                for (int j = 0; j != _order + 1 - k; ++j)
                  {
                    _points[pt](0) = 2.0 * (double)i / (double)_order
                      - 1.0 + (double)k / (double)_order;
                    _points[pt](1) = 2.0 * (double)j / (double)_order
                      - 1.0 + (double)k / (double)_order;
                    _points[pt](2) = (double)k / (double)_order;
                    _weights[pt] = 1.0 / (double)(_order+1) /
                      (double)(_order+2) / (double)(_order+3);
                    pt++;
                  }
              }
          }
	return;
      }



      //---------------------------------------------
      // Unsupported type
    default:
      {
	libMesh::err << "ERROR: Unsupported type: " << type_in << std::endl;
	libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}

} // namespace libMesh
