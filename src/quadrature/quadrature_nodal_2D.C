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
#include "libmesh/quadrature_nodal.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

void QNodal::init_2D(const ElemType, unsigned int)
{
#if LIBMESH_DIM > 1

  switch (_type)
    {

    case QUAD4:
    case QUADSHELL4:
    case TRI3:
    case TRISHELL3:
      {
        QTrap rule(/*dim=*/2, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    case QUAD8:
    case QUADSHELL8:
      {
        // A rule with 8 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(-1,-1), Point(+1,-1), Point(+1,+1), Point(-1,+1),
            Point(0.,-1), Point(+1,0.), Point(0.,+1), Point(-1,0.)
          };

        // vertex (wv), and edge (we) weights are obtained by:
        // 1.) Requiring that they sum to the reference element volume.
        // 2.) Minimizing the Frobenius norm of the difference between
        //     the resulting nodal quadrature (diagonal) mass matrix
        //     and the true mass matrix for the reference element.
        Real wv = Real(19) / 90;
        Real we = Real(71) / 90;

        _weights = {wv, wv, wv, wv, we, we, we, we};

        return;
      }

    case QUAD9:
    case TRI6:
      {
        QSimpson rule(/*dim=*/2, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    default:
      libmesh_error_msg("Element type not supported!:" << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
