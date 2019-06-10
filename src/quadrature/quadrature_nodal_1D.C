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

void QNodal::init_1D(const ElemType, unsigned int)
{
  switch (_type)
    {
    case EDGE2:
      {
        // Nodal quadrature on an Edge2 is QTrap
        QTrap rule(/*dim=*/1, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }
    case EDGE3:
      {
        // Nodal quadrature on an Edge3 is QSimpson
        QSimpson rule(/*dim=*/1, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }
    case EDGE4:
      {
        // The 4-point variant of Simpson's rule. The quadrature
        // points are in the same order as the reference element
        // nodes.
        _points = {Point(-1,0.,0.), Point(+1,0.,0.),
                   Point(-Real(1)/3,0.,0.), Point(Real(1)/3,0.,0.)};
        _weights = {0.25, 0.25, 0.75, 0.75};
        return;
      }
    default:
      libmesh_error_msg("Element type not supported:" << Utility::enum_to_string(_type));
    }
}

} // namespace libMesh
