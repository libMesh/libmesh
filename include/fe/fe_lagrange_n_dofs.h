// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_FE_LAGRANGE_N_DOFS_H
#define LIBMESH_FE_LAGRANGE_N_DOFS_H

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/libmesh.h" // invalid_uint
#include "libmesh/libmesh_device.h"

namespace libMesh
{

namespace FECounts
{

/**
 * \returns The number of Lagrange dofs on an element of type \p t at
 * order \p o, or \p invalid_uint for a combination whose count its
 * element type does not fix -- an order the family does not have, a
 * type that order cannot be built on, or a polytope, whose dofs follow
 * its node count instead.
 *
 * lagrange_n_dofs() in fe_lagrange.C answers the same question and
 * reports those cases as errors; it reads its counts from here, so that
 * device code, which cannot call into the library, gets the same
 * numbers from the same place.
 */
LIBMESH_DEVICE_INLINE constexpr unsigned int
lagrange_n_dofs (const ElemType t, const Order o)
{
  switch (o)
    {
      // lagrange can only be constant on a single node
    case CONSTANT:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          default:
            return invalid_uint;
          }
      }

      // linear Lagrange shape functions
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 2;

          case TRI3:
          case TRISHELL3:
          case TRI3SUBDIVISION:
          case TRI6:
          case TRI7:
            return 3;

          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 4;

          case TET4:
          case TET10:
          case TET14:
            return 4;

          case HEX8:
          case HEX20:
          case HEX27:
            return 8;

          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            return 6;

          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
            return 5;

          case INVALID_ELEM:
            return 0;

          case C0POLYGON:
          case C0POLYHEDRON:
            // A polytope's dof count follows from its node count, which
            // its element type does not fix
            return invalid_uint;

          default:
            return invalid_uint;
          }
      }


      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE3:
            return 3;

          case TRI6:
          case TRI7:
            return 6;

          case QUAD8:
          case QUADSHELL8:
            return 8;

          case QUAD9:
          case QUADSHELL9:
            return 9;

          case TET10:
          case TET14:
            return 10;

          case HEX20:
            return 20;

          case HEX27:
            return 27;

          case PRISM15:
            return 15;

          case PRISM18:
          case PRISM20:
          case PRISM21:
            return 18;

          case PYRAMID13:
            return 13;

          case PYRAMID14:
          case PYRAMID18:
            return 14;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }

    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE4:
            return 4;

          case PRISM20:
            return 20;

          case PRISM21:
            return 21;

          case PYRAMID18:
            return 18;

          case TRI7:
            return 7;

          case TET14:
            return 14;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }

    default:
      return invalid_uint;
    }
}

} // namespace FECounts

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_N_DOFS_H
