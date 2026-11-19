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

#ifndef LIBMESH_FE_MONOMIAL_N_DOFS_H
#define LIBMESH_FE_MONOMIAL_N_DOFS_H

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/libmesh.h" // invalid_uint
#include "libmesh/libmesh_device.h"

namespace libMesh
{

namespace FECounts
{

/**
 * \returns The number of Monomial dofs on an element of type \p t at
 * order \p o, or \p invalid_uint for a combination its element type
 * does not fix.
 *
 * monomial_n_dofs() in fe_monomial.C answers the same question and
 * reports those cases as errors; it reads its counts from here, so that
 * device code, which cannot call into the library, gets the same
 * numbers from the same place.
 *
 * Every count this states is the number of monomials of degree o in the
 * element's dimension, C(o+dim, dim) -- checked for every element type at
 * every order up to sixth.  Writing it as that formula would be a change
 * in its own right, so what is here is the same switch, moved.
 */
LIBMESH_DEVICE_INLINE constexpr unsigned int
monomial_n_dofs (const ElemType t, const Order o)
{
  switch (o)
    {

      // constant shape functions
      // no matter what shape there is only one DOF.
    case CONSTANT:
      return (t != INVALID_ELEM) ? 1 : 0;


      // Discontinuous linear shape functions
      // expressed in the monomials.
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

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 3;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
          case C0POLYHEDRON:
            return 4;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }


      // Discontinuous quadratic shape functions
      // expressed in the monomials.
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 3;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 6;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
          case C0POLYHEDRON:
            return 10;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }


      // Discontinuous cubic shape functions
      // expressed in the monomials.
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 4;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 10;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
          case C0POLYHEDRON:
            return 20;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }


      // Discontinuous quartic shape functions
      // expressed in the monomials.
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 5;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 15;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case C0POLYHEDRON:
            return 35;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }


    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return (order+1);

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return (order+1)*(order+2)/2;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case C0POLYHEDRON:
            return (order+1)*(order+2)*(order+3)/6;

          case INVALID_ELEM:
            return 0;

          default:
            return invalid_uint;
          }
      }
    }
}

} // namespace FECounts

} // namespace libMesh

#endif // LIBMESH_FE_MONOMIAL_N_DOFS_H
