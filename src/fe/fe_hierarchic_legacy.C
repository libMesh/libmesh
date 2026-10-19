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


// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/int_range.h"
#include "libmesh/number_lookups.h"

namespace
{
using namespace libMesh;

// The factor converting the coefficient of a one-dimensional mode of order p from the basis with
// bubbles scaled by 1/p! to the current one: the bubble itself grew by p! s_p
Real mode_ratio(const unsigned int p)
{
  if (p < 2)
    return 1;

  Real factorial = 1;
  for (const auto n : make_range(2u, p + 1))
    factorial *= n;

  return 1 / (factorial * fe_hierarchic_bubble_scaling(p));
}


// The ratio for a triangle shape function: only the edge functions, each a single bubble of the
// edge's order, changed
Real tri_ratio(const unsigned int totalorder, const unsigned int i)
{
  libmesh_assert_less (i, (totalorder+1u)*(totalorder+2u)/2u);

  if (i < 3 || i >= 3u*totalorder)
    return 1;

  // Each edge carries orders 2 through totalorder in turn
  return mode_ratio((i - 3) % (totalorder - 1u) + 2);
}


// The ratio for a hexahedron shape function, whose one-dimensional mode orders cube_indices()
// assigns: vertices, then twelve edges of e = totalorder-1 modes each, then six faces of e*e,
// then the interior
Real hex_ratio(const unsigned int totalorder, const unsigned int i)
{
  const unsigned int e = totalorder - 1u;
  libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u)*(totalorder+1u));

  if (i < 8)
    return 1;

  if (i < 8 + 12*e)
    return mode_ratio((i - 8) % e + 2);

  if (i < 8 + 12*e + 6*e*e)
    {
      // A face mode is a product of two bubbles. The face's orientation decides which axis each
      // lies along, but the product is the same either way
      const unsigned int basisnum = (i - 8 - 12*e) % (e*e);
      return mode_ratio(square_number_row[basisnum] + 2) *
             mode_ratio(square_number_column[basisnum] + 2);
    }

  const unsigned int basisnum = i - 8 - 12*e - 6*e*e;
  return mode_ratio(cube_number_column[basisnum] + 2) *
         mode_ratio(cube_number_row[basisnum] + 2) *
         mode_ratio(cube_number_page[basisnum] + 2);
}


// The ratio for a triangular prism shape function, the product of a triangle function and a
// one-dimensional mode in zeta as prism_indices() lays them out
Real prism_ratio(const unsigned int totalorder, const unsigned int i)
{
  const unsigned int e = totalorder - 1u;
  libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u)*(totalorder+2u)/2u);

  // Vertices
  if (i < 6)
    return 1;

  // Edges 0-2 and 6-8, triangle edges at either end
  if (i < 6 + 3*e)
    return tri_ratio(totalorder, i - 3);

  // Edges 3-5, a zeta bubble at a triangle vertex
  if (i < 6 + 6*e)
    return mode_ratio((i - 6 - 3*e) % e + 2);

  if (i < 6 + 9*e)
    return tri_ratio(totalorder, i - 3 - 6*e);

  // The quadrilateral faces, a triangle edge bubble times a zeta bubble
  if (i < 6 + 9*e + 3*e*e)
    {
      const unsigned int basisnum = (i - 6 - 9*e) % (e*e);
      return mode_ratio(square_number_row[basisnum] + 2) *
             mode_ratio(square_number_column[basisnum] + 2);
    }

  // The triangular faces, which don't contain a bubble
  if (i < 6 + 9*e + 3*e*e + e*(e-1))
    return 1;

  // The interior, a triangle interior function times a zeta bubble
  const unsigned int basisnum = i - 6 - 9*e - 3*e*e - e*(e-1);
  return mode_ratio(prism_number_page[basisnum] + 2);
}


Real scalar_ratio(const Elem & elem,
                  const Order totalorder,
                  const unsigned int i)
{
  switch (elem.type())
    {
    case NODEELEM:
      return 1;

    case EDGE2:
    case EDGE3:
    case EDGE4:
      return mode_ratio(i);

    case TRI3:
    case TRISHELL3:
    case TRI6:
    case TRI7:
      return tri_ratio(totalorder, i);

    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      {
        const auto [i0, i1] = fe_hierarchic_quad_mode_orders(totalorder, i);
        return mode_ratio(i0) * mode_ratio(i1);
      }

    case HEX8:
    case HEX20:
    case HEX27:
      return hex_ratio(totalorder, i);

    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      return prism_ratio(totalorder, i);

    // Edge functions alone carry a bubble
    case TET4:
    case TET10:
    case TET14:
      {
        libmesh_assert_less (i, (totalorder+1u)*(totalorder+2u)*(totalorder+3u)/6u);
        if (i < 4 || i >= 6u*totalorder - 2u)
          return 1;
        return mode_ratio((i - 4) % (totalorder - 1u) + 2);
      }

    default:
      libmesh_error_msg("No HIERARCHIC basis on element type " << Utility::enum_to_string(elem.type()));
    }

  return 1;
}


// The ratio for a SIDE_HIERARCHIC function, a HIERARCHIC function of the side it lives on
Real side_ratio(const Elem & elem,
                const Order totalorder,
                const unsigned int i)
{
  const unsigned int o = totalorder;

  switch (elem.type())
    {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return 1;

    case TRI6:
    case TRI7:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      // Each side is an edge, its functions taken in order 0 through o
      return mode_ratio(i % (o + 1));

    case HEX27:
      {
        // Each side is a quadrilateral. cube_remap() swaps vertex functions with other vertex
        // functions and edge functions with other edge functions of the same order, so the side
        // index gives the orders directly
        const auto [i0, i1] = fe_hierarchic_quad_mode_orders(o, i % ((o+1)*(o+1)));
        return mode_ratio(i0) * mode_ratio(i1);
      }

    case TET14:
      return tri_ratio(o, i % ((o+1)*(o+2)/2));

    case PRISM20:
    case PRISM21:
      {
        const unsigned int dofs_per_quad = (o+1)*(o+1);
        if (i < 3*dofs_per_quad)
          {
            const auto [i0, i1] = fe_hierarchic_quad_mode_orders(o, i % dofs_per_quad);
            return mode_ratio(i0) * mode_ratio(i1);
          }
        return tri_ratio(o, (i - 3*dofs_per_quad) % ((o+1)*(o+2)/2));
      }

    default:
      libmesh_error_msg("No SIDE_HIERARCHIC basis on element type " << Utility::enum_to_string(elem.type()));
    }

  return 1;
}

} // anonymous namespace



namespace libMesh
{

bool fe_hierarchic_bubble_family (const FEFamily family)
{
  switch (family)
    {
    case HIERARCHIC:
    case L2_HIERARCHIC:
    case SIDE_HIERARCHIC:
    case HIERARCHIC_VEC:
    case L2_HIERARCHIC_VEC:
      return true;
    default:
      return false;
    }
}


Real fe_hierarchic_legacy_coefficient_ratio (const FEFamily family,
                                             const Elem & elem,
                                             const Order totalorder,
                                             const unsigned int i)
{
  switch (family)
    {
    case HIERARCHIC:
    case L2_HIERARCHIC:
      return scalar_ratio(elem, totalorder, i);

    case SIDE_HIERARCHIC:
      return side_ratio(elem, totalorder, i);

    // The vector families interleave one scalar function per spatial component
    case HIERARCHIC_VEC:
    case L2_HIERARCHIC_VEC:
      {
        const unsigned int dim = elem.dim();
        return scalar_ratio(elem, totalorder, dim ? i / dim : i);
      }

    default:
      libmesh_error_msg("FE family " << Utility::enum_to_string(family) <<
                        " does not use the HIERARCHIC bubbles");
    }

  return 1;
}

} // namespace libMesh
