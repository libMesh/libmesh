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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe_shim.h"
#include "libmesh/fe_bernstein_impl.h"
#include "libmesh/fe_bernstein_shape_0D_impl.h"
#include "libmesh/fe_bernstein_shape_1D_impl.h"
#include "libmesh/fe_bernstein_shape_2D_impl.h"
#include "libmesh/fe_bernstein_shape_3D_impl.h"

namespace libMesh
{

// ------------------------------------------------------------
// Bernstein-specific implementations, Steffen Petersen 2005

// Anonymous namespace for local helper functions
namespace {

unsigned int bernstein_n_dofs(const ElemType t, const Order o)
{
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
    case EDGE3:
      return (o+1);
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
      {
        if (o == 1)
          return 4;
        else if (o == 2)
          return 8;
        else
          libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for BERNSTEIN FE family!");
      }
    case QUAD9:
      return ((o+1)*(o+1));
    case HEX8:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case HEX20:
      {
        if (o == 1)
          return 8;
        else if (o == 2)
          return 20;
        else
          libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for BERNSTEIN FE family!");
      }
    case HEX27:
      return ((o+1)*(o+1)*(o+1));
    case TRI3:
    case TRISHELL3:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case TRI6:
      return ((o+1)*(o+2)/2);
    case TET4:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case TET10:
      {
        libmesh_assert_less (o, 3);
        return ((o+1)*(o+2)*(o+3)/6);
      }
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for BERNSTEIN FE family!");
    }
} // bernstein_n_dofs()




unsigned int bernstein_n_dofs_at_node(const ElemType t,
                                      const Order o,
                                      const unsigned int n)
{
  switch (t)
    {
    case NODEELEM:
      return 1;

    case EDGE2:
    case EDGE3:
      switch (n)
        {
        case 0:
        case 1:
          return 1;
        case 2:
          libmesh_assert (o>1);
          return (o-1);
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
        }
    case TRI3:
      libmesh_assert_less (n, 3);
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case TRI6:
      switch (n)
        {
        case 0:
        case 1:
        case 2:
          return 1;

        case 3:
        case 4:
        case 5:
          return (o-1);
          // Internal DoFs are associated with the elem, not its nodes
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
        }
    case QUAD4:
      libmesh_assert_less (n, 4);
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case QUAD8:
      libmesh_assert_less (n, 8);
      libmesh_assert_less (o, 3);
      libmesh_fallthrough();
    case QUAD9:
      {
        switch (n)
          {
          case 0:
          case 1:
          case 2:
          case 3:
            return 1;

          case 4:
          case 5:
          case 6:
          case 7:
            return (o-1);

            // Internal DoFs are associated with the elem, not its nodes
          case 8:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD9!");
          }
      }
    case HEX8:
      libmesh_assert_less (n, 8);
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case HEX20:
      libmesh_assert_less (n, 20);
      libmesh_assert_less (o, 3);
      libmesh_fallthrough();
    case HEX27:
      switch (n)
        {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
          return 1;

        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
        case 18:
        case 19:
          return (o-1);

        case 20:
        case 21:
        case 22:
        case 23:
        case 24:
        case 25:
          return ((o-1)*(o-1));

          // Internal DoFs are associated with the elem, not its nodes
        case 26:
          return 0;

        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for HEX27!");
        }
    case TET4:
      libmesh_assert_less (n, 4);
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case TET10:
      libmesh_assert_less (o, 3);
      libmesh_assert_less (n, 10);
      switch (n)
        {
        case 0:
        case 1:
        case 2:
        case 3:
          return 1;

        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
          return (o-1);

        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TET10!");
        }
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for BERNSTEIN FE family!");
    }
} // bernstein_n_dofs_at_node()




unsigned int bernstein_n_dofs_per_elem(const ElemType t, const Order o)
{
  switch (t)
    {
    case NODEELEM:
    case EDGE2:
    case EDGE3:
      return 0;
    case TRI3:
    case TRISHELL3:
    case QUAD4:
    case QUADSHELL4:
      return 0;
    case TRI6:
      return ((o-1)*(o-2)/2);
    case QUAD8:
    case QUADSHELL8:
      if (o <= 2)
        return 0;
      libmesh_fallthrough();
    case QUAD9:
      return ((o-1)*(o-1));
    case HEX8:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case HEX20:
      libmesh_assert_less (o, 3);
      return 0;
    case HEX27:
      return ((o-1)*(o-1)*(o-1));
    case TET4:
      libmesh_assert_less (o, 2);
      libmesh_fallthrough();
    case TET10:
      libmesh_assert_less (o, 3);
      return 0;
    case INVALID_ELEM:
      return 0;
    default:
      libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for BERNSTEIN FE family!");
    }
} // bernstein_n_dofs_per_elem

} // anonymous namespace

// Explicit instantiations for Real
template struct FEShim<0,BERNSTEIN,Real>;
template struct FEShim<1,BERNSTEIN,Real>;
template struct FEShim<2,BERNSTEIN,Real>;
template struct FEShim<3,BERNSTEIN,Real>;

} // namespace libMesh

#endif // LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
