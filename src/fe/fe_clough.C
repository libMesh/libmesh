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

#include "libmesh/fe_shim.h"
#include "libmesh/fe_clough_impl.h"
#include "libmesh/fe_clough_shape_0D_impl.h"
#include "libmesh/fe_clough_shape_3D_impl.h"

namespace libMesh
{

// ------------------------------------------------------------
// Clough-specific implementations

// Anonymous namespace for local helper functions
namespace {



unsigned int clough_n_dofs(const ElemType t, const Order o)
{
  if (t == INVALID_ELEM)
    return 0;

  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      {
        switch (t)
          {
          case TRI6:
            return 9;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
      // Piecewise cubic Clough-Tocher element
    case THIRD:
      {
        switch (t)
          {
          case TRI6:
            return 12;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for CLOUGH FE family!");
    }
} // clough_n_dofs()





unsigned int clough_n_dofs_at_node(const ElemType t,
                                   const Order o,
                                   const unsigned int n)
{
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      {
        switch (t)
          {
            // The 2D Clough-Tocher defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 3;

                case 3:
                case 4:
                case 5:
                  return 0;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
      // The third-order clough shape functions
    case THIRD:
      {
        switch (t)
          {
            // The 2D Clough-Tocher defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 3;

                case 3:
                case 4:
                case 5:
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for CLOUGH FE family!");
    }
} // clough_n_dofs_at_node()




unsigned int clough_n_dofs_per_elem(const ElemType t, const Order o)
{
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      // The third-order Clough-Tocher shape functions
    case THIRD:
      {
        switch (t)
          {
            // The 2D clough defined on a 6-noded triangle
          case TRI6:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for CLOUGH FE family!");
    }
} // clough_n_dofs_per_elem()

} // anonymous

// Explicit instantiations for Real
template struct FEShim<0,CLOUGH,Real>;
template struct FEShim<1,CLOUGH,Real>;
template struct FEShim<2,CLOUGH,Real>;
template struct FEShim<3,CLOUGH,Real>;

} // namespace libMesh
