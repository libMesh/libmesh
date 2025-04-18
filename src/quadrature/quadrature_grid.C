// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// libMesh includes
#include "libmesh/quadrature_grid.h"
#include "libmesh/enum_quadrature_type.h"

namespace libMesh
{

QuadratureType QGrid::type() const
{
  return QGRID;
}

std::unique_ptr<QBase> QGrid::clone() const
{
  return std::make_unique<QGrid>(*this);
}

// See the files:
// quadrature_grid_1D.C
// quadrature_grid_2D.C
// quadrature_grid_3D.C
// for implementation.

}
