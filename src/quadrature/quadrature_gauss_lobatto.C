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


// libMesh includes
#include "libmesh/quadrature_gauss_lobatto.h"
#include "libmesh/enum_quadrature_type.h"

namespace libMesh
{

// See the files:
// quadrature_gauss_lobatto_1D.C
// quadrature_gauss_lobatto_2D.C
// quadrature_gauss_lobatto_3D.C
// for implementation of specific element types.


QGaussLobatto::QGaussLobatto(const unsigned int d,
                             const Order o) : QBase(d,o)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    init(EDGE2);
}


QuadratureType QGaussLobatto::type() const
{
  return QGAUSS_LOBATTO;
}

} // namespace libMesh
