// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_clough.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

// ---------------------------------------
// FEType class members

std::unique_ptr<QBase>
FEType::default_quadrature_rule (const unsigned int dim,
                                 const int extraorder) const
{
  // Clough elements have at least piecewise cubic functions
  if (family == CLOUGH)
    {
      Order o = static_cast<Order>(std::max(static_cast<unsigned int>(this->default_quadrature_order()),
                                            static_cast<unsigned int>(7 + extraorder)));
      return libmesh_make_unique<QClough>(dim, o);
    }

  if (family == SUBDIVISION)
    return libmesh_make_unique<QGauss>(dim, static_cast<Order>(1 + extraorder));

  return libmesh_make_unique<QGauss>(dim, static_cast<Order>(this->default_quadrature_order() + extraorder));
}

} // namespace libMesh
