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
#include "libmesh/quadrature_conical.h"

namespace libMesh
{



void QConical::init_2D(const ElemType type_in,
                       unsigned int p)
{
  switch (type_in)
    {
    case TRI3:
    case TRI6:
      {
        this->conical_product_tri(p);
        return;
      } // end case TRI3, TRI6

      //---------------------------------------------
      // Unsupported element type
    default:
      libmesh_error_msg("ERROR: Unsupported element type: " << type_in);
    } // end switch (type_in)
}

} // namespace libMesh
