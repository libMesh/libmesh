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
#include "libmesh/node_elem.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


Order NodeElem::default_order() const
{
  return FIRST;
}



void NodeElem::connectivity(const unsigned int,
                            const IOPackage,
                            std::vector<dof_id_type> &) const
{
  libmesh_not_implemented();
}

#ifdef LIBMESH_ENABLE_AMR

const float NodeElem::_embedding_matrix[1][1][1] =
  {
    // embedding matrix for child 0
    {
      // 0
      {1.0}, // 0
    }
  };

#endif

} // namespace libMesh
