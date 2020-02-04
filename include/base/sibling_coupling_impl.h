// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_SIBLING_COUPLING_IMPL_H
#define LIBMESH_SIBLING_COUPLING_IMPL_H

// Local Includes
#include "libmesh/sibling_coupling.h"
#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{

template <typename RealType>
void SiblingCouplingTempl<RealType>::operator()
  (const typename MeshBase::const_element_iterator & range_begin,
   const typename MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  LOG_SCOPE("operator()", "SiblingCoupling");

  for (const auto & elem : as_range(range_begin, range_end))
    {
      std::vector<const Elem *> active_siblings;

      const Elem * parent = elem->parent();
      if (!parent)
        continue;

#ifdef LIBMESH_ENABLE_AMR
      parent->active_family_tree(active_siblings);
#endif

      for (const Elem * sibling : active_siblings)
        if (sibling->processor_id() != p)
          coupled_elements.insert
            (std::make_pair(sibling, _dof_coupling));
    }
}


} // namespace libMesh

#endif // LIBMESH_SIBLING_COUPLING_IMPL_H
