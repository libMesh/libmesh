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



#ifndef LIBMESH_SYNC_REFINEMENT_FLAGS_H
#define LIBMESH_SYNC_REFINEMENT_FLAGS_H

// Local includes
#include "libmesh/libmesh_config.h"

// Only compile these functions if the user requests AMR support
#ifdef LIBMESH_ENABLE_AMR

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

// C++ includes
#include <vector>

namespace libMesh {

struct SyncRefinementFlags
{
  typedef unsigned char datum;
  typedef Elem::RefinementState (Elem::*get_a_flag)() const;
  typedef void (Elem::*set_a_flag)(const Elem::RefinementState);

  SyncRefinementFlags(MeshBase & _mesh,
                      get_a_flag _getter,
                      set_a_flag _setter) :
    mesh(_mesh), parallel_consistent(true),
    get_flag(_getter), set_flag(_setter) {}

  MeshBase & mesh;
  bool parallel_consistent;
  get_a_flag get_flag;
  set_a_flag set_flag;
  // References to pointers to member functions segfault?
  // get_a_flag & get_flag;
  // set_a_flag & set_flag;

  // Find the refinement flag on each requested element
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & flags) const
  {
    flags.resize(ids.size());

    for (std::size_t i=0; i != ids.size(); ++i)
      {
        // Look for this element in the mesh
        // We'd better find every element we're asked for
        Elem & elem = mesh.elem_ref(ids[i]);

        // Return the element's refinement flag
        flags[i] = (elem.*get_flag)();
      }
  }

  void act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> & flags)
  {
    for (std::size_t i=0; i != ids.size(); ++i)
      {
        Elem & elem = mesh.elem_ref(ids[i]);

        datum old_flag = (elem.*get_flag)();
        datum new_flag = flags[i];

        if (old_flag != new_flag)
          {
            // It's possible for foreign flags to be (temporarily) more
            // conservative than our own, such as when a refinement in
            // one of the foreign processor's elements is mandated by a
            // refinement in one of our neighboring elements it can see
            // which was mandated by a refinement in one of our
            // neighboring elements it can't see
            // libmesh_assert (!(new_flag != Elem::REFINE &&
            //                   old_flag == Elem::REFINE));
            //
            (elem.*set_flag)
              (static_cast<Elem::RefinementState>(new_flag));
            parallel_consistent = false;
          }
      }
  }
};



} // namespace libMesh


#endif // LIBMESH_ENABLE_AMR

#endif // LIBMESH_SYNC_REFINEMENT_FLAGS_H
