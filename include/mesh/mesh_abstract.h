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

#ifndef LIBMESH_MESH_ABSTRACT_H
#define LIBMESH_MESH_ABSTRACT_H

#include "libmesh/parallel_object.h"

#include <set>
#include <map>

namespace libMesh
{
class GhostingFunctorBase;

/**
 * Abstract structure that allows avoidance of the templated MeshBaseTempl object in algebraic entities
 */
class MeshAbstract : public ParallelObject
{
public:
  MeshAbstract (const Parallel::Communicator & comm_in);

  MeshAbstract(const MeshAbstract &) = default;
  MeshAbstract(MeshAbstract &&) = default;

  MeshAbstract & operator=(const MeshAbstract &) = delete;
  MeshAbstract & operator=(MeshAbstract &&) = delete;

  virtual ~MeshAbstract() = default;

  /**
   * Removes a functor which was previously added to the set of
   * ghosting functors.
   */
  virtual void remove_ghosting_functor(GhostingFunctorBase & ghosting_functor) = 0;

  /**
   * Adds a functor which can specify ghosting requirements for use on
   * distributed meshes.  Multiple ghosting functors can be added; any
   * element which is required by any functor will be ghosted.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the Mesh is destructed.
   */
  virtual void add_ghosting_functor(GhostingFunctorBase & ghosting_functor) = 0;

  /**
   * Adds a functor which can specify ghosting requirements for use on
   * distributed meshes.  Multiple ghosting functors can be added; any
   * element which is required by any functor will be ghosted.
   *
   * GhostingFunctor memory when using this method is managed by the
   * shared_ptr mechanism.
   */
  virtual void add_ghosting_functor(std::shared_ptr<GhostingFunctorBase> ghosting_functor) = 0;

  /**
   * \returns The logical dimension of the mesh; i.e. the manifold
   * dimension of the elements in the mesh.  If we ever support
   * multi-dimensional meshes (e.g. hexes and quads in the same mesh)
   * then this will return the largest such dimension.
   */
  unsigned int mesh_dimension () const;

protected:
  /**
   * We cache the dimension of the elements present in the mesh.
   * So, if we have a mesh with 1D and 2D elements, this structure
   * will contain 1 and 2.
   */
  std::set<unsigned char> _elem_dims;

};

} // namespace libMesh

#endif // LIBMESH_MESH_ABSTRACT_H
