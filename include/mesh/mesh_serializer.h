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


#ifndef LIBMESH_MESH_SERIALIZER_H
#define LIBMESH_MESH_SERIALIZER_H

// Local includes

// C++ includes


namespace libMesh
{
// Forward declarations
class MeshBase;

/**
 * Temporarily serialize a DistributedMesh for non-distributed-mesh
 * capable code paths.  A distributed mesh is allgathered by the
 * MeshSerializer constructor if need_serial is true, in which case
 * remote elements are deleted again by the destructor.
 *
 * Serialization to processor 0 alone (e.g. for serial output from
 * that processor) can be selected in the constructor.
 *
 * If allow_remote_element_removal() is set to true, that will also be
 * temporarily disabled by the serializer, to be reenabled after
 * serializer destruction if so; this allows prepare_for_use() to be
 * called safely from within serialized code.
 *
 * If a mesh is explicitly distributed by a `delete_remote_elements()`
 * call within serialized code, or if allow_remote_element_removal()
 * is explicitly set to true within serialized code, the behavior is
 * undefined.
 *
 * \author Roy Stogner
 * \date 2011-2022
 * \brief Temporarily serializes a DistributedMesh for output.
 */
class MeshSerializer
{
public:
  MeshSerializer(MeshBase & mesh, bool need_serial = true, bool serial_only_needed_on_proc_0 = false);

  ~MeshSerializer();

private:
  /*
   * The mesh which should remain serialized while this serializer
   * exists
   */
  MeshBase & _mesh;

  /*
   * Whether to delete remote elements on the mesh (returning it to a
   * distributed state) when the serializer is destroyed
   */
  bool reparallelize;

  /*
   * Whether to again allow `prepare_for_use` to delete remote
   * elements on the mesh (returning it to a distributed state) after
   * the serializer is destroyed
   */
  bool resume_allow_remote_element_removal;
};

} // namespace libMesh

#endif // LIBMESH_MESH_SERIALIZER_H
