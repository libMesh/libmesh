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

#ifndef LIBMESH_SIDES_TO_ELEM_MAP_H
#define LIBMESH_SIDES_TO_ELEM_MAP_H

// libMesh includes
#include "libmesh/id_types.h" // dof_id_type
#include "libmesh/hashword.h" // Utility::hashword()

// C++ includes
#include <unordered_map>
#include <vector>

namespace libMesh
{

// Forward declarations
class MeshBase;
class Elem;

namespace MeshTools
{

/**
 * This class implements a generalization of the
 * MeshTools::build_nodes_to_elem_map() function, but rather than
 * being a standalone function taking a std::unordered_map argument,
 * we define a class with APIs for both building the underlying
 * map data structure and querying its contents.  This way, we avoid
 * issues with hash collisions resulting from the Elem::low_order_key()
 * function (which could occur with as few as 100k sides when using a
 * 32-bit dof_id_type) in addition to being able to control exactly
 * how the user will interact with the class.
 *
 * \author John W. Peterson
 * \date 2024
 */
class SidesToElemMap
{
public:
  /**
   * Default constructor/destructor
   */
  SidesToElemMap();
  ~SidesToElemMap();

  /**
   * Static build function. Most users will simply call this
   * function to construct the SidesToElemMap object for their
   * Mesh.
   */
  static SidesToElemMap build(const MeshBase & mesh);

  /**
   * Typedef for the iterator type returned by the
   * SidesToeElemMap::get_connected_elems() function.
   */
  typedef std::vector<const Elem *>::const_iterator ElemIter;

  /**
   * Return an iterator pair defining the range of Elems connected to "side"
   * of "elem". The returned list of Elems will also include "elem"
   * itself. Throws an error if the requested (elem, side) combination
   * is not found in the map (this should not happen because every
   * side in the map should be attached to at least one Elem).
   */
  std::pair<ElemIter, ElemIter>
  get_connected_elems(const Elem * elem, unsigned int side) const;

private:

  /**
   * Convenient typedefs for working with std::unordered_map
   */
  typedef std::vector<dof_id_type> Key;
  typedef std::vector<const Elem *> Value;

  struct HashFunction
  {
  public:
    /**
     * The "Hash" template argument.  We just use the same
     * Utility::hashword() function that is used by the various
     * Elem::compute_key() implementations, although here we can
     * assume the sorting has already been done for us. This will be a
     * 32-bit hash when dof_id_type is 32-bits and a 64-bit hash when
     * dof_id_type is 64 bits.
     */
    inline
    std::size_t operator()(const Key & vertex_ids) const
    {
      return cast_int<std::size_t>(Utility::hashword(vertex_ids));
    }

    /**
     * The "KeyEqual" template argument. Calls std::vector::operator==()
     */
    inline
    bool operator()(const Key & lhs, const Key & rhs) const
    {
      return lhs == rhs;
    }
  };

  /**
   * Map from (sorted list of side vertex ids) -> (Elems touching that side)
   * A "side" is uniquely defined by the sorted list of its vertex ids, so this
   * is a good choice to use as a key.
   */
  std::unordered_map<Key, Value, /*Hash*/HashFunction, /*KeyEqual*/HashFunction> _sides_to_elem_map;

  /**
   * Construct a sorted list of vertex ids for the input (elem, side)
   * pair. This is needed in a couple different places, so it makes
   * sense to factor it out. The output array is passed as a reference
   * to facilitate reuse over reallocation.
   */
  void get_sorted_vertex_ids(
    const Elem * elem,
    unsigned int side,
    std::vector<dof_id_type> & sorted_vertex_ids) const;
};

} // namespace MeshTools

} // namespace libMesh

#endif // LIBMESH_SIDES_TO_ELEM_MAP_H
