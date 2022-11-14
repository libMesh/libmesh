// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ELEMCORNER_H
#define LIBMESH_ELEMCORNER_H

#include "libmesh/libmesh_common.h"

#if LIBMESH_DIM > 1

#include "libmesh/elem.h"

#include <utility> // std::pair

namespace libMesh
{

/**
 * The \p ElemCorner is a helper class for defining whether or not
 * _something_ is at an element vertex or edge. In this context,
 * we define corner as a vertex in 2D and as a vertex or an edge
 * in 3D.
 *
 * This class is typically used in intersection algorithms, in which
 * it should be interfaced with by using the \p at_corner(),
 * \p at_edge(), and \p at_vertex() methods.
 *
 * Internally, it keeps track of whether or not at one of these
 * "corners" by storing a set of indicies that denote vertex indices.
 * If the first is valid and the second is not, it implies
 * being at the vertex with said index. If both are valid, it implies
 * being at the edge defined by the verticies with said indices.
 */
class ElemCorner : public std::pair<unsigned short, unsigned short>
{
public:
  /**
   * Default constructor: sets entires to invalid (not at vertex or edge)
   */
  ElemCorner()
    : std::pair<unsigned short, unsigned short>(Elem::invalid_vertex,
                                                Elem::invalid_vertex) { }

  /**
   * Constructor with given vertex indices \p v1 and \p v2.
   */
  ElemCorner(const unsigned short v1, const unsigned short v2)
    : std::pair<unsigned short, unsigned short>(v1, v2) { }

  /**
   * @returns true if at the corner (edge or vertex)
   */
  bool at_corner() const { return first != Elem::invalid_vertex; }

  /**
   * @returns true if the vertex index data is invalid
   */
  bool is_invalid() const
  { return first == Elem::invalid_vertex && second == Elem::invalid_vertex; }

  /**
   * @returns true if at a vertex
   */
  bool at_vertex() const
  { return first != Elem::invalid_vertex && second == Elem::invalid_vertex; }
  /**
   * @returns true if at vertex \p v
   */
  bool at_vertex(const unsigned short v) const
  { return first == v && second == Elem::invalid_vertex; }

  /**
   * Invalidates the current state
   */
  void invalidate()
  {
    first = Elem::invalid_vertex;
    second = Elem::invalid_vertex;
  }

  /**
   * @returns The vertex ID when at a vertex
   */
  unsigned short vertex() const
  { libmesh_assert(at_vertex()); return first; }

  /**
   * Prints the current state (at edge, at vertex, not at either)
   */
  std::string print() const;

  /**
   * Sets the "at vertex" state
   */
  void set_vertex(const unsigned short v)
  {
    libmesh_assert_not_equal_to(v, Elem::invalid_vertex);
    first = v;
    second = Elem::invalid_vertex;
  }

#if LIBMESH_DIM > 2
  /**
   * @returns true if at an edge
   */
  bool at_edge() const
  { return first != Elem::invalid_vertex && second != Elem::invalid_vertex; }
  /**
   * @returns true if at the edge defined by vertices \p v1 and \p v2
   */
  bool at_edge(const unsigned short v1, const unsigned short v2) const
  {
    libmesh_assert_not_equal_to(v1, Elem::invalid_vertex);
    libmesh_assert_not_equal_to(v2, Elem::invalid_vertex);
    return (first == v1 && second == v2) || (first == v2 && second == v1);
  }
  /**
   * @returns true if at the edge with index \p e on element \p
   */
  bool at_edge(const Elem & elem, const unsigned short e) const;

  /**
   * @returns The vertices that contain the edge when at an edge
   */
  const std::pair<unsigned short, unsigned short> & edge_vertices() const
  { libmesh_assert(at_edge()); return *this; }

  /**
   * Sets the "at edge" state
   */
  void set_edge(const unsigned short v1, const unsigned short v2)
  {
    libmesh_assert_not_equal_to(v1, Elem::invalid_vertex);
    libmesh_assert_not_equal_to(v2, Elem::invalid_vertex);
    first = v1;
    second = v2;
  }
  /**
   * Sets the "at edge" state
   */
  void set_edge(const std::pair<unsigned short, unsigned short> & vs)
  { set_edge(vs.first, vs.second); }

  /**
   * @returns The edge when at an edge
   */
  std::unique_ptr<const Elem> build_edge(const Elem & elem) const;
#endif

  /**
   * @returns The vertex point when at a vertex
   */
  const Point & vertex_point(const Elem & elem) const
  { return elem.point(vertex()); }

  /**
   * @returns Whether or not the current state (at vertex/edge) is valid
   * for the given \p elem and \p point. \p tol is used as the tolerance
   * to pass to the equality checks.
   *
   * This ONLY checks for validity when at_corner().
   */
  bool is_valid(const Elem & elem, const Point & point, const Real tol = TOLERANCE) const;
};

} // namespace libMesh

std::ostream & operator<<(std::ostream & os, const libMesh::ElemCorner & elem_corner);

#endif // LIBMESH_DIM > 1

#endif // LIBMESH_ELEMCORNER_H
