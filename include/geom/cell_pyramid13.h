// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_PYRAMID13_H
#define LIBMESH_CELL_PYRAMID13_H

// Local includes
#include "libmesh/cell_pyramid.h"

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p Pyramid13 is an element in 3D composed of 13 nodes, designed
 * to interface with a QUAD8 element on the base and a TRI6 element on
 * each of the triangular faces.  Cubit will generate hybrid meshes
 * with linear pyramids, but as of version 14 will not export
 * quadratic pyramids.  Paraview should support 13-node pyramids...
 *
 * \author John W. Peterson
 * \date 2014
 *
 * The node numbering for the pyramid13 is given below:
 * \verbatim
 * PYRAMID13:
 *                     o 4
 *                   //|\
 *                  // | \
 *                 //  |  \
 *                //   |   \
 *            12 o/    |    o 11
 *              //     |     \
 *             /o 9    o 10   \
 *            //       |       \
 *           //        |        \
 *        3 o/.......o.|........o 2
 *         ./       7  |       /
 *        ./           |      /
 *       ./            |     /
 *      ./             |    /
 *   8 o/              |   o 6
 *    ./               |  /
 *   ./                | /
 *  ./                 |/
 *  o--------o---------o
 *  0        5         1
 *
 * \endverbatim
 */
class Pyramid13 libmesh_final : public Pyramid
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Pyramid13 (Elem * p=libmesh_nullptr) :
    Pyramid(Pyramid13::n_nodes(), p, _nodelinks_data)
  {}

  /**
   * @returns 13.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 13; }

  /**
   * @returns \p PRYAMID13
   */
  virtual ElemType type () const libmesh_override { return PYRAMID13; }

  /**
   * FIXME: we don't yet have a refinement pattern for pyramids...
   * @returns 1
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 1; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override;

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const libmesh_override;

  /**
   * @returns SECOND
   */
  virtual Order default_order() const libmesh_override { return SECOND; }

  /**
   * Builds a \p QUAD8 or \p TRI6 coincident with face i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  virtual UniquePtr<Elem> build_side (const unsigned int i,
                                      bool proxy) const libmesh_override;

  /**
   * Builds a \p EDGE3 coincident with edge i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  virtual UniquePtr<Elem> build_edge (const unsigned int i) const libmesh_override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const libmesh_override;

  /**
   * @returns 2 for all edge nodes
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const libmesh_override;

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const libmesh_override;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[5][8];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[8][3];

  /**
   * Specialization for computing the volume of a Pyramid13.
   */
  virtual Real volume () const libmesh_override;

protected:

  /**
   * Data for links to nodes
   */
  Node * _nodelinks_data[13];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int,
                                  const unsigned int,
                                  const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh


#endif // LIBMESH_CELL_PYRAMID13_H
