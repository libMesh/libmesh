// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_INF_HEX8_H
#define LIBMESH_CELL_INF_HEX8_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf_hex.h"

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p InfHex8 is an infinite element in 3D composed of 8 nodes.
 * It is numbered like this:
 * \verbatim
 * INFHEX8: 7        6                             z^  / y
 *          o        o    closer to infinity        | /
 *          :        |                              |/
 *          :        |                              +----> x
 *     4    :   5    |
 *      o   :    o   |
 *      |   o....|...o 2
 *      |  .3    |  /
 *      | .      | /
 *      |.       |/       base face
 *      o--------o
 *      0        1
 *
 * \endverbatim
 */
class InfHex8 : public InfHex
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfHex8  (Elem* p=NULL);

  /**
   * @returns 8.  The \p InfHex8 has 8 nodes.
   */
  unsigned int n_nodes() const { return 8; }

  /**
   * @returns \p INFHEX8
   */
  ElemType     type() const { return INFHEX8; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const;

  /**
   * @returns FIRST
   */
  Order        default_order() const { return FIRST; }

  /**
   * Returns a \p QUAD4 built coincident with face 0, an \p INFQUAD4
   * built coincident with faces 1 to 4. Note that the \p AutoPtr<Elem>
   * takes care of freeing memory.
   */
  AutoPtr<Elem> build_side (const unsigned int i,
                            bool proxy) const;

  /**
   * Returns an \p EDGE2 built coincident with edges 0 to 3, an \p INFEDGE2
   * built coincident with edges 4 to 7. Note that the \p AutoPtr<Elem>
   * takes care of freeing memory.
   */
  AutoPtr<Elem> build_edge (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type>& conn) const;

  unsigned int vtk_element_type (const unsigned int) const
  { return 12; }

  /**
   * @returns \p true when this element contains the point
   * \p p.  Customized for infinite elements, since knowledge
   * about the envelope can be helpful.
   */
  bool contains_point (const Point& p, Real tol=TOLERANCE) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[5][4];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[8][2];


protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[8];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int i,
                          const unsigned int j,
                          const unsigned int k) const
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[4][8][8];

#endif

};



// ------------------------------------------------------------
// InfHex8 class member functions
inline
InfHex8::InfHex8(Elem* p) :
  InfHex(InfHex8::n_nodes(), p, _nodelinks_data)
{
}


} // namespace libMesh

#endif  // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_CELL_INF_HEX8_H
