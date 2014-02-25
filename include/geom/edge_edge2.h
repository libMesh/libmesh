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



#ifndef LIBMESH_EDGE_EDGE2_H
#define LIBMESH_EDGE_EDGE2_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/edge.h"

// C++ includes
#include <cstddef>

namespace libMesh
{


/**
 * The \p Edge2 is an element in 1D composed of 2 nodes. It is numbered
 * like this:
 *
 \verbatim
 EDGE2: o--------o
 0        1
 \endverbatim
*/

// ------------------------------------------------------------
// Edge class definition
class Edge2 : public Edge
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Edge2 (Elem* p=NULL) :
    Edge(Edge2::n_nodes(), p, _nodelinks_data) {}

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
   * specified edge (i.e. "returns true" in 1D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const;

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const { return true; }

  /**
   * @returns true iff the Lagrange shape functions on this element
   * are linear
   */
  virtual bool is_linear () const { return true; }

  /**
   * @returns \p EDGE2
   */
  ElemType type()  const { return EDGE2; }

  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type>& conn) const;

  /**
   * An optimized method for computing the length of a 2-node edge.
   */
  virtual Real volume () const;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.  This is a finite element.
   */
  bool infinite () const { return false; }

#endif


protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[2];



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
  static const float _embedding_matrix[2][2][2];

#endif
};

} // namespace libMesh


#endif // LIBMESH_EDGE_EDGE2_H
