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


#ifndef LIBMESH_FACE_INF_QUAD4_H
#define LIBMESH_FACE_INF_QUAD4_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/face_inf_quad.h"

// C++ includes
#include <cstddef>

namespace libMesh
{



/**
 * The \p INFQUAD4 is an infinite element in 2D composed of 4 nodes.
 * It is numbered like this:
 \verbatim
 2           3
 INFQUAD4: o           o   closer to infinity
 |           |
 |           |
 |           |
 |           |
 |           |
 o-----------o   base side
 0           1
 \endverbatim
*/
// ------------------------------------------------------------
// InfQuad4 class definition
class InfQuad4 : public InfQuad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  InfQuad4 (Elem* p=NULL) :
    InfQuad(InfQuad4::n_nodes(), p, _nodelinks_data) {}

  /**
   * @returns 4
   */
  unsigned int n_nodes() const { return 4; }

  /**
   * @returns \p INFQUAD4
   */
  ElemType type () const { return INFQUAD4; }

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
   * specified edge (== is_node_on_side in 2D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
  { return this->is_node_on_side(n,e); }

  /**
   * @returns \p FIRST
   */
  Order default_order() const { return FIRST; }

  /**
   * Creates and returns an \p Edge2 for the base side, and an \p InfEdge2 for
   * the sides 1, 2.
   */
  AutoPtr<Elem> build_side (const unsigned int i,
                            bool proxy) const;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type>& conn) const;

  //   void tecplot_connectivity(const unsigned int sf,
  //     std::vector<unsigned int>& conn) const;

  //   void vtk_connectivity(const unsigned int sc,
  // std::vector<unsigned int>*conn = NULL) const;

  //   unsigned int vtk_element_type (const unsigned int) const
  //   { return 9; }

  /**
   * @returns \p true when this element contains the point
   * \p p.  Customized for this \p InfQuad4, since knowledge
   * about the envelope can help avoiding slightly more
   * expensive computations.
   */
  bool contains_point (const Point& p, Real tol=TOLERANCE) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[3][2];


protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[4];



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
  static const float _embedding_matrix[2][4][4];

#endif

};


} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_FACE_INF_QUAD4_H
