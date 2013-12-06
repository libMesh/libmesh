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



#ifndef LIBMESH_CELL_PYRAMID14_H
#define LIBMESH_CELL_PYRAMID14_H

// Local includes
#include "libmesh/cell_pyramid.h"

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p Pyramid14 is an element in 3D composed of 14 nodes, designed
 * to interface with a QUAD9 element on the base and a TRI6 element on
 * each of the triangular faces.  Cubit will generate hybrid meshes
 * with linear pyramids, but as of version 14 will not export
 * quadratic pyramids.  Paraview may support 13-node pyramids, but
 * does not render 14-node pyramids correctly.  So even if this
 * element works in libmesh, we are curently limited in what we can do
 * with it outside the library...
 *
 * The node numbering for the pyramid14 is given below:
   \verbatim
   PYRAMID14:
                       o 4
                     //|\
                    // | \
                   //  |  \
                  //   |   \
              12 o/    |    o 11
                //     |     \
               /o 9    o 10   \
              //       |       \
             //        |        \
          3 o/.......o.|........o 2
           ./       7  |       /
          ./           |      /
         ./            |     /
        ./             |    /
     8 o/       o      |   o 6
      ./        13     |  /
     ./                | /
    ./                 |/
    o--------o---------o
    0        5         1

   \endverbatim
 */

// ------------------------------------------------------------
// Pyramid class definition
class Pyramid14 : public Pyramid
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Pyramid14  (Elem* p=NULL);

  /**
   * @returns 14.
   */
  virtual unsigned int n_nodes() const { return 14; }

  /**
   * @returns \p PRYAMID14
   */
  ElemType type () const { return PYRAMID14; }

  /**
   * FIXME: we don't yet have a refinement pattern for pyramids...
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

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const;

  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  /**
   * Builds a \p QUAD9 or \p TRI6 coincident with face i.
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_side (const unsigned int i,
			    bool proxy) const;

  /**
   * Builds a \p EDGE3 coincident with edge i.
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_edge (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<dof_id_type>& conn) const;

  /**
   * @returns 2 for all edge nodes and 4 for face nodes
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const;

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
							   const unsigned int v) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[5][9];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[8][3];

protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[14];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int,
                          const unsigned int,
                          const unsigned int) const
  { libmesh_error(); return 0.; }

#endif
};



// ------------------------------------------------------------
// Pyramid14 class member functions
inline
Pyramid14::Pyramid14(Elem* p) :
  Pyramid(Pyramid14::n_nodes(), p, _nodelinks_data)
{
}

} // namespace libMesh


#endif // LIBMESH_CELL_PYRAMID14_H
