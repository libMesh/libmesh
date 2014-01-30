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



#ifndef LIBMESH_CELL_PYRAMID5_H
#define LIBMESH_CELL_PYRAMID5_H

// Local includes
#include "libmesh/cell_pyramid.h"

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p Pyramid5 is an element in 3D composed of 5 nodes.
 * It is numbered with a counter-clockwise base like this:
   \verbatim
   PYRAMID5:
             o 4
           //|\
          // | \
         //  |  \
      3 o/...|...o 2
       ./    |  /
      ./     | /
     ./      |/
    o--------o
    0        1

   \endverbatim
 */

// ------------------------------------------------------------
// Pyramid class definition
class Pyramid5 : public Pyramid
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Pyramid5  (Elem* p=NULL);

  /**
   * @returns \p PRYAMID
   */
  ElemType     type () const   { return PYRAMID5; }

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

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const;

  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }

  /**
   * Builds a \p QUAD4 or \p TRI3 built coincident with face i.
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_side (const unsigned int i,
			    bool proxy) const;

  /**
   * Builds a \p EDGE2 built coincident with edge i.
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_edge (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<dof_id_type>& conn) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[5][4];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[8][2];

  /**
   * Specialization for computing the volume of a pyramid.
   */
  virtual Real volume () const;


protected:

  /**
   * Data for links to nodes
   */
  Node* _nodelinks_data[5];



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
// Pyramid5 class member functions
inline
Pyramid5::Pyramid5(Elem* p) :
  Pyramid(Pyramid5::n_nodes(), p, _nodelinks_data)
{
}

} // namespace libMesh


#endif // LIBMESH_CELL_PYRAMID5_H
