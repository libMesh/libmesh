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



#ifndef LIBMESH_CELL_INF_HEX_H
#define LIBMESH_CELL_INF_HEX_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf.h"

namespace libMesh
{




/**
 * The \p InfHex is an element in 3D with 5 sides.
 * The \f$ 6^{th} \f$ side is theoretically located at infinity,
 * and therefore not accounted for.
 * However, one could say that the \f$ 6^{th} \f$ side actually
 * @e does exist in the mesh, since the outer nodes are located
 * at a specific distance from the mesh origin (and therefore
 * define a side).  Still, this face is not to be used!
 */
class InfHex : public InfCell
{
public:

  /**
   * Default infinite brick element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  InfHex(const unsigned int nn, Elem* p, Node** nodelinkdata);

  //   /**
  //    * @returns 4 for the base \p s=0 and 2 for side faces.
  //    */
  //   unsigned int n_children_per_side(const unsigned int s) const;

  /**
   * @returns 5.  Infinite elements have one side less
   * than their conventional counterparts, since one
   * side is supposed to be located at infinity.
   */
  unsigned int n_sides() const { return 5; }

  /**
   * @returns 8.  All infinite hexahedrals (in our
   * setting) have 8 vertices.
   */
  unsigned int n_vertices() const { return 8; }

  /**
   * @returns 8.  All infinite hexahedrals have 8 edges,
   * 4 lying in the base, and 4 perpendicular to the base.
   */
  unsigned int n_edges() const { return 8; }

  /**
   * @returns 5.  All hexahedrals have 5 faces.
   */
  unsigned int n_faces() const { return 5; }

  /**
   * @returns 4
   */
  unsigned int n_children() const { return 4; }

  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const;

  /*
   * @returns true iff the specified edge is on the specified side
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  dof_id_type key (const unsigned int s) const;

  /**
   * @returns a primitive (4-noded) quad or infquad for
   * face i.
   */
  UniquePtr<Elem> side (const unsigned int i) const;

  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  std::pair<Real, Real> qual_bounds (const ElemQuality q) const;



protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem* _elemlinks_data[6+(LIBMESH_DIM>3)];



  /**
   * For higher-order elements, namely \p InfHex16 and
   * \p InfHex18, the matrices for adjacent vertices
   * of second order nodes are quite similar (apart from
   * the face nodes, which are directly handled by \p InfHex18).
   * Therefore hold this matrix here, so that both can re-use
   * this.  Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.
   */
  static const unsigned short int _second_order_adjacent_vertices[8][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[18];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[18];
};



// ------------------------------------------------------------
// InfHex class member functions
inline
InfHex::InfHex(const unsigned int nn, Elem* p, Node** nodelinkdata) :
  InfCell(nn, InfHex::n_sides(), p, _elemlinks_data, nodelinkdata)
{
}


// inline
// unsigned int InfHex::n_children_per_side(const unsigned int s) const
// {
//   libmesh_assert_less (s, this->n_sides());

//   switch (s)
//   {
//     case 0:
//       // every infinite element has 4 children in the base side
//       return 4;

//     default:
//       // on infinite faces (sides), only 2 children exist
//       //
//       // note that the face at infinity is already caught by the libmesh_assertion
//       return 2;
//   }
// }


} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_CELL_INF_HEX_H
